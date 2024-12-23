"""
    KmerDataSet

Struct to hold UInt128 data of a set of fixed-width kmers; optimized for membership
tests with low hit rates. Uses a two-tier lookup structure to optimize the common
case where most membership tests fail.

Implementation
- `lower`: Set of UInt64s containing truncated kmer representations 
- `full`: Set of UInt128s containing complete kmer data.
"""
struct KmerDataSet
    lower::Set{UInt64}
    full::Set{UInt128}

    function KmerDataSet(vals)
        kds = new(
            Set{UInt64}(v % UInt64 for v in vals),
            Set{UInt128}(vals)
        )
        # sets are Dict backed, big size hint seems to improve haskey performance
        dict_size = max(
            length(vals),
            min(2^14, 20 * length(vals))
        )
        # from sizehint! implementation, this should trigger rehashing
        sizehint!(kds.lower.dict, dict_size)
        sizehint!(kds.full.dict, dict_size)
        return kds
    end
end

const ReadPair = Tuple{FASTQRecord, FASTQRecord}

@inline function member(v::UInt128, kds::KmerDataSet)
    return in(v % UInt64, kds.lower) && in(v, kds.full)
end

"""The ordered set of kmers in a DNA sequence."""
function kmer_list(seq::DNASeq, k = 40)
    if BioSequences.hasambiguity(seq)
        error("Cannot construct kmers from a sequence with ambiguity.")
    end
    kmers = Kmer[]
    if length(seq) < k
        return kmers
    end
    push!(kmers, Kmer(@view seq[1:k]))
    last_ind = k
    while last_ind < length(seq)
        last_ind += 1
        push!(kmers, next_kmer(kmers[end], seq[last_ind]))
    end
    return kmers
end

# really just for development, not performant since it forms the query kmerset
# for every read and kmer_list allocates
function filter_reads(reads, query; k::Int = 40, check_rc = true)
    query_kmers = kmer_list(query, k)
    if check_rc
        query_kmers = vcat(query_kmers, reverse_complement.(query_kmers))
    end
    query_set = Set(query_kmers)
    return filter(reads) do read
        any(in(query_set), kmer_list(read, k))
    end
end

"""
Test whether a DNA sequence contains any kmer present in a KmerDataSet. Process sequence in a
single pass, building kmers incrementally and handling ambiguous bases. Return early on first
match found.

Note: Assumes k ≤ 64 and that kmers in kmer_set were constructed using the same k value.
For longer kmers, behavior is undefined.
"""
function _seq_has_kmer(seq::LongDNA{4}, kmer_set::KmerDataSet, k)
    length(seq) < k && return false
    
    # Initialize variables for building kmers incrementally
    empty_bits = 0x80 - ((k % UInt8) << 0x1)
    relevant_bits_mask = typemax(UInt128) >> empty_bits
    kmer_data = zero(UInt128)
    current_bases = 0
    
    # Process all bases in a single pass

    @inbounds for i in 1:length(seq)
        base = seq[i]
        if isambiguous(base)
            current_bases = 0
            kmer_data = zero(UInt128)
            continue
        end
        
        # Always use _next_data - when current_bases is 0 this is equivalent to
        # starting fresh since kmer_data is zero
        kmer_data = _next_data(kmer_data, base, relevant_bits_mask)
        current_bases += 1
            
        # Only check membership once we have a complete kmer
        if current_bases >= k && member(kmer_data, kmer_set)
            return true
        end
    end
    
    return false
end

"""Extract query kmers from a fasta file, validate they have no ambiguous bases, and
return a KmerDataSet for read filtering. Return a pair
`(forward_KmerDataSet, rc_KmerDataSet)`. If `check_rc`, both element are the
same KmerDataSet which contains data for kmers in both forward and RC orientation."""
function _prepare_query_kmers(query_fasta_path, k::Int, check_rc::Bool)
    fwd_kmers = open(FASTAReader, query_fasta_path) do reader
        seqwise_kmers = map(reader) do rec
            seq = sequence(LongDNA{4}, rec)
            if hasambiguity(seq)
                error("Query sequence may not have ambiguous bases.")
            end
            seq_unambig = LongDNA{2}(seq)
            seq_kmers = kmer_list(seq_unambig, k)
            return seq_kmers
        end
        return reduce(vcat, seqwise_kmers)
    end
    if check_rc
        kmers = vcat(fwd_kmers, reverse_complement.(fwd_kmers))
        kds = KmerDataSet(k.data for k in kmers)
        return (kds, kds)
    else # handle fwd and RC orientations separately
        return (
            KmerDataSet(k.data for k in fwd_kmers),
            KmerDataSet(k.data for k in reverse_complement.(fwd_kmers))
        )
    end
end

"""Filter a file of unpaired reads. Needs updating."""
function filter_read_file(fastq_path, query_fasta_path; k::Int = 40, check_rc = true)
    kmer_set = _prepare_query_kmers(query_fasta_path, k, check_rc)[1]

    return open(FASTQReader, fastq_path) do reader
        return _filter_hits(reader, kmer_set, k)
    end
end

function filter_paired_reads(
    reads1_path,
    reads2_path,
    query_fasta_path,
    out1_path,
    out2_path;
    k::Int = 40,
    check_rc = true
)
    fwd_kmerset, rev_kmerset = _prepare_query_kmers(query_fasta_path, k, check_rc)
    hits = Vector{Tuple{FASTQRecord, FASTQRecord}}()
    open(FASTQReader, reads1_path) do reader1
        open(FASTQReader, reads2_path) do reader2
            for (rec1, rec2) in zip(reader1, reader2)
                if (
                    _seq_has_kmer(sequence(LongDNA{4}, rec1), fwd_kmerset, k) ||
                    _seq_has_kmer(sequence(LongDNA{4}, rec2), rev_kmerset, k)
                )
                    push!(hits, (rec1, rec2))
                end
            end
        end
    end
    return write_fastq_pairs(hits, out1_path, out2_path)
end

"""Write pairs of fastq records to a pair of file paths. If either path is Nothing or empty,
return without writing. Return the number of records written."""
function write_fastq_pairs(pairs, path1, path2)
    if isnothing(path1) || isnothing(path2) || isempty(path1) || isempty(path2)
        return 0
    end
    written = 0
    open(FASTQWriter, path1) do writer1
        open(FASTQWriter, path2) do writer2
            for (rec1, rec2) in pairs
                write(writer1, rec1)
                write(writer2, rec2)
                written += 1
            end
        end
    end
    return written
end


"""Filter an interable of FASTQRecords for those that contain a kmer in kmer_set."""
function _filter_hits(records, kmer_set, k)
    hits = FASTQRecord[]
    for rec in records
        if _seq_has_kmer(sequence(LongDNA{4}, rec), kmer_set, k)
            push!(hits, rec)
        end
    end
    return hits
end
