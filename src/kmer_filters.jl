"""Struct to hold UInt128 data of a set of fixed-width kmers; optimized for membership
tests with low hit rates."""
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
            min(2^14, 12 * length(vals))
        )
        sizehint!(kds.lower.dict, dict_size)
        sizehint!(kds.full.dict, dict_size)
        return kds
    end
end

function member(v::UInt128, kds::KmerDataSet)
    return in(v % UInt64, kds.lower) && in(v, kds.full)
end

"""The ordered set of kmers in a DNA sequence."""
function kmer_list(seq::DNASeq, k = 40)
    if BioSequences.hasambiguity(seq)
        error("Cannot construct kmers from a sequence with ambiguity.")
    end
    kmers = Kmer[]
    if length(seq) < 40
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

"""Whether any kmer in `seq` is in `kmer_set`, where `kmer_set` contains just the data
UInt128 of a set of kmers.
!!! `seq` must not contain any ambiguities!"""
function _unambiguous_seq_has_kmer(seq, kmer_set::KmerDataSet, k)
    @debug "In _unambiguous_seq_has_kmer" seq.part # to-do: remove after development! Causes allocation?!
    length(seq) < k && return false
    kmer = Kmer(@view seq[1:k])
    member(kmer.data, kmer_set) && return true
    for next_base in (@view seq[k + 1 : end])
        kmer = next_kmer(kmer, next_base)
        member(kmer.data, kmer_set) && return true
    end
    return false
end

"""Check whether a sequence has a kmer in kmer_set. Do not form/check kmers that would
contain ambiguous bases. For performance, `kmer_set` should contain just the UInt128 data
of the query kmers."""
function _seq_has_kmer(seq::LongDNA{4}, kmer_set::KmerDataSet, k)
    # We cannot form any kmers with ambiguous bases, so we find unambiguous chunks
    # and process each chunk in turn.
    # Plausibly this is not perfect for cache locality: suppose the read is long but
    # unambiguous. Then we walk the read once to confirm it's unambiguous and then again to
    # actually check kmers. Checking for ambiguity right where we form kmers might be better,
    # but definitely messier. Possibly this doesn't matter in practice if seq is a 150bp
    # read and the entire thing easily fits in cache.
    next_ind = 1 # index of the start of the next kmer to check
    while next_ind <= length(seq) - k + 1
        next_ambiguity_ind = findnext(isambiguous, seq, next_ind)
        if isnothing(next_ambiguity_ind) # no more ambiguities
            return _unambiguous_seq_has_kmer((@view seq[next_ind : end]), kmer_set, k)
        elseif next_ambiguity_ind >= next_ind + k
            # the next unambiguous chunk has room to form kmers
            next_chunk_hits = _unambiguous_seq_has_kmer(
                (@view seq[next_ind : (next_ambiguity_ind - 1)]),
                kmer_set,
                k
            )
            next_chunk_hits && return true
        end
        # didn't find a kmer hit, need to keep looking
        next_ind = next_ambiguity_ind + 1
    end
    return false
end

function _prepare_query_kmers(query_fasta_path, k::Int, check_rc::Bool)
    query = open(FASTAReader, query_fasta_path) do reader
        sequence(LongDNA{2}, first(reader))
    end
    if hasambiguity(query)
        error("Query sequence may not have ambiguous bases.")
    end
    query_kmers = kmer_list(query, k)
    if check_rc
        query_kmers = vcat(query_kmers, reverse_complement.(query_kmers))
    end
    return KmerDataSet(k.data for k in query_kmers)
end
function filter_read_file(fastq_path, query_fasta_path; k::Int = 40, check_rc = true)
    kmer_set = _prepare_query_kmers(query_fasta_path, k, check_rc)

    return open(FASTQReader, fastq_path) do reader
        return _filter_hits(reader, kmer_set, k)
    end
end

function filter_paired_reads(reads1_path, reads2_path, query_fasta_path; k::Int = 40, check_rc = true)
    kmer_set = _prepare_query_kmers(query_fasta_path, k, check_rc)
    hits = Vector{Tuple{FASTQRecord, FASTQRecord}}()
    open(FASTQReader, reads1_path) do reader1
        open(FASTQReader, reads2_path) do reader2
            for (rec1, rec2) in zip(reader1, reader2)
                if (
                    _seq_has_kmer(sequence(LongDNA{4}, rec1), kmer_set, k) ||
                    _seq_has_kmer(sequence(LongDNA{4}, rec2), kmer_set, k)
                )
                    push!(hits, (rec1, rec2))
                end
            end
        end
    end
    return hits
end

function _filter_hits(records, kmer_set, k)
    hits = FASTQRecord[]
    for rec in records
        if _seq_has_kmer(sequence(LongDNA{4}, rec), kmer_set, k)
            push!(hits, rec)
        end
    end
    return hits
end
