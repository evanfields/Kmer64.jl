"""Advance a stream to the next record-start position and return this new position.
Exception: if the stream is at position zero and the next character is a '@', ie the first
character of a FASTQ record, return zero."""
function next_rec_start_pos!(io; max_lines = 50)
    init_pos = position(io)
    next_char = peek(io, Char)
    # perhaps we're already at a record start?
    if next_char == '@' && init_pos == 0
        return 0
    end
    for _ in 1:max_lines
        println("Need to read to end of line")
        rest_of_line = readuntil(io, '\n')
        @show rest_of_line
        @show position(io)
        peek(io, Char) == '@' && return position(io)
    end
    error("No @ on a newline found in 50 lines from position $(init_pos)")
end

# todo
function chunk_starts(file, n; min_chunk_bytes = 300_000)
end


"""
    filter_paired_reads_threaded(reads1_path, reads2_path, query_fasta_path, 
                               out1_path, out2_path; k=40, check_rc=true, chunk_size=32)

Parallel version of read pair filtering. Identifies read pairs where either read contains
an exact match to any k-length subsequence from the query sequence.

# Arguments
- `reads1_path`: Path to FASTQ file containing first reads of pairs
- `reads2_path`: Path to FASTQ file containing second reads of pairs
- `query_fasta_path`: Path to FASTA file containing query sequence
- `out1_path`: Output path for filtered first reads
- `out2_path`: Output path for filtered second reads
- `k=40`: Length of kmers to match (must be ≤ 64)
- `check_rc=true`: Whether to also match reverse complements of query kmers
- `chunk_size=32`: Number of read pairs to process in each parallel task

# Returns
Number of read pairs written to output files.

Performance scales with available threads (set via JULIA_NUM_THREADS).
"""
function filter_paired_reads_threaded(
    reads1_path,
    reads2_path,
    query_fasta_path,
    out1_path,
    out2_path;
    k::Int = 40,
    check_rc = true,
    chunk_size = 32
)
    fwd_kmerset, rev_kmerset = _prepare_query_kmers(query_fasta_path, k, check_rc)
    chunkwise_results = Channel{Channel{ReadPair}}(256) # good channel size?
    writer_task = Threads.@spawn write_fastq_pairs(chunkwise_results, out1_path, out2_path)
    open(FASTQReader, reads1_path) do reader1
        open(FASTQReader, reads2_path) do reader2
            # to-do: are we using two threads to read two inputs?
            for chunk in Iterators.partition(zip(reader1, reader2), chunk_size)
                chunk_channel = Channel{ReadPair}(chunk_size)
                Threads.@spawn _filter_chunk!(
                    chunk_channel,
                    chunk,
                    fwd_kmerset,
                    rev_kmerset,
                    k
                )
                put!(chunkwise_results, chunk_channel)
            end
        end
    end
    close(chunkwise_results) # we put! all the chunkwise channels already
    return fetch(writer_task)
end

"""Extract matching read pairs from a ReadPair iterable `read_pairs` and `put!` them in a
Channel, then close the Channel"""
function _filter_chunk!(chunk_channel, read_pairs, fwd_kmerset, rev_kmerset, k)
    for (rec1, rec2) in read_pairs
        if (
            _seq_has_kmer(sequence(LongDNA{4}, rec1), fwd_kmerset, k) ||
            _seq_has_kmer(sequence(LongDNA{4}, rec2), rev_kmerset, k)
        )
            put!(chunk_channel, (rec1, rec2))
        end
    end
    close(chunk_channel)
end

"""Write ReadPairs in a Channel of Channels to a pair of file paths. If either path is
Nothing or empty, return without writing. Return the number of records written."""
function write_fastq_pairs(chunkwise_results::Channel, path1, path2)
    # if we have nothing to write, still consume the channel so it doesn't become blocked
    if isnothing(path1) || isnothing(path2) || isempty(path1) || isempty(path2)
        for chunk_result in chunkwise_results
            for (rec1, rec2) in chunk_result
            end
        end
        return 0
    end
    # we can actually write, do so as results are available in order
    written = 0
    open(FASTQWriter, path1) do writer1
        open(FASTQWriter, path2) do writer2
            for chunk_result in chunkwise_results
                for (rec1, rec2) in chunk_result
                    write(writer1, rec1)
                    write(writer2, rec2)
                    written += 1
                end
            end
        end
    end
    return written
end
