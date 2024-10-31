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

function chunk_starts(file, n; min_chunk_bytes = 300_000)
end