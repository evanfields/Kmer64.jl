"""Partition iterator that @spawn-s its buffer-filling task."""
struct ThreadedPartition{T, S}
    source::T
    chunksize::Int
    buffer::Channel{Vector{S}}
    
    function ThreadedPartition(source, chunksize; items_to_buffer = 32 * chunksize)
        # Type parameters
        T = typeof(source)
        S = eltype(source)
        
        # Create buffered channel for chunks
        buffer_size = max(1, floor(Int, items_to_buffer / chunksize))
        buffer = Channel{Vector{S}}(buffer_size)
        
        # Create iterator
        tp = new{T,S}(source, chunksize, buffer)
        
        # Start background task to fill buffer
        Threads.@spawn begin
            try
                chunk = Vector{S}()
                sizehint!(chunk, chunksize)
                for item in source
                    push!(chunk, item)
                    if length(chunk) >= chunksize
                        put!(buffer, chunk)
                        chunk = Vector{S}()
                        sizehint!(chunk, chunksize)
                    end
                end
                # Handle any remaining items
                if !isempty(chunk)
                    put!(buffer, chunk)
                end
            catch e
                @error "Error in ThreadedPartition" exception=(e, catch_backtrace())
            finally
                close(buffer)
            end
        end
        
        return tp
    end
end

iterate(tp::ThreadedPartition, state = nothing) = iterate(tp.buffer, state)

function try_to_break(path)
    open(path) do io
        stream = FASTX.Automa.TranscodingStreams.NoopStream(io)
        tasks = [Threads.@spawn begin
            while !eof(stream)
                read(stream, 1024)
            end
        end for _ in 1:6]  # Same number of threads that triggered the error
        fetch.(tasks)
    end
end