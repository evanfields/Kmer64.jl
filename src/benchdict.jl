function benchdict(nvals, dictsize; ntest = 1000)
    d = Dict{UInt64, Nothing}()
    sizehint!(d, dictsize)
    for _ in 1:nvals
        d[rand(UInt64)] = nothing
    end
    tests = rand(UInt64, ntest)
    return @elapsed any(haskey(d, t) for t in tests)
end

function benchdict_latesize(nvals, dictsize; ntest = 1000)
    d = Dict{UInt64, Nothing}()
    for _ in 1:nvals
        d[rand(UInt64)] = nothing
    end
    sizehint!(d, dictsize)
    tests = rand(UInt64, ntest)
    return @elapsed any(haskey(d, t) for t in tests)
end

function benchrobindict(nvals, dictsize; ntest = 1000)
    d = RobinDict{UInt64, Nothing}()
    sizehint!(d, dictsize)
    for _ in 1:nvals
        d[rand(UInt64)] = nothing
    end
    tests = rand(UInt64, ntest)
    return @elapsed any(haskey(d, t) for t in tests)
end