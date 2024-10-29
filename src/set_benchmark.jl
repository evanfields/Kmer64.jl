struct DoubleSet
    lower::Set{UInt64}
    full::Set{UInt128}

    function DoubleSet(vals)
        return new(
            Set{UInt64}(v % UInt64 for v in vals),
            Set{UInt128}(vals)
        )
    end
end

function member(v::UInt128, ds::DoubleSet)
    return in(v % UInt64, ds.lower) && in(v, ds.full)
end