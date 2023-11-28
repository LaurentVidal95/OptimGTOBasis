function filter_dual(x::T) where T
    if typeof(x) <: ForwardDiff.Dual
        return filter_dual(x.value)
    end
    x
end
