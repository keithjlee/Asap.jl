# fixed end forces/moments for frame elements

"""
point load at fraction t from beginning, perpendicular to beam in LCS
"""
function frame2d_point(element::Element, P::Union{Int64, Float64}, t::Float64)
    if !(0. <= t <= 1.)
        error("t must be ∈ [0.0, 1.0]")
    end

    L = element.length
    l1 = L * t
    l2 = L - l1

    FSb = P * l2^2 / L^3 * (3l1 + l2)
    FMb = P * l1 * l2^2 / L^2
    FSe = P * l1^2 / L^3 * (l1 + 3l2)
    FMe = - P * l1^2 * l2 / L^2

    return [Fsb, FMb, FSe, FMb]
end

"""
Uniformly distributed load perpendicular to beam in LCS
"""
function frame2d_udl(element::Element, w::Union{Int64, Float64})
    L = element.length

    FSb = w * L / 2
    FMb = w * L^2 / 12
    FSe = w * L / 2
    FMe = - w * L^2 / 12

    return [FSb, FMb, FSe, FMe]
end

"""
Generic truncated UDL on beam
"""
function frame2d_udl2(element::Element, w::Union{Int64, Float64}, t1::Float64, t2::Float64)
    if !(0. <= t1 < 1.) || !(0. <= t2 < 1.)
        error("t1, t2 must be ∈ [0.0, 1.0)")
    end

    L = element.length
    l1 = L * t1
    l2 = L * t2

    FSb = w * L / 2 * (1 - l1 / L^4 * (2L^3 - 2l1^2*L + l1^3) - l2^3/L^4 * (2L - l2))
    FMb = w * L^2 / 12 * (1 - l1^2 / L^4 * (6L^2 - 8l1*L + 3l1^2) - l2^3/L^4 *(4L - 3l2))
    FSe = w * L / 2 * (1 - l1^3/L^4 * (2L - l1) - l2/L^4 * (2L^3 - 2l2^2*L + l2^3))
    FMe = -w * L^2 / 12 * (1 - l1^3/L^4 * (4L - 3l1) - l2^2/L^4 * (6L^2 - 8l2 * L + 3l2^2))

    return [FSb, FMb, FSe, FMe]
end
