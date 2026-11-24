# The Pratt truss shared by the optimization walkthroughs (kN, m): 12 m span,
# 2 m deep, 6 panels, 20 kN at each interior bottom node, planar in XY.
using Asap
function pratt_truss(; L = 12.0, H = 2.0, n = 6, A = 2.0e-3, P = 20.0)
    steel = Material(200e6, 77e6, 7.85, 0.3)
    bar = Section(steel, A)
    b = [Node([L * i / n, 0.0, 0.0], :free, :bottom) for i in 0:n]
    t = [Node([L * i / n, H, 0.0], :free, :top) for i in 1:n-1]
    b[1] = Node([0.0, 0.0, 0.0], Bool[0, 0, 0, 1, 1, 1], :pin)
    b[end] = Node([L, 0.0, 0.0], Bool[1, 0, 0, 1, 1, 1], :roller)
    nodes = vcat(b, t); planarize!(nodes, :XY)
    T(a, c, id) = TrussElement(a, c, bar, id)
    els = AbstractElement{Float64}[[T(b[i], b[i+1], :bottomchord) for i in 1:n]...,
        [T(t[i], t[i+1], :topchord) for i in 1:n-2]..., T(b[1], t[1], :endpost), T(t[end], b[end], :endpost),
        [T(b[i+1], t[i], :vertical) for i in 1:n-1]...,
        [T(t[i], b[i+2], :diagonal) for i in 1:n÷2-1]..., [T(t[i+1], b[i+1], :diagonal) for i in n÷2:n-2]...]
    m = Model(nodes, els, AbstractLoad{Float64}[NodeForce(x, [0.0, -P, 0.0]) for x in b[2:end-1]])
    solve!(m)
    return m, b, t, els
end
