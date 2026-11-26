#=
Enzyme across Asap — one case × mode per process (run via run_enzyme_matrix.sh).

    julia --project=examples examples/ad/enzyme_matrix.jl <case> <mode>

case ∈ truss_area, truss_shape, spaceframe, frame_lineload, frame_joints,
       springs, variable_element, shells, shells_tri, shells_mixed, diaphragm, cases
mode ∈ reverse, forward

Reference: Zygote. Prints one line: `RESULT case mode OK rel_err time_ms first_s`.
Enzyme settings as everywhere in Asap: runtime activity + Const function
annotation (DifferentiationInterface).
=#
using Asap, AsapOptim, LinearAlgebra, Statistics
using Zygote, ChainRulesCore, Enzyme, DifferentiationInterface
import ForwardDiff

const CASE = ARGS[1]
const MODE = ARGS[2]

mat = Material(2.0e8, 7.85, 0.3)

# ── cases: each returns (f, x0) with f: Vector → scalar ─────────────────────
function truss2d()
    sec = Section(mat, 1e-3)
    free(x, y) = Node([x, y, 0.0], [true, true, false, false, false, false])
    nodes = [Node([0.0, 0.0, 0.0], [false, false, false, false, false, false]),
             free(2.0, 0.0), free(4.0, 0.0), free(6.0, 0.0),
             Node([8.0, 0.0, 0.0], [true, false, false, false, false, false]),
             free(2.0, 1.5), free(4.0, 1.5), free(6.0, 1.5)]
    bars = [(1, 2), (2, 3), (3, 4), (4, 5), (6, 7), (7, 8), (1, 6), (6, 2), (2, 7), (7, 3), (3, 8), (8, 4), (8, 5), (6, 3), (3, 8)]
    els = AbstractElement{Float64}[TrussElement(nodes[a], nodes[b], sec) for (a, b) in unique(bars)]
    loads = AbstractLoad{Float64}[NodeForce(nodes[k], [0.0, -10.0, 0.0]) for k in (2, 3, 4)]
    m = Model(nodes, els, loads)
    solve!(m)
    return m, nodes, els
end

function case(name)
    if name == "truss_area"
        m, nodes, els = truss2d()
        p = OptParams(m, AbstractVariable[AreaVariable(el, 1e-3, 1e-5, 1e-2) for el in els])
        return x -> compliance(solve_structure(x, p), p), copy(p.values)
    elseif name == "truss_shape"
        m, nodes, els = truss2d()
        p = OptParams(m, AbstractVariable[SpatialVariable(nodes[k], 0.0, -1.0, 1.0, :Y) for k in 6:8])
        return x -> compliance(solve_structure(x, p), p), copy(p.values) .+ 0.1
    elseif name == "spaceframe"
        sec = Section(mat, 1e-3)
        sf = SpaceFrame(8, 1.0, 8, 1.0, 1.0, sec; support = :corner)
        m = sf.model
        p = OptParams(m, AbstractVariable[AreaVariable(el, 1e-3, 1e-5, 1e-2) for el in m.elements])
        println("spaceframe variables: ", length(p.values))
        return x -> compliance(solve_structure(x, p), p), copy(p.values)
    elseif name == "frame_lineload"
        sec = Section(mat, 1e-2, 1e-4, 5e-5, 1e-6)
        n1 = Node([0.0, 0, 0], :fixed); n2 = Node([0.0, 0, 3], :free); n3 = Node([5.0, 0, 3], :free); n4 = Node([5.0, 0, 0], :fixed)
        beam = FrameElement(n2, n3, sec)
        m = Model([n1, n2, n3, n4], AbstractElement{Float64}[FrameElement(n1, n2, sec), beam, FrameElement(n4, n3, sec)],
            AbstractLoad{Float64}[LineLoad(beam, [0.0, 0.0, -10.0]), NodeForce(n2, [5.0, 0, 0]), SelfWeight(beam)])
        solve!(m)
        p = OptParams(m, AbstractVariable[SpatialVariable(n3, 0.0, -1.0, 1.0, :X), SpatialVariable(n2, 0.0, -1.0, 1.0, :Z),
            AreaVariable(beam, 1e-2, 1e-4, 1e-1)])
        return x -> compliance(solve_structure(x, p), p), copy(p.values) .+ [0.1, 0.1, 0.0]
    elseif name == "frame_joints"
        sec = Section(mat, 1e-2, 1e-4, 5e-5, 1e-6)
        n1 = Node([0.0, 0, 0], :pinned); n2 = Node([0.0, 0, 3], :free); n3 = Node([5.0, 0, 3], :free); n4 = Node([5.0, 0, 0], :pinned)
        beam = FrameElement(n2, n3, sec)
        m = Model([n1, n2, n3, n4], AbstractElement{Float64}[FrameElement(n1, n2, sec), beam, FrameElement(n4, n3, sec)],
            AbstractLoad{Float64}[NodeForce(n2, [5.0, 0, -2.0])])
        solve!(m)
        p = OptParams(m, AbstractVariable[JointVariable(beam, :both, 1e4, 1.0, 1e8)])
        return x -> compliance(solve_structure(x, p), p), copy(p.values)
    elseif name == "springs"
        sec = Section(mat, 1e-2, 1e-4, 5e-5, 1e-6)
        base = Node([0.0, 0, 0], :free); top = Node([0.0, 0, 3], :free); far = Node([4.0, 0, 3], :pinned)
        m = Model([base, top, far], AbstractElement{Float64}[FrameElement(base, top, sec), FrameElement(top, far, sec)],
            AbstractLoad{Float64}[NodeForce(top, [5.0, 0, -3.0])]; springs = [NodalSpring(base, [1e4, 1e4, 1e5, 1e3, 1e3, 1e3])])
        solve!(m)
        p = OptParams(m, AbstractVariable[SpatialVariable(top, 0.0, -1.0, 1.0, :X)])
        return x -> compliance(solve_structure(x, p), p), [0.2]
    elseif name == "variable_element"
        deep = Section(mat, 2e-2, 4e-4, 1e-4, 1e-5); sh = Section(mat, 1e-2, 1e-4, 5e-5, 5e-6)
        n1 = Node([0.0, 0, 0], :fixed); n2 = Node([8.0, 0, 0], :fixed); n3 = Node([4.0, 0, 3.0], :fixed)
        h = VariableElement(n1, n2, AbstractSection{Float64}[deep, sh, deep], [0.3, 0.7])
        m = Model([n1, n2, n3], AbstractElement{Float64}[h, TrussElement(n3, n1, Section(mat, 1e-3))],
            AbstractLoad{Float64}[PointLoad(h, 0.5, [0.0, 0.0, -10.0])])
        solve!(m); st = extract_state(m)
        return x -> compliance(m, ModelState{eltype(x)}(reshape(x, 3, :), st.sections)), vec(copy(st.X))
    elseif name == "shells"
        conc = Material(3.0e4, 2.4e-3, 0.2)
        w = mesh_wall([0.0, 0, 0], [1.0, 0, 0], [0.0, 0, 1], 3.0, 2.0, ShellSection(conc, 0.2); size = 0.5)
        top = [w.grid[i, end] for i in eachindex(w.s)]
        m = Model(w.nodes, AbstractElement{Float64}[w.elements...],
            AbstractLoad{Float64}[[SelfWeight(el; g = [0.0, 0, -9.81]) for el in w.elements]..., [NodeForce(n, [0.5, 0, 0]) for n in top]...])
        solve!(m)
        half = [el for el in w.elements if sum(n -> n.position[1], el.nodes) / 4 < 1.5]
        rest = [el for el in w.elements if !(el in half)]
        p = OptParams(m, AbstractVariable[ThicknessVariable(half, 0.1, 0.5), ThicknessVariable(rest, 0.1, 0.5)])
        return x -> compliance(solve_structure(x, p), p), [0.25, 0.18]
    elseif name == "shells_tri"
        # the same wall with every quad split into two triangles (ANDES-OPT + MITC3+)
        conc = Material(3.0e4, 2.4e-3, 0.2)
        sec = ShellSection(conc, 0.2)
        w = mesh_wall([0.0, 0, 0], [1.0, 0, 0], [0.0, 0, 1], 3.0, 2.0, sec; size = 0.5)
        tris = [t for el in w.elements for t in (ShellElement((el.nodes[1], el.nodes[2], el.nodes[3]), sec),
                                                 ShellElement((el.nodes[1], el.nodes[3], el.nodes[4]), sec))]
        top = [w.grid[i, end] for i in eachindex(w.s)]
        m = Model(w.nodes, AbstractElement{Float64}[tris...],
            AbstractLoad{Float64}[[SelfWeight(el; g = [0.0, 0, -9.81]) for el in tris]..., [NodeForce(n, [0.5, 0, 0]) for n in top]...])
        solve!(m)
        half = [el for el in tris if sum(n -> n.position[1], el.nodes) / 3 < 1.5]
        rest = [el for el in tris if !(el in half)]
        p = OptParams(m, AbstractVariable[ThicknessVariable(half, 0.1, 0.5), ThicknessVariable(rest, 0.1, 0.5)])
        return x -> compliance(solve_structure(x, p), p), [0.25, 0.18]
    elseif name == "shells_mixed"
        # the wall with its right half triangulated: the triangles beside the quads take the
        # quad's drilling lumping on the shared edges (Asap._shell3_K_plain and its rule)
        conc = Material(3.0e4, 2.4e-3, 0.2)
        sec = ShellSection(conc, 0.2)
        w = mesh_wall([0.0, 0, 0], [1.0, 0, 0], [0.0, 0, 1], 3.0, 2.0, sec; size = 0.5)
        els = AbstractElement{Float64}[]
        for el in w.elements
            sum(n -> n.position[1], el.nodes) / 4 > 1.5 ?
                append!(els, (ShellElement((el.nodes[1], el.nodes[2], el.nodes[3]), sec),
                              ShellElement((el.nodes[1], el.nodes[3], el.nodes[4]), sec))) : push!(els, el)
        end
        top = [w.grid[i, end] for i in eachindex(w.s)]
        m = Model(w.nodes, els,
            AbstractLoad{Float64}[[SelfWeight(el; g = [0.0, 0, -9.81]) for el in els]..., [NodeForce(n, [0.5, 0, 0]) for n in top]...])
        solve!(m)
        tris = [el for el in els if el isa ShellElement{Float64,3}]
        quads = [el for el in els if el isa ShellElement{Float64,4}]
        p = OptParams(m, AbstractVariable[ThicknessVariable(quads, 0.1, 0.5), ThicknessVariable(tris, 0.1, 0.5),
            SpatialVariable(w.grid[4, 2], 0.0, -0.2, 0.2, :X)])
        return x -> compliance(solve_structure(x, p), p), [0.25, 0.18, 0.05]
    elseif name == "diaphragm"
        base = [Node([x, y, 0.0], :fixed) for (x, y) in ((0.0, 0.0), (6.0, 0.0), (0.0, 4.0), (6.0, 4.0))]
        top = [Node([x, y, 3.0], :free) for (x, y) in ((0.0, 0.0), (6.0, 0.0), (0.0, 4.0), (6.0, 4.0))]
        cols = AbstractElement{Float64}[FrameElement(base[k], top[k], Section(mat, 1e-2, 1e-4k, 1e-4k, 5e-5)) for k in 1:4]
        d = RigidDiaphragm(top)
        m = Model(vcat(base, top), cols, AbstractLoad{Float64}[NodeForce(d.master, [10.0, -4.0, 0.0])]; diaphragms = [d])
        solve!(m)
        p = OptParams(m, AbstractVariable[SpatialVariable(top[2], 0.0, -0.5, 0.5, :X), SpatialVariable(top[3], 0.0, -0.5, 0.5, :Y)])
        return x -> compliance(solve_structure(x, p), p), [0.1, -0.05]
    elseif name == "cases"
        # pure-path load cases (CaseLoads): a slab of quads and triangles on an edge beam and a
        # column, three cases with design-dependent frame and shell loads; a weighted sum of the
        # per-case compliances, in the node positions
        conc = Material(3.0e4, 2.4e-3, 0.2)
        slab = ShellSection(conc, 0.2)
        P = [0.0 2.0 4.0 0.0 2.0 4.0 0.0 2.0 4.0; 0.0 0.0 0.0 2.0 2.0 2.0 4.0 4.0 4.0]
        sm = shell_mesh(P, [(1, 2, 5, 4), (2, 3, 6, 5), (4, 5, 8), (4, 8, 7), (5, 6, 9), (5, 9, 8)], slab)
        n = sm.nodes
        foreach(k -> fixnode!(n[k], Bool[0, 0, 0, 0, 0, 0]), (1, 3))
        base = Node([4.0, 4.0, -3.0], :fixed)
        sec = Section(mat, 1e-2, 2e-4, 1e-4, 5e-5)
        edge = FrameElement(n[7], n[8], sec); beam = FrameElement(n[8], n[9], sec); col = FrameElement(base, n[9], sec)
        quads = filter(el -> el isa ShellElement{Float64,4}, sm.elements)
        tris = filter(el -> el isa ShellElement{Float64,3}, sm.elements)
        loads = AbstractLoad{Float64}[[SelfWeight(el; g = [0.0, 0, -9.81], case = :dead) for el in sm.elements]...,
            LineLoad(edge, [0.0, 0.0, -5.0]; case = :dead), LineLoad(beam, [0.0, 0.0, -3.0]; case = :dead),
            AreaLoad(quads[1], -2.0; case = :live), AreaLoad(tris[2], -3.0; case = :live),
            PointLoad(edge, 0.3, [0.0, 0.0, -10.0]; case = :live),
            NodeForce(n[9], [6.0, 0.0, 0.0]; case = :wind), LineLoad(col, [1.5, 0.0, 0.0]; case = :wind)]
        m = Model([n; base], AbstractElement{Float64}[sm.elements..., edge, beam, col], loads)
        solve!(m)
        cl = case_loads(m); st = extract_state(m)
        λ = [1.2, 1.6, 0.5]
        x0 = vec(copy(st.X)) .+ 0.02 .* sin.(1:length(st.X))
        return x -> dot(λ, compliance(m, ModelState{eltype(x)}(reshape(x, 3, :), st.sections), cl)), x0
    end
    error("unknown case $name")
end

f, x0 = case(CASE)
# reference gradient: Zygote; where Zygote cannot differentiate the model,
# central finite differences (and say so)
ref = "zygote"
gz = try
    Zygote.gradient(f, x0)[1]
catch err
    global ref = "fd (zygote: " * first(split(sprint(showerror, err), '\n'))[1:min(end, 60)] * ")"
    h = 1e-6 .* max.(abs.(x0), 1.0)
    [(f(x0 .+ h[j] .* (eachindex(x0) .== j)) - f(x0 .- h[j] .* (eachindex(x0) .== j))) / 2h[j] for j in eachindex(x0)]
end
mode = MODE == "reverse" ? Enzyme.set_runtime_activity(Enzyme.Reverse) : Enzyme.set_runtime_activity(Enzyme.Forward)
backend = AutoEnzyme(; mode = mode, function_annotation = Enzyme.Const)
t0 = time()
ge = DifferentiationInterface.gradient(f, backend, x0)
first_s = time() - t0
tms = 1e3 * median([@elapsed DifferentiationInterface.gradient(f, backend, x0) for _ in 1:3])
tz = ref == "zygote" ? 1e3 * median([@elapsed Zygote.gradient(f, x0) for _ in 1:3]) : NaN
println("RESULT $CASE $MODE OK rel_err=$(norm(ge - gz) / norm(gz)) enzyme_ms=$(round(tms, digits = 3)) zygote_ms=$(round(tz, digits = 3)) first_s=$(round(first_s, digits = 1)) n=$(length(x0)) ref=$ref")
