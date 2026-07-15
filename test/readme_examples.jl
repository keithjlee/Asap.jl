# Every code block from README.md, run in sequence — no untested examples.
using Pkg
Pkg.activate(; temp=true, io=devnull)
Pkg.develop([PackageSpec(path="/Users/keithlee/Documents/dev/Asap")]; io=devnull)
Pkg.add("Zygote"; io=devnull)
using Asap, LinearAlgebra

# quick start ---------------------------------------------------------------
steel = Material(200e6, 80.0, 0.3)
column = Section(steel, 1e-2, 8e-5, 3e-5, 5e-7)
n1 = Node([0.0, 0.0, 0.0], :fixed)
n2 = Node([0.0, 0.0, 3.0], :free)
el = FrameElement(n1, n2, column)
model = Model([n1, n2], [el], [NodeForce(n2, [10.0, 0.0, 0.0])])
solve!(model)
@assert displacement(model.results, n2) isa Any
@assert isfinite(moment_z(model, el, 0.0))
println("quickstart OK")

# materials & sections ------------------------------------------------------
steel = Material(200e6, 77e6, 8.0, 0.3)
timber = Material(13.1e6, 560.0, 0.35)
wshape = Section(steel, 1e-2, 8e-5, 3e-5, 5e-7)
bar = Section(steel, 5e-3)
EA(wshape), EIx(wshape)
Ec, Ig, Ag = 30e6, 8e-4, 0.12
cracked_beam = RigiditySection(Ec*Ag, 0.35*Ec*Ig, 0.35*Ec*Ig, 0.10*Ec*Ig, 2.4*Ag)
println("sections OK")

# nodes ----------------------------------------------------------------------
base = Node([0.0, 0.0, 0.0], :fixed)
pin = Node([5.0, 0.0, 0.0], :pinned)
roller = Node([10.0, 0.0, 0.0], [true, true, false, true, true, true])
fixnode!(roller, :zfixed)
println("nodes OK")

# frame elements + releases ---------------------------------------------------
na = Node([0.0, 0.0, 0.0], :fixed); nb = Node([4.0, 0.0, 0.0], :free)
girder = FrameElement(na, nb, wshape, :girder)
rolled = FrameElement(na, nb, wshape; rollangle = 0.0)
hinged = FrameElement(na, nb, wshape; release = :fixedfree)
semi = FrameElement(na, nb, wshape,
    EndConditions(EndSprings(Inf, Inf, 5e4, 5e4), rigid_end()), :connection)
println("frame elements OK")

# mixed truss + frame ----------------------------------------------------------
m1 = Node([0.0, 0.0, 0.0], :fixed)
m2 = Node([3.0, 0.0, 0.0], :free)
m3 = Node([3.0, 0.0, 3.0], :fixed)
beam = FrameElement(m1, m2, wshape, :beam)
tie = TrussElement(m2, m3, bar, :tie)
mixed = Model([m1, m2, m3],
    AbstractElement{Float64}[beam, tie],
    AbstractLoad{Float64}[NodeForce(m2, [0.0, 0.0, -100.0])])
solve!(mixed)
@assert axial_force(mixed.results, tie) > 0
println("mixed OK")

# variable element --------------------------------------------------------------
deep = Section(steel, 2e-2, 4e-4, 1e-4, 2e-6)
mid = Section(steel, 1.5e-2, 2e-4, 6e-5, 1e-6)
shallow = wshape
v1 = Node([0.0, 0.0, 0.0], :fixed); v2 = Node([8.0, 0.0, 0.0], :free)
haunched = VariableElement(v1, v2,
    AbstractSection{Float64}[deep, mid, shallow], [0.25, 0.6])
vmodel = Model([v1, v2], AbstractElement{Float64}[haunched],
    AbstractLoad{Float64}[LineLoad(haunched, [0.0, 0.0, -2.0])])
solve!(vmodel)
@assert isfinite(moment_z(vmodel, haunched, 0.5))
@assert length(element_forces(vmodel.results, haunched, 2)) == 12
push!(vmodel.loads, SelfWeight(haunched; g = [0.0, 0.0, -9.81]))
solve!(vmodel)
println("variable element OK")

# springs -------------------------------------------------------------------------
sb = Node([0.0, 0.0, 0.0], :fixed)
st_ = Node([0.0, 0.0, 3.0], :free)
soil = NodalSpring(st_, [0.0, 0.0, 5e4, 0.0, 0.0, 0.0], :soil)
pad = NodalSpring(st_, 1e5)
smodel = Model([sb, st_], AbstractElement{Float64}[FrameElement(sb, st_, wshape)],
    AbstractLoad{Float64}[NodeForce(st_, [1.0, 0.0, -10.0])]; springs = [soil])
solve!(smodel)
println("springs OK")

# loads ---------------------------------------------------------------------------
NodeForce(m2, [0.0, 0.0, -50.0]; case = :live)
NodeMoment(m2, [0.0, 1e3, 0.0])
LineLoad(beam, [0.0, 0.0, -2.0])
TrapezoidLoad(beam, 0.2, 0.8, 1.0, 3.0, [0.0, 0.0, -1.0])
DistributedLoad(beam, [0.0, 0.3, 0.5], [0.0, 4.0, 0.0], [0.0, 0.0, -1.0])
PointLoad(beam, 0.4, [0.0, 0.0, -10.0])
PointMoment(beam, 0.5, [0.0, 0.0, 800.0])
SelfWeight(beam; g = [0.0, 0.0, -9.81])
println("loads OK")

# recovery ---------------------------------------------------------------------------
rb1 = Node([0.0, 0.0, 0.0], :fixed); rb2 = Node([6.0, 0.0, 0.0], :free)
rbeam = FrameElement(rb1, rb2, wshape, :beam)
rmodel = Model([rb1, rb2], AbstractElement{Float64}[rbeam],
    AbstractLoad{Float64}[LineLoad(rbeam, [0.0, 0.0, -2.0]),
        PointLoad(rbeam, 0.4, [0.0, 0.0, -10.0])])
solve!(rmodel)
moment_z(rmodel, rbeam, 0.5); shear_y(rmodel, rbeam, 0.25)
axial_force(rmodel, rbeam, 0.0); torsion(rmodel, rbeam, 0.5)
f = InternalForces(rmodel, rbeam; resolution = 40)
f.x; f.N; f.Vy; f.Mz; f.Vz; f.My; f.Mx
stt = internal_forces(rmodel, rbeam)
local_displacements(stt, 0.5)
println("recovery OK")

# cases/combos/envelopes ----------------------------------------------------------------
c1 = Node([0.0, 0.0, 0.0], :fixed); c2 = Node([0.0, 0.0, 3.0], :free)
c3 = Node([5.0, 0.0, 3.0], :free); c4 = Node([5.0, 0.0, 0.0], :fixed)
cbeam = FrameElement(c2, c3, wshape, :beam)
cels = AbstractElement{Float64}[FrameElement(c1, c2, wshape), cbeam, FrameElement(c4, c3, wshape)]
cloads = AbstractLoad{Float64}[
    LineLoad(cbeam, [0.0, 0.0, -2.0]; case = :dead),
    PointLoad(cbeam, 0.4, [0.0, 0.0, -10.0]; case = :live),
    NodeForce(c2, [5.0, 0.0, 0.0]; case = :wind)]
cmodel = Model([c1, c2, c3, c4], cels, cloads)
cr = solve_cases!(cmodel)
strength = LoadCombination(:LRFD, [:dead => 1.2, :live => 1.6, :wind => 0.5])
resc = combine(cr, strength)
displacement(resc, c2)
env = envelope(cmodel, cbeam, cr,
    [strength, LoadCombination(:service, [:dead => 1.0, :live => 1.0])])
env.x; env.lo; env.hi
println("cases OK")

# differentiable ----------------------------------------------------------------------
using Zygote
solve!(cmodel)
state = extract_state(cmodel)
g = Zygote.gradient(state.X) do X
    compliance(cmodel, ModelState{Float64}(X, state.sections))
end[1]
@assert size(g) == size(state.X) && any(!iszero, g)
println("differentiable OK")

# generators + geometry extraction -----------------------------------------------------
sec = Section(Steel_kNm, 1e-3, 1e-6, 1e-6, 1e-6)

truss = Warren2D(11, 1.5, 2.0, sec; load = [0.0, -20.0, 0.0])
grid = SpaceFrame(6, 1.2, 6, 1.2, 1.0, sec; support = :corner)
vault = SpaceFrame(6, 1.2, 6, 1.2, 1.0, (u, v) -> 0.5 * sinpi(u) * sinpi(v), sec)
@assert maximum(abs, truss.model.results.u) > 0

fsec = Section(Steel_kNm, 1e-2, 1e-4, 5e-5, 1e-6)
bldg = Frame(2, 6.0, 2, 5.0, 2, 4.0, 2.0, fsec, fsec, fsec, fsec)

gs = XGroundStructure(6.0, 4, 4.0, 3)
candidates = to_truss(gs, sec; load = [1.0, 0.0, 0.0])
@assert candidates isa Model

geo = Geo(truss.model)
@assert length(geo.nodes) == length(truss.model.nodes)

ed = ElementDisplacements(bldg.model.elements[1], bldg.model; resolution = 20)
@assert size(ed.basepositions .+ 100 .* ed.uglobal) == (3, 20)
println("generators OK")

println("ALL README EXAMPLES PASS")
