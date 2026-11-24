#=
REPL display methods for Asap types.

Provides human-readable output when a variable is entered in the REPL.
Each method shows the type name, then aligned `symbol = value (plain-language
explanation)` lines, with dimensions in bracket notation — Asap is
units-agnostic, so dimensions are shown rather than named units.

Follows the FormaSlab.jl ShowMethods pattern. One entry per user-facing type,
grouped by category. Derived quantities that an engineer would reach for
(e.g. a Section's rigidities) are computed and shown inline.
=#

# ─────────────────────────────────────────────────────────────────────────────
# Materials
# ─────────────────────────────────────────────────────────────────────────────

function Base.show(io::IO, ::MIME"text/plain", m::Material)
    println(io, "Material")
    println(io, "  E = $(m.E)  [force/length²]  (Young's modulus — normal stiffness)")
    println(io, "  G = $(m.G)  [force/length²]  (shear modulus)")
    println(io, "  ρ = $(m.ρ)  [mass/length³]   (mass density)")
    print(io,   "  ν = $(m.ν)                   (Poisson's ratio — transverse strain ratio)")
end

Base.show(io::IO, m::Material) = print(io, "Material(E=$(m.E), G=$(m.G), ρ=$(m.ρ), ν=$(m.ν))")

# ─────────────────────────────────────────────────────────────────────────────
# Sections
# ─────────────────────────────────────────────────────────────────────────────

function Base.show(io::IO, ::MIME"text/plain", s::Section)
    println(io, "Section")
    println(io, "  A  = $(s.A)  [length²]  (cross-sectional area)")
    println(io, "  Ix = $(s.Ix)  [length⁴]  (second moment of area, strong axis)")
    println(io, "  Iy = $(s.Iy)  [length⁴]  (second moment of area, weak axis)")
    println(io, "  J  = $(s.J)  [length⁴]  (torsional constant)")
    println(io, "  material: E = $(s.material.E), G = $(s.material.G), ρ = $(s.material.ρ), ν = $(s.material.ν)")
    println(io, "  rigidities (what the stiffness matrix consumes):")
    println(io, "    EA  = $(EA(s))  [force]          (axial)")
    println(io, "    EIx = $(EIx(s))  [force·length²]  (flexural, strong axis)")
    println(io, "    EIy = $(EIy(s))  [force·length²]  (flexural, weak axis)")
    println(io, "    GJ  = $(GJ(s))  [force·length²]  (torsional)")
    print(io,   "    ρA  = $(ρA(s))  [mass/length]    (mass per unit length)")
end

Base.show(io::IO, s::Section) = print(io, "Section(A=$(s.A), Ix=$(s.Ix), Iy=$(s.Iy), J=$(s.J))")

function Base.show(io::IO, ::MIME"text/plain", s::RigiditySection)
    println(io, "RigiditySection  (effective rigidities stored directly)")
    println(io, "  EA  = $(s.EA)  [force]          (axial rigidity)")
    println(io, "  EIx = $(s.EIx)  [force·length²]  (flexural rigidity, strong axis)")
    println(io, "  EIy = $(s.EIy)  [force·length²]  (flexural rigidity, weak axis)")
    println(io, "  GJ  = $(s.GJ)  [force·length²]  (torsional rigidity)")
    print(io,   "  ρA  = $(s.ρA)  [mass/length]    (mass per unit length)")
end

Base.show(io::IO, s::RigiditySection) =
    print(io, "RigiditySection(EA=$(s.EA), EIx=$(s.EIx), EIy=$(s.EIy), GJ=$(s.GJ))")

# ─────────────────────────────────────────────────────────────────────────────
# Nodes
# ─────────────────────────────────────────────────────────────────────────────

# human-readable DOF summary, e.g. "fixed: Tx Ty Tz" or "free"
function _fixity_description(fixity)
    names = ("Tx", "Ty", "Tz", "Rx", "Ry", "Rz")
    fixed = [names[i] for i in 1:6 if !fixity[i]]
    isempty(fixed) && return "free (no supports)"
    length(fixed) == 6 && return "fully fixed (clamped)"
    return "fixed: " * join(fixed, " ") * "  (all others free)"
end

function Base.show(io::IO, ::MIME"text/plain", n::Node)
    println(io, "Node :$(n.id)")
    println(io, "  position = $(collect(n.position))  [length]  (global X, Y, Z)")
    print(io, "  supports = $(_fixity_description(n.fixity))")
end

Base.show(io::IO, n::Node) = print(io, "Node(:$(n.id), $(collect(n.position)))")

# ─────────────────────────────────────────────────────────────────────────────
# End conditions (releases / semi-rigid connections)
# ─────────────────────────────────────────────────────────────────────────────

_spring_word(k) = isinf(k) ? "rigid" : iszero(k) ? "released" : string(k)

function Base.show(io::IO, ::MIME"text/plain", e::EndSprings)
    println(io, "EndSprings  (element-end connection stiffnesses, local axes)")
    println(io, "  kx = $(_spring_word(e.kx))  [force/length]      (axial)")
    println(io, "  kt = $(_spring_word(e.kt))  [force·length/rad]  (torsion — twist about element axis)")
    println(io, "  ky = $(_spring_word(e.ky))  [force·length/rad]  (bending rotation about local y)")
    print(io,   "  kz = $(_spring_word(e.kz))  [force·length/rad]  (bending rotation about local z)")
end

Base.show(io::IO, e::EndSprings) =
    print(io, "EndSprings(kx=$(_spring_word(e.kx)), kt=$(_spring_word(e.kt)), ky=$(_spring_word(e.ky)), kz=$(_spring_word(e.kz)))")

function Base.show(io::IO, ::MIME"text/plain", ec::EndConditions)
    sym = release_symbol(ec)
    header = sym === nothing ? "EndConditions  (semi-rigid)" : "EndConditions  (:$sym)"
    println(io, header)
    println(io, "  start end: $(sprint(show, ec.e1))")
    print(io,   "  far end:   $(sprint(show, ec.e2))")
end

Base.show(io::IO, ec::EndConditions) = begin
    sym = release_symbol(ec)
    sym === nothing ? print(io, "EndConditions($(ec.e1), $(ec.e2))") : print(io, "EndConditions(:$sym)")
end

# ─────────────────────────────────────────────────────────────────────────────
# Elements
# ─────────────────────────────────────────────────────────────────────────────

# shell sections

_fmt3(M) = join((join((string(round(M[i, j]; sigdigits = 5)) for j in axes(M, 2)), "  ") for i in axes(M, 1)), "; ")

function Base.show(io::IO, ::MIME"text/plain", s::ShellSection)
    println(io, "ShellSection  (homogeneous isotropic shell)")
    println(io, "  t = $(s.t)  [length]  (thickness)")
    println(io, "  material: E = $(s.material.E), G = $(s.material.G), ρ = $(s.material.ρ), ν = $(s.material.ν)")
    println(io, "  rigidities (what the stiffness matrix consumes):")
    println(io, "    A  = [$(_fmt3(membrane_rigidity(s)))]  [force/length]   (membrane, N = A·ε)")
    println(io, "    D  = [$(_fmt3(bending_rigidity(s)))]  [force·length]   (bending, M = D·κ)")
    println(io, "    S  = [$(_fmt3(shear_rigidity(s)))]  [force/length]   (transverse shear, Q = S·γ)")
    println(io, "    γt = $(drilling_rigidity(s))  [force/length]   (drilling penalty)")
    print(io,   "    ρt = $(areal_mass(s))  [mass/length²]   (mass per unit area)")
end

Base.show(io::IO, s::ShellSection) = print(io, "ShellSection(t=$(s.t))")

function Base.show(io::IO, ::MIME"text/plain", s::ShellRigiditySection)
    println(io, "ShellRigiditySection  (effective rigidities stored directly)")
    println(io, "  A  = [$(_fmt3(s.A))]  [force/length]   (membrane, N = A·ε)")
    println(io, "  D  = [$(_fmt3(s.D))]  [force·length]   (bending, M = D·κ)")
    println(io, "  S  = [$(_fmt3(s.S))]  [force/length]   (transverse shear, Q = S·γ)")
    println(io, "  γt = $(s.γt)  [force/length]   (drilling penalty)")
    print(io,   "  ρt = $(s.ρt)  [mass/length²]   (mass per unit area)")
end

Base.show(io::IO, s::ShellRigiditySection) = print(io, "ShellRigiditySection(A₁₁=$(s.A[1,1]), D₁₁=$(s.D[1,1]))")

function Base.show(io::IO, ::MIME"text/plain", el::ShellElement{T,N}) where {T,N}
    println(io, "ShellElement{$N} :$(el.id)  (flat shell: drilling membrane + MITC plate)")
    for (k, n) in enumerate(el.nodes)
        println(io, "  node $k: $(collect(n.position))")
    end
    R = local_frame(el)
    println(io, "  area   = $(area(el))  [length²]")
    println(io, "  normal = $(collect(R[3, :]))  (local e₃; right-hand rule over the node order)")
    print(io,   "  section: ", el.section)
end

Base.show(io::IO, el::ShellElement{T,N}) where {T,N} = print(io, "ShellElement{$N}(:$(el.id))")

function Base.show(io::IO, ::MIME"text/plain", el::LinkElement)
    println(io, "LinkElement :$(el.id)  (zero-length 6-spring connection)")
    println(io, "  at $(collect(el.nodeStart.position))")
    println(io, "  k₁, k₂, k₃    = $(collect(el.k[SOneTo(3)]))  [force/length]       (along link axes)")
    print(io,   "  kθ₁, kθ₂, kθ₃ = $(collect(el.k[SVector(4, 5, 6)]))  [force·length/rad]  (about link axes)")
end

Base.show(io::IO, el::LinkElement) = print(io, "LinkElement(:$(el.id))")

function Base.show(io::IO, ::MIME"text/plain", m::WallMesh)
    println(io, "WallMesh  ($(length(m.elements)) shell elements, $(length(m.nodes)) nodes)")
    println(io, "  origin = $(collect(m.origin))  [length]")
    println(io, "  along  = $(collect(m.along))   (wall length direction, local e₁)")
    println(io, "  up     = $(collect(m.up))   (wall height direction)")
    println(io, "  length = $(m.s[end]), height = $(m.z[end])  [length]")
    print(io,   "  grid   = $(length(m.s)) × $(length(m.z)) lines (openings removed)")
end

Base.show(io::IO, m::WallMesh) = print(io, "WallMesh($(length(m.elements)) elements)")

function Base.show(io::IO, ::MIME"text/plain", el::FrameElement)
    sym = release_symbol(el.ends)
    conn = sym === nothing ? "semi-rigid connections" : "connections :$sym"
    println(io, "FrameElement :$(el.id)")
    println(io, "  from $(collect(el.nodeStart.position)) to $(collect(el.nodeEnd.position))")
    println(io, "  length = $(length(el))  [length]")
    println(io, "  rollangle = $(el.rollangle)  [rad]  (section roll about the element axis)")
    println(io, "  $conn")
    println(io, "  section rigidities:")
    println(io, "    EA  = $(EA(el.section)), EIx = $(EIx(el.section)),")
    print(io,   "    EIy = $(EIy(el.section)), GJ = $(GJ(el.section))")
end

Base.show(io::IO, el::FrameElement) =
    print(io, "FrameElement(:$(el.id), L=$(length(el)))")

function Base.show(io::IO, ::MIME"text/plain", el::TrussElement)
    println(io, "TrussElement :$(el.id)  (axial-only two-force member)")
    println(io, "  from $(collect(el.nodeStart.position)) to $(collect(el.nodeEnd.position))")
    println(io, "  length = $(length(el))  [length]")
    print(io,   "  EA = $(EA(el.section))  [force]  (axial rigidity)")
end

Base.show(io::IO, el::TrussElement) =
    print(io, "TrussElement(:$(el.id), L=$(length(el)))")

function Base.show(io::IO, ::MIME"text/plain", el::VariableElement)
    nseg = n_segments(el)
    sym = release_symbol(el.ends)
    conn = sym === nothing ? "semi-rigid outer connections" : "outer connections :$sym"
    println(io, "VariableElement :$(el.id)  ($nseg prismatic segments, one member)")
    println(io, "  from $(collect(el.nodeStart.position)) to $(collect(el.nodeEnd.position))")
    println(io, "  length = $(length(el))  [length];  rollangle = $(el.rollangle)  [rad]")
    println(io, "  $conn; interior joints rigid ($(n_internal_dofs(el)) internal DOFs)")
    ts = segment_fractions(el)
    for s in 1:nseg
        sec = el.sections[s]
        tail = s == nseg ? "" : "\n"
        print(io, "  [$s] t ∈ [$(round(ts[s]; digits=3)), $(round(ts[s+1]; digits=3))]:  " *
                  "EA = $(EA(sec)), EIx = $(EIx(sec)), EIy = $(EIy(sec)), GJ = $(GJ(sec))$tail")
    end
end

Base.show(io::IO, el::VariableElement) =
    print(io, "VariableElement(:$(el.id), $(n_segments(el)) segments, L=$(length(el)))")

# ─────────────────────────────────────────────────────────────────────────────
# Springs
# ─────────────────────────────────────────────────────────────────────────────

function Base.show(io::IO, ::MIME"text/plain", sp::NodalSpring)
    dofnames = ("Tx", "Ty", "Tz", "Rx", "Ry", "Rz")
    active = [i for i in 1:6 if sp.stiffness[i] > 0]
    println(io, "NodalSpring :$(sp.id)  (elastic support at node :$(sp.node.id))")
    if isempty(active)
        print(io, "  no nonzero stiffnesses")
    else
        for (n, i) in enumerate(active)
            unit = i <= 3 ? "[force/length]" : "[force·length/rad]"
            tail = n == length(active) ? "" : "\n"
            print(io, "  k($(dofnames[i])) = $(sp.stiffness[i])  $unit$tail")
        end
    end
end

Base.show(io::IO, sp::NodalSpring) =
    print(io, "NodalSpring(:$(sp.id) at :$(sp.node.id))")

# ─────────────────────────────────────────────────────────────────────────────
# Model & results
# ─────────────────────────────────────────────────────────────────────────────

function Base.show(io::IO, ::MIME"text/plain", d::RigidDiaphragm)
    ax = ("X", "Y", "Z")[d.axis]
    println(io, "RigidDiaphragm :$(d.id)  (rigid in-plane floor, normal $ax)")
    println(io, "  slaves = $(length(d.nodes)) nodes  (in-plane DOFs follow the master rigidly)")
    print(io,   "  master at $(collect(d.master.position))  (apply storey forces/torques here)")
end

Base.show(io::IO, d::RigidDiaphragm) = print(io, "RigidDiaphragm(:$(d.id), $(length(d.nodes)) nodes)")

function Base.show(io::IO, ::MIME"text/plain", m::Model)
    nf = count(el -> el isa FrameElement, m.elements)
    nt = count(el -> el isa TrussElement, m.elements)
    nsh = count(el -> el isa ShellElement, m.elements)
    nl = count(el -> el isa LinkElement, m.elements)
    nother = length(m.elements) - nf - nt - nsh - nl
    kinds = "$nf frame, $nt truss" * (nsh > 0 ? ", $nsh shell" : "") * (nl > 0 ? ", $nl link" : "") *
            (nother > 0 ? ", $nother other" : "")
    println(io, "Model")
    println(io, "  $(length(m.nodes)) nodes, $(length(m.elements)) elements ($kinds)")
    println(io, "  $(length(m.loads)) loads, $(length(m.springs)) nodal springs" *
                (isempty(m.diaphragms) ? "" : ", $(length(m.diaphragms)) rigid diaphragm(s)"))
    if m.cache !== nothing
        p = m.cache.partition
        println(io, "  processed: $(length(p.free)) free / $(length(p.fixed)) fixed / " *
                    "$(length(p.inactive)) inactive DOFs")
    else
        println(io, "  unprocessed (call process! or solve!)")
    end
    print(io, m.results === nothing ? "  unsolved" : "  solved (see model.results)")
end

Base.show(io::IO, m::Model) =
    print(io, "Model($(length(m.nodes)) nodes, $(length(m.elements)) elements)")

function Base.show(io::IO, ::MIME"text/plain", r::LinearResults)
    umax = isempty(r.u) ? 0.0 : maximum(abs, r.u)
    println(io, "LinearResults")
    println(io, "  max |u| = $umax  [length or rad]  (largest displacement component)")
    println(io, "  compliance = $(r.compliance)  [force·length]  (external work uᵀF)")
    print(io,   "  access: displacement(res, node), reaction(res, node), element_forces(res, el)")
end

Base.show(io::IO, r::LinearResults) =
    print(io, "LinearResults(compliance=$(r.compliance))")

# ─────────────────────────────────────────────────────────────────────────────
# Parametric structure generators
# ─────────────────────────────────────────────────────────────────────────────

#=
Generators share one field-introspecting show: the solved model/network is
summarized, scalar generation parameters are listed by name (they mirror the
constructor arguments — see each generator's docstring for meanings), and
bulky index bookkeeping arrays are summarized by size only.
=#
function _show_generator_fields(io::IO, g)
    for f in fieldnames(typeof(g))
        v = getfield(g, f)
        pad = rpad(String(f), 14)
        if v isa Model
            println(io, "  $pad = Model: $(length(v.nodes)) nodes, $(length(v.elements)) elements" *
                        (v.results === nothing ? "" : " (solved)"))
        elseif v isa Network
            println(io, "  $pad = Network: $(length(v.nodes)) nodes, $(length(v.elements)) elements")
        elseif v isa AbstractSection
            println(io, "  $pad = $(nameof(typeof(v))) (EA = $(EA(v)))")
        elseif v isa Real || v isa Symbol
            println(io, "  $pad = $v")
        elseif v isa AbstractVector{<:Real} && length(v) <= 3
            println(io, "  $pad = $v")
        elseif v isa AbstractArray
            println(io, "  $pad = $(join(size(v), "×")) array (index bookkeeping)")
        end
    end
end

function Base.show(io::IO, ::MIME"text/plain", g::AbstractGenerator)
    println(io, nameof(typeof(g)), "  (parametric structure generator)")
    _show_generator_fields(io, g)
    print(io, "  access the solved structure via the model/network field")
end

Base.show(io::IO, g::AbstractGenerator) =
    print(io, nameof(typeof(g)), "(", join(["$f=$(getfield(g, f))"
        for f in fieldnames(typeof(g)) if getfield(g, f) isa Union{Real,Symbol}], ", "), ")")

function Base.show(io::IO, ::MIME"text/plain", g::GroundStructure)
    println(io, nameof(typeof(g)), "  (ground structure — candidate member grid for layout optimization)")
    _show_generator_fields(io, g)
    print(io, "  materialize with to_truss(gs, section) or to_frame(gs, section)")
end

Base.show(io::IO, g::GroundStructure) =
    print(io, nameof(typeof(g)), "(", join(["$f=$(getfield(g, f))"
        for f in fieldnames(typeof(g)) if getfield(g, f) isa Union{Real,Symbol}], ", "), ")")

# ─────────────────────────────────────────────────────────────────────────────
# Plot-ready geometry extraction
# ─────────────────────────────────────────────────────────────────────────────

function Base.show(io::IO, ::MIME"text/plain", g::ModelGeo)
    println(io, "ModelGeo  (plot-ready arrays from a solved frame model)")
    println(io, "  $(length(g.nodes)) nodes, $(length(g.indices)) elements")
    println(io, "  max |u|  = $(maximum(norm.(g.disp)))  [length]         (largest nodal displacement)")
    println(io, "  max |N|  = $(g.max_abs_P)  [force]          (largest axial end force)")
    println(io, "  max |M|  = $(max(g.max_abs_My, g.max_abs_Mz))  [force·length]   (largest bending end moment)")
    print(io,   "  fields: nodes, disp, indices, P, Vy, Vz, Tx, My, Mz, areas, lengths, …")
end

Base.show(io::IO, g::ModelGeo) =
    print(io, "ModelGeo($(length(g.nodes)) nodes, $(length(g.indices)) elements)")

function Base.show(io::IO, ::MIME"text/plain", g::TrussGeo)
    println(io, "TrussGeo  (plot-ready arrays from a solved truss model)")
    println(io, "  $(length(g.nodes)) nodes, $(length(g.indices)) elements")
    println(io, "  max |u| = $(maximum(norm.(g.disp)))  [length]  (largest nodal displacement)")
    println(io, "  max |N| = $(g.max_abs_force)  [force]   (largest axial force)")
    print(io,   "  fields: nodes, disp, indices, forces, areas, lengths, …")
end

Base.show(io::IO, g::TrussGeo) =
    print(io, "TrussGeo($(length(g.nodes)) nodes, $(length(g.indices)) elements)")

function Base.show(io::IO, ::MIME"text/plain", g::NetworkGeo)
    println(io, "NetworkGeo  (plot-ready arrays from a solved FDM network)")
    println(io, "  $(length(g.nodes)) nodes, $(length(g.indices)) elements")
    print(io,   "  fields: nodes, indices, forces, lengths, …")
end

Base.show(io::IO, g::NetworkGeo) =
    print(io, "NetworkGeo($(length(g.nodes)) nodes, $(length(g.indices)) elements)")

function Base.show(io::IO, ::MIME"text/plain", d::ElementDisplacements)
    println(io, "ElementDisplacements  (displaced shape sampled along one member)")
    println(io, "  element    = :$(d.element.id) ($(nameof(typeof(d.element))))")
    println(io, "  resolution = $(d.resolution) stations")
    println(io, "  max |u|    = $(maximum(norm.(eachcol(d.uglobal))))  [length]  (largest sampled displacement)")
    print(io,   "  fields: x [station], ulocal [3×n, LCS], uglobal [3×n, GCS], basepositions [3×n]")
end

Base.show(io::IO, d::ElementDisplacements) =
    print(io, "ElementDisplacements(:$(d.element.id), $(d.resolution) stations)")

# ─────────────────────────────────────────────────────────────────────────────
# FDM (force density method)
# ─────────────────────────────────────────────────────────────────────────────

function Base.show(io::IO, ::MIME"text/plain", n::FDMnode)
    axes = join((n.fixity[a] ? "free" : "FIXED" for a in 1:3), "/")
    println(io, "FDMnode")
    println(io, "  position = $(n.position)  [length]")
    println(io, "  fixity   = x/y/z: $axes  (true = free)")
    print(io,   "  id       = :$(n.id)")
end

Base.show(io::IO, n::FDMnode) =
    print(io, "FDMnode(:$(n.id), $(n.position), $(collect(n.fixity)))")

function Base.show(io::IO, ::MIME"text/plain", e::FDMelement)
    println(io, "FDMelement")
    println(io, "  q  = $(e.q)  [force/length]  (force density)")
    print(io,   "  id = :$(e.id)")
end

Base.show(io::IO, e::FDMelement) = print(io, "FDMelement(:$(e.id), q=$(e.q))")

function Base.show(io::IO, ::MIME"text/plain", l::FDMload)
    println(io, "FDMload")
    print(io, "  force = $(l.force)  [force] at node :$(l.point.id)")
end

Base.show(io::IO, l::FDMload) = print(io, "FDMload($(l.force))")

function Base.show(io::IO, ::MIME"text/plain", nw::Network)
    println(io, "Network  (force density method)")
    println(io, "  $(length(nw.nodes)) nodes, $(length(nw.elements)) elements, $(length(nw.loads)) loads")
    if nw.cache !== nothing
        println(io, "  $(length(nw.cache.N)) fully free / $(length(nw.cache.F)) supported nodes" *
                    (nw.cache.mixed ? "  (per-axis fixity: separable solves)" : ""))
    else
        println(io, "  unprocessed (call process! or solve!)")
    end
    print(io, "  solve!(network) finds the equilibrium geometry in place")
end

Base.show(io::IO, nw::Network) =
    print(io, "Network($(length(nw.nodes)) nodes, $(length(nw.elements)) elements)")
