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

function Base.show(io::IO, ::MIME"text/plain", m::Model)
    nf = count(el -> el isa FrameElement, m.elements)
    nt = count(el -> el isa TrussElement, m.elements)
    println(io, "Model")
    println(io, "  $(length(m.nodes)) nodes, $(length(m.elements)) elements " *
                "($nf frame, $nt truss$(length(m.elements) - nf - nt > 0 ? ", $(length(m.elements) - nf - nt) other" : ""))")
    println(io, "  $(length(m.loads)) loads, $(length(m.springs)) nodal springs")
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
