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
    println(io, "  Ψ = $(el.Ψ)  [rad]  (section roll about the element axis)")
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
