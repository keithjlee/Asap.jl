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
