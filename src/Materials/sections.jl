"""
    AbstractSection{T<:Real}

Supertype of all cross-sections.

A section's entire contract with the analysis layer is the **rigidity
accessor interface** — five functions returning the only quantities a linear
member stiffness matrix (or self-weight computation) ever consumes:

| Accessor | Meaning | Dimension |
|---|---|---|
| [`EA`](@ref)  | axial rigidity                       | [force] |
| [`EIx`](@ref) | flexural rigidity, strong axis       | [force·length²] |
| [`EIy`](@ref) | flexural rigidity, weak axis         | [force·length²] |
| [`GJ`](@ref)  | torsional rigidity                   | [force·length²] |
| [`ρA`](@ref)  | mass per unit length                 | [mass/length] |

Kernels never read struct fields — they call these accessors. Any type
implementing all five is a valid section, which is what lets a classical
geometric [`Section`](@ref) and an effective-stiffness
[`RigiditySection`](@ref) (cracked concrete, homogenized members) be used
interchangeably anywhere.
"""
abstract type AbstractSection{T<:Real} end

"""
    Section{T} <: AbstractSection{T}

A cross-section defined by its geometry and its [`Material`](@ref) — the
standard modeling workflow when section dimensions and material are known
separately.

Rigidities are *derived* on demand: `EA(sec) = sec.material.E * sec.A`, etc.
If you instead know effective rigidities directly (cracked concrete,
transformed sections, optimization parameterizations), use
[`RigiditySection`](@ref).

# Fields
- `material::Material{T}`: the section's material (provides E, G, ρ)
- `A::T`: cross-sectional area [length²]
- `Ix::T`: second moment of area about the **strong** bending axis (local x-axis
  convention of the legacy library, paired with bending in the local x–y
  plane) [length⁴]
- `Iy::T`: second moment of area about the **weak** bending axis [length⁴]
- `J::T`: St. Venant torsional constant [length⁴]

# Constructors
    Section(material, A, Ix, Iy, J)   # full frame section
    Section(material, A)              # axial-only (truss) section: Ix = Iy = J = 0

The axial-only form replaces the legacy `TrussSection`: a truss element only
ever queries `EA`/`ρA`, so zero flexural/torsional properties are simply
never read.

# Examples
```julia-repl
julia> steel = Material(200e6, 77e6, 8.0, 0.3);

julia> w = Section(steel, 1e-2, 8e-5, 3e-5, 5e-7)   # a wide-flange, kN–m

julia> bar = Section(steel, 5e-3)                    # truss bar, axial only
```
"""
struct Section{T} <: AbstractSection{T}
    material::Material{T}
    A::T  # area
    Ix::T # strong-axis second moment of area
    Iy::T # weak-axis second moment of area
    J::T  # torsional constant

    function Section(material::Material{T}, A::Real, Ix::Real, Iy::Real, J::Real) where {T}
        A, Ix, Iy, J = promote(T(A), T(Ix), T(Iy), T(J))
        return new{T}(material, A, Ix, Iy, J)
    end
end

Section(material::Material{T}, A::Real) where {T} = Section(material, A, zero(T), zero(T), zero(T))

"""
    RigiditySection{T} <: AbstractSection{T}

A cross-section defined **directly by its rigidities** — the products the
stiffness matrix actually consumes — rather than by geometry and material.

This is the natural abstraction wherever an "equivalent stiffness" is the
meaningful quantity:

- **Cracked/effective concrete sections**: build from ACI 318-19 stiffness
  modifiers (e.g. EI_eff = 0.35·Ec·Ig for beams, 0.70·Ec·Ig for columns), a
  Branson/Bischoff effective moment of inertia, or a full transformed-section
  analysis. For linear analysis this is *exactly* as expressive as any
  equivalent-E or transformed-I formulation, because only the products enter
  the equations.
- **Homogenized/composite members**, where no single (E, I) pair exists.
- **Optimization parameterizations** where rigidities are the design variables.

Because rigidities are stored, there is no meaningful area or material to
derive mass from — so mass per unit length `ρA` is stored explicitly (do not
back it out of the rigidities).

# Fields
- `EA::T`: axial rigidity — force per unit axial strain [force]
- `EIx::T`: flexural rigidity about the strong axis [force·length²]
- `EIy::T`: flexural rigidity about the weak axis [force·length²]
- `GJ::T`: torsional rigidity [force·length²]
- `ρA::T`: mass per unit length, for self-weight and dynamics [mass/length]

# Examples
```julia-repl
julia> Ec, Ig, A, ρc = 30e6, 8e-4, 0.12, 2.4e3;

julia> cracked_beam = RigiditySection(Ec * A, 0.35 * Ec * Ig, 0.35 * Ec * Ig,
                                      0.1 * Ec * Ig, ρc * A)
```

See also [`Section`](@ref), [`AbstractSection`](@ref).
"""
struct RigiditySection{T} <: AbstractSection{T}
    EA::T  # axial rigidity
    EIx::T # strong-axis flexural rigidity
    EIy::T # weak-axis flexural rigidity
    GJ::T  # torsional rigidity
    ρA::T  # mass per unit length

    function RigiditySection(EA::Real, EIx::Real, EIy::Real, GJ::Real, ρA::Real)
        args = promote(EA, EIx, EIy, GJ, ρA)
        return new{eltype(args)}(args...)
    end
end

# ─────────────────────────────────────────────────────────────────────────────
# The rigidity accessor interface — the section contract consumed by kernels
# ─────────────────────────────────────────────────────────────────────────────

"""
    EA(section) -> T

Axial rigidity of the section: the force required to produce unit axial
strain [force]. `Section` derives it as `E·A`; `RigiditySection` stores it.
"""
EA(s::Section) = s.material.E * s.A
EA(s::RigiditySection) = s.EA

"""
    EIx(section) -> T

Flexural rigidity about the section's **strong** axis [force·length²].
Governs bending in the element's local x–y plane. `Section` derives it as
`E·Ix`; `RigiditySection` stores it.
"""
EIx(s::Section) = s.material.E * s.Ix
EIx(s::RigiditySection) = s.EIx

"""
    EIy(section) -> T

Flexural rigidity about the section's **weak** axis [force·length²].
Governs bending in the element's local x–z plane. `Section` derives it as
`E·Iy`; `RigiditySection` stores it.
"""
EIy(s::Section) = s.material.E * s.Iy
EIy(s::RigiditySection) = s.EIy

"""
    GJ(section) -> T

St. Venant torsional rigidity [force·length²]: torque per unit rate of twist.
`Section` derives it as `G·J`; `RigiditySection` stores it.
"""
GJ(s::Section) = s.material.G * s.J
GJ(s::RigiditySection) = s.GJ

"""
    ρA(section) -> T

Mass per unit length of the member [mass/length]. Multiply by gravitational
acceleration for self-weight per unit length. `Section` derives it as `ρ·A`;
`RigiditySection` stores it explicitly.
"""
ρA(s::Section) = s.material.ρ * s.A
ρA(s::RigiditySection) = s.ρA
