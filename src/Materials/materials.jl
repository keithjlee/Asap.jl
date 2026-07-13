"""
    Material{T<:Real}

A linear-elastic structural material.

A `Material` carries the four scalar properties that linear member analysis
can ever consume: stiffness in normal action (`E`), stiffness in shear (`G`),
mass density (`ρ`), and the ratio that couples them (`ν`). It is immutable —
to analyze a design change, construct a new `Material`.

Asap is units-agnostic: choose any consistent unit system (e.g. kN–m,
N–mm, kip–in) and use it everywhere. Dimensions below are written in
bracket notation.

# Fields
- `E::T`: Young's (elastic) modulus — normal stress per unit normal strain
  [force/length²]
- `G::T`: shear modulus — shear stress per unit shear strain [force/length²]
- `ρ::T`: mass density — used for self-weight and (future) dynamic analysis
  [mass/length³]
- `ν::T`: Poisson's ratio — transverse contraction per unit axial extension
  [unitless]. For isotropic materials `G = E / (2(1 + ν))`.

# Constructors
    Material(E, G, ρ, ν)     # all four properties explicit
    Material(E, ρ, ν)        # isotropic: G derived as E / (2(1 + ν))

Arguments are promoted to a common scalar type, so `Material(200e6, 80, 0.3)`
and `Material(200e6, 80.0, 0.3)` are equivalent.

# Examples
```julia-repl
julia> steel = Material(200e6, 77e6, 8.0, 0.3)      # kN, m: E in kN/m²
julia> concrete = Material(30e6, 2.4, 0.2)           # G derived isotropically
```

See also [`Section`](@ref), [`RigiditySection`](@ref).
"""
struct Material{T<:Real}
    E::T # Young's modulus
    G::T # shear modulus
    ρ::T # mass density
    ν::T # Poisson's ratio

    function Material(E::T, G::T, ρ::T, ν::T) where {T<:Real}
        return new{T}(E, G, ρ, ν)
    end
end

Material(E::Real, G::Real, ρ::Real, ν::Real) = Material(promote(E, G, ρ, ν)...)
Material(E::Real, ρ::Real, ν::Real) = Material(E, E / (2 * (1 + ν)), ρ, ν)

"Structural steel in N–mm units: E = 200e3 N/mm², G = 77e3 N/mm², ρ = 8e-5 kg/mm³ (legacy value), ν = 0.3."
const Steel_Nmm = Material(200e3, 77e3, 8e-5, 0.3)

"Structural steel in kN–m units: E = 200e6 kN/m², G = 77e6 kN/m², ρ = 80 kg/m³ (legacy value), ν = 0.3."
const Steel_kNm = Material(200e6, 77e6, 80.0, 0.3)
