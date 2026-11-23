"""
    AbstractShellSection{T<:Real}

Supertype of all shell (surface) sections — the through-thickness
description of a [`ShellElement`](@ref).

Like [`AbstractSection`](@ref) for members, a shell section's entire contract
with the analysis layer is a **rigidity accessor interface**: the
thickness-integrated constitutive matrices a linear flat-shell stiffness
consumes, plus mass per unit area. Kernels never read struct fields.

| Accessor | Meaning | Dimension |
|---|---|---|
| [`membrane_rigidity`](@ref) | in-plane stress resultant per unit strain, `N = A·ε` (3×3) | [force/length] |
| [`bending_rigidity`](@ref)  | moment per unit curvature, `M = D·κ` (3×3) | [force·length] |
| [`shear_rigidity`](@ref)    | transverse shear per unit shear strain, `Q = S·γ` (2×2) | [force/length] |
| [`drilling_rigidity`](@ref) | penalty tying the drilling rotation to the in-plane rotation of the displacement field | [force/length] |
| [`areal_mass`](@ref)        | mass per unit area | [mass/length²] |

Voigt ordering everywhere: `(xx, yy, xy)` with ENGINEERING shear strain
`γxy = ∂u/∂y + ∂v/∂x` and twist `κxy = ∂βx/∂y + ∂βy/∂x`; transverse shear
`(xz, yz)`. Matrices are expressed in the element's local axes.

No membrane–bending coupling block (`B` in classical lamination theory) is
carried: sections are symmetric through the thickness. Eccentric or layered
sections are out of scope for linear 1.3.
"""
abstract type AbstractShellSection{T<:Real} end

"""
    ShellSection{T} <: AbstractShellSection{T}

A homogeneous, isotropic shell of uniform thickness — the standard section
for concrete walls, slabs, and steel plate.

Rigidities are derived on demand from the material and thickness `t`
(plane-stress Reissner–Mindlin plate theory):

- membrane: `A = t·Q`, with `Q = E/(1−ν²)·[1 ν 0; ν 1 0; 0 0 0] + G·e₃e₃ᵀ`
- bending: `D = t³/12·Q`
- transverse shear: `S = κ·G·t·I₂`, shear correction `κ = 5/6`
- drilling penalty: `γ·t` with `γ = G` (Hughes & Brezzi 1989; Ibrahimbegović,
  Taylor & Wilson 1990 — results are insensitive to γ over several orders
  of magnitude)

The in-plane shear term uses the material's stored `G`, which equals the
isotropic value `E/(2(1+ν))` whenever the material was built with
`Material(E, ρ, ν)`.

# Fields
- `material::Material{T}`: the shell material (E, G, ρ, ν)
- `t::T`: thickness [length]

# Examples
```julia-repl
julia> concrete = Material(4.0e6, 0.15, 0.2);    # kip, ft: E in kip/ft²

julia> wall = ShellSection(concrete, 8 / 12)        # 8 in wall
```

For cracked/effective or orthotropic stiffness use
[`ShellRigiditySection`](@ref), including its stiffness-modifier constructor.
"""
struct ShellSection{T} <: AbstractShellSection{T}
    material::Material{T}
    t::T

    function ShellSection(material::Material{T}, t::Real) where {T}
        Tp = promote_type(T, typeof(t))
        t > 0 || throw(ArgumentError("shell thickness must be positive (got $t)"))
        return new{Tp}(convert(Material{Tp}, material), Tp(t))
    end
end

"""
    ShellRigiditySection{T} <: AbstractShellSection{T}

A shell section defined **directly by its thickness-integrated rigidities**
— the shell analogue of [`RigiditySection`](@ref). Use it for cracked or
effective concrete stiffness, orthotropic or equivalent plates (ribbed,
hollow-core, or composite decks), and optimization parameterizations.

# Fields
- `A::SMatrix{3,3,T}`: membrane rigidity `N = A·ε` [force/length]
- `D::SMatrix{3,3,T}`: bending rigidity `M = D·κ` [force·length]
- `S::SMatrix{2,2,T}`: transverse shear rigidity `Q = S·γ` [force/length]
- `γt::T`: drilling penalty rigidity [force/length] (default: `A[3,3]`, the
  in-plane shear rigidity, i.e. `γ = G`)
- `ρt::T`: mass per unit area [mass/length²]

# Constructors
    ShellRigiditySection(A, D, S, ρt; γt = A[3,3])
    ShellRigiditySection(sec::ShellSection; f11 = 1, f22 = 1, f12 = 1,
                         m11 = 1, m22 = 1, m12 = 1, v13 = 1, v23 = 1, mass = 1)

The second form applies **stiffness modifiers** (the property-modifier
convention of common building-analysis software) to an isotropic section:

- membrane: `A[1,1]·f11`, `A[2,2]·f22`, `A[3,3]·f12`, and the Poisson
  coupling `A[1,2]·√(f11·f22)` (keeps `A` positive definite);
- bending: `D[1,1]·m11`, `D[2,2]·m22`, `D[3,3]·m12`, coupling `·√(m11·m22)`;
- transverse shear: `S[1,1]·v13`, `S[2,2]·v23`;
- the drilling penalty follows the modified in-plane shear rigidity;
- `mass` scales the areal mass.

Direction 1 is the element's local x-axis. A cracked wall under ACI 318-19
Table 6.6.3.1.1(a), for instance, takes `m11 = m22 = m12 = 0.35` (or `0.70`
uncracked); which modifiers a code prescribes is the user's (or a design
kernel's) decision, not Asap's.

# Examples
```julia-repl
julia> wall = ShellSection(Material(4.0e6, 0.15, 0.2), 8 / 12);

julia> cracked = ShellRigiditySection(wall; m11 = 0.35, m22 = 0.35, m12 = 0.35)
```
"""
struct ShellRigiditySection{T} <: AbstractShellSection{T}
    A::SMatrix{3,3,T,9}
    D::SMatrix{3,3,T,9}
    S::SMatrix{2,2,T,4}
    γt::T
    ρt::T
end

function ShellRigiditySection(A::AbstractMatrix, D::AbstractMatrix, S::AbstractMatrix,
    ρt::Real; γt::Real = A[3, 3])
    T = promote_type(eltype(A), eltype(D), eltype(S), typeof(ρt), typeof(γt))
    As, Ds, Ss = SMatrix{3,3,T}(A), SMatrix{3,3,T}(D), SMatrix{2,2,T}(S)
    for (M, name) in ((As, "membrane"), (Ds, "bending"), (Ss, "shear"))
        M ≈ M' || throw(ArgumentError("$name rigidity matrix must be symmetric"))
        isposdef(Symmetric(Matrix(M))) ||
            throw(ArgumentError("$name rigidity matrix must be positive definite"))
    end
    γt > 0 || throw(ArgumentError("drilling rigidity γt must be positive"))
    return ShellRigiditySection{T}(As, Ds, Ss, T(γt), T(ρt))
end

function ShellRigiditySection(sec::ShellSection{T}; f11::Real = 1, f22::Real = 1, f12::Real = 1,
    m11::Real = 1, m22::Real = 1, m12::Real = 1, v13::Real = 1, v23::Real = 1,
    mass::Real = 1) where {T}
    A = membrane_rigidity(sec)
    D = bending_rigidity(sec)
    S = shear_rigidity(sec)
    fA = SMatrix{3,3}(f11, sqrt(f11 * f22), 0, sqrt(f11 * f22), f22, 0, 0, 0, f12)
    fD = SMatrix{3,3}(m11, sqrt(m11 * m22), 0, sqrt(m11 * m22), m22, 0, 0, 0, m12)
    fS = SMatrix{2,2}(v13, 0, 0, v23)
    Am = A .* fA
    return ShellRigiditySection(Am, D .* fD, S .* fS, areal_mass(sec) * mass; γt = Am[3, 3])
end

# ─────────────────────────────────────────────────────────────────────────────
# The shell rigidity accessor interface
# ─────────────────────────────────────────────────────────────────────────────

# plane-stress constitutive matrix of the material [force/length²]
function _plane_stress(m::Material)
    c = m.E / (1 - m.ν^2)
    z = zero(c)
    return SMatrix{3,3}(c, m.ν * c, z, m.ν * c, c, z, z, z, m.G)
end

"""
    membrane_rigidity(section) -> SMatrix{3,3}

Membrane (in-plane) rigidity `A` with `(Nxx, Nyy, Nxy) = A·(εxx, εyy, γxy)`
[force/length]. `ShellSection`: `t·Q` (plane stress).
"""
membrane_rigidity(s::ShellSection) = s.t * _plane_stress(s.material)
membrane_rigidity(s::ShellRigiditySection) = s.A

"""
    bending_rigidity(section) -> SMatrix{3,3}

Bending rigidity `D` with `(Mxx, Myy, Mxy) = D·(κxx, κyy, κxy)`
[force·length]. `ShellSection`: `t³/12·Q`.
"""
bending_rigidity(s::ShellSection) = (s.t^3 / 12) * _plane_stress(s.material)
bending_rigidity(s::ShellRigiditySection) = s.D

"""
    shear_rigidity(section) -> SMatrix{2,2}

Transverse shear rigidity `S` with `(Qx, Qy) = S·(γxz, γyz)` [force/length].
`ShellSection`: `κ·G·t·I` with `κ = 5/6`.
"""
function shear_rigidity(s::ShellSection)
    k = (5 * s.material.G * s.t) / 6
    return SMatrix{2,2}(k, zero(k), zero(k), k)
end
shear_rigidity(s::ShellRigiditySection) = s.S

"""
    drilling_rigidity(section) -> T

Drilling penalty `γ·t` [force/length] of the Hughes–Brezzi / ITW membrane:
weights `∫(skew∇u − ψ)²`, tying the nodal drilling rotation `ψ` to the
in-plane rotation of the displacement field. `ShellSection`: `G·t`.
"""
drilling_rigidity(s::ShellSection) = s.material.G * s.t
drilling_rigidity(s::ShellRigiditySection) = s.γt

"""
    areal_mass(section) -> T

Mass per unit area [mass/length²]. Multiply by gravitational acceleration for
self-weight per unit area. `ShellSection`: `ρ·t`.
"""
areal_mass(s::ShellSection) = s.material.ρ * s.t
areal_mass(s::ShellRigiditySection) = s.ρt

"""
    thickness(section) -> T

Physical thickness [length] (`ShellSection` only — a rigidity section has no
unique thickness).
"""
thickness(s::ShellSection) = s.t
