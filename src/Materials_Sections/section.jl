abstract type AbstractSection end

"""
    Section(A::Float64, E::Float64, G::Float64, Ix::Float64, Iy::Float64, J::Float64, ρ::Float64 = 1.)
    Section(mat::Material, A::Float64, Ix::Float64, Iy::Float64, J::Float64)

A cross section assigned to an element.

# Fields
- `A` Area [Distance²]
- `E` Modulus of Elasticity [Force/Distance²]
- `Ix` Nominal strong moment of inertia [Distance⁴]
- `Iy` Nominal weak moment of inertia [Distance⁴]
- `J` Torsional constant [Distance⁴]
- `ρ=1` Density [Mass/Distance³]
"""
struct Section <: AbstractSection
    A::Float64 # area
    E::Float64 # young's modulus
    G::Float64 # shear modulus
    Ix::Float64 # strong axis I
    Iy::Float64 # weak axis I
    J::Float64 # torsional constant
    ρ::Float64 # density

    function Section(A::Float64, E::Float64, G::Float64, Ix::Float64, Iy::Float64, J::Float64, ρ::Float64 = 1.)
        return new(A, E, G, Ix, Iy, J, ρ)
    end

    function Section(mat::Material, A::Float64, Ix::Float64, Iy::Float64, J::Float64)
        return new(A, mat.E, mat.G, Ix, Iy, J, mat.ρ)
    end

end

"""
    TrussSection(A::Float64, E::Float64)
    TrussSection(mat::Material, A::Float64)

A cross section assigned to a truss element.

# Fields
- `A` Area [Distance²]
- `E` Modulus of Elasticity [Force/Distance²]
- `ρ` Density [Mass/Distance³]
"""
struct TrussSection <: AbstractSection
    A::Float64
    E::Float64
    ρ::Float64

    function TrussSection(A::Float64, E::Float64, ρ = 1.)
        return new(A, E, ρ)
    end

    function TrussSection(mat::Material, A::Float64)
        return new(A, mat.E, mat.ρ)
    end
end
