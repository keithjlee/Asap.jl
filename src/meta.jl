"""
Define a material
"""
struct Material
    E::Float64 #young's modulus
    G::Float64 #shear modulus
    ρ::Float64 #density
    ν::Float64 #poisson's ratio

    function Material(E::Float64, G::Float64, ρ::Float64, ν::Float64)
        return new(E, G, ρ, ν)
    end
end

const Steel_Nmm = Material(200e3, 77e3, 8e-5, 0.3)
const Steel_kNm = Material(200e6, 77e6, 80., 0.3)


abstract type AbstractSection end

"""
Define a full cross section
"""
mutable struct Section <: AbstractSection
    A::Float64 # area
    E::Float64 # young's modulus
    G::Float64 # shear modulus
    Izz::Float64 # strong axis I
    Iyy::Float64 # weak axis I
    J::Float64 # torsional constant
    ρ::Float64 # density

    function Section(A::Float64, E::Float64, G::Float64, Izz::Float64, Iyy::Float64, J::Float64)
        println("Density ρ not provided, defaults to ρ = 1; self-weight analysis will not function")
        return new(A, E, G, Izz, Iyy, J, 1.)
    end

    function Section(mat::Material, A::Float64, Izz::Float64, Iyy::Float64, J::Float64)
        return new(A, mat.E, mat.G, Izz, Iyy, J, mat.ρ)
    end

end

"""
Define a simplified cross section for trusses
"""
mutable struct TrussSection <: AbstractSection
    A::Float64
    E::Float64
    ρ::Float64

    function TrussSection(A::Float64, E::Float64)
        println("Density ρ not provided, defaults to ρ = 1; self-weight analysis will not function")
        return new(A, E, 1.)
    end

    function TrussSection(mat::Material, A::Float64)
        return new(A, mat.E, mat.ρ)
    end
end