abstract type AbstractSection end

"""
Define a full cross section
"""
struct Section <: AbstractSection
    A::Float64 # area
    E::Float64 # young's modulus
    G::Float64 # shear modulus
    Ix::Float64 # strong axis I
    Iy::Float64 # weak axis I
    J::Float64 # torsional constant
    ρ::Float64 # density

    """
    Explicitly define section parameters
    """
    function Section(A::Float64, E::Float64, G::Float64, Ix::Float64, Iy::Float64, J::Float64, ρ::Float64 = 1.)
        # println("Density ρ not provided, defaults to ρ = 1; self-weight analysis will not function")
        return new(A, E, G, Ix, Iy, J, ρ)
    end

    """
    Define a section with a material + geometric parameters
    """
    function Section(mat::Material, A::Float64, Ix::Float64, Iy::Float64, J::Float64)
        return new(A, mat.E, mat.G, Ix, Iy, J, mat.ρ)
    end

end

"""
Define a simplified cross section for trusses
"""
struct TrussSection <: AbstractSection
    A::Float64
    E::Float64
    ρ::Float64

    function TrussSection(A::Float64, E::Float64)
        # println("Density ρ not provided, defaults to ρ = 1; self-weight analysis will not function")
        return new(A, E, 1.)
    end

    function TrussSection(mat::Material, A::Float64)
        return new(A, mat.E, mat.ρ)
    end
end