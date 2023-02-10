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