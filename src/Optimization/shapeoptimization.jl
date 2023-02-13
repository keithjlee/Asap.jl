"""
makeIJ: Make the row and column indices for the global stiffness matrix in sparse CSC format
"""
function makeIJ(model::AbstractModel, n::Integer)
    #initialize
    I = Vector{Int64}()
    J = Vector{Int64}()

    #iterate
    for element in model.elements
        idx = element.globalID
        @inbounds for j = 1:n
            @inbounds for i = 1:n
                push!(I, idx[i])
                push!(J, idx[j])
            end
        end    
    end

    return I, J
end

"""
Function methods for truss and frame models
"""
makeIJ(model::TrussModel) = makeIJ(model, 6)
makeIJ(model::Model) = makeIJ(model, 12)