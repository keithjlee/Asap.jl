using Asap
using Test

@testset "Asap.jl" begin
    # 2D truss test: Example 3.8 from Kassimali "Matrix Analysis of Structures 2e"
    # in kips, ft

    #nodes
    n1 = Node([0., 0.] .* 12, :truss, :fixed)
    n2 = Node([12., 0.] .* 12, :truss, :fixed)
    n3 = Node([24., 0.] .* 12, :truss, :fixed)
    n4 = Node([12., 16.] .* 12, :truss, :free)

    nodes = [n1, n2, n3, n4]

    # elements
    E = 29e3 #ksi

    e1 = Element(nodes, [1,4], E, 8.)
    e2 = Element(nodes, [2,4], E, 6.)
    e3 = Element(nodes, [3,4], E, 8.)

    elements = [e1, e2, e3]

    # load
    l1 = Load(nodes, n4.position, [150., -300.])

    loads = [l1]

    # assemble
    ex38 = Structure(nodes, elements, loads)

    # Analyze
    analyze(ex38)

    # value in textbook
    d_textbook = [0.21552, -0.13995]

    @test d_textbook ≈ round.(n4.disp, digits = 5)


    #2D frame test: Example 7.3 from Ferreira "Matlab Codes for Finite Element Analysis"
    
end
