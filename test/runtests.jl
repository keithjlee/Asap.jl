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
    analyze!(ex38)

    # value in textbook
    d_textbook = [0.21552, -0.13995]

    @test round.(d_textbook, sigdigits = 4) ≈ round.(n4.disp, sigdigits = 4)


    #2D frame test: Example 6.6
    # in kips, in

    # nodes
    n1 = Node([0., 0.] .* 12, :frame, :fixed)
    n2 = Node([10., 20.] .* 12, :frame, :free)
    n3 = Node([30., 20.] .* 12, :frame, :fixed)

    nodes = [n1, n2, n3]

    # elements
    E = 29e3 #ksi
    A = 11.8 # in^2
    I = 310. # in^4

    e1 = Element(nodes, [1,2], E, A, I)
    e2 = Element(nodes, [2,3], E, A, I)

    elements = [e1, e2]

    # loads
    l1 = Load(nodes, n2.position, [0., -60., -750.])

    loads = [l1]

    # assembly + analysis
    ex66 = Structure(nodes, elements, loads)
    analyze!(ex66)

    d_textbook = [0.021302, -0.06732, -0.0025499]

    @test round.(d_textbook, sigdigits = 4) ≈ round.(n2.disp, sigdigits = 4)

    # 3D truss test: Example 8.1
    # in kips, in

    # nodes
    n1 = Node([0., 0., 0.] .* 12, :truss, :fixed)
    n2 = Node([18., 0., 0.] .* 12, :truss, :fixed)
    n3 = Node([12., 16., 0.] .* 12, :truss, :fixed)
    n4 = Node([-6., 16., 0.] .* 12, :truss, :fixed)
    n5 = Node([6., 8., 24.] .* 12, :truss, :free)

    nodes = [n1, n2, n3, n4, n5]

    # elements
    E = 10e3 #ksi
    A = 8.4 #in^2

    e1 = Element(nodes, [1,5], E, A)
    e2 = Element(nodes, [2,5], E, A)
    e3 = Element(nodes, [3,5], E, A)
    e4 = Element(nodes, [4,5], E, A)

    elements = [e1, e2, e3, e4]

    # loads
    l1 = Load(nodes, n5.position, [0., 50., -100.])

    loads = [l1]

    ex81 = Structure(nodes, elements, loads)
    analyze!(ex81)

    d_textbook = [.10913, .57202, -.12104]

    @test round.(d_textbook, sigdigits = 4) ≈ round.(n5.disp, sigdigits = 4)

    # 3D frame test: Example 8.3 Ferreira "MATLAB Codes for Finite Element Analysis"
    # in kN, m

    # nodes
    n1 = Node([0., 0., 0.], :frame, :free)
    n2 = Node([3., 0., 0.], :frame, :fixed)
    n3 = Node([0., 0., -3.], :frame, :fixed)
    n4 = Node([0., -4., 0.], :frame, :fixed)

    nodes = [n1, n2, n3, n4]

    # elements
    E = 210e6
    A = 0.02
    Iy = 10e-5
    Iz = 20e-5
    J = 5e-5
    G = 84e6

    e1 = Element(nodes, [1,2], E, A, G, Iz, Iy, J)
    e2 = Element(nodes, [1,3], E, A, G, Iz, Iy, J)
    e3 = Element(nodes, [1,4], E, A, G, Iz, Iy, J)

    elements = [e1, e2, e3]

    # loads
    l = Load(nodes, n1.position, [-10., 0., 20., 0., 0., 0.])

    loads = [l]

    ex83 = Structure(nodes, elements, loads)
    analyze!(ex83)
    
    d_textbook = [-7.05e-6, -7e-8, 1.418e-5, 1.45e-6, 1.75e-6, 1.14e-6]

    @test round.(d_textbook, digits = 8) ≈ round.(n1.disp, digits = 8)
end
