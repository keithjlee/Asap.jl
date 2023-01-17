using Asap
using Test
using LinearAlgebra

@testset "Asap.jl" begin
    #tol
    tol = 0.1
    # 2D truss test: Example 3.9 from Kassimali "Matrix Analysis of Structures 2e"
    # in kips, ft

    n1 = TrussNode([0., 0., 0.], :fixed)
    n2 = TrussNode([10., 0., 0.], :fixed)
    n3 = TrussNode([0., 8., 0.], :yfree)
    n4 = TrussNode([6., 8., 0.], :free)

    nodes = [n1, n2, n3, n4]
    planarize!(nodes)

    E = 70. #kN/m^2
    A = 4e3 / 1e6 #m^2

    sec = TrussSection(A, E)

    e1 = TrussElement(nodes, [1,3], sec)
    e2 = TrussElement(nodes, [3,4], sec)
    e3 = TrussElement(nodes, [1,4], sec)
    e4 = TrussElement(nodes, [2,3], sec)
    e5 = TrussElement(nodes, [2,4], sec)

    elements = [e1,e2,e3,e4,e5]

    l1 = NodeForce(n3, [0., -400., 0.])
    l2 = NodeForce(n4, [800., -400., 0.])

    loads = [l1, l2]

    model = TrussModel(nodes, elements, loads)
    solve!(model)

    reactions = model.reactions[model.fixedDOFs]
    reactions2d = reactions[reactions .!= 0]

    reactions_textbook = [-.57994,
        320.82,
        -298.39,
        479.17,
        -501.05]

    err = norm(reactions_textbook .- reactions2d)

    @test err <= tol


    #2D frame test: Example 6.6
    # in kips, in

    n1 = Node([0., 0., 0.], :fixed)
    n2 = Node([10., 20., 0.] .* 12, :free)
    n3 = Node([30., 20., 0.] .* 12, :fixed)
    nodes = [n1, n2, n3]
    planarize!(nodes)

    E = 29e3
    A = 11.8
    I = 310.

    sec = Section(A, E, 1., I, I, 1.)

    e1 = Element(nodes, [1,2], sec)
    e1.Ψ = 0.
    e2 = Element(nodes, [2,3], sec)
    e2.Ψ = 0.

    elements = [e1, e2]

    l1 = PointLoad(e1, 0.5, [0., -90., 0.])
    l2 = NodeMoment(n2, [0., 0., -125 * 12.])
    l3 = LineLoad(e2, [0., -1.5/12, 0.])

    loads = [l1, l2, l3]

    model = Model(nodes, elements, loads)
    solve!(model)

    reactions = model.reactions[model.fixedDOFs]
    reactions = reactions[reactions .!= 0]

    reactions_textbook = [30.371,
        102.09,
        1216.,
        -30.372,
        17.913,
        -854.07]

    err = norm(reactions_textbook .- reactions)

    @test err <= tol

    # 3D truss test: Example 8.1
    # in kips, in

    E = 10e3
    A = 8.4
    sec = TrussSection(A, E)

    n1 = TrussNode([-6., 0., 8.] .* 12, :fixed)
    n2 = TrussNode([12., 0., 8.] .* 12, :fixed)
    n3 = TrussNode([6., 0., -8.] .* 12, :fixed)
    n4 = TrussNode([-12., 0., -8.] .* 12, :fixed)
    n5 = TrussNode([0., 24., 0.] .* 12, :free)

    nodes = [n1, n2, n3, n4, n5]

    e1 = TrussElement(nodes, [1,5], sec)
    e2 = TrussElement(nodes, [2,5], sec)
    e3 = TrussElement(nodes, [3,5], sec)
    e4 = TrussElement(nodes, [4,5], sec)

    elements = [e1, e2, e3, e4]

    l1 = NodeForce(n5, [0., -100., -50.])
    loads = [l1]

    model = TrussModel(nodes, elements, loads)
    solve!(model)

    reactions = model.reactions[model.fixedDOFs]

    reactions_textbook = [-5.5581,
        -22.232,
        7.4108,
        1.3838,
        -2.7677,
        0.92255,
        -19.442,
        77.768,
        25.923,
        23.616,
        47.232,
        15.744
        ]

    err = norm(reactions_textbook .- reactions)

    @test err <= tol

    # 3D frame test: Example 8.4
    # in kN, m

    n1 = Node([0., 0., 0.], :free)
    n2 = Node([-240., 0., 0.], :fixed)
    n3 = Node([0., -240., 0.], :fixed)
    n4 = Node([0., 0., -240.], :fixed)

    nodes = [n1, n2, n3, n4]


    E = 29e3
    G = 11.5e3 
    A = 32.9
    Iz = 716.
    Iy = 236.
    J = 15.1

    sec = Section(A, E, G, Iz, Iy, J)

    e1 = Element(nodes, [2,1], sec)
    e1.Ψ = 0.
    e2 = Element(nodes, [3,1], sec)
    e2.Ψ  = pi/2
    e3 = Element(nodes, [4,1], sec)
    e3.Ψ = pi/6

    elements = [e1, e2, e3]

    l1 = LineLoad(e1, [0., -3/12, 0.])
    l2 = NodeMoment(n1, [-150. * 12, 0., 150. * 12])

    loads = [l1, l2]

    model = Model(nodes, elements, loads)
    solve!(model)

    reactions = model.reactions[model.fixedDOFs]

    reactions_textbook = [5.3757,
        44.106,
        -0.74272,
        2.1722,
        58.987,
        2330.5,
        -4.6249,
        11.117,
        -6.4607,
        -515.55,
        -0.76472,
        369.67,
        -0.75082,
        4.7763,
        7.2034,
        -383.5,
        -60.166,
        -4.702]

    err = norm(reactions_textbook .- reactions)

    @test err <= tol
end
