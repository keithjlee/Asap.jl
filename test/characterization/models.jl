# Canonical model builders shared by the fixture generator and the
# characterization tests. These pin the *current* (v0.2.x) behavior of Asap
# ahead of the v1.0 modernization — see docs/MODERNIZATION.md.
#
# Each builder returns an unsolved model. Do not modify existing builders:
# changing a model invalidates its pinned fixtures. Add new builders instead.

using Asap

# 2D truss: Kassimali "Matrix Analysis of Structures 2e" Example 3.9
function char_truss_2d()
    n1 = TrussNode([0.0, 0.0, 0.0], :fixed)
    n2 = TrussNode([10.0, 0.0, 0.0], :fixed)
    n3 = TrussNode([0.0, 8.0, 0.0], :yfree)
    n4 = TrussNode([6.0, 8.0, 0.0], :free)
    nodes = [n1, n2, n3, n4]

    sec = TrussSection(4e3 / 1e6, 70.0)

    elements = [
        TrussElement(n1, n3, sec),
        TrussElement(n3, n4, sec),
        TrussElement(n1, n4, sec),
        TrussElement(n2, n3, sec),
        TrussElement(n2, n4, sec),
    ]

    loads = [
        NodeForce(n3, [0.0, -400.0, 0.0]),
        NodeForce(n4, [800.0, -400.0, 0.0]),
    ]

    model = TrussModel(nodes, elements, loads)
    planarize!(model)
    return model
end

# 2D frame: Kassimali Example 6.6 (kips, in)
function char_frame_2d()
    n1 = Node([0.0, 0.0, 0.0], :fixed)
    n2 = Node([10.0, 20.0, 0.0] .* 12, :free)
    n3 = Node([30.0, 20.0, 0.0] .* 12, :fixed)
    nodes = [n1, n2, n3]

    sec = Section(11.8, 29e3, 1.0, 310.0, 310.0, 1.0)

    e1 = Element(n1, n2, sec)
    e1.Ψ = 0.0
    e2 = Element(n2, n3, sec)
    e2.Ψ = 0.0
    elements = [e1, e2]

    loads = [
        PointLoad(e1, 0.5, [0.0, -90.0, 0.0]),
        NodeMoment(n2, [0.0, 0.0, -125.0 * 12]),
        LineLoad(e2, [0.0, -1.5 / 12, 0.0]),
    ]

    model = Model(nodes, elements, loads)
    planarize!(model)
    return model
end

# 3D truss: Kassimali Example 8.1 (kips, in)
function char_truss_3d()
    sec = TrussSection(8.4, 10e3)

    n1 = TrussNode([-6.0, 0.0, 8.0] .* 12, :fixed)
    n2 = TrussNode([12.0, 0.0, 8.0] .* 12, :fixed)
    n3 = TrussNode([6.0, 0.0, -8.0] .* 12, :fixed)
    n4 = TrussNode([-12.0, 0.0, -8.0] .* 12, :fixed)
    n5 = TrussNode([0.0, 24.0, 0.0] .* 12, :free)
    nodes = [n1, n2, n3, n4, n5]

    elements = [
        TrussElement(n1, n5, sec),
        TrussElement(n2, n5, sec),
        TrussElement(n3, n5, sec),
        TrussElement(n4, n5, sec),
    ]

    loads = [NodeForce(n5, [0.0, -100.0, -50.0])]

    return TrussModel(nodes, elements, loads)
end

# 3D frame with roll angles and torsion: Kassimali Example 8.4 (kips, in)
function char_frame_3d()
    n1 = Node([0.0, 0.0, 0.0], :free)
    n2 = Node([-240.0, 0.0, 0.0], :fixed)
    n3 = Node([0.0, -240.0, 0.0], :fixed)
    n4 = Node([0.0, 0.0, -240.0], :fixed)
    nodes = [n1, n2, n3, n4]

    sec = Section(32.9, 29e3, 11.5e3, 716.0, 236.0, 15.1)

    e1 = Element(n2, n1, sec)
    e1.Ψ = 0.0
    e2 = Element(n3, n1, sec)
    e2.Ψ = pi / 2
    e3 = Element(n4, n1, sec)
    e3.Ψ = pi / 6
    elements = [e1, e2, e3]

    loads = [
        LineLoad(e1, [0.0, -3.0 / 12, 0.0]),
        NodeMoment(n1, [-150.0 * 12, 0.0, 150.0 * 12]),
    ]

    return Model(nodes, elements, loads)
end

# Cantilever with a tip force (12000.org stiffness matrix report, Ex 1)
function char_cantilever_tip()
    n1 = Node([0.0, 0.0, 0.0], :fixed)
    n2 = Node([144.0, 0.0, 0.0], :free)
    sec = Section(1.0, 30e6, 1.0, 57.1, 57.1, 1.0)
    e = Element(n1, n2, sec)
    loads = [NodeForce(n2, [0.0, -400.0, 0.0])]

    model = Model([n1, n2], [e], loads)
    planarize!(model)
    return model
end

# Cantilever with a midspan point load (12000.org Ex 2)
function char_cantilever_midpoint()
    n1 = Node([0.0, 0.0, 0.0], :fixed)
    n2 = Node([144.0, 0.0, 0.0], :free)
    sec = Section(1.0, 30e6, 1.0, 57.1, 57.1, 1.0)
    e = Element(n1, n2, sec)
    loads = [PointLoad(e, 0.5, [0.0, -400.0, 0.0])]

    model = Model([n1, n2], [e], loads)
    planarize!(model)
    return model
end

const CHAR_RELEASES = [:fixedfixed, :fixedfree, :freefixed, :freefree, :joist]

# 3D portal frame whose beam carries the given release. Columns rise in +Z,
# beam spans +X; asymmetric loads (distributed + point + lateral) exercise
# both bending planes, shear, and torsion. Also the canonical model for the
# AsapToolkit InternalForces fixtures used to validate Phase 3 force recovery.
function char_portal_frame(release::Symbol)
    H, L = 3000.0, 5000.0                     # mm-ish proportions

    n1 = Node([0.0, 0.0, 0.0], :fixed)
    n2 = Node([0.0, 0.0, H], :free)
    n3 = Node([L, 0.0, H], :free)
    n4 = Node([L, 0.0, 0.0], :fixed)
    nodes = [n1, n2, n3, n4]

    sec = Section(1e4, 200.0, 77.0, 8e7, 3e7, 5e6)

    col1 = Element(n1, n2, sec, :column)
    beam = Element(n2, n3, sec, :beam; release=release)
    col2 = Element(n4, n3, sec, :column)
    elements = [col1, beam, col2]

    loads = [
        LineLoad(beam, [0.0, 0.0, -2.0]),
        PointLoad(beam, 0.4, [0.0, 3.0, -10.0]),
        NodeForce(n2, [5.0, 0.0, 0.0]),
    ]

    return Model(nodes, elements, loads)
end

# Skewed single element (non-axis-aligned, both nodes fixed) used to pin
# fixed-end force vectors q()/q_local(), local/global stiffness matrices,
# and the transformation matrix per release type.
function char_fef_element(release::Symbol)
    n1 = Node([0.0, 0.0, 0.0], :fixed)
    n2 = Node([3000.0, 1000.0, 2000.0], :fixed)
    sec = Section(1e4, 200.0, 77.0, 8e7, 3e7, 5e6)

    e = Element(n1, n2, sec; release=release)

    loads = [
        LineLoad(e, [0.4, -1.2, 0.9]),
        PointLoad(e, 0.3, [110.0, -70.0, 25.0]),
    ]

    model = Model([n1, n2], [e], loads)
    return model
end

# GravityLoad model — documents the known bug: q_local(::GravityLoad)
# references nonexistent fields, so processing a model with a GravityLoad
# throws. Pinned as @test_broken; when the modernization replaces GravityLoad
# with SelfWeight, this becomes a real test.
function char_gravity_model()
    n1 = Node([0.0, 0.0, 0.0], :fixed)
    n2 = Node([3000.0, 0.0, 0.0], :fixed)
    sec = Section(1e4, 200.0, 77.0, 8e7, 3e7, 5e6, 7.85e-9)
    e = Element(n1, n2, sec)
    loads = [GravityLoad(e, 9810.0)]
    return Model([n1, n2], [e], loads)
end
