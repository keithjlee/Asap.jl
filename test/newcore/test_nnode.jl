# N-node element path (shells S1): the analysis layer must handle elements
# with any number of nodes. A mock 3-node element — three truss bars around a
# triangle — must reproduce the identical structure built from TrussElements,
# through BOTH assembly paths (in-place and pure/AD).

using StaticArrays, Zygote

mutable struct _TriBars{T,S<:Asap.AbstractSection{T}} <: Asap.AbstractElement{T}
    n1::Asap.Node{T}
    n2::Asap.Node{T}
    n3::Asap.Node{T}
    section::S
    id::Symbol
    index::Int
end
_TriBars(a, b, c, sec) = _TriBars(a, b, c, sec, :tribars, 0)

Asap.nodes(el::_TriBars) = (el.n1, el.n2, el.n3)
Asap.dof_signature(::_TriBars) = ntuple(i -> mod1(i, 6) <= 3, 18)

# 6×18 0/1 gather of the translational+rotational slots of nodes a, b
_tri_gather(a::Int, b::Int) = SMatrix{12,18,Float64}(
    [(r <= 6 ? (c == 6 * (a - 1) + r) : (c == 6 * (b - 1) + r - 6)) ? 1.0 : 0.0 for r in 1:12, c in 1:18])
const _TRI_P = (_tri_gather(1, 2), _tri_gather(2, 3), _tri_gather(3, 1))

function _bar12(section, xa, xb)
    k = Asap.truss_stiffness(section, xa, xb)          # 6×6 (translations)
    Z = zero(SMatrix{3,3,eltype(k)})
    # embed translations into a 12×12 two-node block (rotations zero)
    return vcat(hcat(k[1:3, 1:3], Z, k[1:3, 4:6], Z), hcat(Z, Z, Z, Z),
                hcat(k[4:6, 1:3], Z, k[4:6, 4:6], Z), hcat(Z, Z, Z, Z))
end

Asap.stiffness(el::_TriBars, x1, x2, x3) =
    _TRI_P[1]' * _bar12(el.section, x1, x2) * _TRI_P[1] +
    _TRI_P[2]' * _bar12(el.section, x2, x3) * _TRI_P[2] +
    _TRI_P[3]' * _bar12(el.section, x3, x1) * _TRI_P[3]

function _nnode_models()
    mat = Material(200.0, 77.0, 1.0, 0.3)
    sec = Section(mat, 10.0)
    mk() = [Node([0.0, 0.0, 0.0], :pinned), Node([4.0, 0.0, 0.0], :pinned),
            Node([0.0, 4.0, 0.5], :pinned), Node([2.0, 3.0, 0.2], :free),
            Node([2.0, 1.5, 3.0], :free)]
    A, B, E, C, D = mk()
    tri = Model([A, B, E, C, D],
        AbstractElement{Float64}[_TriBars(A, B, C, sec), _TriBars(C, D, E, sec),
                                 TrussElement(B, D, sec), TrussElement(A, D, sec)],
        AbstractLoad{Float64}[NodeForce(C, [1.0, -2.0, -3.0]), NodeForce(D, [0.5, 1.0, -4.0])])
    A, B, E, C, D = mk()
    bars = [(A, B), (B, C), (C, A), (C, D), (D, E), (E, C), (B, D), (A, D)]
    ref = Model([A, B, E, C, D],
        AbstractElement{Float64}[TrussElement(a, b, sec) for (a, b) in bars],
        AbstractLoad{Float64}[NodeForce(C, [1.0, -2.0, -3.0]), NodeForce(D, [0.5, 1.0, -4.0])])
    return tri, ref
end

@testset "N-node element path (S1)" begin
    tri, ref = _nnode_models()
    @test FrameElement{Float64,Section{Float64}} <: AbstractLineElement
    @test TrussElement{Float64,Section{Float64}} <: AbstractLineElement
    @test VariableElement{Float64} <: AbstractLineElement
    @test !(_TriBars{Float64,Section{Float64}} <: AbstractLineElement)

    process!(tri)
    process!(ref)
    el = tri.elements[1]
    @test Asap.element_global_dofs(el) ==
          vcat((6 * (n.index - 1) .+ (1:6) for n in Asap.nodes(el))...)
    g = only(filter(g -> g isa Asap.ElementGroup{<:_TriBars}, tri.cache.groups))
    @test size(g.inodes) == (3, 2)
    @test g.inodes[:, 1] == [n.index for n in Asap.nodes(tri.elements[1])]
    # rotations retired structurally (translation-only signature)
    @test length(tri.cache.partition.free) == length(ref.cache.partition.free) == 6
    @test length(tri.cache.q_local[1][1]) == 18

    rt = solve!(tri).results
    rr = solve!(ref).results
    @test rt.u ≈ rr.u rtol = 1e-12
    @test rt.reactions ≈ rr.reactions rtol = 1e-10 atol = 1e-12
    @test rt.compliance ≈ rr.compliance rtol = 1e-12
    # generic recovery: global nodal forces, 6 per node
    @test length(element_forces(rt, tri.elements[1])) == 18

    # surface-free utilities refuse non-line models instead of misreading fields
    @test_throws ArgumentError connectivity(tri)
    @test_throws MethodError endpoints(tri.elements[1])

    # pure path: parity + gradient through the generic N-node closure
    st, sr = extract_state(tri), extract_state(ref)
    @test compliance(tri, st) ≈ rt.compliance rtol = 1e-12
    gt = Zygote.gradient(X -> compliance(tri, ModelState{Float64}(X, st.sections)), st.X)[1]
    gr = Zygote.gradient(X -> compliance(ref, ModelState{Float64}(X, sr.sections)), sr.X)[1]
    @test gt ≈ gr rtol = 1e-9
end
