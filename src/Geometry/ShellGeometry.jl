"""
    ShellGeo <: AbstractGeo

Plot-ready arrays of the SHELL elements of a solved model (the counterpart of
[`ModelGeo`](@ref), which covers frames and trusses). Matrices are
column-per-item, ready for Makie's `mesh`/`poly`/`arrows`.

# Fields
- `nodes::Matrix{Float64}`: 3 × n_nodes positions (model node order; all
  nodes, so `quads` index directly)
- `disp::Matrix{Float64}`: 3 × n_nodes translational displacements
- `quads::Matrix{Int}`: 4 × n_shells node indices (counter-clockwise); a
  triangle repeats its third node in row 4 (a degenerate quad, which polygon
  plotting draws as the triangle). [`shell_faces`](@ref) gives triangles.
- `elements::Vector{Int}`: model element index of each shell
- `centers`, `normals`: 3 × n_shells element centroids and unit normals
- `frames::Vector{SMatrix{3,3,Float64,9}}`: local axes (rows e₁, e₂, e₃)
- `N`, `M`: 3 × n_shells centre resultants `(xx, yy, xy)` in local axes —
  membrane forces [force/length] and moments [force·length/length]
- `Q`: 2 × n_shells centre transverse shears [force/length]
- `N1`, `N2`: principal membrane forces (N1 ≥ N2), `θ`: angle of N1 from e₁ [rad]
- `thickness::Vector{Float64}`: shell thickness (NaN for rigidity sections)

# Constructor
    ShellGeo(model; results = model.results, factors = nothing)

`results`/`factors` select a case or combination as in
[`shell_resultants`](@ref). Use [`nodal_average`](@ref) to turn any
per-element field into per-node values for smooth contour plots, and
`nodes .+ s .* disp` for a deformed shape scaled by `s`.
"""
struct ShellGeo <: AbstractGeo
    nodes::Matrix{Float64}
    disp::Matrix{Float64}
    quads::Matrix{Int}
    elements::Vector{Int}
    centers::Matrix{Float64}
    normals::Matrix{Float64}
    frames::Vector{SMatrix{3,3,Float64,9}}
    N::Matrix{Float64}
    M::Matrix{Float64}
    Q::Matrix{Float64}
    N1::Vector{Float64}
    N2::Vector{Float64}
    θ::Vector{Float64}
    thickness::Vector{Float64}
end

function ShellGeo(model::Model{Float64}; results = model.results, factors = nothing)
    res = results::LinearResults{Float64}
    nodes_ = reduce(hcat, [Vector(n.position) for n in model.nodes])
    disp = reduce(hcat, [Vector(displacement(res, n)[1:3]) for n in model.nodes])
    idx = [e for (e, el) in enumerate(model.elements) if el isa ShellElement]
    ns = length(idx)
    quads = zeros(Int, 4, ns)
    centers = zeros(3, ns); normals = zeros(3, ns)
    frames = Vector{SMatrix{3,3,Float64,9}}(undef, ns)
    N = zeros(3, ns); M = zeros(3, ns); Q = zeros(2, ns)
    N1 = zeros(ns); N2 = zeros(ns); θ = zeros(ns); thick = fill(NaN, ns)
    for (k, e) in enumerate(idx)
        el = model.elements[e]
        ids = [n.index for n in el.nodes]
        quads[:, k] = length(ids) == 4 ? ids : [ids; ids[3]]
        centers[:, k] = sum(n -> n.position, el.nodes) / length(el.nodes)
        r = shell_resultants(model, el; results = res, factors = factors)
        frames[k] = r.frame
        normals[:, k] = r.frame[3, :]
        N[:, k] = r.N; M[:, k] = r.M; Q[:, k] = r.Q
        N1[k], N2[k], θ[k] = principal_forces(r.N)
        el.section isa ShellSection && (thick[k] = thickness(el.section))
    end
    return ShellGeo(nodes_, disp, quads, idx, centers, normals, frames, N, M, Q, N1, N2, θ, thick)
end

"""
    nodal_average(g::ShellGeo, values::AbstractVector; weight = :count) -> Vector{Float64}

Average a per-shell-element field (length `n_shells`, e.g. `g.N[1, :]`) to
the nodes (length `n_nodes`; 0 at nodes no shell touches) — the usual
smoothing for contour plots. Explicit and opt-in: averaging across
discontinuities (wall/slab junctions, thickness changes) blurs real jumps.

`weight` sets each element's share at a node:
- `:count`: equal shares (a plain mean of the elements meeting there)
- `:area`: proportional to the element's area [length²], so the many small
  elements of a graded region don't outvote the large ones beside it
- `:angle`: proportional to the element's interior angle at the node [rad]
  (independent of element size; on a regular mesh close to `:area`)

Average Cartesian components (`g.N[1, :]`, `g.N[2, :]`, `g.N[3, :]`, in a
common frame) and derive principal forces and directions from the averaged
components ([`principal_forces`](@ref)); never average `N1`, `N2` or `θ`
themselves (principal values are not linear in the field, and angles wrap).
"""
function nodal_average(g::ShellGeo, values::AbstractVector; weight::Symbol = :count)
    length(values) == size(g.quads, 2) || throw(DimensionMismatch("one value per shell element"))
    weight in (:count, :area, :angle) ||
        throw(ArgumentError("weight must be :count, :area or :angle (got :$weight)"))
    acc = zeros(size(g.nodes, 2)); wsum = zeros(size(g.nodes, 2))
    pos(i) = SVector{3}(g.nodes[1, i], g.nodes[2, i], g.nodes[3, i])
    for k in axes(g.quads, 2)
        n = g.quads[4, k] == g.quads[3, k] ? 3 : 4           # triangle: row 4 repeats node 3
        ids = ntuple(a -> g.quads[a, k], 4)
        A = weight === :area ? _polygon_area(map(pos, ids[1:n])) : 1.0
        for a in 1:n
            w = weight === :angle ? _corner_angle(pos(ids[mod1(a - 1, n)]), pos(ids[a]), pos(ids[mod1(a + 1, n)])) : A
            acc[ids[a]] += w * values[k]
            wsum[ids[a]] += w
        end
    end
    return [w > 0 ? s / w : 0.0 for (s, w) in zip(acc, wsum)]
end

# area of a flat (or nearly flat) polygon: half the norm of Σ xₖ × xₖ₊₁
_polygon_area(X) = norm(sum(k -> cross(X[k], X[mod1(k + 1, length(X))]), eachindex(X))) / 2

# interior angle at b of the corner a–b–c [rad]
_corner_angle(a, b, c) = atan(norm(cross(a - b, c - b)), dot(a - b, c - b))

"""
    shell_faces(g::ShellGeo) -> Matrix{Int}

All shell faces as triangles (3 × n_faces node indices, counter-clockwise):
each quadrilateral split along its 1–3 diagonal, each triangle as is. The
element of face `f` is `g.elements[shell_face_owner(g)[f]]`. Ready for
Makie's `mesh(g.nodes', shell_faces(g)')`.
"""
function shell_faces(g::ShellGeo)
    F = Int[]
    for k in axes(g.quads, 2)
        q = view(g.quads, :, k)
        append!(F, (q[1], q[2], q[3]))
        q[4] == q[3] || append!(F, (q[1], q[3], q[4]))
    end
    return reshape(F, 3, :)
end

"""
    shell_face_owner(g::ShellGeo) -> Vector{Int}

For each face of [`shell_faces`](@ref), the column `k` of the ShellGeo arrays
(element-wise fields such as `g.N1[k]`) it belongs to.
"""
function shell_face_owner(g::ShellGeo)
    o = Int[]
    for k in axes(g.quads, 2)
        push!(o, k)
        g.quads[4, k] == g.quads[3, k] || push!(o, k)
    end
    return o
end

"""
    principal_forces(N) -> (N1, N2, θ)

Principal membrane forces of `N = (Nxx, Nyy, Nxy)` [force/length]:
`N1 ≥ N2` and the angle `θ` [rad] of the N1 direction from the local e₁
axis. Compression trajectories (`N2 < 0`) are the strut directions of a
strut-and-tie reading of a wall.
"""
function principal_forces(N::AbstractVector)
    c = (N[1] + N[2]) / 2
    r = sqrt(((N[1] - N[2]) / 2)^2 + N[3]^2)
    return c + r, c - r, atan(2 * N[3], N[1] - N[2]) / 2
end
