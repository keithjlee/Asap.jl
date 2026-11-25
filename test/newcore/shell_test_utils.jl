# included by every shell test file (each also runs standalone); load once
if !isdefined(@__MODULE__, :_SHELL_TEST_UTILS_LOADED)
const _SHELL_TEST_UTILS_LOADED = true

# Shared helpers for the shell test files (membrane-only planar problems).

# structured membrane mesh in the XY plane: `xy(i, j)` gives node (i, j) for
# i ∈ 0:nx, j ∈ 0:ny; elements counter-clockwise. Out-of-plane DOFs
# (uz, θx, θy) restrained everywhere — membrane-only problems. `fix(i, j)`
# returns `nothing` or a (u, v, ψ) restraint triple.
function membrane_mesh(nx, ny, xy, sec; fix = (i, j) -> nothing)
    grid = Matrix{Node{Float64}}(undef, nx + 1, ny + 1)
    for i in 0:nx, j in 0:ny
        f = fix(i, j)
        free = f === nothing ? Bool[1, 1, 0, 0, 0, 1] : Bool[!f[1], !f[2], 0, 0, 0, !f[3]]
        grid[i+1, j+1] = Node(collect(Float64, (xy(i, j)..., 0.0)), free)
    end
    els = AbstractElement{Float64}[ShellElement((grid[i, j], grid[i+1, j], grid[i+1, j+1], grid[i, j+1]), sec)
                                   for i in 1:nx for j in 1:ny]
    return grid, els
end

# consistent nodal loads of an edge traction `q(y)` (quadratic at most) on the
# vertical edge whose node heights are `ys` (linear shape functions; Simpson
# is exact for quadratic × linear)
function _edge_traction_loads(ys, q)
    F = zeros(length(ys))
    for j in 1:length(ys)-1
        ya, yb = ys[j], ys[j+1]
        Na(y) = (yb - y) / (yb - ya)
        simpson(f) = (yb - ya) / 6 * (f(ya) + 4f((ya + yb) / 2) + f(yb))
        F[j] += simpson(y -> q(y) * Na(y))
        F[j+1] += simpson(y -> q(y) * (1 - Na(y)))
    end
    return F
end

# cantilever of length L, depth H, unit thickness (E = 30000, ν = 0.25),
# FULLY clamped at x = 0 (u = v = ψ = 0 at every root node: with a drilling
# element, leaving ψ free lets the "clamped" edge bulge quadratically between
# nodes), parabolic tip shear P = 40. Returns compliance / P (the
# load-weighted mean tip deflection) [length].
function cantilever_wall(nx, ny; L = 48.0, H = 12.0)
    sec = ShellSection(Material(30000.0, 1.0, 0.25), 1.0)
    grid, els = membrane_mesh(nx, ny, (i, j) -> (L * i / nx, H * j / ny), sec;
        fix = (i, j) -> i == 0 ? (true, true, true) : nothing)
    F = _edge_traction_loads([H * j / ny for j in 0:ny], y -> 40 * 1.5 / H * (1 - (2y / H - 1)^2))
    m = Model(vec(grid), els, AbstractLoad{Float64}[NodeForce(grid[nx+1, j+1], [0.0, F[j+1], 0.0]) for j in 0:ny])
    solve!(m)
    return sum(F[j+1] * displacement(m.results, grid[nx+1, j+1])[2] for j in 0:ny) / 40
end

# Cook's membrane (E = 1, ν = 1/3, t = 1; corners (0,0), (48,44), (48,60),
# (0,44)), clamped left edge, uniform shear P = 1 on the right edge; vertical
# displacement at the right-edge MIDPOINT (48, 52) — the point the 23.96
# reference refers to (the top corner converges to ≈ 25.1)
function cook_mid(n)
    sec = ShellSection(Material(1.0, 1.0, 1 / 3), 1.0)
    grid, els = membrane_mesh(n, n, (i, j) -> (48i / n, 44i / n + (44 - 28i / n) * j / n), sec;
        fix = (i, j) -> i == 0 ? (true, true, true) : nothing)
    loads = [NodeForce(grid[n+1, j+1], [0.0, (j == 0 || j == n ? 0.5 : 1.0) / n, 0.0]) for j in 0:n]
    m = Model(vec(grid), els, AbstractLoad{Float64}[loads...])
    solve!(m)
    return displacement(m.results, grid[n+1, n÷2+1])[2]
end

# general structured shell mesh: `X(i, j)` gives the 3D position of node (i, j);
# `fix(i, j)` returns a 6-Bool fixity vector (true = free) or `nothing` (all free)
function grid_shell_mesh(nx, ny, X, sec; fix = (i, j) -> nothing)
    grid = Matrix{Node{Float64}}(undef, nx + 1, ny + 1)
    for i in 0:nx, j in 0:ny
        f = fix(i, j)
        grid[i+1, j+1] = Node(collect(Float64, X(i, j)), f === nothing ? trues(6) |> collect : f)
    end
    els = AbstractElement{Float64}[ShellElement((grid[i, j], grid[i+1, j], grid[i+1, j+1], grid[i, j+1]), sec)
                                   for i in 1:nx for j in 1:ny]
    return grid, els
end

# uniform area load q (global vector) lumped to nodes: q·area/4 per corner — for
# bilinear translations on parallelogram elements this IS the consistent load
function lumped_area_loads(els, q::AbstractVector)
    acc = Dict{Node{Float64},Vector{Float64}}()
    for el in els, n in el.nodes
        acc[n] = get(acc, n, zeros(3)) .+ q .* (Asap.area(el) / 4)
    end
    return AbstractLoad{Float64}[NodeForce(n, f) for (n, f) in acc]
end

# square plate a × a in the XY plane, n × n mesh, pure plate problem (in-plane
# DOFs and drilling restrained); `edge` = :ss (w = 0, rotations free) or
# :clamped. Uniform pressure q (or a centre point load P when `point = true`).
# Returns the centre deflection normalized: w·D/(q a⁴) or w·D/(P a²).
function square_plate(n; a = 1.0, t = 0.01, E = 1.0e4, ν = 0.3, edge = :ss, point = false)
    sec = ShellSection(Material(E, 1.0, ν), t)
    onedge(i, j) = i == 0 || j == 0 || i == n || j == n
    fix(i, j) = onedge(i, j) ? (edge == :ss ? Bool[0, 0, 0, 1, 1, 0] : Bool[0, 0, 0, 0, 0, 0]) :
                Bool[0, 0, 1, 1, 1, 0]
    grid, els = grid_shell_mesh(n, n, (i, j) -> (a * i / n, a * j / n, 0.0), sec; fix = fix)
    loads = point ? AbstractLoad{Float64}[NodeForce(grid[n÷2+1, n÷2+1], [0.0, 0.0, -1.0])] :
            lumped_area_loads(els, [0.0, 0.0, -1.0])
    m = Model(vec(grid), els, loads)
    solve!(m)
    D = E * t^3 / (12 * (1 - ν^2))
    w = -displacement(m.results, grid[n÷2+1, n÷2+1])[3]
    return point ? w * D / a^2 : w * D / a^4
end

# Scordelis–Lo roof (MacNeal & Harder 1985): cylinder R = 25, L = 50, 40° half
# angle, t = 0.25, E = 4.32e8, ν = 0, gravity 90 per unit area; rigid end
# diaphragms. Quarter model, n × n mesh; returns the vertical deflection at the
# midspan of the free edge (reference 0.3024).
function scordelis_lo(n)
    R, L, Φ = 25.0, 50.0, deg2rad(40.0)
    sec = ShellSection(Material(4.32e8, 1.0, 0.0), 0.25)
    X(i, j) = (L / 2 * i / n, R * sin(Φ * j / n), R * cos(Φ * j / n))
    function fix(i, j)
        f = trues(6)
        i == 0 && (f[2] = f[3] = false)                       # diaphragm: uy = uz = 0
        i == n && (f[1] = f[5] = f[6] = false)                # midspan symmetry: ux = θy = θz = 0
        j == 0 && (f[2] = f[4] = f[6] = false)                # crown symmetry: uy = θx = θz = 0
        return collect(f)
    end
    grid, els = grid_shell_mesh(n, n, X, sec; fix = fix)
    m = Model(vec(grid), els, lumped_area_loads(els, [0.0, 0.0, -90.0]))
    solve!(m)
    return -displacement(m.results, grid[n+1, n+1])[3]
end

end # load-once guard
