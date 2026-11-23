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
function edge_loads(ys, q)
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
    F = edge_loads([H * j / ny for j in 0:ny], y -> 40 * 1.5 / H * (1 - (2y / H - 1)^2))
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
