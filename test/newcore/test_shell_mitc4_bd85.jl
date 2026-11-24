# MITC4 against the primary source: Bathe & Dvorkin (1985), "A four-node plate
# bending element based on Mindlin/Reissner plate theory and a mixed
# interpolation", IJNME 21:367–383 (read from the scanned PDF hosted by Bathe's
# group at MIT). Checks: (1) the transverse-shear field of eqs. (8)–(10),
# implemented verbatim in the paper's node numbering, against Asap's
# operator; (2) the Fig. 3 patch (constant bending, shear, twist fields,
# thick and thin); (3) the Fig. 5 cantilever under tip moment — one element
# and two distorted elements reproduce Bernoulli beam theory exactly;
# (4) the Fig. 8 thin-plate (L/h = 1000) accuracy trend.

using StaticArrays, LinearAlgebra
include("shell_test_utils.jl")

@testset "MITC4 vs Bathe & Dvorkin 1985" begin

    @testset "eqs. (8)–(10): shear strains, verbatim" begin
        # paper nodes 1(1,1) 2(−1,1) 3(−1,−1) 4(1,−1) = Asap nodes 3, 4, 1, 2
        xyA = (SVector(0.0, 0.0), SVector(2.3, 0.4), SVector(2.0, 2.6), SVector(-0.3, 1.8))   # Asap order
        P = (xyA[3], xyA[4], xyA[1], xyA[2])                                                   # paper order
        x = getindex.(P, 1); y = getindex.(P, 2)
        Ax = x[1] - x[2] - x[3] + x[4]; Bx = x[1] - x[2] + x[3] - x[4]; Cx = x[1] + x[2] - x[3] - x[4]
        Ay = y[1] - y[2] - y[3] + y[4]; By = y[1] - y[2] + y[3] - y[4]; Cy = y[1] + y[2] - y[3] - y[4]
        dofs = [0.31, -0.7, 0.45, 1.2, 0.05, -0.33, -0.6, 0.9, 0.12, 0.27, -0.44, 0.8]   # (w, θx, θy) × paper node
        w(i) = dofs[3i-2]; θx(i) = dofs[3i-1]; θy(i) = dofs[3i]
        term(i, j) = (w(i) - w(j)) / 2 + (x[i] - x[j]) / 4 * (θy(i) + θy(j)) - (y[i] - y[j]) / 4 * (θx(i) + θx(j))
        function paper_gamma(r, s)
            xr, yr = (Ax + s * Bx) / 4, (Ay + s * By) / 4
            xs, ys = (Cx + r * Bx) / 4, (Cy + r * By) / 4
            detJ = xr * ys - yr * xs
            γrz = sqrt((Cx + r * Bx)^2 + (Cy + r * By)^2) / (8detJ) * ((1 + s) * term(1, 2) + (1 - s) * term(4, 3))   # (9a)
            γsz = sqrt((Ax + s * Bx)^2 + (Ay + s * By)^2) / (8detJ) * ((1 + r) * term(1, 4) + (1 - r) * term(2, 3))   # (9b)
            α, β = atan(yr, xr), atan(ys, xs)                                       # r-, s-axis angles from x
            return SVector(γrz * sin(β) - γsz * sin(α), -γrz * cos(β) + γsz * cos(α))   # (8)
        end
        # the same DOFs in Asap's node order
        perm = (3, 4, 1, 2)   # Asap node a ↦ paper node perm[a]
        ap = SVector{12}(vcat([dofs[3perm[a]-2:3perm[a]] for a in 1:4]...))
        for (r, s) in ((0.0, 0.0), (0.37, -0.81), (-0.9, 0.4), (1 / sqrt(3), -1 / sqrt(3)))
            @test Asap._mitc_gamma(xyA, ap, r, s) ≈ paper_gamma(r, s) rtol = 1e-12
        end
    end

    @testset "Fig. 3 patch: constant bending, shear and twist (thick and thin)" begin
        # the paper's patch: 10 × 10 square, interior nodes (2,2), (8,3), (8,7), (4,7)
        P = [SVector(0.0, 0.0), SVector(10.0, 0.0), SVector(10.0, 10.0), SVector(0.0, 10.0),
             SVector(2.0, 2.0), SVector(8.0, 3.0), SVector(8.0, 7.0), SVector(4.0, 7.0)]
        quads = [(1, 2, 6, 5), (2, 3, 7, 6), (3, 4, 8, 7), (4, 1, 5, 8), (5, 6, 7, 8)]
        fields = (
            bending = p -> [(p[1]^2 - 0.4 * p[2]^2) / 2, -0.4 * p[2], -p[1]],   # Kirchhoff, κ constant
            twist = p -> [p[1] * p[2], p[1], -p[2]],                           # w = xy (θx = w,y, θy = −w,x)
            shear = p -> [0.3 * p[1] - 0.2 * p[2], 0.0, 0.0])                   # constant γ = (0.3, −0.2)
        for t in (1.0, 0.001), (name, field) in pairs(fields)
            sec = ShellSection(Material(2.1e6, 1.0, 0.3), t)
            K = zeros(24, 24)
            for q in quads
                Ke = Asap.mitc4_plate_stiffness(Tuple(P[i] for i in q), bending_rigidity(sec), shear_rigidity(sec))
                d = vcat((3 * (i - 1) .+ (1:3) for i in q)...)
                K[d, d] .+= Ke
            end
            U = vcat(field.(P)...)
            # the paper's SHEARING case restrains θx = θy = 0 at ALL nodes (a constant
            # transverse shear is not a Mindlin equilibrium state otherwise: constant
            # Q needs a moment gradient) — only the interior w are unknown there
            inner = name === :shear ? [13, 16, 19, 22] : collect(13:24)
            known = setdiff(1:24, inner)
            Ui = -K[inner, inner] \ (K[inner, known] * U[known])
            @test Ui ≈ U[inner] rtol = 1e-8 atol = 1e-12
        end
    end

    @testset "Fig. 5: cantilever, tip moment — Bernoulli beam theory exactly" begin
        # E = 2.1e6, ν = 0.3, L = 100, width 1; plate DOFs only
        # the paper's t = 0.1 (L/t = 1000) is exact up to round-off of a stiffness
        # matrix with condition ≈ 1e10 (shear vs bending of one 100 × 0.1 element):
        # checked at 1e-5 there and at 1e-7 for t = 1 (4e-11 for one element, 4e-9 distorted)
        E, ν, L, b, M = 2.1e6, 0.3, 100.0, 1.0, 1.0
        for (t, tol) in ((1.0, 1e-7), (0.1, 1e-5))
        sec = ShellSection(Material(E, 1.0, ν), t)
        EI = E * b * t^3 / 12            # free edges: anticlastic curvature ⇒ beam theory with E
        function cantilever(xs_bottom, xs_top)
            # plate DOFs (w, θx, θy) free; root: w = 0 and bending rotation θy = 0
            fix(x) = x == 0 ? Bool[0, 0, 0, 1, 0, 0] : Bool[0, 0, 1, 1, 1, 0]
            nb = [Node([xb, 0.0, 0.0], fix(xb)) for xb in xs_bottom]
            nt = [Node([xt, b, 0.0], fix(xt)) for xt in xs_top]
            els = [ShellElement((nb[k], nb[k+1], nt[k+1], nt[k]), sec) for k in 1:length(nb)-1]
            loads = AbstractLoad{Float64}[NodeMoment(nb[end], [0.0, -M / 2, 0.0]), NodeMoment(nt[end], [0.0, -M / 2, 0.0])]
            m = Model(vcat(nb, nt), AbstractElement{Float64}[els...], loads)
            solve!(m)
            return displacement(m.results, nb[end])[3], displacement(m.results, nb[end])[5]
        end
        for (xb, xt) in (([0.0, L], [0.0, L]), ([0.0, 0.4L, L], [0.0, 0.6L, L]))   # one element; two distorted (±0.1L)
            wtip, θtip = cantilever(xb, xt)
            @test wtip ≈ M * L^2 / (2EI) rtol = tol
            @test abs(θtip) ≈ M * L / EI rtol = tol
        end
        end
    end

    @testset "Fig. 8: thin plates (L/h = 1000), quarter-plate meshes N = 2, 4, 6" begin
        # the paper's errors in centre deflection are a few percent at N = 2 and
        # vanish by N = 6 (read from its plots); a full n × n mesh ≡ quarter N = n/2
        errs = [square_plate(2N; t = 1e-3) / 0.00406 - 1 for N in (2, 4, 6)]
        @test abs(errs[1]) < 0.03
        @test abs(errs[3]) < 0.005
        @test abs(errs[3]) < abs(errs[2]) < abs(errs[1])
    end
end
