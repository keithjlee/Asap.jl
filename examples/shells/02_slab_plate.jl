#=
Walkthrough 2 — a slab, standalone (plate bending)

Shift+enter through the blocks. Units: kN, m.

The same ShellElement carries out-of-plane load as a Reissner–Mindlin plate
(MITC4 — no shear locking, thin or thick). Two cases:

A. a simply supported square slab under uniform pressure, checked against
   Timoshenko & Woinowsky-Krieger (centre deflection 0.00406·qa⁴/D, centre
   moment 0.0479·qa² at ν = 0.3, rescaled for our ν). Those are THIN-plate (Kirchhoff) values: a real
   200 mm slab over 6 m (a/t = 30) is a little more flexible, because a
   Reissner–Mindlin plate also deflects in transverse shear;
B. the same slab carried only by point supports (four corner columns and a
   central one) — a flat plate, where moments concentrate in column strips.
=#

begin
    using Asap, LinearAlgebra
    include(joinpath(@__DIR__, "..", "figures.jl")); set_figdir(@__DIR__)
end

#=
## A — simply supported slab vs closed form

`shell_mesh`-style structured slab built directly with `ShellElement`s; the
edge nodes get w = 0 with rotations free (soft simple support), and the
in-plane DOFs are restrained so this is a pure plate problem.
=#
begin
    a, t, q = 6.0, 0.20, 5.0                       # span [m], thickness [m], pressure [kPa]
    concrete = Material(30.0e6, 2.4, 0.2)

    function square_slab(n; support = :edges, t = t)
        sec = ShellSection(concrete, t)
        grid = Matrix{Node{Float64}}(undef, n + 1, n + 1)
        for i in 0:n, j in 0:n
            onedge = i == 0 || j == 0 || i == n || j == n
            corner_or_mid = (i in (0, n) && j in (0, n)) || (i == n ÷ 2 && j == n ÷ 2)
            supported = support === :edges ? onedge : corner_or_mid
            # free: in-plane DOFs + drilling restrained (pure plate); w restrained where supported
            fix = Bool[0, 0, !supported, 1, 1, 0]
            grid[i+1, j+1] = Node([a * i / n, a * j / n, 0.0], fix)
        end
        els = [ShellElement((grid[i, j], grid[i+1, j], grid[i+1, j+1], grid[i, j+1]), sec)
               for i in 1:n for j in 1:n]
        loads = AbstractLoad{Float64}[AreaLoad(el, [0.0, 0.0, -q]) for el in els]
        m = Model(vec(grid), AbstractElement{Float64}[els...], loads)
        solve!(m)
        return m, grid, els
    end

    Dt(t) = 30.0e6 * t^3 / (12 * (1 - 0.2^2))
    # thin slab (a/t = 300): the Kirchhoff closed form applies
    tthin = a / 300
    mA, gridA, elsA = square_slab(16; t = tthin)
    centre = gridA[9, 9]
    wc = -displacement(mA.results, centre)[3]
    # moments: sample at an element CENTRE (superconvergent) — with an odd
    # mesh, the central element's centre is the slab centre
    mO, gridO, elsO = square_slab(17; t = tthin)
    elc = only(filter(e -> norm(sum(n -> n.position, e.nodes) / 4 - [a / 2, a / 2, 0]) < 1e-9, elsO))
    Mc = shell_resultants(mO, elc).M
    # Timoshenko's table is for ν = 0.3. By symmetry, the centre moment is
    # M = (1 + ν)·D·|w,xx|, and D·w does not depend on ν (SS edges), so for
    # this ν = 0.2 concrete the reference is 0.0479·1.2/1.3. The deflection
    # coefficient 0.00406 holds for any ν.
    Mref = 0.0479 * (1 + 0.2) / (1 + 0.3)
    checkA = (w_ratio = wc * Dt(tthin) / (q * a^4) / 0.00406, M_ratio = abs(Mc[1]) / (q * a^2) / Mref)
    @show checkA
end

#=
Both ratios ≈ 1. Deflection and bending-moment fields (plan view):
=#
begin
    mA, gridA, elsA = square_slab(16)          # the real 200 mm slab for the field plots
    f = fig_slab(mA; field = :w, title = "SS slab: deflection")
    savefig("02_slab_ss_w.pdf", f)
end

begin
    f = fig_slab(mA; field = :Mxx, title = rich("SS SLAB: M", subscript("xx")))
    savefig("02_slab_ss_Mxx.pdf", f)
end

#=
## Convergence of the closed-form coefficient

The thin slab converges to the Kirchhoff value; the real a/t = 30 slab
settles ≈ 3% above it — the transverse-shear flexibility a thin-plate
element would miss (and MITC4 captures without locking in the thin limit).
=#
begin
    ns = [4, 8, 16, 24]
    coeff(n, t) = (r = square_slab(n; t = t); -displacement(r[1].results, r[2][n÷2+1, n÷2+1])[3] * Dt(t) / (q * a^4))
    thin = [coeff(n, a / 300) for n in ns]
    real_ = [coeff(n, t) for n in ns]
    f = fig_series(ns, [("thin, {a/t} = 300", thin), ("200 {mm} slab, {a/t} = 30", real_)];
        xlabel = "elements per side", ylabel = "{w·D/(q a⁴)}", ref = 0.00406,
        reflabel = "Kirchhoff 0.00406", title = "centre deflection", legendpos = :rb)
    savefig("02_slab_convergence.pdf", f)
    @show thin real_
end

#=
## B — flat plate on point supports

Supports only at the four corners and the centre: the slab now spans to
points, and the hogging moment concentrates over the central column.
=#
begin
    mB, gridB, elsB = square_slab(24; support = :points)
    gB = ShellGeo(mB)
    flatplate = (max_deflection = maximum(gB.disp[3, :] .|> abs), Mxx_range = extrema(gB.M[1, :]))
    @show flatplate
end

begin
    f = fig_slab(mB; geo = gB, field = :w, title = "flat plate: deflection")
    savefig("02_flatplate_w.pdf", f)
end

begin
    f = fig_slab(mB; geo = gB, field = :Mxx, title = rich("FLAT PLATE: M", subscript("xx")))
    savefig("02_flatplate_Mxx.pdf", f)
end

#=
## Numbers for the report

`reports/S2_slab_plate.pdf` prints these (see `figures.jl`, `report_data`).
=#
begin
    report_data("S2";
        span = a, thickness = t, pressure = q, E = concrete.E, nu = concrete.ν,
        thin_thickness = tthin, slenderness_thin = a / tthin, slenderness_real = a / t,
        D_thin = Dt(tthin), D_real = Dt(t),
        n_w = 16, n_M = 17,
        w_centre_thin_mm = 1000 * wc, w_coeff = wc * Dt(tthin) / (q * a^4), w_coeff_ref = 0.00406,
        M_centre = Mc[1], M_coeff = abs(Mc[1]) / (q * a^2), M_coeff_ref_nu03 = 0.0479, M_coeff_ref = Mref,
        checkA = checkA,
        w_centre_real_mm = -1000 * displacement(mA.results, gridA[9, 9])[3],
        convergence = (ns = ns, thin = thin, real = real_,
                       thin_ratio = thin ./ 0.00406, real_ratio = real_ ./ 0.00406),
        flatplate = (n = 24, max_deflection_mm = 1000 * flatplate.max_deflection,
                     Mxx_min = flatplate.Mxx_range[1], Mxx_max = flatplate.Mxx_range[2]))
end
