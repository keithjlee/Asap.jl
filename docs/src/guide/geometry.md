# Plot-ready geometry extraction

```@setup geo
using Asap
sec  = Section(Steel_kNm, 1e-3, 1e-6, 1e-6, 1e-6)
fsec = Section(Steel_kNm, 1e-2, 1e-4, 5e-5, 1e-6)
truss = Warren2D(11, 1.5, 2.0, sec; load = [0.0, -20.0, 0.0])
bldg  = Frame(2, 6.0, 2, 5.0, 2, 4.0, 2.0, fsec, fsec, fsec, fsec)
```

[`Geo`](@ref)`(model)` (or `Geo(network)`) flattens a solved structure into plain arrays for plotting — node positions, displaced positions, element connectivity, and per-element force ranges with their maxima. [`ElementDisplacements`](@ref)`(element, model)` samples the exact deflected curve along a member.

The examples reuse the `truss` (a solved `Warren2D`) and `bldg` (a solved `Frame`) from [Structure generators](generators.md):

```@example geo
geo = Geo(truss.model)          # geo.nodes, geo.disp, geo.indices, geo.P, …
```

```@example geo
ed = ElementDisplacements(bldg.model.elements[1], bldg.model; resolution = 20)
ed.basepositions .+ 100 .* ed.uglobal   # 3×20 displaced curve, scaled ×100
```
