# Parametric structure generators

```@setup generators
using Asap
```

A family of one-line generators (absorbed from AsapToolkit) builds, loads, and solves common structural typologies — handy for testing, benchmarking, and optimization studies:

```@example generators
sec = Section(Steel_kNm, 1e-3, 1e-6, 1e-6, 1e-6)

truss = Warren2D(11, 1.5, 2.0, sec; load = [0.0, -20.0, 0.0])   # 2D Warren arch
grid  = SpaceFrame(6, 1.2, 6, 1.2, 1.0, sec; support = :corner) # double-layer grid

# variable-depth spaceframe: any callable surface(u, v) on [0,1]² sets the top layer
vault = SpaceFrame(6, 1.2, 6, 1.2, 1.0, (u, v) -> 0.5 * sinpi(u) * sinpi(v), sec)

# every generator returns a solved model plus its generation parameters
maximum(abs, truss.model.results.u)
```

```@example generators
# 3D building frame: columns, primary beams, joists (with releases), braces
fsec = Section(Steel_kNm, 1e-2, 1e-4, 5e-5, 1e-6)
bldg = Frame(2, 6.0, 2, 5.0, 2, 4.0, 2.0, fsec, fsec, fsec, fsec)
```

```@example generators
# ground structures for layout optimization: dense candidate grids → models
gs = XGroundStructure(6.0, 4, 4.0, 3)
candidates = to_truss(gs, sec; load = [1.0, 0.0, 0.0])
```

The full set: [`Warren2D`](@ref), [`Pratt2D`](@ref), [`BakerTruss`](@ref), [`TrussFrame`](@ref), [`SpaceFrame`](@ref), [`SpaceFrameBeam`](@ref), [`Frame`](@ref), [`GridFrame`](@ref), [`GridNetwork`](@ref) (FDM), and the [`XGroundStructure`](@ref) / [`DenseGroundStructure`](@ref) / [`BoundedGroundStructure`](@ref) family with [`to_truss`](@ref) / [`to_frame`](@ref).
