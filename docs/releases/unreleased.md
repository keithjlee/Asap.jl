# Unreleased (after v1.4.0)

## Added
- `mesh_polygon(…; segments)`: interior lines the element edges must follow, given as point pairs
  `(a, b)` (section cuts through the material, interior load or bearing lines). Their nodes come from
  `nodes_on_segment`, so `section_cut` integrates the exact force across any chosen line. A segment
  may end on the outline, a hole or another segment; a segment leaving the material, or crossing an
  edge or another segment, is refused.
- Walkthrough `examples/shells/07_dapped_end.jl`: a dapped-end beam. Bearing pads as pressures with
  rigid-body supports only, exact interface cut forces, the nib tie and the hanger as integrated cut
  forces, the crack direction at the re-entrant corner, and a mesh study.
