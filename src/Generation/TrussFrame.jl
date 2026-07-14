"""
    TrussFrame <: AbstractGenerator

A 2D trussed tower/frame in the XY plane: `nx × ny` braced bays stacked
vertically, generated and solved on construction.

# Fields
- `model::Model{Float64}` solved structural model
- `nx::Integer` number of bays in x
- `dx::Real` bay width [length]
- `ny::Integer` number of stories in y
- `dy::Real` story height [length]
- `section::AbstractSection` cross section applied to all members
- `single_base::Bool` support only the central base node (true) or all base nodes (false)
"""
struct TrussFrame <: AbstractGenerator
    model::Model{Float64}
    nx::Integer
    dx::Real
    ny::Integer
    dy::Real
    section::AbstractSection
    single_base::Bool
    loaded_nodes::Symbol

    function TrussFrame(nx::Integer, dx::Real, ny::Integer, dy::Real, section::AbstractSection, load = [0., -10., 0.]; single_base = false, loaded_nodes = :row2)

        #checks
        @assert nx % 2 == 0 "nx must be even"
        @assert in(loaded_nodes, [:col1, :col2, :col3, :col4, :row1, :row2]) "loaded_nodes must be one of: :col1, :col2, :col3, :col4, :row1, :row2"
    
        #left most columns
        col1 = [Node([0., dy * n, 0.], :free, :col1) for n = 0:ny]

        #pin support
        fixnode!(col1[1], :pinned)
        col1[1].id = :support
        
        #col2 
        col2 = [Node([dx, dy * n, 0.], :free, :col2) for n = 0:ny]

        #pin support
        fixnode!(col2[1], :pinned)
        col2[1].id = :support
        
        #col1
        col3 = [Node([(nx-1) * dx, dy * n, 0.], :free, :col3) for n = 0:ny]

        #pin support
        fixnode!(col3[1], :pinned)
        col3[1].id = :support
        
        #col2 
        col4 = [Node([nx * dx, dy * n, 0.], :free, :col4) for n = 0:ny]

        #pin support
        fixnode!(col4[1], :pinned)
        col4[1].id = :support
        
        #row1
        row1 = [Node([dx + dx * n, ny * dy, 0.], :free, :row1) for n = 1:nx-3]
        
        #row 2
        row2 = [Node([dx + dx * n, (ny - 1) * dy, 0.], :free, :row2) for n = 1:nx-3]
        
        #column chords
        e_c1 = single_base ? [TrussElement(col1[[i, i+1]]..., section, :col1chord) for i = 2:ny] : [TrussElement(col1[[i, i+1]]..., section, :col1chord) for i = 1:ny]
        
        e_c2 = [TrussElement(col2[[i, i+1]]..., section, :col2chord) for i = 1:ny]
    
        e_c3 = [TrussElement(col3[[i, i+1]]..., section, :col3chord) for i = 1:ny]

        e_c4 = single_base ? [TrussElement(col4[[i, i+1]]..., section, :col4chord) for i = 2:ny] : [TrussElement(col4[[i, i+1]]..., section, :col4chord) for i = 1:ny]
        
        #beam chords
        e_r1 = [TrussElement(row1[[i, i+1]]..., section, :row1chord) for i = 1:nx-4]
        push!(e_r1, TrussElement(col2[end], row1[1], section, :row1chord))
        push!(e_r1, TrussElement(row1[end], col3[end], section, :row1chord))
        
        e_r2 = [TrussElement(row2[[i, i+1]]..., section, :row2chord) for i = 1:nx-4]
        push!(e_r2, TrussElement(col2[end-1], row2[1], section, :row2chord))
        push!(e_r2, TrussElement(row2[end], col3[end-1], section, :row2chord))
        
        #column webs (left)
        w_c12d = [TrussElement(col2[i], col1[i+1], section, :col12web) for i = 1:ny] #diagonals
        w_c12h = [TrussElement(c1, c2, section, :col12web) for (c1, c2) in zip(col1[2:end], col2[2:end])]
        
        w_c12 = [w_c12d; w_c12h]
        
        #column webs (right)
        w_c34d = [TrussElement(col3[i], col4[i+1], section, :col34web) for i = 1:ny]
        w_c34h = [TrussElement(c3, c4, section, :col34web) for (c3, c4) in zip(col3[2:end], col4[2:end])]
        
        w_c34 = [w_c34d; w_c34h]
        
        #beam webs
        w_rv = [TrussElement(r1, r2, section, :rowweb) for (r1, r2) in zip(row1, row2)]
        w_rd1 = [TrussElement(row1[i], row2[i+1], section, :rowweb) for i = 1:Int(floor((nx-3)/2))]
        push!(w_rd1, TrussElement(col2[end], row2[1], section, :rowweb))
        w_rd2 = [TrussElement(row2[i], row1[i+1], section, :rowweb) for i = Int(ceil((nx-3)/2)):nx-4]
        push!(w_rd2, TrussElement(col3[end], row2[end], section, :rowweb))
        
        w_r = [w_rv; w_rd1; w_rd2]
        
        nodes = [col1; col2; col3; col4; row1; row2]
        elements = [e_c1; e_c2; e_c3; e_c4; e_r1; e_r2; w_c12; w_c34; w_r]
        loads = [NodeForce(node, load) for node in nodes[loaded_nodes]]
        
        model = Model(nodes, elements, loads)
        planarize!(model)
        solve!(model)

        new(
            model,
            nx,
            dx,
            ny,
            dy,
            section,
            single_base,
            loaded_nodes
        )
    end

end