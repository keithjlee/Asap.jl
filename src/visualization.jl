```
Functions and packages for visualizing structural analysis
```

# define custom colors for tension/compression
gray1 = colorant"#B2B2B2"
const blue = colorant"#3EA8DE"
const pink = colorant"#FF7BAC"
const red = colorant"#EA2127"
const redWhiteBlue = cgrad([pink, :white, blue])
const green  = colorant"#48B674"
const gray2blue = cgrad([gray1, blue])

# custom arrowheads for support conditions
const pin = load(joinpath(@__DIR__, "customObjects", "pinSupport.obj"))
const fix = load(joinpath(@__DIR__, "customObjects", "fixSupport.obj"))

function fixedArrow(N)
    if N == 2
        return '■' 
    else
        return fix
    end
end

function pinArrow(N)
    if N == 2
        return '▲'
    else
        return pin
    end
end

# converts loads to proper position + vector for plotting
function loadConverter(structure::Structure; scaleFactor = :auto)
    loads = vcat([load.load for load in structure.loads]...)
    
    denom = maximum(abs.(loads))

    if scaleFactor == :auto
        lengths = [element.length for element in structure.elements]
        scaleFactor = mean(lengths) * 0.3
    end

    positions = [load.position for load in structure.loads]
    
    if structure.dims == 2
        loads = [load.load[1:2] ./ denom .* scaleFactor for load in structure.loads]
        return [[p[1] for p in positions], [p[2] for p in positions], [l[1] for l in loads], [l[2] for l in loads]]
        # return [[p[1] for p in positions], [p[2] for p in positions], [0.0 for p in positions], [l[1] for l in loads], [l[2] for l in loads], [0.0 for l in loads]]
    else
        loads = [load.load[1:3] ./ denom .* scaleFactor for load in structure.loads]
        return [[p[1] for p in positions], [p[2] for p in positions], [p[3] for p in positions], [l[1] for l in loads], [l[2] for l in loads], [l[3] for l in loads]]
    end
end

# converting supports to plottable vectors (to be optimized)
function pinDofConverter(structure::Structure)
    positions = [node.position for node in structure.nodes]
    dofs = [node.DOFS for node in structure.nodes]

    l = length(positions)

    x = Float64[]
    y = Float64[]
    z = Float64[]
    u = Float64[]
    v = Float64[]
    w = Float64[]

    if structure.dims == 2
        for i = 1:l
            if all(dofs[i])
                continue
            else
                if dofs[i][1] == false
                    push!(x, positions[i][1] - 1)
                    push!(y, positions[i][2])
                    push!(z, 0)
                    push!(u, 1)
                    push!(v, 0)
                    push!(w, 0)
                end
                if dofs[i][2] == false
                    push!(x, positions[i][1])
                    push!(y, positions[i][2] - 1)
                    push!(z, 0)
                    push!(u, 0)
                    push!(v, 1)
                    push!(w, 0)
                end
            end
        end

        return x, y, z, u, v, w
    else
        for i = 1:l
            if all(dofs[i])
                continue
            else
                if dofs[i][1] == false
                    push!(x, positions[i][1] - 1)
                    push!(y, positions[i][2])
                    push!(z, positions[i][3])
                    push!(u, 1)
                    push!(v, 0)
                    push!(w, 0)
                end
                if dofs[i][2] == false
                    push!(x, positions[i][1])
                    push!(y, positions[i][2] - 1)
                    push!(z, positions[i][3])
                    push!(u, 0)
                    push!(v, 1)
                    push!(w, 0)
                end
                if dofs[i][3] == false
                    push!(x, positions[i][1])
                    push!(y, positions[i][2])
                    push!(z, positions[i][3] - 1)
                    push!(u, 0)
                    push!(v, 0)
                    push!(w, 1)
                end
            end
        end

        return x, y, z, u, v, w
    end
end

function fixDofConverter(structure::Structure)
    positions = [node.position for node in structure.nodes]
    dofs = [node.DOFS for node in structure.nodes]

    l = length(positions)

    x = Float64[]
    y = Float64[]
    z = Float64[]
    u = Float64[]
    v = Float64[]
    w = Float64[]

    if structure.dims == 2
        for i = 1:l
            if all(dofs[i]) || length(dofs[i]) == 2
                continue
            else
                if dofs[i][3] == false
                    push!(x, positions[i][1])
                    push!(y, positions[i][2])
                    push!(z, -1)
                    push!(u, 0)
                    push!(v, 0)
                    push!(w, 1)
                end
            end
        end

        return x, y, z, u, v, w
    else
        for i = 1:l
            if all(dofs[i]) || length(dofs[i]) == 3
                continue
            else
                if dofs[i][4] == false
                    push!(x, positions[i][1] - 1)
                    push!(y, positions[i][2])
                    push!(z, positions[i][3])
                    push!(u, 1)
                    push!(v, 0)
                    push!(w, 0)
                end
                if dofs[i][5] == false
                    push!(x, positions[i][1])
                    push!(y, positions[i][2] - 1)
                    push!(z, positions[i][3])
                    push!(u, 0)
                    push!(v, 1)
                    push!(w, 0)
                end
                if dofs[i][6] == false
                    push!(x, positions[i][1])
                    push!(y, positions[i][2])
                    push!(z, positions[i][3] - 1)
                    push!(u, 0)
                    push!(v, 0)
                    push!(w, 1)
                end
            end
        end

        return x, y, z, u, v, w
    end

end

function reactionConverter(structure::Structure)
    if structure.dims == 2
        if structure.elements[1].type == :truss
            gap = 2
        else
            gap = 3
        end
        group = 1
        reactions = Vec3.([[structure.reactions[i:i+group]; 0.0] for i = 1:gap:structure.nDOFS])
    else
        if structure.elements[1].type == :truss
            gap = 3
        else
            gap = 6
        end
        group = 2
        reactions = Vec3.([structure.reactions[i:i+group] for i = 1:gap:structure.nDOFS])
    end
    reactions ./= maximum(norm.(reactions))
    reactions .*= 0.5mean([element.length for element in structure.elements])
    return reactions
end

######################

function axo(geo::Geometry; 
    cmforce = redWhiteBlue,
    mode = :undisplaced,
    forces = true,
    scale = false,
    elev = pi/6,
    az = pi/4,
    lw = 2,
    lcoverride = nothing,
    bgc = :transparent,
    arrowSize = 0.15,
    res = (3000,3000))

    fig = Figure(resolution = res, backgroundcolor = bgc)

    if scale == false
        ax1 = Axis3(fig[1,1],
            backgroundcolor = bgc,
            elevation = elev,
            azimuth = az,
            aspect = :data)
        hidedecorations!(ax1)
        hidespines!(ax1)
    else
        ax1 = Axis3(fig[1,1],
            backgroundcolor = bgc,
            elevation = elev,
            azimuth = az,
            aspect = :data)
    end

    if mode == :undisplaced
        els = geo.elements
        normalizedAreas = [e.A for e in geo.structure.elements]
        normalizedAreas ./= maximum(normalizedAreas)
        thickness = lw .* normalizedAreas
        lc = - ones(geo.structure.nElements)
        cm = :grays
        cr = (-1,1)
    else
        els = geo.displacedElements
        peakforce = maximum(abs.(geo.axialForce))
        thickness = 5 .* lw .* abs.(geo.axialForce) ./ peakforce
        lc = geo.axialForce
        cm = cmforce
        cr = (-peakforce, peakforce)
    end

    # override linecolor if desired
    if !isnothing(lcoverride)
        lc = lcoverride
    end

    linesegments!(vcat(els...),
        color = lc,
        colormap = cm,
        colorrange = cr,
        linewidth = thickness)

    # arrows!(geo.pinDOFS...,
    #     arrowhead = pin,
    #     arrowcolor = (:black,0.6),
    #     arrowsize = 2 * arrowSize,
    #     linewidth = 0)

    if mode == :undisplaced && forces
        arrows!(geo.loads...,
            color = (red, 0.5),
            arrowsize = arrowSize,
            linewidth = arrowSize / 3)
    end

    ax = Axis3(fig[2,1],
        backgroundcolor = bgc,
        elevation = 0,
        azimuth = 0,
        aspect = :data)

    hidespines!(ax)
    hidedecorations!(ax)

    linesegments!(vcat(els...),
        color = lc,
        colormap = cm,
        colorrange = cr,
        linewidth = thickness)

    # arrows!(geo.pinDOFS...,
    #     arrowhead = pin,
    #     arrowcolor = (:black,0.6),
    #     arrowsize = 2 * arrowSize,
    #     linewidth = 0)

    ax3 = Axis3(fig[2,2],
        backgroundcolor = bgc,
        elevation = 0,
        azimuth = pi/2,
        aspect = :data)

    hidedecorations!(ax3)
    hidespines!(ax3)

    linesegments!(vcat(els...),
        color = lc,
        colormap = cm,
        colorrange = cr,
        linewidth = thickness)

    # arrows!(geo.pinDOFS...,
    #     arrowhead = pin,
    #     arrowcolor = (:black,0.6),
    #     arrowsize = 2 * arrowSize,
    #     linewidth = 0)

    ax4 = Axis3(fig[1,2],
        backgroundcolor = bgc,
        elevation = pi/2,
        azimuth = pi/2,
        aspect = :data)
        
    hidedecorations!(ax4)
    hidespines!(ax4)

    linesegments!(vcat(els...),
        color = lc,
        colormap = cm,
        colorrange = cr,
        linewidth = thickness)

    # arrows!(geo.pinDOFS...,
    #     arrowhead = pin,
    #     arrowcolor = (:black,0.6),
    #     arrowsize = 2 * arrowSize,
    #     linewidth = 0)

    return fig
end

#######################

function structurePlot(geo::Geometry; 
        res = (1200,800),
        lineSize = 2,
        lineColor = :black, 
        backgroundColor = :transparent,
        colorMap = redWhiteBlue,
        nodeColor = green,
        elementColor = :black,
        lineSizeRange = 0.1:0.1:10,
        arrowSize = 0.3, 
        arrowSizeRange = 0.05:0.05:5,
        textSize = 1,
        textSizeRange = 0.1:0.1:10,
        scaleToForce = false,
        constructionPath = false)
    
    # Initialize
    # fig, ax = scatter(geo.nodes, 
    #     markersize = 0,
    #     figure = (resolution = res, 
    #     backgroundcolor = backgroundColor),
    #     show_axis = false)

    fig = Figure(resolution = res)
    ax = Axis3(fig[1,1],
        aspect = :data)

    # Alternative method using Axis3... provides orthographic viewmode
    # fig = Figure(resolution = res)

    # ax, nds = meshscatter(fig[1,1],
    #     geo.nodes, markersize = 0, 
    #     axis = (type = Axis3, perspective = 0))
    # hidedecorations!(ax)
    # hidespines!(ax)

    # create colorrange
    colorLimit = maximum(abs.(geo.axialForce))


    # Create adjustment sliders
    lsgrid = labelslidergrid!(
        fig,
        ["Line Scale", "Arrow Scale", "Text Size"],
        [lineSizeRange, arrowSizeRange, textSizeRange];
        width = Auto(),
        # height = Auto(),
        # tellheight = false
    )

    areas = geo.areas ./ maximum(geo.areas)
    # Variables for sliders
    element_thickness = lift(s-> lineSize .* areas .* s, lsgrid.sliders[1].value)
    if scaleToForce
        displaced_thickness = lift(s -> lineSize .* abs.(geo.axialForce) ./ maximum(geo.axialForce) .* s, lsgrid.sliders[1].value)
    else
        displaced_thickness = lift(s-> 2 .* areas .* s, lsgrid.sliders[1].value)
    end
    arrow_thickness = lift(s-> s / 3, lsgrid.sliders[2].value)
    arrow_size = lift(s-> s, lsgrid.sliders[2].value)
    ts = lift(s-> 10 * s, lsgrid.sliders[3].value)

    if constructionPath
        path_thickness = lift(s -> 4 .* s, lsgrid.sliders[1].value)
    end

    # Set defaults
    set_close_to!(lsgrid.sliders[1], lineSize)
    set_close_to!(lsgrid.sliders[2], arrowSize)
    set_close_to!(lsgrid.sliders[3], textSize)

    # Add sliders
    fig[2,1] = lsgrid.layout

    # Variables for toggles
    if constructionPath
        toggles = [Toggle(fig, active = active) for active in [false, true, false, true, false]]
        labels = [Label(fig, x)  for x in ["Displaced", "Loads", "Labels", "DOFs", "Path"]]
        opposite = lift(t -> !t[], toggles[1].active)
    else
        toggles = [Toggle(fig, active = active) for active in [false, true, false, true]]
        labels = [Label(fig, x)  for x in ["Displaced", "Loads", "Labels", "DOFs"]]
        opposite = lift(t -> !t[], toggles[1].active)
    end

    # Add toggles
    fig[1,2] = grid!(hcat(toggles, labels), tellheight = false)

    # Plot features
    bars = linesegments!(vcat(geo.elements...), 
        color = lineColor, 
        linewidth= element_thickness, 
        visible = opposite)
    loadArrows = arrows!(geo.loads..., 
        color = (red, 0.5), 
        arrowsize = arrow_size, 
        linewidth = arrow_thickness)
    reactionArrows = arrows!(geo.nodes .- geo.reactions, geo.reactions,
        color = keithGreen,
        arrowsize = arrow_size,
        linewidth = arrow_thickness)
    nodeLabels = text!(geo.nodeLabels, position = geo.nodes, 
        align = (:center, :bottom), 
        color = nodeColor, 
        textsize=  ts)
    elementLabels = text!(geo.elementLabels, position = geo.midpoints, 
        color = elementColor, 
        textsize = ts)
    displacedShape = linesegments!(vcat(geo.displacedElements...), 
        linewidth = displaced_thickness, 
        color = geo.axialForce, 
        colormap = colorMap,
        colorrange = (-colorLimit, colorLimit),
        visible = false)
    pinDofs = arrows!(geo.pinDOFS..., 
        arrowhead = pin, 
        arrowsize = arrow_size, 
        arrowcolor = (:black,0.6), 
        linewidth = 0)
    fixDofs = arrows!(geo.fixDOFS..., 
        arrowhead = fix, 
        arrowsize = arrow_size, 
        arrowcolor = (:black,0.6), 
        linewidth = 0)

    # plot node construction path (only relevant if assembly sequence is of interest and/or list of nodes was ordered)
    if constructionPath
        cp = vcat([geo.nodes[1]], [[geo.nodes[i], geo.nodes[i]] for i = 2:length(geo.nodes)-1]..., [geo.nodes[end]])
        path = linesegments!(cp, 
            color = 0:length(cp)-1, 
            colormap = gray2blue, 
            linewidth = path_thickness)
        connect!(path.visible, toggles[5].active)
    end

    # Connect toggles to features
    connect!(displacedShape.visible, toggles[1].active)
    connect!(loadArrows.visible, toggles[2].active)
    connect!(reactionArrows.visible, toggles[2].active)
    connect!(nodeLabels.visible, toggles[3].active)
    connect!(elementLabels.visible, toggles[3].active)
    connect!(pinDofs.visible, toggles[4].active)
    connect!(fixDofs.visible, toggles[4].active)
    
    # show figure
    display(fig)
    return fig
end

# Alternative method of Geometry structure is not created
function structurePlot(structure::Structure)
    geo = Geometry(structure)
    structurePlot(geo)
end

# Plotting before analysis (useful when creating structure piece by piece to check for errors)
function prePlot(structure::Structure; lineSize = 2,
    lineColor = :black, 
    lineSizeRange = 0.1:0.1:10,
    arrowSize = 0.15, 
    arrowSizeRange = 0.05:0.05:5,
    textSize = 1,
    textSizeRange = 0.1:0.1:10,
    scaleToForce = false)

    # initialize geometry
    geo = Geometry(structure)

    if isdefined(structure, :loads)
        # Initialize
        fig, ax = scatter(geo.nodes, 
        figure = (resolution = (1200,800),),
        show_axis = false)

        # Create adjustment sliders
        lsgrid = labelslidergrid!(
            fig,
            ["Line Scale", "Arrow Scale", "Text Size"],
            [lineSizeRange, arrowSizeRange, textSizeRange];
            width = Auto(),
            # height = Auto(),
            # tellheight = false
        )

        # Variables for sliders
        element_thickness = lift(s-> geo.areas .* s, lsgrid.sliders[1].value)
        if scaleToForce
            displaced_thickness = lift(s -> lineSize .* abs.(geo.axialForce) ./ maximum(geo.axialForce) .* s, lsgrid.sliders[1].value)
        else
            displaced_thickness = lift(s-> 2 .* geo.areas .* s, lsgrid.sliders[1].value)
        end
        arrow_thickness = lift(s-> s / 3, lsgrid.sliders[2].value)
        arrow_size = lift(s-> s, lsgrid.sliders[2].value)
        ts = lift(s-> 10 * s, lsgrid.sliders[3].value)

        # Set defaults
        set_close_to!(lsgrid.sliders[1], lineSize)
        set_close_to!(lsgrid.sliders[2], arrowSize)
        set_close_to!(lsgrid.sliders[3], textSize)

        # Add sliders
        fig[2,1] = lsgrid.layout

        # Variables for toggles
        toggles = [Toggle(fig, active = active) for active in [true, false, true]]
        labels = [Label(fig, x)  for x in ["Loads", "Labels", "DOFs"]]

        # Add toggles
        fig[1,2] = grid!(hcat(toggles, labels), tellheight = false)

        # Plot features
        bars = linesegments!(vcat(geo.elements...), 
            color = lineColor, 
            linewidth= element_thickness)
        loadArrows = arrows!(geo.loads..., 
            color = (red, 0.5), 
            arrowsize = arrow_size, 
            linewidth = arrow_thickness)
        nodeLabels = text!(geo.nodeLabels, position = geo.nodes, align = (:center, :bottom), 
            color = green, 
            textsize=  ts)
        elementLabels = text!(geo.elementLabels, 
            position = geo.midpoints, 
            color = :black, 
            textsize = ts)
        pinDofs = arrows!(geo.pinDOFS..., 
            arrowhead = pin, 
            arrowsize = arrow_size, 
            arrowcolor = (:black, 0.8), 
            linewidth = 0)
        fixDofs = arrows!(geo.fixDOFS..., 
            arrowhead = fix, 
            arrowsize = arrow_size, 
            arrowcolor = (:black, 0.8), 
            linewidth = 0)

        # Connect toggles to features
        connect!(loadArrows.visible, toggles[1].active)
        connect!(nodeLabels.visible, toggles[2].active)
        connect!(elementLabels.visible, toggles[2].active)
        connect!(pinDofs.visible, toggles[3].active)
        connect!(fixDofs.visible, toggles[3].active)


        display(fig)
    else
        # Initialize
        fig, ax = scatter(geo.nodes, 
        figure = (resolution = (1200,800),),
        show_axis = false)

        # Create adjustment sliders
        lsgrid = labelslidergrid!(
            fig,
            ["Line Scale", "DOF Scale", "Text Size"],
            [lineSizeRange, arrowSizeRange, textSizeRange];
            width = Auto(),
            # height = Auto(),
            # tellheight = false
        )

        # Variables for sliders
        element_thickness = lift(s-> geo.areas .* s, lsgrid.sliders[1].value)
        arrow_size = lift(s-> s, lsgrid.sliders[2].value)
        ts = lift(s-> 10 * s, lsgrid.sliders[3].value)

        # Set defaults
        set_close_to!(lsgrid.sliders[1], lineSize)
        set_close_to!(lsgrid.sliders[2], arrowSize)
        set_close_to!(lsgrid.sliders[3], textSize)

        # Add sliders
        fig[2,1] = lsgrid.layout

        # Variables for toggles
        toggles = [Toggle(fig, active = active) for active in [false, true]]
        labels = [Label(fig, x)  for x in ["Labels", "DOFs"]]

        # Add toggles
        fig[1,2] = grid!(hcat(toggles, labels), tellheight = false)

        # Plot features
        bars = linesegments!(vcat(geo.elements...), 
            color = lineColor, 
            linewidth= element_thickness, 
            visible = opposite)
        nodeLabels = text!(geo.nodeLabels, 
            position = geo.nodes, 
            align = (:center, :bottom), 
            color = green, 
            textsize=  ts)
        elementLabels = text!(geo.elementLabels, 
            position = geo.midpoints, 
            color = :black, 
            textsize = ts)
        pinDofs = arrows!(geo.pinDOFS..., 
            arrowhead = pin, 
            arrowsize = arrow_size, 
            arrowcolor = (:black, 0.8), 
            linewidth = 0)
        fixDofs = arrows!(geo.fixDOFS..., 
            arrowhead = fix, 
            arrowsize = arrow_size, 
            arrowcolor = (:black, 0.8), 
            linewidth = 0)

        # Connect toggles to features
        connect!(nodeLabels.visible, toggles[1].active)
        connect!(elementLabels.visible, toggles[1].active)
        connect!(pinDofs.visible, toggles[2].active)
        connect!(fixDofs.visible, toggles[2].active)


        display(fig)
        
    end
end


function plot(geo::Geometry;
        showaxis = true)

    fig = Figure()

    if length(geo.nodes[1]) == 2
        ax = Axis(fig[1,1],
            aspect = DataAspect())
    else
        ax = Axis3(fig[1,1],
            aspect = :data)
    end

    # hide axes
    if !showaxis
        hidedecorations!(ax)
        hidespines!(ax)
    end

    colorLimit = maximum(abs.(geo.axialForce))

    # Slider grid
    lsgrid = SliderGrid(fig[2,1],
        (label = "Line Scale", range = 0.1:0.1:10, startvalue = 1.),
        (label = "Object Scale", range = 0.1:0.1:10, startvalue = 1.),
        (label = "Text Scale", range = 0.1:0.1:10, startvalue = 1.))

    # Slider values
    areas = geo.areas ./ maximum(geo.areas)
    lw = lift(lsgrid.sliders[1].value) do v
        areas .* v
    end

    forces = geo.axialForce ./ maximum(abs.(geo.axialForce))
    lw2 = lift(lsgrid.sliders[1].value) do v
        2 .* forces .* v
    end

    objectSize1 = lift(lsgrid.sliders[2].value) do v
        v
    end

    objectSize2 = lift(lsgrid.sliders[2].value) do v
        v / 3
    end

    textSize = lift(lsgrid.sliders[3].value) do v
        v * 10
    end

    #toggles
    toggles = [Toggle(fig, active = active) for active in [false, true, false, true]]
    labels = [Label(fig, "Displaced"),
        Label(fig, "Loads"),
        Label(fig, "Labels"),
        Label(fig, "DOFs")]

    opp = lift(t -> !t[], toggles[1].active)
    
    fig[1,2] = grid!(hcat(toggles, labels), tellheight = false)

    # plot features
    bars = linesegments!(vcat(geo.elements...),
        color = :black,
        linewidth = lw,
        visible = opp)

    loadArrows = arrows!(geo.loads...,
        color = (red, 0.5),
        arrowsize = objectSize1,
        linewidth = objectSize2)

    connect!(loadArrows.visible, toggles[2].active)

    reactionArrows = arrows!(geo.nodes .- geo.reactions, geo.reactions,
        color = green,
        arrowsize = objectSize1,
        linewidth = objectSize2)

    connect!(reactionArrows.visible, toggles[2].active)

    nodeLabels = text!(geo.nodeLabels, position = geo.nodes,
        align = (:center, :bottom),
        color = green,
        textsize = textSize)

    connect!(nodeLabels.visible, toggles[3].active)

    elementLabels = text!(geo.elementLabels, position = geo.midpoints,
        color = :gray,
        textsize = textSize)

    connect!(elementLabels.visible, toggles[3].active)

    displacedBars = linesegments!(vcat(geo.displacedElements...),
        linewidth = lw2,
        color = geo.axialForce,
        colormap = redWhiteBlue,
        colorrange = (-colorLimit, colorLimit))

    connect!(displacedBars.visible, toggles[1].active)

    pinDofs = arrows!(geo.pinDOFS..., 
        arrowhead = pin, 
        arrowsize = objectSize1, 
        arrowcolor = (:black,0.6), 
        linewidth = 0)

    connect!(pinDofs.visible, toggles[4].active)

    fixDofs = arrows!(geo.fixDOFS..., 
        arrowhead = fix, 
        arrowsize = objectSize1, 
        arrowcolor = (:black,0.6), 
        linewidth = 0)

    connect!(fixDofs.visible, toggles[4].active)

    display(fig)
    return fig
end