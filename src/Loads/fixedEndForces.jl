"""
fixed end forces for point load
"""
function Qlocal(load::PointLoad)
    #values
    LCS = load.element.LCS
    l = load.element.length

    # load vector in LCS
    plocal = load.element.R[1:3, 1:3] * load.value .* LCS

    #axial end forces
    ax1 = ax2 = - dot(plocal[1], LCS[1]) / 2

    #position of load from start node
    a = load.position * l
    b = l - a #remainder

    # perpendicular load in local Y
    py = - dot(plocal[2], LCS[2])

    # moments in local Z
    mz1 = py * b^2 * a / l^2 
    # mz2 = -py * a^2 * b / l^2
    mz2 = -py * a^2 * b / l^2

    # shear in local Y
    # vy2 = (py * a - mz1 - mz2) / l
    # vy2 = (py * a - mz1 - mz2) / l
    # vy1 = py - vy2
    vy1 = py * b^2 / l^3 * (3a + b)
    vy2 = py * a^2 / l^3 * (a + 3b)

    # perpendicular load in local Z
    pz = - dot(plocal[3], LCS[3])

    # moments in local Y
    my1 = -pz * b^2 * a / l^2
    my2 = pz * a^2 * b / l^2

    # shear in local Z
    # vz2 = (pz * a - my1 - my2) / l
    # vz2 = (pz * a - my1 + my2) / l
    # vz1 = pz - vz2
    vz1 = pz * b^2 / l^3 * (3a + b)
    vz2 = pz * a^2 / l^3 * (a + 3b)

    return [ax1, vy1, vz1, 0., my1, mz1, ax2, vy2, vz2, 0., my2, mz2]
end

"""
fixed end forces for distributed load
"""
function Qlocal(load::LineLoad)
    LCS = load.element.LCS
    l = load.element.length

    # load vector in LCS
    plocal = load.element.R[1:3, 1:3] * load.value .* LCS

    #axial end forces
    ax1 = ax2 = - dot(plocal[1], LCS[1]) * l / 2

    #perpendicular load in local y
    py = - dot(plocal[2], LCS[2])

    vy1 = vy2 = py * l / 2 #shears in Y
    mz1 = py * l^2 / 12 #moment 1 in Z
    mz2 = -mz1 #moment 2 in Z


    # perpendicular load in local z
    pz = -dot(plocal[3], LCS[3])

    vz1 = vz2 = pz * l / 2 #shears
    my1 = -pz * l^2 / 12 #moment 1 in Y
    my2 = -my1 #moment 2 in y

    return [ax1, vy1, vz1, 0., my1, mz1, ax2, vy2, vz2, 0., my2, mz2]
end

"""
fixed end forces for gravity load (special case of line load)
"""
function Qlocal(load::GravityLoad)

    LCS = load.element.LCS
    l = load.element.length
    value = [0., 0., -1.] .* load.element.ρ .* load.element.section.A

    # load vector in LCS
    plocal = element.R[1:3, 1:3] * value .* LCS

    #axial end forces
    ax1 = ax2 = - dot(plocal[1], LCS[1]) * l / 2

    #perpendicular load in local y
    py = - dot(plocal[2], LCS[2])

    vy1 = vy2 = py * l / 2 #shears in Y
    mz1 = py * l^2 / 12 #moment 1 in Z
    mz2 = -mz1 #moment 2 in Z

    # perpendicular load in local z
    pz = -dot(plocal[3], LCS[3])

    vz1 = vz2 = pz * l / 2 #shears
    my1 = -pz * l^2 / 12 #moment 1 in Y
    my2 = -my1 #moment 2 in y

    return [ax1, vy1, vz1, 0., my1, mz1, ax2, vy2, vz2, 0., my2, mz2]
end

"""
Free-fixed Qf
"""
function Q_freefixed(load::Load)
    #length of element
    factor = 3 / 2 / load.element.length
    #fixed end components
    FAb, FSby, FSbz, FTb, FMby, FMbz, FAe, FSey, FSez, FTe, FMey, FMez = Qlocal(load)

    #modified fixed end forces
    return [FAb, #axial beginning
        FSby - factor*FMbz, #shear local Y
        FSbz + factor*FMby, #shear local Z
        0, #torsion beginning
        0, #moment local Y
        0, #moment local Z
        FAe, # axial end
        FSey + factor*FMbz, #shear local Y
        FSez - factor*FMby, #shear local Z
        FTb + FTe, #torsion end
        FMey - 1/2*FMby, #moment local Y
        FMez - 1/2*FMbz] #moment local Z
end

"""
Fixed-free Qf
"""
function Q_fixedfree(load::Load)
    #length of element
    factor = 3 / 2 / load.element.length
    #fixed end components
    FAb, FSby, FSbz, FTb, FMby, FMbz, FAe, FSey, FSez, FTe, FMey, FMez = Qlocal(load)

    #modified fixed end forces
    return [FAb, 
        FSby - factor*FMez,
        FSbz + factor*FMey,
        FTb + FTe,
        FMby - 1/2 * FMey,
        FMbz - 1/2 * FMez,
        FAe,
        FSey + factor*FMez,
        FSez - factor*FMey,
        0,
        0,
        0]
end

"""
Free-free Qf
"""
function Q_freefree(load::Load)
    #length of element
    factor = 1 / load.element.length
    #fixed end components
    FAb, FSby, FSbz, FTb, FMby, FMbz, FAe, FSey, FSez, FTe, FMey, FMez = Qlocal(load)

    #modified fixed end forces
    return [FAb, 
        FSby - factor*(FMbz + FMez),
        FSbz + factor*(FMby + FMey),
        FTb,
        0,
        0,
        FAe,
        FSey + factor*(FMbz + FMez),
        FSez - factor*(FMby + FMey),
        FTe,
        0,
        0]
end

function Q_joist(load::Load)
    #fixed end components
    FAb, FSby, FSbz, FTb, FMby, FMbz, FAe, FSey, FSez, FTe, FMey, FMez = Qlocal(load)

    #modified fixed end forces
    return [FAb, 
        FSby,
        FSbz,
        FTb,
        0,
        0,
        FAe,
        FSey,
        FSez,
        FTe,
        0,
        0]
end

"""
Fixed-fixed Qf
"""
function Q_fixedfixed(load::Load)
    return Qlocal(load)
end

"""
Map of release to proper Q function
"""
qDict = Dict(:fixedfixed => Q_fixedfixed,
    :freefixed => Q_freefixed,
    :fixedfree => Q_fixedfree,
    :freefree => Q_freefree,
    :joist => Q_joist)

"""
Generate the fixed-end forces for a given load type
"""
function Q(load::Load)
    #appropriate function
    qFunction = qDict[load.element.release]
    return qFunction(load)
end