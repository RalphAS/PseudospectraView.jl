module PSApp


# WARNING: some comments are obsolete!
# This is an incomplete refactor of a 3-year old code, a Julian dinosaur.

using QML
# Python must not manage a gui:
ENV["MPLBACKEND"] = "Agg"
using Plots
using Pseudospectra
PseudospectraPlots = Pseudospectra.PseudospectraPlots
PlotsGUIState = Pseudospectra.PseudospectraPlots.PlotsGUIState
using PseudospectraQML

using Formatting
PSAData = PseudospectraQML.PSAData

verbosity=2

# shorthand
PSA = Pseudospectra

using Printf: @sprintf
using SparseArrays: issparse
# for schur!
using LinearAlgebra

################################################################
# Global data
# Since the App runs in global scope, we may as well use global vars
# where it's convenient.
# A few more are created by big function calls below.

maindisplay = nothing
mainpltobj = nothing

otherdisplay = nothing
pltobj2 = nothing

# the instance of PSA opts maintained by the GUI
guiopts = Dict{Symbol,Any}()

ps_data = nothing

firstcall = true

@enum DirVsIter dvi_unset dvi_direct dvi_iter
dvi = dvi_unset


# flags indicating whether items changed by GUI require action
need_recomp = false
need_redraw = false

# In order to force updates of the GUI window, we use tasks
# and a QTimer. For now we use produce/consume logic.
compute_chnl = nothing
running = false

initmsg() = "Welcome to PseudospectraView"

const axname2idx = Dict(:xMin=>1,:xMax=>2,:yMin=>3,:yMax=>4)

tmpdata = nothing

# Until we can get a plot location from the pointer and GUI canvas,
# we use a dialog which saves its result here:
zselected = NaN+0.0im

################################################################
# QMLFunctions: how the GUI talks to Julia
# Apparently qmlfunctions need to be defined in global scope
# and QML.exec() needs to run in global scope to see them.

# These are called by the App code, apparently only on setup or resize
# Display is generally forced via drawcmd().

function init_backend(width::Float64, height::Float64)
    if width < 5 || height < 5
        return
    end
    # FIXME: allow other backends for Plots.
    pyplot(size=(Int64(round(width)),Int64(round(height))))
    return
end

function mainplot(d::JuliaDisplay)
    global maindisplay
    maindisplay = d
    if mainpltobj != nothing
        display(d, mainpltobj)
    end
    return
end

function otherplot(d::JuliaDisplay)
    global otherdisplay
    otherdisplay = d
    if pltobj2 != nothing
        display(d, pltobj2)
    end
    return
end

qmlfunction("mainplot", mainplot)
qmlfunction("init_backend", init_backend)
qmlfunction("otherplot", otherplot)

# Functions called by various GUI controls

# setoptbydlg() is called when a GUI dialog accepted a value from user.
# Named for principal use case, but can also invoke a compute/plot routine.
function setoptbydlg(optkey::AbstractString,val::AbstractString)
    k = Symbol(optkey)
    global need_recomp
#    println("setting $optkey by dialog to $val")
    ao = get(ps_data.ps_dict,:arpack_opts,nothing)
    curzoom = ps_data.zoom_list[ps_data.zoom_pos]
    ao_current = true
    if (optkey[1:6] == "Arpack") && ao == nothing
        @warn "no arpack opts in ps_data"
        return
    end
    local n,x

    if k == :ArpackMaxiter
        try
            n = parse(Int,val)
        catch JE
            @warn "Invalid entry for ARPACK maxiter"
            return
        end
        if ao.maxiter != n
            ao.maxiter = n
            need_recomp = true
            ao_current = false
        end
    elseif k == :ArpackNcv
        try
            n = parse(Int,val)
        catch JE
            @warn "Invalid entry for ARPACK ncv"
            return
        end
        if ao.ncv != n
            ao.ncv = n
            need_recomp = true
            ao_current = false
        end
    elseif k == :ArpackTol
        try
            x = parse(Float64,val)
        catch JE
            @warn "Invalid entry for ARPACK tol"
            return
        end
        if !isapprox(ao.tol,x)
            ao.tol = x
            need_recomp = true
            ao_current = false
        end
    elseif k == :AbscissaEps
        try
            x = parse(Float64,val)
        catch JE
            @warn "Invalid entry for ϵ"
            return
        end
        if x < 0
            @warn "Invalid entry for ϵ"
        end
        abscissacb(x)
    elseif k == :RadiusEps
        try
            x = parse(Float64,val)
        catch JE
            @warn "Invalid entry for ϵ"
            return
        end
        if x < 0
            @warn "Invalid entry for ϵ"
        end
        radiuscb(x)
    elseif k == :ProjLev
        try
            x = parse(Float64,val)
        catch JE
            @warn "Invalid entry for Projection Level"
            return
        end
        if !isapprox(curzoom.proj_lev, x)
            curzoom.proj_lev = x
            curzoom.computed = false
            ps_data.ps_dict[:proj_valid] = false
            need_recomp = true
        end
    else
        @warn "setoptbydlg: Unexpected key $k"
    end
    # TODO:  (upstream also has v0, v0Cmd)
    if !ao_current
        setarpackgui(ao)
    end
    return
end

function setoptbykey(optkey::AbstractString,val::AbstractString)
    k = Symbol(optkey);
    global need_recomp
    global need_redraw
    goable = false
    local x,n
    if ps_data == nothing
        # FIXME: should we populate opts instead?
        (k == :which) && return # known case
        @warn "setopt invoked before pop for $optkey"
        return
    end
    ps_dict = ps_data.ps_dict
    curzoom = ps_data.zoom_list[ps_data.zoom_pos]
    # Treat ARPACK params shown on main GUI window specially:
    # following upstream we expect user to verify them before "Update".
    if k == :which
        if haskey(ps_dict,:arpack_opts)
            if ps_dict[:arpack_opts].which != Symbol(val)
                ps_dict[:arpack_opts].which = Symbol(val)
                need_recomp = true
            end
        end
    elseif k == :kArpack
        if haskey(ps_dict,:arpack_opts)
            try
                n = parse(Int,val)
            catch
                @warn "Invalid entry for ARPACK k"
                return
            end
            if (n < 1) || (n > size(ps_data.input_matrix,1)-2)
                @warn "Invalid entry for ARPACK k"
            end
            if ps_dict[:arpack_opts].nev != n
                ps_dict[:arpack_opts].nev = n
                ps_dict[:proj_valid] = false
                need_recomp = true
            end
        end
    elseif k == :npts
        try
            n = parse(Int,val)
        catch
            @warn "Invalid entry for npts"
            return
        end
        if n < 5
            @warn "npts must be at least 5; got $n"
            return
        end
        curzoom.npts = n
        need_recomp = true
        goable = true
    elseif k ∈ keys(axname2idx)
        try
            x = parse(Float64,val)
        catch
            @warn "Invalid entry for $k"
            return
        end
        isempty(curzoom.ax) && (curzoom.ax = fill(NaN,(4,)))
        curzoom.ax[axname2idx[k]] = x
        need_recomp = true
        goable = !any(isnan.(curzoom.ax))
    elseif k ∈ [:first,:last,:step]
        try
            x = parse(Float64,val)
        catch
            @warn "Invalid entry for $k"
            return
        end
        setfield!(curzoom.levels,k,x)
        curzoom.autolev = false
        need_redraw = true
        goable = true
    elseif k ∈ [:showimagax, :showunitcircle]
        guiopts[k] = parse(Bool,val)
        need_redraw = true
        goable = true
    elseif k == :showFov
        tf = (val == "true")
        if get(guiopts,:showfov,false) != tf
            need_redraw = true
            goable = true
        end
        guiopts[:showfov] = tf
    else
        @warn "setoptbykey: unrecognized key $k"
    end

    if goable
        @emit enableGo(Int32(1))
    end
    return
end

function smartlevels()
    global need_redraw
    if ps_data == nothing
        @warn "smartlev invoked before population of ps_data"
        return
    end
    curzoom = ps_data.zoom_list[ps_data.zoom_pos]
    # FIXME: needs much more logic for zooming
    t_levels,err = PSA.recalc_levels(curzoom.Z,curzoom.ax)
    if err == -1
        printwarning("Range is too small: no contours to plot. Refine grid or zoom out.")
        return
    elseif err == -2
        printwarning("Matrix is too non-normal: resolvent norm → ∞. Zoom out!")
        return
    end
    lev = PSA.LevelDesc(t_levels)
    curzoom.levels = lev
    curzoom.autolev = true
    @emit setLevels(lev.first,lev.last,lev.step)
    need_redraw = true
    @emit enableGo(Int32(1))
    nothing
end

function myproduce(chnl::Nothing,val)
    println("? cant produce to empty channel; val=",val)
end

myproduce(chnl,val) = put!(chnl,val)

function compute(chnl)
    # println("compute arg is ",typeof(chnl))
    myproduce(chnl,2)
    global running
    running = true
    printinfo("calculating...")
    PSA.driver!(ps_data,guiopts,gs,myprintln=printinfo)
    setedittext(ps_data.zoom_list[ps_data.zoom_pos])
    running = false
    myproduce(chnl,-1)
    nothing
end

# This code runs in the main task, invoked by the QTimer.
# Its role is to receive kicks from the compute task which
# provoke GUI redraws.
function monitor()
#        println("monitor: attempting to take")
    status = take!(compute_chnl)
    if status <= 0
        @emit runTimer(Int32(0))
        printinfo("Done.")
        if compute_chnl == nothing
            println("done but no channel?")
        end
        global compute_chnl = nothing
    end
end

function go()
    global need_recomp
    global need_redraw
    if ps_data == nothing
        @warn "go invoked before population"
        return
    end
    @emit showPlot(Int32(1))
    global compute_chnl
    if need_recomp
#            println("opening channel")
        compute_chnl = Channel(compute)
        @emit runTimer(Int32(1))
        need_recomp = false
        need_redraw = false
    elseif need_redraw
        printinfo("redrawing")
        PseudospectraPlots.redrawcontour(gs,ps_data,guiopts)
        need_redraw = false
    end
#    ax = (mainpltobj != nothing) ? PSA.getxylims(mainpltobj) : zeros(0)
    setedittext(ps_data.zoom_list[ps_data.zoom_pos])
    printinfo("ready")
    @emit enableGo(Int32(0))
end

function setselectedz(zrs::AbstractString,zis::AbstractString)
    global zselected
    zr,zi = NaN,NaN
    try
        zr = parse(Float64,zrs)
        zi = parse(Float64,zis)
    catch
        printwarning("invalid z")
        zr,zi = NaN,NaN
    end
    zselected = zr + zi*im
#    notify(dialogcv)
    nothing # returns to QML-land, which doesn't like Complex
end

function loadmtx(varname::AbstractString,myexpr::AbstractString)
    #  any existing plots are obsolete, so clear
    global mainplotobj, plotobj2
    mainplotobj = nothing
    plotobj2 = nothing
    @emit showPlot(Int32(0))
    @emit showPlot2(Int32(0))

    # this should really be done by the GUI, but this is easier:
    if isempty(varname)
        varname = "A"
    end

    global firstcall
    if firstcall
        printinfo("Please be patient, plotting code needs to be compiled...")
    end

    try
        eval(PSAData,parse("$(varname) = $(myexpr)"))
        eval(Main.PSApp,
             parse("ps_data = Pseudospectra.new_matrix(PSAData.$(varname),guiopts)"))
    catch JE
        printwarning("New Matrix failed. See console.")
        @warn "loadmtx failed, exception was $JE"
        return nothing
    end
    global guiopts
    printbanner("Working on $(varname) = $(myexpr)")
    PSAData.clearopts() # delete obsolete matrix-specific stuff
    guiopts = PSA.fillopts(gs,PSAData.getopts())
    refresh()
    firstcall = false
    (verbosity > 1) && println("matrix loaded")
    nothing
end

function savedata(varname::AbstractString,mykey)
    mymodname = current_module()
    global tmpdata
    # FIXME: test for valid variable name
    if mykey == "portrait"
        tmpdata = ps_data.zoom_list[ps_data.zoom_pos]
    elseif mykey == "ps_data"
        tmpdata = ps_data
    elseif mykey == "zpsradius"
        tmpdata = ps_data.ps_dict[:zpsradius]
    else
        @warn "unrecognized save key $mykey"
        return
    end
    try
        eval(Main,parse("$(varname) = deepcopy($(mymodname).tmpdata)"))
    catch JE
        @warn "savedata failed, exception was $JE"
    end
    tmpdata = nothing
end

function use_eigs(x::Int32)
    global dvi
    direct = (x == 1)
    ps_dict = ps_data.ps_dict
    ps_dict[:direct] = direct
    m,n = size(ps_data.input_matrix)
    if direct
        dvi = dvi_direct
        if !isempty(get(ps_dict,:schur_mtx,zeros(0)))
            # revert to input matrix if dense square
            ps_data.matrix = ps_dict[:schur_mtx]
            ps_dict[:ews] = ps_dict[:orig_ews]
            ps_dict[:ew_estimates] = false
        elseif (m==n) && issparse(ps_data.input_matrix)
            # if sparse square, forget ARPACK compression
            ps_data.matrix = ps_data.input_matrix
        end
        # if there's a unitary matrix from Schur decomposition, use that
        # otherwise use the one give by the user, if any
        if !isempty(get(ps_dict,:schur_unitary_mtx,zeros(0)))
            ps_data.unitary_mtx = (ps_data.input_unitary_mtx
                                     * ps_dict[:schur_unitary_mtx])
        else
            ps_data.unitary_mtx = ps_data.input_unitary_mtx
        end
        ps_dict[:proj_valid] = false
        # if we're reverting to a square matrix, forget ARPACK projection
        ss = size(ps_data.matrix)
        if ss[1] == ss[2]
            ps_dict[:isHessenberg] = false
        end
        if (ps_data.zoom_list[ps_data.zoom_pos].computed
            && !isempty(get(ps_dict,:schur_mtx,zeros(0))))
            @emit toggleFoV(Int32(1))
        end
    else
        dvi = dvi_iter
        if !haskey(ps_dict,:arpack_opts)
            ps_dict[:arpack_opts] = PSA.ArpackOptions{eltype(ps_data.matrix)}()
            @warn "setting default ARPACK options (from GUI)"
        end
        @emit toggleFoV(Int32(0))
    end
    global need_recomp
    need_recomp = true
    @emit enableGo(Int32(1))
end

function mainzoom(inorout::Int32, xpxl::Float64, ypxl::Float64, yref::Float64)
    zooming_in = (inorout == 0)
    if !isa(Plots.backend(),Plots.PyPlotBackend)
        printwarning("zooming requires PyPlot backend")
        return nothing
    end
    fig = gs.mainph.o
    ax_plt = fig.axes[1]
    inv = ax_plt.transData.inverted()
    xy = inv.transform((xpxl,yref+1-ypxl))
    # println("zooming "*ifelse(zooming_in,"in","out")*
    #         " at $xy for $yref $ypxl")
    z = xy[1]+1im*xy[2]
    ax = PseudospectraPlots.getxylims(gs.mainph)
    if zooming_in
        res = zoomin!(ps_data,z,ax)
    else
        res = zoomout!(ps_data,z,ax,include_fov=get(guiopts,:showfov,false))
    end
    if res < 0
        printwarning("unable to zoom based on requested point")
    elseif res == 0
        # redraw existing portrait
        printinfo("redrawing...")
        Pseudospectra.redrawcontour(gs, ps_data, guiopts)
        printinfo("done")
    else
        # compute and draw new portrait
        printinfo("recomputing...")
        driver!(ps_data, guiopts, gs)
        printinfo("done")
    end
    setedittext(ps_data.zoom_list[ps_data.zoom_pos])
    nothing
end

qmlfunction("loadmtx", loadmtx)
qmlfunction("savedata", savedata)
qmlfunction("initmsg", initmsg)
qmlfunction("go", go)
qmlfunction("use_eigs", use_eigs)
qmlfunction("setoptbykey", setoptbykey)
qmlfunction("setoptbydlg", setoptbydlg)
qmlfunction("setselectedz", setselectedz)
qmlfunction("monitor", monitor)
qmlfunction("mainzoom", mainzoom)

#=
# This looks clever but deadlocks
function zdlg()
    @emit getZ("Select z for pseudomode")
end

function guizdlg(gs; what="a z-value")
    global dialogcv
    global zselected
    t = Task(zdlg)
    schedule(t)
    wait(dialogcv)
    println("dialog selection: $zselected")
    zselected
end
=#
function getguiz(gs; what="a z-value")
    zselected
end

function abscissacb(epsilon)
    f,z = PSA.psa_abscissa(ps_data.input_matrix,epsilon)
    print("abscissa(ϵ = $epsilon): $f at ",z)
    msg = "abscissa(ϵ = $epsilon) = $f at $(z[1])"
    printinfo(msg)
    ps_data.ps_dict[:zpsabscissa] = z
    @emit saveVariable("Save point(s) at pseudospectral abscissa","zpsabscissa")
end

function radiuscb(epsilon)
    f,z = PSA.psa_radius(ps_data.input_matrix,epsilon)
    print("radius(ϵ = $epsilon): $f at ",z)
    msg = "radius(ϵ = $epsilon) = $f at $(z[1])"
    printinfo(msg)
    ps_data.ps_dict[:zpsradius] = z
    @emit saveVariable("Save point(s) at pseudospectral radius","zpsradius")
end

function ploteigv()
    if isnan(zselected)
        printwarning("Select z first.")
    else
        # TODO: wrap in gui stupefy/waken
        PSA.modeplot(ps_data,gs,0,zgetter=getguiz)
        @emit showPlot2(Int32(1))
    end
end
function plotpmode()
    if isnan(zselected)
        printwarning("Select z first.")
    else
        # TODO: wrap in gui stupefy/waken
        PSA.modeplot(ps_data,gs,1,zgetter=getguiz)
        @emit showPlot2(Int32(1))
    end
end

function mtxpwrcb(nmax)
    # This breaks easily, so we wrap it up in cotton wool.

    # FIXME: this is really ugly. Surely there is a better way?
    flag = false
    local JE
    try
        PSA.mtxpowersplot(gs,ps_data,nmax,gradual=false)
    catch JE
        flag = true
    end
    if flag
        printwarning("Matrix power comp/plot failed. See console.")
        rethrow(JE)
    else
        @emit showPlot2(Int32(1))
    end
end

function mtxexpcb(dt,nmax)
    flag = false
    local JE
    try
        PSA.mtxexpsplot(gs,ps_data,dt,nmax, gradual=false)
    catch JE
        flag = true
    end
    if flag
        printwarning("Matrix exp comp/plot failed. See console.")
        rethrow(JE)
    else
        @emit showPlot2(Int32(1))
    end
end

function plotps3d()
    PSA.surfplot(gs, ps_data, guiopts)
    @emit showPlot2(Int32(1))
end

qmlfunction("ploteigv", ploteigv)
qmlfunction("plotpmode", plotpmode)
qmlfunction("mtxpwrcb", mtxpwrcb)
qmlfunction("mtxexpcb", mtxexpcb)
qmlfunction("plotps3d", plotps3d)
qmlfunction("smartlevels", smartlevels)

################################################################
# functions to be invoked from Julialand

# drawing command, invoked as a callback inside PSA
function drawcmd(gs::PlotsGUIState,ph,id)
    global mainpltobj
    global pltobj2
    if id == 1
        mainpltobj = ph
        if maindisplay != nothing
            display(maindisplay, ph)
        end
        if running
            myproduce(compute_chnl,1)
            sleep(0.01)
            myproduce(compute_chnl,1)
        end
    else
        pltobj2 = ph
        if otherdisplay != nothing
            display(otherdisplay, ph)
        end
    end
end

"""
format a pair of reals so that they are distinguishable (if distinct)
and short (if possible)
"""
function formatpair(x1::Real,x2::Real)
    dx = abs(x2-x1)
    if !isfinite(dx) || (dx == 0) # punt
        s1 = @sprintf("%g",x1)
        s2 = @sprintf("%g",x2)
    else
        xscale = max(abs(x1),abs(x2))
        ndig = 2+ceil(Int,log10(xscale/dx))
        # this gave world-age problems (revisit if interpolation yields bad fmt)
        # "%." * @sprintf("%dg",ndig)
        fmt = "%.$(ndig)g"
        s1 = sprintf1(fmt,x1)
        s2 = sprintf1(fmt,x2)
    end
    return [s1,s2]
end

# This is a remedy for the annoying behavior of
# Qt toString(): 1/10 -> 0.100000000001 etc.
formataxislims(ax) = vcat(formatpair(ax[1],ax[2]), formatpair(ax[3],ax[4]))

# fills in the main text fields
function setedittext(zoom)
    ax = zoom.ax
    if !isempty(ax)
        @emit setAxisLims(formataxislims(ax)...)
    else
        @emit setAxisLims("NaN","NaN","NaN","NaN")
    end
    if zoom.levels.isunif
        @emit setLevels(zoom.levels.first,zoom.levels.last,zoom.levels.step)
    else
        # println("Not setting levels fields; zoom is")
        # println(zoom)
        @emit setLevels(NaN,NaN,NaN)
    end
    @emit setGrid(Int32(zoom.npts))
    if dvi != dvi_unset
        @emit toggleDirect(dvi == dvi_direct ? Int32(1) : Int32(0))
    end
    global need_recomp
    global need_redraw
    need_recomp = false
    need_redraw = false
end

# Warning: must be coordinated w/ ListModel in QML file
const which2idx = Dict{Symbol,Int}(:LM => 0, :LR => 1, :SR => 2,
                                   :LI => 3, :SI => 4)
function setarpackgui(ao)
    @emit setArpackOpts(Int32(ao.nev),
                        Int32(ao.ncv),
                        Int32(ao.maxiter),
                        ao.tol,
                        Int32(which2idx[ao.which]))
end

"""
place some text in the main App banner
"""
function printbanner(msg::AbstractString)
    @emit setMainMsg(msg,Int32(0))
end

"""
place some text into the informational text area
"""
function printinfo(msg::AbstractString)
    @emit setMainMsg(msg,Int32(1))
    if running
        myproduce(compute_chnl,1)
        sleep(0.01)
        myproduce(compute_chnl,1)
    end
end

"""
place some text into the informational text area
"""
function printwarning(msg::AbstractString)
    @emit setMainMsg(msg,Int32(2))
    if running
        myproduce(compute_chnl,1)
        sleep(0.01)
        myproduce(compute_chnl,1)
    end
end

"""
called after a new matrix is available in ps_data
"""
function refresh()
    global mainplotobj,plotobj2,need_recomp
    mainplotobj,plotobj2 = nothing,nothing
    # Fixme: clear canvasses
    m,n = size(ps_data.input_matrix)
    zoom = ps_data.zoom_list[ps_data.zoom_pos]
    ps_dict = ps_data.ps_dict
    ax = zoom.ax
    if ps_dict[:direct]
        # FIXME: need a function to check whether ew's will be available
        if ((issparse(ps_data.input_matrix)
             && get(ps_dict,:sparse_direct,false))
            || (m != n)
            || length(methods(schur!,(typeof(ps_data.input_matrix),))) == 0
            )
            # rectangular, sparse-direct, or just plain weird
            @emit toggleAutoAx(Int32(-1))
            if isempty(ax)
                # flag need to get ax
                printinfo("Axis limits must be specified to proceed")
                # FIXME: should wait until they are ready
                @emit enableGo(Int32(0))
                need_recomp = true
            end
        else
            # should be ready for this
            @emit toggleDirect(Int32(1))
            # FIXME: wrap this in QTimer-monitored task
            if isempty(ax) && !isempty(get(ps_dict,:ews,[]))
                PSA.origplot!(ps_data,guiopts,gs)
            else
                PSA.driver!(ps_data,guiopts,gs)
            end
            @emit showPlot(Int32(1))
            @emit enableGo(Int32(0))
        end
    else # iterative
        @emit toggleDirect(Int32(0))
        if !haskey(ps_dict,:arpack_opts)
            # FIXME: warn via gui and bail
            throw(ArgumentError("ps_data has no arpack opts"))
        end
        setarpackgui(ps_dict[:arpack_opts])
        # If the user specified (some) ARPACK options, we might be ready to go.
        # Nah: esp. w/o a stop button, premature start is more likely and
        # more troublesome.
        #=
        if haskey(guiopts,:arpack_opts)
            # FIXME: couldn't clear, so may get junk while starting up
            @emit showPlot(Int32(1))
            PSA.driver!(ps_data,guiopts,gs)
            @emit enableGo(Int32(0))
        else
        =#
        begin
            printinfo("Set/check ARPACK options before proceeding")
            need_recomp = true
            # CHECKME: should we wait until options are ready?
            @emit enableGo(Int32(1))
        end
    end
    setedittext(ps_data.zoom_list[ps_data.zoom_pos])
end


################################################################
# This is the main App execution sequence.
# QML seems to require that key data structures be in global scope.

gs = PlotsGUIState(nothing,0,drawcmd)

# Load up and check data first so if PSA throws an exception
# we are not inchoate.

if !isempty(PSAData.getdefaultmtx())
    (verbosity > 0) && println("Using existing matrix from PSAData")
    # guiopts is empty until now
    merge!(guiopts,PSA.fillopts(gs,PSAData.opts))
    ps_data = PSA.new_matrix(PSAData.getdefaultmtx(),guiopts)
    # forget opts that should not apply to newly ingested data
    # viz. when new_matrix() is invoked from GUI.
    PSAData.clearopts(false)
end

# Now start up the engine

qml_file = joinpath(dirname(@__FILE__), "psagui.qml")

load(qml_file, timer = QTimer())

if !isempty(PSAData.getdefaultmtx())
    refresh()
end

# Run the application
QML.exec()

end
