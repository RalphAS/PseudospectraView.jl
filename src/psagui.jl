module PSApp

# WARNING: some comments are obsolete!
# Even worse, some control flow is probably obsolete too.

function cprint(s::String)
    n = length(s)+1
    fid = 2
    ccall(:write,Clong,(Cint,Cstring,Clong),fid,s,n)
end

# needed for use of JuliaPaintedItem
ENV["QSG_RENDER_LOOP"] = "basic"

using QML
using Observables

using ..PseudospectraQML: _verbosity, _pbackend
verbosity() = _verbosity[]

# the App will store return values here:
const results = Dict{Symbol,Any}()

# which entry in QML stacks to use for graphs
live_pane() = _pbackend[] == :pyplot ? 2 : 1

using Plots

using CxxWrap

# Other backends are not yet really functional, so try to force this for now.
if _pbackend[] == :pyplot
    # Python must not manage a gui
    ENV["MPLBACKEND"] = "Agg"
    pyplot()
else
    # Experimental!
    # @warn "GR backend is experimental!"
    ENV["GKSwstype"] = "use_default"
    ENV["GKS_WSTYPE"] = 381
    ENV["GKS_QT_VERSION"] = 5
    gr(show=false)
end

(verbosity() > 1) && println("Plots loaded, backend is $(typeof(Plots.backend())).")

using Pseudospectra
# shorthand
PSA = Pseudospectra

PseudospectraPlots = Pseudospectra.PseudospectraPlots
PlotsGUIState = Pseudospectra.PseudospectraPlots.PlotsGUIState
function drawcmd end # forward decl
const gs = PlotsGUIState(nothing,0,drawcmd)


using PseudospectraQML

PSAData = PseudospectraQML.PSAData

using Formatting

using Printf: @sprintf
using SparseArrays: issparse
# for schur!
using LinearAlgebra


################################################################
# Global data
# Since the App runs in global scope, and dependencies make this godawful
# slow anyway, we may as well use global vars where it's convenient.
# A few more are created by big function calls below.

maindisplay = nothing
mainpltobj = nothing
const mainsize = [0,0]

otherdisplay = nothing
pltobj2 = nothing
const othersize = [0,0]

# the instance of PSA opts maintained by the GUI
guiopts = Dict{Symbol,Any}()

ps_data = nothing

firstcall = true

@enum DirVsIter dvi_unset dvi_direct dvi_iter
dvi = dvi_unset

const timer = QTimer()

const cpaint = Channel{Int}(1)
const mcounter = Ref(0)

# flags indicating whether items changed by GUI require action
need_recomp = false
need_redraw = false

# In order to force updates of the GUI window, we use tasks
# and a QTimer. For now we use produce/consume logic.
compute_chnl = nothing
const running = Ref(false)

const axname2idx = Dict(:xMin=>1,:xMax=>2,:yMin=>3,:yMax=>4)

tmpdata = nothing

# Until we can get a plot location from the pointer and GUI canvas,
# we use a dialog which saves its result here:
zselected = NaN+0.0im

################################################################
# Global state
# Observables for data exchange between GUI and Julia

const o_banner_msg = Observable("Welcome to PseudospectraQML")
const o_info_msg = Observable("")
const computation = Observable{Any}(nothing)
const computation_key = Observable(Cint(0))
const o_progress = Observable(0.0)
const o_gridpts = Observable("")
const o_xmin = Observable("")
const o_xmax = Observable("")
const o_ymin = Observable("")
const o_ymax = Observable("")
const o_firstlev = Observable("")
const o_lastlev = Observable("")
const o_levstep = Observable("")
const o_fov = Observable(Cint(0))
const o_imaxis = Observable(Cint(0))
const o_unitcircle = Observable(Cint(0))
const o_arpack_nev = Observable("6")
const o_arpack_ncv = Observable("8")
const o_arpack_which = Observable("LM")
const o_ticks = Observable(Cint(0)) # Number of times the timer has ticked
const o_surf_ok = Observable(_pbackend[] == :pyplot ? Cint(1) : Cint(0))

propmap = JuliaPropertyMap(
    "computationkey" => computation_key,
    "gridpts" => o_gridpts,
    "xmin" => o_xmin,
    "xmax" => o_xmax,
    "ymin" => o_ymin,
    "ymax" => o_ymax,
    "ticks" => o_ticks,
    "banner_msg" => o_banner_msg,
    "info_msg" => o_info_msg,
    "firstlev" => o_firstlev,
    "lastlev" => o_lastlev,
    "levstep" => o_levstep,
    "fov" => o_fov,
    "imag_axis" => o_imaxis,
    "unit_circle" => o_unitcircle,
    "arpack_nev" => o_arpack_nev,
    "arpack_ncv" => o_arpack_ncv,
    "arpack_which" => o_arpack_which,
    "surface_ok" => o_surf_ok,
)

mutable struct Monitor
    status::Int # 0:initialized 1:running 2:interrupted 3:done 4:error
    key::Int32
    function Monitor()
        me = new(0,computation_key[])
        on(computation_key) do key
            me.key = key
            me.status = 0
        end
        return me
    end
end
isfinished(c::Monitor) = c.status >= 3
struct Runner
    channel::Channel
    ready::Bool
end
function compute_chan(channel::Channel)
    verbosity() > 1 && println("starting computation")
    c = Monitor()
    while !isfinished(c)
        run!(c)
        if isready(channel)
            cprint("blocking in cc\n")
        end
        cprint("putting $(c.status) in cc\n")
        put!(channel, c.status)
        cprint("continuing in cc\n")
        # if c.status == 3
        #    QML.stop(timer)
        # end
    end
    cprint("cc done\n")
    put!(channel, 4)
    cprint("cc exiting\n")
    nothing
end

const pvals = [0.0, 0.5, 0.99, 1.0, -1.0]
function step(r::Runner)
    # failed attempt to avoid deadlock by letting another tick go by
  # if isopen(r.channel)
    # NOTE: if we don't block on r.channel, everything waits in limbo
    if !isready(r.channel)
        cprint("blocking in step\n")
    end
    s = take!(r.channel)
    cprint("taking $s in step\n")
    verbosity() > 1 && println("step got $s")
#    if s == 3
#        QML.stop(timer)
#    end
    o_progress[] = pvals[s]
    @emit updateMainPlot()
  # end
end

on(o_ticks) do t
    cprint("tick\n")
    if t > 0
        step(computation[])
    end
#    if !isready(cpaint)
#        verbosity() > 1 && println("ticker puts to cpaint")
#        put!(cpaint, -1)
#    end
end

function run!(c::Monitor)
    c.status = 1
    running[] = true
    (verbosity() > 1) && println("calculating...")
    printinfo("calculating...")
    PSA.driver!(ps_data,guiopts,gs,myprintln=printinfo)
    setedittext(ps_data.zoom_list[ps_data.zoom_pos])
    running[] = false
    (verbosity() > 1) && println("calculation done.")
    printinfo("done.")
    c.status = 3
end

function setup_computation(key)
    computation_key[] = key
    o_progress[] = 0.0
    computation[] = Runner(Channel(compute_chan), false)
end

on(o_progress) do p
    if p >= 1
        # this acts as finalizer for Runner
        println("ticks: ", o_ticks[])
        QML.stop(timer)
        o_ticks[] = Cint(0)
        computation[] = nothing
    end
end

on(o_gridpts) do s
    local n
    n = tryparse(Int, s)
    if n === nothing
        @warn "Invalid entry for npts"
        return
    end
    global need_recomp
    global need_redraw
    if ! (ps_data === nothing)
        curzoom = ps_data.zoom_list[ps_data.zoom_pos]
        old_n = curzoom.npts
        if n >= 5
            if n != old_n
                curzoom.npts = n
                need_recomp = true
                need_redraw = true
                @emit enableGo(Int32(1))
            end
        else
            # this should be obviated by QML logic.
            @warn "npts must be at least 5; got $n"
            o_gridpts[] = max(curzoom.npts, 5)
        end
    end
end

function update_ax(k, s::AbstractString)
    global need_recomp
    global need_redraw
    local x
    x = tryparse(Float64,s)
    if x === nothing
        @warn "Invalid entry for $k"
        return
    end
    (ps_data === nothing) && return
    curzoom = ps_data.zoom_list[ps_data.zoom_pos]
    isempty(curzoom.ax) && (curzoom.ax = fill(NaN,(4,)))
    curzoom.ax[axname2idx[k]] = x
    need_recomp = true
    if !any(isnan.(curzoom.ax)) # goable
        @emit enableGo(Int32(1))
    end
end

on(o_xmin) do s
    update_ax(:xMin, s)
end
on(o_xmax) do s
    update_ax(:xMax, s)
end
on(o_ymin) do s
    update_ax(:yMin, s)
end
on(o_ymax) do s
    update_ax(:yMax, s)
end

function update_levels(k, s)
    global need_redraw
    x = tryparse(Float64, s)
    if x === nothing
        @warn "Invalid entry for contour levels $k"
        return
    end
    curzoom = ps_data.zoom_list[ps_data.zoom_pos]
    setfield!(curzoom.levels, k, x)
    curzoom.autolev = false
    need_redraw = true
    @emit enableGo(Int32(1))
end
on(o_firstlev) do s
    update_levels(:first, s)
end
on(o_lastlev) do s
    update_levels(:last, s)
end
on(o_levstep) do s
    update_levels(:step, s)
end

on(o_arpack_nev) do s
    global need_recomp
    ps_dict = ps_data.ps_dict
    if haskey(ps_dict,:arpack_opts)
        n = tryparse(Int, s)
        if n === nothing
            @warn "Invalid entry for ARPACK nev"
            return
        end
        if (n < 1) || (n > size(ps_data.input_matrix,1)-2)
            @warn "Invalid entry for ARPACK nev"
        end
        if ps_dict[:arpack_opts].nev != n
            ps_dict[:arpack_opts].nev = n
            ps_dict[:proj_valid] = false
            need_recomp = true
            @emit enableGo(Int32(1))
        end
    end
end

on(o_arpack_ncv) do s
    global need_recomp
    ps_dict = ps_data.ps_dict
    if haskey(ps_dict,:arpack_opts)
        n = tryparse(Int, s)
        if n === nothing
            @warn "Invalid entry for ARPACK ncv"
            return
        end
        if (n < 1) || (n > size(ps_data.input_matrix,1)-2)
            @warn "Invalid entry for ARPACK ncv"
        end
        if ps_dict[:arpack_opts].ncv != n
            ps_dict[:arpack_opts].ncv = n
            ps_dict[:proj_valid] = false
            need_recomp = true
            @emit enableGo(Int32(1))
        end
    end
end

on(o_arpack_which) do s
    (verbosity() > 1) && println("attempting to set which to $s")
    global need_recomp
    ps_dict = ps_data.ps_dict
    if s ∈ ["LM", "SM", "LR", "SR", "LI", "SI"]
        if haskey(ps_dict,:arpack_opts)
            if ps_dict[:arpack_opts].which != Symbol(s)
                ps_dict[:arpack_opts].which = Symbol(s)
                need_recomp = true
                @emit enableGo(Int32(1))
            end
        end
    else
        @warn "Invalid entry for ARPACK \"which\""
    end
end

on(o_fov) do s
    tf = Bool(s)
    global need_redraw
    if get(guiopts,:showfov,false) != tf
        need_redraw = true
        @emit enableGo(Int32(1))
    end
    guiopts[:showfov] = tf
end

on(o_imaxis) do s
    tf = Bool(s)
    global need_redraw
    if get(guiopts,:showimagax,false) != tf
        need_redraw = true
        @emit enableGo(Int32(1))
    end
    guiopts[:showimagax] = tf
end

on(o_unitcircle) do s
    tf = Bool(s)
    global need_redraw
    if get(guiopts,:showunitcircle,false) != tf
        need_redraw = true
        @emit enableGo(Int32(1))
    end
    guiopts[:showunitcircle] = tf
end

################################################################
# QMLFunctions: how the GUI talks to Julia
# Apparently qmlfunctions need to be defined in global scope
# and QML.exec() needs to run in global scope to see them.

# These are called by the App code, apparently only on setup or resize
# Display is generally forced via drawcmd().

# somebody doesn't manage dimensions correctly
const pyplot_fudge = 0.9

function init_backend(width::Float64, height::Float64, idx)
    if width < 5 || height < 5
        return
    end
    if idx == 1
        mainsize .= round.(Int,[width,height])
    else
        othersize .= round.(Int,[width,height])
    end
    if isa(Plots.backend(),Plots.PyPlotBackend)
        pyplot(size=(round(Int64,pyplot_fudge*width),
                     round(Int64,height)))
    elseif isa(Plots.backend(),Plots.GRBackend)
        @warn "Support for GR backend is incomplete; expect trouble." maxlog=1
        # CHECKME: GR may need a vertical fudge factor
        gr(size=(round(Int64,width),round(Int64,height)))
        Plots.GR.inline()
    else
        @error "unsupported Plots backend"
    end
    return
end

# This is supposed to help with z-picking. It doesn't, yet.
function set_paint_size(width::Float64, height::Float64, idx)
    if width < 5 || height < 5
        return
    end
    if idx == 1
        mainsize .= round.(Int,[width,height])
    else
        othersize .= round.(Int,[width,height])
    end
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
qmlfunction("size_paint", set_paint_size)

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
        n = tryparse(Int,val)
        if n === nothing
            @warn "Invalid entry for ARPACK maxiter"
            return
        end
        if ao.maxiter != n
            ao.maxiter = n
            need_recomp = true
            ao_current = false
        end
    elseif k == :ArpackNcv
        n = tryparse(Int,val)
        if n === nothing
            @warn "Invalid entry for ARPACK ncv"
            return
        end
        if ao.ncv != n
            ao.ncv = n
            need_recomp = true
            ao_current = false
        end
    elseif k == :ArpackTol
        x = tryparse(Float64,val)
        if x === nothing
            @warn "Invalid entry for ARPACK tol"
            return
        end
        if !isapprox(ao.tol,x)
            ao.tol = x
            need_recomp = true
            ao_current = false
        end
    elseif k == :AbscissaEps
        x = tryparse(Float64,val)
        if x ===  nothing
            @warn "Invalid entry for ϵ"
            return
        end
        if x < 0
            @warn "Invalid entry for ϵ"
        end
        abscissacb(x)
    elseif k == :RadiusEps
        x = tryparse(Float64,val)
        if x === nothing
            @warn "Invalid entry for ϵ"
            return
        end
        if x < 0
            @warn "Invalid entry for ϵ"
        end
        radiuscb(x)
    elseif k == :ProjLev
        x = tryparse(Float64,val)
        if x === nothing
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

function go()
    global need_recomp
    global need_redraw
    if ps_data == nothing
        @warn "go invoked before population"
        return
    end
    @emit showPlot(Int32(0), Int32(live_pane()))
    if need_recomp
        (verbosity() > 1) && println("opening channel for computation task")
        @emit startTimer(Int32(10))
        setup_computation(0)
        need_recomp = false
        need_redraw = false
    elseif need_redraw
        printinfo("redrawing")
        (verbosity() > 1) && println("redrawing")
        Pseudospectra.redrawcontour(gs,ps_data,guiopts)
        need_redraw = false
    end
#    ax = (mainpltobj != nothing) ? PSA.getxylims(mainpltobj) : zeros(0)
    setedittext(ps_data.zoom_list[ps_data.zoom_pos])
    printinfo("ready")
    @emit enableGo(Int32(0))
end


# assorted callbacks

function setselectedz(zrs::AbstractString,zis::AbstractString)
    global zselected
    zr,zi = NaN,NaN
    zr = tryparse(Float64,zrs)
    zi = tryparse(Float64,zis)
    if zr === nothing || zi === nothing
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
    @emit showPlot(Int32(0), Int32(0))
    @emit showPlot(Int32(1), Int32(0))

    # this should really be done by the GUI, but this is easier:
    if isempty(varname)
        varname = "A"
    end

    global firstcall
    if firstcall
        printinfo("Please be patient, plotting code needs to be compiled...")
    end

    try
        Core.eval(PSAData, Meta.parse("$(varname) = $(myexpr)"))
        Core.eval(PseudospectraQML.PSApp,
             Meta.parse("ps_data = Pseudospectra.new_matrix(PSAData.$(varname),guiopts)"))
        (verbosity() > 1) && println("Matrix has been digested.")
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
    (verbosity() > 1) && println("matrix loaded")
    nothing
end

function savedata(varname::AbstractString,mykey)
    # FIXME: test for valid variable name
    if mykey == "portrait"
        tmpdata = ps_data.zoom_list[ps_data.zoom_pos]
    elseif mykey == "ps_data"
        tmpdata = ps_data
    elseif mykey == "zpsradius"
        tmpdata = ps_data.ps_dict[:zpsradius]
    elseif mykey == "zpsabscissa"
        tmpdata = ps_data.ps_dict[:zpsabscissa]
    else
        @warn "unrecognized save key $mykey"
        return
    end
    try
        results[Symbol(varname)] = deepcopy(tmpdata)
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
            # CHECKME: what's the reasoning here?
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

function zmoused(leftright::Int32, xpxl::Float64, ypxl::Float64, yref::Float64)
    if isa(Plots.backend(),Plots.PyPlotBackend)
        fig = gs.mainph.o
        ax_plt = fig.axes[1]
        inv = ax_plt.transData.inverted()
        xscaled = xpxl * pyplot_fudge
        xy = inv.transform((xscaled,yref+1-ypxl))
    elseif isa(Plots.backend(),Plots.GRBackend)
        fig = Plots.GR.gcf()
        nc, nr = mainsize # fig[:size]
        # println("GR thinks size is $nc $nr")
        fac = (nc > nr) ? (1.0 / nc) : (1.0 / nr)
        xndc = fac * xpxl
        yndc = fac * (yref+1-ypxl)
        xy = Plots.GR.ndctowc(xndc, yndc)
    else
        printwarning("zooming or point selection not implemented for $(typeof(Plots.backend()))")
        return nothing
    end
    (verbosity() > 1) && println("selected $xy from yvals $yref $ypxl")
    z = xy[1]+1im*xy[2]
    if true
        zooming_in = (leftright == 0)
        mainzoom(zooming_in, z)
    else
    end
end

function mainzoom(zooming_in::Bool, z)
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
qmlfunction("go", go)
qmlfunction("use_eigs", use_eigs)
qmlfunction("setoptbydlg", setoptbydlg)
qmlfunction("setselectedz", setselectedz)
# qmlfunction("monitor", monitor)
qmlfunction("zmoused", zmoused)

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
    println("abscissa(ϵ = $epsilon): $f at ",z)
    msg = "abscissa(ϵ = $epsilon) = $f at $(z[1])"
    printinfo(msg)
    ps_data.ps_dict[:zpsabscissa] = (ϵ = epsilon, x = f, z = z)
    @emit saveVariable("Save point(s) at pseudospectral abscissa","zpsabscissa")
end

function radiuscb(epsilon)
    f,z = PSA.psa_radius(ps_data.input_matrix,epsilon)
    println("radius(ϵ = $epsilon): $f at ",z)
    msg = "radius(ϵ = $epsilon) = $f at $(z[1])"
    printinfo(msg)
    ps_data.ps_dict[:zpsradius] = (ϵ = epsilon, r = f, z = z)
    @emit saveVariable("Save point(s) at pseudospectral radius","zpsradius")
end

function ploteigv()
    if isnan(zselected)
        printwarning("Select z first.")
    else
        # TODO: wrap in gui stupefy/waken
        PSA.modeplot(ps_data, 0; gs=gs, zgetter=getguiz)
        @emit showPlot(Int32(1), Int32(live_pane()))
    end
end
function plotpmode()
    if isnan(zselected)
        printwarning("Select z first.")
    else
        # TODO: wrap in gui stupefy/waken
        PSA.modeplot(ps_data, 1; gs=gs, zgetter=getguiz)
        @emit showPlot(Int32(1), Int32(live_pane()))
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
        @emit showPlot(Int32(1), Int32(live_pane()))
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
        @emit showPlot(Int32(1), Int32(live_pane()))
    end
end

function plotps3d()
    PSA.surfplot(gs, ps_data, guiopts)
    @emit showPlot(Int32(1), Int32(live_pane()))
end

qmlfunction("ploteigv", ploteigv)
qmlfunction("plotpmode", plotpmode)
qmlfunction("mtxpwrcb", mtxpwrcb)
qmlfunction("mtxexpcb", mtxexpcb)
qmlfunction("plotps3d", plotps3d)
qmlfunction("smartlevels", smartlevels)

# these will be wrapped below
function paint_main(p::CxxPtr{QPainter}, item::CxxPtr{JuliaPaintedItem})
    (_pbackend[] == :gr) || return
    ENV["GKS_CONID"] = split(repr(p.cpp_object), "@")[2]
    dev = device(p[])[]
    r = effectiveDevicePixelRatio(window(item[])[])
    w, h = width(dev) / r, height(dev) / r
    if running[]
        c = mcounter[] + 1
        mcounter[] = c
#        if !isready(cpaint)
#            verbosity() > 1 && println("painter puts to cpaint")
#            put!(cpaint, c)
#       end
        # @emit stopTimer()
    end
    if (gs.mainph === nothing)
        (verbosity() > 1) && println("null paint_main called")
        return
    end
    (verbosity() > 2) && println("paint_main called")
    if mainsize[1] > 0
        plot!(gs.mainph, size=(w,h)) # Tuple(mainsize))
        display(gs.mainph)
    end
    return
end
function paint_second(p::CxxPtr{QPainter}, item::CxxPtr{JuliaPaintedItem})
    (_pbackend[] == :gr) || return
    ENV["GKS_CONID"] = split(repr(p.cpp_object), "@")[2]
    dev = device(p[])[]
    r = effectiveDevicePixelRatio(window(item[])[])
    w, h = width(dev) / r, height(dev) / r
    if (gs.ph2 === nothing)
        (verbosity() > 1) && println("null paint_second called")
        return
    end
    (verbosity() > 1) && println("paint_second called")
    if othersize[1] > 0
        plot!(gs.ph2, size=(w,h))
        #@warn "suppressing display"
        px = plot(gs.ph2[1],layout=(1,1))
        display(px)
    end
    return
end

################################################################
# functions to be invoked from Julialand

# drawing command, invoked as a callback inside PSA
function drawcmd(gs::PlotsGUIState,ph,id)
    global mainpltobj
    global pltobj2
    if id == 1
        if _pbackend[] == :pyplot
            mainpltobj = ph
            if maindisplay != nothing
                display(maindisplay, ph)
            end
        else
            if running[]
                if computation[] === nothing
                    cprint("? running w/o computation\n")
                    # @emit dieDieDie()
                else
                    if isready(computation[].channel)
                        cprint("blocking in drawcmd\n")
                    end
                    cprint("putting 2 in drawcmd\n")
                    put!(computation[].channel, 2)
                    cprint("continuing in drawcmd\n")
                end
            else
                @emit updateMainPlot()
            end
        end
    else
        if _pbackend[] == :pyplot
            pltobj2 = ph
            if otherdisplay != nothing
                display(otherdisplay, ph)
            end
        else
            gs.ph2 = ph
            @emit updatePlot2()
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
    o_banner_msg[] = msg
end

"""
place some text into the informational text area
"""
function printinfo(msg::AbstractString)
    o_info_msg[] = msg
    # TODO: set an info color or icon
end

"""
place a warning into the informational text area
"""
function printwarning(msg::AbstractString)
    o_info_msg[] = msg
    # TODO: set a warning color or icon
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
            @emit showPlot(Int32(0), Int32(live_pane()))
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
            @emit showPlot(Int32(0), Int32(1))
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

# Load up and check data first so if PSA throws an exception
# we are not inchoate.

if !isempty(PSAData.getdefaultmtx())
    (verbosity() > 0) && println("Using existing matrix from PSAData")
    # guiopts is empty until now
    merge!(guiopts,PSA.fillopts(gs,PSAData.opts))
    ps_data = PSA.new_matrix(PSAData.getdefaultmtx(),guiopts)
    (verbosity() > 1) && println("Matrix has been digested.")
    # forget opts that should not apply to newly ingested data
    # viz. when new_matrix() is invoked from GUI.
    PSAData.clearopts(false)
end

# Now start up the engine

qml_file = joinpath(dirname(@__FILE__), "psagui.qml")

load(qml_file, timer=timer, observables=propmap,
     paint_main_wrapped = @safe_cfunction(paint_main, Cvoid,
                                          (CxxPtr{QPainter}, CxxPtr{JuliaPaintedItem})),
     paint_second_wrapped = @safe_cfunction(paint_second, Cvoid,
                                          (CxxPtr{QPainter}, CxxPtr{JuliaPaintedItem}))
     )

if !isempty(PSAData.getdefaultmtx())
    refresh()
end

function runme()
    @emit setInterval(Cint(100))
    # Run the application
    QML.exec()
    return results
end

end
