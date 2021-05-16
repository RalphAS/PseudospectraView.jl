module PseudospectraView

using Plots
using Pseudospectra
using QML
using Qt5QuickControls_jll

export psagui

# CHECKME: why do QML examples use dirname(@__FILE__) instead?
const pvsrcdir = realpath(@__DIR__)

# this is for tracking GUI-related activity, distinct from PSA option.
# 0: quiet, 1: reassure us that something is going on, 2: help the developer
const _verbosity = Ref(0)

# backend for Plots
const _pbackend = Ref(:unset)
const _supported_backends = [:pyplot, :gr]

include("psadata.jl")

function psagui(; backend::Symbol=:pyplot)
    (backend in _supported_backends) || throw(ArgumentError("unsupported Plots backend"))
    _pbackend[] = backend
    include(joinpath(pvsrcdir,"psagui.jl"))
    results = Base.invokelatest(PSApp.runme)
    return (length(keys(results)) == 0) ? nothing : results
end

"""
    psagui()
    psagui(A[,opts])

Start the QML GUI for Pseudospectra, optionally loading matrix `A` first.
`opts` is a `Dict{Symbol,Any}` of Pseudospectra options (many of which may
be reset by GUI objects).
Returns a `Dict` with variables saved by the GUI, or `nothing`.
"""
function psagui(A::AbstractMatrix,
                opts::Dict{Symbol,Any}=Dict{Symbol,Any}();
                backend::Symbol=:pyplot)
    (backend in _supported_backends) || throw(ArgumentError("unsupported Plots backend"))
    _pbackend[] = backend
    Base.invokelatest(PSAData.setdefaultmatrix,A)
    isempty(opts) || Base.invokelatest(PSAData.addopts,opts)
    include(joinpath(pvsrcdir,"psagui.jl"))
    results = Base.invokelatest(PSApp.runme)
    return (length(keys(results)) == 0) ? nothing : results
end

end # module
