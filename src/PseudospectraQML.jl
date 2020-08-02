module PseudospectraQML

using Plots
using Pseudospectra

export psagui

# CHECKME: why do QML examples use dirname(@__FILE__) instead?
const pvsrcdir = realpath(@__DIR__)

function psagui()
    if !isdefined(PseudospectraQML,:PSAData)
        include(joinpath(pvsrcdir,"psadata.jl"))
    end
    include(joinpath(pvsrcdir,"psagui.jl"))
end

"""
    psagui()
    psagui(A[,opts])

Start the QML GUI for Pseudospectra, optionally loading matrix `A` first.
"""
function psagui(A::AbstractMatrix,opts::Dict{Symbol,Any}=Dict{Symbol,Any}())
    # Yes, wise reader, this function violates several rules.
    # Suggestions for effective alternatives are welcome.

    # !isdefined(:PSAData) is inadequate for full reload.
    if !isdefined(PseudospectraQML,:PSAData)
        include(joinpath(pvsrcdir,"psadata.jl"))
    end
    Base.invokelatest(PSAData.setdefaultmatrix,A)
    isempty(opts) || Base.invokelatest(PSAData.addopts,opts)
    include(joinpath(pvsrcdir,"psagui.jl"))
    results = Base.invokelatest(PSApp.runme)
    return (length(keys(results)) == 0) ? nothing : results
end

end # module
