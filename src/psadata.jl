"""
This is a hack to work around scoping issues in QML.
We put user options and data in the global scope of a known module
so we can run the GUI code as a script after they are defined.
"""
module PSAData
import Pseudospectra

# allow shorthand
PSA = Pseudospectra

# The opts variable has two uses:
# (a) options for use w/ a particular matrix, when it is provided to the
#     GUI-invoking function
# (b) general preferences, e.g. graphics attributes
opts = Dict{Symbol,Any}()

Adefault = zeros(0,0)

# FIXME: need to analyze full catalog
const mtxspecific = [:ax,:direct,:levels,:sparse_direct,:proj_lev,:npts,
                     :real_matrix]

"""
    addopts(newopts::Dict{Symbol,Any})

store options to be used with cached matrices
"""
function addopts(newopts::Dict{Symbol,Any})
    global opts
    for k in keys(newopts)
        opts[k] = newopts[k]
    end
end
function clearopts(everything::Bool=false)
    global opts
    for k in keys(opts)
        if everything || (k ∈ mtxspecific)
            delete!(opts,k)
        end
    end
end
"""
    setdefaultmatrix(Anew::AbstractMatrix)

Cache a matrix to be analyzed by PseudospectraQML.psagui().
"""
function setdefaultmatrix(Anew::AbstractMatrix)
    global Adefault
    Adefault = Anew
end

function getdefaultmtx()
    Adefault
end

function getopts()
    opts
end

end # PSAData
