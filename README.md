# PseudospectraQML.jl

![Lifecycle](https://img.shields.io/badge/lifecycle-experimental-orange.svg)<!--
![Lifecycle](https://img.shields.io/badge/lifecycle-maturing-blue.svg)
![Lifecycle](https://img.shields.io/badge/lifecycle-stable-green.svg)
![Lifecycle](https://img.shields.io/badge/lifecycle-retired-orange.svg)
![Lifecycle](https://img.shields.io/badge/lifecycle-archived-red.svg)
![Lifecycle](https://img.shields.io/badge/lifecycle-dormant-blue.svg) -->
[![Build Status](https://travis-ci.org/RalphAS/PseudospectraQML.jl.svg?branch=master)](https://travis-ci.org/RalphAS/PseudospectraQML.jl)
[![codecov.io](http://codecov.io/github/RalphAS/PseudospectraQML.jl/coverage.svg?branch=master)](http://codecov.io/github/RalphAS/PseudospectraQML.jl?branch=master)
<!--
[![Documentation](https://img.shields.io/badge/docs-stable-blue.svg)](https://RalphAS.github.io/PseudospectraQML.jl/stable)
[![Documentation](https://img.shields.io/badge/docs-master-blue.svg)](https://RalphAS.github.io/PseudospectraQML.jl/dev)
-->

# Introduction
PseudospectraQML is a Julia/Qt GUI for analyzing nonsymmetric matrices by
displaying eigenvalues and pseudospectra, along with related properties.
It is built as a driver for [Pseudospectra.jl](https://github.com/RalphAS/Pseudospectra.jl).

# Installation

The current version of this package requires Plots.jl with the PyPlot backend.

The graphical user interface is based on
[QML.jl](https://github.com/barche/QML.jl). Install that and make sure it
works first.

Install [Pseudospectra.jl](https://github.com/RalphAS/Pseudospectra.jl); this
must be done manually, pending registration.

Finally, add this repo ("https://github.com/RalphAS/PseudospectraQML.jl")
with the package manager.

# Basics
Invoke the App with

```julia
using QML
using PseudospectraQML
A=your_matrix_generating_function()
psagui(A)
```

(Alternatively one can invoke `psagui()` with no arguments, and
run the `File/New Matrix` menu command.)

Please be patient; the infrastructure takes a while to initialize.
An app window should appear, and a spectral portrait (i.e. a contour
plot of base-10 logarithms of inverse resolvent norms) should be shown
in the main display.

The default mesh is very coarse; adjust the grid size and axis limits to
taste and press the "Update" button to recompute.

One can zoom in by clicking the left mouse button with the pointer in the
main plot area. (Zoom out with the right button.) Note that this usually
involves recomputation.

To display an eigenmode or pseudo-eigenmode, use the "Select z" button to
specify the value, then the appropriate mode button.

## FIXME: walk-throughs with screenshots needed here

Most GUI behavior mimics EigTool.
Much of [the EigTool documentation](http://www.cs.ox.ac.uk/pseudospectra/eigtool/documentation/index.html)
also applies to PseudospectraQML.

# Saving figures
To save a plot for external use, the recommended procedure is as follows:
1. save the current data as a `Portrait` type, using the "File" menu
2. exit the App (Alt-F4 works for me)
3. use fields of the saved variable to make a plot
E.g., if you saved the `Portrait` as variable `baz`,
```julia
using Plots
contour(baz.x,baz.y,baz.Z)
png("baz")
```

# Striking differences from EigTool
* One must use the "Select z" command button **before** plotting an
  eigenmode/pseudomode or transients.
* The "pause/resume/stop" logic is not implemented.

# Other notable peculiarities
* At present, there is a use dichotomy: in a given session, one may
  use **either** a direct interface to Pseudospectra.jl with interactive
  graphics **or** the GUI. (This is because the plotter must use a
  non-interactive backend to prepare images for QML, and coordination of
  backend switching with the Pseudospectra API is only a dream.)
* The App is in its own module, which is reloaded for each instance.
* PseudospectraQML uses a separate module (PSAData) to store
  some variables.

# Known issues
* The GUI is finicky.
  * Text fields may not be properly converted to associated variables. Using
    the Tab key to check in the value seems most reliable.
  * Coordinate conversion for zooming is flaky. I find that resizing the
    whole app window, then forcing a redraw, makes it more accurate.
  * Obscure callback problems will arise occasionally.
* Some warnings and errors provoked from the GUI are reported in the console.
* Some buttons & menu items remain enabled when incompatible with state.
* The control flow is very convoluted.
  * Zooming and the axis-setting widgets are not coordinated (yet).
  * Deadlocks or crashes are possible. File an issue if you encounter them.

# Caveat
The current translator/author is a novice at GUI construction, and generally
prefers command-line interfaces. Hence this part of the project may languish
unless others show interest.

# Disclaimer
This software is provided by the copyright holders and contributors "as is" and
any express or implied warranties, including, but not limited to, the implied
warranties of merchantability and fitness for a particular purpose are
disclaimed. In no event shall the copyright holder or contributors be liable for
any direct, indirect, incidental, special, exemplary, or consequential damages
(including, but not limited to, procurement of substitute goods or services;
loss of use, data, or profits; or business interruption) however caused and
on any theory of liability, whether in contract, strict liability, or tort
(including negligence or otherwise) arising in any way out of the use of this
software, even if advised of the possibility of such damage.

# Note on licensing
The source files are under the BSD license, in accordance with derivation
from EigTool. However, an executable built with the GUI
would include Qt code, which has additional restrictions.
