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
It is largely based on the famous
[EigTool](http://www.cs.ox.ac.uk/pseudospectra/eigtool/) from Oxford.

It uses core computational routines from a separate package,
[Pseudospectra.jl](https://github.com/RalphAS/Pseudospectra.jl).


# Warning

This package is experimental and may misbehave (i.e., some remaining
bugs may deadlock or crash a Julia session), so save any valuable
results before running the app.  It is quite complicated and has some
tricky dependencies, so it takes a while to load. (Also, see "Known issues"
below.) That said, most of the basic functionality works.
Interested users are encouraged to try it out and report problems:
please provide details of how to reproduce the latter.

# Installation and direct dependencies

The graphical user interface is based on [QML.jl](https://github.com/barche/QML.jl).
Install that and make sure it works first.

The current version of this package requires Plots.jl. The PyPlot
backend for Plots is supported; support for the GR backend is work
in progress.

Install [Pseudospectra.jl](https://github.com/RalphAS/Pseudospectra.jl); this
must be done manually, pending registration.

Finally, add this repo ("https://github.com/RalphAS/PseudospectraQML.jl")
with the package manager.

# Basics
Invoke the App with

```julia
using QML
using PseudospectraQML
A = your_matrix_generating_function()
saved_vars = psagui(A)
```

(Alternatively one can invoke `psagui()` with no arguments, and
run the `File/New Matrix` menu command, but that operates in a
context which may not have what you need for your matrix.)

Please be patient; the infrastructure takes a while to initialize.
An app window should appear. If the matrix is not large, a spectral portrait
(i.e. a contour plot of base-10 logarithms of inverse resolvent norms) should
be shown in the main display area. If the matrix is large, one can
select `Direct` if one is patient, or select parameters for the
iterative scheme in ARPACK to compute pseudospectra for a reduced
form, then press `Update` for computation and display.

The default mesh is very coarse; adjust the grid size and axis limits to
taste and press the `Update` button to recompute.

One can zoom in by clicking the left mouse button with the pointer in the
main plot area. (Zoom out with the right button.) Note that this usually
involves recomputation. ("Left" and "right" here refer to default button
assigments.)

To display an eigenmode or pseudo-eigenmode, click the `Select z` button
which brings up a dialog to specify the value, then click
the appropriate mode button.

## TODO: walk-throughs with screenshots needed here

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
saved_vars = psagui(A)
# run the app, then close the window

using Plots
gr() # or some other backend - see below
baz = saved_vars[:baz]
contour(baz.x,baz.y,baz.Z)
png("baz")
```

# Striking differences from EigTool
* One must use the "Select z" command button **before** plotting an
  eigenmode/pseudomode or transients.
* The "pause/resume/stop" logic is not implemented.

# Other notable peculiarities
* At present, using Pseudospectra.jl with command-line graphics calls
  and with the QML App in the same session requires the
  use of a `Plots` backend **other than** PyPlot for the command-line work.
  (This is because the plotter must use a
  non-interactive interface to prepare images for QML)
* The App is in its own module, which is reloaded for each instance.
* PseudospectraQML uses a separate module (PSAData) to store
  some variables.

# Known issues
* The GUI is finicky.
  * Text fields may not be properly converted to associated variables unless
    one uses the Tab key to check in the value.
  * Coordinate conversion for zooming is flaky. Sometimes resizing the
    whole app window, then forcing a redraw, makes it more accurate.
  * Obscure callback problems arise occasionally, some harmless,
	others not.
* Some informational messages are not correctly synchronized.
* Some warnings and errors provoked from the GUI are reported in the console.
* Some buttons & menu items remain enabled when incompatible with state.
* The control flow is very convoluted.
  * Zooming and the axis-setting widgets are not always coordinated.
  * Deadlocks or crashes are possible. File an issue if you encounter them.

# Caveat
The current translator/author is a novice at GUI construction, and generally
prefers command-line interfaces. Hence this project may languish
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
