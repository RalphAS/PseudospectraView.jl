# PseudospectraQML.jl

![Lifecycle](https://img.shields.io/badge/lifecycle-experimental-orange.svg)<!--
![Lifecycle](https://img.shields.io/badge/lifecycle-maturing-blue.svg)
![Lifecycle](https://img.shields.io/badge/lifecycle-stable-green.svg)
![Lifecycle](https://img.shields.io/badge/lifecycle-retired-orange.svg)
![Lifecycle](https://img.shields.io/badge/lifecycle-archived-red.svg)
![Lifecycle](https://img.shields.io/badge/lifecycle-dormant-blue.svg) -->
<!--
[![GitHub CI Build Status](https://github.com/RalphAS/PseudospectraQML.jl/workflows/CI/badge.svg)](https://github.com/RalphAS/PseudospectraQML.jl/actions)
[![codecov.io](http://codecov.io/github/RalphAS/PseudospectraQML.jl/coverage.svg?branch=master)](http://codecov.io/github/RalphAS/PseudospectraQML.jl?branch=master)
-->
<!--
[![Documentation](https://img.shields.io/badge/docs-stable-blue.svg)](https://RalphAS.github.io/PseudospectraQML.jl/stable)
[![Documentation](https://img.shields.io/badge/docs-master-blue.svg)](https://RalphAS.github.io/PseudospectraQML.jl/dev)
-->

# Introduction

PseudospectraQML implements a Julia/Qt GUI for analyzing nonsymmetric matrices by
displaying eigenvalues and pseudospectra, along with related properties.
It is largely based on the famous
[EigTool](http://www.cs.ox.ac.uk/pseudospectra/eigtool/) from Oxford.

It uses core computational routines from a separate package,
[Pseudospectra.jl](https://github.com/RalphAS/Pseudospectra.jl).

**Please check the [Warnings](README.md#warnings) section below to avoid nasty surprises.**

# Installation and direct dependencies

The graphical user interface is based on [QML.jl](https://github.com/barche/QML.jl).
Install that and make sure it works first.

The current version of this package requires Plots.jl. The GR and PyPlot
backends for Plots are supported.

Install [Pseudospectra.jl](https://github.com/RalphAS/Pseudospectra.jl); this
package is in the General registry.

Finally, add this repo ("https://github.com/RalphAS/PseudospectraQML.jl")
with the package manager.

# Basics
Invoke the App with

```julia
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
If you use the PyPlot backend for Plots, saving from the "File" menu is possible.

Otherwise, to save a plot for external use, the recommended procedure is as follows:
1. save the current data as a `Portrait` type, using the "File" menu
2. exit the App (Alt-F4 works for me)
3. use fields of the saved variable to make a plot
E.g., if you saved the `Portrait` as variable `baz`,
```julia
saved_vars = psagui(A; backend=:gr)
# run the app, then close the window

using Plots
pyplot() # you must change the backend - see below
baz = saved_vars[:baz]
contour(baz.x,baz.y,baz.Z)
png("baz")
```

# Striking differences from EigTool
* One must use the "Select z" command button **before** plotting an
  eigenmode/pseudomode or transients.
* The "pause/resume/stop" logic is not implemented.

# Other notable peculiarities
* The App is in its own module, which is reloaded for each instance.
* PseudospectraQML uses yet another module (PSAData) to store
  some variables.

# Warnings

The graphics backends used by the GUI do not reliably allow for hot-swapping of
online destinations; if you want to generate plots with other interfaces in the same
Julia session, you should use a different backend to do so.

This package is experimental and may misbehave (i.e., some remaining
bugs may deadlock or crash a Julia session), so save any valuable
results before running the app.  It is quite complicated and has some
tricky dependencies, so it takes a while to load.

# Other known issues

* The GUI is finicky.
  * Text fields may not be properly converted to associated variables unless
    one uses the Tab key to check in the value.
  * Coordinate conversion for zooming is flaky. Sometimes resizing the
    whole app window, then forcing a redraw, makes it more accurate.
  * Obscure callback (ok, "signalling") problems arise occasionally; some are harmless,
	others not.
* When secondary plots are generated under the GR backend, the primary (portrait)
  plot may be disrupted. Tinkering with mesh, contour, or axis limits settings
  and clicking "Update" usually fixes this.
* Sometimes the PyPlot graphics don't fit in their assigned canvas.
* GR often fails to display legends, titles, etc.
* PyPlot doesn't produce animated plots that I could show in the App.
* GR doesn't produce surface plots that I could show in the App.

## Just to be fair, a few for which I am responsible:
* Informational messages are not correctly synchronized.
* Some warnings and errors provoked from the GUI are only reported in the console.
* Some buttons & menu items may remain enabled when incompatible with the state.
* Some exceptions are not caught, and will return to the REPL, leaving the GUI in an
  unstable state.


# Disclaimers
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
