/* -*-javascript-*- */
import QtQuick 2.0
import QtQuick.Controls 1.0
import QtQuick.Layouts 1.3
import QtQuick.Dialogs 1.2
import org.julialang 1.0

/*
 * WARNING: some comments are obsolete!
 * This is an incomplete refactor of a 3-year old code, a Julian dinosaur.
 */

/*
 * Outline:
 *  connections
 *  dialog boxes
 *  menu bar
 *  main body (ColumnLayout)
 *    functions
 *    signals
 *    header (rectangles)
 *    graphics row layout
 *    axis, contours, etc. row layout
 *    footer (rectangles)
 *
 */

ApplicationWindow {
    id: appEps
    title: "Pseudospectrum Application"
    width: 960
    height: 640
    visible: true

    // Set up timer connection
    Connections {
	target: timer
	function onTimeout() { Julia.monitor(); }
    }

    // diaglog box for setting a parameter
    Dialog {
	id: setParamDlg
        property string curparam
	visible: false
	modality: Qt.WindowModal
	title: "Set a parameter"
	standardButtons: StandardButton.Ok | StandardButton.Cancel
	onAccepted: Julia.setoptbydlg(curparam, paramDlgParam.text)
	ColumnLayout {
	    spacing:6
	    Text {
		id: setParamDlgText
		text: "author: please rewrite setParamDlgText"
	    }
	    TextField {
		id: paramDlgParam
		Layout.alignment: Qt.AlignCenter
		Layout.minimumWidth: 150
		Layout.preferredWidth: 200
		placeholderText: qsTr("?")
		onAccepted: Julia.setoptbydlg(setParamDlg.curparam, paramDlgParam.text)
	    }
	}
    }

    // dialog box for 
    Dialog {
	id: newMtxDlg
	visible: false
	modality: Qt.WindowModal
	title: "Load a New Matrix"
	standardButtons: StandardButton.Ok | StandardButton.Cancel
	onAccepted: Julia.loadmtx(newMtxName.text,newMtxExpr.text)
	ColumnLayout {
	    spacing:6
	    Text {
		id: newDlgText
		text: "Load a new matrix into EPS-App"
	    }
	    Text {
		font.bold: true;
		text: "The expression will be evaluated,\nand the variable will reside\nin the global scope of the PSAData module."
	    }
	    Text {
		text: "Variable name for the new matrix"
	    }
            TextField {
		id: newMtxName
		Layout.alignment: Qt.AlignCenter
		Layout.minimumWidth: 150
		Layout.preferredWidth: 200
		placeholderText: qsTr("A")
            }
	    Text {
		text: "Expression to construct the matrix:"
	    }
            TextField {
		id: newMtxExpr
		Layout.alignment: Qt.AlignCenter
		Layout.minimumWidth: 300
		Layout.preferredWidth: 400
		placeholderText: qsTr("Julia expression")
            }
	}
    }

    Dialog {
	id: saveDlg
	visible: false
        property string curkey
	modality: Qt.WindowModal
	title: "Save Data to Main"
	standardButtons: StandardButton.Save | StandardButton.Cancel
	onAccepted: Julia.savedata(saveVarName.text,curkey)
	ColumnLayout {
	    spacing:6
	    Text {
		id: saveDlgText
		text: "author: please rewrite saveDlgText"
	    }
	    Text {
		text: "Specify a variable name for saved data"
	    }
	    Text {
		font.bold: true;
		text: "This variable will be in the global scope of the Main module."
	    }

            TextField {
		id: saveVarName
		Layout.alignment: Qt.AlignCenter
		Layout.minimumWidth: 150
		Layout.preferredWidth: 200
		placeholderText: qsTr("foobar")
            }
	}
    }

    Dialog {
	id: pickZDlg
	visible: false
	modality: Qt.WindowModal
	title: "Select a Point"
	standardButtons: StandardButton.Ok | StandardButton.Cancel
	onAccepted: Julia.setselectedz(zReal.text,zImag.text)
	ColumnLayout {
	    spacing:6
	    Text {
		id: pickZDlgText
		text: "author: please rewrite pickZDlgText"
	    }
	    RowLayout {
		spacing: 6
		Text {
		    text: "Re(z): "
		}
		TextField {
		id: zReal
		Layout.alignment: Qt.AlignCenter
		Layout.minimumWidth: 150
		Layout.preferredWidth: 200
		placeholderText: qsTr("0.0")
		}
		Text {
		    text: "Im(z): "
		}
		TextField {
		id: zImag
		Layout.alignment: Qt.AlignCenter
		Layout.minimumWidth: 150
		Layout.preferredWidth: 200
		placeholderText: qsTr("0.0")
		}
	    }
	}
    }

    menuBar: MenuBar {
	Menu {
            title: "\&File"
	    MenuItem {
		text: "\&New Matrix"
		onTriggered: {
			newMtxDlg.visible = true
		}
		enabled: true
	    }
/*
            // currently failing
	    MenuItem {
		text: "\&Quit"
		shortcut: StandardKey.Quit
		onTriggered: Qt.Quit()
	    }
*/
/*
            MenuItem {
            text: "Open"
            shortcut: "Ctrl+O"
            onTriggered: { some_function(); }
            enabled: false
            }
            MenuSeparator { }
*/
            /* can also have Menu for submenu */
            Menu {
		title: "E\&xport..."
		MenuItem {
		    text: "Eigen\&values"
		    /* onTriggered: */
		    enabled: false
		}
		MenuItem {
		    text: "\&Portrait data"
		    onTriggered: {
			saveDlgText.text = "Save spectral portrait data (Portrait)"
			saveDlg.curkey = "portrait"
			saveDlg.visible = true
		    }
		    enabled: true
		}
		MenuItem {
		    text: "\&Everything"
		    onTriggered: {
			saveDlgText.text = "Save full data structure (PSAStruct)"
			saveDlg.curkey = "ps_data"
			saveDlg.visible = true
		    }
		    enabled: true
		}
            }
	} // end of File
	Menu {
            title: "\&Numbers"
            MenuItem {
		text: "Pseudospectral \&abscissa"
		onTriggered: {
		    setParamDlgText.text = "Set epsilon for abscissa";
		    setParamDlg.curparam = "AbscissaEps";
		    paramDlgParam.text = "0";
		    setParamDlg.visible = true;
		}
		enabled: true /* not needed if true */
            }
            MenuItem {
		text: "Pseudospectral \&radius"
		onTriggered: {
		    setParamDlgText.text = "Set epsilon for radius";
		    setParamDlg.curparam = "RadiusEps";
		    paramDlgParam.text = "0";
		    setParamDlg.visible = true;
		}
		enabled: true /* not needed if true */
            }
	    // "Departure from normality"
	    // "Numerical abscissa"
	    // "Display points"
	    // "Eigenvector matrix condition nr."
	} // end of Numbers
	Menu {
	    title: "\&Transients"
	    MenuItem {
		id: mtxPowers
		text: "Matrix \&powers"
		onTriggered: {
		    // FIXME: get App property var
		    Julia.mtxpwrcb(50)
		}
	    }
	    MenuItem {
		id: mtxExps
		text: "Matrix \&exponentials"
		onTriggered: {
		    // FIXME: get App property var
		    Julia.mtxexpcb(0.1,50)
		}
	    }
	    MenuItem {
		text: "\&No. of steps for transient plots"
		onTriggered: {
		    setParamDlgText.text = "Set step count for transient plots";
		    setParamDlg.curparam = "TransientNmax";
		    paramDlgParam.text = "50";
		    setParamDlg.visible = true;
		}
		enabled: false /* not needed if true */
	    }
	    // ET has additional menuitems:
	    // "Compute a bound"
	    // "Best estimate lower bound"
	} // end of Transients
	Menu {
            title: "E\&xtras"
            MenuItem {
		id: showImagAxis
		text: "Display \&imaginary axis"
		checkable: true
		checked: false
		onToggled: {
		    if (showImagAxis.checked)
			Julia.setoptbykey("showimagax","true");
		    else
			Julia.setoptbykey("showimagax","false");
		}
            }
            MenuItem {
		id: showUnitCircle
		text: "Display \&unit circle"
		checkable: true
		checked: false
		onToggled: {
		    if (showUnitCircle.checked)
			Julia.setoptbykey("showunitcircle","true");
		    else
			Julia.setoptbykey("showunitcircle","false");
		}
            }
            MenuItem {
		text: "\&Projection level"
		onTriggered: {
		    setParamDlgText.text = "Set projection level";
		    setParamDlg.curparam = "ProjLev";
		    paramDlgParam.text = "0";
		    setParamDlg.visible = true;
		}
		//enabled: false /* not needed if true */
            }
	} // end of Extras
	Menu {
	    title: "\&ARPACK/eigs"
	    MenuItem {
		id: arpackMaxiter
		text: "\&Maxiter"
		onTriggered: {
		    setParamDlgText.text = "Set Maxiter for ARPACK/Eigs";
		    setParamDlg.curparam = "ArpackMaxiter";
		    paramDlgParam.text = root.arpMaxiter;
		    setParamDlg.visible = true;
		}
	    }
	    MenuItem {
		id: arpackNcv
		text: "\&NCV"
		onTriggered: {
		    setParamDlgText.text = "Set No. Krylov vecs for ARPACK/Eigs";
		    setParamDlg.curparam = "ArpackNcv";
		    paramDlgParam.text = root.arpNcv;
		    setParamDlg.visible = true;
		}
	    }
	    MenuItem {
		id: arpackTol
		text: "\&Tol"
		onTriggered: {
		    setParamDlgText.text = "Set Tolerance for ARPACK/Eigs";
		    setParamDlg.curparam = "ArpackTol";
		    paramDlgParam.text = root.arpTol;
		    setParamDlg.visible = true;
		}
	    }
	} // end of ARPACK/eigs menu
	// TODO:
	/*
	Menu {
	    title: "\&Demos"
	} // end of Demos menu
	Menu {
	    title: "\&Help"
	    MenuItem {
		id: aboutMI
		text: "About"
		onToggled: {
		}
	    }
	} // end of Help menu
	*/
    } /* end of menubar */

    ColumnLayout {
	id: root
	spacing: 6
	anchors.fill: parent
	/* Layout.fillHeight: true */
	property string arpMaxiter
	property string arpNcv
	property string arpWhichIdx
	property string arpTol

	function do_plot()
	{
            if(jdisp === null)
		return;
            Julia.mainplot(jdisp)
	}

	function init_and_plot()
	{
            if(jdisp === null)
		return;
            Julia.init_backend(jdisp.width, jdisp.height);
            do_plot();
	}

	function do_plot2()
	{
            if(jdisp2 === null)
		return;

            Julia.otherplot(jdisp2)
	}

	function init_and_plot2()
	{
            if(jdisp2 === null)
		return;

            Julia.init_backend(jdisp2.width, jdisp2.height);
            do_plot2();
	}

	JuliaSignals {
	    /*
	     * JuliaSignals allow Julia code to call into the GUI thread
	     */
	    /*
	    signal redrawAppW()
	    onRedrawAppW: {
		appEps.show()
	    }
	    */
            signal setMainMsg(var msg, var which) // string, int
            onSetMainMsg: {
		if (which == 0) {
		    juliaMainMsg.text = msg;
		    juliaMainMsg.visible = true;
		}
		if (which == 1) {
		    juliaInfoMsg.text = msg;
		    juliaInfoMsg.color = "black";
		    juliaInfoMsg.visible = true;
		}
		if (which == 2) {
		    juliaInfoMsg.text = msg;
		    juliaInfoMsg.color = "red";
		    juliaInfoMsg.visible = true;
		}
	    }
	    /*
	     * use signals to enable/disable assoc. items, set fields, etc.
	     */
            signal enableGo(var goval) // int
	    onEnableGo: {
		goBtn.enabled = (goval > 0);
	    }

            signal toggleFoV(var x) // int
	    onToggleFoV: {
		fovCkBox.checked = (x > 0);
		fovCkBox.enabled = (x >= 0);
	    }

	    signal toggleAutoAx(var x) // int
	    onToggleAutoAx: {
		autoaxisCkBox.checked = (x > 0);
		autoaxisCkBox.enabled = (x >= 0);
	    }

	    signal toggleDirect(var x) // int
	    onToggleDirect: {
		directBtn.checked = (x > 0);
		eigsBtn.checked = (x == 0);
	    }

	    signal showPlot(var flag) // int
	    onShowPlot: {
		jdispstack.currentIndex = flag;
	    }

	    signal showPlot2(var flag) // int
	    onShowPlot2: {
		jdisp2stack.currentIndex = flag;
	    }
/*
	    signal clearPlot(int idnum)
	    onClearPlot: {
		if (idnum == 1)
		    Julia.cleardisplay(jdisp);
		else
		    Julia.cleardisplay(jdisp2);
	    }
*/

	    /* toString() is inadequate for these, so preformat them */
	    // args are strings
	    signal setAxisLims(var x0, var x1, var y0, var y1)
	    onSetAxisLims: {
		xMin.text = x0;
		xMax.text = x1;
		yMin.text = y0;
		yMax.text = y1;
	    }
	    // args are double
	    signal setLevels(var x0, var x1, var x2)
	    onSetLevels: {
		ctrMin.text = x0.toString();
		ctrMax.text = x1.toString();
		ctrStep.text = x2.toString();
	    }
	    signal setGrid(var x0) // int
	    onSetGrid: {
		gridSize.text = x0.toString();
	    }
	    // int, int, int, double, int
	    signal setArpackOpts(var nev, var ncv, var maxiter, var tol,
				 var whichidx)
	    onSetArpackOpts: {
		root.arpMaxiter = maxiter.toString();
		kArpack.text = nev.toString();
		root.arpNcv = ncv.toString();
		root.arpTol = tol.toString();
		whichArpackBox.currentIndex = whichidx;
	    }
	    signal saveVariable(var dlgtext, var keytext) // strings
	    onSaveVariable: {
		saveDlgText.text = dlgtext;
		saveDlg.curkey = keytext;
		saveDlg.visible = true;
	    }

/*
// deadlocks if invoked in a task
	    signal getZ(string str)
	    onGetZ: {
		pickZDlgText.text = str;
		pickZDlg.visible = true;
	    }
*/
	    signal runTimer(var x) // int
	    onRunTimer: {
		if (x > 0) {
		    timer.interval = 10;
		    timer.start();
		}
		else
		    timer.stop();
	    }
	}

/*
 * Here begins the layout of the main GUI.
 * The first few rectangles are the header.
 */
	Rectangle {
            color: "black"
            Layout.preferredHeight: 1
            Layout.fillWidth: true
	}
	Rectangle {
            id: mainmsgRect
            color: "lightsteelblue"
            Layout.fillWidth: true
            Layout.preferredHeight: 20
            /* Layout.margins: 5 */
            Layout.alignment: Qt.AlignTop

            Text {
		id: juliaMainMsg
		Layout.minimumHeight: 30
		anchors.left: parent.left
		anchors.leftMargin: 6
		anchors.verticalCenter: parent.verticalCenter
		font.pointSize: 20; font.bold: true;
		text: Julia.initmsg() /* only used for initialization */
            }
	}

	Rectangle {
            color: "black"
            Layout.preferredHeight: 1
            Layout.fillWidth: true
	}

	/*
	 * The main body is a row: stack-column-stack.
	 */
	RowLayout {
            Layout.fillWidth: true
            spacing: 6

	    /* This is a workaround allowing us to hide the plots when
	     * they are invalid (e.g. from an old matrix).
	     * QML.jl doesn't wrap a clear() method for JuliaDisplay
	     * so we swap with a blank rectangle via stack indexing.
	     * This has the nice side-effect of handling app-resizes.
	     */
	    StackLayout {
		id:jdispstack
		Layout.alignment: Qt.AlignLeft
		Layout.preferredWidth: 400
		Layout.preferredHeight: 380

		Rectangle {
		    id: jdrect
		    implicitWidth: 400
		    implicitHeight: 380
		}

		JuliaDisplay {
		    id: jdisp
		    implicitWidth: 400
		    implicitHeight: 380
		    onHeightChanged: root.init_and_plot()
		    onWidthChanged: root.init_and_plot()
		    // apparently wanted inside the display
		    MouseArea {
			id: jdispev
			anchors.fill:parent
			acceptedButtons: Qt.LeftButton | Qt.RightButton
			onClicked: {
			    if (mouse.button == Qt.LeftButton)
				Julia.mainzoom(0, mouse.x, mouse.y, jdisp.height);
			    else
				Julia.mainzoom(1, mouse.x, mouse.y, jdisp.height);
			}
		    }
		}
	    }

            ColumnLayout {
		/*          Layout.fillHeight: true */
		Layout.alignment: Qt.AlignTop
		spacing: 6

		GroupBox { /* Rectangle failed */
		    title: "Commands:"
		    id: goRect
		    Layout.alignment: Qt.AlignTop
		    ColumnLayout {
			Button {
			    id: goBtn
			    Layout.alignment: Qt.AlignCenter
			    // apparently colors need button styles
			    // color: "green"
			    text: "Update"
			    onClicked: Julia.go()
			    enabled: false
			}
			/*
			// placeholder pending a really nasty coding session
			Button {
			id: pauseBtn
			Layout.alignment: Qt.AlignCenter
			// color: "yellow"
			text: "Pause"
			onClicked: Julia.psapause()
			}
			*/
			CheckBox {
			    id: fovCkBox
			    text: "Show numerical range"
			    checked: false
			    onClicked: {
				if (fovCkBox.checked)
				    Julia.setoptbykey("showFov","true");
				else
				    Julia.setoptbykey("showFov","false");
			    }
			}
			Button {
			    Layout.alignment: Qt.AlignCenter
			    text: "Select z"
			    onClicked: {
				pickZDlg.visible = true;
			    }
			}
			Button {
			    Layout.alignment: Qt.AlignCenter
			    text: "Eigenmode"
			    onClicked: {
				Julia.ploteigv()
			    }
			}
			Button {
			    Layout.alignment: Qt.AlignCenter
			    text: "Pseudomode"
			    onClicked: Julia.plotpmode()
			}
			Button {
			    Layout.alignment: Qt.AlignCenter
			    text: "3D Plot"
			    onClicked: Julia.plotps3d()
			}
			/*
			  Button {
			  Layout.alignment: Qt.AlignCenter
			  text: "Quit"
			  onClicked: Julia.quitgui()
			  }
			*/
		    }
		}

		GroupBox {
		    title: "Mesh:"
		    Layout.alignment: Qt.AlignBottom
		    RowLayout {
			Text {
			    text: "Grid Size: "
			}
			TextField {
			    id: gridSize
			    Layout.alignment: Qt.AlignCenter
			    Layout.minimumWidth: 30
			    Layout.preferredWidth: 40
			    placeholderText: qsTr("0")
			    validator: IntValidator {}
			    onTextChanged: Julia.setoptbykey("npts",gridSize.text)
			}
		    }
		}
            }

	    StackLayout {
		id:jdisp2stack
		Layout.alignment: Qt.AlignRight
		Layout.preferredWidth: 400
		Layout.preferredHeight: 380

		Rectangle {
		    implicitWidth: 400
		    implicitHeight: 380
		}

		JuliaDisplay {
		    id: jdisp2
		    implicitWidth: 400
		    implicitHeight: 380
		    onHeightChanged: root.init_and_plot2()
		    onWidthChanged: root.init_and_plot2()
		}
	    }

	} /* main middle row */

	RowLayout {

            Layout.alignment: Qt.AlignBottom
            Layout.margins: 5
            spacing: 6

            GroupBox {
		title: "Figure Axes: "
		Layout.alignment: Qt.AlignTop
		ColumnLayout{
		    Layout.alignment: Qt.AlignCenter
		    /*
		      ExclusiveGroup { id:axlims }
		      RadioButton {
		      text: "Auto"
		      checked: true
		      exclusiveGroup: axlims
		      onClicked: Julia.setoptbykey("autolims","true")
		      }
		      RadioButton {
		      text: "Take from plot"
		      exclusiveGroup: axlims
		      onClicked: Julia.set_axes_fig()
		      }
		      RadioButton {
		      text: "Set from here"
		      exclusiveGroup: axlims
		      onClicked: Julia.set_axes(xMin,yMin,xMax,yMax)
		      }
		    */
		    CheckBox {
			id: autoaxisCkBox
			text: "Auto"
			checked: true
		    }
		    RowLayout {
			Layout.alignment: Qt.AlignCenter
			Text {
			    text: "Y max: "
			}
			TextField {
			    id: yMax
			    enabled: { !autoaxisCkBox.checked; }
			    Layout.alignment: Qt.AlignCenter
			    Layout.minimumWidth: 20
			    Layout.preferredWidth: 70
			    placeholderText: qsTr("0")
			    validator: DoubleValidator {}
			    onEditingFinished: Julia.setoptbykey("yMax",yMax.text)
			}
		    }

		    RowLayout {
			Text {
			    text: "X min: "
			}
			TextField {
			    id: xMin
			    enabled: { !autoaxisCkBox.checked; }
			    Layout.alignment: Qt.AlignCenter
			    Layout.minimumWidth: 20
			    Layout.preferredWidth: 70
			    placeholderText: qsTr("0")
			    validator: DoubleValidator {}
			    onEditingFinished: Julia.setoptbykey("xMin",xMin.text)
			}
			Text {
			    text: " X max: " /* note cheap kerning */
			}
			TextField {
			    id: xMax
			    enabled: { !autoaxisCkBox.checked; }
			    Layout.alignment: Qt.AlignCenter
			    Layout.minimumWidth: 20
			    Layout.preferredWidth: 70
			    placeholderText: qsTr("0")
			    validator: DoubleValidator {}
			    onEditingFinished: Julia.setoptbykey("xMax",xMax.text)
			}
		    }

		    RowLayout {
			Layout.alignment: Qt.AlignCenter
			Text {
			    text: "Y min: "
			}
			TextField {
			    id: yMin
			    enabled: { !autoaxisCkBox.checked; }
			    Layout.alignment: Qt.AlignCenter
			    Layout.minimumWidth: 20
			    Layout.preferredWidth: 70
			    placeholderText: qsTr("0")
			    validator: DoubleValidator {}
			    onEditingFinished: Julia.setoptbykey("yMin",yMin.text)
			}
		    }

		}
            }

            GroupBox {
		title: "Direct/Iterative: "
		Layout.alignment: Qt.AlignTop
		ColumnLayout{
		    ExclusiveGroup { id:dirvsiter }
		    RadioButton {
			id: directBtn
			text: "Direct"
			checked: true
			exclusiveGroup: dirvsiter
			onClicked: Julia.use_eigs(0)
		    }
		    RadioButton {
			id: eigsBtn
			text: "ARPACK/eigs"
			exclusiveGroup: dirvsiter
			onClicked: Julia.use_eigs(1)
		    }
		    RowLayout {
			Text {
			    text: "No. eigs: "
			}
			TextField {
			    id: kArpack
			    enabled: eigsBtn.checked
			    Layout.alignment: Qt.AlignCenter
			    Layout.preferredWidth: 40
			    placeholderText: qsTr("6")
			    inputMethodHints: Qt.ImhDigitsOnly
			    onEditingFinished: Julia.setoptbykey("kArpack",kArpack.text)
			}
		    }
		    ComboBox {
			id: whichArpackBox
			currentIndex: 0
			property string curkey
			enabled: eigsBtn.checked

			model: ListModel {
			    id: arpackcb

			    // Warning: order must agree w/ which2idx
			    // in Julia code

			    ListElement {
				text: "Large Mod."
				key: "LM"
			    }
			    /*
			    // not implemented upstream either
			    ListElement {
			    text: "Small Mod."
			    key: "SM"
			    }
			    */
			    ListElement {
				text: "Large real"
				key: "LR"
			    }
			    ListElement {
				text: "Small real"
				key: "SR"
			    }
			    ListElement {
				text: "Large imag."
				key: "LI"
			    }
			    ListElement {
				text: "Small imag."
				key: "SI"
			    }
			}
			onCurrentIndexChanged: {
			    whichArpackBox.curkey = arpackcb.get(currentIndex).key
			    Julia.setoptbykey("which",whichArpackBox.curkey);
			}
		    }
		}
            }

            GroupBox {
		Layout.alignment: Qt.AlignTop
		title: "Contour Levels:"
		ColumnLayout{
		    Button {
			Layout.alignment: Qt.AlignCenter
			text: "Smart Levels"
			onClicked: Julia.smartlevels()
		    }
		    RowLayout {
			Layout.alignment: Qt.AlignCenter
			Text {
			    text: "log10(max): "
			}
			TextField {
			    id: ctrMax
			    Layout.alignment: Qt.AlignCenter
			    Layout.minimumWidth: 20
			    Layout.preferredWidth: 60
			    placeholderText: qsTr("0")
			    inputMethodHints: Qt.ImhFormattedNumbersOnly
			    onEditingFinished: Julia.setoptbykey("last",ctrMax.text)
			}
		    }
		    RowLayout {
			Layout.alignment: Qt.AlignCenter
			Text {
			    text: "log10(min): "
			}
			TextField {
			    id: ctrMin
			    Layout.alignment: Qt.AlignCenter
			    Layout.minimumWidth: 20
			    Layout.preferredWidth: 60
			    placeholderText: qsTr("0")
			    inputMethodHints: Qt.ImhFormattedNumbersOnly
			    onEditingFinished: Julia.setoptbykey("first",ctrMin.text)
			}
		    }
		    RowLayout {
			Layout.alignment: Qt.AlignCenter
			Text {
			    text: "Step size: "
			}
			TextField {
			    id: ctrStep
			    Layout.alignment: Qt.AlignCenter
			    Layout.minimumWidth: 20
			    Layout.preferredWidth: 60
			    placeholderText: qsTr("0")
			    inputMethodHints: Qt.ImhFormattedNumbersOnly
			    onEditingFinished: Julia.setoptbykey("step",ctrStep.text)
			}
		    }
		}
            } /* contour level group box */

	} /* main bottom row */
	Rectangle {
            color: "black"
            Layout.preferredHeight: 1
            Layout.fillWidth: true
	}
	Rectangle {
            id: infomsgRect
            color: "lightsteelblue"
            Layout.fillWidth: true
            Layout.preferredHeight: 16
            /* Layout.margins: 5 */
            Layout.alignment: Qt.AlignTop

            Text {
		id: juliaInfoMsg
		Layout.minimumHeight: 16
		anchors.left: parent.left
		anchors.leftMargin: 6
		anchors.verticalCenter: parent.verticalCenter
		font.pointSize: 12;
		// we tried a context property here, but updates are erratic
		text: "watch for news here"
            }
	}

    } /* main column */
} /* app */
