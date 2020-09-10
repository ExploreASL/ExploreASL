from PySide2.QtWidgets import *
from PySide2.QtGui import *
from PySide2.QtCore import *
from xASL_GUI_HelperClasses import DandD_Graphing_ListWidget2LineEdit, DandD_FileExplorer2LineEdit, \
    DandD_FileExplorerFile2LineEdit
from xASL_GUI_FacetPlot import xASL_GUI_FacetGridOrganizer
from xASL_GUI_PlotLabels import xASL_GUI_PlotLabels
from xASL_GUI_Graph_Subsetter import xASL_GUI_Subsetter
from xASL_GUI_Graph_Loader import xASL_GUI_Data_Loader
import os
import sys
import json
from glob import glob
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar
from matplotlib.figure import Figure
import seaborn as sns
from pprint import pprint

pd.set_option("display.max_columns", 10)
pd.set_option("display.width", 400)


# noinspection PyCallingNonCallable
class xASL_PostProc(QMainWindow):
    """
    Central Graphing Widget which will act as a "central hub" for other widgets to be placed in
    """
    def __init__(self, parent_win=None):
        # Parent window is fed into the constructor to allow for communication with parent window devices
        super().__init__(parent=parent_win)

        # Window Size and initial visual setup
        self.setMinimumSize(1920, 1000)
        self.cw = QWidget(self)
        self.setCentralWidget(self.cw)
        self.mainlay = QVBoxLayout()
        self.cw.setLayout(self.mainlay)
        self.setWindowTitle("Explore ASL - Post Processing Visualization")
        self.setWindowIcon(QIcon(os.path.join(os.getcwd(), "media", "ExploreASL_logo.png")))

        self.canvas_generate(None)

        # Central Classes
        self.subsetter = xASL_GUI_Subsetter(self)
        self.loader = xASL_GUI_Data_Loader(self)
        self.plotlabels = xASL_GUI_PlotLabels(self)

        # Initialize blank givens
        self.fig_manager = None
        self.fig_wid = None
        self.axes_wid = None
        self.loaded_wide_data = pd.DataFrame()
        self.loaded_long_data = pd.DataFrame()
        self.plotstylenames = ['seaborn-ticks', 'seaborn-white', 'seaborn-whitegrid', 'seaborn-dark',
                               'seaborn-darkgrid', 'seaborn-paper', 'seaborn-poster', 'seaborn-talk',
                               'bmh', 'classic', 'dark_background', 'fivethirtyeight', 'ggplot']
        # Main Widgets setup
        self.UI_Setup_Docker()

    def UI_Setup_Docker(self):
        self.dock = QDockWidget("Data Visualization Settings", self.cw)
        self.dock.setMinimumWidth(480)
        self.dock.setFeatures(QDockWidget.AllDockWidgetFeatures)
        self.addDockWidget(Qt.LeftDockWidgetArea, self.dock)
        self.cont_maindock = QWidget(self.dock)
        self.vlay_maindock = QVBoxLayout(self.cont_maindock)
        self.dock.setWidget(self.cont_maindock)

        # Set up the directory settings - what analysis folder
        self.grp_directories = QGroupBox("Directory settings", self.cont_maindock)
        self.formlay_directories = QFormLayout(self.grp_directories)
        self.le_analysis_dir = DandD_FileExplorer2LineEdit(self.grp_directories)
        self.le_analysis_dir.setText(self.parent().config["DefaultRootDir"])
        self.le_demographics_file = DandD_FileExplorerFile2LineEdit([".tsv", ".csv", ".xlsx"], self.grp_directories)
        self.le_demographics_file.setPlaceholderText("Drag & Drap a supporting .tsv/.csv/.xlsx file")
        self.cmb_atlas_selection = QComboBox(self.grp_directories)
        self.cmb_atlas_selection.addItems(["MNI", "Hammers"])
        self.cmb_pvc_selection = QComboBox(self.grp_directories)
        self.cmb_pvc_selection.addItems(["Without Partial Volume Correction", "With Partial Volume Correction"])
        self.cmb_stats_selection = QComboBox(self.grp_directories)
        self.cmb_stats_selection.addItems(["Mean", "Median", "Coefficient of Variation"])
        self.btn_subset_data = QPushButton("Subset Data", self.grp_directories, clicked=self.subsetter.show)
        self.btn_subset_data.setEnabled(False)
        self.btn_load_in_data = QPushButton("Load Data", self.grp_directories, clicked=self.loader.load_exploreasl_data)
        self.formlay_directories.addRow("Analysis Directory", self.le_analysis_dir)
        self.formlay_directories.addRow("Ancillary Study Dataframe", self.le_demographics_file)
        self.formlay_directories.addRow("Which Atlas to Utilize", self.cmb_atlas_selection)
        self.formlay_directories.addRow("Which Partial-Volume Stats to View", self.cmb_pvc_selection)
        self.formlay_directories.addRow("Which Statistic to View", self.cmb_stats_selection)
        self.formlay_directories.addRow(self.btn_subset_data)
        self.formlay_directories.addRow(self.btn_load_in_data)

        # Connect the appropriate lineedits to the subsetter class
        self.le_analysis_dir.textChanged.connect(self.subsetter.clear_contents)
        self.le_demographics_file.textChanged.connect(self.subsetter.clear_contents)

        # Setup the main Variable Viewer
        self.grp_varview = QGroupBox("Variables", self.cont_maindock)
        self.vlay_varview = QVBoxLayout(self.grp_varview)
        self.lst_varview = QListWidget(self.grp_varview)
        self.lst_varview.setFixedHeight(250)
        self.lst_varview.setDragEnabled(True)
        self.vlay_varview.addWidget(self.lst_varview)

        # Setup the start of Plotting Settings
        self.grp_pltsettings = QGroupBox("Plotting Settings", self.cont_maindock)
        self.vlay_pltsettings = QVBoxLayout(self.grp_pltsettings)
        self.cmb_figuretypeselection = QComboBox(self.grp_pltsettings)
        self.cmb_figuretypeselection.addItems(["Select an option", "Facet Grid"])
        self.cmb_figuretypeselection.setEnabled(False)
        self.tab_pltsettings = QTabWidget(self.grp_pltsettings)

        self.cont_pltsettings_univlvlparms = QWidget(self.grp_pltsettings)
        self.cont_pltsettings_figlvlparms = QWidget(self.grp_pltsettings)
        self.cont_pltsettings_axeslvlparms = QWidget(self.grp_pltsettings)

        self.vlay_pltsetting_figlvlparms = QVBoxLayout(self.cont_pltsettings_figlvlparms)
        self.vlay_pltsettings_axeslvlparms = QVBoxLayout(self.cont_pltsettings_axeslvlparms)

        self.tab_pltsettings.addTab(self.cont_pltsettings_univlvlparms, "Common Parameters")
        self.tab_pltsettings.addTab(self.cont_pltsettings_figlvlparms, "Figure Parameters")
        self.tab_pltsettings.addTab(self.cont_pltsettings_axeslvlparms, "Axes Parameters")
        self.tab_pltsettings.setMaximumHeight(400)

        self.vlay_pltsettings.addWidget(self.cmb_figuretypeselection)
        self.vlay_pltsettings.addWidget(self.tab_pltsettings)

        # Add the groups to the main vertical layout
        self.vlay_maindock.addWidget(self.grp_directories)
        self.vlay_maindock.addWidget(self.grp_varview)
        self.vlay_maindock.addWidget(self.grp_pltsettings)
        self.vlay_maindock.addStretch()

        # This has to go at the end due to order of initialization
        self.UI_Setup_CommonParameters()
        self.cmb_figuretypeselection.currentTextChanged.connect(self.set_current_figure_manager)

    # These parameters are universal to all plots and can therefore be set up here
    def UI_Setup_CommonParameters(self):
        # Format layout must be the first item declared
        self.formlay_commonparms = QFormLayout(self.cont_pltsettings_univlvlparms)

        # Set up the overall plot style
        self.cmb_plotstyle = QComboBox(self.cont_pltsettings_univlvlparms)
        self.cmb_plotstyle.addItems(self.plotstylenames)
        self.cmb_plotstyle.currentTextChanged.connect(self.plotupdate_plotstyle)
        self.formlay_commonparms.addRow("Overall Plot Style", self.cmb_plotstyle)

        # Set up the padding
        self.cont_pltsettings_univlvlparms.setLayout(self.formlay_commonparms)
        self.spinbox_toppad = QDoubleSpinBox(self.cont_pltsettings_univlvlparms)
        self.spinbox_bottompad = QDoubleSpinBox(self.cont_pltsettings_univlvlparms)
        self.spinbox_leftpad = QDoubleSpinBox(self.cont_pltsettings_univlvlparms)
        self.spinbox_rightpad = QDoubleSpinBox(self.cont_pltsettings_univlvlparms)
        self.spinbox_wspace = QDoubleSpinBox(self.cont_pltsettings_univlvlparms)
        self.spinbox_hspace = QDoubleSpinBox(self.cont_pltsettings_univlvlparms)
        for description, widget, default in zip(["Top Padding", "Bottom Padding", "Left Padding", "Right Padding",
                                                 "Horizontal Inter-Facet Spacing", "Vertical Inter-Facet Spacing"],
                                                [self.spinbox_toppad, self.spinbox_bottompad, self.spinbox_leftpad,
                                                 self.spinbox_rightpad, self.spinbox_wspace, self.spinbox_hspace],
                                                [0.91, 0.057, 0.055, 0.95, 0.2, 0.2]):
            widget.setValue(default)
            widget.setRange(0, 1)
            widget.setSingleStep(0.01)
            widget.valueChanged.connect(self.plotupdate_padding)
            self.formlay_commonparms.addRow(description, widget)

        # Set up the text options and connect the plotlabels signals to the appropriate function
        self.btn_showplotlabels = QPushButton("Show Plot Label Settings", self.cont_pltsettings_univlvlparms,
                                              clicked=self.plotlabels.show)
        self.formlay_commonparms.addRow(self.btn_showplotlabels)
        self.plotlabels.signal_xaxislabel_changed.connect(self.plotupdate_xlabel)
        self.plotlabels.signal_yaxislabel_changed.connect(self.plotupdate_ylabel)
        self.plotlabels.signal_title_changed.connect(self.plotupdate_title)

    # When called, this updates the general plot style occuring for figures
    def plotupdate_plotstyle(self):
        plt.style.use(self.cmb_plotstyle.currentText())
        print(f"Using style: {self.cmb_plotstyle.currentText()}")
        # First clear the canvas
        self.canvas_full_clear()
        # Then generate the new canvas with the updated facetgrid and figure
        self.canvas_generate("Facet Grid Figure Call")
        # Then see if the axes can be updated, and if so, update them
        self.plotupdate_facetgrid_axes()
        self.plotupdate_padding()


    @Slot(dict)
    def plotupdate_xlabel(self, xlabel_kwargs):
        for ax in self.mainfig.axes:
            plt.sca(ax)
            plt.xlabel(**xlabel_kwargs)
        self.canvas.draw()

    @Slot(dict)
    def plotupdate_ylabel(self, ylabel_kwargs):
        for ax in self.mainfig.axes:
            plt.sca(ax)
            plt.ylabel(**ylabel_kwargs)
        self.canvas.draw()

    @Slot(dict)
    def plotupdate_title(self, title_kwargs):
        plt.title(**title_kwargs)
        self.canvas.draw()

    # When called, this updates the padding to whatever is the current value in the padding doublespinboxes
    def plotupdate_padding(self):
        print("Adjusting padding")
        plt.subplots_adjust(left=self.spinbox_leftpad.value(),
                            bottom=self.spinbox_bottompad.value(),
                            right=self.spinbox_rightpad.value(),
                            top=self.spinbox_toppad.value(),
                            wspace=self.spinbox_wspace.value(),
                            hspace=self.spinbox_hspace.value())
        self.canvas.draw()

    # Controls the current type of figure and its arguments in the Figure Parameters Tab
    # by necessity, it must also clear the contents of whatever is in the axes table
    def set_current_figure_manager(self, fig_type):
        if fig_type == "Facet Grid":
            # Must first clear the current canvas and set one in that is prepared with a facet grid
            self.canvas_full_clear()
            plt.style.use(self.cmb_plotstyle.currentText())
            print(f"Using style: {self.cmb_plotstyle.currentText()}")
            self.canvas_generate("Facet Grid Default")
            self.plotupdate_padding()
            self.fig_manager = xASL_GUI_FacetGridOrganizer(self)
            self.fig_manager.change_figparms_updateplot.connect(self.plotupdate_facetgrid_figurecall_plot)
            self.fig_manager.change_axesparms_updateplot.connect(self.plotupdate_facetgrid_axescall_plot)
            self.fig_manager.change_legendparms_updateplot.connect(self.plotupdate_facetgrid_legend)
            self.fig_manager.change_axesparms_widget.connect(self.update_axesparms)

            # Clear and replace the figure parameters widget for the figure parameters tab. The new figure parameters
            # widget is specific to the currently-assigned manager
            self.clear_figparms()
            self.fig_wid = self.fig_manager.cont_figparms
            self.vlay_pltsetting_figlvlparms.addWidget(self.fig_wid)

            # Must also clear the axes widget for the axes parameters tab
            self.clear_axesparms()

        else:
            self.canvas_full_clear()
            self.canvas_generate(None)
            self.plotupdate_padding()

            # Delete the figure manager just in case
            if self.fig_manager:
                del self.fig_manager
                self.fig_manager = None

            # Clear the figure widget tab
            self.clear_figparms()
            self.clear_axesparms()

    # Convenience Function - clears the Figure Parameters Tab
    def clear_figparms(self):
        if self.fig_wid is not None:
            self.vlay_pltsetting_figlvlparms.removeWidget(self.fig_wid)
            self.fig_wid.setParent(None)
            del self.fig_wid
            self.fig_wid = None

    # Convenience Function - clears the Axes Parameters Tab
    def clear_axesparms(self):
        print("Clearing Axes Parms")
        if self.axes_wid is not None:
            self.vlay_pltsettings_axeslvlparms.removeWidget(self.axes_wid)
            self.axes_wid.setParent(None)
            del self.axes_wid
            self.axes_wid = None

    # Convenience Function - updates the "Axes Parameters" tab with the new tab that must have been prepared
    @Slot()
    def update_axesparms(self):
        print("update_axesparms received a signal")
        self.axes_wid = self.fig_manager.cont_axesparms
        self.vlay_pltsettings_axeslvlparms.addWidget(self.axes_wid)

    # Convenience Function - clears the canvas
    def canvas_full_clear(self):
        plt.clf()
        plt.close(self.mainfig)
        self.mainlay.removeWidget(self.canvas)
        self.mainlay.removeWidget(self.nav)
        self.canvas.setParent(None)
        self.nav.setParent(None)
        del self.nav
        del self.mainfig
        del self.canvas

    # Called after a clear to restore a blank/default canvas with the appropriate figure manager
    def canvas_generate(self, action):

        if action is None:
            self.mainfig = Figure()
        # In "Facet Grid Default" a simple placeholder canvas is set
        elif action == 'Facet Grid Default':
            self.grid = sns.FacetGrid(data=self.loader.long_data)
            self.mainfig = self.grid.fig
        # In "Facet Grid Figure Call", a change in the figure parms is forcing an update
        elif action == "Facet Grid Figure Call":
            # For false call with style change
            if self.cmb_figuretypeselection == 'Select an option':
                return

            # First create the constructor, then feed it into the Facetgrid class
            constructor = {}
            for kwarg, call in self.fig_manager.fig_kwargs.items():
                if call() == "":
                    constructor[kwarg] = None
                else:
                    constructor[kwarg] = call()

            # Create the FacetGrid; take note that if the grid type is a regression plot, then the palette argument
            # is fed into the grid call instead of the axes contructor
            if not self.cmb_figuretypeselection.currentText() == "Regression Plot":
                self.grid = sns.FacetGrid(data=self.loader.long_data, **constructor)
                self.mainfig = self.grid.fig
            else:
                self.grid = sns.FacetGrid(data=self.loader.long_data,
                                          palette=self.fig_manager.cmb_palette.currentText(),
                                          **constructor)
                self.mainfig = self.grid.fig

        self.canvas = FigureCanvas(self.mainfig)
        self.nav = NavigationToolbar(self.canvas, self.canvas)
        self.mainlay.addWidget(self.nav)
        self.mainlay.addWidget(self.canvas)

    # Called by updates from the Figure Parameters
    @Slot()
    def plotupdate_facetgrid_figurecall_plot(self):
        # First clear the canvas
        self.canvas_full_clear()
        # Then generate the new canvas with the updated facetgrid and figure
        self.canvas_generate("Facet Grid Figure Call")
        # Then see if the axes can be updated, and if so, update them
        self.plotupdate_facetgrid_axes()


    # Called by updates from the Axes Parameters
    @Slot()
    def plotupdate_facetgrid_axescall_plot(self):
        self.plotupdate_facetgrid_axes()


    # Convenience function
    def plotupdate_facetgrid_axes(self):
        # If x or y axes arguments are still in their non-callable form (i.e string type), don't proceed
        if any([isinstance(self.fig_manager.axes_arg_x, str), isinstance(self.fig_manager.axes_arg_y, str)]):
            return

        # If x and y axes arguments are blank, don't proceed
        if any([self.fig_manager.axes_arg_x() == '', self.fig_manager.axes_arg_y() == '']):
            return

        # Account for a user typing the palette name, accidentally triggering this function, don't proceed
        if any([not self.fig_manager.cmb_palette,
                self.fig_manager.cmb_palette.currentText() not in self.fig_manager.palettenames]):
            return

        # Clear all axes
        axes = self.grid.axes
        for ax in axes.flat:
            ax.clear()

        # Otherwise, proceed
        func = self.fig_manager.plotting_func
        x, y, hue = self.fig_manager.axes_arg_x(), self.fig_manager.axes_arg_y(), self.fig_manager.axes_arg_hue()

        axes_constructor = {}
        for kwarg, call in self.fig_manager.axes_kwargs.items():
            if call() == "":
                axes_constructor[kwarg] = None
            else:
                axes_constructor[kwarg] = call()

        # Account for a user not selecting a palette
        if axes_constructor['palette'] in ['', "None", "Default Blue", "No Palette"]:
            axes_constructor['palette'] = None

        print("AXES CONSTRUCTOR:")
        pprint(axes_constructor)
        if hue == '':
            self.grid = self.grid.map(func, x, y, data=self.loader.long_data, **axes_constructor)
        else:
            self.grid = self.grid.map(func, x, y, hue, data=self.loader.long_data, **axes_constructor)

        # Redraw the legend
        self.plotupdate_facetgrid_legend()

    # Called by updates to the legend parameters of an axes within a FacetGrid
    def plotupdate_facetgrid_legend(self):
        if all([len(self.fig_manager.legend_kwargs) > 0,  # kwargs for the legend function must exist
                self.fig_manager.chk_showlegend.isChecked(),  # the option to show the legend must be there
                self.fig_manager.axes_arg_hue() != ''  # the hue argument cannot be empty
                ]):
            plt.legend(**self.fig_manager.legend_kwargs)
        else:
            plt.legend("", frameon=False)

        # This contains the canvas.draw() call within it
        self.plotupdate_padding()
