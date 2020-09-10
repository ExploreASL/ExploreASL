from PySide2.QtWidgets import *
from PySide2.QtGui import *
from PySide2.QtCore import *
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar
from matplotlib.figure import Figure
import seaborn as sns
from xASL_GUI_HelperClasses import DandD_FileExplorer2LineEdit, DandD_Graphing_ListWidget2LineEdit
from xASL_GUI_HelperFuncs import connect_widget_to_signal
import json
import sys
import os
from pprint import pprint


# noinspection PyCallingNonCallable
class xASL_GUI_FacetGridOrganizer(QWidget):
    """
    Class for generating the widgets that will go into the Figure Parms and Axes Parms of the main PostProc.
    Also, prepares arguments to feed into the constructor of the plotting. The goal is to keep the parameter generation
    vs plotting and data loading separate
    """
    change_axesparms_widget = Signal()  # This signal is for updating the mainwidget to change the Axes Parameters tab
    change_figparms_updateplot = Signal()  # This signal is for updating the plot due to a change in a figure parms
    change_axesparms_updateplot = Signal()  # This signal is for updating the plot due to a change in an axes parms
    change_legendparms_updateplot = Signal()  # This signal is for updating the plot due to a change in legend parms

    def __init__(self, parent):
        super(xASL_GUI_FacetGridOrganizer, self).__init__(parent=parent)
        self.parent_cw = parent  # Parent at the time of initialization, which will be the PostProc Widget
        self.manager_type = "Facet Grid"

        # Prepare arguments for legend widget
        self.legend_widget = xASL_GUI_FacetLegend(self)
        self.legend_widget.process_legendkwargs.connect(self.legend_kwarg_change_updateplot)
        self.legend_kwargs = {}

        # Prepare the figure parameters
        self.UI_Setup_FigureParms()

        # Prepare misc attributes and default arguments
        self.all_dtypes = ["object", "category", "bool", "float", "float16", "float32", "float64", "int", "int8",
                           "int16", "int32", "int64"]
        self.numeric_dtypes = ["float", "float16", "float32", "float64", "int", "int8", "int16", "int32", "int64"]
        self.palettenames = ["None", "Default Blue", "No Palette", 'Accent', 'Accent_r', 'Blues', 'Blues_r', 'BrBG',
                             'BrBG_r', 'BuGn', 'BuGn_r', 'BuPu',
                             'BuPu_r', 'CMRmap', 'CMRmap_r', 'Dark2', 'Dark2_r', 'GnBu', 'GnBu_r', 'Greens', 'Greens_r',
                             'Greys', 'Greys_r', 'OrRd', 'OrRd_r', 'Oranges', 'Oranges_r', 'PRGn', 'PRGn_r', 'Paired',
                             'Paired_r', 'Pastel1', 'Pastel1_r', 'Pastel2', 'Pastel2_r', 'PiYG', 'PiYG_r', 'PuBu',
                             'PuBuGn', 'PuBuGn_r', 'PuBu_r', 'PuOr', 'PuOr_r', 'PuRd', 'PuRd_r', 'Purples', 'Purples_r',
                             'RdBu', 'RdBu_r', 'RdGy', 'RdGy_r', 'RdPu', 'RdPu_r', 'RdYlBu', 'RdYlBu_r', 'RdYlGn',
                             'RdYlGn_r', 'Reds', 'Reds_r', 'Set1', 'Set1_r', 'Set2', 'Set2_r', 'Set3', 'Set3_r',
                             'Spectral', 'Spectral_r', 'Wistia', 'Wistia_r', 'YlGn', 'YlGnBu', 'YlGnBu_r', 'YlGn_r',
                             'YlOrBr', 'YlOrBr_r', 'YlOrRd', 'YlOrRd_r', 'afmhot', 'afmhot_r', 'autumn', 'autumn_r',
                             'binary', 'binary_r', 'bone', 'bone_r', 'brg', 'brg_r', 'bwr', 'bwr_r', 'cividis',
                             'cividis_r', 'cool', 'cool_r', 'coolwarm', 'coolwarm_r', 'copper', 'copper_r', 'cubehelix',
                             'cubehelix_r', 'flag', 'flag_r', 'gist_earth', 'gist_earth_r', 'gist_gray', 'gist_gray_r',
                             'gist_heat', 'gist_heat_r', 'gist_ncar', 'gist_ncar_r', 'gist_rainbow', 'gist_rainbow_r',
                             'gist_stern', 'gist_stern_r', 'gist_yarg', 'gist_yarg_r', 'gnuplot', 'gnuplot2',
                             'gnuplot2_r', 'gnuplot_r', 'gray', 'gray_r', 'hot', 'hot_r', 'hsv', 'hsv_r', 'icefire',
                             'icefire_r', 'inferno', 'inferno_r', 'jet', 'jet_r', 'magma', 'magma_r', 'mako', 'mako_r',
                             'nipy_spectral', 'nipy_spectral_r', 'ocean', 'ocean_r', 'pink', 'pink_r', 'plasma',
                             'plasma_r', 'prism', 'prism_r', 'rainbow', 'rainbow_r', 'rocket', 'rocket_r', 'seismic',
                             'seismic_r', 'spring', 'spring_r', 'summer', 'summer_r', 'tab10', 'tab10_r', 'tab20',
                             'tab20_r', 'tab20b', 'tab20b_r', 'tab20c', 'tab20c_r', 'terrain', 'terrain_r', 'twilight',
                             'twilight_r', 'twilight_shifted', 'twilight_shifted_r', 'viridis', 'viridis_r', 'vlag',
                             'vlag_r', 'winter', 'winter_r']
        self.default_palette_idx = self.palettenames.index("Set1")

        # Some default values for important variables
        self.axes_arg_x = ''  # Set to '' instead of None because the latter is not supported by lineedits
        self.axes_arg_y = ''
        self.axes_arg_hue = ''

    def sendSignal_figparms_updateplot(self):
        self.change_figparms_updateplot.emit()

    def sendSignal_axesparms_updateplot(self):
        self.change_axesparms_updateplot.emit()

    @Slot(dict)
    # Slot to process the legend_kwargs and then inform the main widget to update the canvas to reflect the change
    def legend_kwarg_change_updateplot(self, legend_kwargs: dict):
        print("legend_kwarg_change_updateplot received a signal")
        self.legend_kwargs = {}
        if "bbox_x" in legend_kwargs.keys():
            bbox_x = legend_kwargs["bbox_x"]()
            bbox_y = legend_kwargs["bbox_y"]()
            del legend_kwargs["bbox_x"]
            del legend_kwargs["bbox_y"]
            self.legend_kwargs["bbox_to_anchor"] = (bbox_x, bbox_y)
        for key, call in legend_kwargs.items():
            self.legend_kwargs[key] = call()

        self.change_legendparms_updateplot.emit()

    # Sets up the widgets that will be contained within the "Figure Parameters" Tab
    def UI_Setup_FigureParms(self):
        self.le_row = DandD_Graphing_ListWidget2LineEdit(self.parent_cw, ["object", "category"])
        self.le_row.setPlaceholderText("Drag & Drop variable to define facet rows")
        self.le_col = DandD_Graphing_ListWidget2LineEdit(self.parent_cw, ["object", "category"])
        self.le_col.setPlaceholderText("Drag & Drop variable to define facet columns")
        self.chk_sharex = QCheckBox(checked=True)
        self.chk_sharey = QCheckBox(checked=True)
        self.spinbox_facetheight = QDoubleSpinBox(maximum=10, minimum=1, value=5, singleStep=0.1)
        self.spinbox_aspect = QDoubleSpinBox(maximum=2, minimum=0.5, value=1, singleStep=0.1)
        self.chk_legend_out = QCheckBox(checked=True)
        self.chk_despine = QCheckBox(checked=True)
        self.chk_margin_titles = QCheckBox(checked=True)

        self.fig_kwargs = {"row": self.le_row.text,
                           "col": self.le_col.text,
                           "sharex": self.chk_sharex.isChecked,
                           "sharey": self.chk_sharey.isChecked,
                           "height": self.spinbox_facetheight.value,
                           "aspect": self.spinbox_aspect.value,
                           "legend_out": self.chk_legend_out.isChecked,
                           "despine": self.chk_despine.isChecked,
                           "margin_titles": self.chk_margin_titles.isChecked
                           }
        self.cont_figparms = QWidget()
        self.formlay_figparms = QFormLayout(self.cont_figparms)
        for description, widget in zip(["Row Grouping Variable", "Column Grouping Variable", "Share x axes?",
                                        "Share y axes?", "Height of each facet (inches)", "Aspect Ratio",
                                        "Draw legend outside grid?", "Remove top and right spines?", "Margin titles"],
                                       [self.le_row, self.le_col, self.chk_sharex, self.chk_sharey,
                                        self.spinbox_facetheight, self.spinbox_aspect,
                                        self.chk_legend_out, self.chk_despine, self.chk_margin_titles]):
            self.formlay_figparms.addRow(description, widget)
            connect_widget_to_signal(widget, self.sendSignal_figparms_updateplot)

        self.cmb_axestype = QComboBox()
        self.cmb_axestype.addItems(["Select a plot type", "Point Plot", "Bar Plot", "Strip Plot", "Swarm Plot",
                                    "Box Plot", "Violin Plot", "Boxen Plot", "Scatter Plot", "Line Plot"])
        self.cmb_axestype.currentTextChanged.connect(self.UI_Setup_AxesParms)
        self.formlay_figparms.addRow(self.cmb_axestype)

    # Sets up the widgets that will be contained within the "Axes Parameters" Tab
    def UI_Setup_AxesParms(self, plot_type):
        # Start off by removing the previous axes
        self.parent_cw.clear_axesparms()
        print(f"Selected {plot_type} as the Axes Type")

        # These are always a given
        if self.cmb_axestype.currentText() in ["Point Plot", "Bar Plot", "Strip Plot", "Swarm Plot",
                                               "Box Plot", "Violin Plot", "Boxen Plot"]:
            self.le_x = DandD_Graphing_ListWidget2LineEdit(self.parent_cw, self.all_dtypes)
            self.le_x.setPlaceholderText("Drag & Drop the X-axis Variable")
            self.le_y = DandD_Graphing_ListWidget2LineEdit(self.parent_cw, self.all_dtypes)
            self.le_y.setPlaceholderText("Drag & Drop the Y-axis Variable")

        elif self.cmb_axestype.currentText() in ["Scatter Plot", "Line Plot",
                                                 "Regression Plot", "Residuals Plot"]:
            self.le_x = DandD_Graphing_ListWidget2LineEdit(self.parent_cw, self.numeric_dtypes)
            self.le_x.setPlaceholderText("Drag & Drop the X-axis Variable")
            self.le_y = DandD_Graphing_ListWidget2LineEdit(self.parent_cw, self.numeric_dtypes)
            self.le_y.setPlaceholderText("Drag & Drop the Y-axis Variable")

        else:
            print("THIS SHOULD NEVER PRINT!!!!!")
            return

        self.le_hue = DandD_Graphing_ListWidget2LineEdit(self.parent_cw, ["object", "category"])
        self.le_hue.setPlaceholderText("Drag & Drop the Hue Grouping Variable")
        self.axes_arg_x = self.le_x.text
        self.axes_arg_y = self.le_y.text
        self.axes_arg_hue = self.le_hue.text

        self.le_x.textChanged.connect(self.sendSignal_axesparms_updateplot)
        self.le_y.textChanged.connect(self.sendSignal_axesparms_updateplot)
        self.le_hue.textChanged.connect(self.sendSignal_axesparms_updateplot)

        if plot_type == "Point Plot":
            self.plotting_func = sns.pointplot
            self.ci = QDoubleSpinBox(maximum=100, minimum=0, value=95, singleStep=1)
            self.dodge = QCheckBox(checked=True)
            self.join = QCheckBox(checked=True)
            self.errwidth = QDoubleSpinBox(maximum=2, minimum=0, value=2, singleStep=0.05)
            self.capsize = QDoubleSpinBox(maximum=1, minimum=0, value=0.05, singleStep=0.005)
            self.cmb_palette = self.UI_Setup_PaletteCombobox()
            self.axes_kwargs = {"ci": self.ci.value,
                                "dodge": self.dodge.isChecked,
                                "join": self.join.isChecked,
                                "errwidth": self.errwidth.value,
                                "capsize": self.capsize.value,
                                "palette": self.cmb_palette.currentText
                                }
            self.cont_axesparms = QWidget()
            self.formlay_axesparms = QFormLayout(self.cont_axesparms)
            for description, widget in zip(["X Axis Variable", "Y Axis variable", "Hue Grouping Variable",
                                            "Confidence interval", "Offset points?", "Join points?",
                                            "Errorbar thickness", "Cap width", "Palette"],
                                           [self.le_x, self.le_y, self.le_hue, self.ci, self.dodge, self.join,
                                            self.errwidth, self.capsize, self.cmb_palette]):
                self.formlay_axesparms.addRow(description, widget)
                connect_widget_to_signal(widget, self.sendSignal_axesparms_updateplot)

            # Send the signal for this widget group to be added to the main Axes Parameters tab
            self.change_axesparms_widget.emit()

        elif plot_type == "Bar Plot":
            self.plotting_func = sns.barplot
            self.ci = QDoubleSpinBox(maximum=100, minimum=0, value=95, singleStep=1)
            self.dodge = QCheckBox(checked=True)
            self.errwidth = QDoubleSpinBox(maximum=3, minimum=0, value=1.55, singleStep=0.05)
            self.capsize = QDoubleSpinBox(maximum=1, minimum=0, value=0.04, singleStep=0.01)
            self.cmb_palette = QComboBox()
            self.cmb_palette.addItems(self.palettenames)
            self.axes_kwargs = {"ci": self.ci.value,
                                "dodge": self.dodge.isChecked,
                                "errwidth": self.errwidth.value,
                                "capsize": self.capsize.value,
                                "palette": self.cmb_palette.currentText
                                }
            self.cont_axesparms = QWidget()
            self.formlay_axesparms = QFormLayout(self.cont_axesparms)
            for description, widget in zip(["X Axis Variable", "Y Axis variable", "Hue Grouping Variable",
                                            "Confidence interval", "Offset when hue nesting?",
                                            "Errorbar thickness", "Cap width", "Palette"],
                                           [self.le_x, self.le_y, self.le_hue, self.ci, self.dodge,
                                            self.errwidth, self.capsize, self.cmb_palette]):
                self.formlay_axesparms.addRow(description, widget)
                connect_widget_to_signal(widget, self.sendSignal_axesparms_updateplot)

            # Send the signal for this widget group to be added to the main Axes Parameters tab
            self.change_axesparms_widget.emit()

        elif plot_type == "Strip Plot" or plot_type == "Swarm Plot":
            if plot_type == "Strip Plot":
                self.plotting_func = sns.stripplot
            else:
                self.plotting_func = sns.swarmplot
            self.dodge = QCheckBox(checked=True)
            self.size = QDoubleSpinBox(maximum=10, minimum=1, value=4.5, singleStep=0.1)
            self.linewidth = QDoubleSpinBox(maximum=2, minimum=0, value=0.1, singleStep=0.1)
            self.cmb_palette = self.UI_Setup_PaletteCombobox()
            self.axes_kwargs = {"dodge": self.dodge.isChecked,
                                "size": self.size.value,
                                "linewidth": self.linewidth.value,
                                "palette": self.cmb_palette.currentText
                                }
            self.cont_axesparms = QWidget()
            self.formlay_axesparms = QFormLayout(self.cont_axesparms)
            for description, widget in zip(["X Axis Variable", "Y Axis variable", "Hue Grouping Variable",
                                            "Separate groupings when hue nesting?", "Marker size", "Marker edge width",
                                            "Palette"],
                                           [self.le_x, self.le_y, self.le_hue, self.dodge, self.size, self.linewidth,
                                            self.cmb_palette]):
                self.formlay_axesparms.addRow(description, widget)
                connect_widget_to_signal(widget, self.sendSignal_axesparms_updateplot)

            # Send the signal for this widget group to be added to the main Axes Parameters tab
            self.change_axesparms_widget.emit()

        elif plot_type == "Box Plot":
            self.plotting_func = sns.boxplot
            self.barwidth = QDoubleSpinBox(maximum=1, minimum=0, value=0.7, singleStep=0.05)
            self.dodge = QCheckBox(checked=True)
            self.fliersize = QDoubleSpinBox(maximum=10, minimum=0, value=5, singleStep=0.25)
            self.linewidth = QDoubleSpinBox(maximum=5, minimum=0, value=1.55, singleStep=0.05)
            self.whis = QDoubleSpinBox(maximum=3, minimum=0, value=1.5, singleStep=0.1)
            self.cmb_palette = self.UI_Setup_PaletteCombobox()
            self.axes_kwargs = {"width": self.barwidth.value,
                                "dodge": self.dodge.isChecked,
                                "fliersize": self.fliersize.value,
                                "linewidth": self.linewidth.value,
                                "palette": self.cmb_palette.currentText,
                                "whis": self.whis.value
                                }
            self.cont_axesparms = QWidget()
            self.formlay_axesparms = QFormLayout(self.cont_axesparms)
            for description, widget in zip(["X Axis Variable", "Y Axis variable", "Hue Grouping Variable",
                                            "Bar width when not hue nesting?", "Shift bars when hue nesting?",
                                            "Marker size of outlier observations", "Box linewidth",
                                            "Palette",
                                            "Fraction of IQR to identify inliers"],
                                           [self.le_x, self.le_y, self.le_hue, self.barwidth, self.dodge,
                                            self.fliersize, self.linewidth, self.cmb_palette, self.whis]):
                self.formlay_axesparms.addRow(description, widget)
                connect_widget_to_signal(widget, self.sendSignal_axesparms_updateplot)

            # Send the signal for this widget group to be added to the main Axes Parameters tab
            self.change_axesparms_widget.emit()

        elif plot_type == "Violin Plot":
            self.plotting_func = sns.violinplot
            self.kernalbwalgo = QComboBox()
            self.kernalbwalgo.addItems(["scott", "silverman"])
            self.scale = QComboBox()
            self.scale.addItems(["area", "count", "width"])
            self.scale_hue = QCheckBox(checked=True)
            self.dodge = QCheckBox(checked=True)
            self.linewidth = QDoubleSpinBox(maximum=5, minimum=0, value=1.8, singleStep=0.1)
            self.cmb_palette = self.UI_Setup_PaletteCombobox()
            self.axes_kwargs = {"bw": self.kernalbwalgo.currentText,
                                "scale": self.scale.currentText,
                                "scale_hue": self.scale_hue.isChecked,
                                "dodge": self.dodge.isChecked,
                                "linewidth": self.linewidth.value,
                                "palette": self.cmb_palette.currentText
                                }
            self.cont_axesparms = QWidget()
            self.formlay_axesparms = QFormLayout(self.cont_axesparms)
            for description, widget in zip(["X Axis Variable", "Y Axis variable", "Hue Grouping Variable",
                                            "Kernel algorithm", "Basis of scaling width",
                                            "Recalculate scaling when hue nesting?", "Shift plots when hue nesting?",
                                            "Linewidth of outlines", "Palette"],
                                           [self.le_x, self.le_y, self.le_hue, self.kernalbwalgo, self.scale,
                                            self.scale_hue, self.dodge, self.linewidth, self.cmb_palette]):
                self.formlay_axesparms.addRow(description, widget)
                connect_widget_to_signal(widget, self.sendSignal_axesparms_updateplot)

            # Send the signal for this widget group to be added to the main Axes Parameters tab
            self.change_axesparms_widget.emit()

        elif plot_type == "Boxen Plot":
            self.plotting_func = sns.boxenplot
            self.barwidth = QDoubleSpinBox(maximum=1, minimum=0, value=0.8, singleStep=0.05)
            self.dodge = QCheckBox(checked=True)
            self.linewidth = QDoubleSpinBox(maximum=5, minimum=0, value=0.7, singleStep=0.1)
            self.outlier_prop = QDoubleSpinBox(maximum=1, minimum=0, value=0.007, singleStep=0.001)
            self.show_outliers = QCheckBox(checked=True)
            self.cmb_palette = self.UI_Setup_PaletteCombobox()
            self.axes_kwargs = {"width": self.barwidth.value,
                                "dodge": self.dodge.isChecked,
                                "linewidth": self.linewidth.value,
                                "outlier_prop": self.outlier_prop.value,
                                "showfliers": self.show_outliers.isChecked,
                                "palette": self.cmb_palette.currentText
                                }
            self.cont_axesparms = QWidget()
            self.formlay_axesparms = QFormLayout(self.cont_axesparms)
            for description, widget in zip(["X Axis Variable", "Y Axis variable", "Hue Grouping Variable",
                                            "Relative bar width", "Shift plots when hue nesting?",
                                            "Linewidth of the boxes", "Expected proportion of outliers",
                                            "Show outliers?", "Palette"],
                                           [self.le_x, self.le_y, self.le_hue, self.barwidth, self.dodge,
                                            self.linewidth, self.outlier_prop, self.show_outliers, self.cmb_palette]):
                self.formlay_axesparms.addRow(description, widget)
                connect_widget_to_signal(widget, self.sendSignal_axesparms_updateplot)

            # Send the signal for this widget group to be added to the main Axes Parameters tab
            self.change_axesparms_widget.emit()

        elif plot_type == "Scatter Plot":
            self.plotting_func = sns.scatterplot
            self.size_grouper = DandD_Graphing_ListWidget2LineEdit(self.parent_cw, self.all_dtypes)
            self.size_grouper.setPlaceholderText("(Optional) Drag & Drop the Markersize-Grouping Variable")
            self.style_grouper = DandD_Graphing_ListWidget2LineEdit(self.parent_cw, ["object", "category"])
            self.style_grouper.setPlaceholderText("(Optional) Drag & Drop the Markerstyle-Grouping Variable")
            self.cmb_palette = self.UI_Setup_PaletteCombobox()
            self.spinbox_markersize = QDoubleSpinBox(maximum=100, minimum=0, value=40, singleStep=1)
            self.axes_kwargs = {"size": self.size_grouper.text,
                                "style": self.style_grouper.text,
                                "palette": self.cmb_palette.currentText,
                                "s": self.spinbox_markersize.value
                                }
            self.cont_axesparms = QWidget()
            self.formlay_axesparms = QFormLayout(self.cont_axesparms)
            for description, widget in zip(["X Axis Variable", "Y Axis variable", "Hue Grouping Variable",
                                            "Size Grouping Variable", "Marker Style Grouping Variable", "Marker size",
                                            "Palette"],
                                           [self.le_x, self.le_y, self.le_hue, self.size_grouper, self.style_grouper,
                                            self.spinbox_markersize, self.cmb_palette]):
                self.formlay_axesparms.addRow(description, widget)
                connect_widget_to_signal(widget, self.sendSignal_axesparms_updateplot)

            # Send the signal for this widget group to be added to the main Axes Parameters tab
            self.change_axesparms_widget.emit()

        elif plot_type == "Line Plot":
            self.plotting_func = sns.lineplot
            self.size_grouper = DandD_Graphing_ListWidget2LineEdit(self.parent_cw, self.all_dtypes)
            self.size_grouper.setPlaceholderText("(Optional) Drag & Drop the Markersize-Grouping Variable")
            self.style_grouper = DandD_Graphing_ListWidget2LineEdit(self.parent_cw, ["object", "category"])
            self.style_grouper.setPlaceholderText("(Optional) Drag & Drop the Markerstyle-Grouping Variable")
            self.sort = QCheckBox(checked=True)
            self.dashes = QCheckBox(checked=False)
            self.err_style = QComboBox()
            self.err_style.addItems(["band", "bars"])
            self.cmb_palette = self.UI_Setup_PaletteCombobox()
            self.axes_kwargs = {"size": self.size_grouper.text,
                                "style": self.style_grouper.text,
                                "sort": self.sort.isChecked,
                                "dashes": self.dashes.isChecked,
                                "err_style": self.err_style.currentText,
                                "palette": self.cmb_palette.currentText
                                }
            self.cont_axesparms = QWidget()
            self.formlay_axesparms = QFormLayout(self.cont_axesparms)
            for description, widget in zip(["X Axis Variable", "Y Axis variable", "Hue Grouping Variable",
                                            "Size Grouping Variable", "Marker Style Grouping Variable",
                                            "Pre-sort data first?",
                                            "Use dashes instead of solid lines?", "Style of errorbars", "Palette"],
                                           [self.le_x, self.le_y, self.le_hue, self.size_grouper, self.style_grouper,
                                            self.sort, self.dashes, self.err_style, self.cmb_palette]):
                self.formlay_axesparms.addRow(description, widget)
                connect_widget_to_signal(widget, self.sendSignal_axesparms_updateplot)

            # Send the signal for this widget group to be added to the main Axes Parameters tab
            self.change_axesparms_widget.emit()

        else:
            print("THIS SHOULD NEVER PRINT. AN IMPOSSIBLE AXES TYPE WAS SELECTED")

        # Setup the wigets associated with connecting to the legend parameters each time
        self.UI_Setup_Btn2LegendParms()

        # Regardless of what type of plot was selected, the current axes should be cleared
        axes = self.parent_cw.grid.axes
        for ax in axes.flat:
            ax.clear()
        self.parent_cw.canvas.draw()

    # Convenience function to generate the combobox that will respond to palette changes
    def UI_Setup_PaletteCombobox(self):
        cmb_palette = QComboBox()
        cmb_palette.addItems(self.palettenames)
        cmb_palette.setMaxVisibleItems(10)
        cmb_palette.setEditable(True)
        cmb_palette.setCurrentIndex(self.default_palette_idx)
        if len(self.parent_cw.loader.long_data) > 2000:
            cmb_palette.setAutoCompletion(False)
        return cmb_palette

    # Convenience function to generate the axes-level widgets that will grant access to legend parameters
    def UI_Setup_Btn2LegendParms(self):
        self.chk_showlegend = QCheckBox(clicked=self.legend_widget.sendSignal_legendparms_updateplot)
        self.btn_legendparms = QPushButton("Change Legend Parameters", clicked=self.legend_widget.show)
        self.formlay_axesparms.addRow("Show legend in figure?", self.chk_showlegend)
        self.formlay_axesparms.addRow(self.btn_legendparms)


class xASL_GUI_FacetLegend(QWidget):
    process_legendkwargs = Signal(dict)

    def __init__(self, parent=None):
        super().__init__(parent)
        self.parent = parent  # The Facet Grid Manager
        self.setWindowFlag(Qt.Window)
        self.setWindowTitle("FacetGrid Legend Settings")
        self.setMinimumSize(300, 480)
        self.legend_kwargs = {}
        self.formlay_legend = QFormLayout(self)
        self.cmb_location = QComboBox()
        self.cmb_location.addItems(["best", "upper right", "upper left", "lower left", "lower right", "center left",
                                    "center right", "lower center", "upper center"])
        self.chk_manual_loc = QCheckBox(checked=False)
        self.spinbox_legend_x = QDoubleSpinBox(maximum=2, minimum=0, value=1.05, singleStep=0.01)
        self.spinbox_legend_y = QDoubleSpinBox(maximum=20, minimum=0, value=1, singleStep=0.01)
        self.cmb_fontsize = QComboBox()
        self.cmb_fontsize.addItems(["xx-small", "x-small", "small", "medium", "large", "x-large", "xx-large"])
        self.cmb_fontsize.setCurrentIndex(3)
        self.chk_frameon = QCheckBox(checked=False)
        self.chk_fancybox = QCheckBox(checked=False)
        self.chk_shadow = QCheckBox(checked=False)
        self.spinbox_markersize = QDoubleSpinBox(maximum=2, minimum=0, value=1, singleStep=0.01)
        self.chk_markerfirst = QCheckBox(checked=True)

        for widget, description in zip([self.cmb_location, self.chk_manual_loc, self.spinbox_legend_x,
                                        self.spinbox_legend_y, self.cmb_fontsize, self.chk_frameon, self.chk_fancybox,
                                        self.chk_shadow, self.spinbox_markersize, self.chk_markerfirst],
                                       ["Legend location", "Manually set location?", "Manual X position",
                                        "Manual Y Position", "Legend fontsize", "Legend frame is on?",
                                        "Legend frame has rounded corners?", "Legend frame has shadow?",
                                        "Marker size", "Marker first before label?"]):
            self.formlay_legend.addRow(description, widget)
            connect_widget_to_signal(widget, self.sendSignal_legendparms_updateplot)

    def sendSignal_legendparms_updateplot(self):
        self.prepare_legend_kwargs()
        self.process_legendkwargs.emit(self.legend_kwargs)

    def prepare_legend_kwargs(self):
        self.legend_kwargs = {
            "loc": self.cmb_location.currentText,
            "fontsize": self.cmb_fontsize.currentText,
            "frameon": self.chk_frameon.isChecked,
            "fancybox": self.chk_fancybox.isChecked,
            "shadow": self.chk_shadow.isChecked,
            "markerscale": self.spinbox_markersize.value,
            "markerfirst": self.chk_markerfirst.isChecked
        }
        if self.chk_manual_loc.isChecked():
            self.legend_kwargs["bbox_x"] = self.spinbox_legend_x.value
            self.legend_kwargs["bbox_y"] = self.spinbox_legend_y.value
