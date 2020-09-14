from PySide2.QtWidgets import *
from PySide2.QtGui import *
from PySide2.QtCore import *
import seaborn as sns
from ExploreASL_GUI.xASL_GUI_HelperClasses import DandD_Graphing_ListWidget2LineEdit
from ExploreASL_GUI.xASL_GUI_HelperFuncs import connect_widget_to_signal
import json
import os
from ExploreASL_GUI.xASL_GUI_Graph_FacetArtist import xASL_GUI_FacetArtist
from ExploreASL_GUI.xASL_GUI_Graph_FacetLabels import xASL_GUI_FacetLabels
# from pprint import pprint


class xASL_GUI_FacetManager(QWidget):
    signal_axesparms_updateplot = Signal()
    signal_figparms_updateplot = Signal()
    signal_paddingparms_updateplot = Signal()

    def __init__(self, parent):
        super().__init__(parent=parent)
        self.parent_cw = parent
        self.artist: xASL_GUI_FacetArtist = self.parent_cw.fig_artist
        self.mainlay = QVBoxLayout(self)

        ##########################################
        # Widget Classes for more complex settings
        ##########################################
        # Prepare arguments for legend widget
        self.legend_widget = xASL_GUI_FacetLegend(self)
        self.labels_widget = xASL_GUI_FacetLabels(self)

        # Set up variables that the artist will reference when updating the plot
        self.axes_arg_x = ''  # Set to '' instead of None because the latter is not supported by lineedits
        self.axes_arg_y = ''
        self.axes_arg_hue = ''
        self.legend_kwargs = {}
        self.padding_kwargs = {}
        self.axes_widget = None

        with open(os.path.join(self.parent_cw.config["ProjectDir"], "JSON_LOGIC",
                               "GraphingParameters.json")) as freader:
            parms = json.load(freader)
            self.all_dtypes = parms["all_dtypes"]
            self.numeric_dtypes = parms["numeric_dtypes"]
            self.categorical_dtypes = ["object", "category"]
            self.palettenames = parms["palettenames"]
            self.default_palette_idx = self.palettenames.index("Set1")

        # Setup the widget
        self.UI_Setup_Tabs()
        self.UI_Setup_CommonParameters()
        self.UI_Setup_FigureParms()

    def UI_Setup_Tabs(self):
        self.tab_pltsettings = QTabWidget()
        self.cont_commonparms, self.cont_figparms, self.cont_axesparms = QWidget(), QWidget(), QWidget()
        self.tab_pltsettings.addTab(self.cont_commonparms, "Common Parameters")
        self.tab_pltsettings.addTab(self.cont_figparms, "Facet-Level Parameters")
        self.tab_pltsettings.addTab(self.cont_axesparms, "Plot-Level Parameters")
        self.tab_pltsettings.setMaximumHeight(400)

        # Also define the formlayouts here and in advance
        self.formlay_commonparms = QFormLayout(self.cont_commonparms)
        self.formlay_figparms = QFormLayout(self.cont_figparms)
        self.vlay_axesparms = QVBoxLayout(self.cont_axesparms)
        self.mainlay.addWidget(self.tab_pltsettings)

    # Set up connections with an artist
    def UI_Setup_ConnectManager2Artist(self):
        self.artist: xASL_GUI_FacetArtist = self.parent_cw.fig_artist
        if self.artist is not None:
            self.cmb_axestype.currentTextChanged.connect(self.UI_Setup_AxesParms)
            # Connect Facet Signals to the artist's Slots
            self.signal_figparms_updateplot.connect(self.artist.plotupdate_figurecall)
            self.signal_axesparms_updateplot.connect(self.artist.plotupdate_axescall)
            self.signal_paddingparms_updateplot.connect(self.artist.plotupdate_paddingcall)
            # Connect the legend_widget's signals to the artist's Slots
            self.legend_widget.signal_legendcall_updateplot.connect(self.artist.plotupdate_legendcall)
            # Connect the label_widget's signals to the artist's Slots
            self.labels_widget.signal_xaxislabel_plotupdate.connect(self.artist.plotupdate_xlabel)
            self.labels_widget.signal_yaxislabel_plotupdate.connect(self.artist.plotupdate_ylabel)
            self.labels_widget.signal_title_plotupdate.connect(self.artist.plotupdate_title)

            # Perform some basic adjustments first as well, such as padding
            self.sendSignal_plotupdate_paddingcall()

        else:
            print("UI_Setup_Connections the artist had a value of None")

    def UI_Setup_CommonParameters(self):
        # Set up the padding
        self.spinbox_toppad = QDoubleSpinBox()
        self.spinbox_bottompad = QDoubleSpinBox()
        self.spinbox_leftpad = QDoubleSpinBox()
        self.spinbox_rightpad = QDoubleSpinBox()
        self.spinbox_wspace = QDoubleSpinBox()
        self.spinbox_hspace = QDoubleSpinBox()
        for description, widget, default in zip(["Top Padding", "Bottom Padding", "Left Padding", "Right Padding",
                                                 "Horizontal Inter-Facet Spacing", "Vertical Inter-Facet Spacing"],
                                                [self.spinbox_toppad, self.spinbox_bottompad, self.spinbox_leftpad,
                                                 self.spinbox_rightpad, self.spinbox_wspace, self.spinbox_hspace],
                                                [0.91, 0.057, 0.055,
                                                 0.95, 0.2, 0.2]):
            widget.setValue(default)
            widget.setRange(0, 1)
            widget.setSingleStep(0.01)
            widget.valueChanged.connect(self.sendSignal_plotupdate_paddingcall)
            self.formlay_commonparms.addRow(description, widget)

        # Set up the text options and connect the plotlabels signals to the appropriate function
        self.btn_showplotlabels = QPushButton("Show Plot Label Settings", clicked=self.labels_widget.show)
        self.formlay_commonparms.addRow(self.btn_showplotlabels)

    def UI_Setup_FigureParms(self):
        # Define Widgets
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

        # Generate Facet Figure-level Mappings
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

        # Attach defined widgets to the form layout
        for description, widget in zip(["Row Grouping Variable", "Column Grouping Variable", "Share x axes?",
                                        "Share y axes?", "Height of each facet (inches)", "Aspect Ratio",
                                        "Draw legend outside grid?", "Remove top and right spines?", "Margin titles"],
                                       [self.le_row, self.le_col, self.chk_sharex, self.chk_sharey,
                                        self.spinbox_facetheight, self.spinbox_aspect,
                                        self.chk_legend_out, self.chk_despine, self.chk_margin_titles]):
            self.formlay_figparms.addRow(description, widget)
            connect_widget_to_signal(widget, self.sendSignal_plotupdate_figurecall)

        # Finally add the combobox for selecting the axes type. The connect signal will be established later in
        # UI_Setup_ConnectManager2Artist to ensure that the artist exists
        self.cmb_axestype = QComboBox()
        self.cmb_axestype.addItems(["Select a plot type", "Point Plot", "Strip Plot", "Swarm Plot", "Box Plot",
                                    "Violin Plot", "Boxen Plot", "Scatter Plot"])  # "Bar Plot", "Line Plot"
        font = QFont()
        font.setPointSize(12)
        self.cmb_axestype.setFont(font)
        self.formlay_figparms.addRow(self.cmb_axestype)

    # Convenience Function - clears the Axes Parameters Tab
    # noinspection PyTypeChecker
    def clear_axesparms(self):
        print("Clearing Axes Parms")
        if self.axes_widget is not None:
            self.vlay_axesparms.removeWidget(self.subcont_axesparms)
            self.subcont_axesparms.setParent(None)
            del self.subcont_axesparms
            self.axes_widget = None

    # Sets up the widgets that will be contained within the "Axes Parameters" Tab
    def UI_Setup_AxesParms(self, plot_type):
        self.clear_axesparms()
        print(f"Selected {plot_type} as the Axes Type")

        # These are always a given
        if self.cmb_axestype.currentText() in ["Point Plot", "Bar Plot", "Strip Plot", "Swarm Plot",
                                               "Box Plot", "Violin Plot", "Boxen Plot"]:
            self.le_x = DandD_Graphing_ListWidget2LineEdit(self.parent_cw, self.categorical_dtypes)
            self.le_x.setPlaceholderText("Drag & Drop the X-axis Variable")
            self.le_y = DandD_Graphing_ListWidget2LineEdit(self.parent_cw, self.numeric_dtypes)
            self.le_y.setPlaceholderText("Drag & Drop the Y-axis Variable")

        elif self.cmb_axestype.currentText() in ["Scatter Plot", "Line Plot",
                                                 "Regression Plot", "Residuals Plot"]:
            self.le_x = DandD_Graphing_ListWidget2LineEdit(self.parent_cw, self.numeric_dtypes)
            self.le_x.setPlaceholderText("Drag & Drop the X-axis Variable")
            self.le_y = DandD_Graphing_ListWidget2LineEdit(self.parent_cw, self.numeric_dtypes)
            self.le_y.setPlaceholderText("Drag & Drop the Y-axis Variable")
        else:
            pass

        self.le_hue = DandD_Graphing_ListWidget2LineEdit(self.parent_cw, ["object", "category"])
        self.le_hue.setPlaceholderText("Drag & Drop the Hue Grouping Variable")
        self.axes_arg_x = self.le_x.text
        self.axes_arg_y = self.le_y.text
        self.axes_arg_hue = self.le_hue.text

        if plot_type == "Point Plot":
            # Define row widgets
            self.plotting_func = sns.pointplot
            self.ci = QDoubleSpinBox(maximum=100, minimum=0, value=95, singleStep=1)
            self.dodge = QCheckBox(checked=True)
            self.join = QCheckBox(checked=True)
            self.errwidth = QDoubleSpinBox(maximum=2, minimum=0, value=2, singleStep=0.05)
            self.capsize = QDoubleSpinBox(maximum=1, minimum=0, value=0.05, singleStep=0.005)
            self.cmb_palette = self.UI_Setup_PaletteCombobox()
            self.axes_kwargs = self.UI_Setup_AxesMappings(["ci", "dodge", "join", "errwidth", "capsize", "palette"],
                                                          [self.ci, self.dodge, self.join, self.errwidth, self.capsize,
                                                           self.cmb_palette])
            # Create the underlying widget container and form layout
            self.UI_Setup_AxesParms_Subcontainer()
            # Add widgets to form layout and connect to appropriate function/signal
            for description, widget in zip(["X Axis Variable", "Y Axis variable", "Hue Grouping Variable",
                                            "Confidence interval", "Offset points?", "Join points?",
                                            "Errorbar thickness", "Cap width", "Palette"],
                                           [self.le_x, self.le_y, self.le_hue, self.ci, self.dodge, self.join,
                                            self.errwidth, self.capsize, self.cmb_palette]):
                self.formlay_axesparms.addRow(description, widget)
                connect_widget_to_signal(widget, self.sendSignal_plotupdate_axescall)

        elif plot_type == "Strip Plot" or plot_type == "Swarm Plot":
            if plot_type == "Strip Plot":
                self.plotting_func = sns.stripplot
            else:
                self.plotting_func = sns.swarmplot
            self.dodge = QCheckBox(checked=True)
            self.size = QDoubleSpinBox(maximum=10, minimum=1, value=4.5, singleStep=0.1)
            self.linewidth = QDoubleSpinBox(maximum=2, minimum=0, value=0.1, singleStep=0.1)
            self.cmb_palette = self.UI_Setup_PaletteCombobox()
            self.axes_kwargs = self.UI_Setup_AxesMappings(["dodge", "size", "linewidth", "palette"],
                                                          [self.dodge, self.size, self.linewidth, self.cmb_palette])
            # Create the underlying widget container and form layout
            self.UI_Setup_AxesParms_Subcontainer()
            # Add widgets to form layout and connect to appropriate function/signal
            for description, widget in zip(["X Axis Variable", "Y Axis variable", "Hue Grouping Variable",
                                            "Separate groupings when hue nesting?", "Marker size", "Marker edge width",
                                            "Palette"],
                                           [self.le_x, self.le_y, self.le_hue, self.dodge, self.size, self.linewidth,
                                            self.cmb_palette]):
                self.formlay_axesparms.addRow(description, widget)
                connect_widget_to_signal(widget, self.sendSignal_plotupdate_axescall)

        elif plot_type == "Box Plot":
            self.plotting_func = sns.boxplot
            self.barwidth = QDoubleSpinBox(maximum=1, minimum=0, value=0.7, singleStep=0.05)
            self.dodge = QCheckBox(checked=True)
            self.fliersize = QDoubleSpinBox(maximum=10, minimum=0, value=5, singleStep=0.25)
            self.linewidth = QDoubleSpinBox(maximum=5, minimum=0, value=1.55, singleStep=0.05)
            self.whis = QDoubleSpinBox(maximum=3, minimum=0, value=1.5, singleStep=0.1)
            self.cmb_palette = self.UI_Setup_PaletteCombobox()
            self.axes_kwargs = self.UI_Setup_AxesMappings(["width", "dodge", "fliersize",
                                                           "linewidth", "palette", "whis"],
                                                          [self.barwidth, self.dodge, self.fliersize,
                                                           self.linewidth, self.cmb_palette, self.whis])
            # Create the underlying widget container and form layout
            self.UI_Setup_AxesParms_Subcontainer()
            # Add widgets to form layout and connect to appropriate function/signal
            for description, widget in zip(["X Axis Variable", "Y Axis variable", "Hue Grouping Variable",
                                            "Bar width when not hue nesting?", "Shift bars when hue nesting?",
                                            "Marker size of outlier observations", "Box linewidth",
                                            "Palette",
                                            "Fraction of IQR to identify inliers"],
                                           [self.le_x, self.le_y, self.le_hue, self.barwidth, self.dodge,
                                            self.fliersize, self.linewidth, self.cmb_palette, self.whis]):
                self.formlay_axesparms.addRow(description, widget)
                connect_widget_to_signal(widget, self.sendSignal_plotupdate_axescall)

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
            self.axes_kwargs = self.UI_Setup_AxesMappings(["bw", "scale", "scale_hue", "dodge", "linewidth", "palette"],
                                                          [self.kernalbwalgo, self.scale, self.scale_hue, self.dodge,
                                                           self.linewidth,
                                                           self.cmb_palette])
            # Create the underlying widget container and form layout
            self.UI_Setup_AxesParms_Subcontainer()
            # Add widgets to form layout and connect to appropriate function/signal
            for description, widget in zip(["X Axis Variable", "Y Axis variable", "Hue Grouping Variable",
                                            "Kernel algorithm", "Basis of scaling width",
                                            "Recalculate scaling when hue nesting?", "Shift plots when hue nesting?",
                                            "Linewidth of outlines", "Palette"],
                                           [self.le_x, self.le_y, self.le_hue, self.kernalbwalgo, self.scale,
                                            self.scale_hue, self.dodge, self.linewidth, self.cmb_palette]):
                self.formlay_axesparms.addRow(description, widget)
                connect_widget_to_signal(widget, self.sendSignal_plotupdate_axescall)

        elif plot_type == "Boxen Plot":
            self.plotting_func = sns.boxenplot
            self.barwidth = QDoubleSpinBox(maximum=1, minimum=0, value=0.8, singleStep=0.05)
            self.dodge = QCheckBox(checked=True)
            self.linewidth = QDoubleSpinBox(maximum=5, minimum=0, value=0.7, singleStep=0.1)
            self.outlier_prop = QDoubleSpinBox(maximum=1, minimum=0, value=0.007, singleStep=0.001)
            self.show_outliers = QCheckBox(checked=True)
            self.cmb_palette = self.UI_Setup_PaletteCombobox()
            self.axes_kwargs = self.UI_Setup_AxesMappings(["width", "dodge", "linewidth", "outlier_prop",
                                                           "showfliers", "palette"],
                                                          [self.barwidth, self.dodge, self.linewidth, self.outlier_prop,
                                                           self.show_outliers, self.cmb_palette])
            # Create the underlying widget container and form layout
            self.UI_Setup_AxesParms_Subcontainer()
            # Add widgets to form layout and connect to appropriate function/signal
            for description, widget in zip(["X Axis Variable", "Y Axis variable", "Hue Grouping Variable",
                                            "Relative bar width", "Shift plots when hue nesting?",
                                            "Linewidth of the boxes", "Expected proportion of outliers",
                                            "Show outliers?", "Palette"],
                                           [self.le_x, self.le_y, self.le_hue, self.barwidth, self.dodge,
                                            self.linewidth, self.outlier_prop, self.show_outliers, self.cmb_palette]):
                self.formlay_axesparms.addRow(description, widget)
                connect_widget_to_signal(widget, self.sendSignal_plotupdate_axescall)

        elif plot_type == "Scatter Plot":
            self.plotting_func = sns.scatterplot
            self.size_grouper = DandD_Graphing_ListWidget2LineEdit(self.parent_cw, self.all_dtypes)
            self.size_grouper.setPlaceholderText("(Optional) Drag & Drop the Markersize-Grouping Variable")
            self.style_grouper = DandD_Graphing_ListWidget2LineEdit(self.parent_cw, ["object", "category"])
            self.style_grouper.setPlaceholderText("(Optional) Drag & Drop the Markerstyle-Grouping Variable")
            self.cmb_palette = self.UI_Setup_PaletteCombobox()
            self.spinbox_markersize = QDoubleSpinBox(maximum=100, minimum=0, value=40, singleStep=1)
            self.axes_kwargs = self.UI_Setup_AxesMappings(["size", "style", "palette", "s"],
                                                          [self.size_grouper, self.style_grouper,
                                                           self.cmb_palette, self.spinbox_markersize])
            # Create the underlying widget container and form layout
            self.UI_Setup_AxesParms_Subcontainer()
            # Add widgets to form layout and connect to appropriate function/signal
            for description, widget in zip(["X Axis Variable", "Y Axis variable", "Hue Grouping Variable",
                                            "Size Grouping Variable", "Marker Style Grouping Variable", "Marker size",
                                            "Palette"],
                                           [self.le_x, self.le_y, self.le_hue, self.size_grouper, self.style_grouper,
                                            self.spinbox_markersize, self.cmb_palette]):
                self.formlay_axesparms.addRow(description, widget)
                connect_widget_to_signal(widget, self.sendSignal_plotupdate_axescall)

        # Setup the wigets associated with connecting to the legend parameters each time
        if plot_type != "Select a plot type":
            self.UI_Setup_Btn2LegendParms()

        # Regardless of what type of plot was selected, the current axes should be cleared
        self.artist.clear_axes()
        self.artist.canvas.draw()

    ##########################################################
    # Functions for SetupUI convenience and codeline reduction
    ##########################################################
    def UI_Setup_AxesParms_Subcontainer(self):
        self.subcont_axesparms = QWidget()
        self.axes_widget = True
        self.vlay_axesparms.addWidget(self.subcont_axesparms)
        self.formlay_axesparms = QFormLayout(self.subcont_axesparms)

    @staticmethod
    def UI_Setup_AxesMappings(keywords: list, widgets: list):
        mapping = {}
        for keyword, widget in zip(keywords, widgets):
            if isinstance(widget, (QDoubleSpinBox, QSpinBox)):
                mapping[keyword] = widget.value
            elif isinstance(widget, QComboBox):
                mapping[keyword] = widget.currentText
            elif isinstance(widget, QCheckBox):
                mapping[keyword] = widget.isChecked
        return mapping

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

    #############################################################################
    # Functions Designed to send signals out and sometimes set up parameter dicts
    #############################################################################

    def sendSignal_plotupdate_paddingcall(self):
        self.padding_kwargs = {"left": self.spinbox_leftpad.value(),
                               "bottom": self.spinbox_bottompad.value(),
                               "right": self.spinbox_rightpad.value(),
                               "top": self.spinbox_toppad.value(),
                               "wspace": self.spinbox_wspace.value(),
                               "hspace": self.spinbox_hspace.value()}
        self.signal_paddingparms_updateplot.emit()

    def sendSignal_plotupdate_axescall(self):
        self.signal_axesparms_updateplot.emit()

    def sendSignal_plotupdate_figurecall(self):
        self.signal_figparms_updateplot.emit()

    def on_subset(self):
        self.sendSignal_plotupdate_figurecall()


class xASL_GUI_FacetLegend(QWidget):
    signal_legendcall_updateplot = Signal(dict)

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

        self.legend_kwargs = {
            "loc": self.cmb_location.currentText,
            "fontsize": self.cmb_fontsize.currentText,
            "frameon": self.chk_frameon.isChecked,
            "fancybox": self.chk_fancybox.isChecked,
            "shadow": self.chk_shadow.isChecked,
            "markerscale": self.spinbox_markersize.value,
            "markerfirst": self.chk_markerfirst.isChecked
        }

    def sendSignal_legendparms_updateplot(self):
        constructor = {key: call() for key, call in self.legend_kwargs.items()}
        if self.chk_manual_loc.isChecked():
            constructor["bbox_to_anchor"] = (self.spinbox_legend_x.value(), self.spinbox_legend_y.value())
        self.signal_legendcall_updateplot.emit(constructor)
