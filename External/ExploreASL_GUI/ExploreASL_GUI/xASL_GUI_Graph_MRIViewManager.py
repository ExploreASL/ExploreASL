from PySide2.QtWidgets import *
from PySide2.QtGui import *
from PySide2.QtCore import *
import seaborn as sns
from ExploreASL_GUI.xASL_GUI_Graph_MRIViewArtist import xASL_GUI_MRIViewArtist
from ExploreASL_GUI.xASL_GUI_HelperClasses import DandD_Graphing_ListWidget2LineEdit
from ExploreASL_GUI.xASL_GUI_HelperFuncs import connect_widget_to_signal
import json
import os


class xASL_GUI_MRIViewManager(QWidget):

    signal_axesparms_updateplot = Signal()

    def __init__(self, parent):
        super().__init__(parent=parent)
        self.parent_cw = parent
        self.artist: xASL_GUI_MRIViewArtist = self.parent_cw.fig_artist
        self.mainlay = QVBoxLayout(self)
        self.spacer = QSpacerItem(0, 1, QSizePolicy.Preferred, QSizePolicy.Expanding)

        with open(os.path.join(self.parent_cw.config["ProjectDir"],
                               "JSON_LOGIC",
                               "GraphingParameters.json")) as freader:
            parms = json.load(freader)
            self.all_dtypes = parms["all_dtypes"]
            self.numeric_dtypes = parms["numeric_dtypes"]
            self.categorical_dtypes = ["object", "category"]
            self.palettenames = parms["palettenames"]
            self.default_palette_idx = self.palettenames.index("Set1")

        # Key Variables
        self.analysis_dir = self.parent_cw.le_analysis_dir.text()
        try:
            with open(os.path.join(self.analysis_dir, "DataPar.json")) as parms_reader:
                parms = json.load(parms_reader)
                self.subject_regex = parms["subject_regexp"].strip("$^")

            self.axes_arg_x = ''  # Set to '' instead of None because the latter is not supported by lineedits
            self.axes_arg_y = ''
            self.axes_arg_hue = ''
            self.axes_widget = None

            self.subject_sessions = self.parent_cw.loader.loaded_wide_data["SUBJECT"].tolist()

            self.UI_Setup_Tabs()
            self.UI_Setup_MRIViewParameters()
            self.error_init = False

        except FileNotFoundError:
            QMessageBox.warning(self.parent_cw,
                                "Unable to load Plot & MRI Viewer",
                                "No DataPar.json was detected in the analysis directory, which is required by the MRI"
                                " Viewer to be able to select subject images",
                                QMessageBox.Ok)
            self.error_init = True
        except KeyError:
            QMessageBox.warning(self.parent_cw,
                                "Unable to load Plot & MRI Viewer",
                                "The DataPar.json file either did not contain the appropriate regex key or the regex"
                                " value was corrupt",
                                QMessageBox.Ok)
            self.error_init = True

    def UI_Setup_Tabs(self):
        self.tab_pltsettings = QTabWidget()
        self.cont_scatterparms, self.cont_mriparms = QWidget(), QWidget()
        self.tab_pltsettings.addTab(self.cont_scatterparms, "Scatterplot Parameters")
        self.tab_pltsettings.addTab(self.cont_mriparms, "MRI Viewer Parameters")

        self.vlay_axesparms = QVBoxLayout(self.cont_scatterparms)
        self.formlay_mriparms = QFormLayout(self.cont_mriparms)
        self.mainlay.addWidget(self.tab_pltsettings)

        # Might move this elsewhere; but for now define the axes type selection here
        self.cmb_axestype = QComboBox()
        self.cmb_axestype.addItems(["Select a plot type", "Strip Plot", "Scatter Plot"])
        self.cmb_axestype.currentTextChanged.connect(self.UI_Setup_AxesParms)
        self.vlay_axesparms.addWidget(self.cmb_axestype)
        self.vlay_axesparms.addSpacerItem(self.spacer)

    def UI_Setup_ConnectManager2Artist(self):
        self.artist: xASL_GUI_MRIViewArtist = self.parent_cw.fig_artist
        if self.artist is not None:
            self.cmb_selectsubject.currentTextChanged.connect(self.artist.switch_subject)
            self.slider_axialslice.valueChanged.connect(self.artist.plotupdate_axialslice)
            self.slider_coronalslice.valueChanged.connect(self.artist.plotupdate_coronalslice)
            self.slider_sagittalslice.valueChanged.connect(self.artist.plotupdate_sagittalslice)
            self.signal_axesparms_updateplot.connect(self.artist.plotupdate_axescall)
        else:
            print("UI_Setup_Connections; the artist had a value of None")

    def UI_Setup_AxesParms(self, plot_type):
        self.clear_axesparms()
        print(f"Selected {plot_type} as the Axes Type")

        if plot_type == "Strip Plot":
            self.plotting_func = sns.stripplot
            self.le_x = DandD_Graphing_ListWidget2LineEdit(self.parent_cw, self.categorical_dtypes)
            self.le_x.setPlaceholderText("Drag & Drop the X-axis Variable")
            self.le_y = DandD_Graphing_ListWidget2LineEdit(self.parent_cw, self.numeric_dtypes)
            self.le_y.setPlaceholderText("Drag & Drop the Y-axis Variable")
            self.le_hue = DandD_Graphing_ListWidget2LineEdit(self.parent_cw, ["object", "category"])
            self.le_hue.setPlaceholderText("Drag & Drop the Hue Grouping Variable")
            self.axes_arg_x = self.le_x.text
            self.axes_arg_y = self.le_y.text
            self.axes_arg_hue = self.le_hue.text

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

            self.vlay_axesparms.removeItem(self.spacer)

        elif plot_type == "Scatter Plot":
            self.plotting_func = sns.scatterplot
            self.le_x = DandD_Graphing_ListWidget2LineEdit(self.parent_cw, self.numeric_dtypes)
            self.le_x.setPlaceholderText("Drag & Drop the X-axis Variable")
            self.le_y = DandD_Graphing_ListWidget2LineEdit(self.parent_cw, self.numeric_dtypes)
            self.le_y.setPlaceholderText("Drag & Drop the Y-axis Variable")
            self.le_hue = DandD_Graphing_ListWidget2LineEdit(self.parent_cw, ["object", "category"])
            self.le_hue.setPlaceholderText("Drag & Drop the Hue Grouping Variable")
            self.axes_arg_x = self.le_x.text
            self.axes_arg_y = self.le_y.text
            self.axes_arg_hue = self.le_hue.text

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

            self.vlay_axesparms.removeItem(self.spacer)

        else:
            self.axes_arg_x = ''
            self.axes_arg_y = ''
            self.axes_arg_hue = ''

        # Regardless of what type of plot was selected, the current axes should be cleared
        self.artist.clear_axes()
        self.artist.plotting_canvas.draw()
        self.cmb_selectsubject.setCurrentIndex(0)

    def UI_Setup_MRIViewParameters(self):
        self.cmb_selectsubject = QComboBox()
        self.cmb_selectsubject.addItems(["Select a subject and run"] + self.subject_sessions)

        self.slider_axialslice = QSlider(Qt.Horizontal, minimum=0, maximum=120, singleStep=1, value=50)
        self.slider_coronalslice = QSlider(Qt.Horizontal, minimum=0, maximum=144, singleStep=1, value=50)
        self.slider_sagittalslice = QSlider(Qt.Horizontal, minimum=0, maximum=120, singleStep=1, value=50)

        self.formlay_mriparms.addRow("Subject_Run Selection", self.cmb_selectsubject)
        self.formlay_mriparms.addRow("Axial Slice", self.slider_axialslice)
        self.formlay_mriparms.addRow("Coronal Slice", self.slider_coronalslice)
        self.formlay_mriparms.addRow("Sagittal Slice", self.slider_sagittalslice)

    # Convenience Function - clears the Axes Parameters Tab
    # noinspection PyTypeChecker
    def clear_axesparms(self):
        print("Clearing Axes Parms")
        if self.axes_widget is not None:
            self.vlay_axesparms.removeWidget(self.subcont_axesparms)
            self.vlay_axesparms.addSpacerItem(self.spacer)
            self.subcont_axesparms.setParent(None)
            del self.subcont_axesparms
            self.axes_widget = None

    #############################################################################
    # Functions Designed to send signals out and sometimes set up parameter dicts
    #############################################################################
    def sendSignal_plotupdate_axescall(self):
        """
        On a change in any of the Axes Parameters, a signal is sent out to update the artist
        """
        self.signal_axesparms_updateplot.emit()

    def on_subset(self):
        """
        On a change in the data being subset, a signal is sent out to update the artist
        """
        self.signal_axesparms_updateplot.emit()

    ##########################################################
    # Functions for SetupUI convenience and codeline reduction
    ##########################################################
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
