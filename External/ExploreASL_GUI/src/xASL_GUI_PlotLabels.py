from PySide2.QtWidgets import *
from PySide2.QtGui import *
from PySide2.QtCore import *
from src.xASL_GUI_HelperFuncs_WidgetFuncs import connect_widget_to_signal, disconnect_widget_and_reset
from src.xASL_GUI_HelperFuncs_WidgetFuncs import set_formlay_options
from platform import system


class xASL_GUI_PlotLabels(QWidget):
    """
    This class is responsible for preparing the widget that will be associated with altering plot settings such as
    x-axis label, y-axis label, title, etc.
    """

    signal_xaxislabel_changed = Signal(dict)
    signal_yaxislabel_changed = Signal(dict)
    signal_title_changed = Signal(dict)

    def __init__(self, parent=None):
        super().__init__(parent)
        self.parent = parent  # The Facet Grid Manager
        self.setWindowFlag(Qt.Window)
        self.setWindowTitle("Axis Labels Settings")
        self.setMinimumSize(300, 480)
        self.mainlay = QVBoxLayout(self)
        self.grp_xaxislabel = QGroupBox(title="X-Axis Label Settings")
        self.grp_yaxislabel = QGroupBox(title="Y-Axis Label Settings")
        self.grp_title = QGroupBox(title="Title Settings")
        self.formlay_xaxislabel = QFormLayout(self.grp_xaxislabel)
        self.formlay_yaxislabel = QFormLayout(self.grp_yaxislabel)
        self.formlay_title = QFormLayout(self.grp_title)

        self.Setup_UI_XAxisLabelParms()
        self.Setup_UI_YAxisLabelParms()
        self.Setup_UI_TitleParms()

        self.mainlay.addWidget(self.grp_title)
        self.mainlay.addWidget(self.grp_xaxislabel)
        self.mainlay.addWidget(self.grp_yaxislabel)

        # Additional MacOS actions
        if system() == "Darwin":
            set_formlay_options(self.formlay_title)
            set_formlay_options(self.formlay_xaxislabel)
            set_formlay_options(self.formlay_yaxislabel)

    def Setup_UI_XAxisLabelParms(self):
        # X-Axis Label Properties
        self.le_xaxislabeltext = QLineEdit()
        self.cmb_xaxislabelsize = QComboBox()
        self.cmb_xaxislabelsize.addItems(['xx-small', 'x-small', 'small', 'medium', 'large', 'x-large', 'xx-large'])
        self.cmb_xaxislabelsize.setCurrentIndex(3)
        self.spinbox_xaxislabelrot = QDoubleSpinBox(maximum=90, minimum=0, value=0, singleStep=1)
        self.cmb_xaxislabelweight = QComboBox()
        self.cmb_xaxislabelweight.addItems(["normal", "bold", "black"])
        self.cmb_xaxislabelweight.setCurrentIndex(0)
        self.xaxis_widgets = [self.le_xaxislabeltext, self.cmb_xaxislabelsize, self.spinbox_xaxislabelrot,
                              self.cmb_xaxislabelweight]
        self.xaxis_defaults = ["", 3, 0, 0]
        # Connect the widgets to their functions and add to the appropriate form layout
        for description, widget in zip(["Text", "Fontsize", "Rotation", "Fontweight"], self.xaxis_widgets):
            connect_widget_to_signal(widget, self.sendSignal_xaxislabel_updateplot)
            self.formlay_xaxislabel.addRow(description, widget)

    def Setup_UI_YAxisLabelParms(self):
        # Y-Axis Label Properties
        self.le_yaxislabeltext = QLineEdit()
        self.cmb_yaxislabelsize = QComboBox()
        self.cmb_yaxislabelsize.addItems(['xx-small', 'x-small', 'small', 'medium', 'large', 'x-large', 'xx-large'])
        self.cmb_yaxislabelsize.setCurrentIndex(3)
        self.spinbox_yaxislabelrot = QDoubleSpinBox(maximum=90, minimum=0, value=90, singleStep=1)
        self.cmb_yaxislabelweight = QComboBox()
        self.cmb_yaxislabelweight.addItems(["normal", "bold", "black"])
        self.cmb_yaxislabelweight.setCurrentIndex(0)
        self.yaxis_widgets = [self.le_yaxislabeltext, self.cmb_yaxislabelsize, self.spinbox_yaxislabelrot,
                              self.cmb_yaxislabelweight]
        self.yaxis_defaults = ["", 3, 90, 0]
        # Connect the widgets to their functions and add to the appropriate form layout
        for description, widget in zip(["Text", "Fontsize", "Rotation", "Fontweight"], self.yaxis_widgets):
            connect_widget_to_signal(widget, self.sendSignal_yaxislabel_updateplot)
            self.formlay_yaxislabel.addRow(description, widget)

    def Setup_UI_TitleParms(self):
        # Title Properties
        self.le_titletext = QLineEdit()
        self.cmb_titleloc = QComboBox()
        self.cmb_titleloc.addItems(['center', 'left', 'right'])
        self.cmb_titlesize = QComboBox()
        self.cmb_titlesize.addItems(['xx-small', 'x-small', 'small', 'medium', 'large', 'x-large', 'xx-large'])
        self.cmb_titlesize.setCurrentIndex(3)
        self.cmb_titleweight = QComboBox()
        self.cmb_titleweight.addItems(["normal", "bold", "black"])
        self.cmb_titleweight.setCurrentIndex(0)
        self.title_widgets = [self.le_titletext, self.cmb_titleloc, self.cmb_titlesize, self.cmb_titleweight]
        self.title_defaults = ["", 0, 3, 0]
        # Connect the widgets to their functions and add to the appropriate form layout
        for description, widget in zip(["Text", "Location", "Fontsize", "Fontweight"], self.title_widgets):
            connect_widget_to_signal(widget, self.sendSignal_title_updateplot)
            self.formlay_title.addRow(description, widget)

    def reset_labels(self):
        """
        This function resets every widget back to their defaults and then reconnects the widget to that signal.
        This avoids overhead of constantly updating a canvas as widgets are reset back to their defaults in the
        event of a new figure being implemented
        """
        for widget, signal, default in zip(self.xaxis_widgets + self.yaxis_widgets + self.title_widgets,
                                           [self.sendSignal_xaxislabel_updateplot] * 4 +
                                           [self.sendSignal_yaxislabel_updateplot] * 4 +
                                           [self.sendSignal_title_updateplot] * 4,
                                           self.xaxis_defaults + self.yaxis_defaults + self.title_defaults):
            disconnect_widget_and_reset(widget=widget, target_signal=signal, default=default)
            connect_widget_to_signal(widget=widget, target_signal=signal)

    def sendSignal_xaxislabel_updateplot(self):
        """
        This function prepares the arguments necessary for updating the x-axis label and emits them
        """
        xaxis_kwargs = {"xlabel": self.le_xaxislabeltext.text(),
                        "fontsize": self.cmb_xaxislabelsize.currentText(),
                        "fontweight": self.cmb_xaxislabelweight.currentText(),
                        "rotation": self.spinbox_xaxislabelrot.value()}
        self.signal_xaxislabel_changed.emit(xaxis_kwargs)

    def sendSignal_yaxislabel_updateplot(self):
        """
        This function prepares the arguments necessary for updating the y-axis label
        """
        xaxis_kwargs = {"ylabel": self.le_yaxislabeltext.text(),
                        "fontsize": self.cmb_yaxislabelsize.currentText(),
                        "fontweight": self.cmb_yaxislabelweight.currentText(),
                        "rotation": self.spinbox_yaxislabelrot.value()}
        self.signal_yaxislabel_changed.emit(xaxis_kwargs)

    def sendSignal_title_updateplot(self):
        """
        This function prepares the arguments necessary for updating the title
        """
        title_kwargs = {"label": self.le_titletext.text(),
                        "fontdict": {
                            "fontsize": self.cmb_titlesize.currentText(),
                            "fontweight": self.cmb_titleweight.currentText()
                        },
                        "loc": self.cmb_titleloc.currentText()
                        }
        self.signal_title_changed.emit(title_kwargs)
        print("sendSignal_title_update")
