import os
from PySide2.QtGui import QIcon
from PySide2.QtCore import QSize
from PySide2.QtWidgets import QDoubleSpinBox, QSpinBox, QComboBox, QLineEdit, QCheckBox
from ExploreASL_GUI.xASL_GUI_HelperClasses import DandD_Graphing_ListWidget2LineEdit


def set_widget_icon(widget, config: dict, icon_name: str, size: tuple = None):
    """
    Convenience function for setting a widget to contain an icon of a particular size
    :param widget: the widget for which an icon should be set
    :param config: the config instance, so that the appropriate filepath stored may be accessed
    :param icon_name: the basename of the icon
    :param size: tuple of width by height, in pixels
    """
    icon_path = os.path.join(config["ProjectDir"], "media", icon_name)
    widget.setIcon(QIcon(icon_path))
    if size is not None:
        widget.setIconSize(QSize(*size))


def connect_widget_to_signal(widget, target_signal):
    """
    Convenience function for connecting a widget to a signal; useful in for loops
    :param widget: the widget to connect
    :param target_signal: the signal to connect to
    """
    if isinstance(widget, QComboBox):
        widget.currentTextChanged.connect(target_signal)
    elif isinstance(widget, (QSpinBox, QDoubleSpinBox)):
        widget.valueChanged.connect(target_signal)
    elif isinstance(widget, (DandD_Graphing_ListWidget2LineEdit, QLineEdit)):
        widget.textChanged.connect(target_signal)
    elif isinstance(widget, QCheckBox):
        widget.clicked.connect(target_signal)
    else:
        print(f'{widget} could not be connected')


def disconnect_widget_and_reset(widget, target_signal, default):
    """
    Convenience function for disconnecting a widget from a signal and resetting the widget back to a default value
    without triggering those previous connections
    :param widget: the widget to disconnect
    :param target_signal: the signal to disconnect from
    :param default: the default value to change to after wards
    """
    if isinstance(widget, QComboBox):
        widget.currentTextChanged.disconnect(target_signal)
        widget.setCurrentIndex(default)
    elif isinstance(widget, (QSpinBox, QDoubleSpinBox)):
        widget.valueChanged.disconnect(target_signal)
        widget.setValue(default)
    elif isinstance(widget, (DandD_Graphing_ListWidget2LineEdit, QLineEdit)):
        widget.textChanged.disconnect(target_signal)
        widget.setText(default)
    elif isinstance(widget, QCheckBox):
        widget.clicked.disconnect(target_signal)
        widget.setChecked(default)
    else:
        print(f'{widget} could not be connected')
