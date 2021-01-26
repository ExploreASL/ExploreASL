import os
from PySide2.QtGui import QIcon, Qt
from PySide2.QtCore import QSize
from PySide2.QtWidgets import (QDoubleSpinBox, QSpinBox, QComboBox, QLineEdit, QCheckBox, QHBoxLayout, QPushButton,
                               QVBoxLayout, QScrollArea, QWidget, QFormLayout)
from src.xASL_GUI_HelperClasses import (DandD_Graphing_ListWidget2LineEdit, DandD_FileExplorer2LineEdit)
from typing import Tuple, Union, Any


def set_formlay_options(formlay: QFormLayout,
                        field_growth: Union[str, QFormLayout.FieldGrowthPolicy] = QFormLayout.ExpandingFieldsGrow,
                        formside_alignment: Union[Tuple[str, str],
                                                  Tuple[Qt.Alignment, Qt.Alignment]] = (Qt.AlignLeft, Qt.AlignTop),
                        labelside_alignment: Union[str, Qt.Alignment] = Qt.AlignLeft,
                        row_wrap_policy: Union[str, QFormLayout.RowWrapPolicy] = QFormLayout.WrapLongRows,
                        vertical_spacing: int = None,
                        horizontal_spacing: int = None):
    """
    Convenience function for setting the options for a QFormLayout

    :param formlay: the QFormLayout widget to be altered
    :param field_growth: string or FieldGrowthPolicy. For string options, acceptable strings are:
    "at_size_hint", "expanding_fields_grow", or "all_nonfixed_grow"
    :param formside_alignment: string or Qt.Alignment. For string options, acceptable strings are:
    "left", "right", "top", or "bottom".
    :param labelside_alignment: string or Qt.Alignment. For string options, acceptable strings are:
    "left", "right", "top" or "bottom".
    :param row_wrap_policy: string or QFormLayout.RowWrapPolicy. For string options, acceptable strings are:
    "dont_wrap" (Fields are always beside Labels); "wrap_long" (Labels column has enough spacing to accomodate
    the widest label); "wrap_all" (Labels are above their Fields)
    :param vertical_spacing: the vertical spacing between rows, in pixels
    :param horizontal_spacing: the horizontal spacing within a row, in pixels
    """
    field_dict = {"at_size_hint": QFormLayout.FieldsStayAtSizeHint,
                  "expanding_fields_grow": QFormLayout.ExpandingFieldsGrow,
                  "all_nonfixed_grow": QFormLayout.AllNonFixedFieldsGrow}
    align_dict = {"left": Qt.AlignLeft, "right": Qt.AlignRight, "top": Qt.AlignTop, "bottom": Qt.AlignBottom}
    wrap_dict = {"dont_wrap": QFormLayout.DontWrapRows,
                 "wrap_long": QFormLayout.WrapLongRows,
                 "wrap_all": QFormLayout.WrapAllRows}

    if isinstance(field_growth, str):
        formlay.setFieldGrowthPolicy(field_dict[field_growth])
    else:
        formlay.setFieldGrowthPolicy(field_growth)

    if isinstance(formside_alignment[0], str):
        formlay.setFormAlignment(align_dict[formside_alignment[0]] | align_dict[formside_alignment[1]])
    else:
        formlay.setFormAlignment(formside_alignment[0] | formside_alignment[1])

    if isinstance(labelside_alignment, str):
        formlay.setLabelAlignment(align_dict[labelside_alignment])
    else:
        formlay.setLabelAlignment(labelside_alignment)

    if isinstance(row_wrap_policy, str):
        formlay.setRowWrapPolicy(wrap_dict[row_wrap_policy])
    else:
        formlay.setRowWrapPolicy(row_wrap_policy)

    if vertical_spacing is not None:
        formlay.setVerticalSpacing(vertical_spacing)
    if horizontal_spacing is not None:
        formlay.setHorizontalSpacing(horizontal_spacing)

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


def make_droppable_clearable_le(le_connect_to=None, btn_connect_to=None, default='', **kwargs):
    """
    Function for creating a typical QLineEdit-QPushButton pair encapsulated within a QHboxLayout.

    :param le_connect_to: the function that the lineedit's textChanged signal should connect to, if any
    :param btn_connect_to: the function that the pushbutton's clicked signal should connect to, if any
    :param default: the default text that should be present in the lineedit
    :param kwargs: additional keywords that are fed into DandD_FileExplorer2LineEdit
    """
    hlay = QHBoxLayout()
    le = DandD_FileExplorer2LineEdit(**kwargs)
    le.setText(default)
    le.setClearButtonEnabled(True)
    if le_connect_to is not None:
        le.textChanged.connect(le_connect_to)
    btn = QPushButton("...", )
    if btn_connect_to is not None:
        btn.clicked.connect(btn_connect_to)
    hlay.addWidget(le)
    hlay.addWidget(btn)
    return hlay, le, btn


def make_scrollbar_area(parent, orientation: str = "v",
                        margins: Tuple[int, int, int, int] = (0, 0, 0, 0)) -> Tuple[QVBoxLayout, QScrollArea, QWidget]:
    """
    Function for creating a QVBoxLayout within which is placed a QScrollArea set to a QWidget container

    :param parent: The parent widget of the layout to be returned
    :param orientation: The type of layout to return - "v" or "h" for QVBoxLayout or QHBoxLayout, respectively. Default
    is vertical.
    :param margins: The margins to set around the layout. Tuple of (left, top, right, bottom) as integers
    """
    o_dict = {"v": QVBoxLayout, "h": QHBoxLayout}
    vlay, scrollarea, container = o_dict[orientation](parent), QScrollArea(), QWidget()
    vlay.setContentsMargins(*margins)
    scrollarea.setWidget(container)
    scrollarea.setWidgetResizable(True)
    vlay.addWidget(scrollarea)
    return vlay, scrollarea, container
