from platform import system
from PySide2.QtWidgets import QLineEdit


def set_os_dependent_text(linedit: QLineEdit, config_ossystem: str = '', text_to_set: str = ''):
    if config_ossystem == '' or config_ossystem:
        config_ossystem = system()

    if config_ossystem == "Windows":
        linedit.setText(text_to_set.replace('/', '\\'))
    else:
        linedit.setText(text_to_set.replace('\\', '/'))
