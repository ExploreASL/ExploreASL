from PySide2.QtWidgets import QLineEdit, QAbstractItemView, QListWidget
from PySide2.QtCore import Qt, Signal, QModelIndex
import os
from ExploreASL_GUI.xASL_GUI_HelperFuncs_StringOps import set_os_dependent_text


class DandD_Graphing_ListWidget2LineEdit(QLineEdit):
    """
    Modified QLineEdit to support accepting text drops from a QListWidget or QAbstractItemModel derivatives
    """

    def __init__(self, postproc_widget, dtype_list, parent=None):
        # print(dtype_list)
        super().__init__(parent)
        self.permitted_dtypes = dtype_list
        self.setAcceptDrops(True)
        self.PostProc_widget = postproc_widget

    def dragEnterEvent(self, event) -> None:
        # print(event.mimeData().formats())
        if event.mimeData().hasFormat("application/x-qabstractitemmodeldatalist"):
            event.accept()
        else:
            event.ignore()
            print("dragEnterEvent ignored")

    def dragMoveEvent(self, event) -> None:
        # print(event.mimeData().formats())
        if event.mimeData().hasFormat("application/x-qabstractitemmodeldatalist"):
            event.setDropAction(Qt.CopyAction)
            event.accept()
        else:
            event.ignore()

    def dropEvent(self, event) -> None:
        if all([event.mimeData().hasFormat("application/x-qabstractitemmodeldatalist"),
                isinstance(event.source(), QAbstractItemView)]):
            # Get the text that has been dragged
            ix: QModelIndex = event.source().currentIndex()
            colname = ix.data()

            # Get the current data types of the loaded data
            data_dtypes = {col: str(name) for col, name in
                           self.PostProc_widget.loader.long_data.dtypes.to_dict().items()}
            # print(f"The dtypes of the current dataframe are: {data_dtypes}")
            # print(f"For column: {colname}, the detected datatype was: {data_dtypes[colname]}")
            # print(f"The permitted dtypes for this widget are: {self.permitted_dtypes}")
            if data_dtypes[colname] in self.permitted_dtypes:
                event.accept()
                self.setText(colname)
            else:
                event.ignore()
        else:
            event.ignore()


class DandD_FileExplorer2LineEdit(QLineEdit):
    """
    Modified QLineEdit to support accepting text drops from a file explorer

    Signature: DandD_FileExplorer2LineEdit(parent, acceptable_path_type, supported_extensions)
        - parent = the parent widget
        - acceptable_path_type = the type of filepath it will accept.
            - "File",
            - "Directory"
            - "Both"
        - supported_extensions = a list of the .filetypes that this lineedit will accept (i.e .csv , .txt , etc.),

    """

    def __init__(self, parent=None, acceptable_path_type: str = "Both", supported_extensions: list = None):
        super().__init__(parent)
        self.setAcceptDrops(True)
        self.path_type = acceptable_path_type
        if supported_extensions is None or supported_extensions == "All":
            self.supported_ext = []
        else:
            self.supported_ext = supported_extensions

        # Quality Control
        if acceptable_path_type not in ["File", "Directory", "Both"]:
            raise ValueError("acceptable_path_type must be one of: 'File', 'Directory', or 'Both'")

    def dragEnterEvent(self, event) -> None:
        if event.mimeData().hasUrls():
            event.accept()
        else:
            event.ignore()

    def dragMoveEvent(self, event) -> None:
        if event.mimeData().hasUrls():
            event.setDropAction(Qt.CopyAction)
            event.accept()
        else:
            event.ignore()

    def dropEvent(self, event) -> None:
        if event.mimeData().hasUrls():
            event.accept()
            for url in event.mimeData().urls():
                if url.isLocalFile():
                    path_string = str(url.toLocalFile())

                    # Scenario 1: Accept all filetypes
                    if self.path_type == "Both" and os.path.exists(path_string):
                        set_os_dependent_text(linedit=self,
                                              config_ossystem='',
                                              text_to_set=path_string)
                        return

                    # Scenario 2: Accept only Files
                    elif self.path_type == "File" and os.path.isfile(path_string):
                        if len(self.supported_ext) == 0:
                            set_os_dependent_text(linedit=self,
                                                  config_ossystem='',
                                                  text_to_set=path_string)
                            return
                        else:
                            if os.path.splitext(path_string)[1] in self.supported_ext:
                                set_os_dependent_text(linedit=self,
                                                      config_ossystem='',
                                                      text_to_set=path_string)
                                return
                            else:
                                event.ignore()

                    # Scenario 3: Accept only Direcories
                    elif self.path_type == "Directory" and os.path.isdir(path_string):
                        set_os_dependent_text(linedit=self,
                                              config_ossystem='',
                                              text_to_set=path_string)
                        return

                    else:
                        event.ignore()
        else:
            event.ignore()


class DandD_FileExplorer2ListWidget(QListWidget):
    """
    Class meant to accept MULTIPLE directory inputs and add them to the underlying QListWidget
    """
    alert_regex = Signal()

    def __init__(self, parent=None):
        super().__init__(parent)
        self.setAcceptDrops(True)

    def dragEnterEvent(self, event) -> None:
        if event.mimeData().hasUrls():
            event.accept()
        else:
            event.ignore()

    def dragMoveEvent(self, event) -> None:
        if event.mimeData().hasUrls():
            event.setDropAction(Qt.CopyAction)
            event.accept()
        else:
            event.ignore()

    def dropEvent(self, event) -> None:
        if event.mimeData().hasUrls():
            event.accept()
            subject_directories = []
            for url in event.mimeData().urls():
                if url.isLocalFile():
                    subject_directories.append(str(url.toLocalFile()))
            # Only filter for the basenames of directories; also, avoid bringing in unnecessary directories like lock
            basenames = [os.path.basename(directory) for directory in subject_directories if os.path.isdir(directory)
                         and directory not in ['lock', "Population"]]
            # Avoid double-dipping the names
            current_names = [self.item(idx).text() for idx in range(self.count())]
            filtered_basenames = [name for name in basenames if name not in current_names]
            self.addItems(filtered_basenames)
            self.alert_regex.emit()
        else:
            event.ignore()
