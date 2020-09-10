from PySide2.QtWidgets import QWidget, QLabel, QComboBox, QFormLayout, QVBoxLayout
from PySide2.QtGui import Qt, QFont


# noinspection PyAttributeOutsideInit
class xASL_GUI_Subsetter(QWidget):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.setWindowFlag(Qt.Window)
        self.setWindowTitle("Subset the data")
        self.font_format = QFont()
        self.font_format.setPointSize(16)
        self.parent_cw = parent

        # This will be used to keep track of columns that have already been added to the subset form layout
        self.current_rows = {}

        # Main Setup
        self.Setup_UI_MainWidgets()

    def Setup_UI_MainWidgets(self):
        self.mainlay = QVBoxLayout(self)
        self.formlay_headers = QFormLayout()
        self.formlay_subsets = QFormLayout()
        self.lab_inform_reload = QLabel("Changes will not take place until you re-load the data")
        self.lab_inform_reload.setFont(self.font_format)
        self.formlay_headers.addRow(QLabel("Column Names"), QLabel("\t\tSubset on"))

        self.mainlay.addLayout(self.formlay_headers)
        self.mainlay.addLayout(self.formlay_subsets)
        self.mainlay.addWidget(self.lab_inform_reload)

    # Updates the current fields that are permitted to be subset.
    def update_subsetable_fields(self):
        df = self.parent_cw.loader.loaded_long_data
        colnames = df.columns
        for name in colnames:
            # Skip any column based on criteria
            if any([name in self.current_rows.keys(),  # skip if already encountered
                    str(df[name].dtype) not in ["object", "category"],  # skip if not categorical
                    name in ["SUBJECT", "SubjectNList"]  # skip some unnecessary defaults with too many categories
                    ]):
                continue

            # Otherwise, create a combobox, add it to the form layout, and add the method to the current text as a
            # value, as this will be used to
            cmb = QComboBox()
            cmb.addItems(["Select a subset"] + df[name].unique().tolist())
            self.formlay_subsets.addRow(name, cmb)
            self.current_rows[name] = cmb.currentText

    # Convenience function for purging everything in the event that either the analysis directory changes or the
    # supplementary data file provided is changed
    def clear_contents(self):
        # Clear the dict
        self.current_rows.clear()
        # Clear the rows of the format layout
        for ii in range(self.formlay_subsets.rowCount()):
            self.formlay_subsets.removeRow(0)

    # The active function that will subset the data; takes in a dataframe and subsets it
    def subset_data(self, long_df):
        for key, call in self.current_rows.items():
            if call() != "Select a subset":
                print(f"Subsetting {key} on {call()}")
                long_df = long_df.loc[long_df[key] == call(), :]

        return long_df
