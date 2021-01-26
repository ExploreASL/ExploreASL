from PySide2.QtWidgets import (QWidget, QLabel, QComboBox, QFormLayout, QVBoxLayout, QPushButton,
                               QScrollArea, QSizePolicy)
from PySide2.QtGui import Qt, QFont
from PySide2.QtCore import Signal, Slot
import pandas as pd
import numpy as np
from src.xASL_GUI_HelperFuncs_WidgetFuncs import set_formlay_options
from platform import system


# noinspection PyAttributeOutsideInit
class xASL_GUI_Subsetter(QWidget):

    def __init__(self, parent=None):
        super().__init__(parent)
        self.setWindowFlag(Qt.Window)
        self.setWindowTitle("Subset the data")
        self.font_format = QFont()
        self.font_format.setPointSize(16)
        self.parent_cw = parent
        self.setMinimumSize(400, 200)

        # This will be used to keep track of columns that have already been added to the subset form layout
        self.current_subsettable_fields = {}
        self.do_not_add = ['', np.nan, 'nan']
        # Main Setup
        self.Setup_UI_MainWidgets()

        # Additional MacOS actions
        if system() == "Darwin":
            set_formlay_options(self.formlay_subsets)

    def Setup_UI_MainWidgets(self):
        self.mainlay = QVBoxLayout(self)
        self.scroll_subsetter = QScrollArea()
        self.cont_subsetterrows = QWidget()
        self.scroll_subsetter.setWidget(self.cont_subsetterrows)
        self.scroll_subsetter.setWidgetResizable(True)
        self.formlay_subsets = QFormLayout(self.cont_subsetterrows)
        self.formlay_subsets.addRow(QLabel("Column Names"), QLabel("Subset on"))
        self.btn_subset_data = QPushButton("Subset the data", clicked=self.call_subset_data)

        self.mainlay.addWidget(self.scroll_subsetter)
        self.mainlay.addWidget(self.btn_subset_data)

    # Updates the current fields that are permitted to be subset.
    def update_subsetable_fields_on_load(self, df):
        """
        Updates the known fields that may be subsetted when new data is loaded in. This function is necessary, as its
        behavior differs during loading versus intercepting messages from the dtype indicator
        @param df: The dataframe being loaded in
        """
        # Always start off with a clearing of the contents
        self.clear_contents()

        colnames = df.columns
        for colname in colnames:

            # Skip any column based on criteria
            if any([colname in self.current_subsettable_fields.keys(),  # skip if already encountered
                    str(df[colname].dtype) not in ["object", "category"],  # skip if not categorical
                    colname in ["SUBJECT", "SubjectNList", "MeanMotion"]  # skip some unnecessary defaults
                    ]):
                continue

            # Otherwise, create a combobox, add it to the form layout, and add the method to the current text as a
            # value, as this will be used to
            cmb = Subsetter_QCombobox(column_name=colname)
            cmb.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Preferred)
            values = self.parent_cw.loader.loaded_long_data[colname].unique()
            to_add = [value for value in values if value not in self.do_not_add]
            cmb.addItems(["Select a subset"] + to_add)
            self.formlay_subsets.addRow(colname, cmb)
            self.current_subsettable_fields[colname] = cmb

    def add_subsetable_field(self, colname):
        """
        Adds a combobox to the subsetter widget form layout. May be called either during a load or during a dtype
        update.
        @param colname: the name of the column to add
        """
        cmb = Subsetter_QCombobox(column_name=colname)
        values = self.parent_cw.loader.long_data[colname].unique()
        to_add = [value for value in values if value not in self.do_not_add]
        cmb.addItems(["Select a subset"] + to_add)
        self.formlay_subsets.addRow(colname, cmb)
        self.current_subsettable_fields[colname] = cmb

    def remove_subsetable_field(self, colname: str):
        """
        Removes a combobox from the subsetter widget form layout. Called when there is a categorical --> numerical
        change.
        """
        # First remove it from the form layout
        widget = self.current_subsettable_fields[colname]
        row, _ = self.formlay_subsets.getWidgetPosition(widget)
        self.formlay_subsets.removeRow(row)

        # Also remove it from the dict holding the comboboxes
        del self.current_subsettable_fields[colname], widget, row

    @Slot(str, str)
    def add_or_remove_subsetable_field(self, colname: str, new_dtype: str):
        """
        Slot called by a dtype change by the dtype indicator
        """
        print(f"Subsetter add_or_remove_subsetable_field received a signal to alter column {colname} as it became a "
              f"dtype: {new_dtype}")
        # Case 1: dtype switched to is categorical; therefore add it to the form layout
        if new_dtype == "categorical":
            self.add_subsetable_field(colname=colname)

        elif new_dtype == "numerical":
            try:
                is_default = self.current_subsettable_fields[colname].currentText() == "Select a subset"
                # A KeyError means that the user accidentally wishes to revert a decision to convert from numerical to
                # categorical
            except KeyError:
                print(f"Received a KeyError for column: {colname}.\n"
                      f"Check if the following keys match the listed rows in the subsetter:\n"
                      f"{self.current_subsettable_fields.keys()}")
                return

            self.remove_subsetable_field(colname=colname)

            # If the removed field was the default "Select a subset option", don't do anything.
            if is_default:
                return
            # Otherwise, the long_data must be updated (and the figure manager as well)
            else:
                self.subset_data(already_converted=True)

        else:
            raise ValueError(f"add_or_remove_subsetable_field received an inappropriate dtype descriptor.\n"
                             f"The new_dtype indicated was: {new_dtype}")

    def clear_contents(self):
        """
        Convenience function for purging everything in the event that either the analysis directory changes or the
        supplementary data file provided is changed.
        """
        # Clear the dict
        self.current_subsettable_fields.clear()
        # Clear the rows of the format layout
        for ii in range(self.formlay_subsets.rowCount() - 1):
            self.formlay_subsets.removeRow(1)

    def subset_data_on_load(self, long_df: pd.DataFrame):
        """
        Takes in a dataframe & subsets it. This is required due to the fact that the loader's self.long_data_orig
        may not exist at this point
        :param long_df: A dataframe in long format. This gets subset repeatedly according to the comboboxes present
        that do not have 'Select a subset' as their current option
        """
        for key, cmb in self.current_subsettable_fields.items():
            if cmb.currentText() != "Select a subset":
                print(f"Subsetting {key} on {cmb.currentText()}")
                long_df = long_df.loc[long_df[key] == cmb.currentText(), :]

        return long_df

    def call_subset_data(self):
        self.subset_data()

    def subset_data(self, already_converted=False):
        """
        Updates the long_data variable based on the long_data_orig and calls the figure managers "on_subset" callable.
        """
        long_df = self.parent_cw.loader.long_data_orig

        # Part 1: Apply any dtype transformations as appropriate
        if len(self.parent_cw.dtype_indicator.covariate_cols) > 0 and not already_converted:
            print("CONVERTING WITHIN SUBSET_DATA")
            self.parent_cw.loader.long_data = self.parent_cw.dtype_indicator.update_dataframe_dtypes(long_df)

        # Part 2: update the loader's long_data variable by subsetting the long_data_orig
        for key, cmb in self.current_subsettable_fields.items():
            if cmb.currentText() != "Select a subset":
                long_df = long_df.loc[long_df[key] == cmb.currentText(), :]
        self.parent_cw.loader.long_data = long_df

        del long_df

        # Part 3: have the current figure manager update its artist with respect to the changes
        if self.parent_cw.fig_manager is not None:
            self.parent_cw.fig_manager.on_subset()

        print(f"ON SUBSET the shape of the dataframe is {self.parent_cw.loader.long_data.shape}")


class Subsetter_QCombobox(QComboBox):

    def __init__(self, column_name, parent=None):
        super(Subsetter_QCombobox, self).__init__(parent=parent)
        self.associated_column = column_name


class xASL_GUI_Datatype_Indicator(QWidget):
    """
    Class meant to assist users with defining the datatype of their covariates
    """

    signal_inform_subsetter = Signal(str)  # Informs the subsetter of whether the dtype change was cat2num or num2cat

    def __init__(self, parent=None):
        super().__init__(parent)
        self.setWindowFlag(Qt.Window)
        self.setWindowTitle("Clarify Covariate Datatypes")
        self.setMinimumWidth(400)
        self.parent_cw = parent

        self.covariate_cols = {}
        self.off_limits = {"Side of the Brain", "CBF", "Anatomical Area", "MeanMotion", "LongitudinalTimePoint",
                           "AcquisitionTime", "GM_vol", "WM_vol", "CSF_vol", "GM_ICVRatio", "GMWM_ICVRatio",
                           "WMH_vol", "WMH_count", "SUBJECT", "Site", "SubjectNList"}
        self.dtypes_to_general = {"object": "categorical", "category": "categorical", "bool": "categorical",
                                  "float": "numerical", "float16": "numerical", "float32": "numerical",
                                  "float64": "numerical", "int": "numerical", "int8": "numerical", "int16": "numerical",
                                  "int32": "numerical", "int64": "numerical"}
        # Main Setup
        self.Setup_UI_MainWidgets()

    def Setup_UI_MainWidgets(self):
        self.mainlay = QVBoxLayout(self)
        self.scroll_dtype_indicator = QScrollArea()
        self.cont_indicate_dtype_rows = QWidget()
        self.scroll_dtype_indicator.setWidget(self.cont_indicate_dtype_rows)
        self.scroll_dtype_indicator.setWidgetResizable(True)

        self.formlay_dtypes = QFormLayout(self.cont_indicate_dtype_rows)
        self.formlay_dtypes.addRow(QLabel("Column Names"), QLabel("Datatype"))
        self.mainlay.addWidget(self.scroll_dtype_indicator)

    def update_known_covariates(self, df):
        """
        Updates the known covariates and their datatypes during a dataload
        """
        # Always start off with a clearing of the contents
        self.clear_contents()
        # If the lineedit for the covariates dataframe file is empty during a dataload, clear everything
        if self.parent_cw.le_metadata.text() == '':
            return

        # Otherwise, update the widget
        for colname in df.columns.tolist():
            if colname in self.off_limits:
                continue
            else:
                translated_dtype = self.dtypes_to_general[str(df[colname].dtype)]
                cmb = DtypeIndicator_QCombobox(column_name=colname)
                cmb.setCurrentIndex(cmb.findText(translated_dtype))
                self.formlay_dtypes.addRow(colname, cmb)
                self.covariate_cols[colname] = cmb

    def update_dataframe_dtypes(self, df):
        """
        Updates a dataframe to have any column listed as a covariate be automatically converted to the indicated dtype
        in the combobox within the dtype_indicator widget
        @param df: the dataframe whose columns may or may not be converted
        @return: the processed dataframe
        """
        current_covs = set(self.covariate_cols.keys())
        for colname in df.columns.tolist():
            if colname in current_covs:
                cmb_text = self.covariate_cols[colname].currentText()
                newtype = {"numerical": "float32", "categorical": "category"}[cmb_text]

                if cmb_text == "numerical":
                    df[colname] = df[colname].values.astype(np.float)
                elif cmb_text == "categorical":
                    df[colname] = df[colname].values.astype(np.str)
                else:
                    raise ValueError(f"update_dataframe_dtypes inferred an incorrect combobox "
                                     f"option at colname {colname}.")

                df[colname] = df[colname].astype(newtype)

        return df

    # Convenience function for purging everything in the event that either the analysis directory changes or the
    # supplementary data file provided is changed
    def clear_contents(self):
        # Clear the dicts
        self.covariate_cols.clear()
        # Clear the rows of the format layout
        for ii in range(self.formlay_dtypes.rowCount() - 1):
            self.formlay_dtypes.removeRow(1)


class DtypeIndicator_QCombobox(QComboBox):
    signal_sendupdateddtype = Signal(str, str)  # Signals the colname, general type (numerical, categorical)

    def __init__(self, column_name, parent=None):
        super(DtypeIndicator_QCombobox, self).__init__(parent=parent)
        self.associated_column = column_name
        self.addItems(["numerical", "categorical"])

    def activate(self):
        self.currentTextChanged.connect(self.inform_update)

    def deactivate(self):
        self.currentTextChanged.disconnect(self.inform_update)

    def inform_update(self, new_text):
        print(f"The combobox associated with variable {self.associated_column} is sending signal {new_text}")
        self.signal_sendupdateddtype.emit(self.associated_column, new_text)
