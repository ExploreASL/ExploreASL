from PySide2.QtWidgets import *
import pandas as pd
import os
from glob import glob


class xASL_GUI_Data_Loader:
    """
    Class specifically dedicated for loading in the ExploreASL data from its Stats directory as well as loading in any
    ancillary data the user may provide
    """

    def __init__(self, parent):
        self.parent_cw = parent
        self.loaded_wide_data = pd.DataFrame()
        self.loaded_long_data = pd.DataFrame()
        self.dtype_guide = {"SUBJECT": "object",
                            "LongitudinalTimePoint": "category",
                            "SubjectNList": "category",
                            "Site": "category",
                            "AcquisitionTime": "float64",
                            "GM_vol": "float64",
                            "WM_vol": "float64",
                            "CSF_vol": "float64",
                            "GM_ICVRatio": "float64",
                            "GMWM_ICVRatio": "float64"}

    def load_exploreasl_data(self):
        # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        # First Section - Load in the ExploreASL Stats directory data
        if os.path.exists(os.path.join(self.parent_cw.le_analysis_dir.text(), "Population", "Stats")):
            stats_dir = os.path.join(self.parent_cw.le_analysis_dir.text(), "Population", "Stats")
            atlas = {"MNI": "MNI_structural", "Hammers": "Hammers"}[self.parent_cw.cmb_atlas_selection.currentText()]
            pvc = {"With Partial Volume Correction": "PVC2",
                   "Without Partial Volume Correction": "PVC0"}[self.parent_cw.cmb_pvc_selection.currentText()]
            stat = {"Mean": "mean", "Median": "median",
                    "Coefficient of Variation": "CoV"}[self.parent_cw.cmb_stats_selection.currentText()]
            # Get the relevant files
            gm_file = glob(os.path.join(stats_dir, f'{stat}_*_TotalGM*{pvc}.tsv'))
            wm_file = glob(os.path.join(stats_dir, f'{stat}_*_DeepWM*{pvc}.tsv'))
            atlas_file = glob(os.path.join(stats_dir, f'{stat}_*_{atlas}*{pvc}.tsv'))
            # Exit if not all files can be found
            for file in [gm_file, wm_file, atlas_file]:
                if len(file) == 0:
                    return
            # Clearing of appropriate widgets to accomodate new data
            self.parent_cw.lst_varview.clear()
            # Extract each as a dataframe and merge them
            dfs = []
            for file in [gm_file, wm_file, atlas_file]:
                df = pd.read_csv(file[0], sep='\t')
                df.drop(0, axis=0, inplace=True)  # First row is unnecessary
                df = df.loc[:, [col for col in df.columns if "Unnamed" not in col]]
                dfs.append(df)
            df: pd.DataFrame = pd.concat(dfs, axis=1)
            df = df.T.drop_duplicates().T

            # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            # Second Section - Fix the ExploreASL native data dtypes
            for col in df.columns:
                if col in self.dtype_guide.keys():
                    df[col] = df[col].astype(self.dtype_guide[col])
                else:
                    df[col] = df[col].astype("float64")
            self.loaded_wide_data = df
            self.backup_data = df.copy()

            # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            # Third Section - If there is any ancillary data specified, load it in
            if all([os.path.exists(self.parent_cw.le_demographics_file.text()),
                    os.path.isfile(self.parent_cw.le_demographics_file.text()),
                    os.path.splitext(self.parent_cw.le_demographics_file.text())[1] in [".tsv", ".csv", ".xlsx"]
                    ]):
                result = self.load_ancillary_data(df)
                if result is not None:
                    self.loaded_wide_data = result
                # If the merging failed, default to just using the ExploreASL datasets. In a future update, add some
                # sort of user feedback that this went wrong
                else:
                    self.loaded_wide_data = self.backup_data

            # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            # Fourth Section - Convert the wide format data into a long format
            vars_to_keep_constant = [col for col in self.loaded_wide_data.columns if not any([col.endswith("_B"),
                                                                                              col.endswith("_L"),
                                                                                              col.endswith("_R")])]
            vars_to_melt = [col for col in self.loaded_wide_data.columns if col not in vars_to_keep_constant]
            self.loaded_long_data = self.loaded_wide_data.melt(id_vars=vars_to_keep_constant,
                                                               value_vars=vars_to_melt,
                                                               var_name="Atlas Location",
                                                               value_name="CBF")
            self.loaded_long_data["CBF"] = self.loaded_long_data["CBF"].astype("float64")
            atlas_location = self.loaded_long_data.pop("Atlas Location")
            atlas_loc_df: pd.DataFrame = atlas_location.str.extract("(.*)_(B|L|R)", expand=True)
            atlas_loc_df.rename(columns={0: "Anatomical Area", 1: "Side of the Brain"}, inplace=True)
            atlas_loc_df["Side of the Brain"] = atlas_loc_df["Side of the Brain"].apply(lambda x: {"B": "Bilateral",
                                                                                                   "R": "Right",
                                                                                                   "L": "Left"}[x])
            atlas_loc_df = atlas_loc_df.astype("category")
            self.loaded_long_data: pd.DataFrame = pd.concat([self.loaded_long_data, atlas_loc_df], axis=1)
            self.loaded_long_data = self.loaded_long_data.infer_objects()
            self.current_dtypes = self.loaded_long_data.dtypes
            self.current_dtypes = {col: str(str_name) for col, str_name in
                                   zip(self.current_dtypes.index, self.current_dtypes.values)}
            self.parent_cw.lst_varview.addItems(self.loaded_long_data.columns.tolist())

            # The user may have loaded in new data and the subsetter's fields should reflect that
            self.parent_cw.subsetter.update_subsetable_fields()

            # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            # Fifth Section - Subset the data accordingly if the criteria is set
            self.loaded_long_data = self.parent_cw.subsetter.subset_data(self.loaded_long_data)

            # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            # Sixth Section - Housekeeping and Finishing touches
            # Alter this when Section 5 is completed; long_data is the "good copy" of the data that will be plotted
            self.long_data = self.loaded_long_data
            self.parent_cw.cmb_figuretypeselection.setEnabled(
                True)  # Data is loaded; figure selection settings can be enabled
            self.parent_cw.btn_subset_data.setEnabled(True)  # Data is loaded; subsetting is allowed

            # In case any of this was done again (data was already loaded once before), we must account for what may
            # have already been plotted or set; everything must be cleared. This should be as easy as setting the
            # figureselection to the first index, as plots & settings can only exist if its current index is non-zero,
            # and setting it to zero has the benefit of clearing everything else already
            if self.parent_cw.cmb_figuretypeselection.currentIndex() != 0:
                self.parent_cw.cmb_figuretypeselection.setCurrentIndex(0)

    def load_ancillary_data(self, exasl_df):
        # Load in the other dataframe, with flexibility for filetype
        file = self.parent_cw.le_demographics_file.text()
        filetype = os.path.splitext(file)[1]
        if filetype == '.tsv':
            demo_df = pd.read_csv(file, sep='\t')
        elif filetype == '.csv':
            demo_df = pd.read_csv(file)
        elif filetype == '.xlsx':
            demo_df = pd.read_excel(file)
        else:
            print("An unsupported filetype was given")
            QMessageBox().warning(self.parent_cw,
                                  "Unsupported File Type",
                                  "This program only accepts the following filetypes:\n"
                                  "Comma-separated values (*.csv)\n"
                                  "Tab-separated values (*.tsv)\n"
                                  "Microsoft Excel spreadsheets (*.xlsx)",
                                  QMessageBox.Ok)
            return None
        # Abort if the pertinent "SUBJECT" column is not in the read columns. In a future update, add support for user
        # specification of which column to interpret as the SUBJECT column
        if "SUBJECT" not in demo_df.columns:
            return None
        merged = pd.merge(left=demo_df, right=exasl_df, how='inner', on='SUBJECT', sort=True)
        if len(merged) == 0:
            return None
        sub_in_merge, sub_in_demo, sub_in_exasl = set(merged["SUBJECT"].tolist()), set(
            demo_df["SUBJECT"].tolist()), set(exasl_df["SUBJECT"].tolist())
        diff_in_demo = sub_in_demo.difference(sub_in_merge)
        diff_in_exasl = sub_in_exasl.difference(sub_in_merge)
        if any([len(diff_in_demo) > 0, len(diff_in_exasl) > 0]):
            QMessageBox().information(self.parent_cw,
                                      "Merge successful, but differences were found:\n",
                                      f"You provided a file with {len(sub_in_demo)} subjects.\n"
                                      f"ExploreASL's output had {len(sub_in_exasl)} subjects.\n"
                                      f"During the merge {len(diff_in_demo)} subjects present in the file "
                                      f"had to be excluded.\n"
                                      f"During the merge {len(diff_in_exasl)} subjects present in ExploreASL's output "
                                      f"had to be excluded",
                                      QMessageBox.Ok)
        return merged
