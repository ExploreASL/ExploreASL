from pathlib import Path
from json import load, dump, JSONDecodeError
from typing import Union, List, Any
import pandas as pd
from numpy import isnan


########################################################################################################################
# PREFACE
# This module contains functions meant to ease with the operation of directories. Merging directories; establishing
# symlinks; coordinated changing for filedata; etc.
# Current Main Functions:
#       - alter_sidecars ; using either a csv dataframe or a list of subjects + key + value ; alter the json sidecars
#       in a given study
########################################################################################################################
def robust_read_csv(df_path: Union[Path, str], **kwargs):
    """
    Reads in a dataframe, accounting for common file extensions

    :param df_path: Path object to the
    :param kwargs: keyword arguments to feed the expected pandas import function
    :return: the loaded in dataframe
    """
    if isinstance(df_path, str):
        df_path = Path(df_path)

    if df_path.suffix == ".csv":
        df = pd.read_csv(df_path, sep=",", **kwargs)
    elif df_path.suffix == ".tsv":
        df = pd.read_csv(df_path, sep="\t", **kwargs)
    elif df_path.suffix == ".xlsx":
        df = pd.read_excel(df_path, **kwargs)
    else:
        return None
    for column in df.columns:
        if column.startswith("Unnamed:"):
            df.drop(column, axis=1, inplace=True)
    return df


def interpret_value(value: Any):
    if isinstance(value, (int, float, list, dict, set)):
        return value
    if value == "":  # Case: empty string; just return it
        return ""
    elif value.isdigit():  # Case: something numerical; return the appropriate number type
        try:
            return int(value)
        except ValueError:
            return float(value)
    elif value.startswith("[") and value.endswith("]"):
        splitter = ", " if ", " in value else ","
        to_return = [interpret_value(sub_val) for sub_val in value.strip("[]").split(splitter)]  # Case: list;
        if len(to_return) == 1:
            return to_return[0]
        else:
            return to_return
    elif value.lower() in ["true", "t", "yes", "y"]:  # Case: it's a positive boolean
        return True
    elif value.lower() in ["false", "f", "no", "n"]:  # Case: it's a negative boolean
        return False
    else:
        return value


def alter_json_sidecar(json_path: Union[Path, str], action: str, key: str, value: Any = None):
    """
    Changes a key within the json sidecars to either have a specific key removed altogether or its value changed
    :param json_path: path to the json sidecar file
    :param action: a string denoting the action to take ("remove" to remove the key; "alter" to alter the value)
    :param key: which key to remove or change
    :param value: if altering a key, what new value it should take on
    """
    json_path = Path(json_path)
    if any([not json_path.exists(), key is None]):
        return
    try:
        with open(json_path, "r") as sidecar_reader:
            sidecar_data = load(sidecar_reader)
            if action.lower() in ["remove", "purge", "delete"]:
                del sidecar_data[key]
            else:
                sidecar_data[key] = value
    except KeyError as key_err:
        print(f"Encountered a KeyError with key {key}.\n{key_err}")
        return
    except JSONDecodeError as json_err:
        print(f"Encountered a JSONDecodeError with with file {json_path}\n{json_err}")
        return
    with open(json_path, "w") as sidecar_writer:
        dump(sidecar_data, sidecar_writer, indent=3)


def alter_sidecars(root_dir: Union[str, Path], subjects: Union[List[str], str, Path],
                   which_scan: str, action: str, key: str = None, value: Any = None):
    """
    Changes the json sidecars of specified subjects in a study directory

    :param root_dir: The root directory from which subjects will be searched for
    :param subjects: Either a list of subject names or a csv file whose first column SUBJECT is the list of subjects
    and whose second column whose name is the Key to alter/remove and whose values are new values to implement, if any
    :param which_scan: one of "ASL", "T1" or "M0"; defines which sidecars will be targetted
    :param action: one of "remove" or "alter"; defines whether the key in the sidecar is changed or deleted
    :param key: the name of the sidecar key to change
    :param value: the new value the sidecar should take on, if any
    """
    # Defensive Programming
    root_dir = Path(root_dir)
    if not root_dir.exists():
        raise ValueError("root_dir argument was provided a filepath that does not exist")
    if which_scan.lower() not in ["asl", "t1", "m0"]:
        raise ValueError("which_scan argument must be one of 'ASL' or 'T1' or 'M0'")
    else:
        scan_translator = {"asl": "*ASL*.json", "t1": "*T1*.json", "m0": "*M0*.json"}

    # Get a dict whose keys are subject names and whose values are a dict of sidecar key and new value
    if isinstance(subjects, (str, Path)):
        df = robust_read_csv(subjects)
        df.set_index(df.columns[0], inplace=True)
        iter_dict: dict = df.T.to_dict()
    elif isinstance(subjects, list):
        iter_dict: dict = {subject: {key: value} for subject in subjects}
    else:
        raise TypeError("The subjects argument did not fit into any of the expected types.\n"
                        "Please provide either a file with a column of names and a column of values\n"
                        "or provide a list of strings representing subjects as well as the appropriate\n"
                        "key and values that will be applied to all")

    for subject, key_val_dict in iter_dict.items():
        try:
            subject_path = next(root_dir.rglob(subject))
            if not subject_path.is_dir():
                continue
            print(subject_path)
        except StopIteration:
            continue
        for key, value in key_val_dict.items():
            interpreted_value = interpret_value(value)
            try:
                json_path = next(subject_path.rglob(scan_translator[which_scan.lower()]))
            except StopIteration:
                continue
            if action.lower() in ["remove", "purge", "delete"]:
                alter_json_sidecar(json_path=json_path, action=action, key=key, value=interpreted_value)
            else:
                if isinstance(interpreted_value, (list, dict, str)):
                    alter_json_sidecar(json_path=json_path, action=action, key=key, value=interpreted_value)
                elif not isnan(interpreted_value):
                    alter_json_sidecar(json_path=json_path, action=action, key=key, value=interpreted_value)
                else:
                    continue


def merge_directories(roots: List[Union[Path, str]],
                      merge_root: Union[Path, str],
                      symbolic: bool = True,
                      overwrite: bool = False):
    """
    Merges several root directories into a single study symbolically or by copying filepaths over

    :param roots: list of study root paths (i.e. /home/mpasternak/GENFI_utils/Philips_2D_EPI/analysis)
    :param merge_root: path to the new root that will contain all the merged data
    :param symbolic: whether to make the new paths symbolic (True) or real copies (False); default is True
    :param overwrite: whether to overwrite existent paths downstream from merge_root; default is False
    """

    def make_symlinks(parentpath: Path, current_root: str, target_root: str, overwrite_links: bool = overwrite):
        """
        Walks recursively down a filepath, creating symlinks as needed based on a reference "other" root directory
        """
        for childpath in parentpath.iterdir():
            symlink_path = Path(str(childpath).replace(current_root, target_root))
            if all([childpath.exists(), childpath.is_dir()]):
                symlink_path.mkdir(parents=True, exist_ok=True)
                make_symlinks(childpath, current_root, target_root, overwrite_links)
            elif all([childpath.exists(), childpath.is_file(), overwrite_links]):
                if symlink_path.exists():
                    symlink_path.unlink(missing_ok=True)
                symlink_path.symlink_to(childpath)
            elif all([childpath.exists(), childpath.is_file(), not overwrite_links]):
                if not symlink_path.exists():
                    symlink_path.symlink_to(childpath)
            else:
                continue

    def make_reallinks(parentpath: Path, current_root: str, target_root: str, overwrite_links: bool = overwrite):
        """
        Walks recursively down a filepath, copying files & dirs as needed based on a reference "other" root directory
        """
        for childpath in parentpath.iterdir():
            reallink_path = Path(str(childpath).replace(current_root, target_root))
            if all([childpath.exists(), childpath.is_dir()]):
                reallink_path.mkdir(parents=True, exist_ok=True)
                make_reallinks(childpath, current_root, target_root, overwrite_links)
            elif all([childpath.exists(), childpath.is_file(), overwrite_links]):
                if reallink_path.exists():
                    reallink_path.unlink(missing_ok=True)
                copyfile(src=childpath, dst=reallink_path)
            elif all([childpath.exists(), childpath.is_file(), not overwrite_links]):
                if not reallink_path.exists():
                    copyfile(src=childpath, dst=reallink_path)
            else:
                continue

    # Account for typing
    if isinstance(merge_root, str):
        merge_root = Path(merge_root)
    if not merge_root.exists():
        merge_root.mkdir(parents=True)

    # Iterate over the analysis directories
    for root in roots:
        if isinstance(root, str):
            root = Path(root)
        if symbolic:
            make_symlinks(parentpath=root, current_root=str(root),
                          target_root=str(merge_root), overwrite_links=symbolic)
        else:
            make_reallinks(parentpath=root, current_root=str(root),
                           target_root=str(merge_root), overwrite_links=symbolic)
