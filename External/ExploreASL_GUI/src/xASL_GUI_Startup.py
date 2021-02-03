from src.xASL_GUI_MainWin import xASL_MainWin
from PySide2.QtWidgets import QApplication, QMessageBox, QWidget, QDialog, QDialogButtonBox
from PySide2.QtGui import QIcon
from platform import system
from shutil import which
from glob import iglob
from os import chdir
from pathlib import Path
import sys
from json import load
from typing import Union
import subprocess
import re


def get_local_matlab() -> (bool, Union[str, None]):
    matlab_cmd_path = which("matlab")
    matlab_ver_regex = re.compile(r"R\d{4}[ab]")
    # Get version #
    matlab_ver = None

    # Is on PATH
    if matlab_cmd_path is not None:

        match = matlab_ver_regex.search(matlab_cmd_path)

        # If the version number was found
        if match:
            print(f"shutil method was a success and located: {match.group()}")
            matlab_ver = match.group()
            return matlab_ver, matlab_cmd_path

        # Otherwise,
        # For Linux/MacOS with promising root
        if system() == "Linux" and '/usr/' in matlab_cmd_path:
            print(f"User clearly has the matlab command in {matlab_cmd_path}, but the version number could not be "
                  f"ascertained. Attempting to locate around '/usr/local/")
            for search_pattern in ["/usr/local/matlab**/bin", "/usr/local/**/MATLAB/*/bin"]:
                try:
                    local_result = next(iglob(search_pattern, recursive=True))
                    local_match = matlab_ver_regex.search(local_result[0])
                    if local_match:
                        matlab_ver = local_match.group()
                        return matlab_ver, matlab_cmd_path
                except StopIteration:
                    continue

        # If no luck so far, resort to using subprocess since matlab is on PATH
        if matlab_ver is None:
            print("Version was not readily visible in PATH. Attempting backup subprocess method to extract version")
            result = subprocess.run(["matlab", "-nosplash", "-nodesktop", "-batch", "matlabroot"],
                                    capture_output=True, text=True)
            match = matlab_ver_regex.search(result.stdout)
            if result.returncode == 0 and match:
                matlab_ver = match.group()
            return matlab_ver, matlab_cmd_path

        # Windows or no lucky search pattern
        return matlab_ver, matlab_cmd_path

    # Not on PATH
    else:
        if system() != "Darwin":
            return matlab_ver, matlab_cmd_path
        # MacOS, default installation to Applications seems to avoid adding MATLAB to PATH. Look for it in applications
        applications_path = Path("/Applications").resolve()
        try:
            matlab_cmd_path = str(next(applications_path.rglob("bin/matlab")))
            matlab_ver = matlab_ver_regex.search(matlab_cmd_path).group()
            return matlab_ver, matlab_cmd_path
        except (StopIteration, AttributeError):
            return matlab_ver, matlab_cmd_path


def startup():
    app = QApplication(sys.argv)
    screen = app.primaryScreen()
    screen_size = screen.availableSize()
    project_dir = Path(__file__).resolve().parent.parent
    print(f"Launching script at: {Path(__file__)} ")
    print(f"Project Directory is: {project_dir}")

    # Get the appropriate default style based on the user's operating system
    app.setStyle("Fusion") if system() in ["Windows", "Linux"] else app.setStyle("macintosh")

    # Ensure essential directories exist
    for essential_dir in ["JSON_LOGIC", "media", "External"]:
        if not (project_dir / essential_dir).exists():
            QMessageBox().warning(QWidget(),
                                  "No JSON_LOGIC directory found",
                                  f"The program directory structure is compromised. No JSON_LOGIC directory was located"
                                  f"in {project_dir}",
                                  QMessageBox.Ok)
            sys.exit(1)

    # Check if the master config file exists; if it doesn't, the app will initialize one on the first startup
    if (project_dir / "JSON_LOGIC" / "ExploreASL_GUI_masterconfig.json").exists():
        print("Loading masterconfig file.")
        with open(project_dir / "JSON_LOGIC" / "ExploreASL_GUI_masterconfig.json") as master_config_reader:
            master_config = load(master_config_reader)
        # Update the ProjectDir and ScriptsDir variables in the event the user moves the location of this folder
        # First, make sure the Startup.py is located in the src folder
        if project_dir != master_config["ProjectDir"]:
            master_config["ProjectDir"] = str(project_dir)
            master_config["ScriptsDir"] = str(project_dir / "src")

    # Otherwise, this is a first time startup and additional things need to be checked
    else:
        master_config = {"ExploreASLRoot": "",  # The filepath to the ExploreASL directory
                         "DefaultRootDir": str(Path.home()),  # The default root for the navigator to watch from
                         "ScriptsDir": str(project_dir / "src"),  # The location of where this script is launched from
                         "ProjectDir": str(project_dir),  # The location of the src main dir
                         "Platform": f"{system()}",
                         "ScreenSize": (screen_size.width(), screen_size.height()),  # Screen dimensions
                         "DeveloperMode": True}  # Whether to launch the app in developer mode or not

        # We must also check for the MATLAB version present on the machine
        desc = "Is a standard MATLAB program installed on this machine?"
        check_for_local = QMessageBox.question(QWidget(), "MATLAB Detection", desc,
                                               (QMessageBox.Yes | QMessageBox.No | QMessageBox.Cancel))
        if check_for_local == QMessageBox.Cancel:
            sys.exit(0)

        if check_for_local == QMessageBox.Yes:
            version, cmd_path = get_local_matlab()
            master_config["MATLAB_VER"] = version
            master_config["MATLAB_CMD_PATH"] = cmd_path
            if cmd_path is None:
                QMessageBox.warning(QWidget(), "No MATLAB Command Found",
                                    "The matlab command could not be found on this local system. ASL analysis "
                                    "cannot be locally performed without the GUI knowing from where MATLAB can be "
                                    "launched.", QMessageBox.Ok)
            elif cmd_path is not None and version is None:
                QMessageBox.warning(QWidget(), "No MATLAB Version Found",
                                    f"The matlab command was found at:\n{str(cmd_path)}.\nHowever, the version could "
                                    f"not be determined. The GUI requires knowledge of the MATLAB version being run "
                                    f"in order to process data", QMessageBox.Ok)
            else:
                QMessageBox.information(QWidget(), "Local MATLAB Located && Version discerned",
                                        f"Detected the matlab path to be: {cmd_path}\n"
                                        f"Detected the matlab version to be: {version}", QMessageBox.Ok)
        else:
            QMessageBox.information(QWidget(), "No off-local support at the current time",
                                    "The current version of ExploreASL_GUI does not offer support for compiled or "
                                    "virtual MATLAB runtimes.", QMessageBox.Ok)

        # Assuming the above was successful, dcm2niix may not have executable permission; add execute permissions
        dcm2niix_dir = project_dir / "External" / "DCM2NIIX" / f"DCM2NIIX_{system()}"
        dcm2niix_file = next(dcm2niix_dir.glob("dcm2niix*"))
        stat = oct(dcm2niix_file.stat().st_mode)
        if not stat.endswith("775"):
            dcm2niix_file.chmod(0o775)
        else:
            print(f"dcm2niix already has execute permissions")

    # If all was successful, launch the GUI
    app.setWindowIcon(QIcon(str(project_dir / "media" / "ExploreASL_logo.ico")))
    chdir(project_dir / "src")

    # Memory cleanup
    del project_dir

    main_win = xASL_MainWin(master_config)
    main_win.show()
    sys.exit(app.exec_())


if __name__ == '__main__':
    startup()
