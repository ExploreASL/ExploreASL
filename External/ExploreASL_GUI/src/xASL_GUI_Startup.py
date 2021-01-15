from src.xASL_GUI_MainWin import xASL_MainWin
from PySide2.QtWidgets import QApplication, QMessageBox, QWidget, QDialog
from PySide2.QtGui import QIcon
from platform import system
from shutil import which
from glob import glob
from os import chdir
from pathlib import Path
import sys
from json import load
import subprocess
import re


def startup():
    print(f"Launching script at: {Path(__file__)} ")
    app = QApplication(sys.argv)
    project_dir = Path(__file__).resolve().parent.parent
    print(f"Project Directory is: {project_dir}")
    # Get the appropriate default style based on the user's operating system
    if system() in ["Windows", "Linux"]:  # Windows
        app.setStyle("Fusion")
    else:  # Mac
        app.setStyle("macintosh")

    if not (project_dir / "JSON_LOGIC").exists():
        QMessageBox().warning(QWidget(),
                              "No JSON_LOGIC directory found",
                              f"The program directory structure is compromised. No JSON_LOGIC directory was located"
                              f"in {project_dir}",
                              QMessageBox.Ok)
        sys.exit(1)

    # Get the screen credentials
    screen = app.primaryScreen()
    screen_size = screen.availableSize()

    if (project_dir / "JSON_LOGIC" / "ErrorsListing.json").exists():
        with open(project_dir / "JSON_LOGIC" / "ErrorsListing.json") as startup_errs_reader:
            startup_errs = load(startup_errs_reader)
    else:
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
        # Two common possibilities on Linux:
        # 1) the PATH registered is /usr/local/MATLAB/[VERSION]/bin
        # 2) the PATH registered is /usr/bin
        # Try option 1 first to get the version, opt for looking around /usr/local if Option 2 comes around
        print("First time startup. Please be patient as the matlab version is detected\n"
              "Attempting shutil method first")

        # Is matlab even on PATH?
        result = which("matlab")
        if result is not None:
            # If matlab was found on PATH, see if the version number can be extracted from it as per possibility #1
            # On Windows, this should always work. On Linux...*sigh*
            regex = re.compile(r".*R\d{4}[ab]")
            match = regex.search(result)
            if match:
                master_config["MATLABROOT"] = match.group()
                print(f"shutil method was a success and located: {match.group()}")

            # Not a success, see if the version number can be extracted from possibility #2
            elif not match and '/usr/' in result:
                version_is_located = False
                print(f"User clearly has the matlab command in {result}, but the version number could not be "
                      f"ascertained. Attempting to locate around '/usr/local/")
                for search_pattern in ["/usr/local/matlab**/bin", "/usr/local/**/MATLAB/*/bin"]:
                    local_result = glob(search_pattern, recursive=True)
                    if len(local_result) != 0:
                        local_match = regex.search(local_result[0])
                        if local_match:
                            print(f"Located MATLAB version number: {local_match.group()}")
                            master_config["MATLABROOT"] = local_match.group()
                            version_is_located = True
                            break

                if not version_is_located:
                    QMessageBox().warning(QWidget(), startup_errs["NoMATLABVer"][0],
                                          startup_errs["NoMATLABVer"][1], QMessageBox.Ok)
                    sys.exit(1)

            # Neither a success with possibility #1 or #2
            else:
                QMessageBox().warning(QWidget(), startup_errs["NoMATLABVer"][0],
                                      startup_errs["NoMATLABVer"][1], QMessageBox.Ok)
                sys.exit(1)

            # Assuming the above was successful, dcm2niix may not have executable permission; add execute permissions
            dcm2niix_dir = project_dir / "External" / "DCM2NIIX" / f"DCM2NIIX_{system()}"
            dcm2niix_file = next(dcm2niix_dir.glob("dcm2niix*"))
            stat = oct(dcm2niix_file.stat().st_mode)
            if not stat.endswith("775"):
                dcm2niix_file.chmod(0o775)
            else:
                print(f"dcm2niix already has execute permissions")

        # Otherwise, try the old subprocess methods
        else:
            print("shutil method unsuccessful in locating a path. Attempting backup subprocess method")
            regex = re.compile("'(.*)'")
            result = subprocess.run(["matlab", "-nosplash", "-nodesktop", "-batch", "matlabroot"],
                                    capture_output=True, text=True)
            if result.returncode == 0:
                match = regex.search(result.stdout)
                if match:
                    master_config["MATLABROOT"] = match.group(1)
                    print(f"subprocess method was a success and located: {match.group(1)}")
                else:
                    QMessageBox().warning(QWidget(), startup_errs["NoMATLAB"][0],
                                          startup_errs["NoMATLAB"][1], QMessageBox.Ok)
                    sys.exit(1)
            else:
                QMessageBox().warning(QWidget(), startup_errs["NoMATLAB"][0],
                                      startup_errs["NoMATLAB"][1], QMessageBox.Ok)
                sys.exit(1)

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
