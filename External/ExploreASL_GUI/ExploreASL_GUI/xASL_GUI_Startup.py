from ExploreASL_GUI.xASL_GUI_MainWin import xASL_MainWin
from PySide2.QtWidgets import QApplication, QMessageBox, QWidget
from PySide2.QtGui import QIcon
from platform import system
from shutil import which
import os
import sys
import json
import subprocess
import re


def startup():
    app = QApplication(sys.argv)
    current_dir = os.getcwd()
    project_dir = current_dir
    json_logic_dir = os.path.join(current_dir, "JSON_LOGIC")
    scripts_dir = os.path.join(current_dir, "ExploreASL_GUI")
    regex = re.compile("'(.*)'")
    # Get the appropriate default style based on the user's operating system
    if system() == "Windows":  # Windows
        app.setStyle("Fusion")
    elif system() == "Darwin":  # Mac
        app.setStyle("macintosh")
    elif system() == "Linux":  # Linux
        app.setStyle("Fusion")
    else:
        QMessageBox().warning(QWidget(),
                              "Unsupported Operating System",
                              f"The operating system identified as: {system()} "
                              f"is not compatible with this program. "
                              f"Please run this program on any of the following:\n"
                              f"- Windows 10\n-Mactintosh\n-Linux using Ubuntu Kernel",
                              QMessageBox.Ok)
        sys.exit()

    if not os.path.exists(json_logic_dir):
        QMessageBox().warning(QWidget(),
                              "No JSON_LOGIC directory found",
                              f"The program directory structure is compromised. No JSON_LOGIC directory was located"
                              f"in {project_dir}",
                              QMessageBox.Ok)
        sys.exit()

    # Get the screen credentials
    screen = app.primaryScreen()
    screen_size = screen.availableSize()

    # Check if the master config file exists; if it doesn't, the app will initialize one on the first startup
    if os.path.exists(os.path.join(json_logic_dir, "ExploreASL_GUI_masterconfig.json")):
        print("Loading masterconfig file.")
        with open(os.path.join(json_logic_dir, "ExploreASL_GUI_masterconfig.json")) as master_config_reader:
            master_config = json.load(master_config_reader)

    # Otherwise, this is a first time startup and additional things need to be checked
    else:
        master_config = {"ExploreASLRoot": "",  # The filepath to the ExploreASL directory
                         "DefaultRootDir": current_dir,  # The default root for the navigator to watch from
                         "ScriptsDir": scripts_dir,  # The location of where this script is launched from
                         "ProjectDir": project_dir,  # The location of the ExploreASL_GUI main dir
                         "Platform": f"{system()}",
                         "ScreenSize": (screen_size.width(), screen_size.height()),  # Screen dimensions
                         "DeveloperMode": True  # Whether to launch the app in developer mode or not
                         }

        # We must also check for the MATLAB version present on the machine
        # First try the faster shutil.which method
        print("First time startup. Please be patient as the matlab version is detected")
        print("Attempting shutil method first")
        result = which("matlab")
        if result is not None:
            regex = re.compile(r".*R\d{4}[ab]")
            match = regex.search(result)
            if match:
                master_config["MATLABROOT"] = match.group()
                print(f"shutil method was a success and located: {match.group()}")
            else:
                QMessageBox().warning(QWidget(),
                                      "No MATLAB directory found",
                                      "No path to the MATLAB root directory could be located on this device. "
                                      "If MATLAB is installed on this device and this message is displaying, "
                                      "please contact your system administration and check whether MATLAB is "
                                      "listed in your system's PATH variable.",
                                      QMessageBox.Ok)
                sys.exit()
        # Otherwise, try the old subprocess methods
        else:
            print("shutil method unsuccessful in locating a path. Attempting backup subprocess method")
            result = subprocess.run(["matlab", "-nosplash", "-nodesktop", "-batch", "matlabroot"],
                                    capture_output=True, text=True)
            if result.returncode == 0:
                match = regex.search(result.stdout)
                if match:
                    master_config["MATLABROOT"] = match.group(1)
                    print(f"subprocess method was a success and located: {match.group(1)}")
                else:
                    QMessageBox().warning(QWidget(),
                                          "No MATLAB directory found",
                                          "No path to the MATLAB root directory could be located on this device. "
                                          "If MATLAB is installed on this device and this message is displaying, "
                                          "please contact your system administration and check whether MATLAB is "
                                          "listed in your system's PATH variable.",
                                          QMessageBox.Ok)
                    sys.exit()
            else:
                QMessageBox().warning(QWidget(),
                                      "No MATLAB directory found",
                                      "No path to the MATLAB root directory could be located on this device. "
                                      "If MATLAB 2019 or later is installed on this device and this message is still "
                                      "displaying, please contact your system administration and check whether MATLAB"
                                      "is in your system's PATH variable. Otherwise, if your MATLAB version is not "
                                      "2019a or later, then this GUI is incompatible with your version of MATLAB",
                                      QMessageBox.Ok)
                sys.exit()

    # Memory cleanup
    del regex, current_dir

    # If all was successful, launch the GUI
    app.setWindowIcon(QIcon(os.path.join(master_config["ProjectDir"], "media", "ExploreASL_logo.ico")))
    os.chdir(scripts_dir)
    main_win = xASL_MainWin(master_config)
    main_win.show()
    sys.exit(app.exec_())


if __name__ == '__main__':
    startup()
