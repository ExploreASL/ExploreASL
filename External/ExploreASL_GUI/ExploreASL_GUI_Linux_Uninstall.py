import subprocess as sp
from getpass import getpass
import os
import sys
from platform import system
from glob import glob


def Linux_Uninstall():
    print("Implementing the Linux Install")
    # Get the root directory and the location of the bash directory
    try:
        if any([os.getcwd().endswith("ExploreASL_GUI/ExploreASL_GUI"),
                os.getcwd().endswith("ExploreASL_GUI-master/ExploreASL_GUI")]):
            explore_asl_root = os.path.dirname(os.getcwd())
            if not os.path.isdir(explore_asl_root):
                return 1
            else:
                os.chdir(explore_asl_root)
        else:
            uninstaller_name = "ExploreASL_GUI_Linux_Uninstall.py"
            explore_asl_root = os.path.dirname(glob(os.path.join(os.getcwd(), "**", uninstaller_name),
                                                    recursive=True)[0])
            if not os.path.isdir(explore_asl_root):
                sys.exit(1)
            else:
                os.chdir(explore_asl_root)
    except IndexError:
        sys.exit(1)
    print(f"{explore_asl_root=}")

    # Get rid of the .desktop file
    password = getpass("To remove the .desktop file from /usr/local/applications and all the ExploreASL_GUI "
                       "directories, your sudo password is required: ")
    proc = sp.Popen(f"sudo -S rm /usr/share/applications/ExploreASL_GUI.desktop; sudo rm -rf {explore_asl_root}",
                    shell=True, stdin=sp.PIPE, stdout=sp.PIPE, stderr=sp.PIPE, universal_newlines=True)
    proc.communicate(password)
    print(f"Result code of removing the desktop file from /usr/local/applications: "
          f"{proc.returncode} (0 means successful)")
    if proc.returncode != 0:
        print(f"{proc.returncode=}")
        sys.exit(proc.returncode)
    else:
        print("ExploreASL_GUI has successfully uninstalled. Exiting.")
        sys.exit(0)


def MacOS_Uninstall():
    pass


if __name__ == '__main__':
    if system() == "Linux":
        Linux_Uninstall()
        sys.exit(0)
    elif system() == "Darwin":
        MacOS_Uninstall()
        sys.exit(0)
    else:
        print("Neither Linux nor Mac System described")
        sys.exit(1)
