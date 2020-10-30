import os
import sys
import re
import subprocess as sp
from glob import glob
from platform import system
import shutil
from getpass import getpass


def Linux_Install():
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
            installer_name = "ExploreASL_GUI_Linux_Install.py"
            explore_asl_root = os.path.dirname(glob(os.path.join(os.getcwd(), "**", installer_name), recursive=True)[0])
            if not os.path.isdir(explore_asl_root):
                sys.exit(1)
            else:
                os.chdir(explore_asl_root)
    except IndexError:
        sys.exit(1)
    print(f"{explore_asl_root=}")

    # Make the venv directory if it does not exist
    if not os.path.isdir(os.path.join(explore_asl_root, "venv")):
        print("Making the venv directory...")
        os.mkdir(os.path.join(explore_asl_root, 'venv'))

    # Create the venv
    packages_path = os.path.join(explore_asl_root, "venv", "lib", "python3.8", "site-packages")
    if not os.path.exists(packages_path):
        venvcreate_result = sp.run(f"python3.8 -m venv {explore_asl_root}/venv; ", shell=True, text=True)
        print(f"Return code of creating a new environment: {venvcreate_result.returncode=} (0 means successful)")
        if venvcreate_result.returncode != 0:
            sys.exit(venvcreate_result.returncode)
    else:
        print("Venv was already created. Continuing...")

    # Activate the venv, then install the packages
    bin_path = os.path.join(explore_asl_root, "venv", "bin")
    activate_command = f". {explore_asl_root}/venv/bin/activate"
    packages_install_command = f"pip3 install -r {explore_asl_root}/requirements.txt"
    if not os.path.isdir(bin_path):
        print("The venv/bin path does not exist!")
        sys.exit(1)
    if not os.path.isfile(os.path.join(explore_asl_root, "requirements.txt")):
        print(f"Could not find the requirements.txt file within {explore_asl_root}")
        sys.exit(1)
    if not os.path.exists(os.path.join(explore_asl_root, "venv", "lib", "python3.8", "site-packages", "nilearn")):
        venvactandinstall_result = sp.run(f"{activate_command} && {packages_install_command}", shell=True, text=True)
        print(f"Result code of activating the env and installing packages {venvactandinstall_result.returncode=} "
              f"(0 means successful)")
        if venvactandinstall_result.returncode != 0:
            print(venvactandinstall_result.stderr)
            sys.exit(venvactandinstall_result.returncode)
    else:
        print("Packages already exist in venv site-packages. Continuing...")

    # Prepare the shell launcher
    shell_launcher = os.path.join(explore_asl_root, "xASL_GUI_run.sh")
    if not os.path.exists(shell_launcher):
        with open(shell_launcher, 'w') as shell_writer:
            to_write = f"""#!/bin/bash
. {explore_asl_root}/venv/bin/activate
cd {explore_asl_root}
python3.8 {explore_asl_root}/xASL_GUI_run.py
"""
            shell_writer.write(to_write)
        os.chmod(shell_launcher, 0o744)
    else:
        print(f"{shell_launcher} already exists. Continuing...")

    # Prepare the .desktop file
    desktop_file = os.path.join(explore_asl_root, "ExploreASL_GUI.desktop")
    if not os.path.exists(desktop_file):
        with open(desktop_file, "w") as desktop_writer:
            to_write = f"""[Desktop Entry]
Version=0.2.5
Name=ExploreASL_GUI
Comment=ExploreASL_GUI is an accompanying graphical user interface for ExploreASL workflows
Exec={shell_launcher}
Icon={explore_asl_root}/media/ExploreASL_logo.png
Terminal=false
Type=Application
Categories=Utility;Development;
"""
            desktop_writer.write(to_write)
    else:
        print(f"{desktop_file} already exists. Continuing...")

    # Move the .desktop file into /usr/share/applications
    password = getpass("To move a .desktop file to /usr/local/applications , we will only as for your "
                       "\nsudo password this once: ")
    proc = sp.Popen(f"sudo -S mv {desktop_file} /usr/share/applications/ExploreASL_GUI.desktop", shell=True,
                    stdin=sp.PIPE, stdout=sp.PIPE, stderr=sp.PIPE, universal_newlines=True)
    proc.communicate(password)
    print(f"Result code of moving the desktop file to /usr/local/applications: {proc.returncode} (0 means successful)")
    if proc.returncode != 0:
        print(f"{proc.returncode=}")
        sys.exit(proc.returncode)

    print("Installation complete. If you are seeing this message, everything should be a success")
    sys.exit(0)


def MacOS_Install():
    sys.exit(1)


if __name__ == '__main__':
    if system() == "Linux":
        Linux_Install()
    elif system() == "Darwin":
        MacOS_Install()
    else:
        print("Neither Linux nor Mac System described")
        sys.exit(1)
