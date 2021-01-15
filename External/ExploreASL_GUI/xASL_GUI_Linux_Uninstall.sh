#!/bin/bash
EXPLOREASLROOT="$( cd "$(dirname "$0")" >/dev/null 2>&1 ; pwd -P )"
echo "Sudo password required to remove ExploreASL_GUI.desktop from /usr/share/applications"
rm -rf $EXPLOREASLROOT || exit 1
echo "Successfully removed the contents of {$EXPLOREASLROOT}"
sudo rm /usr/share/applications/ExploreASL_GUI.desktop || exit 1
echo "Removed ExploreASL_GUI.desktop file successfully. Uninstall complete."
