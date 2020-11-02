## Summary
The purpose of this folder is to demonstrate how to schedule xASL jobs on a linux cluster with remote access (SSH). To this end, credentials have to be provided and a task has to be set up with a scheduler.

## Requirements
Execution has been tested under Ubuntu 20.04.1 LTS using python 3.8.5 with the packages listed in requirements.txt. The execution relies on the following linux packages:
"""
sudo apt install sshpass
sudo apt install vpnc
sudo apt install python3-pip
sudo pip3 install keyrings.alt
"""

## Execution
To set up a scheduled task, please follow the steps below:
1; Define your credentials by typing
"""
python3 define_credentials.py
"""
into the terminal and fill in your details. These details can be checked and removed any time with the check_credentials.py and remove_credentials.py scripts which can be executed similarly. (Selecting "all" in the prompt of remove_credentials.py will delete every confidential information. Furthermore, the passwords stored by the keyring python module cannot be accesed by other users even with root permission.)

2; If VPN connection is required then copy the vpn_sample.conf file and edit it to fill in your details.
"""
cp vpn_sample.conf .my_vpn.conf
nano .my_vpn.conf
sudo chown root .my_vpn.conf
sudo chgrp root .my_vpn.conf
"""
**Warning: setting up a VPN require sudo rights! Storing credentials in the *.conf file is insecure**

3; On the cluster, commands included in the cluster_commands.sh file will be executed. This should include submitting a job file, such as the xASL_job.run file which is processed by the slurm load manager. Please edit the xASL_job.run file to specify the expected execution time and potentially an email address for updates related to the job.

4; Run "python3 job_submission.py" to test that everything is set up correctly. This script copies data to the user's home directory on the cluster. If this is not suitable then the code should be edited.

5; New tasks can be scheduled by typing
"""
crontab -e
"""
in a linux terminal. This leads to an editable file specifying scheduled tasks. Adding, for instance, the following line means that the job_submission.py script is executed at noon every sunday.
> 0 12 * * 7 python3 /[path to working directory]/job_submission.py


### ToDo
- [ ] Find a more secure solution to store VPN login details
- [ ] Verification on test data

(The runner4scheduler.sh file became useless.)
