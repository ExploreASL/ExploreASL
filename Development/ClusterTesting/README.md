## Summary
The purpose of this folder is to demonstrate how to schedule xASL jobs on a linux cluster with remote access (SSH). To this end, credentials have to be provided and a task has to be set up with a scheduler.

## Requirements
Execution has been tested under Ubuntu 20.04.1 LTS using python 3.8.5 with the packages listed in requirements.txt. The required linux packages can be installed with:

```
sudo apt install sshpass
sudo apt install vpnc
sudo apt install python3-pip
sudo pip3 install --upgrade pip
sudo pip3 install keyring --upgrade
sudo pip3 install keyrings.alt
```

## Execution
To set up a scheduled task, please follow the steps below:

1; Define your credentials by typing
```
python3 define_credentials.py
```
into the terminal and fill in your details. These details can be checked, modified, and removed any time with the check_credentials.py, modify_credentials.py, and the remove_credentials.py scripts which can be executed similarly. (Selecting "all" in the prompt of remove_credentials.py will delete every confidential information. Passwords stored by the keyring python module cannot be directly accesed by other users.)

**Warning: setting up a VPN require sudo rights! Storing credentials in the *.conf file is insecure**

2; On the cluster, commands included in the cluster_commands.sh file will be executed. This includes submitting the xASL_job.run file which should be created by copying the xASL_job_template.run and will be processed by the slurm load manager. The submitted job executes the xASL_test.m MATLAB file. Please edit the xASL_job.run file to specify the expected (overestimated) execution time and potentially an email address for updates related to the job.

3; Run
```
python3 job_submission.py
```
to test that everything is set up correctly. This script copies data to the user's home directory on the cluster. If this is not suitable then the code should be edited.

4; New tasks can be scheduled by typing
```
crontab -e
```
in a linux terminal. This leads to an editable file specifying scheduled tasks. Adding, for instance, the following line means that the job_submission.py script is executed at noon every sunday.
```
0 12 * * 7 python3 /[path to working directory]/job_submission.py
```

### ToDo
- [ ] Find more secure solution to store credentials
- [ ] Verification on test data
