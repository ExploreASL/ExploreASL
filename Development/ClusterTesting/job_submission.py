import getpass
import keyring
import os
import sys

abspath = os.path.abspath(__file__)
dname = os.path.dirname(abspath)
os.chdir(dname)

# check whether files containing vpn and arcus username exist
username_file = "./.login_srnm"
if os.path.exists(username_file):
    f = open(username_file, "r")
    username = f.read()
    f.close()
else:
    print('PERMISSION DENIED.\n\n')
    sys.exit()
        
cluster_ddrss  = keyring.get_password('ssh_address',username)
cluster_psswrd = keyring.get_password('ssh_password',username)
sudo_psswrd    = keyring.get_password('sudo_password',username)
vpn_address    = keyring.get_password('vpn_address',username)
vpn_username   = keyring.get_password('vpn_username',username)
vpn_password   = keyring.get_password('vpn_password',username)
vpn_ID         = keyring.get_password('vpn_ID',username)
vpn_secret     = keyring.get_password('vpn_secret',username)


# check whether user have the required passwords stored by keyring
if (cluster_psswrd==None or cluster_ddrss==None):
    print('PERMISSION DENIED.\n\n')
    sys.exit()
else:
    
    if sudo_psswrd==None:
        print('VPN connection cannot be used because root password has not been specified.')
        
        print('\n Attempting SSH connection. \n')
        # ssh job submission
        ssh_copy_command  = 'sshpass -p ' + cluster_psswrd + ' scp ./cluster_commands.sh ' + username + '@' + cluster_ddrss + ':~/cluster_commands.sh'
        os.system(ssh_copy_command)
        ssh_copy_command  = 'sshpass -p ' + cluster_psswrd + ' scp ./xASL_job.run ' + username + '@' + cluster_ddrss + ':~/xASL_job.run'
        os.system(ssh_copy_command)
        ssh_copy_command  = 'sshpass -p ' + cluster_psswrd + ' scp ./xASL_test.m ' + username + '@' + cluster_ddrss + ':~/xASL_test.m'
        os.system(ssh_copy_command)
        ssh_login_command = 'sshpass -p ' + cluster_psswrd + ' ssh ' + username + '@' + cluster_ddrss + ' "bash cluster_commands.sh; exit;"'
        os.system(ssh_login_command)
    else:
        print('\nCreate VPN configuration file \n')
        os.system('echo "' + sudo_psswrd + '" | sudo -S rm -rf ./.my_vpn.conf')
        vpn_file = open('./.my_vpn.conf', "w")
        vpn_file.write('IPSec gateway '+vpn_address+'\n')
        vpn_file.write('IPSec ID '+vpn_ID+'\n')
        vpn_file.write('IPSec secret '+vpn_secret+'\n')
        vpn_file.write('Xauth username '+vpn_username+'\n')
        vpn_file.write('Xauth password '+vpn_password+'\n')
        vpn_file.close()
        
        os.system('echo "' + sudo_psswrd + '" | sudo -S chown root .my_vpn.conf')
        os.system('echo "' + sudo_psswrd + '" | sudo -S chgrp root .my_vpn.conf')
        
        print('Attempting VPN connection. \n')
        # vpn connect
        vpn_login_command = 'echo "'+ sudo_psswrd +'" | sudo -S vpnc-connect ./.my_vpn.conf'
        os.system(vpn_login_command)
        
        os.system('echo "' + sudo_psswrd + '" | sudo -S rm -rf ./.my_vpn.conf')
        print('\nVPN configuration file deleted \n')
        
        print('\nAttempting SSH connection. \n')
        # ssh job submission
        ssh_copy_command  = 'sshpass -p ' + cluster_psswrd + ' scp ./cluster_commands.sh ' + username + '@' + cluster_ddrss + ':~/cluster_commands.sh'
        os.system(ssh_copy_command)
        ssh_copy_command  = 'sshpass -p ' + cluster_psswrd + ' scp ./xASL_job.run ' + username + '@' + cluster_ddrss + ':~/xASL_job.run'
        os.system(ssh_copy_command)
        ssh_copy_command  = 'sshpass -p ' + cluster_psswrd + ' scp ./xASL_test.m ' + username + '@' + cluster_ddrss + ':~/xASL_test.m'
        os.system(ssh_copy_command)
        ssh_login_command = 'sshpass -p ' + cluster_psswrd + ' ssh ' + username + '@' + cluster_ddrss + ' "bash cluster_commands.sh; exit;"'
        os.system(ssh_login_command)
        
        
        print('\nDisconnect VPN\n')
        # vpn disconnet
        os.system('echo "' + sudo_psswrd + '" | sudo -S vpnc-disconnect')

