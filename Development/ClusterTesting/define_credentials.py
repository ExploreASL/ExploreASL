import keyring
import getpass

username = input('Define username for password protection:\t')
cluster_address = input('Define cluster SSH address:\t\t\t')
cluster_password = getpass.getpass('Define cluster SSH password:\t\t\t')

myFile = open(".login_srnm", "w")
myFile.write(username)
myFile.close()

service_id = 'ssh_address'
keyring.set_password(service_id, username, cluster_address)
service_id = 'ssh_password'
keyring.set_password(service_id, username, cluster_password)

print('\nIs VPN access required to log in to the defined cluster?')
answer = input("Yes or No [Y/N]:\t\t\t\t")

if answer.lower() =='yes' or answer.lower() =='y':
    print('Additional credentials are required for VPN connection.')
    
    sudo_password = getpass.getpass('Provide sudo password for local machine:\t')
    service_id = 'sudo_password'
    keyring.set_password(service_id, username, sudo_password)
    
    print('\nThe following information regarding the VPN connection will be temporarily stored in ".my_vpn.conf" until the job submission is complete.')
    
    vpn_address = input('Define VPN address:\t\t\t\t')
    vpn_username = input('Define VPN username:\t\t\t\t')
    vpn_password = getpass.getpass('Define VPN password:\t\t\t\t')
    vpn_ID = getpass.getpass('Define VPN identification:\t\t\t')
    vpn_secret = getpass.getpass('Define VPN secret word:\t\t\t\t')

    
    service_id = 'vpn_address'
    keyring.set_password(service_id, username, vpn_address)
    service_id = 'vpn_username'
    keyring.set_password(service_id, username, vpn_username)
    service_id = 'vpn_password'
    keyring.set_password(service_id, username, vpn_password)
    service_id = 'vpn_ID'
    keyring.set_password(service_id, username, vpn_ID)
    service_id = 'vpn_secret'
    keyring.set_password(service_id, username, vpn_secret)

#    print('Then use the chgrp and chown commands to ensure that .my_vpn.conf is owned by root.')
else:
    print('VPN credentials are not registered.')

# delete password
# keyring.delete_password(service_id, username)
# retrieve password
# keyring.get_password(service_id, username)

