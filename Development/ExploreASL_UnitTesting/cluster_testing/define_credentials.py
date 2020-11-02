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
    
    print('\nPlease edit the vpn_sample.conf file so that it includes the required login details and save it as ".my_vpn.conf".')
    print('Then use the chgrp and chown commands to ensure that .my_vpn.conf is owned by root.')
else:
    print('VPN credentials are not registered.')

# delete password
# keyring.delete_password(service_id, username)
# retrieve password
# keyring.get_password(service_id, username)

