import keyring
import getpass
import os

answer = 'yes'

while answer.lower() == 'yes' or answer.lower() == 'y':
    answer = input('Do you wish to delete any password/address? [yes/no/all]:\t\t\t')
    if answer.lower() == 'yes' or answer.lower() == 'y':
        print('Options:\n[ssh_address/ssh_password/sudo_password/vpn_address/vpn_username/vpn_password/vpn_ID/vpn_secret]\n')
        service_id = input('Select an option:\t\t\t\t')
        username = input("What is the username corresponding to the password you wish to delete?\t")
        keyring.delete_password(service_id, username)
    if answer.lower() == 'all':
        username = input("What is the username corresponding to the passwords you wish to delete?\t")
        for i in ['ssh_address','ssh_password','sudo_password','vpn_address','vpn_username','vpn_password','vpn_ID','vpn_secret']:
            try:
                keyring.delete_password(i, username)
            except:
                print(i+' probably has been already removed.')
        os.system("bash clear.sh")