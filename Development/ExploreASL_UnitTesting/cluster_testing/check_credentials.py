import keyring
import getpass

answer = 'yes'

while answer.lower() == 'yes' or answer.lower() == 'y':
    answer = input('Do you wish to check any password/address? [yes/no]:\t\t\t')
    if answer.lower() == 'yes' or answer.lower() == 'y':
        service_id = input('[ssh_address/ssh_password/sudo_password]:\t\t\t\t')
        username = input("What is the username corresponding to the password you wish to check?\t")
        print( keyring.get_password(service_id, username) )
        
import readline
readline.clear_history()

print('Terminal history has been ereased.')