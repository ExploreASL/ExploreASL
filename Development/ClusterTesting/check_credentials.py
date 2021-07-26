import keyring
import getpass

answer = 'yes'

while answer.lower() == 'yes' or answer.lower() == 'y':
    answer = input('Do you wish to check a stored keyword? [yes/no]:\t\t\t')
    if answer.lower() == 'yes' or answer.lower() == 'y':
        print('Options: [ssh_address/ssh_password/sudo_password/vpn_address/vpn_username/vpn_password/vpn_ID/vpn_secret]\n')
        service_id = input('Select an option:\t\t\t\t\t\t\t')
        username = input("What is the username corresponding to the keyword you wish to check?\t")
        keyword = keyring.get_password(service_id, username)
        if keyword==None: print('Keyword does not exist.')
        else: print( 'The requested keyword is:\t\t\t\t\t\t'+keyword )
    elif answer.lower() != 'no' or answer.lower() == 'n':
        answer = 'yes'
        
import readline
readline.clear_history()

print('Terminal history has been ereased.')