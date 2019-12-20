import logging

def init():
    global args
    global available_types
    global pack_output
    global logger

    available_types = ['csv', 'json', 'file', 'nii.img', 'nii.hdr']


init()
