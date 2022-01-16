#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun  1 19:01:47 2021

Personal Archive of generic useful functions

@author: aman
"""

def dictfetchall(cursor):
    "Return all rows from a cursors as a dict"
    columns = [col[0] for col in cursor.description]
    return [
        dict(zip(columns, row))
        for row in cursor.fetchall()
    ]

def pretty_plot(plt):
    "Define Global Plotting Parameters"
    plt.style.use('seaborn')
    plt.rcParams["figure.dpi"] = 600
    plt.rcParams["font.family"] = 'Times New Roman'
    plt.rcParams['xtick.major.size'] = 4
    plt.rcParams['ytick.major.size'] = 4

def ssh(user,host,jump=True,jumphost=None):
    """Pass username and hostname as parameters for ssh"""
    # SSH Client Setup
    import paramiko,os
    pkey_path = os.path.expanduser('~/.ssh/id_rsa') #path to local private key
    port = 22 #default
    jhost = None

    try:
        #Setup a ssh connection using prior login inofrmation

        ssh = paramiko.SSHClient() #SSHClient Object
        k = paramiko.RSAKey.from_private_key_file(pkey_path) #Import/Decode private key
        ssh.set_missing_host_key_policy(paramiko.AutoAddPolicy()) #Add host if unrecognized
        ssh.connect(hostname=host, username=user, pkey=k) #Establish ssh connection

        if jump:
            #Setup channel to jump primary host server
            vmchannel = ssh.get_transport().open_channel("direct-tcpip",(jumphost,port),(host,port))
            #Setup another ssh client to final host
            jhost = paramiko.SSHClient()
            jhost.set_missing_host_key_policy(paramiko.AutoAddPolicy())
            jhost.connect(jumphost, username=user, sock=vmchannel) #Pass channel to socket

    except Exception as e:
        print(str(e))

    return ssh,jhost

def goto(target):
    """Returns a string that mimics goto alias
    that finds/navigates to a given target"""

    string = 'f(){{ cd "$(find /nfs/recons4/CTIOPI/regions -type d -name {} | head -1)";\
            unset -f f;}}; f;'.format(target)

    return string

def postgres():

    import psycopg2
    #Setup for local PostgreSQL server connection
    try:
        db = psycopg2.connect("dbname='postgres'")
    except:
        print("I am unable to connect to the database")

    return db