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
    plt.rcParams["font.family"] = 'Times New Roman'
    #plt.rcParams["font.size"] = 15
    #plt.figure(figsize=(5,3))
    #plt.rcParams["figure.dpi"] = 300
    plt.style.use('seaborn')