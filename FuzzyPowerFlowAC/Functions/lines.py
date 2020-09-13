# -*- coding: utf-8 -*-
"""
Created on Sat Sep 12 12:48:15 2020

@author: Sofia
"""
import numpy as np

def R(lines, bus , data_lines):
    r1_lines = np.zeros([bus,bus])
    for i in range(lines):
        x = int(data_lines[i,0])-1
        y = int(data_lines[i,1])-1
        r1_lines[x,y] = data_lines[i,2]
        r1_lines[y,x] = r1_lines[x,y]
    return r1_lines

def X(lines, bus , data_lines):
    x1_lines = np.zeros([bus,bus])
    for i in range(lines):
        x = int(data_lines[i,0])-1
        y = int(data_lines[i,1])-1
        x1_lines[x,y] = data_lines[i,3]
        x1_lines[y,x] = x1_lines[x,y]
    return x1_lines

def Z(lines, bus ,data_lines):
    z1_lines=np.zeros([bus,bus],dtype = complex)
    for i in range(lines):
        x = int(data_lines[i,0])-1
        y = int(data_lines[i,1])-1
        a = data_lines[i,2]
        b = data_lines[i,3]
        z1_lines[x,y] = complex(a,b)
        z1_lines[y,x] = z1_lines[x,y]
    return z1_lines

def Y(lines,bus , z_lines, data_lines):
    y1_lines = np.zeros([bus,bus],dtype = complex)
    for i in range(lines):
        x = int(data_lines[i,0])-1
        y = int(data_lines[i,1])-1
        y1_lines[x,y] = 1/(z_lines[x,y])
        y1_lines[y,x] = 1/(z_lines[y,x])
        
    return y1_lines
    
