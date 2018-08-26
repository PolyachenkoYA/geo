#!/usr/bin/python3 -tt
# Copyright 2010 Google Inc.
# Licensed under the Apache License, Version 2.0
# http://www.apache.org/licenses/LICENSE-2.0

# Google's Python Class
# http://code.google.com/edu/languages/google-python-class/

import sys
import numpy as np
import math
import os
import matplotlib.pyplot as plt
import copy as cp

#import mylib_molecules as my

def load_file(file_name):
    return (open(file_name, 'rU').read()).split('\n')

def read_params(filename):
    params = load_file(filename)
    params = params[0:2]
    kys = params[0].split()
    vals = params[1].split()
    
    if(len(kys) != len(vals)):
        print('wrong params file:\n' + param_path)
        sys.exit(1)
    
    params = {}
    for i in range(len(kys)):
        params[kys[i]] = float(vals[i])
        
    params['Nrays'] = int(params['Nrays'])
        
    return params    
    
def read_res(filename, params):
    data = load_file(filename)
    head_data = data[0].split()
    Ndet = int(head_data[0])
    Nbars = int(head_data[1])
    data = data[1:]
    
    res_data = []
    x = np.zeros(Nbars)
    y = np.zeros(Nbars)
    for i_det in range(Ndet):
        for i in range(Nbars):
            xy_data = data[i_det*Nbars + i].split()
            x[i] = float(xy_data[0])
            y[i] = float(xy_data[1])
        res_data.append(cp.deepcopy([x,y])) # res_data[i_det][0\1 <=> x\y][i]
        
    return Ndet, Nbars, res_data

# Define a main() function that prints a little greeting.
def main():
    args = sys.argv[1:]
    argc_min = 1
    if len(args) < argc_min:
        print('usage: ./draw_res.py    model_name    [keys]')
        sys.exit(1)
        
    model_name = args[0]
    params = read_params(model_name + '.prm')
    Ndet, Nbars, data = read_res(model_name + '.res', params)
    
    fig_c = -1
    fig = []
    draw_on_screen = 1
    #draw_on_screen = my.find_key(keys, 'keep')
    if(not draw_on_screen):
        plt.switch_backend('Agg') 
            
    for i_det in range(Ndet):
        fig_c += 1
        fig.append(plt.figure(fig_c))
        x = data[i_det][0][:]
        y = data[i_det][1][:]
        plt.plot(x, y)
        plt.xlabel('time')
        plt.ylabel('A')
        plt.grid(True)
    
        if(draw_on_screen):
            fig[fig_c].show()
        #path = os.path.join(graph_dir, 'Energy_' + time_gaps_str + '.png')
        #fig[fig_c].savefig(path)    
        
    if(draw_on_screen):
        input()
    
        
# This is the standard boilerplate that calls the main() function.
if __name__ == '__main__':
  main()
