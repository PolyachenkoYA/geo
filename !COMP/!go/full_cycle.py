#!/usr/bin/python3 -tt
# Copyright 2010 Google Inc.
# Licensed under the Apache License, Version 2.0
# http://www.apache.org/licenses/LICENSE-2.0

# Google's Python Class
# http://code.google.com/edu/languages/google-python-class/

import sys
import numpy as np
import math as mth
import os
import matplotlib.pyplot as plt

all_proc_flags = '-mv'

def find_key(keys, key0):
    for _key in keys:
        if(_key == key0):
            return 1
    
    return 0

def run_it(command):
    print(command)
    return os.system(command) == 0

# Define a main() function that prints a little greeting.
def main():
    args = sys.argv[1:]
    argc_min = 1
    if len(args) < argc_min:
        print('usage:\n./full_proc.py      model_name     [keys]                                                    ')
        sys.exit(1)
        
    model_name = args[0]
    if(len(args) > 1):       
        keys = (all_proc_flags if (args[1] == '-all') else args[1])
        #if(args[1] == '-all'):
        #    keys = all_proc_flags
        #else
        #    keys = args[1]
    else:
        keys = all_proc_flags
    keys = keys.split('-')
        
    #extra_args_str = ''
    #for i in range(1,len(args)):
    #   extra_args_str += (args[i] + ' ')
        
    #if(my.find_key(keys, 'gen')):
    #    if(not my.run_it('./gen ' + model_name)):
    #        return

    if(not run_it('./comp_geo ' + model_name + '            ')):
        return

    file_to_copy = model_name + '.log'
    run_it('mv -v ' + file_to_copy + ' ' + os.path.join('./', model_name, file_to_copy) +  '                                ')
    file_to_copy = model_name + '.res'
    run_it('mv -v ' + file_to_copy + ' ' + os.path.join('./', model_name, file_to_copy) +  '                                ')
#    file_to_copy = model_name + '.prm'
#    run_it('cp -v ' + file_to_copy + ' ' + os.path.join('./', model_name, file_to_copy) +  '                                ')
    file_to_copy = model_name + '.mat'
    run_it('cp -v ' + file_to_copy + ' ' + os.path.join('./', model_name, file_to_copy) +  '                                ')
    file_to_copy = model_name + '.tet'
    run_it('cp -v ' + file_to_copy + ' ' + os.path.join('./', model_name, file_to_copy) +  '                                ')
            
    if(find_key(keys, 'mv')):
        command_to_run = 'cd ../../' + '            '
        print(command_to_run)
        os.chdir('../../')
    
        res_path = os.path.join('./', 'RES', 'DATA')
        run_it('mv ' + os.path.join('./', '!COMP', '!go', model_name) + ' ' + res_path + '            ')
        
        command_to_run = 'cd ' + res_path + '            '
        print(command_to_run)
        os.chdir(res_path)   
        
        #my.run_it('./full_post_proc.py ' + model_name + ' ' + extra_args_str)

# This is the standard boilerplate that calls the main() function.
if __name__ == '__main__':
  main()
