# from __future__ import print_function
import numpy as np
import scipy
import os
import sys
import glob
import time
import datetime as dt
# import progressbar as pb
import argparse
import subprocess as sp
import multiprocessing

def process_pair(value):
    process_time=time.time()
    sp.check_call(value,stderr=sp.STDOUT,shell=True)
    runtime=time.time()-process_time
    return runtime

if __name__ == '__main__':

    parser = argparse.ArgumentParser( \
        description="""run commands in command_file on num_simultaneous_processes processors""",
        epilog='>>  <<',
        formatter_class=argparse.RawDescriptionHelpFormatter)


    parser.add_argument('-num_simultaneous_processes', 
                        action='store', 
                        type=int, 
                        default=5,
                        help='number of simultaneous processes (cores) to use [%(default)s]')
    parser.add_argument('pycorr_command_file', 
                        action='store', 
                        type=str, 
                        default='pycorr_command_file.txt',
                        help='text file full of pycorr commands to run [%(default)s]')
    args = parser.parse_args()


    with open(args.pycorr_command_file,'r') as pyc_cf:
        process_list = pyc_cf.read().splitlines()

    print('found %d commands to run'%(len(process_list)))

    start_time=time.time()

    num_simultaneous_processes=args.num_simultaneous_processes

    pool = multiprocessing.Pool(num_simultaneous_processes)
    results = []
    r = pool.map_async(process_pair, process_list, callback=results.append)
    r.wait() # Wait on the results
    print(results)


    print('%f seconds for %d processes - %f seconds/process'%(time.time()-start_time,len(process_list),(time.time()-start_time)/len(process_list)))