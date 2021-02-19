#!/bin/python

# runs a test to check the floor of VLCT 


import os.path
import sys
import shutil
import subprocess

import numpy as np
import yt

yt.mylog.setLevel(30) # set yt log level to "WARNING"

from testing_utils import prep_cur_dir, EnzoEWrapper

def run_tests(executable):

    temp = 'input/vlct/passive_scalar_floor/{}'
    wrapper = EnzoEWrapper(executable,temp)

    wrapper('passive_scalar_floor_test.in')


def check_passive_scalar_floor(cycle0_output, cycle1_output,
                               passive_scalar_density_floor):

    for fname, floor_application in [(cycle0_output, False),
                                     (cycle1_output, True)]:
        ds = yt.load(fname)
        ad = ds.all_data()
        # these are hard-coded expectations:
        for field, val in [('density', 1.0), ('velocity_x', 1.0),
                           ('velocity_y', 0.0), ('velocity_z', 0.0),
                           ('total_energy',  4.0)]:
            # I'm not quite sure I understand why the simulation values don't
            # exactly match in the output values
            if val == 0.0:
                diff = np.abs((ad[field].in_cgs().v))
            else:
                diff = np.abs((ad[field].in_cgs().v - val)/val)
            
            if (diff > 5e-16).any():
                msg = ("FAILED: the \"{}\" field is expected to have a "
                       "uniform value of {} (in cgs units). This is not the "
                       "case for {}.")
                print(msg.format(field,val,fname))
                return False


        if floor_application:
            if (ad['red'].v != passive_scalar_density_floor).any():
                msg = ("FAILED: the \"red\" passive scalar field is expected "
                       "to have a uniform value of {} g/cm^3 after the floor "
                       "for passive scalars is applied (in {})")
                print(msg.format(passive_scalar_density_floor, fname))
                return False
        else:
            if (ad['red'].v != 0.0).any():
                msg = ("FAILED: the \"red\" passive scalar field is expected "
                       "to have a uniform value of 0.0 g/cm^3 before the "
                       "floor for passive scalars is applied (in {})")
                print(msg.format(fname))
                return False
    
    # check that density, velocity_x, velocity_y, velocity_z, and total energy
    # are all unchanged
    return True

def analyze_tests():
    r = []
    r.append(check_passive_scalar_floor(
        'passive_floor_0/passive_floor_0.block_list',
        'passive_floor_1/passive_floor_1.block_list',
        passive_scalar_density_floor = 0.25
    ))
    n_passed = np.sum(r)
    n_tests = len(r)
    print("{:d} Tests passed out of {:d} Tests.".format(n_passed,n_tests))

    return n_passed == n_tests

def cleanup():
    dir_names = ['passive_floor_0',
                 'passive_floor_1']
    for dir_name in dir_names:
        if os.path.isdir(dir_name):
            shutil.rmtree(dir_name)

if __name__ == '__main__':
    executable = 'bin/enzo-p'

    # this script can either be called from the base repository or from
    # the subdirectory: input/vlct
    prep_cur_dir(executable)

    # run the tests
    tests_complete = run_tests(executable)

    # analyze the tests
    tests_passed = analyze_tests()

    # cleanup the tests
    cleanup()

    if tests_passed:
        sys.exit(0)
    else:
        sys.exit(3)
