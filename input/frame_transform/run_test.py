import glob
import os.path
import shutil
import subprocess
import sys

import h5py
import numpy as np
import yt

yt.funcs.mylog.setLevel(40)

def _load_block_list(block_list):
    out = []
    with open(block_list, 'r') as f:
        for line in f:
            out.append(line.strip().split(' '))
    return out

def _load_block_prop(fname,block_name):
    with h5py.File(fname, 'r') as f:
        assert block_name in f.keys()
        block = f[block_name]
        
        frame_velocity = block.attrs.get(
            'frame_velocity',
            np.array([np.nan, np.nan, np.nan]))
        origin_offset = block.attrs.get(
            'last_updated_origin_offset',
            np.array([np.nan, np.nan, np.nan]))
        last_updated_origin_offset = block.attrs.get(
            'last_updated_origin_offset',
            np.array([np.nan, np.nan, np.nan]))
        last_frame_update_time = block.attrs.get(
            'last_frame_update_time', [np.nan])[0]
        return (frame_velocity, last_updated_origin_offset,
                last_frame_update_time)

def load_block_frame_prop(blocklist,load_all = False):
    if load_all:
        slc = slice(None)
    else:
        slc = slice(0,1)

    vel_l = []
    last_updated_origin_offset_l = []
    last_frame_update_time_l = []

    block_basename_pairs = _load_block_list(blocklist)
    dirname = os.path.dirname(blocklist)

    for block_name, basename in block_basename_pairs[slc]:
        fname = os.path.join(dirname, basename)
        tmp = _load_block_prop(fname,block_name)
        vel_l.append(tmp[0])
        last_updated_origin_offset_l.append(tmp[1])
        last_frame_update_time_l.append(tmp[2])
    if load_all:
        return (np.array(vel_l),
                np.array(last_updated_origin_offset_l),
                np.array(last_frame_update_time_l))
    else:
        return (vel_l[0], 
                last_updated_origin_offset_l[0],
                last_frame_update_time_l[0])
    
def load_block_frame_series(blocklist_l, 
                            load_all_blocks = False):
    out_vel_l = []
    out_last_updated_offset_l = []
    out_last_update_time_l = []
    for blocklist in blocklist_l:
        tmp = load_block_frame_prop(blocklist,
                                    load_all_blocks)
        out_vel_l.append(tmp[0])
        out_last_updated_offset_l.append(tmp[1])
        out_last_update_time_l.append(tmp[2])
    return (np.array(out_vel_l),
            np.array(out_last_updated_offset_l),
            np.array(out_last_update_time_l))

def check_equally_spaced_times(name, snaps,times,t_spacing, start_t = 0.,
                               rtol=1.e-14, atol=1e-12):
    # utility function used when the timestep is expected to be constant
    # returns true when the times are consistent with expectations and False otherwise
    assert sorted(snaps) == snaps
    
    expected_times = []
    cur_t = start_t
    cur_index = 0
    for cur_snap in range(0, snaps[-1] + 1):
        if snaps[cur_index] == cur_snap:
            expected_times.append(cur_t)
            cur_index+=1
        cur_t += t_spacing
    assert len(expected_times) == len(snaps)
    
    return np.allclose(times, expected_times, rtol = rtol, atol = atol)


class TrackerParam:
    def __init__(self, field, cut_region_str):
        self.field = field
        self.cut_region_str = cut_region_str
    def _get_data(self,ds):
        if self.cut_region_str is None:
            return ds.all_data()
        else:
            return ds.all_data().cut_region([self.cut_region_str])

    def calc_average_velocity(self, ds, axis):
        # in reality, we should compute the values for 
        # each block. Order the blocks and then sum the values
        assert axis >=0 and axis <3
        vname = 'velocity_' + 'xyz'[axis]
        data = self._get_data(ds)
        
        total = data[self.field].sum()
        if total.v == 0:
            return None
        weighted_sum = (data[vname]*data[self.field]).sum()
        return (weighted_sum/total).to('code_velocity').v
    
    def calc_min_velocitiy(self,ds,axis):
        assert axis >=0 and axis <3
        vname = 'velocity_' + 'xyz'[axis]
        data = self._get_data(ds)
        return (data[vname].min().to('code_velocity').v)


def get_velocity_prop(fnames, tracker_param, axis = None):

    bulk_velocities = []
    min_velocities = []
    times = []
    for fname in fnames:
        ds = yt.load(fname)
        times.append(ds.current_time.to('code_time').v)
        for l, method in [(bulk_velocities,'calc_average_velocity'),
                          (min_velocities,'calc_min_velocitiy')]:
            if axis is None:
                ax_l = [0,1,2]
            else:
                ax_l = [axis]
            cur_vel = []
            for ax in ax_l:
                bulk_velocity_comp = getattr(tracker_param,method)(ds,ax)
                if bulk_velocity_comp is None:
                    cur_vel.append(np.nan)
                else:
                    cur_vel.append(bulk_velocity_comp)
            if len(cur_vel) == 1:
                l.append(cur_vel[0])
            else:
                l.append(cur_vel)
    return np.array(times), np.array(bulk_velocities), np.array(min_velocities)



def check_bulk_velocity_residual(sim_name, axis_str, snaps, times, 
                                 bulk_velocities, atol):
    # returns true if it meets expectations. Otherwise, it returns False
    # snaps, times, bulk_velocity should only include values when the
    # velocity is updated.
    if isinstance(atol,(list,tuple)):
        assert len(atol) == len(snaps)
        atol_l = atol
    else:
        atol_l = [atol for i in snaps]
        
    for snap, time, bulk_velocity, cur_atol in zip(snaps, times, 
                                                   bulk_velocities, atol_l):
        if np.isnan(cur_atol) or np.isnan(bulk_velocity):
            continue
        elif np.abs(bulk_velocity) > np.abs(cur_atol):
            _msg = ('FAILED: The residual velocity-{axis_str} for snapshot '
                    '{snap} (t = {time}) of the {sim_name} simulation is '
                    '{residual}. It is expected to be smaller than {cur_atol}')
            print(_msg.format(axis_str = axis_str, snap = snap, time = time,
                              sim_name = sim_name, residual=bulk_velocity,
                              cur_atol = cur_atol))
            return False
    return True

def load_relevant_frame_transform_data(sim_dir_pattern, axis, 
                                       density_thresh=16.):
    tracker_param = TrackerParam(
        'density','obj["density"] >= {}'.format(float(density_thresh)))
    sim_dirs = sorted(glob.glob(sim_dir_pattern))
    fnames = [_prep_sim_fname(d) for d in sim_dirs]
    
    snaps = [int(e.split('_')[-1]) for e in sim_dirs]
    vels,last_updated_offsets, last_update_times \
    = load_block_frame_series(fnames, load_all_blocks = False)
    times, bulk_velocities, min_velocities \
        = get_velocity_prop(fnames, tracker_param, axis = axis)
    return (fnames, snaps, times, vels, last_updated_offsets, last_update_times, 
            bulk_velocities, min_velocities)


def _prep_sim_fname(sim_dir):
    return os.path.join(
        sim_dir, os.path.basename(sim_dir) + '.block_list')

def _get_update_cycles(last_update_times,
                       cycles, times):
    if (np.diff(cycles) != 1.).any():
        raise ValueError(
            'The cycle numbers, {!r}, are not contiguous'.\
            format(cycles))
    
    # whenever a frame update occurs the last_update_time is set
    # to the time for the following cycle
    update_occurred = (last_update_times == times)[1:]
    # in other words, we know that updates occurred during 
    # during the cycle snaps[:-1][update_occurred].
    update_cycles = np.array(cycles)[:-1][update_occurred]
    # we won't be able to check if an update occured in 
    # update_cycles.max()
    # note that we will only see the updates in snapshot output at
    # the cycle immediately following the update (update_cycles+1)
    return update_cycles

def _check_update_cycles(expected_update_cycles, update_cycles):
    expect = np.sort(expected_update_cycles)
    actual = np.sort(update_cycles)
    
    if (actual == expect).all():
        return True
    else:
        _msg = ('FAILED: expected updates during cycles {!r}. '
                'Updates actually occurred during cycles {!r}')
        print(_msg.format(expect,actual))
        return False
    
    
def check_origin_offset(update_cycles,times, frame_velocities,
                        last_updated_offset, atol=0., rtol= 1.e-14):
    expected_last_origin_offset = np.array([0., 0., 0.])
    last_update_time = np.finfo('float64').min
    current_frame_velocity = np.array([0., 0., 0.])
    for cycle in range(len(times)):
        if not np.allclose(last_updated_offset[cycle],
                           expected_last_origin_offset,
                           atol = atol, rtol = rtol):
            _msg = ('FAILED: During cycle {:d}, the expected last '
                    'origin offset is {!r}. The actual '
                    'last_updated_origin_offset is {!r}')
            print(_msg.format(cycle,expected_last_origin_offset,
                              last_updated_offset[cycle]))
            return False
        if cycle in update_cycles:
            assert cycle+1 < len(times)
            next_cycle_time = times[cycle+1]
            if last_update_time != np.finfo('float64').min:
                dt = next_cycle_time - last_update_time
                expected_last_origin_offset += current_frame_velocity * dt
            current_frame_velocity = frame_velocities[cycle+1]
            last_update_time = next_cycle_time
    return True


def check_velocity_properties(update_type, update_cycles, frame_velocities,
                              bulk_velocities, min_velocities,
                              atol = 1.e-16):
    update_types = ['weighted_average', 'min', 'min_zero_floor']
    assert update_type in update_types
    
    if update_type in update_types[1:]:
        w = (np.array(update_cycles)+1,)
        actual_mins = (min_velocities + frame_velocities[:,0])[w]
        if not ((actual_mins<0.).any() and (actual_mins>0.).any()):
            _msg = ('The update_type: "{}" can\'t be fully tested because'
                    'there needs to be at least one update cycle when '
                    'the min velocity negative and another where its positive')
            raise ValueError(_msg.format(update_type))
            
    orig_frame_velocity = 0.

    for cycle in range(len(frame_velocities)):
        if cycle not in update_cycles:
            continue
        # remember the effects of the frame transform is seen
        # in the snapshot immediately after the update (cycle+1)
        if update_type == 'weighted_average':
            residual = np.abs(bulk_velocities[cycle + 1])
            if residual > atol:
                _msg = ('FAILED: The residual weighted velocity for snapshot '
                        'after the frame transform in cycle {cycle} of the is '
                        '{residual}. It is expected to be smaller than {atol}')
                print(_msg.format(cycle = cycle,residual = residual, atol = atol))
                return False
        elif update_type == 'min':
            if (min_velocities[cycle+1] != 0.):
                _msg = ('FAILED: The minimum velocity for the snapshot saved '
                        'after the frame transform in cycle {cycle} is {min_v}.'
                        'It is expected to be 0.')
                print(_msg.format(cycle = cycle,
                                  min_v = repr(min_velocities[cycle+1])))
                return False
        elif update_type == 'min_zero_floor':
            if (min_velocities[cycle+1] > 0.):
                _msg = ('FAILED: The minimum velocity for the snapshot saved '
                        'after the frame transform in cycle {cycle} is {min_v}.'
                        'It is expected to be 0 or negative.')
                print(_msg.format(cycle = cycle,
                                  min_v = repr(min_velocities[cycle+1])))
                return False
    return True

def check_frame_transform_sim(template, expected_update_cycles, 
                              update_type = 'weighted_average',
                              density_thresh = 16.):
    """
    This should probably be refactored.
    """

    fnames, snaps, times, vels, last_updated_offsets, last_update_times, \
    bulk_velocities, min_velocities \
        = load_relevant_frame_transform_data(template, 0,
                                             density_thresh = density_thresh)
    assert snaps[0] == 0
    if vels[0][0] != 0.:
        print('FAILED: The frame velocity at snapshot 0 should be 0')
        return False
    if np.amax(expected_update_cycles) >= snaps[-1]:
        raise ValueError(
            'Because we can only determine whether a frame transform '
            'occurs in a given cycle if we have the snapshot from the '
            'next cycle, the max expected_update_cycle cannot exceed '
            'the maximum snapshot number ({})'.format(snaps[-1]))

    update_cycles = _get_update_cycles(
        last_update_times = last_update_times,
        cycles = snaps, times = times)
    if not _check_update_cycles(expected_update_cycles, 
                                update_cycles):
        return False
    if not check_origin_offset(update_cycles, times, vels,
                               last_updated_offsets, atol=0., rtol= 1.e-14):
        return False
    # finally, let's check the bulk velocity properties
    if not check_velocity_properties(update_type, update_cycles, vels,
                                     bulk_velocities, min_velocities,
                                     atol = 1.e-16):
        return False

    return True


def run_tests(executable):
    
    temp = 'input/frame_transform/{}.in'
    
    for elem in ['update_every_cycle',
                 'update_every_third_cycle',
                 'update_every_third_cycle_by_time',
                 'update_every_1.5_cycles_by_time',
                 'min-update_every_cycle',
                 'min_zero_floor-every_cycle']:
        input_file = temp.format(elem)
        assert os.path.isfile(input_file)
        subprocess.call(executable.strip() + ' ' + input_file,
                        shell=True)

def analyze_tests():
    # Since we are only running a simulation for an incredibly short period of 
    # time, and fixing the timestep per cycle, a more robust test would be to
    # run a simulation without frame tracking and then predict all of the 
    # simulation properties once frame tracking is included. Doing this would
    # also enable us to ensure that inflow conditions are updated properly

    # For simplicity, the timestep for each cycle is fixed to 1/512 (so that 
    # we don't need to worry about floating point errors). In effect the time
    # at the start of the cycle (when output is written to disk) is simply 
    # (cycle number)/512.
    r = []

    # the first 4 tests test out frame updates with a incremental velocity 
    # changes set by the weighted average. The difference is the schedule 
    # type that they use. We could probably replace 3 of them with unit tests
    r.append(check_frame_transform_sim(
        'frame_transform-every_cycle_*', list(range(16)), 
        update_type = 'weighted_average', density_thresh = 16.))

    # arbitrarily, we decided to make the following tests start at cycle 1
    r.append(check_frame_transform_sim(
        'frame_transform-every_third_cycle_[0-9]*', [1,4,7,10,13], 
        update_type = 'weighted_average', density_thresh = 16.))
    r.append(check_frame_transform_sim(
        'frame_transform-every_third_cycle_time_*', [1,4,7,10,13],
        update_type = 'weighted_average', density_thresh = 16.))

    # now for the test of updates every 1.5 cycles (using a "minimum_time"
    # based schedule)
    # the scheduled minimum times are:
    #     t = [0, 1.5, 3, 4.5, 6, 7.5, 9, 10.5, 12, 13.5, 15, 16.5, ...]/512
    # The first simulation times that occur after each of the above are:
    #     t = [0,  2,  3,   5, 6,   8, 9,   11, 12,   14, 15,   17, ...]/512
    # The cycle number during which each update occurs is 512 times each 
    # of the above. We drop all cycle numbers >= 12 since cycle 12 is the 
    # last time we save a snapshot and we can only tell if a frame 
    # transformation occured in a given cycle if we have the snapshot from
    # the following cycle.
    three_half_update_cycles = [0, 2, 3, 5, 6, 8, 9, 11, 12, 14, 15]
    r.append(check_frame_transform_sim(
        'frame_transform-every_1.5_cycles_*', three_half_update_cycles,
        update_type = 'weighted_average', density_thresh = 16.))

    r.append(check_frame_transform_sim(
        'frame_transform-min-every_cycle_*', list(range(16)),
        update_type = 'min', density_thresh = 15.999999))

    r.append(check_frame_transform_sim(
        'frame_transform-min_zero_floor-every_cycle_*', list(range(16)),
        update_type = 'min_zero_floor', density_thresh = 15.999999))

    print('{} tests passed out of {}'.format(np.sum(r), len(r)))
    if np.sum(r) == len(r):
        return True
    return False


def cleanup():
    dir_templates = ['frame_transform-every_cycle_*',
                     'frame_transform-every_third_cycle_[0-9]*',
                     'frame_transform-every_third_cycle_time_*',
                     'frame_transform-every_1.5_cycles_*',
                     'frame_transform-min-every_cycle_*',
                     'frame_transform-min_zero_floor-every_cycle_*']
    for template in dir_templates:
        for dir_name in glob.glob(template):
            print(dir_name)
            if os.path.isdir(dir_name):
                shutil.rmtree(dir_name)
if __name__ == '__main__':
    run_tests('bin/enzo-p')
    tests_passed = analyze_tests()
    cleanup()

    if tests_passed:
        sys.exit(0)
    else:
        sys.exit(3)
