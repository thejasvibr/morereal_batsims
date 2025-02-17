'''
From scratch bat-realistic simulated trajectories
=================================================
'''
import argparse
import matplotlib.pyplot as plt
import pyroomacoustics as pra
import numpy as np 
import scipy.spatial as spl
import scipy.signal as signal 
import scipy.interpolate as si
choose = np.random.choice

import pandas as pd
#%%
parse_roomdims = lambda X: np.fromstring(X[1:-1], dtype=np.float64, sep=',')
#%%

args = argparse.ArgumentParser()
args.add_argument('-seed', type=int)
args.add_argument('-room_dims', type=parse_roomdims, default=[7,5,3], help='LxBxH of the room')
args.add_argument('-ipi', type=float, default=100, help='interpulse interval in milliseconds')
args.add_argument('-firstcall_delay', type=float, default=0,
                  help='Time at which the first call is emitted in seconds')
param = args.parse_args()

np.random.seed(param.seed)
#%%
array_geom = np.array(([3, 4.9, 1.5],
                      [2.5, 4.9, 1],
                      [2, 4.9, 1.5],
                      [1.5, 4.9,1],
                      [1.0, 4.9, 1.5],
                      [0.01, 2.5, 1.5],
                      [0.01, 3.25, 1.5],
                      [0.01, 4, 2.0],
                      ))

ipi = param.ipi
room_dims = param.room_dims # l,b, h

fs = 192000
ray_tracing=False
rt60_tgt = 0.2  # seconds
e_absorption, ref_order = pra.inverse_sabine(rt60_tgt, room_dims)
ref_order = 1

room = pra.ShoeBox(
    room_dims, fs=fs, materials=pra.Material(e_absorption),
    max_order=ref_order
    )

#%%
# set duration of bat flight and other parameters


duration_range = 1.5
z_min, z_max = np.percentile(array_geom[:,2], [0,100])
xmin, xmax = room.get_bbox()[0,:]
ymin, ymax = room.get_bbox()[1,:]
x_range = np.arange(xmin+0.01, xmax-0.01, 0.01)
y_range = np.arange(ymin+0.01, ymax-0.01, 0.01)
zminmax_range = np.arange(z_min, z_max, 0.1)

n_interpolation_pts = 10
vmin, vmax = 3, 6
n_pts = 4
#%%
bad_speedprofile = True
jj = 0
while bad_speedprofile:
    jj += 1 
    if np.remainder(jj,1000)==0:
        print(f'In the {jj}th iteration...still trying to find a good random traj')
    # choose 4 z points that are within the min-max of the microphone values
    z_coods = np.sort(choose(zminmax_range, n_pts) )
    # choose 4 x-y points that are within the 2D outline of the room
    # super-inefficient, but it works in a loop perhaps
    x_coods = choose(x_range, n_pts)
    y_coods = choose(y_range, n_pts)
    time_at_points = np.linspace(0,1,n_pts)
    
    
    time_span = np.linspace(0,1,n_interpolation_pts)
    seed_xyz = np.array((x_coods, y_coods, z_coods)).T
    spline = [si.interp1d(time_at_points, seed_xyz[:,i], 'cubic') for i in range(3)]
    interp_spline = np.zeros((time_span.size, 3))
    
    for i in range(3):
        interp_spline[:,i] = spline[i](time_span)
    
    # check that the points when connected by a 3D spline don't have an average
    # speed that is 3-5 m/s
    distmat = spl.distance_matrix(interp_spline, interp_spline)
    interpoint_distances = []
    for i,j in zip(range(1,distmat.shape[0]), range(distmat.shape[0]-1)):
        interpoint_distances.append(distmat[i,j])
    speed_profile = np.array(interpoint_distances)/(1/n_interpolation_pts)

    # if yes, then TICK  
    speed_in_limits =  np.all(np.logical_and(speed_profile>=vmin, speed_profile<=vmax))
    if speed_in_limits:
        t_highresolution = np.linspace(0,1,1000)
        traj_highresolution = np.zeros((t_highresolution.size,3))
        for i in range(3):
            traj_highresolution[:,i] = spline[i](t_highresolution)
        xyz_inroom = [np.all(np.logical_and(traj_highresolution[:,i]>=0, traj_highresolution[:,i]<=room_dims[i])) for i in range(3)]
    
    
        if np.all(xyz_inroom):
            bad_speedprofile = False
                              
                          



#%%
# Now choose the emission points to be spaced evenly according to ipi set
t_xyz = np.column_stack((t_highresolution, traj_highresolution))
emision_inds = range(t_highresolution.size)[::int(ipi)]
t_xyz_emission = np.column_stack((t_xyz, np.tile(0,t_xyz.shape[0])))
t_xyz_emission[emision_inds,-1] = 1 

only_emission_points = t_xyz_emission[t_xyz_emission[:,-1]==1,:]
emission_points = pd.DataFrame(data={'t': only_emission_points[:,0],
                                     'x':only_emission_points[:,1],
                                     'y': only_emission_points[:,2],
                                     'z': only_emission_points[:,3],
                                     })

# design the synthetic bat call
call_durn = float(choose(np.linspace(5e-3,7e-3,5),1))
minf, maxf = float(choose([15000,20000,22000],1)), float(choose([88000, 89000, 92000],1))
t_call = np.linspace(0,call_durn, int(fs*call_durn))
call_type = str(choose(['logarithmic','linear','hyperbolic'], 1)[0])
batcall = signal.chirp(t_call, maxf, t_call[-1], minf,call_type)
batcall *= signal.tukey(batcall.size, 0.1)
batcall *= 0.1
#batcall = np.random.normal(0,0.1, batcall.size)

# Now filter out only the emission points and generate the simulated audio

for rownum, row in emission_points.iterrows():
    t,x,y,z = row
    t = t+ param.firstcall_delay
    call_position = np.array([x, y, z])
    if room.is_inside(call_position):
        room.add_source(position=call_position.tolist(),
                        signal=batcall, delay=t)
    else:
        raise ValueError(f'Position {call_position} is not inside')

room.add_microphone_array(array_geom.T)
print('...computing RIR...')
room.compute_rir()
print('room simulation started...')
room.simulate()


print('room simulation ended...')
sim_audio = room.mic_array.signals.T

#%%
    
#%%
plt.figure()
plt.subplot(111, projection='3d')
plt.plot(emission_points['x'], emission_points['y'], emission_points['z'],'*')
plt.plot(array_geom[:,0], array_geom[:,1], array_geom[:,2],'*')
plt.gca().set_xlim(0,room_dims[0])
plt.gca().set_ylim(0,room_dims[1])
plt.gca().set_zlim(0,room_dims[2])


plt.savefig(f'trajectory_mics_seed-{param.seed}.png')

#%%
dd_empoints = spl.distance_matrix(emission_points.loc[:,'x':'z'], emission_points.loc[:,'x':'z'])
distance_travelled = [dd_empoints[i,j] for i,j in zip(range(dd_empoints.shape[0]), range(1,dd_empoints.shape[0]))]
speed = np.array(distance_travelled)/(param.ipi*1e-3)


#%% SAving all simulation related parameters
room.mic_array.to_wav(
    f"single_bat_{param.seed}.wav",
    norm=True,
    bitdepth=np.int16,
)

# save the array geometry
pd.DataFrame(array_geom, columns=['x','y','z']).to_csv(f'single_bat_{param.seed}.csv')

# save the flight path with the emission points and times 
pd.DataFrame(t_xyz_emission, index=range(t_xyz_emission.shape[0]), columns='t,x,y,z,is_emissionpoint'.split(',')).to_csv(f'flight_and_call_groundtruth_{param.seed}.csv')


