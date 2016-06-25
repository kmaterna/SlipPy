#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt
import slippy.patch
import slippy.basis
import slippy.io
import slippy.bm
import slippy.inversion
import matplotlib.pyplot as plt
import logging
import modest
logging.basicConfig(level=logging.INFO)
np.random.seed(1)

params = {'strike':70.0,
          'dip':45.0,
          'length':200000.0,
          'width':60000.0,
          'position':[-84.2,43.3,0.0],
          'Nlength':40,
          'Nwidth':20,
          'basis':[[1.0,1.0,0.0],[1.0,-1.0,0.0]],
          'penalty':0.01}

slippy.inversion.main(params,
                      gps_input_file='synthetic_gps.txt',
                      gps_output_file='predicted_gps.txt',
                      slip_output_file='predicted_slip.txt')
modest.summary()                   

pred_pos_geo,pred_disp,pred_sigma = slippy.io.read_gps_data('predicted_gps.txt')
bm = slippy.bm.create_default_basemap(pred_pos_geo[:,0],pred_pos_geo[:,1])
pred_pos_cart = slippy.bm.geodetic_to_cartesian(pred_pos_geo,bm)

obs_pos_geo,obs_disp,obs_sigma = slippy.io.read_gps_data('synthetic_gps.txt')
obs_pos_cart = slippy.bm.geodetic_to_cartesian(obs_pos_geo,bm)

input = slippy.io.read_slip_data('predicted_slip.txt')
patch_pos_geo = input[0]
patch_pos_cart = slippy.bm.geodetic_to_cartesian(patch_pos_geo,bm)
patch_strike = input[1]
patch_dip = input[2]
patch_length = input[3]
patch_width = input[4]
slip = input[5]
Ps = [slippy.patch.Patch(p,l,w,s,d) for p,l,w,s,d in zip(patch_pos_cart,
                                                         patch_length,
                                                         patch_width,
                                                         patch_strike,
                                                         patch_dip)]
fig,ax = plt.subplots()
ax.set_title('left lateral')
bm.drawcoastlines(ax=ax)
q = ax.quiver(obs_pos_cart[:,0],obs_pos_cart[:,1],
              obs_disp[:,0],obs_disp[:,1],
              zorder=1,color='k',scale=1.0)
ax.quiver(pred_pos_cart[:,0],pred_pos_cart[:,1],
          pred_disp[:,0],pred_disp[:,1],
          zorder=1,color='m',scale=1.0)
          
ax.quiverkey(q,0.8,0.2,0.05,'0.05 m')
ps = slippy.patch.draw_patches(Ps,colors=slip[:,0],ax=ax,edgecolor='none',zorder=0)
fig.colorbar(ps,ax=ax)
fig,ax = plt.subplots()

ax.set_title('thrust')
bm.drawcoastlines(ax=ax)
q = ax.quiver(obs_pos_cart[:,0],obs_pos_cart[:,1],
              obs_disp[:,0],obs_disp[:,1],
              zorder=1,color='k',scale=1.0)
ax.quiver(pred_pos_cart[:,0],pred_pos_cart[:,1],
          pred_disp[:,0],pred_disp[:,1],
          zorder=1,color='m',scale=1.0)
          
ax.quiverkey(q,0.8,0.2,0.05,'0.05 m')
ps = slippy.patch.draw_patches(Ps,colors=slip[:,1],ax=ax,edgecolor='none',zorder=0)
fig.colorbar(ps,ax=ax)
plt.show()
quit()
