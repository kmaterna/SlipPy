#!/usr/bin/env python
from slippy.okada import dislocation
import numpy as np

def build_system_matrix(pos,patches,disp_directions,slip_directions,Nleveling=0,leveling=False, leveling_offset_sign=1):
  ''' 
  builds the system matrix 

  Parameters
  ----------
    pos : (N,3) array of surface observation points
    
    patches : (M,) list of Patch instances
    
    disp_direction : (N,3) array
      displacement direction
      
    slip_direction : (M,3) array
      slip directions

    Nleveling : Integer, how many of the observations are leveling 
      (assumed to be at the bottom of the vector)
  
  Returns
  -------
    out : (N,M) array of dislocation greens functions

  '''
  pos = np.asarray(pos)
  slip_directions = np.asarray(slip_directions)
  disp_directions = np.asarray(disp_directions)
  ifleveling = Nleveling>0 or leveling
  
  G = np.zeros((len(pos),len(patches)+ifleveling))  # one extra column if leveling

  for i,p in enumerate(patches):
    top_center = p.patch_to_user([0.5,1.0,0.0])
    disp,derr = dislocation(pos,slip_directions[i],top_center,
                            p.length,p.width,p.strike,p.dip)
    disp = np.einsum('...j,...j',disp,disp_directions)
    G[:,i] = disp

  # Add an extra column containing the offset parameter for leveling
  if ifleveling:
    vector_of_ones = np.zeros(np.shape(G[:,0]));
    for i in range(Nleveling):
      vector_of_ones[i]=leveling_offset_sign;
    vector_of_ones=np.flipud(vector_of_ones);
    G[:,len(patches)]=vector_of_ones  
                                
  return G  
  