#!/usr/bin/env python
import slippy.io
import slippy.basis
import slippy.patch
import slippy.gbuild
import slippy.tikhonov
import matplotlib.pyplot as plt
import numpy as np
import scipy.optimize
import scipy.linalg
import sys

def reg_nnls(G,L,d):
  dext = np.concatenate((d,np.zeros(L.shape[0])))
  Gext = np.vstack((G,L))
  return scipy.optimize.nnls(Gext,dext)[0]

def main(config):
  ### load in all data
  ###################################################################
  if config["plotter"] == "basemap":
    import slippy.bm as plotting_library
  else:
    import slippy.xyz2geo as plotting_library

  # Repackage into a list of faults (each fault being a dict)
  fault_list = [];
  for key in config["faults"].keys():
    fault_segment = {
      "strike":config["faults"][key]["strike"],
      "dip":config["faults"][key]["dip"],
      "length":config["faults"][key]["length"],
      "width":config["faults"][key]["width"],
      "seg_pos_geo":config["faults"][key]["position"],
      "Nlength":config["faults"][key]["Nlength"],
      "Nwidth":config["faults"][key]["Nwidth"],
      "slip_basis":config["faults"][key]["basis"],
      "penalty":config["faults"][key]["penalty"]};
    fault_list.append(fault_segment)

  gps_input_file = config['gps_input_file']  
  insar_input_file = config['insar_input_file']
  leveling_input_file = config['leveling_input_file']  
  gps_output_file = config['gps_output_file']  
  insar_output_file = config['insar_output_file']  
  leveling_output_file = config['leveling_output_file']
  slip_output_file = config['slip_output_file']
  plotter = config['plotter']
  
  # The overall objects
  obs_disp_f = np.zeros((0,))
  obs_sigma_f = np.zeros((0,))
  obs_pos_geo_f = np.zeros((0,3))  # going to contain gps and insar obs
  obs_basis_f = np.zeros((0,3))
  obs_pos_total = np.zeros((0,3))  # for the basemap

  if gps_input_file is not None:
    gps_input = slippy.io.read_gps_data(gps_input_file)    
    Ngps = len(gps_input[0])
    obs_gps_pos_geo = gps_input[0]
    obs_gps_disp = gps_input[1]
    obs_gps_sigma = gps_input[2]
    obs_gps_basis = slippy.basis.cardinal_basis((Ngps,3))
    
    obs_disp_fi = obs_gps_disp.reshape((Ngps*3,))
    obs_sigma_fi = obs_gps_sigma.reshape((Ngps*3,))
    obs_basis_fi = obs_gps_basis.reshape((Ngps*3,3))
    obs_pos_geo_fi = obs_gps_pos_geo[:,None,:].repeat(3,axis=1).reshape((Ngps*3,3))
    
    obs_disp_f = np.concatenate((obs_disp_f,obs_disp_fi),axis=0)
    obs_sigma_f = np.concatenate((obs_sigma_f,obs_sigma_fi),axis=0)
    obs_basis_f = np.concatenate((obs_basis_f,obs_basis_fi),axis=0)    
    obs_pos_geo_f = np.concatenate((obs_pos_geo_f,obs_pos_geo_fi),axis=0)
    obs_pos_total = np.concatenate((obs_pos_total,obs_pos_geo_fi),axis=0)
    
  else:
    obs_gps_pos_geo = np.zeros((0,3))
    obs_gps_disp = np.zeros((0,3))
    obs_gps_sigma = np.zeros((0,3))
    obs_gps_basis = np.zeros((0,3,3))
    Ngps = 0

  if insar_input_file is not None:
    insar_input = slippy.io.read_insar_data(insar_input_file)    
    Ninsar = len(insar_input[0])
    obs_insar_pos_geo = insar_input[0]
    obs_insar_disp = insar_input[1]
    obs_insar_sigma = insar_input[2]
    obs_insar_basis = insar_input[3]

    obs_disp_f = np.concatenate((obs_disp_f,obs_insar_disp),axis=0)
    obs_sigma_f = np.concatenate((obs_sigma_f,obs_insar_sigma),axis=0)
    obs_basis_f = np.concatenate((obs_basis_f,obs_insar_basis),axis=0)    
    obs_pos_geo_f = np.concatenate((obs_pos_geo_f,obs_insar_pos_geo),axis=0)
    obs_pos_total = np.concatenate((obs_pos_total,obs_insar_pos_geo),axis=0)
  
  else:
    obs_insar_pos_geo = np.zeros((0,3))
    obs_insar_disp = np.zeros((0,))
    obs_insar_sigma = np.zeros((0,))
    obs_insar_basis = np.zeros((0,3))
    Ninsar = 0

  if leveling_input_file is not None: 
    leveling_input = slippy.io.read_insar_data(leveling_input_file)
    Nleveling = len(leveling_input[0])
    obs_leveling_pos_geo = leveling_input[0]
    obs_leveling_disp = leveling_input[1]
    obs_leveling_sigma = leveling_input[2]
    obs_leveling_basis = leveling_input[3]

    obs_pos_total = np.concatenate((obs_pos_total,obs_leveling_pos_geo),axis=0)
    obs_disp_f = np.concatenate((obs_disp_f,obs_leveling_disp),axis=0)
    obs_sigma_f = np.concatenate((obs_sigma_f,obs_leveling_sigma),axis=0)   
    obs_basis_f = np.concatenate((obs_basis_f,obs_leveling_basis),axis=0)    
    obs_pos_geo_f = np.concatenate((obs_pos_geo_f,obs_leveling_pos_geo),axis=0)
  
  else:
    obs_leveling_pos_geo = np.zeros((0,3))
    obs_leveling_disp = np.zeros((0,))
    obs_leveling_sigma = np.zeros((0,))
    obs_leveling_basis = np.zeros((0,3))
    Nleveling = 0    

  if gps_output_file is None:
    gps_output_file = sys.stdout

  if insar_output_file is None:
    insar_output_file = sys.stdout

  if leveling_output_file is None:
    leveling_output_file = sys.stdout 

  if slip_output_file is None:
    slip_output_file = sys.stdout

  # ###################################################################
  ### set up basemap for calculation
  ### discretize the fault segments
  ### create slip basis vectors for each patch  
  ### build regularization matrix
  ###################################################################
  patches = [];  # a growing list of fault patches. 
  slip_basis_f = np.zeros((0,3)) # a growing list of basis functions for slip patches
  total_fault_slip_basis=[];
  patches_f = [];
  L_array = [];
  Ns_total=0;
  
  # Set up the map for the calculation
  bm = plotting_library.create_default_basemap(obs_pos_total[:,0],obs_pos_total[:,1]) 
  obs_pos_cart_f = plotting_library.geodetic_to_cartesian(obs_pos_geo_f,bm)  
  
  # Fault processing
  for fault in fault_list:  
    # Convert fault to cartesian coordinates  
    fault["seg_pos_cart"] = plotting_library.geodetic_to_cartesian(fault["seg_pos_geo"],bm)

    # Discretize fault segment
    seg = slippy.patch.Patch(fault["seg_pos_cart"],
                             fault["length"],fault["width"],
                             fault["strike"],fault["dip"])
    single_fault_patches = np.array(seg.discretize(fault["Nlength"],fault["Nwidth"]))
    Ns = len(single_fault_patches);

    # Create slip basis vectors
    Ds = len(fault["slip_basis"])  # the number of basis vectors for this slip patch. 
    single_fault_slip_basis = np.array([fault["slip_basis"] for j in range(Ns)])  
    if total_fault_slip_basis==[]:
      total_fault_slip_basis = single_fault_slip_basis;
    else:
      total_fault_slip_basis=np.concatenate((total_fault_slip_basis, single_fault_slip_basis),axis=0)
    single_fault_silp_basis_f = single_fault_slip_basis.reshape((Ns*Ds,3))

    # Packaging of slip_basis_f
    single_fault_patches_f = single_fault_patches[:,None].repeat(Ds,axis=1).reshape((Ns*Ds,))
    patches=np.concatenate((patches,single_fault_patches),axis=0)
    patches_f=np.concatenate((patches_f,single_fault_patches_f),axis=0)
    slip_basis_f=np.concatenate((slip_basis_f, single_fault_silp_basis_f),axis=0);

    ### build regularization matrix
    L = np.zeros((0,Ns*Ds))
    indices = np.arange(Ns*Ds).reshape((Ns,Ds))
    for i in range(Ds): 
      connectivity = indices[:,i].reshape((fault["Nlength"],fault["Nwidth"]))
      Li = slippy.tikhonov.tikhonov_matrix(connectivity,2,column_no=Ns*Ds)
      L = np.vstack((Li,L))

    L *= fault["penalty"] 
    L_array.append(L)
    Ns_total = Ns_total+Ns

  # Build System Matrix
  G = slippy.gbuild.build_system_matrix(obs_pos_cart_f, 
                                        patches_f,
                                        obs_basis_f,
                                        slip_basis_f, 
                                        Nleveling) 
  

  if Nleveling>0:  # IF LEVELING: 
    ### build larger regularization matrix
    L_array.append([0]);
    L = scipy.linalg.block_diag(*L_array)   

    ### weigh system matrix and data by the uncertainty
    ###################################################################  
    G /= obs_sigma_f[:,None]
    obs_disp_f /= obs_sigma_f

    ### estimate slip and compute predicted displacement
    #####################################################################
    slip_f = reg_nnls(G,L,obs_disp_f)
    pred_disp_f = G.dot(slip_f)*obs_sigma_f 

    slip_f = slip_f[0:-1];  # LEVELING: Will ignore the last model parameter, which is the leveling offset
    slip = slip_f.reshape((Ns_total,Ds))  # THIS ASSUMES ALL FAULTS HAVE THE SAME NUMBER OF BASIS VECTORS
    cardinal_slip = slippy.basis.cardinal_components(slip,total_fault_slip_basis)

    # split predicted displacements into insar and GPS component  # LEVELING: Things just got more complicated
    pred_disp_f_gps = pred_disp_f[:3*Ngps]
    pred_disp_gps = pred_disp_f_gps.reshape((Ngps,3))
    pred_disp_insar = pred_disp_f[3*Ngps:3*Ngps+Ninsar]
    pred_disp_leveling = pred_disp_f[3*Ngps+Ninsar:];

  if Nleveling==0:
    # build regularization matrix
    L = scipy.linalg.block_diag(*L_array) 

    ### weigh system matrix and data by the uncertainty
    ###################################################################
    G /= obs_sigma_f[:,None]
    obs_disp_f /= obs_sigma_f

    ### estimate slip and compute predicted displacement
    #####################################################################
    slip_f = reg_nnls(G,L,obs_disp_f)
    pred_disp_f = G.dot(slip_f)*obs_sigma_f 
    slip = slip_f.reshape((Ns_total,Ds))  # THIS ASSUMES ALL FAULTS HAVE THE SAME NUMBER OF BASIS VECTORS
    cardinal_slip = slippy.basis.cardinal_components(slip,total_fault_slip_basis)

    # split predicted displacements into insar and GPS component 
    pred_disp_f_gps = pred_disp_f[:3*Ngps]
    pred_disp_gps = pred_disp_f_gps.reshape((Ngps,3))
    pred_disp_insar = pred_disp_f[3*Ngps:]  
    pred_disp_leveling = pred_disp_f[-1:-1]; # no leveling; padding with nothing

  ### get slip patch data for outputs
  #####################################################################
  patches_pos_cart =[i.patch_to_user([0.5,1.0,0.0]) for i in patches]
  patches_pos_geo = plotting_library.cartesian_to_geodetic(patches_pos_cart,bm)
  patches_strike = [i.strike for i in patches]
  patches_dip = [i.dip for i in patches]
  patches_length = [i.length for i in patches]
  patches_width = [i.width for i in patches]

  ### write output
  #####################################################################
  slippy.io.write_slip_data(patches_pos_geo,
                            patches_strike,patches_dip,
                            patches_length,patches_width,
                            cardinal_slip,slip_output_file)

  slippy.io.write_gps_data(obs_gps_pos_geo, 
                           pred_disp_gps,0.0*pred_disp_gps,
                           gps_output_file)

  slippy.io.write_insar_data(obs_insar_pos_geo,
                           pred_disp_insar,0.0*pred_disp_insar,
                           obs_insar_basis, 
                           insar_output_file)

  slippy.io.write_insar_data(obs_leveling_pos_geo,
                           pred_disp_leveling,0.0*pred_disp_leveling,
                           obs_leveling_basis, 
                           leveling_output_file)

  return
