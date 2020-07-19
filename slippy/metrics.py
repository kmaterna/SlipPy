#!/usr/bin/env python
# After you've done an inversion, what are the results? 
# How much moment is created, and how big is the misfit? 
# Writes into a summary text file, and prints to screen. 
# Useful for L-curve analysis.

import numpy as np 
import sys
import json
import slippy.io


# -------- INPUT FUNCTIONS ----------- # 
def welcome_and_parse(argv):
	print("Metrics for this inversion:")
	if len(argv)<2:
		print("Error! Please provide the name of a config json. Exiting. "); sys.exit(0);
	else:
		configname = argv[1];
	return configname;

def parse_json(configname):
	config_file = open(configname,'r');
	config = json.load(config_file);
	for key in config:  # adding the output directory onto the output files, for ease of use. 
		if "output" in key and key != "output_dir":
			config[key]=config["output_dir"]+config[key];
	config["summary_file"]=config["output_dir"]+"summary_stats.txt";  # Creates an output file
	return config;


# -------- COMPUTE FUNCTIONS ----------- # 
def get_slip_moment(slip_filename):
	# From the inversion results, what is the moment of the slip distribution? 
	moment_total = 0;
	mu=30e9;  # Pa, assumed. 
	length, width, leftlat, thrust, tensile = np.loadtxt(slip_filename,skiprows=1,unpack=True,usecols=(5,6,7,8,9));
	for i in range(len(length)):
		slip = np.sqrt(leftlat[i]*leftlat[i] + thrust[i]*thrust[i]);
		area = length[i]*width[i]; # m^2
		momenti = moment_from_muad(mu, area, slip);
		moment_total = moment_total+momenti;
	mw = mw_from_moment(moment_total);
	print("Calculating metrics for inversion results.");
	print("Calculating moment from %s" % slip_filename);
	return [moment_total, mw];

def get_total_misfit(config):
	lev_misfit, lev_norm_misfit, insar_misfit, insar_norm_misfit, gps_misfit, gps_norm_misfit = None, None, None, None, None, None;
	lev_npts, insar_npts, gps_npts = None, None, None;
	if "observed_leveling_file" in config.keys():
		[lev_misfit, lev_norm_misfit, lev_npts] = get_misfit_leveling(config["observed_leveling_file"],config["predicted_leveling_file"]);
	if "observed_insar_file" in config.keys():
		[insar_misfit, insar_norm_misfit, insar_npts] = get_misfit_insar(config["observed_insar_file"],config["predicted_insar_file"]);
	if "observed_gps_file" in config.keys():
		[gps_misfit, gps_norm_misfit, gps_npts] = get_misfit_gps(config["observed_gps_file"],config["predicted_gps_file"]);
	return [gps_misfit, gps_norm_misfit, gps_npts, insar_misfit, insar_norm_misfit, insar_npts, lev_misfit, lev_norm_misfit, lev_npts];

def get_misfit_gps(obs_file, pred_file):
	# Misfit from each data pair (GPS, UAVSAR, Leveling, S1, TSX)
	# Want in both absolute numbers and relative to the respective uncertainties. 
	gps_input = slippy.io.read_gps_data(obs_file); 
	gps_pred = slippy.io.read_gps_data(pred_file); 
	abs_misfit = np.abs(gps_input[1]-gps_pred[1]);
	norm_misfit = np.divide(abs_misfit, gps_input[2]);  # divide by sigma
	mean_average_misfit = np.nanmean(abs_misfit);
	mean_norm_average_misfit = np.nanmean(norm_misfit);
	npts = len(gps_input[1]);
	return [mean_average_misfit, mean_norm_average_misfit, npts];

def get_misfit_insar(obs_file, pred_file):
	# Misfit from each data pair (GPS, UAVSAR, Leveling, S1, TSX)
	# Want in both absolute numbers and relative to the respective uncertainties. 
	insar_input = slippy.io.read_insar_data(obs_file)
	insar_pred = slippy.io.read_insar_data(pred_file)
	abs_misfit = np.abs(insar_input[1]-insar_pred[1]);
	norm_misfit = np.divide(abs_misfit, insar_input[2]);  # divide by sigma	
	mean_average_misfit = np.nanmean(abs_misfit);
	mean_norm_average_misfit = np.nanmean(norm_misfit);
	npts = len(insar_input[1]);
	return [mean_average_misfit, mean_norm_average_misfit, npts];

def get_misfit_leveling(obs_file, pred_file):
	# Misfit from each data pair (GPS, UAVSAR, Leveling, S1, TSX)
	# Want in both absolute numbers and relative to the respective uncertainties. 
	leveling_input = slippy.io.read_insar_data(obs_file)
	leveling_pred = slippy.io.read_insar_data(pred_file)
	abs_misfit = np.abs(leveling_input[1]-leveling_pred[1]);
	norm_misfit = np.divide(abs_misfit, leveling_input[2]);  # divide by sigma	
	mean_average_misfit = np.nanmean(abs_misfit);
	mean_norm_average_misfit = np.nanmean(norm_misfit);
	npts = len(leveling_input[1]);
	return [mean_average_misfit, mean_norm_average_misfit, npts];


# -------- WRITE FUNCTIONS ----------- # 
def write_outputs(moments, misfits, outfile):
	ofile=open(outfile,'w');  # cleaning the file from last times

	# moments
	scinot = "{:e}".format(moments[0]);
	print("Total Slip Moment is %s N-m, equivalent to mw=%f \n" % (scinot, moments[1]) );
	ofile.write("Total Slip Moment is %s N-m, equivalent to mw=%f \n" % (scinot, moments[1]) );

	if misfits[0] is not None:
		print("Average GPS misfit: %f mm" % (1000*misfits[0]) );
		print("Average normalized GPS misfit: %f sigma \n" % (misfits[1]) );
		ofile.write("Average GPS misfit: %f mm\n" % (1000*misfits[0]) );
		ofile.write("Average normalized GPS misfit: %f sigma \n" % (misfits[1]) );
		ofile.write("GPS npts: %d \n" % misfits[2]);
	
	if misfits[3] is not None:
		print("Average InSAR misfit: %f mm" % (1000*misfits[3]) );
		print("Average normalized InSAR misfit: %f sigma \n" % (misfits[4]) );
		ofile.write("Average InSAR misfit: %f mm\n" % (1000*misfits[3]) );
		ofile.write("Average normalized InSAR misfit: %f sigma \n" % (misfits[4]) );
		ofile.write("InSAR npts: %d \n" % misfits[5]);

	if misfits[6] is not None:
		print("Average Leveling misfit: %f mm" % (1000*misfits[6]) );
		print("Average normalized Leveling misfit: %f sigma \n" % (misfits[7]) );	
		ofile.write("Average Leveling misfit: %f mm\n" % (1000*misfits[6]) );
		ofile.write("Average normalized Leveling misfit: %f sigma \n" % (misfits[7]) );
		ofile.write("Leveling npts: %d \n" % misfits[8]);
	
	ofile.close();
	return;

# -------- MATH FUNCTIONS ----------- # 
def moment_from_muad(mu, A, d):
	# moment = mu * A * d
	return mu*A*d;

def mw_from_moment(moment):
	# Takes newton meters, returns a moment magnitude
	moment = moment*1e7;
	mw = (2/3)*np.log10(moment) - 10.7
	return mw;


# -------- ACCESS FUNCTIONS ----------- # 
def main_function(config):
	misfits = get_total_misfit(config);
	moments = get_slip_moment(config["slip_output_file"]);
	write_outputs(moments, misfits, config["summary_file"]);
	return;


if __name__=="__main__":
	configfile = welcome_and_parse(sys.argv);  # expects a config file passed in by argument
	# main_function(configfile=configfile);


