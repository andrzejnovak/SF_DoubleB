import os, sys, errno
import copy
import numpy as np
import pandas as pd
import time
from get_SF_cfit import runSF_x
import pickle
pd.set_option('display.max_columns', 500)
pd.set_option('display.width', 1000)


import argparse
parser = argparse.ArgumentParser(description='Process some integers.')
parser.add_argument('--load', action='store_true', default=False, help="Use this if you've ran the first part, but crashed when computing the final systematics, it saves the intermediate steps.")
parser.add_argument('-n', '--name', default='DDBvL', help="Name of the tagger, to load WPs etc...")
parser.add_argument('-o', '--outname', default=None, help="Name to append to out files...")
args = parser.parse_args()
if args.outname == None:
	args.outname = args.name

def make_dirs(dirname):
    """
    Ensure that a named directory exists; if it does not, attempt to create it.
    """
    try:
        os.makedirs(dirname)
    except OSError, e:
        if e.errno != errno.EEXIST:
            raise

###################################
###################################
###################################
# Prepare for outside systematic (SF) calculations

# Specify up and down shape systematics
templates = ['b_up', 'b_down', 'c_up', 'c_down', 'cfromg_up', 'cfromg_down', 'l_up', 'l_down']

#Make a dictionary to record intermediate SF values
SF_dict_empty = {
	'SF_b_down' : [],
	'SF_b_up' : [], 
	'SF_c_down' : [], 
	'SF_c_up' : [], 
	'SF_cfromg_down' : [], 
	'SF_cfromg_up' : [], 
	'SF_l_down' : [], 
	'SF_l_up' : [], 
	'SF_5_temp' : [], 
	'SF_JP' : [] }

# Calculate up/down shape SFs and store to dictionary
def step1(templates=templates, WP=None):
	glue=True; addSYS=True; calcSYS=False; # Glue templates to match final SF, add systematics to correlate, don't actually calcualate errors.
	print WP
	for n, template in enumerate(templates):		
		for m, pt_bin in enumerate(pt_bins):
			start = time.time()
			if int(pt_bin.split("to")[-1]) <= 350:
				file_name = r1+name1+"_"+WP+"_"+template+root
				SF, pars, chi2 = runSF_x(file_name, pt_bins, m, WP, glue=glue,  addSYS=addSYS, calcSYS=calcSYS)
			else:
				file_name = r+name2+"_"+template+root		
				print file_name
				SF, pars, chi2 = runSF_x(file_name, pt_bins, m, WP, glue=glue, addSYS=addSYS, calcSYS=calcSYS, LTSV=LTSV)
				#SF, pars, chi2 = runSF_x(file_name, pt_bins, m, WP, merge=merge, glue=glue, addSYS=addSYS, calcSYS=calcSYS)
			
			print "		", SF	
			print "		Time to run: ", np.round((time.time() - start)/60, 2), "min"

			SF_dicts['SF_dict_{}'.format(WPs.index(WP))]['SF_'+template].append(float(str(SF)))
	return

# Calculate systematic of fitting all tempaltes separately regardeless of similar shapes
def step2(WP=None):
	glue=False; addSYS=True; calcSYS=False; 
	for m, pt_bin in enumerate(pt_bins):
		start = time.time()
		if int(pt_bin.split("to")[-1]) <= 350:
			file_name = r+name1+root
			SF, pars, chi2 = runSF_x(file_name, pt_bins, m, WP, glue=glue,  addSYS=addSYS, calcSYS=calcSYS, LTSV=LTSV)
		else: 
			#file_name = r+name2+WP+root
			file_name = r+name2+root
			SF, pars, chi2 = runSF_x(file_name, pt_bins, m, WP, glue=glue,  addSYS=addSYS, calcSYS=calcSYS, LTSV=LTSV)
		
		print "		", SF	
		print "      Time to run: ", np.round((time.time() - start)/60, 2), "min"
		
		SF_dicts['SF_dict_{}'.format(WPs.index(WP))]['SF_5_temp'].append(float(str(SF)))
	return 

# Calculate MCJP Calib - different input files, data is modified
def step2_1(WP=None):
	glue=True; addSYS=True;  calcSYS=False
	for m, pt_bin in enumerate(pt_bins):
		start = time.time()
		if int(pt_bin.split("to")[-1]) <= 350:
			file_name = rjp+"Run2017BCDEF_ReReco_QCDMuonEnriched_AK4DiJet170_Pt250_Final_DoubleMuonTaggedFatJets_histograms_btagval_allVars_ptReweighted_dataUseMCJPcalib_SysMerged_SFtemplates_"+WP+root
			SF, pars, chi2 = runSF_x(file_name, pt_bins, m, WP, glue=glue,  addSYS=addSYS, calcSYS=calcSYS)
			
		else:
			file_name = rjp+"Run2017BCDEF_ReReco_QCDMuonEnriched_AK8Jet300orAK4Jet300_Pt350_Final_DoubleMuonTaggedFatJets_histograms_btagval_allVars_ptReweighted_dataUseMCJPcalib_SysMerged_SFtemplates_"+WP+root
			SF, pars, chi2 = runSF_x(file_name, pt_bins, m, WP, glue=glue, addSYS=addSYS, calcSYS=calcSYS)
		
		print SF				
		print "Time to run: ", np.round((time.time() - start)/60, 2), "min"

		SF_dicts['SF_dict_{}'.format(WPs.index(WP))]['SF_JP'].append(float(str(SF)))

# Calculate actual SF values, include all systematics
def step3(WP=None, SF_dict=SF_dict_empty):
	glue=True; calcSYS=True; addSYS=True
	# Make some arrays to catch outputs and log them
	SFs, sigma_stats, syst_ups, syst_downs = [], [], [], []
	errors_all, variances_all = [], []
	SF_SV, SF_JP, preeff, tageff = [], [], [], []
	for m, pt_bin in enumerate(pt_bins):
		start = time.time()
		if int(pt_bin.split("to")[-1]) <= 350:
			file_name = r1+name1+"_"+WP+root
			file_name = r+name2+root
			SF, sigma_stat, syst_up, syst_down, variances_names, errors, variances, nom_pars, chi2 = runSF_x(file_name, 
				pt_bins, m, WP, glue=glue,  addSYS=addSYS, calcSYS=calcSYS, SF_dict=SF_dict, LTSV=LTSV)
		else:
			file_name = r+name2+root
			print file_name
			SF_dict = SF_dicts['SF_dict_{}'.format(WPs.index(WP))]
			SF, sigma_stat, syst_up, syst_down, variances_names, errors, variances, nom_pars, chi2 = runSF_x(file_name, 
				pt_bins, m, WP, glue=glue,  addSYS=addSYS, calcSYS=calcSYS, SF_dict=SF_dict, LTSV=LTSV) 
		
		errors_all.append(errors)
		variances_all.append(variances)
		SFs.append(np.round(SF, 2))
		sigma_stats.append(np.round(sigma_stat, 3))
		syst_ups.append(np.round(syst_up, 3))
		syst_downs.append(np.round(syst_down, 3))
		#pars = [par, par_tag, parJPtagged,par_tagJPtagged, parSV, par_tagSV]
		if LTSV: 
			SF_SV.append(nom_pars[5]/nom_pars[4])
			SF_JP.append(nom_pars[2]/nom_pars[0])
			tageff.append(nom_pars[3]/nom_pars[2])
			preeff.append(nom_pars[1]/nom_pars[0])
		else: 
			SF_SV.append(0.)
			SF_JP.append(0.)
			tageff.append(0)
			preeff.append(0)
		
		
		
		
		print "Time to run: ", np.round((time.time() - start)/60, 2), "min"
	return SFs, sigma_stats, syst_ups, syst_downs, variances_names, errors_all, variances_all, SF_SV, preeff, tageff, SF_JP


if __name__ == "__main__":
	# Specify input files
	LTSV = True
	root = '.root'

	r = 'col3/'
	name2 = "collated_normRun2016_"+args.name

	# Read WPs from JSON appropriate to tagger name
	import json
	WPs, WPs_value = [], []
	with open('DDX.json') as json_file:
		data = json.load(json_file)[args.name]
		for key, value in sorted(data.iteritems(), key=lambda (k,v): (v,k)):
			WPs.append(args.name+str(key))
			WPs_value.append(value)

	# Specify binning
	#pt_bins = ['pt250to350','pt350to430', 'pt430to2000']
	pt_bins = ['pt350to430', 'pt430to2000']

	SF_dicts = {}
	for i in range(len(WPs)):
		SF_dicts['SF_dict_{}'.format(i)] = copy.deepcopy(SF_dict_empty)

	# Load intermediate results - steps: 1,2, 2_1
	if args.load:
		SF_dicts = {}
		with open('SF_dicts.pkl'.format(i), 'r') as f:
			SF_dicts = pickle.load(f)

	# If running from start - run all pre SFs
	else:
		pass
		#for WP in WPs: M = step2(WP=WP)
		#sys.exit()
		for WP in WPs: M = step1(WP=WP)
		for WP in WPs: M = step2(WP=WP)
		#for WP in WPs: M = step2_1(WP=WP)

	print SF_dicts
	import pprint
	pprint.pprint(SF_dicts)

	with open('SF_dicts.pkl'.format(i), 'w') as f:
		pickle.dump(SF_dicts, f)

	# Run final SF calculation
	headers = ["Systematic"]
	arrays = []
	reprint = []
	parinfo = []
	for i, WP in enumerate(WPs): 
		WP_dict = SF_dicts['SF_dict_{}'.format(i)]
		SFs, sigma_stats, syst_ups, syst_downs, variances_names, errors_allpt, variances_allpt, SF_SV, preeff, tageff, SF_JP = step3(WP=WP, SF_dict=WP_dict)	
		#SFs, sigma_stats, syst_ups, syst_downs, SF_SV, preeff, tageff = [0]*len(pt_bins), [0]*len(pt_bins), [0]*len(pt_bins),[0]*len(pt_bins), [0]*len(pt_bins), [0]*len(pt_bins), [0]*len(pt_bins) # temp
		#SFs = step2(WP=WP) # temp 
		#reprint.append([WP[len("DoubleB"):], SFs, sigma_stats, syst_ups, syst_downs, SF_SV, preeff, tageff])
		tagger = args.name
		reprint.append([WP[len(tagger):], SFs, sigma_stats, syst_ups, syst_downs, SF_SV, SF_JP,preeff, tageff])
		for pt in pt_bins:
			#headers.append(WP[len("DoubleB"):]+pt)
			headers.append(WP[len(tagger):]+pt)

		names = variances_names
		for ar in errors_allpt:
			arrays.append(ar)	


	# Store outputs in data frames
	results = []
	for row in reprint:
		results.append(np.array(row[1:]))

	df = pd.DataFrame(np.array(arrays).transpose(), columns=headers[1:], index=names)
	#print df
	print headers
	df2 = pd.DataFrame(np.concatenate((results), axis=1), columns=headers[1:], index=['SF', 'Stat', 'Sys up', 'Sys down', 'SF_SV', 'SF_JP', 'preeff', 'tageff'] )
	df2.loc['Combined up'] = np.sqrt(df2.loc['Stat'].values**2 + df2.loc['Sys up'].values**2)
	df2.loc['Combined down'] = np.sqrt(df2.loc['Stat'].values**2 + df2.loc['Sys down'].values**2)
	print df2

	df.to_csv('DF_sys{}.csv'.format(args.outname))
	df2.to_csv('DF_res{}.csv'.format(args.outname))

