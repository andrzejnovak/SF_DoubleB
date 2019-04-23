import os, sys, errno
import copy
import numpy as np
import pandas as pd
import time
import pprint
import ROOT as R
from get_SF_cfit import runSF_x
import pickle
pd.set_option('display.max_columns', 500)
pd.set_option('display.width', 1000)


import argparse
parser = argparse.ArgumentParser(description='Process some integers.')
parser.add_argument('--load', action='store_true', default=False, help="Use this if you've ran the first part, but crashed when computing the final systematics, it saves the intermediate steps.")
parser.add_argument('--load2', action='store_true', default=False, help="")
parser.add_argument('--test', action='store_true', default=False, help="")
parser.add_argument('--nomonly', action='store_true', default=False, help="")
parser.add_argument('-n', '--name', default='DDBvL', help="Name of the tagger, to load WPs etc...")
parser.add_argument('-y', '--year', default='2016', help="Year of the run...")
parser.add_argument('-o', '--outname', default=None, help="Name to append to out files...")
args = parser.parse_args()
if args.outname == None:
	args.outname = args.name

if "DDCvL" in args.name:
	ccSignal = True
else: ccSignal = False

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


# Get nominal SF
def step0(WP=None, pt_bin_ix=0, systname=None, addSYS=True):
	glue=True; calcSYS=False; # Glue templates to match final SF, add systematics to correlate, don't actually calcualate errors.
	
	pt_bin = pt_bins[pt_bin_ix]
	start = time.time()
	if int(pt_bin.split("to")[-1]) <= 350:
		file_name = r+name1+root
		SF, pars, chi2 = runSF_x(file_name, pt_bins, m, WP, glue=glue,  addSYS=addSYS, calcSYS=calcSYS, systlist=syst_list, ccSignal=ccSignal)
	else:
		file_name = r+name2+root
		SF, pars, chi2 = runSF_x(file_name,
			pt_bins, pt_bin_ix, WP, glue=glue,  addSYS=addSYS, calcSYS=calcSYS, LTSV=LTSV, systlist=syst_list, systname=systname, ccSignal=ccSignal) 

	#print "		Time to run: ", np.round((time.time() - start)/60, 2), "min"
	return SF

# Calculate up/down shape SFs and store to dictionary	
def step1(templates=templates, WP=None):
	glue=True; addSYS=True; calcSYS=False; # Glue templates to match final SF, add systematics to correlate, don't actually calcualate errors.
	print WP
	for n, template in enumerate(templates):		
		for m, pt_bin in enumerate(pt_bins):
			start = time.time()
			if int(pt_bin.split("to")[-1]) <= 350:
				file_name = r1+name1+"_"+WP+"_"+template+root
				SF, pars, chi2 = runSF_x(file_name, pt_bins, m, WP, glue=glue,  addSYS=addSYS, calcSYS=calcSYS, systlist=syst_list)
			else:
				file_name = r+name2+"_"+template+root		
				print file_name
				SF, pars, chi2 = runSF_x(file_name, pt_bins, m, WP, glue=glue, addSYS=addSYS, calcSYS=calcSYS, LTSV=LTSV, systlist=syst_list)
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
			SF, pars, chi2 = runSF_x(file_name, pt_bins, m, WP, glue=glue,  addSYS=addSYS, calcSYS=calcSYS, LTSV=LTSV, systlist=syst_list)
		else: 
			#file_name = r+name2+WP+root
			file_name = r+name2+root
			SF, pars, chi2 = runSF_x(file_name, pt_bins, m, WP, glue=glue,  addSYS=addSYS, calcSYS=calcSYS, LTSV=LTSV, systlist=syst_list)
		
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
			SF, pars, chi2 = runSF_x(file_name, pt_bins, m, WP, glue=glue,  addSYS=addSYS, calcSYS=calcSYS, systlist=syst_list)
			
		else:
			file_name = rjp+"Run2017BCDEF_ReReco_QCDMuonEnriched_AK8Jet300orAK4Jet300_Pt350_Final_DoubleMuonTaggedFatJets_histograms_btagval_allVars_ptReweighted_dataUseMCJPcalib_SysMerged_SFtemplates_"+WP+root
			SF, pars, chi2 = runSF_x(file_name, pt_bins, m, WP, glue=glue, addSYS=addSYS, calcSYS=calcSYS, systlist=syst_list)
		
		print SF				
		print "Time to run: ", np.round((time.time() - start)/60, 2), "min"

		SF_dicts['SF_dict_{}'.format(WPs.index(WP))]['SF_JP'].append(float(str(SF)))

# Calculate actual SF values, include all systematics
def step3(WP=None, pt_bin_ix=0, SF_dict=SF_dict_empty):
	glue=True; calcSYS=True; addSYS=True

	pt_bin = pt_bins[pt_bin_ix]
	start = time.time()
	if int(pt_bin.split("to")[-1]) <= 350:
		file_name = r1+name1+"_"+WP+root
		file_name = r+name2+root
		SF, sigma_stat, syst_up, syst_down, variances_names, errors, variances, nom_pars, chi2 = runSF_x(file_name, 
			pt_bins, pt_bin_ix, WP, glue=glue,  addSYS=addSYS, calcSYS=calcSYS, SF_dict=SF_dict, LTSV=LTSV, systlist=syst_list)
	else:
		file_name = r+name2+root
		SF_dict = SF_dicts['SF_dict_{}'.format(WPs.index(WP))]
		SF, sigma_stat, syst_up, syst_down, variances_names, errors, variances, nom_pars, chi2 = runSF_x(file_name, 
			pt_bins, pt_bin_ix, WP, glue=glue,  addSYS=addSYS, calcSYS=calcSYS, SF_dict=SF_dict, LTSV=LTSV, systlist=syst_list) 
		
	if LTSV: 
		SF_SV = (nom_pars[5]/nom_pars[4])
		SF_JP = (nom_pars[2]/nom_pars[0])
		tageff = (nom_pars[3]/nom_pars[2])
		preeff = (nom_pars[1]/nom_pars[0])
	else: 
		SF_SV = 0.
		SF_JP = 0.
		tageff = 0.
		preeff = 0.	
		
	print "Time to run: ", np.round((time.time() - start)/60, 2), "min"
	return SF, sigma_stat, syst_up, syst_down, variances_names, errors, variances, SF_SV, preeff, tageff, SF_JP


if __name__ == "__main__":
	# Specify input files
	LTSV = True
	root = '.root'

	#r = 'col4/'
	r = 'col_fin450/'
	name2 = "collated_norm"+args.name

	# Prepare for outside systematic (SF) calculations
	all_syst_list = ["JES","BFRAG", "CFRAG","K0L", "PU", "CDFRAG", "NTRACKS", "GCC", "GBB"]
	#all_syst_list = ["BFRAG", "CFRAG","K0L", "CDFRAG", "NTRACKS", "GCC", "GBB"]
	_f = R.TFile(r+name2+root)
	sys_avail = set([_k.GetName().split("_")[-1].replace("up","").replace("down", "").replace("_", "") for _k in _f.GetListOfKeys()])
	syst_list = list(set(all_syst_list) & set(sys_avail))	

	# Read WPs from JSON appropriate to tagger name
	import json
	WPs, WPs_value = [], []
	tagger = args.name.split("_")[-1]
	with open('DDX.json') as json_file:
		data = json.load(json_file)[tagger]
		for key, value in sorted(data.iteritems(), key=lambda (k,v): (v,k)):
			WPs.append(tagger+str(key))
			WPs_value.append(value)

	# Specify binning
	#pt_bins = ['pt350to430', 'pt430to2000', 'pt350to2000']
	pt_bins = ['pt350to450', 'pt450to2000', 'pt350to2000']

	SF_dicts = {}
	for i in range(len(WPs)):
		SF_dicts['SF_dict_{}'.format(i)] = copy.deepcopy(SF_dict_empty)

	
	if args.test:
		testSFs = {}
		for ix, pt in enumerate(pt_bins):
			for i, WP in enumerate(WPs): 
				testSFs[pt+WP] = {}
				for _sys in [None]+syst_list:
					if _sys==None: 
						SF_nom = step0(WP=WP, pt_bin_ix=ix, systname=_sys, addSYS=False)	
						sys_n = "Nominal"
					else:
						SF_nom = step0(WP=WP, pt_bin_ix=ix, systname=_sys, addSYS=True)	
						sys_n = _sys
					testSFs[pt+WP][sys_n] = SF_nom
					print pt, WP, _sys, SF_nom
		pprint.pprint(testSFs)
		sys.exit()

	# Load intermediate results - steps: 1,2, 2_1
	if args.load:
		SF_dicts = {}
		with open('SF_dicts_{}.pkl'.format(args.name), 'r') as f:
			SF_dicts = pickle.load(f)

	# If running from start - run all pre SFs
	else:
		if not args.nomonly:
			pass
			#for WP in WPs: M = step1(WP=WP)
			#for WP in WPs: M = step2(WP=WP)
			#for WP in WPs: M = step2_1(WP=WP)

	pprint.pprint(SF_dicts)

	with open('SF_dicts_{}.pkl'.format(args.name), 'w') as f:
		pickle.dump(SF_dicts, f)

	# Run final SF calculation
	headers = ["Systematic"]
	arrays = []
	reprint = []
	parinfo = []

	ptwp = [WP[len(tagger):]+pt for pt in pt_bins for WP in WPs]
	df = pd.DataFrame()
	df2 = pd.DataFrame(
		index=ptwp,
		columns=['SF', 'SF_nom', 'Stat', 'Sys up', 'Sys down', 'SF_SV', 'SF_JP', 'preeff', 'tageff', 'Combined up', 'Combined down'] )
	
	for ix, pt in enumerate(pt_bins):
		for i, WP in enumerate(WPs): 	
			row_name = WP[len(tagger):]+pt
			WP_dict = SF_dicts['SF_dict_{}'.format(i)]

			if args.load2:
				df2 = pd.read_csv('DF_res{}.csv'.format(args.name), index_col=0)
				df = pd.read_csv('DF_sys{}.csv'.format(args.name), index_col=0)

			if np.isnan(df2.loc[row_name, 'SF']) :
				print df2
				if args.nomonly:
					SF_nom = step0(WP=WP, pt_bin_ix=ix)	
					SF, sigma_stat, syst_up, syst_down, SF_SV, SF_JP, preeff, tageff = 0,0,0,0,0,0,0,0
				else:
					SF_nom = step0(WP=WP, pt_bin_ix=ix)	
					SF, sigma_stat, syst_up, syst_down, variances_names, errors, variances, SF_SV, preeff, tageff, SF_JP = step3(WP=WP, pt_bin_ix=ix, SF_dict=WP_dict)	

					df = df.append(dict(zip(variances_names, errors)), ignore_index=True)
					df.to_csv('DF_sys{}.csv'.format(args.name))

				df2.loc[row_name] = [SF, SF_nom, sigma_stat, syst_up, syst_down, SF_SV, SF_JP, preeff, tageff, 
								np.sqrt(sigma_stat**2 + syst_up**2),
								np.sqrt(sigma_stat**2 + syst_down**2)
				]
			 				
				df2.to_csv('DF_res{}.csv'.format(args.name))
				
			else:
				continue
			
			#if args.load2:
			# 	df2 = pd.read_csv('DF_res{}.csv'.format(args.outname))
			print df
			print df2
			
			#df.to_csv('DF_sys{}.csv'.format(args.outname))


