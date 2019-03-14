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
args = parser.parse_args()

def make_dirs(dirname):
    """
    Ensure that a named directory exists; if it does not, attempt to create it.
    """
    try:
        os.makedirs(dirname)
    except OSError, e:
        if e.errno != errno.EEXIST:
            raise

# Probably needs to be hot coded, when a fit fails due to empty template
# Specify bin and WP when to use merged templates, returns input to run_SF_x
def when_to_merge(WP, bins, bin_idx):
	print "NBINS: ", len(bins)
	merge = False
	# if len(bins) == 4:
	# 	if WP == 'DoubleBM2' and bin_idx in [1]: merge=2
	# 	elif WP == 'DoubleBM2' and bin_idx in [2]: merge=True
	# 	elif WP == 'DoubleBH' and bin_idx in [1]: merge=2
	# 	elif WP == 'DoubleBH' and bin_idx in [0]: merge=True
	# 	else: merge = False
	# if len(bins) == 3:
	# 	if WP == 'DoubleBM2' and bin_idx in [1]: merge=True
	# 	elif WP == 'DoubleBH' and bin_idx in [0,1]: merge=True
	# 	else: merge = False
	# if len(bins) == 2:
	# 	if WP == 'DoubleBM2' and bin_idx in []: merge=True
	# 	elif WP == 'DoubleBH' and bin_idx in [0]: merge=True
	# 	else: merge = False
	if len(bins) == 3:
		if WP == 'DeepDoubleXH' and bin_idx in [0]: merge=True
		#elif WP == 'DoubleBH' and bin_idx in [0,1]: merge=True
		else: merge = False

	return merge

###################################
###################################
###################################
# Prepare for outside systematic (SF) calculations

# Specify up and down shape systematics
templates = ['b_0p5', 'b_1p5', 'c_0p5', 'c_1p5', 'cfromg_0p5', 'cfromg_1p5', 'l_0p5', 'l_1p5']

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

# Make a copy for each WP
SF_dict_DoubleBL = copy.deepcopy(SF_dict_empty)
SF_dict_DoubleBM1 = copy.deepcopy(SF_dict_empty)
SF_dict_DoubleBM2 = copy.deepcopy(SF_dict_empty)
SF_dict_DoubleBT = copy.deepcopy(SF_dict_empty)
SF_dict_DDBvLL = copy.deepcopy(SF_dict_empty)
SF_dict_DDBvLM1 = copy.deepcopy(SF_dict_empty)
SF_dict_DDBvLM2 = copy.deepcopy(SF_dict_empty)
SF_dict_DDBvLT1 = copy.deepcopy(SF_dict_empty)
SF_dict_DDBvLT2 = copy.deepcopy(SF_dict_empty)

# Calculate up/down shape SFs and store to dictionary
WP = "BL"
def step1(templates=templates, WP=WP):
	glue=True; addSYS=True; calcSYS=False; # Glue templates to match final SF, add systematics to correlate, don't actually calcualate errors.
	print WP
	for n, template in enumerate(templates):		
		for m, pt_bin in enumerate(pt_bins):
			start = time.time()
			print WP, pt_bin, template
			merge = False
			merge = when_to_merge(WP, pt_bins, m)
			if int(pt_bin.split("to")[-1]) <= 350:
				file_name = r1+name1+"_"+WP+"_"+template+root
				SF, pars, chi2 = runSF_x(file_name, pt_bins, m, WP, merge=merge, glue=glue,  addSYS=addSYS, calcSYS=calcSYS)
			else:
				#file_name = r1+name2+"_"+WP+"_"+template+root			
				file_name = "DeepDoubleX_PtMinIs350_plain_test_dataUseMCJPcalib_SysMerged_SFtemplates_DDX_H_WP.root"			
				SF, pars, chi2 = runSF_x(file_name, pt_bins, m, WP, merge=merge, glue=glue, addSYS=addSYS, calcSYS=calcSYS, SV=False)
				#SF, pars, chi2 = runSF_x(file_name, pt_bins, m, WP, merge=merge, glue=glue, addSYS=addSYS, calcSYS=calcSYS)
			
			print "		", SF	
			print "		Time to run: ", np.round((time.time() - start)/60, 2), "min"

			eval("SF_dict_"+WP+"['SF_'+template.replace('0p5', 'down').replace('1p5', 'up')].append(float("+str(SF)+"))")
	return

# Calculate systematic of fitting all tempaltes separately regardeless of similar shapes
def step2(WP=WP):
	#glue=False; addSYS=True; calcSYS=False; 
	glue=True; addSYS=False; calcSYS=False; # just a quick test
	SFs = [] # temp
	for m, pt_bin in enumerate(pt_bins):
		start = time.time()
		print WP, pt_bin
		merge=False
		merge = when_to_merge(WP, pt_bins, m)
		if int(pt_bin.split("to")[-1]) <= 350:
			#file_name = r+name1+WP+root
			file_name = r+name1+root
			# LTSVIs350_dataUseMCJPcalib_SysMerged_SFtemplates_DeepDoubleXH.root
			SF, pars, chi2 = runSF_x(file_name, pt_bins, m, WP, merge=merge, glue=glue,  addSYS=addSYS, calcSYS=calcSYS, LTSV=False)
		else: 
			#file_name = r+name2+WP+root
			file_name = r+name2+root
			print file_name
			SF, pars, chi2 = runSF_x(file_name, pt_bins, m, WP, merge=merge, glue=glue,  addSYS=addSYS, calcSYS=calcSYS, LTSV=False)
		
		print "		", SF	
		print "      Time to run: ", np.round((time.time() - start)/60, 2), "min"
		#make_dirs('results5temp/'+str(WP)+str(pt_bin))
		#os.system('cp -r pics/* results5temp/'+str(WP)+str(pt_bin))

		eval("SF_dict_"+WP+"['SF_5_temp'].append(float("+str(SF)+"))")
		SFs.append(SF)
	return SFs

# Calculate MCJP Calib - different input files, data is modified
def step2_1(WP=WP):
	glue=True; addSYS=False;  calcSYS=False
	for m, pt_bin in enumerate(pt_bins):
		start = time.time()
		print WP, pt_bin
		merge=False
		merge = when_to_merge(WP, pt_bins, m)
		if int(pt_bin.split("to")[-1]) <= 350:
			file_name = rjp+"Run2017BCDEF_ReReco_QCDMuonEnriched_AK4DiJet170_Pt250_Final_DoubleMuonTaggedFatJets_histograms_btagval_allVars_ptReweighted_dataUseMCJPcalib_SysMerged_SFtemplates_"+WP+root
			SF, pars, chi2 = runSF_x(file_name, pt_bins, m, WP, merge=merge,glue=glue,  addSYS=addSYS, calcSYS=calcSYS)
			
		else:
			file_name = rjp+"Run2017BCDEF_ReReco_QCDMuonEnriched_AK8Jet300orAK4Jet300_Pt350_Final_DoubleMuonTaggedFatJets_histograms_btagval_allVars_ptReweighted_dataUseMCJPcalib_SysMerged_SFtemplates_"+WP+root
			SF, pars, chi2 = runSF_x(file_name, pt_bins, m, WP, merge=merge, glue=glue, addSYS=addSYS, calcSYS=calcSYS)
		
		print SF				
		print "Time to run: ", np.round((time.time() - start)/60, 2), "min"

		eval("SF_dict_"+WP+"['SF_JP'].append(float("+str(SF)+"))")

# Calculate actual SF values, include all systematics
def step3(WP=WP, SF_dict=SF_dict_empty):
	glue=True; calcSYS=True; addSYS=True
	# Make some arrays to catch outputs and log them
	SFs, sigma_stats, syst_ups, syst_downs = [], [], [], []
	errors_all, variances_all = [], []
	SF_SV, preeff, tageff = [], [], []
	for m, pt_bin in enumerate(pt_bins):
		start = time.time()
		print WP, m, pt_bin
		merge = False
		merge = when_to_merge(WP, pt_bins, m)
		#print merge
		if int(pt_bin.split("to")[-1]) <= 350:
			file_name = r1+name1+"_"+WP+root
			SF, sigma_stat, syst_up, syst_down, variances_names, errors, variances, nom_pars, chi2 = runSF_x(file_name, pt_bins, m, WP, merge=merge, glue=glue,  addSYS=addSYS, calcSYS=calcSYS, SF_dict=SF_dict) 
		else:
			file_name = r1+name2+"_"+WP+root
			SF, sigma_stat, syst_up, syst_down, variances_names, errors, variances, nom_pars, chi2 = runSF_x(file_name, pt_bins, m, WP, merge=merge, glue=glue,  addSYS=addSYS, calcSYS=calcSYS, SF_dict=SF_dict) 
		
		errors_all.append(errors)
		variances_all.append(variances)
		SFs.append(np.round(SF, 2))
		sigma_stats.append(np.round(sigma_stat, 3))
		syst_ups.append(np.round(syst_up, 3))
		syst_downs.append(np.round(syst_down, 3))
		SF_SV.append(nom_pars[5]/nom_pars[4])
		preeff.append(nom_pars[1]/nom_pars[0])
		tageff.append(nom_pars[3]/nom_pars[2])
		
		print "Time to run: ", np.round((time.time() - start)/60, 2), "min"
	return SFs, sigma_stats, syst_ups, syst_downs, variances_names, errors_all, variances_all, SF_SV, preeff, tageff



# Specify input files
r1 = 'May17single/'
rjp = 'May17JPsingle3/'
name1 = 'Run2017BCDEF_ReReco_QCDMuonEnriched_AK4DiJet170_Pt250_Final_DoubleMuonTaggedFatJets_histograms_btagval_allVars_ptReweighted_SysMerged_SFtemplates'
name2 = 'Run2017BCDEF_ReReco_QCDMuonEnriched_AK8Jet300orAK4Jet300_Pt350_Final_DoubleMuonTaggedFatJets_histograms_btagval_allVars_ptReweighted_SysMerged_SFtemplates'
WP = "BL"
root = '.root'
r = 'March19/'
#name1 = "LTSVIs250to350_dataUseMCJPcalib_SysMerged_SFtemplates_"
#name2 = "LTSVIs350_dataUseMCJPcalib_SysMerged_SFtemplates_"

#name2 = "collatedDDBvL"
#name2 = "collatedDDBvLv3"
#name2 = "collatedv3"
name2 = "collatedDoubleBsmall"


# Specify binning
#pt_bins = ['pt250to350', 'pt350to430', 'pt430to2000']
#pt_bins = ['pt250to300', 'pt300to350',  'pt350to450', 'pt450to2000']
#pt_bins = ['pt250to350', 'pt350to450', 'pt450to2000']
#pt_bins = ['pt250to350', 'pt350to2000']
#pt_bins = ['pt250to350','pt350to430', 'pt430to2000']
pt_bins = ['pt350to430', 'pt430to2000']

# Specify WPs to run
WPs = ['DoubleBL', 'DoubleBM1', 'DoubleBM2', 'DoubleBT']
#WPs = ['DDBvLL','DDBvLM1', 'DDBvLM2', 'DDBvLT1', 'DDBvLT2']
WP_dicts = [SF_dict_DoubleBL, SF_dict_DoubleBM1, SF_dict_DoubleBM2, SF_dict_DoubleBT]
#WP_dicts = [SF_dict_DDBvLL,SF_dict_DDBvLM1,SF_dict_DDBvLM2,SF_dict_DDBvLT1,SF_dict_DDBvLT2 ]

# Load intermediate results - steps: 1,2, 2_1
if args.load:
	with open('pkl1.pkl', 'r') as f:
	    SF_dict_DoubleBL = pickle.load(f)
	with open('pkl2.pkl', 'r') as f:
	    SF_dict_DoubleBM1 = pickle.load(f)
	with open('pkl3.pkl', 'r') as f:
	    SF_dict_DoubleBM2 = pickle.load(f)
	with open('pkl4.pkl', 'r') as f:
	    SF_dict_DoubleBT = pickle.load(f)

	WP_dicts = [SF_dict_DoubleBL, SF_dict_DoubleBM1, SF_dict_DoubleBM2, SF_dict_DoubleBT]

# If running from start - run all pre SFs
else:
	pass
	#for WP in WPs: M = step2(WP=WP)
	#sys.exit()
	#for WP in WPs: M = step1(WP=WP)
	#for WP in WPs: M = step2(WP=WP)
	#for WP in WPs: M = step2_1(WP=WP)


with open('pkl1.pkl', 'w') as f:
    pickle.dump(SF_dict_DoubleBL, f)
with open('pkl2.pkl', 'w') as f:
    pickle.dump(SF_dict_DoubleBM1, f)
with open('pkl3.pkl', 'w') as f:
    pickle.dump(SF_dict_DoubleBM2, f)
with open('pkl4.pkl', 'w') as f:
    pickle.dump(SF_dict_DoubleBT, f)


# Run final SF calculation
headers = ["Systematic"]
arrays = []
reprint = []
parinfo = []
for WP, WP_dict in zip(WPs, WP_dicts): 
	print WP
	#SFs, sigma_stats, syst_ups, syst_downs, variances_names, errors_allpt, variances_allpt, SF_SV, preeff, tageff = step3(WP=WP, SF_dict=WP_dict)	
	SFs, sigma_stats, syst_ups, syst_downs, SF_SV, preeff, tageff = [0]*len(pt_bins), [0]*len(pt_bins), [0]*len(pt_bins),[0]*len(pt_bins), [0]*len(pt_bins), [0]*len(pt_bins), [0]*len(pt_bins) # temp
	SFs = step2(WP=WP) # temp 
	

	#reprint.append([WP[len("DoubleB"):], SFs, sigma_stats, syst_ups, syst_downs, SF_SV, preeff, tageff])
	reprint.append([WP[len("DeepDoubleX"):], SFs, sigma_stats, syst_ups, syst_downs, SF_SV, preeff, tageff])
	for pt in pt_bins:
		#headers.append(WP[len("DoubleB"):]+pt)
		headers.append(WP[len("DeepDoubleX"):]+pt)

	#names = variances_names
	#for ar in errors_allpt:
	#	arrays.append(ar)

# Store outputs in data frames
results = []
for row in reprint:
	results.append(np.array(row[1:]))

#df = pd.DataFrame(np.array(arrays).transpose(), columns=headers[1:], index=names)
#print df
df2 = pd.DataFrame(np.concatenate((results), axis=1), columns=headers[1:], index=['SF', 'Stat', 'Sys up', 'Sys down', 'SF_SV', 'preeff', 'tageff'] )
df2.loc['Combined up'] = np.sqrt(df2.loc['Stat'].values**2 + df2.loc['Sys up'].values**2)
df2.loc['Combined down'] = np.sqrt(df2.loc['Stat'].values**2 + df2.loc['Sys down'].values**2)
print df2

#df.to_csv('DF_sys.csv')
df2.to_csv('DF_res.csv')
