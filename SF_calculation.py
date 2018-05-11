import os, sys, errno
import copy
import numpy as np
import subprocess
import pandas as pd
import time
#from test_alice import runSF_x
from get_SF_cfit import runSF_x
from pprint import pprint
import pickle
pd.set_option('display.max_columns', 500)
pd.set_option('display.width', 1000)

def make_dirs(dirname):
    """
    Ensure that a named directory exists; if it does not, attempt to create it.
    """
    try:
        os.makedirs(dirname)
    except OSError, e:
        if e.errno != errno.EEXIST:
            raise

import argparse
parser = argparse.ArgumentParser(description='Process some integers.')
parser.add_argument('--load', action='store_true', default=False, help="Use this if you've ran the first part, but crashed when computing the final systematics, it saves the intermediate steps.")
args = parser.parse_args()

# Names for 2018 Spring run, incomplete 2017 dataset
#r1 = '2018filesSyst/'
#r2 = '2018filesSyst/'
#JP_r = '2018filesSystJP/'
#name = 'Run2017BCDEF_ReReco_QCDMuonEnriched_AK8Jet300orAK4Jet300_Pt350_Final_DoubleMuonTaggedFatJets_histograms_btagval_v2_ptReweighted_SysMerged_SFtemplates'
#name2 = 'Run2017BCDEF_ReReco_QCDMuonEnriched_AK4DiJet170_Pt250_Final_DoubleMuonTaggedFatJets_histograms_btagval_v1_v2_ptReweighted_SysMerged_SFtemplates'
# Names for2018 Spring run, complete 2017 dataset (almost)
#r1 = 'Mar15-2018/SFtemplates/'
#r2 = 'Mar15-2018/SFtemplates/'
#P_r = 'Mar15-2018/SFtemplates_dataUseMCJPCalib/'
#name = 'Run2017BCDEF_ReReco_QCDMuonEnriched_AK8Jet300orAK4Jet300_Pt350_Final_DoubleMuonTaggedFatJets_histograms_btagval_v3_ptReweighted_SysMerged_SFtemplates'
#name2 = 'Run2017BCDEF_ReReco_QCDMuonEnriched_AK4DiJet170_Pt250_Final_DoubleMuonTaggedFatJets_histograms_btagval_v1_v3_ptReweighted_SysMerged_SFtemplates'
# Names for SV try
#r1 = 'normfile/'
#r1 = 'singlefilenorm/'
#r1 = 'April27single/'
r1 = 'May9single/'
rjp = 'May9JPsingle/'
#r1 = 'April27forcenorm/'
name1 = 'Run2017BCDEF_ReReco_QCDMuonEnriched_AK4DiJet170_Pt250_Final_DoubleMuonTaggedFatJets_histograms_btagval_allVars_ptReweighted_SysMerged_SFtemplates'
name2 = 'Run2017BCDEF_ReReco_QCDMuonEnriched_AK8Jet300orAK4Jet300_Pt350_Final_DoubleMuonTaggedFatJets_histograms_btagval_allVars_ptReweighted_SysMerged_SFtemplates'

WP = "BL"
root = '.root'

templates = ['b_0p5', 'b_1p5', 'c_0p5', 'c_1p5', 'cfromg_0p5', 'cfromg_1p5', 'l_0p5', 'l_1p5']
pt_bins = ['pt250to350', 'pt350to430', 'pt430to2000']

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
	'SF_JP' : []
}

SF_dict_DoubleBL = copy.deepcopy(SF_dict_empty)
SF_dict_DoubleBM1 = copy.deepcopy(SF_dict_empty)
SF_dict_DoubleBM2 = copy.deepcopy(SF_dict_empty)
SF_dict_DoubleBH = copy.deepcopy(SF_dict_empty)

def step1(templates=templates, WP=WP):
	M = []
	glue=True; inclSYS=False
	print WP
	for n, template in enumerate(templates):		
		for m, pt_bin in enumerate(pt_bins):
			start = time.time()
			print WP, pt_bin, template
			merge=False
			if WP == 'DoubleBM2' and m == 1: merge=True
			if WP == 'DoubleBH' and m in [0,1]: merge=True
			if m == 0: 
				file_name = r1+name1+"_"+WP+"_"+template+root
				SF, pars = runSF_x(file_name, pt_bin, WP, merge=merge, glue=glue, inclSYS=inclSYS)
			else:
				file_name = r1+name2+"_"+WP+"_"+template+root			
				SF, pars = runSF_x(file_name, pt_bin, WP, merge=merge, glue=glue, inclSYS=inclSYS)
			
			print "		", SF	
			print "		Time to run: ", np.round((time.time() - start)/60, 2), "min"

			eval("SF_dict_"+WP+"['SF_'+template.replace('0p5', 'down').replace('1p5', 'up')].append(float("+str(SF)+"))")
	return M

def step2(WP=WP):
	bin_pars = []
	bin_nams = []
	glue=False;  inclSYS=False
	M = []
	print WP
	for m, pt_bin in enumerate(pt_bins):
		start = time.time()
		print pt_bin
		merge=False
		if WP == 'DoubleBM2' and m == 1: merge=True
		if WP == 'DoubleBH' and m in [0,1]: merge=True
		if merge == True :print "MERGE"
		if m == 0: 
			file_name = r1+name1+"_"+WP+root
			SF, pars = runSF_x(file_name, pt_bin, WP, merge=merge, glue=glue, inclSYS=inclSYS)		
		else: 
			file_name = r1+name2+"_"+WP+root
			SF, pars = runSF_x(file_name, pt_bin, WP, merge=merge, glue=glue, inclSYS=inclSYS)
		
		bin_nams.append(WP+pt_bin)
		bin_pars.append(pars)
		print "Time to run: ", np.round((time.time() - start)/60, 2), "min"
		print WP, pt_bin
		make_dirs('results5temp/'+str(WP)+str(pt_bin))
		os.system('cp -r pics/* results5temp/'+str(WP)+str(pt_bin))

		eval("SF_dict_"+WP+"['SF_5_temp'].append(float("+str(SF)+"))")

	df = pd.DataFrame(data=bin_pars, columns=["parJP", "parJP_tag", "parJPpass", "parJPpass_tag", "parSV", "parSV_tag"], index=bin_nams)
	print df


def step2_1(WP=WP):
	glue=False;  inclSYS=False
	M = []
	print WP
	for m, pt_bin in enumerate(pt_bins):
		start = time.time()
		print pt_bin
		merge=False
		if WP == 'DoubleBM2' and m == 1: merge=True
		if WP == 'DoubleBH' and m in [0,1]: merge=True
		if merge == True :print "MERGE"
		if m == 0: 
			file_name = rjp+"Run2017BCDEF_ReReco_QCDMuonEnriched_AK4DiJet170_Pt250_Final_DoubleMuonTaggedFatJets_histograms_btagval_allVars_ptReweighted_dataUseMCJPcalib_SysMerged_SFtemplates_"+WP+root
			SF, pars = runSF_x(file_name, pt_bin, WP, merge=merge,glue=glue, inclSYS=inclSYS)
			
		else:
			file_name = rjp+"Run2017BCDEF_ReReco_QCDMuonEnriched_AK8Jet300orAK4Jet300_Pt350_Final_DoubleMuonTaggedFatJets_histograms_btagval_allVars_ptReweighted_dataUseMCJPcalib_SysMerged_SFtemplates_"+WP+root
			SF, pars = runSF_x(file_name, pt_bin, WP, merge=merge, glue=glue, inclSYS=inclSYS)
		
		print SF				
		print "Time to run: ", np.round((time.time() - start)/60, 2), "min"

		eval("SF_dict_"+WP+"['SF_JP'].append(float("+str(SF)+"))")


def step3(WP=WP, SF_dict=SF_dict_empty):
	glue=True;  inclSYS=True
	SFs, sigma_stats, syst_ups, syst_downs = [], [], [], []
	errors_all, variances_all = [], []
	for m, pt_bin in enumerate(pt_bins):
		start = time.time()
		print m, pt_bin
		merge = False
		if WP == 'DoubleBM2' and m == 1: merge=True
		if WP == 'DoubleBH' and m in [0,1]: merge=True
		if merge == True :print "MERGE"
		if m == 0: 
			#file_name = r2+name2+"_"+WP+root
			file_name = r1+name1+"_"+WP+root
			print file_name		
			SF, sigma_stat, syst_up, syst_down, variances_names, errors, variances = runSF_x(file_name, pt_bin, WP, merge=merge, glue=glue, inclSYS=inclSYS, SF_dict=SF_dict) 
		else:
			file_name = r1+name2+"_"+WP+root
			print file_name
			SF, sigma_stat, syst_up, syst_down, variances_names, errors, variances = runSF_x(file_name, pt_bin, WP, merge=merge, glue=glue, inclSYS=inclSYS, SF_dict=SF_dict) 
		
		errors_all.append(errors)
		variances_all.append(variances)
		SFs.append(np.round(SF, 2))
		sigma_stats.append(np.round(sigma_stat, 3))
		syst_ups.append(np.round(syst_up, 3))
		syst_downs.append(np.round(syst_down, 3))
		
		#make_dirs('results/'+str(WP)+str(pt_bin))
		#os.system('cp -r pics/* results/'+str(WP)+str(pt_bin))
		print "Time to run: ", np.round((time.time() - start)/60, 2), "min"
	return SFs, sigma_stats, syst_ups, syst_downs, variances_names, errors_all, variances_all

WPs = ['DoubleBL', 'DoubleBM1', 'DoubleBM2', 'DoubleBH']
WP_dicts = [SF_dict_DoubleBL, SF_dict_DoubleBM1, SF_dict_DoubleBM2, SF_dict_DoubleBH]

import shelve
if args.load:
	WP_dicts = []
	shelf = shelve.open("temp.shlf", flag='r')
	for key in shelf.keys():
	    WP_dicts.append(shelf[key])
	shelf.close()

else:
	#pass
	for WP in WPs: M = step1(WP=WP)
	for WP in WPs: M = step2(WP=WP)
	for WP in WPs: M = step2_1(WP=WP)

for i in range(len(WPs)):
	pprint(WP_dicts[i])


shelf = shelve.open("temp.shlf", flag="c")
for i in range(len(WPs)):
	shelf[WPs[i]] = WP_dicts[i]
shelf.close()


print WP_dicts

bins = ["pt1", "pt2", "pt3"]
headers = ["Systematic"]
arrays = []
for WP, WP_dict in zip(WPs, WP_dicts): 
	#if WP != "DoubleBL": continue
	SFs, sigma_stats, syst_ups, syst_downs, variances_names, errors_allpt, variances_allpt = step3(WP=WP, SF_dict=WP_dict)	
	print [WP[len("DoubleB"):], SFs, sigma_stats, syst_ups, syst_downs]
	
	for pt in bins:
		headers.append(WP[len("DoubleB"):]+pt)

	names = variances_names
	for ar in errors_allpt:
		arrays.append(ar)
		print ar
		print len(names), len(ar)
import csv
from prettytable import PrettyTable
f = open('Test.csv', 'w')
writer = csv.writer(f, delimiter=',', quotechar='|', quoting=csv.QUOTE_MINIMAL, lineterminator='\n')
writer.writerow(headers)

print headers
t = PrettyTable(headers)
for i in range(len(arrays[0])):
	row = [names[i]]
	csvrow = [names[i]]
	for errors in arrays:
		val = np.round(errors[i],4)
		row.append("{:.2E}".format(val))
		csvrow.append(val)
	writer.writerow(csvrow)
	t.add_row(row)
print t
f.close()

table_txt = t.get_string()
with open('Test.txt','w') as file:
    file.write(table_txt)
