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


r1 = 'May17single/'
rjp = 'May17JPsingle3/'
name1 = 'Run2017BCDEF_ReReco_QCDMuonEnriched_AK4DiJet170_Pt250_Final_DoubleMuonTaggedFatJets_histograms_btagval_allVars_ptReweighted_SysMerged_SFtemplates'
name2 = 'Run2017BCDEF_ReReco_QCDMuonEnriched_AK8Jet300orAK4Jet300_Pt350_Final_DoubleMuonTaggedFatJets_histograms_btagval_allVars_ptReweighted_SysMerged_SFtemplates'

WP = "BL"
root = '.root'

templates = ['b_0p5', 'b_1p5', 'c_0p5', 'c_1p5', 'cfromg_0p5', 'cfromg_1p5', 'l_0p5', 'l_1p5']
#pt_bins = ['pt250to350', 'pt350to430', 'pt430to2000']
pt_bins = ['pt250to300', 'pt300to350',  'pt350to450', 'pt450to2000']
#pt_bins = ['pt250to350', 'pt350to450', 'pt450to2000']
#pt_bins = ['pt250to350', 'pt350to2000']

def step0(pt_bins, WP=WP, sys=None, SV=True):
	glue=True; calcSYS=False; addSYS=False
	SFs = []
	chi2sS = []
	for m, pt_bin in enumerate(pt_bins):
		start = time.time()
		print WP, m, pt_bin
		merge = False
		if WP == 'DoubleBM2' and m in [1]: merge=2
		if WP == 'DoubleBM2' and m in [2]: merge=True
		if WP == 'DoubleBH' and m in [1]: merge=2
		if WP == 'DoubleBH' and m in [0,2]: merge=True
		if m in [0,1]: 
			file_name = r1+name1+"_"+WP+root
			SF, pars, chi2s = runSF_x(file_name, pt_bins, m, WP, merge=merge, glue=glue, addSYS=addSYS, calcSYS=calcSYS, systname=sys, SV=SV)
		else: 
			file_name = r1+name2+"_"+WP+root
			SF, pars, chi2s = runSF_x(file_name, pt_bins, m, WP, merge=merge, glue=glue, addSYS=addSYS, calcSYS=calcSYS, systname=sys, SV=SV)
		#SF = 1
		#chi2s = [1,1,1,1,1,1]
		chi2sS.append(chi2s)
		SFs.append(SF)
		print "Time to run: ", np.round((time.time() - start)/60, 2), "min"
	return SFs, chi2sS

WPs = ['DoubleBL', 'DoubleBM1', 'DoubleBM2', 'DoubleBH']
"""
# Quick Checks
allSFs = []
for WP in WPs: 
	SFsWP, chi2s = step0(pt_bins,WP=WP)
	allSFs.append(SFsWP)
df = pd.DataFrame(allSFs)
print df


# Comparing syst influence JP only
allSFs = []
systs = ["JES", "NTRACKS", "BFRAG", "CFRAG", "CD", "K0L", "PU"]
for sys in systs:
	row = []
	headers = []
	for WP in WPs: 
		SFsWP, chi2s = step0(pt_bins, 	WP=WP, sys=sys, SV=False)
		row += SFsWP
		for m, pt in enumerate(pt_bins):
			headers.append(WP[len("DoubleB"):]+'pt'+str(m))
	allSFs.append(row)

df = pd.DataFrame(allSFs, index=systs, columns=headers)
print df
df.to_csv('DF_JPsyseffects.csv')
"""

import ROOT
def count(pt_bins, wp=WP):
	count_along_wp = []
	for m, pt in enumerate(pt_bins):
		print WP, m, pt
		if m in [0,1]: 
			file_name = r1+name1+"_"+WP+root
		else: 
			file_name = r1+name2+"_"+WP+root
		F = ROOT.TFile.Open(file_name)
		#keys = F.GetListOfKeys()
		#names = [key.GetName() for key in keys]
		var = "tau1VertexMassCorr"
		hname_data = "UNWEIGHTED__DATA__FatJet_"+var+"_all_"+pt+"_data_opt"
		hname_data_tag = "UNWEIGHTED__DATA__FatJet_"+var+"_"+wp+"pass_"+pt+"_data_opt"
		hname_data_untag = "UNWEIGHTED__DATA__FatJet_"+var+"_"+wp+"fail_"+pt+"_data_opt"
		
		var = "JP"
		jhname_data = "UNWEIGHTED__DATA__FatJet_"+var+"_all_"+pt+"_data_opt"
		jhname_data_tag = "UNWEIGHTED__DATA__FatJet_"+var+"hasSV"+"_all_"+pt+"_data_opt"
		jhname_data_untag = "UNWEIGHTED__DATA__FatJet_"+var+"noSV"+"_all_"+pt+"_data_opt"
		jthname_data = "UNWEIGHTED__DATA__FatJet_"+var+"_"+wp+"pass_"+pt+"_data_opt"
		jthname_data_tag = "UNWEIGHTED__DATA__FatJet_"+var+"hasSV_"+wp+"pass_"+pt+"_data_opt"
		jthname_data_untag = "UNWEIGHTED__DATA__FatJet_"+var+"noSV_"+wp+"pass_"+pt+"_data_opt"

		a = F.Get(hname_data).Integral()
		b = F.Get(hname_data_tag).Integral()
		c = F.Get(jhname_data).Integral()
		d = F.Get(jhname_data_tag).Integral()
		e = F.Get(jthname_data).Integral()
		f = F.Get(jthname_data_tag).Integral()
		#count_along_wp.append([a,b,c,d,e,f])
		count_along_wp.append([c,d,e,f,a,b])

	return count_along_wp

allNs = []
binnings = [['pt250to300', 'pt300to350', 'pt350to450', 'pt450to2000']]#, ['pt250to350', 'pt350to450', 'pt450to2000'],['pt250to350', 'pt350to2000']]
headers = []
for pt_bins in binnings:
	row = []
	for WP in WPs: 
		count_along_wp = count(pt_bins, wp=WP)
		row += count_along_wp
		for m, pt in enumerate(pt_bins):
			headers.append(WP[len("DoubleB"):].replace("H", "T")+str(pt))
	allNs += row

print allNs
df3 = pd.DataFrame(np.array(allNs).transpose(), index=["JP", "JP_hasSV", "JPpass", "JPpass_hasSV", "SV", "SV_tag"], columns=headers)
print df3
df3.to_csv('DF_counts4bins.csv')

# Comparing syst influence
"""
allSFs = []
systs = ["JES", "NTRACKS", "BFRAG", "CFRAG", "CD", "K0L", "PU"]
for sys in systs:
	row = []
	headers = []
	for WP in WPs: 
		SFsWP, chi2s = step0(pt_bins, 	WP=WP, sys=sys)
		row += SFsWP
		for m, pt in enumerate(pt_bins):
			headers.append(WP[len("DoubleB"):]+'pt'+str(m))
	allSFs.append(row)

df = pd.DataFrame(allSFs, index=systs, columns=headers)
print df
df.to_csv('DF_syseffects.csv')


# Binning chi2s
allChis = []
#binnings = [['pt250to350', 'pt350to430', 'pt430to2000'], ['pt250to350', 'pt350to450', 'pt450to2000'],['pt250to350', 'pt350to2000']]
binnings = [['pt250to300', 'pt300to350', 'pt350to450', 'pt450to2000']]#, ['pt250to350', 'pt350to450', 'pt450to2000'],['pt250to350', 'pt350to2000']]
headers = []
for pt_bins in binnings:
	row = []
	for WP in WPs: 
		SFsWP, chi2s = step0(pt_bins,WP=WP)
		print chi2s
		row += chi2s
		for m, pt in enumerate(pt_bins):
			headers.append(WP[len("DoubleB"):].replace("H", "T")+str(pt))
	allChis += row

df2 = pd.DataFrame(np.array(allChis).transpose(), index=["JP", "JP_hasSV", "JPpass", "JPpass_hasSV", "SV", "SV_tag"], columns=headers)
print df2
df2.to_csv('DF_chi2s4bins.csv')

"""