#from __future__ import unicode_literals
from cfit import *
import numpy as np

def runSF_x(file, pt, tag, var = "JP", merge=False, glue=True, inclSYS=True, SF_dict={}):
	N_digits_SF = 3
	N_digits_err = 4
	if pt == "pt250to350": ptbin = 0
	if pt == "pt350to430": ptbin = 1
	if pt == "pt430to2000": ptbin = 2
	cf = cfit(var + " discriminator")
	cf.SetVerbose(0)
	cf.ProducePlots(True)
	cf.SetLegendHeader(pt+" "+tag)
	cf.SetOptimization(OPT_MORPH_SGN_SIGMA)
	cf.SetCovarianceMode(COV_MAX)
	cf.SetMorphing(OPTMORPH_CUTOFF,0.5)
	cf.SetInputFile(file)

	if inclSYS:
		#if merge:
		#systlist = ["JES", "NTRACKS", "BFRAG", "CFRAG", "CD"]
		#else:
		systlist = ["JES", "NTRACKS", "BFRAG", "CFRAG", "CD", "K0L", "PU"]
		for sysName in systlist:
			# Flipped down up in Alice's code
			cf.AddSys(sysName, "_"+sysName+"down" ,"_"+sysName+"up")
			#cf.AddSys(sysName, "_"+sysName+"up" ,"_"+sysName+"down") #-- this one is correct

	cf.SetMatrixName("matrix_"+pt)
	cf.SetMatrixOption("WRITE")

	# Set DATA names
	hname_data = "UNWEIGHTED__DATA__FatJet_"+var+"_all_"+pt+"_data_opt"
	hname_data_tag = "UNWEIGHTED__DATA__FatJet_"+var+"_"+tag+"pass_"+pt+"_data_opt"
	hname_data_untag = "UNWEIGHTED__DATA__FatJet_"+var+"_"+tag+"fail_"+pt+"_data_opt"

	cf.SetData(hname_data)   
	cf.SetDataTag(hname_data_tag)
	cf.SetDataUntag(hname_data_untag)	

	# Hist name setup
	def_name_qcd = "UNWEIGHTED__QCDMuEnr__FatJet_"
	hname_bfromg = def_name_qcd+var+"_all_"+pt+"_bfromg_opt"
	hname_b = def_name_qcd+var+"_all_"+pt+"_b_opt"
	hname_cfromg = def_name_qcd+var+"_all_"+pt+"_cfromg_opt"
	hname_c = def_name_qcd+var+"_all_"+pt+"_c_opt"
	hname_l = def_name_qcd+var+"_all_"+pt+"_l_opt"
	hname_b_cfromg = def_name_qcd+var+"_all_"+pt+"_b_cfromg_opt"
	hname_c_l = def_name_qcd+var+"_all_"+pt+"_c_l_opt"

	hname_bfromg_tag = hname_bfromg.replace("_all_", "_"+tag+"pass_")
	hname_b_tag = hname_b.replace("_all_", "_"+tag+"pass_")
	hname_cfromg_tag = hname_cfromg.replace("_all_", "_"+tag+"pass_")
	hname_c_tag = hname_c.replace("_all_", "_"+tag+"pass_")
	hname_l_tag = hname_l.replace("_all_", "_"+tag+"pass_")
	hname_b_cfromg_tag = hname_b_cfromg.replace("_all_", "_"+tag+"pass_")
	hname_c_l_tag = hname_c_l.replace("_all_", "_"+tag+"pass_")

	hname_bfromg_untag = hname_bfromg.replace("_all_", "_"+tag+"fail_")
	hname_b_untag = hname_b.replace("_all_", "_"+tag+"fail_")
	hname_cfromg_untag = hname_cfromg.replace("_all_", "_"+tag+"fail_")
	hname_c_untag = hname_c.replace("_all_", "_"+tag+"fail_")
	hname_l_untag = hname_l.replace("_all_", "_"+tag+"fail_")
	hname_b_cfromg_untag = hname_b_cfromg.replace("_all_", "_"+tag+"fail_")
	hname_c_l_untag = hname_c_l.replace("_all_", "_"+tag+"fail_")

	# Mergning for low count templates:
	if merge: tempNs = ["g #rightarrow b#bar{b}", "b + g #rightarrow c#bar{c}", "c + dusg"]
	else: 	tempNs = ["g #rightarrow b#bar{b}", "b", "g #rightarrow c#bar{c}", "c", "dusg"]
	def add_templates(glue, merge):
		# AddTemplate(label, name in input file, color)
		if merge:			
			cf.AddTemplate(tempNs[0],hname_bfromg,	65)
			cf.AddTemplate(tempNs[1],hname_b_cfromg,628)
			cf.AddTemplate(tempNs[2],hname_c_l,	597)
			#if glue:
			#	cf.GlueTemplates(tempNs[1:],"other flavours",28);
		else:
			cf.AddTemplate(tempNs[0],hname_bfromg,	65)
			cf.AddTemplate(tempNs[1],hname_b,		213)
			cf.AddTemplate(tempNs[2],hname_cfromg,	208)
			cf.AddTemplate(tempNs[3],hname_c, 		206)
			cf.AddTemplate(tempNs[4],hname_l, 		212)
			if glue:
				cf.GlueTemplates([tempNs[1], tempNs[2]],"b + g #rightarrow c#bar{c}",628);
				cf.GlueTemplates([tempNs[3], tempNs[4]],"c + dusg",597)

		if merge:
			cf.AddTemplateTag(tempNs[0],hname_bfromg_tag,	65)
			cf.AddTemplateTag(tempNs[1],hname_b_cfromg_tag,628)
			cf.AddTemplateTag(tempNs[2],hname_c_l_tag,	597)
			#if glue:
			#	cf.GlueTemplatesTag(tempNs[1:],"other flavours",28);
		else:
			cf.AddTemplateTag(tempNs[0],hname_bfromg_tag,	65)
			cf.AddTemplateTag(tempNs[1],hname_b_tag,		213)
			cf.AddTemplateTag(tempNs[2],hname_cfromg_tag,	208)
			cf.AddTemplateTag(tempNs[3],hname_c_tag,		206)
			cf.AddTemplateTag(tempNs[4],hname_l_tag,		212)
			if glue:
				cf.GlueTemplatesTag(tempNs[1:],"other flavours",28);

		if merge:
			cf.AddTemplateUntag(tempNs[0],hname_bfromg_untag,	65)
			cf.AddTemplateUntag(tempNs[1],hname_b_cfromg_untag,628)
			cf.AddTemplateUntag(tempNs[2],hname_c_l_untag,	597)
		else:			
			cf.AddTemplateUntag(tempNs[0],hname_bfromg_untag,	7)
			cf.AddTemplateUntag(tempNs[1],hname_b_untag,		213)
			cf.AddTemplateUntag(tempNs[2],hname_cfromg_untag,	42)
			cf.AddTemplateUntag(tempNs[3],hname_c_untag,		8)
			cf.AddTemplateUntag(tempNs[4],hname_l_untag,		4)
		
	add_templates(glue, merge)

	def getSF(cf, tempNs, sysVar="", statVar=0):
		par, err = [], []
		par_tag, err_tag = [], []

		cf.SetSysVariation(sysVar)
		if statVar > 0: cf.SetStatVariation(statVar)
		cf.Run() # Run untag
		for i in range(cf.GetNPar()):
			par.append(cf.GetPar(i))
			err.append(cf.GetParErr(i))

		chi2 = cf.GetChisq()
		ndata = cf.GetNData()
		nmc1 = cf.GetNTemplate(tempNs[0]) #Signal
		nmc = 0. #Total
		for i, name in enumerate(tempNs): nmc += cf.GetNTemplate(name)
		
		cf.SetSysVariation(sysVar)
		if statVar > 0: cf.SetStatVariation(statVar)
		cf.Run("tag")  # Run tag
		for i in range(cf.GetNPar()):
			par_tag.append(cf.GetPar(i))
			err_tag.append(cf.GetParErr(i))

		chi2_tag = cf.GetChisq();
		ndata_tag = cf.GetNData();
		nmc1_tag = cf.GetNTemplate(tempNs[0]); #Signal
		nmc_tag = 0. #Total
		for i, name in enumerate(tempNs): nmc_tag += cf.GetNTemplate(name)

		fr = nmc1/nmc;
		fr_tag = nmc1_tag/nmc_tag;	   
		effMC = nmc1_tag/nmc1
		#effDATA = nmc_tag/nmc*par_tag[0]/par[0]*fr_tag/fr		
		#effDATA =  nmc1_tag/nmc1 * par_tag[0]/par[0]
		effDATA = effMC*par_tag[0]/par[0]

		sf = effDATA/effMC
		return sf

	#if inclSYS==False:
	SF = getSF(cf, tempNs, sysVar="")
	#print "SF =", SF

	if inclSYS==True:
		print "Calculating Errors"
		sigma_stat = 0
   		sigma_syst_up = 0
		sigma_syst_down = 0
		variances_names = []
		errors = []
		
		#Calculate systematics 1	
		ups, downs = ["_"+f+"up" for f in systlist], ["_"+f+"down" for f in systlist]
		systVarlist = list(sum(zip(ups, downs+[0]), ()))
		if len(systVarlist) != len(ups+downs): 
			print "Incorrect number of systs considered"
		
		cf.SetMatrixOption("READ")
		cf.ProducePlots(False)
		for i, sysVar in enumerate(systVarlist):
			SF_sys_i = getSF(cf, tempNs, sysVar=sysVar)
			print "SF", sysVar, SF, SF_sys_i, (SF - SF_sys_i)**2
			delta = SF - SF_sys_i
			sigma2 = delta**2
			errors.append(delta)
			variances_names.append(sysVar[1:])
			if delta < 0:
				sigma_syst_up += sigma2
			else:
				sigma_syst_down += sigma2

		# Calculate statistical errs
		cf.SetMatrixOption("READ")
		cf.ProducePlots(False)
		stat_SFs = []
		for i in range(100):
			SF_stat_i = getSF(cf, tempNs, sysVar="", statVar=667+i)
			stat_SFs.append(SF_stat_i)

		sf_stat_mean = sum(stat_SFs)/100.
		for SF_i in stat_SFs:
			sigma_stat += (sf_stat_mean - SF_i)**2/100.
		sigma_stat = np.sqrt(sigma_stat);

		# 100 -105 systematics missing
		for i in range(100,105):
			cf.SetMatrixOption("WRITE")
			cf.ProducePlots(False)
		
			if i == 100: cf.SetOptimization(OPT_NOCORR);
			if i == 101: cf.SetMorphing(OPTMORPH_CUTOFF,0.25);
			if i == 102: cf.SetMorphing(OPTMORPH_CUTOFF,0.75);
			if i == 103: cf.SetMorphing(OPTMORPH_GEOMETRIC);
			if i == 104:
				cf = cfit("")
				cf.SetVerbose(0)
				cf.ProducePlots(False)
				cf.SetOptimization(OPT_MORPH_SGN_SIGMA)
				cf.SetCovarianceMode(COV_MAX)
				cf.SetMorphing(OPTMORPH_CUTOFF,0.5)
			    			     
				cf.SetInputFile(file)
		   
				mname = "matrix_"+pt;
				cf.SetMatrixName(mname)
				cf.SetMatrixOption("WRITE")

				cf.SetData(hname_data)   
				cf.SetDataTag(hname_data_tag)
				cf.SetDataUntag(hname_data_untag)

				add_templates(glue, merge)

			SF_sys2_i = getSF(cf, tempNs, sysVar="", statVar=-1)
			cf.SetOptimization(OPT_MORPH_SGN_SIGMA)
			cf.SetMorphing(OPTMORPH_CUTOFF,0.5)

			delta = SF - SF_sys2_i
			sigma2 = delta**2
			errors.append(delta)
			variances_names.append("Morph"+str(i))
			print SF_sys2_i, delta
			sigma_syst_up += sigma2
			sigma_syst_down += sigma2

		# Calculate systematics 2
		"""for SF_i in SF_list:
			if SF_i == []: continue
			if SF < SF_i[ptbin]:
				sigma_syst_up += (SF - SF_i[ptbin])*(SF - SF_i[ptbin])
			else:
				sigma_syst_down += (SF - SF_i[ptbin])*(SF - SF_i[ptbin])
		"""
		for SF_i in SF_dict.keys():
			if SF_dict[SF_i] == []: continue
			error = SF - SF_dict[SF_i][ptbin]
			sigma2 = error**2
			if SF < SF_dict[SF_i][ptbin]:
				sigma_syst_up += sigma2
			else:
				sigma_syst_down += sigma2
			errors.append(error)
			variances_names.append(SF_i)

		syst_up = np.sqrt(sigma_syst_up)
		syst_down = np.sqrt(sigma_syst_down)

		print "Err_stat =", sigma_stat
		print "Err_syst_up =", syst_up
		print "Err_syst_down =", syst_down
		#print variances_names
		#errors = [np.sqrt(var) for var in variances]
		variances = [err**2 for err in errors]
		#print errors

	if inclSYS:
		print np.round(sigma_stat, N_digits_err), np.round(syst_up, N_digits_err), np.round(syst_down, N_digits_err)
		#return SF, np.round(sigma_stat, N_digits_err), np.round(syst_up, N_digits_err), np.round(syst_down, N_digits_err)
		return SF, sigma_stat, syst_up, syst_down, variances_names, errors, variances
	else:
		return SF

if __name__ == "__main__":
	fname = "/eos/user/a/anovak/SF/CMSSW_7_6_4/src/CFIT/SFcalculation/2018filesSyst/Run2017BCDEF_ReReco_QCDMuonEnriched_AK8Jet300orAK4Jet300_Pt350_Final_DoubleMuonTaggedFatJets_histograms_btagval_v2_ptReweighted_SysMerged_SFtemplates_DoubleBL.root"

	fname = "/eos/user/a/anovak/SF/CMSSW_7_6_4/src/CFIT/SFcalculation/2018filesSyst/Run2017BCDEF_ReReco_QCDMuonEnriched_AK4DiJet170_Pt250_Final_DoubleMuonTaggedFatJets_histograms_btagval_v1_v2_ptReweighted_SysMerged_SFtemplates_DoubleBL.root"

	fname2 = fname.replace(".root", "_b_0p5.root")

	SF_dict = {
	'SF_b_down' : [1.02, 0.98, 1.0],
	'SF_b_up' : [1.02, 0.99, 1.0], 
	'SF_c_down' : [1.02, 0.98, 1.0], 
	'SF_c_up' : [1.02, 0.99, 1.01], 
	'SF_cfromg_down' : [1.03, 1.01, 1.01], 
	'SF_cfromg_up' : [1.02, 0.98, 1.0], 
	'SF_l_down' : [1.02, 0.98, 1.0], 
	'SF_l_up' : [1.02, 0.99, 1.0], 
	'SF_5_temp' : [1.0, 1.06, 0.98], 
	'SF_JP' : []
	}

	tag = "DoubleBL"
	#tag = "DoubleBM2"
	pt = "pt250to350"
	#pt = "pt350to430"
	#pt = "pt430to2000"

	SF, sigma_stat, syst_up, syst_down, variances_names, errors, variances = runSF_x(fname, pt, tag, glue=True, inclSYS=True, SF_dict=SF_dict) 
	#SF = runSF_x(fname, pt, tag, glue=False, inclSYS=False, SF_dict=SF_dict) 

	from prettytable import PrettyTable
	t = PrettyTable(["WP+pt", "Sigma", "Sigma^2"])
	variances = ["{:.2E}".format(var) for var in variances]
	errors = ["{:.2E}".format(err) for err in errors]
	for name, var, err in zip(variances_names, errors, variances):
		t.add_row([name, err, var])
	print t