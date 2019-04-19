#from __future__ import unicode_literals
from cfit import *
import ROOT as r
import numpy as np
import os

def runSF_x(file, bins, pt_bin, wp, 
			glue=True, #glue to bind separate templates and vary them together
			addSYS=False, calcSYS=True, SF_dict={},
			ccSignal = False,
			LTSV=True, # Use SV mass method instead of JP only
			systlist = ["BFRAG", "CFRAG","K0L", "PU"],
			debug=False,
			systname=None): #Add only one template (for debugging)

	stat_n = 100 #100 #Set how many time to run statistical variation
	pt = bins[pt_bin] # Pick a pt bin

	cfJP = cfit("JP discriminator")
	cfSV = cfit("SV discriminator")
	cfJPtag = cfit("JPtagged discriminator")

	cfJP.SetMatrixName("matrices/matrix_"+pt+'_'+wp)
	cfJP.SetMatrixOption("WRITE")
	cfSV.SetMatrixName("matrices/SVmatrix_"+pt+'_'+wp)
	cfSV.SetMatrixOption("WRITE")
	cfJPtag.SetMatrixName("matrices/JPtagmatrix_"+pt+'_'+wp)
	cfJPtag.SetMatrixOption("WRITE")

	# Settings
	for cf in [cfJP, cfSV, cfJPtag]:
		cf.SetVerbose(0)
		cf.ProducePlots(True)
		make_dirs('pics/'+wp+'/'+cf.GetRun().split(" ")[0])
		cf.SetPicsDir('pics/'+wp+'/'+cf.GetRun().split(" ")[0])
		cf.SetBatch(True)
		cf.SetLegendHeader(pt+" "+wp)
		cf.SetOptimization(OPT_MORPH_SGN_SIGMA)
		cf.SetCovarianceMode(COV_MAX)
		cf.SetMorphing(OPTMORPH_CUTOFF,0.5)
		cf.SetInputFile(file)
	

	#if inclSYS:
	if addSYS:
		if systname != None:
			systlist = [systname]
		else:
			#systlist = ["JES", "NTRACKS", "BFRAG", "CFRAG", "CD", "K0L", "PU"]
			systlist = systlist
		print "Added systs", systlist
		for sysName in systlist:
			for cf in [cfJP, cfSV, cfJPtag]:
				cf.AddSys(sysName, "_"+sysName+"up" ,"_"+sysName+"down")
			

	# Mergning for low count templates:
	# if merge == 1: tempNs = ["g #rightarrow b#bar{b}", "b + g #rightarrow c#bar{c}", "c + dusg"]
	# elif merge == 2: tempNs = ["g #rightarrow b#bar{b}", "b + g #rightarrow c#bar{c} + c + dusg"]
	# else: 	tempNs = ["g #rightarrow b#bar{b}", "b", "g #rightarrow c#bar{c}", "c", "dusg"]
	def add_templates(cf, glue, var = "JP", LTSV=LTSV, tag=False, ccSignal=False):
		# Set DATA names
		hname_data = "UNWEIGHTED__DATA__FatJet_"+var+"_all_"+pt+"_data_opt"
		hname_data_tag = "UNWEIGHTED__DATA__FatJet_"+var+"_"+wp+"pass_"+pt+"_data_opt"
		hname_data_untag = "UNWEIGHTED__DATA__FatJet_"+var+"_"+wp+"fail_"+pt+"_data_opt"
		if LTSV and not tag:
			hname_data = "UNWEIGHTED__DATA__FatJet_"+var+"_all_"+pt+"_data_opt"
			hname_data_tag = "UNWEIGHTED__DATA__FatJet_"+var+"hasSV"+"_all_"+pt+"_data_opt"
			hname_data_untag = "UNWEIGHTED__DATA__FatJet_"+var+"noSV"+"_all_"+pt+"_data_opt"
		elif LTSV and tag:
			hname_data = "UNWEIGHTED__DATA__FatJet_"+var+"_"+wp+"pass_"+pt+"_data_opt"
			hname_data_tag = "UNWEIGHTED__DATA__FatJet_"+var+"hasSV_"+wp+"pass_"+pt+"_data_opt"
			hname_data_untag = "UNWEIGHTED__DATA__FatJet_"+var+"noSV_"+wp+"pass_"+pt+"_data_opt"

		cf.SetData(hname_data)   
		cf.SetDataTag(hname_data_tag)
		cf.SetDataUntag(hname_data_untag)	

		# Hist name setup	
		## Nominal names	
		def_name_qcd = "UNWEIGHTED__QCDMu+__FatJet_"
		hname_bfromg = 	def_name_qcd+	var+"_all_"+	pt+"_bfromg_opt"
		hname_b = 		def_name_qcd+	var+"_all_"+	pt+"_b_opt"
		hname_cfromg = 	def_name_qcd+	var+"_all_"+	pt+"_cfromg_opt"
		hname_c = 		def_name_qcd+	var+"_all_"+	pt+"_c_opt"
		hname_l = 		def_name_qcd+	var+"_all_"+	pt+"_l_opt"
		hname_b_cfromg =def_name_qcd+	var+"_all_"+	pt+"_b_cfromg_opt"
		hname_c_l = 	def_name_qcd+	var+"_all_"+	pt+"_c_l_opt"
		hname_bkg = 	def_name_qcd+	var+"_all_"+	pt+"_b_cfromg_c_l_opt"

		# Make a list
		hists_nom = [hname_bfromg, hname_b, hname_cfromg, hname_c, hname_l, hname_b_cfromg, hname_c_l, hname_bkg]
		# Make lists for tag and untag
		hists_tag = [] 
		hists_untag = [] 
		hists_nom_temp = []
		## Tag and untag
		### Regular mode
		if not LTSV : 
			for i, hist_nom in enumerate(hists_nom):
				hist = hist_nom.replace("_all_", "_"+wp+"pass_")
				hists_tag.append(hist)
			for i, hist_nom in enumerate(hists_nom):
				hist = hist_nom.replace("_all_", "_"+wp+"fail_")
				hists_untag.append(hist)
		### LTSV mode all			
		elif LTSV and not tag:
			for i, hist_nom in enumerate(hists_nom):
				hist = hist_nom.replace(var, var+"hasSV")
				hists_tag.append(hist)
			for i, hist_nom in enumerate(hists_nom):
				hist = hist_nom.replace(var, var+"noSV")
				hists_untag.append(hist)

		### LTSV mode tagged		
		elif LTSV and tag:
			for i, hist_nom in enumerate(hists_nom):
				hist = hist_nom.replace(var, var+"hasSV").replace("_all_", "_"+wp+"pass_")
				hists_tag.append(hist)
			for i, hist_nom in enumerate(hists_nom):
				hist = hist_nom.replace(var, var+"noSV").replace("_all_", "_"+wp+"pass_")
				hists_untag.append(hist)
			for i, hist_nom in enumerate(hists_nom):
				hist = hist_nom.replace("_all_", "_"+wp+"pass_")
				hists_nom_temp.append(hist)
			hists_nom = hists_nom_temp
		else: 
			print "Really should not happen"

		# Check if histograms are non-empty
		f = r.TFile(file)

		if addSYS:
			_ups, _downs = ["_"+_f+"up" for _f in systlist], ["_"+_f+"down" for _f in systlist]
			_systVarlist = [""] + list(sum(zip(_ups, _downs+[0]), ()))
		else: _systVarlist = [""]
		
		empty_hist_id_tot = np.array([False] * len(hists_nom))
		for _syst in _systVarlist:
			empty_hist_ix = []			
			for i, (_nom, _tag, _untag) in enumerate(zip(hists_nom, hists_tag, hists_untag)):
				#print f.Get(_nom), f.Get(_tag), f.Get(_untag)
				#print i, _nom, f.Get(_nom).Integral(), f.Get(_tag).Integral(), f.Get(_untag).Integral()	
				try:
					empty_hist_ix.append(any(_int == 0. for _int in [f.Get(_nom+_syst).Integral(), f.Get(_tag+_syst).Integral(), f.Get(_untag+_syst).Integral()]))
				except:
					print [f.Get(_nom+_syst), f.Get(_tag+_syst), f.Get(_untag+_syst)]
					import sys
					sys.exit()
			empty_hist_id_tot = (empty_hist_id_tot ^ empty_hist_ix)
	
		if debug: print "Empty templates :", empty_hist_id_tot
		if any(empty_hist_id_tot[5:7]) : 
			if debug: print "Dynamic merge 2", var, tag
			merge = 2
		elif any(empty_hist_id_tot[0:5]): 
			if debug: print "Dynamic merge 1", var, tag
			merge = 1
		else:
			if debug: print "No merging required"
			merge=0

		if merge == 1: tempNs = ["g #rightarrow b#bar{b}", "b + g #rightarrow c#bar{c}", "c + dusg"]
		elif merge == 2: tempNs = ["g #rightarrow b#bar{b}", "b + g #rightarrow c#bar{c} + c + dusg"]
		else: 	tempNs = ["g #rightarrow b#bar{b}", "b", "g #rightarrow c#bar{c}", "c", "dusg"]

		# AddTemplate(label, name in input file, color)
		## Nominal
		if merge==1:	
			cf.AddTemplate(tempNs[0], hists_nom[0], 65)
			cf.AddTemplate(tempNs[1], hists_nom[5], 628)
			cf.AddTemplate(tempNs[2], hists_nom[6],	597)
			if glue:
				if not LTSV:
					pass 
				else:
					cf.GlueTemplates(tempNs[1:],"other flavours",28);
		elif merge==2:	
			cf.AddTemplate(tempNs[0], hists_nom[0], 65)
			cf.AddTemplate(tempNs[1], hists_nom[7], 28)
		else:
			cf.AddTemplate(tempNs[0], hists_nom[0],	65)
			cf.AddTemplate(tempNs[1], hists_nom[1],	213)
			cf.AddTemplate(tempNs[2], hists_nom[2],	208)
			cf.AddTemplate(tempNs[3], hists_nom[3], 206)
			cf.AddTemplate(tempNs[4], hists_nom[4], 212)
			if glue:
				if not LTSV:
					if not ccSignal: cf.GlueTemplates([tempNs[1], tempNs[2]],"b + g #rightarrow c#bar{c}",628);
					else: cf.GlueTemplates([tempNs[0], tempNs[1]],"b + g #rightarrow b#bar{b}",628);
					cf.GlueTemplates([tempNs[3], tempNs[4]],"c + dusg",597)
				else:
					cf.GlueTemplates(tempNs[1:],"other flavours",28);
		## Taggged templates
		if merge==1:
			cf.AddTemplateTag(tempNs[0], hists_tag[0],		 65)
			cf.AddTemplateTag(tempNs[1], hists_tag[5], 		628)
			cf.AddTemplateTag(tempNs[2], hists_tag[6], 		597)
			if glue:
				if not LTSV:
					pass
				else:
					cf.GlueTemplatesTag(tempNs[1:],"other flavours",28);
		elif merge==2:	
			cf.AddTemplateTag(tempNs[0], hists_tag[0], 65)
			cf.AddTemplateTag(tempNs[1], hists_tag[7], 28)
		else:
			cf.AddTemplateTag(tempNs[0], hists_tag[0],		65)
			cf.AddTemplateTag(tempNs[1], hists_tag[1],		213)
			cf.AddTemplateTag(tempNs[2], hists_tag[2],		208)
			cf.AddTemplateTag(tempNs[3], hists_tag[3],		206)
			cf.AddTemplateTag(tempNs[4], hists_tag[4],		212)
			if glue:
				if not LTSV:
					if not ccSignal: cf.GlueTemplates([tempNs[1], tempNs[2]],"b + g #rightarrow c#bar{c}",628);
					else: cf.GlueTemplates([tempNs[0], tempNs[1]],"b + g #rightarrow b#bar{b}",628);
					cf.GlueTemplatesTag([tempNs[3], tempNs[4]],"c + dusg",597)
				else:
					cf.GlueTemplatesTag(tempNs[1:],"other flavours",28);
		## Untag templates
		if merge==1:
			cf.AddTemplateUntag(tempNs[0], hists_untag[0],	65)
			cf.AddTemplateUntag(tempNs[1], hists_untag[5],	628)
			cf.AddTemplateUntag(tempNs[2], hists_untag[6],	597)
		elif merge==2:	
			cf.AddTemplateUntag(tempNs[0], hists_untag[0], 65)
			cf.AddTemplateUntag(tempNs[1], hists_untag[7], 28)
		else:			
			cf.AddTemplateUntag(tempNs[0], hists_untag[0],	7)
			cf.AddTemplateUntag(tempNs[1], hists_untag[1],	213)
			cf.AddTemplateUntag(tempNs[2], hists_untag[2],	42)
			cf.AddTemplateUntag(tempNs[3], hists_untag[3],	8)
			cf.AddTemplateUntag(tempNs[4], hists_untag[4],	4)

		return tempNs
	
	if not LTSV:
		tempNs = add_templates(cfJP, glue, var="JP", ccSignal=ccSignal)
	else:
		tempNs = add_templates(cfSV, glue, LTSV=False, var="tau1VertexMassCorr", ccSignal=ccSignal)
		tempNsJPnotag = add_templates(cfJP, glue, var="JP", LTSV=True, tag=False, ccSignal=ccSignal)
		tempNsJPtag = add_templates(cfJPtag, glue, var="JP", LTSV=True, tag=True, ccSignal=ccSignal)
		
	def getpars(cf, tempNs, sysVar, statVar, ccSignal=False): 
		# Like CFIT example
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
		effDATA = effMC*par_tag[0]/par[0]
		sf = effDATA/effMC
		if debug: print "{: <24} : EffMC = {} , EffData = {} ".format(cf.GetRun().split("_")[0] ,  effMC, effDATA)
		if debug: print "          Total MC pre-tag / post-tag = {}, signal only = {}".format(nmc, nmc_tag)
		if debug: print "         Signal MC pre-tag / post-tag = {}, signal only = {}".format(nmc1_tag, nmc1_tag)
 
		# r_name = cf.GetRun().split(" ")[0]
		# print r_name
		# try: os.mkdir('out_pics/')
		# except: pass
		# try: os.mkdir('out_pics/{}'.format(r_name))
		# except: pass
		# os.system('mv pics/*.eps out_pics/{}/'.format(r_name))

		if not ccSignal:
			return nmc1, nmc, nmc1_tag, nmc_tag, par[0], par_tag[0], chi2, chi2_tag
		else: 
			return nmc1, nmc, nmc1_tag, nmc_tag, par[2], par_tag[2], chi2, chi2_tag

	def getSF(cfmain, tempNs, cfJPall=None, cfJPtag=None, sysVar="", statVar=-1, LTSV=LTSV):

		if not LTSV :
			nmc1, nmc, nmc1_tag, nmc_tag, par, par_tag, chi2, chi2_tag = getpars(cfmain, tempNs, sysVar=sysVar, statVar=statVar)
			fr = nmc1/nmc;
			fr_tag = nmc1_tag/nmc_tag;	   
			effMC = nmc1_tag/nmc1
			effDATA = nmc_tag/nmc*par_tag/par*fr_tag/fr		
			
			sf = effDATA/effMC
			#SF = sf
			SF = par_tag/par
			pars = [par, par_tag]
			chi2s = [chi2, chi2_tag]

		else:	
			#print "Fitting JPall"
			nmc1, nmc, nmc1_tag, nmc_tag, par, par_tag, chi2, chi2_tag = getpars(cfJPall[0], cfJPall[1], sysVar, statVar)
			#print "Fitting JPtagged" 
			nmc1_tagged, nmc, nmc1_tagged_SV, nmc_tag, parJPtagged, par_tagJPtagged, chi2tag, chi2_tagtag = getpars(cfJPtag[0], cfJPtag[1], sysVar, statVar)			  
			#print "Fitting SVmass"
			nmc1_SV_SVfit, nmc, nmc1_tagged_SV_SVfit, nmc_tag, parSV, par_tagSV, chi2SV, chi2_tagSV = getpars(cfmain, tempNs, sysVar, statVar)
		
			# Resulting SF
			SF = par_tagSV*parJPtagged*par_tag/(par_tagJPtagged*parSV*par)
			pars = [par, par_tag, parJPtagged,par_tagJPtagged, parSV, par_tagSV]
			chi2s = [chi2, chi2_tag, chi2tag, chi2_tagtag, chi2SV, chi2_tagSV]

		return SF, pars, chi2s

	if not LTSV:
		SF, nom_pars, chi2s = getSF(cfJP, tempNs, sysVar="", LTSV=False)
	else:
		SF, nom_pars, chi2s = getSF(cfSV, tempNs, cfJPall=[cfJP,tempNsJPnotag], cfJPtag=[cfJPtag,tempNsJPtag], sysVar="", LTSV=LTSV)
	if debug: print "SF =", SF

	if calcSYS:
		print "Calculating Errors"
		for cf in [cfJP, cfSV, cfJPtag]:
			cf.ProducePlots(False)
			cf.SetMatrixOption("READ")
		sigma_stat = 0
   		sigma_syst_up = 0
		sigma_syst_down = 0
		variances_names = []
		errors = []
		
		#Calculate systematics 1 (Within input file)
		ups, downs = ["_"+f+"up" for f in systlist], ["_"+f+"down" for f in systlist]
		systVarlist = list(sum(zip(ups, downs+[0]), ()))
		if len(systVarlist) != len(ups+downs): 
			print "Incorrect number of systs considered"
		
		for i, sysVar in enumerate(systVarlist):
			if not LTSV:
				SF_sys_i, pars,ch = getSF(cfJP, tempNs, sysVar=sysVar)
			else:
				SF_sys_i, pars, ch = getSF(cfSV, tempNs, cfJPall=[cfJP,tempNsJPnotag], cfJPtag=[cfJPtag, tempNsJPtag], sysVar=sysVar, LTSV=LTSV)
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
		stat_SFs = []
		for i in range(stat_n):
			if not LTSV:
				SF_stat_i, pars, ch = getSF(cfJP, tempNs, sysVar="", statVar=667+i)
			else:
				SF_stat_i, pars, ch = getSF(cfSV, tempNs, cfJPall=[cfJP,tempNsJPnotag], cfJPtag=[cfJPtag, tempNsJPtag], sysVar="", statVar=667+i, LTSV=LTSV)
			stat_SFs.append(SF_stat_i)

		sf_stat_mean = sum(stat_SFs)/float(stat_n)
		for SF_i in stat_SFs:
			sigma_stat += (sf_stat_mean - SF_i)**2/float(stat_n)
		sigma_stat = np.sqrt(sigma_stat);

		# 100 -105 systematics  - Mophing systematics from using CFIT
		for i in range(100,105):
			for cf in [cfJP, cfSV, cfJPtag]:
				cf.ProducePlots(False)
				cf.SetMatrixOption("WRITE")
	
				if i == 100: cf.SetOptimization(OPT_NOCORR);
				if i == 101: cf.SetMorphing(OPTMORPH_CUTOFF,0.25);
				if i == 102: cf.SetMorphing(OPTMORPH_CUTOFF,0.75);
				if i == 103: cf.SetMorphing(OPTMORPH_GEOMETRIC);
				if i == 104:
					#continue
					cf = cfit("")
					cf.SetVerbose(0)
					cf.ProducePlots(False)
					cf.SetOptimization(OPT_MORPH_SGN_SIGMA)
					cf.SetCovarianceMode(COV_MAX)
					cf.SetMorphing(OPTMORPH_CUTOFF,0.5)
				    			     
					cf.SetInputFile(file)
			   
					mname = "matrix_"+pt;
					cf.SetMatrixName(mname+'_'+wp)
					cf.SetMatrixOption("WRITE")

					#cf.SetData(hname_data)   
					#cf.SetDataTag(hname_data_tag)
					#cf.SetDataUntag(hname_data_untag)

					#add_templates(cf, glue, merge)

			if not LTSV:
				SF_sys2_i, par, ch = getSF(cfJP, tempNs, sysVar="", statVar=-1)
			else:
				SF_sys2_i, par, ch = getSF(cfSV, tempNs, cfJPall=[cfJP,tempNsJPnotag], cfJPtag=[cfJPtag, tempNsJPtag], sysVar="", statVar=-1, LTSV=LTSV)
			cf.SetOptimization(OPT_MORPH_SGN_SIGMA)
			cf.SetMorphing(OPTMORPH_CUTOFF,0.5)

			delta = SF - SF_sys2_i
			sigma2 = delta**2
			errors.append(delta)
			variances_names.append("Morph"+str(i))
			print SF_sys2_i, delta
			if SF < SF_sys2_i:
				sigma_syst_up += sigma2
			else:
				sigma_syst_down += sigma2

		# Calculate systematics 2 (from SF values calculated separately like MCJPCalib)
		
		for SF_i in SF_dict.keys():
			if SF_dict[SF_i] == []: continue
			error = SF - SF_dict[SF_i][pt_bin]
			sigma2 = error**2
			print  SF_i, error
			if SF < SF_dict[SF_i][pt_bin]:
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
		variances = [err**2 for err in errors]
		#print errors

	if calcSYS:
		#print np.round(sigma_stat, N_digits_err), np.round(syst_up, N_digits_err), np.round(syst_down, N_digits_err)
		#return SF, np.round(sigma_stat, N_digits_err), np.round(syst_up, N_digits_err), np.round(syst_down, N_digits_err)
		return SF, sigma_stat, syst_up, syst_down, variances_names, errors, variances, nom_pars, chi2s
	else:
		return SF, nom_pars, chi2s

def make_dirs(dirname):
	import os, errno
	try:
		os.makedirs(dirname)
	except OSError, e:
		if e.errno != errno.EEXIST:
			raise

if __name__ == "__main__":
	
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

	WP = "DDBvLM2"
	pt_bins = ['pt350to2000', 'pt430to2000']
	m = 0

	file_name = 'col3/collated_normRun2016_DoubleB.root'
	file_name = 'col3/collatedRun2016_DoubleB.root'
	
	glue=False; addSYS=True; calcSYS=True
	#glue=False; addSYS=False; merge=False; calcSYS=False
	#SF, pars, chi2s = runSF_x(file_name, pt_bins, m, WP, glue=glue, addSYS=addSYS, calcSYS=calcSYS, systname=None, LTSV=True)
	
	# SF, sigma_stat, syst_up, syst_down, variances_names, errors, variances, nom_pars, chi2 = runSF_x(file_name, 
	# 	pt_bins, m, WP, glue=glue, addSYS=addSYS, calcSYS=calcSYS, systname=None, LTSV=True)
	print "DoubleB 16"
	file_name = 'col_fin450/collated_normRun2016_DoubleB.root'
	SF, pars, chi2s = runSF_x(file_name, pt_bins, m, "DoubleBM1", glue=glue, addSYS=False, calcSYS=False, systname=None, LTSV=True, ccSignal=False, debug=True)
	print "DoubleB 17"
	file_name = 'col_fin450/collated_normRun2017_DoubleB.root'
	SF, pars, chi2s = runSF_x(file_name, pt_bins, m, "DoubleBM1", glue=glue, addSYS=False, calcSYS=False, systname=None, LTSV=True, ccSignal=False, debug=True)
	
	# file_name = 'colnom/collated_normRun2016_DDBvL.root'
	# print "DDBvL"
	# SF, pars, chi2s = runSF_x(file_name, pt_bins, m, "DDBvLM2", glue=glue, addSYS=False, calcSYS=False, systname=None, LTSV=True, ccSignal=False, debug=True)
	# print "DeepAK8"
	# file_name = 'colnom/collated_normRun2016_DeepAK8ZHbb.root'
	# SF, pars, chi2s = runSF_x(file_name, pt_bins, m, "DeepAK8ZHbbM2", glue=glue, addSYS=False, calcSYS=False, systname=None, LTSV=True	, ccSignal=False, debug=True)

	# Threading example ( Doesn't run faster though, likely due to the wrapped cfit class)
	# import threading
	# threading.Thread( target=runSF_x, 
	# 	args=(file_name, pt_bins, m, WP,), 
	# 	kwargs={'glue':glue, 'addSYS':addSYS, 'calcSYS':calcSYS, 'systname':None, 'LTSV':True},
	# 	name='thread_function').start()	
	

