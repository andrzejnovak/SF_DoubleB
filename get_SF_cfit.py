#from __future__ import unicode_literals
from cfit import *
import numpy as np

def runSF_x(file, pt, wp, merge=False, glue=True, inclSYS=True, SF_dict={}, SV=True):
	N_digits_SF = 3
	N_digits_err = 4
	if pt == "pt250to350": pt_bin = 0
	if pt == "pt350to430": pt_bin = 1
	if pt == "pt430to2000": pt_bin = 2

	cfJP = cfit("JP discriminator")
	cfSV = cfit("SV discriminator")
	cfJPtag = cfit("JPtagged discriminator")

	cfJP.SetMatrixName("matrix_"+pt)
	cfJP.SetMatrixOption("WRITE")
	cfSV.SetMatrixName("SVmatrix_"+pt)
	cfSV.SetMatrixOption("WRITE")
	cfJPtag.SetMatrixName("JPtagmatrix_"+pt)
	cfJPtag.SetMatrixOption("WRITE")

	# Settings
	for cf in [cfJP, cfSV, cfJPtag]:
		cf.SetVerbose(0)
		cf.ProducePlots(True)
		cf.SetLegendHeader(pt+" "+wp)
		cf.SetOptimization(OPT_MORPH_SGN_SIGMA)
		cf.SetCovarianceMode(COV_MAX)
		cf.SetMorphing(OPTMORPH_CUTOFF,0.5)
		cf.SetInputFile(file)
	

	if inclSYS:
		#systlist = ["JES", "NTRACKS", "BFRAG", "CFRAG", "CD", "K0L", "PU"]
		systlist = ["JES",  "BFRAG", "CFRAG", "CD", "K0L", "PU"]
		for sysName in systlist:
			# Flipped down up in Alice's code
			#pass
			cfJP.AddSys(sysName, "_"+sysName+"down" ,"_"+sysName+"up")
			cfSV.AddSys(sysName, "_"+sysName+"down" ,"_"+sysName+"up")
			cfJPtag.AddSys(sysName, "_"+sysName+"down" ,"_"+sysName+"up")
			#cf.AddSys(sysName, "_"+sysName+"up" ,"_"+sysName+"down") #-- this one is correct

	# Mergning for low count templates:
	if merge: tempNs = ["g #rightarrow b#bar{b}", "b + g #rightarrow c#bar{c}", "c + dusg"]
	else: 	tempNs = ["g #rightarrow b#bar{b}", "b", "g #rightarrow c#bar{c}", "c", "dusg"]
	def add_templates(cf, glue, merge, var = "JP", SV=SV, tag=False):
		# Set DATA names
		hname_data = "UNWEIGHTED__DATA__FatJet_"+var+"_all_"+pt+"_data_opt"
		hname_data_tag = "UNWEIGHTED__DATA__FatJet_"+var+"_"+wp+"pass_"+pt+"_data_opt"
		hname_data_untag = "UNWEIGHTED__DATA__FatJet_"+var+"_"+wp+"fail_"+pt+"_data_opt"
		if SV and not tag:
			hname_data = "UNWEIGHTED__DATA__FatJet_"+var+"_all_"+pt+"_data_opt"
			hname_data_tag = "UNWEIGHTED__DATA__FatJet_"+var+"hasSV"+"_all_"+pt+"_data_opt"
			hname_data_untag = "UNWEIGHTED__DATA__FatJet_"+var+"noSV"+"_all_"+pt+"_data_opt"
		elif SV and tag:
			hname_data = "UNWEIGHTED__DATA__FatJet_"+var+"_"+wp+"pass_"+pt+"_data_opt"
			hname_data_tag = "UNWEIGHTED__DATA__FatJet_"+var+"hasSV_"+wp+"pass_"+pt+"_data_opt"
			hname_data_untag = "UNWEIGHTED__DATA__FatJet_"+var+"noSV_"+wp+"pass_"+pt+"_data_opt"

		cf.SetData(hname_data)   
		cf.SetDataTag(hname_data_tag)
		cf.SetDataUntag(hname_data_untag)	

		# Hist name setup	
		## Nominal names	
		def_name_qcd = "UNWEIGHTED__QCDMuEnr__FatJet_"
		hname_bfromg = 	def_name_qcd+	var+"_all_"+	pt+"_bfromg_opt"
		hname_b = 		def_name_qcd+	var+"_all_"+	pt+"_b_opt"
		hname_cfromg = 	def_name_qcd+	var+"_all_"+	pt+"_cfromg_opt"
		hname_c = 		def_name_qcd+	var+"_all_"+	pt+"_c_opt"
		hname_l = 		def_name_qcd+	var+"_all_"+	pt+"_l_opt"
		hname_b_cfromg =def_name_qcd+	var+"_all_"+	pt+"_b_cfromg_opt"
		hname_c_l = 	def_name_qcd+	var+"_all_"+	pt+"_c_l_opt"

		# Make a list
		hists_nom = [hname_bfromg, hname_b, hname_cfromg, hname_c, hname_l, hname_b_cfromg, hname_c_l]
		# Make lists for tag and untag
		#0 					 1 				2 					3				4 				5 						6
		#hname_bfromg_tag,	 hname_b_tag, 	hname_cfromg_tag, 	hname_c_tag, 	hname_l_tag, 	hname_b_cfromg_tag, 	hname_c_l_tag 
		#hname_bfromg_untag, hname_b_untag, hname_cfromg_untag, hname_c_untag, 	hname_l_untag, 	hname_b_cfromg_untag, 	hname_c_l_untag 
		hists_tag = [] 
		hists_untag = [] 
		hists_nom_temp = []
		## Tag and untag
		### Regular mode
		if not SV : 
			for i, hist_nom in enumerate(hists_nom):
				hist = hist_nom.replace("_all_", "_"+wp+"pass_")
				hists_tag.append(hist)
			for i, hist_nom in enumerate(hists_nom):
				hist = hist_nom.replace("_all_", "_"+wp+"fail_")
				hists_untag.append(hist)
		### SV mode all			
		elif SV and not tag:
			for i, hist_nom in enumerate(hists_nom):
				hist = hist_nom.replace(var, var+"hasSV")
				hists_tag.append(hist)
			for i, hist_nom in enumerate(hists_nom):
				hist = hist_nom.replace(var, var+"noSV")
				hists_untag.append(hist)

		### SV mode tagged		
		elif SV and tag:
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

		#import sys
		#sys.exit()
		# AddTemplate(label, name in input file, color)
		if merge:	
			cf.AddTemplate(tempNs[0], hists_nom[0]	,65)
			cf.AddTemplate(tempNs[1], hists_nom[5],628)
			cf.AddTemplate(tempNs[2], hists_nom[6],	597)
			if glue and SV:
				cf.GlueTemplates(tempNs[1:],"other flavours",28);
		else:
			cf.AddTemplate(tempNs[0], hists_nom[0],	65)
			cf.AddTemplate(tempNs[1], hists_nom[1],		213)
			cf.AddTemplate(tempNs[2], hists_nom[2],	208)
			cf.AddTemplate(tempNs[3], hists_nom[3], 		206)
			cf.AddTemplate(tempNs[4], hists_nom[4], 		212)
			if glue:
				if not SV:
					cf.GlueTemplates([tempNs[1], tempNs[2]],"b + g #rightarrow c#bar{c}",628);
					cf.GlueTemplates([tempNs[3], tempNs[4]],"c + dusg",597)
				else:
					cf.GlueTemplates(tempNs[1:],"other flavours",28);

		if merge:
			cf.AddTemplateTag(tempNs[0], hists_tag[0],		 65)
			cf.AddTemplateTag(tempNs[1], hists_tag[5], 		628)
			cf.AddTemplateTag(tempNs[2], hists_tag[6], 		597)
			if glue and SV:
				cf.GlueTemplatesTag(tempNs[1:],"other flavours",28);
		else:
			cf.AddTemplateTag(tempNs[0], hists_tag[0],		65)
			cf.AddTemplateTag(tempNs[1], hists_tag[1],		213)
			cf.AddTemplateTag(tempNs[2], hists_tag[2],		208)
			cf.AddTemplateTag(tempNs[3], hists_tag[3],		206)
			cf.AddTemplateTag(tempNs[4], hists_tag[4],		212)
			if glue:
				if not SV:
					cf.GlueTemplatesTag([tempNs[1], tempNs[2]],"b + g #rightarrow c#bar{c}",628);
					cf.GlueTemplatesTag([tempNs[3], tempNs[4]],"c + dusg",597)
				else:
					cf.GlueTemplatesTag(tempNs[1:],"other flavours",28);

		if merge:
			cf.AddTemplateUntag(tempNs[0], hists_untag[0],	65)
			cf.AddTemplateUntag(tempNs[1], hists_untag[5],	628)
			cf.AddTemplateUntag(tempNs[2], hists_untag[6],	597)
		else:			
			cf.AddTemplateUntag(tempNs[0], hists_untag[0],	7)
			cf.AddTemplateUntag(tempNs[1], hists_untag[1],	213)
			cf.AddTemplateUntag(tempNs[2], hists_untag[2],	42)
			cf.AddTemplateUntag(tempNs[3], hists_untag[3],	8)
			cf.AddTemplateUntag(tempNs[4], hists_untag[4],	4)
	
	if not SV:
		#add_templates(cf, glue, merge, var="JP")
		pass
	else:
		add_templates(cfSV, glue, merge, SV=False, var="tau1VertexMassCorr")
		add_templates(cfJP, glue, merge, var="JP", SV=True, tag=False)
		add_templates(cfJPtag, glue, merge, var="JP", SV=True, tag=True)
		

	#def getpars(cf, tempNs, sysVar="", statVar=-1):
	def getpars(cf, tempNs, sysVar, statVar): # TO make it work with multiprocessing
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

		#print chi2
		#print chi2_tag

		fr = nmc1/nmc;
		fr_tag = nmc1_tag/nmc_tag;	   
		effMC = nmc1_tag/nmc1
		#effDATA = nmc_tag/nmc*par_tag[0]/par[0]*fr_tag/fr		
		#effDATA =  nmc1_tag/nmc1 * par_tag[0]/par[0]
		effDATA = effMC*par_tag[0]/par[0]

		sf = effDATA/effMC
		return nmc1, nmc, nmc1_tag, nmc_tag, par[0], par_tag[0]

	def getSF(cfmain, tempNs, cfJPall=None, cfJPtag=None, sysVar="", statVar=-1, SV=SV):

		if not SV :
			if cfJPtag!=None and cfJPall!=None: 
				raise "Something is wrong"
			nmc1, nmc, nmc1_tag, nmc_tag, par, par_tag = getpars(cfmain, tempNs, sysVar=sysVar, statVar=statVar)
			fr = nmc1/nmc;
			fr_tag = nmc1_tag/nmc_tag;	   
			effMC = nmc1_tag/nmc1
			#effDATA = nmc_tag/nmc*par_tag[0]/par[0]*fr_tag/fr		
			#effDATA =  nmc1_tag/nmc1 * par_tag[0]/par[0]
			effDATA = effMC*par_tag/par

			sf = effDATA/effMC
			SF = sf
			#SF = par_tag/par

		else:	
			import os

			#try:
			#	os.system('rm -r JPfit')
			#	os.system('rm -r JPtagfit')
			#	os.system('rm -r SV')
			#except:
			#	print "dirs have not yet been created, all good"
			#CFJP		
			#print "Fitting JPall"
			nmc1, nmc, nmc1_tag, nmc_tag, par, par_tag = getpars(cfJPall, tempNs, sysVar, statVar)
			#ndata1 = nmc1*par			 
			#nmc1_SV = nmc1_tag
			#ndata1_SV = nmc1_SV*par_tag
			#effMC_SV = nmc1_SV/nmc1
			#effDATA_SV = ndata1_SV/ndata1;
			#os.mkdir('JPfit')
			#os.system('mv pics/*png JPfit')
			
			#print "Fitting JPtagged" 
			nmc1_tagged, nmc, nmc1_tagged_SV, nmc_tag, parJPtagged, par_tagJPtagged = getpars(cfJPtag, tempNs, sysVar, statVar)			  
			#ndata1_tagged = nmc1_tagged*parJPtagged		   
			#ndata1_tagged_SV = nmc1_tagged_SV*par_tagJPtagged
			#effMC_SV_tagged = nmc1_tagged_SV/nmc1_tagged;
			#effDATA_SV_tagged = ndata1_tagged_SV/ndata1_tagged;
			#os.mkdir('JPtagfit')
			#os.system('mv pics/*png JPtagfit')
					
			#print "Fitting SVmass"
			nmc1_SV_SVfit, nmc, nmc1_tagged_SV_SVfit, nmc_tag, parSV, par_tagSV = getpars(cfmain, tempNs, sysVar, statVar)
			#ndata1_SV_SVfit    = nmc1_SV_SVfit*parSV
			#ndata1_tagged_SV_SVfit = nmc1_tagged_SV_SVfit*par_tagSV
			#nmc1_tag_estim = nmc1_tagged_SV_SVfit/effMC_SV_tagged;
			#nmc1_estim     = nmc1_SV_SVfit/effMC_SV;  
			#effMC = nmc1_tag_estim/nmc1_estim;
			#os.mkdir('SV')
			#os.system('mv pics/*png SV')
			
			#data1_tag_estim = ndata1_tagged_SV_SVfit/effDATA_SV_tagged;
			#ndata1_estim     = ndata1_SV_SVfit/effDATA_SV; 
			#effDATA = ndata1_tag_estim/ndata1_estim;
			#SF = effDATA/effMC; 

			#import multiprocessing as mp
			#print "MULTI"
			#from pathos.multiprocessing import ProcessingPool as Pool
			#p = Pool(4)
			#all_pars = (p.map(getpars, [[cfJPall, tempNs, sysVar, statVar], [cfJPtag, tempNs, sysVar, statVar], [cfmain, tempNs, sysVar, statVar]]))
			#print all_pars

			SF = par_tagSV*parJPtagged*par_tag/(par_tagJPtagged*parSV*par)
			#print "parSV_tag", par_tagSV
			#print "parJPtag", parJPtagged
			#print "parJP_tag", par_tag
			#print "parJPtag_tag", par_tagJPtagged
			#print "parSV", parSV
			#print "parJP", par

			pars = [par, par_tag, parJPtagged,par_tagJPtagged, parSV, par_tagSV]

			#os.system('mv JPfit pics/')
			#os.system('mv JPtagfit pics/')
			#os.system('mv SV pics/')
		
		return SF, pars

	#if inclSYS==False:
	if not SV:
		SF = getSF(cf, tempNs, sysVar="", SV=False)
	else:
		SF, pars = getSF(cfSV, tempNs, cfJPall=cfJP, cfJPtag=cfJPtag, sysVar="", SV=SV)
	print "SF =", SF

	if inclSYS==True:
		print "Calculating Errors"
		for cf in [cfJP, cfSV, cfJPtag]:
			cf.ProducePlots(False)
			cf.SetMatrixOption("READ")
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
		
		for i, sysVar in enumerate(systVarlist):
			if not SV:
				SF_sys_i = getSF(cf, tempNs, sysVar=sysVar)
			else:
				SF_sys_i, pars = getSF(cfSV, tempNs, cfJPall=cfJP, cfJPtag=cfJPtag, sysVar=sysVar, SV=SV)
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
		for i in range(20):
			if not SV:
				SF_stat_i, pars = getSF(cf, tempNs, sysVar="", statVar=667+i)
			else:
				SF_stat_i, pars = getSF(cfSV, tempNs, cfJPall=cf, cfJPtag=cfJPtag, sysVar="", statVar=667+i, SV=SV)
			stat_SFs.append(SF_stat_i)

		sf_stat_mean = sum(stat_SFs)/20.
		for SF_i in stat_SFs:
			sigma_stat += (sf_stat_mean - SF_i)**2/20.
		sigma_stat = np.sqrt(sigma_stat);

		# 100 -105 systematics missing
		for i in range(100,105):
			for cf in [cfJP, cfSV, cfJPtag]:
				cf.ProducePlots(False)
				cf.SetMatrixOption("WRITE")
	
			if i == 100: cf.SetOptimization(OPT_NOCORR);
			if i == 101: cf.SetMorphing(OPTMORPH_CUTOFF,0.25);
			if i == 102: cf.SetMorphing(OPTMORPH_CUTOFF,0.75);
			if i == 103: cf.SetMorphing(OPTMORPH_GEOMETRIC);
			if i == 104:
				continue
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

				#cf.SetData(hname_data)   
				#cf.SetDataTag(hname_data_tag)
				#cf.SetDataUntag(hname_data_untag)

				#add_templates(cf, glue, merge)

			if not SV:
				SF_sys2_i = getSF(cf, tempNs, sysVar="", statVar=-1)
			else:
				SF_sys2_i, par = getSF(cfSV, tempNs, cfJPall=cfJP, cfJPtag=cfJPtag, sysVar="", statVar=-1, SV=SV)
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
			error = SF - SF_dict[SF_i][pt_bin]
			sigma2 = error**2
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
		#print variances_names
		#errors = [np.sqrt(var) for var in variances]
		variances = [err**2 for err in errors]
		#print errors

	if inclSYS:
		print np.round(sigma_stat, N_digits_err), np.round(syst_up, N_digits_err), np.round(syst_down, N_digits_err)
		#return SF, np.round(sigma_stat, N_digits_err), np.round(syst_up, N_digits_err), np.round(syst_down, N_digits_err)
		return SF, sigma_stat, syst_up, syst_down, variances_names, errors, variances
	else:
		return SF, pars

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

	WP = "DoubleBH"
	pt_bin = "pt430to2000"


	#SF, sigma_stat, syst_up, syst_down, variances_names, errors, variances = runSF_x(fname, pt, tag, glue=True, inclSYS=True, SF_dict=SF_dict) 
	#SF = runSF_x(fname, pt, wp, glue=False, inclSYS=False, SF_dict=SF_dict) 
	r = 'singlefilenorm/'
	#r = 'normfile/'
	file_name = r+'Run2017BCDEF_ReReco_QCDMuonEnriched_AK8Jet300orAK4Jet300_Pt350_Final_DoubleMuonTaggedFatJets_histograms_btagval_allVars_ptReweighted_SysMerged_SFtemplates_DoubleBH.root'
	glue=True;  inclSYS=False; merge=False
	SF, pars = runSF_x(file_name, pt_bin, WP, merge=merge, glue=glue, inclSYS=inclSYS)		
	