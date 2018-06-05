import ROOT

r = 'May17/'
#r = 'May17JP/'

suff = "_ptReweighted_SysMerged_SFtemplates_"  
#suff = "_ptReweighted_dataUseMCJPcalib_SysMerged_SFtemplates_"
wp = "DoubleBL"
root = ".root"

import glob
import os
for trig in [trig1, trig2]:
	for wp in ["DoubleBL", "DoubleBM1", "DoubleBM2", "DoubleBH"]:#, "JPnoSV", "JPhasSV", "SVmass"]:
		for updown in ["",'_b_0p5', '_b_1p5', '_c_0p5', '_c_1p5', '_cfromg_0p5', '_cfromg_1p5', '_l_0p5', '_l_1p5']:
			var = "*"
			tohadd = []
			#print r+trig+"*"+var+suff+wp+updown+"_ADDBINNING"+root
			#print 'Run2017BCDEF_ReReco_QCDMuonEnriched_AK8Jet300orAK4Jet300_Pt350_Final_DoubleMuonTaggedFatJets_histograms_btagval_v4_SVmass_ptReweighted_dataUseMCJPcalib_SysMerged_SFtemplates_DoubleBM1_ADDBINNING.root'

			tohadd += glob.glob(r+trig+"*"+var+suff+wp+updown+"_ADDBINNING"+root)
			#tohadd += glob.glob(r+trig+"*"+var+suff+wp+updown+root)
			print tohadd
		#print r+trig+"*"+var+suff+wp+root

			try:
				os.system("hadd -f May17single/"+trig+"allVars"+suff+wp+updown+root+" "+tohadd[0]+" "+tohadd[1]+" "+tohadd[2]+" "+tohadd[3])
			except:
				pass
