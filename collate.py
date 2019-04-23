from __future__ import print_function
import os, sys
import ROOT as R
import numpy as np

import argparse
parser = argparse.ArgumentParser(description='Process some integers.')
parser.add_argument('-n', '--name', default='Run2017_DDBvL', help="Name of the tagger, to load WPs etc...")
parser.add_argument('--debug', action="store_true", default=False, help="")
args = parser.parse_args()
print("Run name (--name): {}".format(args.name))

din = 'plots_final5'
#pt_bins = ['pt350to430', 'pt430to2000', 'pt350to2000']
pt_bins = ['pt350to450', 'pt450to2000', 'pt350to2000']
names = [f.strip(".root") for f in os.listdir(din) if args.name in f]
tagger = args.name.split("_")[-1]

root_file_dict = {}
for name in names:
	root_file_dict[name] = R.TFile(os.path.join(din, name+'.root'))
systnames = [n.replace(args.name, "").strip("_") for n in names] # Autmated available systnames
print("Systematics found: {}".format(systnames))
print("Additional systematics: {}".format("JES, GCC, GBB"))

# Read tagger values
import json
WPs, WPs_value = [], []
with open('DDX.json') as json_file:
	data = json.load(json_file)[tagger]
	for key, value in sorted(data.iteritems(), key=lambda (k,v): (v,k)):
		WPs.append(str(key))
		WPs_value.append(value)

def collate_systematics(run_name=args.name,
						rfiles = root_file_dict,
						tagger = tagger,
						WPs = WPs, 
						shapein=False,
						debug=args.debug):		
	fout = R.TFile("col_fin450/collated{}.root".format(run_name), 'RECREATE')
	fout.cd()
	out_names = ["UNWEIGHTED__DATA__FatJet_", "UNWEIGHTED__QCDMu+__FatJet_"] 
	in_names = ["DATA__FatJet_", "QCDMu+__FatJet_"]
	for iQD, QD in enumerate(["data", "qcd"]):
		for LTSVvar in ["JP", "JPhasSV", "JPnoSV", "tau1VertexMassCorr"]:
			for ptbin in pt_bins:
				for WP in WPs:
					for passfail in ["all", "{}_pass", "{}_fail"]:
						if passfail == "all" and WPs.index(WP) > 0: continue
						passfail = passfail.format(tagger+"_"+WP)
						if debug: print ("    "+in_names[0]+LTSVvar+"_"+passfail+"_"+ptbin+"_"+"<flavor>")
						# Clone default histos
						for flavor in ["data", "b", "bfromg", "c", "cfromg", "l"]:
							if QD == "data" and not flavor == "data": continue
							if QD == "qcd" and flavor == "data": continue
							for systname in systnames:
								if len(systname) < 1: # No syst case					
									hin = rfiles[args.name].Get(in_names[iQD]+LTSVvar+"_"+passfail+"_"+ptbin+"_"+flavor)
									hout_name = out_names[iQD]+LTSVvar+"_"+passfail.replace("_", "")+"_"+ptbin+"_"+flavor+"_opt"									
									try: hout = hin.Clone(hout_name)
									except: print(in_names[iQD]+LTSVvar+"_"+passfail+"_"+ptbin+"_"+flavor)
									hout.Write()
								elif not systname.endswith("UP") and not systname.endswith("DOWN"): # one sided syst case
									# Copy syst as up
									hin = rfiles[run_name+"_"+systname].Get(in_names[iQD]+LTSVvar+"_"+passfail+"_"+ptbin+"_"+flavor)
									hout_name = out_names[iQD]+LTSVvar+"_"+passfail.replace("_", "")+"_"+ptbin+"_"+flavor+"_opt_"+systname+"up"
									hout = hin.Clone(hout_name)
									hout.Write()
									# Set down to nominal
									hin = rfiles[run_name].Get(in_names[iQD]+LTSVvar+"_"+passfail+"_"+ptbin+"_"+flavor)
									hout_name = out_names[iQD]+LTSVvar+"_"+passfail.replace("_", "")+"_"+ptbin+"_"+flavor+"_opt_"+systname+"down"
									hout = hin.Clone(hout_name)
									hout.Write()
								else: # Do two sided separately
									if systname.endswith("UP"): systupdown = systname.replace("UP", "up")
									elif systname.endswith("DOWN"): systupdown = systname.replace("DOWN", "down")
									hin = rfiles[run_name+"_"+systname].Get(in_names[iQD]+LTSVvar+"_"+passfail+"_"+ptbin+"_"+flavor)
									hout_name = out_names[iQD]+LTSVvar+"_"+passfail.replace("_", "")+"_"+ptbin+"_"+flavor+"_opt_"+systupdown
									hout = hin.Clone(hout_name)
									hout.Write()

							for jes in ["JESup", "JESdown"]:	
								hin = rfiles[args.name].Get(in_names[iQD]+LTSVvar+"_"+passfail+"_"+ptbin+"_"+jes+"_"+flavor)
								hout_name = out_names[iQD]+LTSVvar+"_"+passfail.replace("_", "")+"_"+ptbin+"_"+flavor+"_opt_"+jes
								hout = hin.Clone(hout_name)
								hout.Write()

							for systname, shape_name in zip(["GBB", "GCC"], ['bfromg','cfromg']):	
								for mult, updown in zip([1.5, 0.5], ["up", "down"]):
									hin = rfiles[args.name].Get(in_names[iQD]+LTSVvar+"_"+passfail+"_"+ptbin+"_"+flavor)
									hout_name = out_names[iQD]+LTSVvar+"_"+passfail.replace("_", "")+"_"+ptbin+"_"+flavor+"_opt_"+systname+updown
									hout = hin.Clone(hout_name)
									if flavor == shape_name: hout.Scale(mult) 
									hout.Write()

						# Merge some QCD templates
						if QD == "data": continue
						for merge, temps in zip (["b_cfromg", "b_bfromg", "c_l", "b_cfromg_c_l", "b_bfromg_c_l"], 
												[["b","cfromg"], ["b","bfromg"], ["c","l"], ["b", "c", "cfromg", "l"], ["b", "c", "bfromg", "l"]]):
							for systname in systnames:
								if len(systname) < 1: # No syst case									
									hin = rfiles[run_name].Get(in_names[iQD]+LTSVvar+"_"+passfail+"_"+ptbin+"_"+temps[0])
									hout_name = out_names[iQD]+LTSVvar+"_"+passfail.replace("_","")+"_"+ptbin+"_"+merge+"_opt"
									hout = hin.Clone(hout_name)
									for temp in temps[1:]:
										hout.Add(rfiles[run_name].Get(in_names[iQD]+LTSVvar+"_"+passfail+"_"+ptbin+"_"+temp))
									hout.Write()
								elif not systname.endswith("UP") and not systname.endswith("DOWN"): # one sided syst case
									# Copy syst as up
									hin = rfiles[run_name+"_"+systname].Get(in_names[iQD]+LTSVvar+"_"+passfail+"_"+ptbin+"_"+temps[0])
									hout_name = out_names[iQD]+LTSVvar+"_"+passfail.replace("_","")+"_"+ptbin+"_"+merge+"_opt_"+systname+"up"
									hout = hin.Clone(hout_name)
									for temp in temps[1:]:
										hout.Add(rfiles[run_name+"_"+systname].Get(in_names[iQD]+LTSVvar+"_"+passfail+"_"+ptbin+"_"+temp))
									hout.Write()
									# Set down to nominal
									hin = rfiles[run_name].Get(in_names[iQD]+LTSVvar+"_"+passfail+"_"+ptbin+"_"+temps[0])
									hout_name = out_names[iQD]+LTSVvar+"_"+passfail.replace("_","")+"_"+ptbin+"_"+merge+"_opt_"+systname+"down"
									hout = hin.Clone(hout_name)
									for temp in temps[1:]:
										hout.Add(rfiles[run_name].Get(in_names[iQD]+LTSVvar+"_"+passfail+"_"+ptbin+"_"+temp))
									hout.Write()
								else: # Do two sided syst one by one 									
									hin = rfiles[run_name+"_"+systname].Get(in_names[iQD]+LTSVvar+"_"+passfail+"_"+ptbin+"_"+temps[0])
									if systname.endswith("UP"): systupdown = systname.replace("UP", "up")
									elif systname.endswith("DOWN"): systupdown = systname.replace("DOWN", "down")
									hout_name = out_names[iQD]+LTSVvar+"_"+passfail.replace("_", "")+"_"+ptbin+"_"+merge+"_opt_"+systupdown
									hout = hin.Clone(hout_name)
									for temp in temps[1:]:
										hout.Add(rfiles[run_name+"_"+systname].Get(in_names[iQD]+LTSVvar+"_"+passfail+"_"+ptbin+"_"+temp))
									hout.Write()
							for jes in ["JESup", "JESdown"]:							
								hin = rfiles[run_name].Get(in_names[iQD]+LTSVvar+"_"+passfail+"_"+ptbin+"_"+jes+"_"+flavor)
								hout_name = out_names[iQD]+LTSVvar+"_"+passfail.replace("_", "")+"_"+ptbin+"_"+merge+"_opt_"+jes
								hout = hin.Clone(hout_name)
								for temp in temps[1:]:									
									hout.Add(rfiles[run_name].Get(in_names[iQD]+LTSVvar+"_"+passfail+"_"+ptbin+"_"+jes+"_"+temp))
								hout.Write()

							for systname, shape_name in zip(["GBB", "GCC"], ['bfromg','cfromg']):
								for mult, updown in zip([1.5, 0.5], ["up", "down"]):
									hin = rfiles[args.name].Get(in_names[iQD]+LTSVvar+"_"+passfail+"_"+ptbin+"_"+temps[0])
									hout_name = out_names[iQD]+LTSVvar+"_"+passfail.replace("_", "")+"_"+ptbin+"_"+merge+"_opt_"+systname+updown
									hout = hin.Clone(hout_name)
									if temps[0] == shape_name: hout.Scale(mult)
									for temp in temps[1:]:
										htemp = rfiles[args.name].Get(in_names[iQD]+LTSVvar+"_"+passfail+"_"+ptbin+"_"+temp)	
										if temp == shape_name: htemp.Scale(mult)				
										hout.Add(htemp)
									hout.Write()

	fout.Close()
	print("Collated syst to {}".format(fout))

def produce_scalevars(	scaletemplate = None,
						updown = "up",
						tagger = tagger,
						run_name=args.name,
						debug=args.debug,
						keepsys=True
						):
	
	fin = R.TFile("col_fin450/collated{}.root".format(run_name))
	if scaletemplate != None:
		fout = R.TFile("col_fin450/collated_norm{}.root".format(run_name+"_"+scaletemplate+"_"+updown), 'RECREATE')
	else:
		fout = R.TFile("col_fin450/collated_norm{}.root".format(run_name), 'RECREATE')
	fout.cd()
	
	for LTSVvar in ["JP", "JPhasSV", "JPnoSV", "tau1VertexMassCorr"]:
		for ptbin in pt_bins:
			for WP in WPs:
				for passfail in ["all", "{}_pass", "{}_fail"]:
					if passfail == "all" and WPs.index(WP) > 0: continue
					passfail = passfail.format(tagger+"_"+WP)
					if keepsys: 
						systlist = [""]+[s.replace("UP","").replace("DOWN", "")+"up" for s in systnames if len(s) > 0]+[s.replace("UP","").replace("DOWN","")+"down" for s in systnames if len(s) > 0] # Must include up/down
						systlist = systlist + ["JESup", "JESdown", "GBBdown", "GCCdown", "GBBup", "GCCup"]
					else: systlist = [""]
					for systname in systlist :
						if len(systname) > 0: systname = "_"+systname
						if debug: print(systname)
						hout_name = "{}"+LTSVvar+"_"+passfail.replace("_", "")+"_"+ptbin+"_{}_opt{}"								
						if debug: print(hout_name.format("UNWEIGHTED__DATA__FatJet_", "data", systname))
						int_data = fin.Get(hout_name.format("UNWEIGHTED__DATA__FatJet_", "data", systname)).Integral()
						int_flavs, int_old = [], []
						for flavor in ["b", "bfromg", "c", "cfromg", "l"]:
							int_h = fin.Get(hout_name.format("UNWEIGHTED__QCDMu+__FatJet_", flavor, systname)).Integral()
							if flavor == scaletemplate and updown == "up": 
								int_flavs.append(int_h*1.5)
							elif flavor == scaletemplate and updown == "down": 
								int_flavs.append(int_h*0.5)
							else:
								int_flavs.append(int_h)
							int_old.append(int_h) # Debug
						int_qcd = np.sum(int_flavs)
						norm = 1#int_data/(int_qcd+0.0000000000000000000000000000000001
						int_flavs_corr = []
						hout = fin.Get(hout_name.format("UNWEIGHTED__DATA__FatJet_", "data", systname))
						hout.Write()
						for flavor in ["b", "bfromg", "c", "cfromg", "l", "b_cfromg", "c_l", "b_cfromg_c_l"]:
							hout = fin.Get(hout_name.format("UNWEIGHTED__QCDMu+__FatJet_", flavor, systname))
							if scaletemplate!= None and scaletemplate in flavor.split("_"):
								if updown == "up": 
									hout.Scale(1.5)							
								elif updown == "down": 
									hout.Scale(0.5)
							else:
								pass
							hout.Scale(norm)
							int_flavs_corr.append(hout.Integral())
							hout.Write()
						#for merge, temps in zip (["b_cfromg", "c_l", "b_cfromg_c_l"], [["b","cfromg"], ["c","l"], ["b", "c", "cfromg", "l"]]):

	fout.Close()
	print("Normed data/MC templates, while scaling {} 50% {}.".format(scaletemplate, updown))
	print("Output written to: {}".format(fout))

if __name__ == "__main__":
	collate_systematics()
	produce_scalevars(keepsys=True) # Produce norms
	sys.exit()
	for ud in ["up", "down"]:
		for temp in ['b', 'bfromg', 'c', 'cfromg', 'l']:
			produce_scalevars(updown=ud, scaletemplate=temp)
