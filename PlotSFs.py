import ROOT
from ROOT import *
from ROOT import gROOT
from array import array
import numpy as np

def SFComp(WP, SFs, Stat, Sysup, Sysdown):
    # Settings
    # Bins
    bins = [250, 350, 430, 840]       
    # Inputs
    #SFs, Stat, Sysup, Sysdown = [1, 0.99, 1.01], [0.014, 0.032, 0.023], [0.015, 0.03, 0.018], [0.013, 0.017, 0.01]  #BL

    # Bin center/size calc
    bin_tuples = [bins[max(i, 0):i + 2] for i in range(0, len(bins)-1)] # Rewrite as tuples
    bin_center, bin_size = [], []  
    for i in bin_tuples:
        cent = (i[0]+i[1])/2
        size = cent - i[0]
        bin_center.append(cent); bin_size.append(size)
    
    def error_bounds(bins, SFs, Stat, Sysup, Sysdown, smoothing = 1):
        # Calculate a smoothened set of error values for plotting a continuous band
        new_bins = []
        for i in range(len(bins)-1):
            new_bins.append(bins[i])
            for j in range(1,smoothing):
                new_bins.append(bins[i] + (bins[i+1]-bins[i])/float(smoothing)*j)
        new_bins.append(bins[-1])
        mids = [(new_bins[i]+new_bins[i+1])/2 for i in range(0, len(new_bins)-1)]
        xs = [x for temps in zip(new_bins, mids) for x in temps]
        xs.append(new_bins[-1])
        # New SF nominal values
        new_SFs = []
        for i in range(len(SFs)-1):
            new_SFs.append(SFs[i])
            for j in range(1,smoothing):
                new_SFs.append(SFs[i] + (SFs[i+1]-SFs[i])/float(smoothing)*j)
        new_SFs.append(SFs[-1])
        mid_SFs = [(new_SFs[i]+new_SFs[i+1])/2 for i in range(0, len(new_SFs)-1)]
        ys = [x for temps in zip(new_SFs, mid_SFs) for x in temps]
        ys.append(new_SFs[-1])
        ys = [SFs[0]]*smoothing+ys+[SFs[-1]]*smoothing

        # New error values        
        def sumunderroot(x):
            sig2 = 0
            for ix in x: sig2+=ix**2
            return np.sqrt(sig2)
        eyu = [sumunderroot(x)+SFs[i] for i,x in enumerate(zip(Stat, Sysup))]
        eyd = [SFs[i]-sumunderroot(x) for i,x in enumerate(zip(Stat, Sysdown))]
        new_EysU = []
        for i in range(len(eyu)-1):
            new_EysU.append(eyu[i])
            for j in range(1,smoothing):
                new_EysU.append(eyu[i] + (eyu[i+1]-eyu[i])/float(smoothing)*j)
        new_EysU.append(eyu[-1])
        mid_eyu = [(new_EysU[i]+new_EysU[i+1])/2 for i in range(0, len(new_EysU)-1)]
        eysU = [x for temps in zip(new_EysU, mid_eyu) for x in temps]
        eysU.append(new_EysU[-1])
        eysU = [eyu[0]]*smoothing+eysU+[eyu[-1]]*smoothing        
        new_EysD = []
        for i in range(len(eyd)-1):
            new_EysD.append(eyd[i])
            for j in range(1,smoothing):
                new_EysD.append(eyd[i] + (eyd[i+1]-eyd[i])/float(smoothing)*j)
        new_EysD.append(eyd[-1])
        mid_eyd = [(new_EysD[i]+new_EysD[i+1])/2 for i in range(0, len(new_EysD)-1)]
        eysD = [x for temps in zip(new_EysD, mid_eyd) for x in temps]
        eysD.append(new_EysD[-1])
        eysD = [eyd[0]]*smoothing+eysD+[eyd[-1]]*smoothing

        if smoothing == 1:
            xs = bins
            ys = SFs
            eysU = eyu
            eysD = eyd

        return xs, ys, eysU, eysD

    xs, ys, eyu, eyd = error_bounds(bins, SFs, Stat, Sysup, Sysdown)
    shapes = []
    #print xs
    for i in range(len(xs)*4):
        shapes.append(ROOT.TGraph())
    for shape in shapes:
        ROOT.SetOwnership(shape, False)
    for i in range(len(xs)-1):
        #print i
        for version in range(4):
            #print shapes[i+(len(xs))*version]
            #print eyu[i]
            shapes[i+(len(xs))*version].SetPoint((version+0)%4, xs[i],   eyu[i])
            shapes[i+(len(xs))*version].SetPoint((version+1)%4, xs[i+1],   eyu[i])
            shapes[i+(len(xs))*version].SetPoint((version+2)%4, xs[i+1],   eyd[i])
            shapes[i+(len(xs))*version].SetPoint((version+3)%4, xs[i],   eyd[i])
 
    for shape in shapes:
        shape.SetFillColor(kBlue-4)
        shape.SetFillStyle(3004)
 
    # Put into arrays for ROOT
    xs = array('f', xs)
    ys = array('f', ys)
    #exl = array('f', exl)
    #exr = array('f', exr)
    eyu = array('f', eyu)
    eyd = array('f', eyd)
    SFs = array('f',SFs)   
    sys_up = array('f', Sysup) 
    sys_down = array('f',Sysdown)  
    stat = array('f',Stat)
    bin_center = array('f',bin_center)
    bin_size = array('f',bin_size) 

    #print xs
    #print ys
    #print exl, exr
    #print eyu, eyd

    # Style
    gROOT.SetBatch()    
    gROOT.SetStyle("Plain")
    gStyle.SetHistTopMargin(0)
    gStyle.SetOptStat(0)       
    # Set custom style and make activate
    style = set_style()
    style.cd()
      
    # Start plotting
    c1 = ROOT.TCanvas("c1","c1",0,0,600,500)
    c1.Draw()
    c1.cd()
    c1.SetLogx(0)
    c1.SetGrid(1)

    plot = ROOT.TGraphAsymmErrors(len(bin_center), bin_center, SFs, bin_size, bin_size, sys_down, sys_up)   
    
    plot.SetMarkerSize(1)
    plot.SetMarkerStyle(2)
    plot.SetLineWidth(3)
    plot.SetMarkerColor(46)
    plot.SetLineColor(46)
    plot.SetTitle("")
    plot.GetXaxis().SetMoreLogLabels()
    plot.GetXaxis().SetNoExponent()
    plot.SetMinimum(0.6)
    plot.SetMaximum(1.4)
    plot.GetXaxis().SetTitle("p_{T} [GeV]")
    plot.GetYaxis().SetTitle("SF_{double b}")

    # Draw
    plot.Draw("AP")
    for shape in shapes:
        shape.Draw("f same")    

    #leg = ROOT.TLegend(0.35,0.35,0.70,0.20)
    leg = ROOT.TLegend(0.55,0.7,0.9,0.9)
    leg.SetFillColor(253)    
    leg.SetBorderSize(0)
    leg.AddEntry(plot, "SF (systematics)", "PL")
    leg.AddEntry(shape, "Stat+syst", "F")
    tag = ROOT.TLatex()
    tag.SetNDC();
    tag.SetTextAlign(22);
    tag.SetTextFont(63);
    tag.SetTextSizePixels(30);
    #tag.DrawLatex(0.55,0.82,"DoubleB"+WP)
    tag.DrawLatex(0.35,0.82,"DoubleB"+WP)
    tag2 = ROOT.TLatex()
    tag2.SetNDC();
    tag2.SetTextAlign(22);
    tag2.SetTextFont(63);
    tag2.SetTextSizePixels(14);
    #tag2.DrawLatex(0.55,0.76,"(Preliminary)")
    tag2.DrawLatex(0.35,0.76,"(Preliminary)")
    leg.Draw()       
    
    c1.RedrawAxis("g") 

    c1.Print("pics/SFComp_cMVAv2"+WP+".png")    
    c1.Clear()  

    #WP, SFs, Stat, Sysup, Sysdown
    for bin in range(len(SFs)):
        print WP, "bin", bin,  SFs[bin], "+/-", np.sqrt(Stat[bin]**2+Sysup[bin]**2), "/", np.sqrt(Stat[bin]**2+Sysdown[bin]**2)


def set_style():
    style = TStyle("STYLE","User style")
    icol=0

    style.SetFrameBorderMode(icol)
    style.SetFrameFillColor(icol)
    style.SetCanvasBorderMode(icol)
    style.SetCanvasColor(icol)
    style.SetPadBorderMode(icol)
    style.SetPadColor(icol)
    style.SetStatColor(icol)
    
    style.SetPaperSize(20,26) 
    style.SetPadTopMargin(0.05)
    style.SetPadRightMargin(0.05)
    style.SetPadBottomMargin(0.16)
    style.SetPadLeftMargin(0.16)

    style.SetTitleXOffset(1.4)
    style.SetTitleYOffset(1.4)

    font=42; # Helvetica
    tsize=0.05;
    style.SetTextFont(font)

    style.SetTextSize(tsize);
    style.SetLabelFont(font,"x")
    style.SetTitleFont(font,"x")
    style.SetLabelFont(font,"y")
    style.SetTitleFont(font,"y")
    style.SetLabelFont(font,"z")
    style.SetTitleFont(font,"z")  
    style.SetLabelSize(tsize,"x")
    style.SetTitleSize(tsize,"x")
    style.SetLabelSize(tsize,"y")
    style.SetTitleSize(tsize,"y")
    style.SetLabelSize(tsize,"z")
    style.SetTitleSize(tsize,"z")

    style.SetMarkerStyle(20)
    style.SetMarkerSize(1.2)
    #style.SetHistLineWidth(2.);
    style.SetLineStyleString(2,"[12 12]") 
    # get rid of X error bars 
    #style.SetErrorX(0.001);
    # get rid of error bar caps
    style.SetEndErrorSize(0.)

    # do not display any of the standard histogram decorations
    style.SetOptTitle(0)
    #style.SetOptStat(1111);
    style.SetOptStat(0)
    #style.SetOptFit(1111);
    style.SetOptFit(0)

    #put tick marks on top and RHS of plots
    style.SetPadTickX(1)
    style.SetPadTickY(1)
    return style


if __name__ == "__main__":
    var = 1
    # DoubleBL
    #SFs, Stat, Sysup, Sysdown = [1, 0.99, 1.01], [0.014, 0.032, 0.023], [0.015, 0.03, 0.018], [0.013, 0.017, 0.01] 
    #DoubleBM1
    #SFs, Stat, Sysup, Sysdown = [0.95, 0.97, 0.95], [0.016, 0.038, 0.027], [0.011, 0.022, 0.015], [0.01, 0.004, 0.01] 
    #DoubleBM2
    #SFs, Stat, Sysup, Sysdown = [0.9, 0.92, 0.88], [0.019, 0.048, 0.033], [0.013, 0.018, 0.011 ], [0., 0.011, 0.01 ] 
    #DoubleBH
    #SFs, Stat, Sysup, Sysdown = [0.82, 0.84, 0.78], [0.028, 0.091, 0.035], [0.,  0.027, 0.01], [0.003, 0.017, 0.] 

    # 2018 STAT
    #List of lists
    #BL, SFs, Stat, Sysup, Sysdown
    #BL = ["L", [1.02, 0.99, 1], [0.017, 0.032, 0.022], [0.01, 0.079, 0.013], [0.02, 0.01, 0.]]
    #BM1
    #BM1 = ["M1", [0.97, 0.94, 0.97], [0.02, 0.037, 0.025], [0.011,0.059,0.013],[0.017, 0.004, 0.014]]
    #BM2
    #BM2 = ["M2", [0.91, 0.92, 0.9], [0.022, 0.042, 0.032], [0.009, 0.043, 0.009],[0.012, 0.013, 0.01]]
    #BT
    #BT = ["T", [0.84, 0.84, 0.79], [0.026,0.051,0.034], [0.011, 0.029, 0.01], [0.006, 0.02, 0.011]]
    
    """
    # Partial 2017 Dataset, Full systematics
    BL = ['L', [1.01, 0.98999999999999999, 1.0], [0.016, 0.034000000000000002, 0.023], [0.025000000000000001, 0.10000000000000001, 0.019], [0.029000000000000001, 0.037999999999999999, 0.043999999999999997]]
    BM1 = ['M1', [0.95999999999999996, 0.94999999999999996, 0.96999999999999997], [0.019, 0.037999999999999999, 0.027], [0.048000000000000001, 0.072999999999999995, 0.017000000000000001], [0.024, 0.040000000000000001, 0.037999999999999999]]
    BM2 = ['M2', [0.90000000000000002, 0.92000000000000004, 0.90000000000000002], [0.021999999999999999, 0.043999999999999997, 0.032000000000000001], [0.033000000000000002, 0.058999999999999997, 0.041000000000000002], [0.02, 0.017999999999999999, 0.052999999999999999]]
    BH = ['H', [0.82999999999999996, 0.83999999999999997, 0.79000000000000004], [0.024, 0.052999999999999999, 0.035000000000000003], [0.031, 0.034000000000000002, 0.010999999999999999], [0.019, 0.032000000000000001, 0.021000000000000001]]
    
    """
    # Full (almost) dataset
    """
    BL = ['L', [1.0, 0.97999999999999998, 1.01], [0.017000000000000001, 0.027, 0.021999999999999999], [0.057000000000000002, 0.109, 0.021999999999999999], [0.033000000000000002, 0.032000000000000001, 0.059999999999999998]]

    BM1 = ['M1', [0.96999999999999997, 0.95999999999999996, 0.96999999999999997], [0.019, 0.031, 0.025999999999999999], [0.029999999999999999, 0.023, 0.027], [0.025000000000000001, 0.02, 0.041000000000000002]]

    BM2 = ['M2', [0.90000000000000002, 0.900000000000000021000000000000003, 0.87], [0.028000000000000001, 0.042000000000000003, 0.031], [0.021000000000000001, 0.021999999999999999, 0.023], [0.027, 0.023, 0.045999999999999999]]

    BH = ['H', [0.83999999999999997, 0.83999999999999997, 0.78000000000000003], [0.024, 0.049000000000000002, 0.029999999999999999], [0.014999999999999999, 0.042000000000000003, 0.021999999999999999], [0.016, 0.039, 0.028000000000000001]]
    """
    """
    # SV syst - NTracks, CFRAG
    BL = ['L', [0.95999999999999996, 1.05, 0.92000000000000004], [0.027, 0.065000000000000002, 0.040000000000000001], [0.014, 0.016, 0.029000000000000001], [0.012999999999999999, 0.064000000000000001, 0.014]]
    BM1 = ['M1', [0.93000000000000005, 1.01, 0.83999999999999997], [0.033000000000000002, 0.070999999999999994, 0.042999999999999997], [0.014, 0.019, 0.036999999999999998], [0.012999999999999999, 0.062, 0.012999999999999999]]
    BM2 = ['M2', [0.84999999999999998, 0.93000000000000005, 0.78000000000000003], [0.033000000000000002, 0.067000000000000004, 0.044999999999999998], [0.01, 0.017000000000000001, 0.037999999999999999], [0.016, 0.065000000000000002, 0.016]]
    BH = ['H', [0.78000000000000003, 0.84999999999999998, 0.67000000000000004], [0.037999999999999999, 0.086999999999999994, 0.044999999999999998], [0.012999999999999999, 0.033000000000000002, 0.019], [0.021000000000000001, 0.050999999999999997, 0.02]]
    
    # SV - all syst
    BL = ['L', [0.92000000000000004, 1.0600000000000001, 0.93000000000000005], [0.048000000000000001, 0.070999999999999994, 0.043999999999999997], [0.10299999999999999, 0.021999999999999999, 0.023], [0.023, 0.087999999999999995, 0.029000000000000001]]
    BM1 = ['M1', [0.87, 1.03, 0.84999999999999998], [0.050000000000000003, 0.079000000000000001, 0.047], [0.17999999999999999, 0.031, 0.043999999999999997], [0.024, 0.126, 0.02]]
    BM2 = ['M2', [0.79000000000000004, 0.95999999999999996, 0.79000000000000004], [0.051999999999999998, 0.075999999999999998, 0.049000000000000002], [0.17399999999999999, 0.047, 0.035999999999999997], [0.021999999999999999, 0.156, 0.033000000000000002]]
    BH = ['T', [0.71999999999999997, 0.90000000000000002, 0.67000000000000004], [0.051999999999999998, 0.098000000000000004, 0.047], [0.184, 0.051999999999999998, 0.021999999999999999], [0.037999999999999999, 0.157, 0.027]]
    """
    # SV adjusted shape systs
    BL = ['L', [0.92000000000000004, 1.0600000000000001, 0.93000000000000005], [0.048000000000000001, 0.070999999999999994, 0.043999999999999997], [0.042999999999999997, 0.025000000000000001, 0.032000000000000001], [0.028000000000000001, 0.055, 0.021000000000000001]]
    BM1 = ['M1', [0.87, 1.03, 0.84999999999999998], [0.050000000000000003, 0.079000000000000001, 0.047], [0.069000000000000006, 0.032000000000000001, 0.044999999999999998], [0.028000000000000001, 0.066000000000000003, 0.02]]
    BM2 = ['M2', [0.79000000000000004, 0.95999999999999996, 0.79000000000000004], [0.051999999999999998, 0.075999999999999998, 0.049000000000000002], [0.073999999999999996, 0.060999999999999999, 0.042999999999999997], [0.021999999999999999, 0.067000000000000004, 0.023]]
    BH = ['H', [0.71999999999999997, 0.90000000000000002, 0.67000000000000004], [0.051999999999999998, 0.098000000000000004, 0.047], [0.070000000000000007, 0.063, 0.029000000000000001], [0.037999999999999999, 0.065000000000000002, 0.02]]    
    #SFComp("T", SFs, Stat, Sysup, Sysdown)
    for lis in [BL, BM1, BM2, BH]:
        SFComp(lis[0], lis[1], lis[2], lis[3], lis[4])
        for i, p in enumerate(["250-350","350-430","430-840"]):
            #a,b,c = lis[0], lis[1][i], "+/-", np.sqrt(lis[2][i]**2+lis[3][i]**2), "/", np.sqrt(lis[2][i]**2+lis[4][i]**2)
            a,b,c = lis[1][i], np.sqrt(lis[2][i]**2+lis[3][i]**2), np.sqrt(lis[2][i]**2+lis[4][i]**2)
            print "|||||"+ p+" |"+str(round(a,2)) +"&plusmn;"+str(round(b,2))+"/"+str(round(c,2))+"|"
    	
