//----------------------Version 1.0-------------------------//
//Plotting Macro for Mu -> Tau Fake Rates study
//Author: Yiwen Wen
//DESY
//-----------------------------------------------------------//
#include "HttStylesNew.cc"
#include "TColor.h"
#include <math.h>
void ProduceDatacardInputs_MuTauFR_m_vis_Scale_ZMM(
                                 TString varName = "m_vis_reso_scale_up",
                                 TString uncertSuffix = "reso_Up",
                                 TString etaSuffix = "Lt0p4",
                                 TString wp = "taubyLooseDeepTau2017v2p1VSmu",
                                 float xLower = 70,
                                 float xUpper = 120,
                                 int numberofbins = 10,
                                 bool passProbe = true,
                                 bool deepTau = true,
                                 TString extraCuts = ""
                                 )
{
	SetStyle();
	TString samples[1] =
                    {
                        "NTuples_DYJetsToLL_M-50"// (0)Drell-Yan Z->MuMu
                    };

	float xsec[1] = {
                        6077.22,// (0)Drell-Yan Z->MuMu
                    };
    TString suffix("");
    if(passProbe)
        suffix="_pass";
    else
        suffix="_fail";
    
    
	float lumi = 59740;
	//Deal with bins and binning
	float xMin = xLower;
	float xMax = xUpper;
	int nBins = numberofbins;
	float bins[100];
	float binWidth = (xMax-xMin)/float(nBins);
        for (int iB=0; iB<=nBins; ++iB)
        	bins[iB] = xMin + float(iB)*binWidth;

	int nSamples = 1;
	TH1D * hist[1];
    TString whatTauid = "DecayModeFindingProbe>0.5 &&taubyTightIsolationMVArun2017v2DBoldDMwLT2017";
    if (deepTau)
    {
        whatTauid = "DecayModeFindingNewProbe>0.5 && taubyMediumDeepTau2017v2p1VSjet";
    }

    // *******************************
    // ***** DY+Jets samples *******
    // *******************************
    TH1D * histZMM[9];
    
    TString dyRefSamples[5] = {"NTuples_DYJetsToLL_M-50",
        "NTuples_DY1JetsToLL_M-50",
        "NTuples_DY2JetsToLL_M-50",
        "NTuples_DY3JetsToLL_M-50",
        "NTuples_DY4JetsToLL_M-50"
    };
    double dyRefXSec[5] = {6077.22,
        1.137*877.8,
        1.137*304.4,
        1.137*111.5,
        1.137*44.03};//2017 values
    double dyRefEvents[5] = {97800939,34859434,9790490,6897933,4346952};
    
    for (int iDY=0; iDY<5; ++iDY) {
        TFile * file = new TFile("/home/cmsdas/public/store/TausShortExercise/"+dyRefSamples[iDY]+".root");
        TH1D * histWeightsH = (TH1D*)file->Get("histWeightsH");
        dyRefEvents[iDY] = histWeightsH->GetSumOfWeights();
    }
    TString dySampleNames[9] = {"NTuples_DYJetsToLL_M-50",
			     "NTuples_DYJetsToLL_M-50",
			     "NTuples_DYJetsToLL_M-50",
			     "NTuples_DYJetsToLL_M-50",
			     "NTuples_DYJetsToLL_M-50",
        "NTuples_DY1JetsToLL_M-50",
        "NTuples_DY2JetsToLL_M-50",
        "NTuples_DY3JetsToLL_M-50",
        "NTuples_DY4JetsToLL_M-50"
    };
    double dyNorm[9];
    dyNorm[0] = lumi*dyRefXSec[0]/dyRefEvents[0];
    dyNorm[1] = lumi/(dyRefEvents[0]/dyRefXSec[0]+dyRefEvents[1]/dyRefXSec[1]);
    dyNorm[2] = lumi/(dyRefEvents[0]/dyRefXSec[0]+dyRefEvents[2]/dyRefXSec[2]);
    dyNorm[3] = lumi/(dyRefEvents[0]/dyRefXSec[0]+dyRefEvents[3]/dyRefXSec[3]);
    dyNorm[4] = lumi/(dyRefEvents[0]/dyRefXSec[0]+dyRefEvents[4]/dyRefXSec[4]);
    dyNorm[5] = lumi/(dyRefEvents[0]/dyRefXSec[0]+dyRefEvents[1]/dyRefXSec[1]);
    dyNorm[6] = lumi/(dyRefEvents[0]/dyRefXSec[0]+dyRefEvents[2]/dyRefXSec[2]);
    dyNorm[7] = lumi/(dyRefEvents[0]/dyRefXSec[0]+dyRefEvents[3]/dyRefXSec[3]);
    dyNorm[8] = lumi/(dyRefEvents[0]/dyRefXSec[0]+dyRefEvents[4]/dyRefXSec[4]);
    
    TString npartonCuts[9] = {"&&(npartons==0||npartons>4)",
        "&&npartons==1",
        "&&npartons==2",
        "&&npartons==3",
        "&&npartons==4",
        "",
        "",
        "",
        ""
    };
    TString cutsZMM[9];
    for (int iDY=0; iDY<9; ++iDY)
    {
        cutsZMM[iDY]   = "zptweight*puweight*effweight*mcweight*(mt_1<40&&os>0.5&& gen_match_1 == 2 && gen_match_2 == 2 &&iso_1<0.15"+npartonCuts[iDY]+")";
        
    }
    if(passProbe)
    {
        for (int i=0; i<9; ++i)
        {
            cutsZMM[i] ="("+wp+">0.5&&"+whatTauid+">0.5)*"+cutsZMM[i];
        }
        
    }
    else
    {
        for (int i=0; i<9; ++i)
        {
            cutsZMM[i] ="("+wp+"<0.5&&"+whatTauid+">0.5)*"+cutsZMM[i];
            
        }
    }
    
    if(etaSuffix=="Lt0p4")
    {
        for (int i=0; i<9; ++i)
        {
            cutsZMM[i] = "(fabs(EtaProbe)<0.4)*"+cutsZMM[i];
        }
    }
    
    if(etaSuffix=="0p4to0p8")
    {
        for (int i=0; i<9; ++i)
        {
            cutsZMM[i] = "(fabs(EtaProbe)<0.8 && fabs(EtaProbe)>0.4)*"+cutsZMM[i];
        }
    }
    
    if(etaSuffix=="0p8to1p2")
    {
        for (int i=0; i<9; ++i)
        {
            cutsZMM[i] = "(fabs(EtaProbe)<1.2 && fabs(EtaProbe)>0.8)*"+cutsZMM[i];
        }
    }
    
    if(etaSuffix=="1p2to1p7")
    {
        for (int i=0; i<9; ++i)
        {
            cutsZMM[i] = "(fabs(EtaProbe)<1.7 && fabs(EtaProbe)>1.2)*"+cutsZMM[i];
        }
    }
    if(etaSuffix=="Gt1p7")
    {
        for (int i=0; i<9; ++i)
        {
            cutsZMM[i] = "(fabs(EtaProbe)>1.7)*"+cutsZMM[i];
        }
    }
    // filling histograms for DYJets samples
    for (int i=0; i<9; ++i) { // run over samples
        
        TFile * file = new TFile("/home/cmsdas/public/store/TausShortExercise/"+dySampleNames[i]+".root");
        TTree * tree = (TTree*)file->Get("MuTauFR");
        double normaliza = dyNorm[i];
        TString histNameZMM = dySampleNames[i] + "_"+varName+"_ZMM";
        histZMM[i] = new TH1D(histNameZMM,"",nBins,xMin,xMax);
        histZMM[i]->Sumw2();
        
        tree->Draw(varName+">>"+histNameZMM,cutsZMM[i] + extraCuts);
        
        for (int iB=1; iB<=nBins; ++iB) {
            double x = histZMM[i]->GetBinContent(iB);
            double e = histZMM[i]->GetBinError(iB);
            histZMM[i]->SetBinContent(iB,normaliza*x);
            histZMM[i]->SetBinError(iB,normaliza*e);
        }
        std::cout << dySampleNames[i] << " -> DY ZMM = " << histZMM[i]->GetEntries() << " : " << histZMM[i]->GetSumOfWeights() << std::endl;
    }
    
    hist[0]  = histZMM[0];
    for (int iDY=1; iDY<9; ++iDY)
    {
        hist[0]->Add(hist[0],histZMM[iDY]);
    }
    
    std::cout << "end: DY stitching" << std::endl;
    
    TFile * file = new TFile("./shapes/MuTauFR"+wp+etaSuffix+suffix+"_"+varName+".root","recreate");
    file->mkdir("MuTauFR"+suffix);
    file->cd("MuTauFR"+suffix);
    TH1D * ZMM = (TH1D*)hist[0]->Clone("ZMM_"+uncertSuffix);

    float totZMM = ZMM->GetSumOfWeights();

    std::cout << "ZMM : " << ZMM->GetSumOfWeights() << std::endl;
    
    file->Write();
    file->Close();


}

