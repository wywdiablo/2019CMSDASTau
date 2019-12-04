//----------------------Version 1.2-------------------------//
//Producing Shapes Macro for CMS Mu -> Tau Fake Rates Study
//Author: Yiwen Wen
//DESY
//----------------------------------------------------------//
//v1.1: ZJ and ZTT template
//v1.2: Stitched WJets samples
//v1.3: Apply veto events with bad/duplicate muon 
//----------------------------------------------------------//
#include "HttStylesNew.cc"
#include "TColor.h"
#include <math.h>
void ProduceDatacardInputs_MuTauFR(
                                 TString varName = "m_vis",
                                 TString etaSuffix = "Lt0p4",
                                 TString wp = "taubyLooseDeepTau2017v2p1VSmu",
                                 float xLower = 70,
                                 float xUpper = 120,
                                 int numberofbins = 10,
                                 bool passProbe = true,
                                 bool deepTau = true,
                                 TString extraCuts = "*tauidweight"
                                 )
{
	SetStyle();
    TString samples[11] =
    {
        "NTuples_SingleMuon_Run2018ABCD",//(0)data
        "NTuples_TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8",//(1)TTbar to 2L2Nu
        "NTuples_TTToHadronic_TuneCP5_13TeV-powheg-pythia8", //(2) TTbar to hadronic
        "NTuples_TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8", //(3) TTbar to semileptonic
        "NTuples_WW_TuneCP5_13TeV-pythia8",// (4)WW
        "NTuples_WZ_TuneCP5_13TeV-pythia8",// (5)WZ
        "NTuples_ZZ_TuneCP5_13TeV-pythia8",// (6)ZZ
        "NTuples_ST_tW_antitop_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8", // (7) SingleTop tW tbar
        "NTuples_ST_tW_top_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8", // (8) SingleTop tW t
        "NTuples_ST_t-channel_antitop_4f_InclusiveDecays_TuneCP5_13TeV-powheg-madspin-pythia8",// (9) SingleTop t antitop
        "NTuples_ST_t-channel_top_4f_InclusiveDecays_TuneCP5_13TeV-powheg-madspin-pythia8"// (10) SingleTop t top
    };
    
    float xsec[11] = {
        1, // (0)data
        88.29, // (1)TT to 2L2Nu
        365.35, // (2)TT to hadronic
        377.96, // (3)TT to semileptonic
        63.21,  // WW (4)
        22.82,   // WZ (5)
        10.32,  // ZZ (6)
        35.06,   // ST_tW_antitop_5f_inclusiveDecays (7)
        35.06,   // ST_tW_top_5f_inclusiveDecays (8)
        80.95,   // ST_t_antitop (9)
        136.02   // ST_t_top (10)
    };
    /*
    double eventNum[16] = {
        1,
        1, //TT to 2L2Nu
        1, //TT to hadronic
        1, //TT to semileptonic
        7791498, //WW
        3928630, //WZ
        1949768, //ZZ
        1, //ST_tW_antitop
        1, //ST_tW_top
        1, //ST_t_antitop
        1, //ST_t_top
    };*/
    
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

	int nSamples = 11;
	TH1D * hist[16];
	TH1D * histSS[16];

    
    TString whatTauid = "DecayModeFindingProbe>0.5 &&taubyTightIsolationMVArun2017v2DBoldDMwLT2017";
    if (deepTau)
    {
        whatTauid = "DecayModeFindingNewProbe>0.5 && taubyMediumDeepTau2017v2p1VSjet";
    }
	//inistiating cuts
	TString cuts[15];
	TString cutsSS[15];
    
    TString WJetsweight="0.950578*";
    TString qcdweight ="1.06*";
    
	for (int i=0; i<nSamples; ++i)
	{
        cuts[i] = "puweight*effweight*mcweight*(mt_1<40&&os>0.5&&iso_1<0.15)";
        cutsSS[i] = qcdweight+"puweight*effweight*mcweight*(mt_1<40&&os<0.5&&iso_1<0.15)";
  	}

    cuts[0] = "(mt_1<40 && os>0.5 &&iso_1<0.15)";
    cutsSS[0] = qcdweight+"(mt_1<40 && os<0.5 && iso_1<0.15)";

    if(passProbe)
    {
        for (int i=0; i<nSamples; ++i)
        {
            cuts[i] ="("+wp+">0.5&&"+whatTauid+">0.5)*"+cuts[i];
            cutsSS[i] ="("+wp+">0.5&&"+whatTauid+">0.5)*"+cutsSS[i];
        }
        
    }
    else
    {
        for (int i=0; i<nSamples; ++i)
        {
            cuts[i] ="("+wp+"<0.5 &&"+whatTauid+">0.5)*"+cuts[i];
            cutsSS[i] ="("+wp+"<0.5 &&"+whatTauid+">0.5)*"+cutsSS[i];
        }
    }
    
    if(etaSuffix=="Lt0p4")
    {
        for (int i=0; i<nSamples; ++i)
        {
            cuts[i] = "(fabs(EtaProbe)<0.4)*"+cuts[i];
            cutsSS[i] = "(fabs(EtaProbe)<0.4)*"+cutsSS[i];
        }
    }
    
    if(etaSuffix=="0p4to0p8")
    {
        for (int i=0; i<nSamples; ++i)
        {
            cuts[i] = "(fabs(EtaProbe)<0.8 && fabs(EtaProbe)>0.4)*"+cuts[i];
            cutsSS[i] = "(fabs(EtaProbe)<0.8 && fabs(EtaProbe)>0.4)*"+cutsSS[i];
        }
    }
    
    if(etaSuffix=="0p8to1p2")
    {
        for (int i=0; i<nSamples; ++i)
        {
            cuts[i] = "(fabs(EtaProbe)<1.2 && fabs(EtaProbe)>0.8)*"+cuts[i];
            cutsSS[i] = "(fabs(EtaProbe)<1.2 && fabs(EtaProbe)>0.8)*"+cutsSS[i];
        }
    }
    
    if(etaSuffix=="1p2to1p7")
    {
        for (int i=0; i<nSamples; ++i)
        {
            cuts[i] ="(fabs(EtaProbe)<1.7 && fabs(EtaProbe)>1.2)*"+cuts[i];
            cutsSS[i] ="(fabs(EtaProbe)<1.7 && fabs(EtaProbe)>1.2)*"+cutsSS[i];
        }
    }
    if(etaSuffix=="Gt1p7")
    {
        for (int i=0; i<nSamples; ++i)
        {
            cuts[i] ="(fabs(EtaProbe)>1.7)*"+cuts[i];
            cutsSS[i] ="(fabs(EtaProbe)>1.7)*"+cutsSS[i];
        }
    }
    
	for (int i=0; i<nSamples; ++i)
	{
		std::cout <<i<<":"+samples[i]<<std::endl;
		TFile * file = new TFile("./NTuples/"+samples[i]+".root");
		TH1D * histWeightsH = (TH1D*)file->Get("histWeightsH");
		TTree * tree = (TTree*)file->Get("MuTauFR");
        double normaliza = xsec[i]*lumi/histWeightsH->GetSumOfWeights();
		TString histName = samples[i] + "_"+varName;
        TString histNameSS = samples[i] + "_"+varName+"_ss";
        hist[i] = new TH1D(histName,"",nBins,xMin,xMax);
		hist[i]->Sumw2();
        histSS[i] = new TH1D(histNameSS,"",nBins,xMin,xMax);
        histSS[i]->Sumw2();
        tree->Draw(varName+">>"+histName,cuts[i]+extraCuts);
        tree->Draw(varName+">>"+histNameSS,cutsSS[i]+extraCuts);
		if(i > 0)
		{	
			for (int iB=1; iB<=nBins; ++iB) 
			{
                double x = hist[i]->GetBinContent(iB);
                double e = hist[i]->GetBinError(iB);
                hist[i]->SetBinContent(iB,normaliza*x);
                hist[i]->SetBinError(iB,normaliza*e);
                double xSS = histSS[i]->GetBinContent(iB);
                double eSS = histSS[i]->GetBinError(iB);
                histSS[i]->SetBinContent(iB,normaliza*xSS);
                histSS[i]->SetBinError(iB,normaliza*eSS);
            }
		}	
	}
    // *******************************
    // ***** DY+Jets samples *******
    // *******************************
    TH1D * histZMM[9];
    TH1D * histZMMSS[9];
    TH1D * histZJ[9];
    TH1D * histZJSS[9];
    TH1D * histZTTml[9];
    TH1D * histZTTmlSS[9];
    TH1D * histZTTmt[9];
    TH1D * histZTTmtSS[9];
    
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
        TFile * file = new TFile("./NTuples/"+dyRefSamples[iDY]+".root");
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
    TString cutsZMMSS[9];
    TString cutsZJ[9];
    TString cutsZJSS[9];
    TString cutsZTTml[9];
    TString cutsZTTmlSS[9];
    TString cutsZTTmt[9];
    TString cutsZTTmtSS[9];
    for (int iDY=0; iDY<9; ++iDY)
    {
        cutsZMM[iDY]   = "zptweight*puweight*effweight*mcweight*(mt_1<40&&os>0.5&&gen_match_1 == 2 && gen_match_2 == 2  &&iso_1<0.15"+npartonCuts[iDY]+")";
        cutsZMMSS[iDY] = qcdweight+"zptweight*puweight*effweight*mcweight*(mt_1<40&&os<0.5&&gen_match_1 == 2 && gen_match_2 == 2  &&iso_1<0.15"+npartonCuts[iDY]+")";
        cutsZJ[iDY]   = "zptweight*puweight*effweight*mcweight*(mt_1<40&&os>0.5&& gen_match_2 == 6 &&iso_1<0.15"+npartonCuts[iDY]+")";
        cutsZJSS[iDY] = qcdweight+"zptweight*puweight*effweight*mcweight*(mt_1<40&&os<0.5&& gen_match_2 == 6 &&iso_1<0.15"+npartonCuts[iDY]+")";
        cutsZTTml[iDY]   = "zptweight*puweight*effweight*mcweight*(mt_1<40&&os>0.5&& gen_match_1 == 4 && (gen_match_2 == 3 || gen_match_2 == 4) &&iso_1<0.15"+npartonCuts[iDY]+")";
        cutsZTTmlSS[iDY] = qcdweight+"zptweight*puweight*effweight*mcweight*(mt_1<40&&os<0.5&& gen_match_1 == 4 && (gen_match_2 == 3 || gen_match_2 == 4) &&iso_1<0.15"+npartonCuts[iDY]+")";
        cutsZTTmt[iDY]   = "zptweight*puweight*effweight*mcweight*(mt_1<40&&os>0.5&& gen_match_1 == 4 && gen_match_2 == 5 &&iso_1<0.15"+npartonCuts[iDY]+")";
        cutsZTTmtSS[iDY] = qcdweight+"zptweight*puweight*effweight*mcweight*(mt_1<40&&os<0.5&& gen_match_1 == 4 && gen_match_2 == 5 &&iso_1<0.15"+npartonCuts[iDY]+")";
        
    }
    if(passProbe)
    {
        for (int i=0; i<9; ++i)
        {
            cutsZMM[i] ="("+wp+">0.5&&"+whatTauid+">0.5)*"+cutsZMM[i];
            cutsZMMSS[i] ="("+wp+">0.5&&"+whatTauid+">0.5)*"+cutsZMMSS[i];
            cutsZJ[i] ="("+wp+">0.5&&"+whatTauid+">0.5)*"+cutsZJ[i];
            cutsZJSS[i] ="("+wp+">0.5&&"+whatTauid+">0.5)*"+cutsZJSS[i];
            cutsZTTml[i] ="("+wp+">0.5&&"+whatTauid+">0.5)*"+cutsZTTml[i];
            cutsZTTmlSS[i] ="("+wp+">0.5&&"+whatTauid+">0.5)*"+cutsZTTmlSS[i];
            cutsZTTmt[i] ="("+wp+">0.5&&"+whatTauid+">0.5)*"+cutsZTTmt[i];
            cutsZTTmtSS[i] ="("+wp+">0.5&&"+whatTauid+">0.5)*"+cutsZTTmtSS[i];
        }
        
    }
    else
    {
        for (int i=0; i<9; ++i)
        {
            cutsZMM[i] ="("+wp+"<0.5&&"+whatTauid+">0.5)*"+cutsZMM[i];
            cutsZMMSS[i] ="("+wp+"<0.5&&"+whatTauid+">0.5)*"+cutsZMMSS[i];
            cutsZJ[i] ="("+wp+"<0.5&&"+whatTauid+">0.5)*"+cutsZJ[i];
            cutsZJSS[i] ="("+wp+"<0.5&&"+whatTauid+">0.5)*"+cutsZJSS[i];
            cutsZTTml[i] ="("+wp+"<0.5&&"+whatTauid+">0.5)*"+cutsZTTml[i];
            cutsZTTmlSS[i] ="("+wp+"<0.5&&"+whatTauid+">0.5)*"+cutsZTTmlSS[i];
            cutsZTTmt[i] ="("+wp+"<0.5&&"+whatTauid+">0.5)*"+cutsZTTmt[i];
            cutsZTTmtSS[i] ="("+wp+"<0.5&&"+whatTauid+">0.5)*"+cutsZTTmtSS[i];
            
        }
    }
    
    if(etaSuffix=="Lt0p4")
    {
        for (int i=0; i<9; ++i)
        {
            cutsZMM[i] = "(fabs(EtaProbe)<0.4)*"+cutsZMM[i];
            cutsZMMSS[i] = "(fabs(EtaProbe)<0.4)*"+cutsZMMSS[i];
            cutsZJ[i] = "(fabs(EtaProbe)<0.4)*"+cutsZJ[i];
            cutsZJSS[i] = "(fabs(EtaProbe)<0.4)*"+cutsZJSS[i];
            cutsZTTml[i] = "(fabs(EtaProbe)<0.4)*"+cutsZTTml[i];
            cutsZTTmlSS[i] = "(fabs(EtaProbe)<0.4)*"+cutsZTTmlSS[i];
            cutsZTTmt[i] = "(fabs(EtaProbe)<0.4)*"+cutsZTTmt[i];
            cutsZTTmtSS[i] = "(fabs(EtaProbe)<0.4)*"+cutsZTTmtSS[i];
        }
    }
    
    if(etaSuffix=="0p4to0p8")
    {
        for (int i=0; i<9; ++i)
        {
            cutsZMM[i] = "(fabs(EtaProbe)<0.8 && fabs(EtaProbe)>0.4)*"+cutsZMM[i];
            cutsZMMSS[i] = "(fabs(EtaProbe)<0.8 && fabs(EtaProbe)>0.4)*"+cutsZMMSS[i];
            cutsZJ[i] = "(fabs(EtaProbe)<0.8 && fabs(EtaProbe)>0.4)*"+cutsZJ[i];
            cutsZJSS[i] = "(fabs(EtaProbe)<0.8 && fabs(EtaProbe)>0.4)*"+cutsZJSS[i];
            cutsZTTml[i] = "(fabs(EtaProbe)<0.8 && fabs(EtaProbe)>0.4)*"+cutsZTTml[i];
            cutsZTTmlSS[i] = "(fabs(EtaProbe)<0.8 && fabs(EtaProbe)>0.4)*"+cutsZTTmlSS[i];
            cutsZTTmt[i] = "(fabs(EtaProbe)<0.8 && fabs(EtaProbe)>0.4)*"+cutsZTTmt[i];
            cutsZTTmtSS[i] = "(fabs(EtaProbe)<0.8 && fabs(EtaProbe)>0.4)*"+cutsZTTmtSS[i];
        }
    }
    
    if(etaSuffix=="0p8to1p2")
    {
        for (int i=0; i<9; ++i)
        {
            cutsZMM[i] = "(fabs(EtaProbe)<1.2 && fabs(EtaProbe)>0.8)*"+cutsZMM[i];
            cutsZMMSS[i] = "(fabs(EtaProbe)<1.2 && fabs(EtaProbe)>0.8)*"+cutsZMMSS[i];
            cutsZJ[i] = "(fabs(EtaProbe)<1.2 && fabs(EtaProbe)>0.8)*"+cutsZJ[i];
            cutsZJSS[i] = "(fabs(EtaProbe)<1.2 && fabs(EtaProbe)>0.8)*"+cutsZJSS[i];
            cutsZTTml[i] = "(fabs(EtaProbe)<1.2 && fabs(EtaProbe)>0.8)*"+cutsZTTml[i];
            cutsZTTmlSS[i] = "(fabs(EtaProbe)<1.2 && fabs(EtaProbe)>0.8)*"+cutsZTTmlSS[i];
            cutsZTTmt[i] = "(fabs(EtaProbe)<1.2 && fabs(EtaProbe)>0.8)*"+cutsZTTmt[i];
            cutsZTTmtSS[i] = "(fabs(EtaProbe)<1.2 && fabs(EtaProbe)>0.8)*"+cutsZTTmtSS[i];
        }
    }
    
    if(etaSuffix=="1p2to1p7")
    {
        for (int i=0; i<9; ++i)
        {
            cutsZMM[i] = "(fabs(EtaProbe)<1.7 && fabs(EtaProbe)>1.2)*"+cutsZMM[i];
            cutsZMMSS[i] = "(fabs(EtaProbe)<1.7 && fabs(EtaProbe)>1.2)*"+cutsZMMSS[i];
            cutsZJ[i] = "(fabs(EtaProbe)<1.7 && fabs(EtaProbe)>1.2)*"+cutsZJ[i];
            cutsZJSS[i] = "(fabs(EtaProbe)<1.7 && fabs(EtaProbe)>1.2)*"+cutsZJSS[i];
            cutsZTTml[i] = "(fabs(EtaProbe)<1.7 && fabs(EtaProbe)>1.2)*"+cutsZTTml[i];
            cutsZTTmlSS[i] = "(fabs(EtaProbe)<1.7 && fabs(EtaProbe)>1.2)*"+cutsZTTmlSS[i];
            cutsZTTmt[i] = "(fabs(EtaProbe)<1.7 && fabs(EtaProbe)>1.2)*"+cutsZTTmt[i];
            cutsZTTmtSS[i] = "(fabs(EtaProbe)<1.7 && fabs(EtaProbe)>1.2)*"+cutsZTTmtSS[i];
        }
    }
    if(etaSuffix=="Gt1p7")
    {
        for (int i=0; i<9; ++i)
        {
            cutsZMM[i] = "(fabs(EtaProbe)>1.7)*"+cutsZMM[i];
            cutsZMMSS[i] = "(fabs(EtaProbe)>1.7)*"+cutsZMMSS[i];
            cutsZJ[i] = "(fabs(EtaProbe)>1.7)*"+cutsZJ[i];
            cutsZJSS[i] = "(fabs(EtaProbe)>1.7)*"+cutsZJSS[i];
            cutsZTTml[i] = "(fabs(EtaProbe)>1.7)*"+cutsZTTml[i];
            cutsZTTmlSS[i] = "(fabs(EtaProbe)>1.7)*"+cutsZTTmlSS[i];
            cutsZTTmt[i] = "(fabs(EtaProbe)>1.7)*"+cutsZTTmt[i];
            cutsZTTmtSS[i] = "(fabs(EtaProbe)>1.7)*"+cutsZTTmtSS[i];
        }
    }
    // filling histograms for DYJets samples
    for (int i=0; i<9; ++i) { // run over samples
        
        TFile * file = new TFile("./NTuples/"+dySampleNames[i]+".root");
        TTree * tree = (TTree*)file->Get("MuTauFR");
        double normaliza = dyNorm[i];
        TString histNameZMM = dySampleNames[i] + "_"+varName+"_ZMM";
        TString histNameZMMSS = dySampleNames[i] + "_"+varName+"_ZMMSS";
        TString histNameZJ = dySampleNames[i] + "_"+varName+"_ZJ";
        TString histNameZJSS = dySampleNames[i] + "_"+varName+"_ZJSS";
        TString histNameZTTml = dySampleNames[i] + "_"+varName+"_ZTTml";
        TString histNameZTTmlSS = dySampleNames[i] + "_"+varName+"_ZTTmlSS";
        TString histNameZTTmt = dySampleNames[i] + "_"+varName+"_ZTTmt";
        TString histNameZTTmtSS = dySampleNames[i] + "_"+varName+"_ZTTmtSS";
        
        histZMM[i] = new TH1D(histNameZMM,"",nBins,xMin,xMax);
        histZMM[i]->Sumw2();
        histZMMSS[i] = new TH1D(histNameZMMSS,"",nBins,xMin,xMax);
        histZMMSS[i]->Sumw2();
        histZJ[i] = new TH1D(histNameZJ,"",nBins,xMin,xMax);
        histZJ[i]->Sumw2();
        histZJSS[i] = new TH1D(histNameZJSS,"",nBins,xMin,xMax);
        histZJSS[i]->Sumw2();
        histZTTml[i] = new TH1D(histNameZTTml,"",nBins,xMin,xMax);
        histZTTml[i]->Sumw2();
        histZTTmlSS[i] = new TH1D(histNameZTTmlSS,"",nBins,xMin,xMax);
        histZTTmlSS[i]->Sumw2();
        histZTTmt[i] = new TH1D(histNameZTTmt,"",nBins,xMin,xMax);
        histZTTmt[i]->Sumw2();
        histZTTmtSS[i] = new TH1D(histNameZTTmtSS,"",nBins,xMin,xMax);
        histZTTmtSS[i]->Sumw2();
        
        tree->Draw(varName+">>"+histNameZMM,cutsZMM[i] + extraCuts);
        tree->Draw(varName+">>"+histNameZMMSS,cutsZMMSS[i] + extraCuts);
        tree->Draw(varName+">>"+histNameZJ,cutsZJ[i] + extraCuts);
        tree->Draw(varName+">>"+histNameZJSS,cutsZJSS[i] + extraCuts);
        tree->Draw(varName+">>"+histNameZTTml,cutsZTTml[i] + extraCuts);
        tree->Draw(varName+">>"+histNameZTTmlSS,cutsZTTmlSS[i] + extraCuts);
        tree->Draw(varName+">>"+histNameZTTmt,cutsZTTmt[i] + extraCuts);
        tree->Draw(varName+">>"+histNameZTTmtSS,cutsZTTmtSS[i] + extraCuts);
        
        for (int iB=1; iB<=nBins; ++iB) {
            
            double x = histZMM[i]->GetBinContent(iB);
            double e = histZMM[i]->GetBinError(iB);
            histZMM[i]->SetBinContent(iB,normaliza*x);
            histZMM[i]->SetBinError(iB,normaliza*e);
            x = histZMMSS[i]->GetBinContent(iB);
            e = histZMMSS[i]->GetBinError(iB);
            histZMMSS[i]->SetBinContent(iB,normaliza*x);
            histZMMSS[i]->SetBinError(iB,normaliza*e);
            
            x = histZJ[i]->GetBinContent(iB);
            e = histZJ[i]->GetBinError(iB);
            histZJ[i]->SetBinContent(iB,normaliza*x);
            histZJ[i]->SetBinError(iB,normaliza*e);
            x = histZJSS[i]->GetBinContent(iB);
            e = histZJSS[i]->GetBinError(iB);
            histZJSS[i]->SetBinContent(iB,normaliza*x);
            histZJSS[i]->SetBinError(iB,normaliza*e);
            
            x = histZTTml[i]->GetBinContent(iB);
            e = histZTTml[i]->GetBinError(iB);
            histZTTml[i]->SetBinContent(iB,normaliza*x);
            histZTTml[i]->SetBinError(iB,normaliza*e);
            x = histZTTmlSS[i]->GetBinContent(iB);
            e = histZTTmlSS[i]->GetBinError(iB);
            histZTTmlSS[i]->SetBinContent(iB,normaliza*x);
            histZTTmlSS[i]->SetBinError(iB,normaliza*e);
            
            x = histZTTmt[i]->GetBinContent(iB);
            e = histZTTmt[i]->GetBinError(iB);
            histZTTmt[i]->SetBinContent(iB,normaliza*x);
            histZTTmt[i]->SetBinError(iB,normaliza*e);
            x = histZTTmtSS[i]->GetBinContent(iB);
            e = histZTTmtSS[i]->GetBinError(iB);
            histZTTmtSS[i]->SetBinContent(iB,normaliza*x);
            histZTTmtSS[i]->SetBinError(iB,normaliza*e);
        }
        std::cout << dySampleNames[i] << " -> DY ZMM = " << histZMM[i]->GetEntries() << " : " << histZMM[i]->GetSumOfWeights() << std::endl;
        std::cout << dySampleNames[i] << " -> DY ZJ = " << histZJ[i]->GetEntries() << " : " << histZJ[i]->GetSumOfWeights() << std::endl;
        std::cout << dySampleNames[i] << " -> DY ZTTml = " << histZTTml[i]->GetEntries() << " : " << histZTTml[i]->GetSumOfWeights() << std::endl;
        std::cout << dySampleNames[i] << " -> DY ZTTmt = " << histZTTmt[i]->GetEntries() << " : " << histZTTmt[i]->GetSumOfWeights() << std::endl;
    }
    
    hist[11]   = histZMM[0];
    histSS[11] = histZMMSS[0];
    hist[12]   = histZJ[0];
    histSS[12] = histZJSS[0];
    hist[13]   = histZTTml[0];
    histSS[13] = histZTTmlSS[0];
    hist[14]   = histZTTmt[0];
    histSS[14] = histZTTmtSS[0];
    for (int iDY=1; iDY<9; ++iDY)
    {
        hist[11]->Add(hist[11],histZMM[iDY]);
        histSS[11]->Add(histSS[11],histZMMSS[iDY]);
        hist[12]->Add(hist[12],histZJ[iDY]);
        histSS[12]->Add(histSS[12],histZJSS[iDY]);
        hist[13]->Add(hist[13],histZTTml[iDY]);
        histSS[13]->Add(histSS[13],histZTTmlSS[iDY]);
        hist[14]->Add(hist[14],histZTTmt[iDY]);
        histSS[14]->Add(histSS[14],histZTTmtSS[iDY]);
    }
    
    std::cout << "end: DY stitching" << std::endl;
    
    // *******************************
    // ***** W+Jets samples *******
    // *******************************
    TH1D * histW[9];
    TH1D * histWSS[9];
    TString wRefSamples[5] = {"NTuples_WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8",
        "NTuples_W1JetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8",
        "NTuples_W2JetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8",
        "NTuples_W3JetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8",
        "NTuples_W4JetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8"
    };
    
    double wRefXSec[5] = {61527,
        1.162*8104.0,
        1.162*2793.0,
        1.162*992.5,
        1.162*544.3}; //2017 values
    
    double wRefEvents[5] = {74635450,54988117,32368249,19700377,11333705};
    
    for (int iDY=0; iDY<5; ++iDY) {
        TFile * file = new TFile("./NTuples/"+wRefSamples[iDY]+".root");
        TH1D * histWeightsH = (TH1D*)file->Get("histWeightsH");
        wRefEvents[iDY] = histWeightsH->GetSumOfWeights();
    }
    TString wSampleNames[9] = {"NTuples_WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8",
			     "NTuples_WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8",
			     "NTuples_WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8",
			     "NTuples_WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8",
			     "NTuples_WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8",
        "NTuples_W1JetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8",
        "NTuples_W2JetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8",
        "NTuples_W3JetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8",
        "NTuples_W4JetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8"
    };
    
    double wNorm[9];
    wNorm[0] = lumi*wRefXSec[0]/wRefEvents[0];
    wNorm[1] = lumi/(wRefEvents[0]/wRefXSec[0]+wRefEvents[1]/wRefXSec[1]);
    wNorm[2] = lumi/(wRefEvents[0]/wRefXSec[0]+wRefEvents[2]/wRefXSec[2]);
    wNorm[3] = lumi/(wRefEvents[0]/wRefXSec[0]+wRefEvents[3]/wRefXSec[3]);
    wNorm[4] = lumi/(wRefEvents[0]/wRefXSec[0]+wRefEvents[4]/wRefXSec[4]);
    wNorm[5] = lumi/(wRefEvents[0]/wRefXSec[0]+wRefEvents[1]/wRefXSec[1]);
    wNorm[6] = lumi/(wRefEvents[0]/wRefXSec[0]+wRefEvents[2]/wRefXSec[2]);
    wNorm[7] = lumi/(wRefEvents[0]/wRefXSec[0]+wRefEvents[3]/wRefXSec[3]);
    wNorm[8] = lumi/(wRefEvents[0]/wRefXSec[0]+wRefEvents[4]/wRefXSec[4]);
    
    TString cutsW[9];
    TString cutsWSS[9];
    for (int iDY=0; iDY<9; ++iDY)
    {
        cutsW[iDY]   = WJetsweight+"puweight*effweight*mcweight*(os>0.5&&mt_1<40&&iso_1<0.15"+npartonCuts[iDY]+")";
        cutsWSS[iDY] = WJetsweight+qcdweight+"puweight*effweight*mcweight*(os<0.5&&mt_1<40&&iso_1<0.15"+npartonCuts[iDY]+")";
    }
    
    if(passProbe)
    {
        for (int i=0; i<9; ++i)
        {
            cutsW[i] ="("+wp+">0.5&&"+whatTauid+">0.5)*"+cutsW[i];
            cutsWSS[i] ="("+wp+">0.5&&"+whatTauid+">0.5)*"+cutsWSS[i];
        }
        
    }
    else
    {
        for (int i=0; i<9; ++i)
        {
            cutsW[i] ="("+wp+"<0.5 &&"+whatTauid+">0.5)*"+cutsW[i];
            cutsWSS[i] ="("+wp+"<0.5 && "+whatTauid+">0.5)*"+cutsWSS[i];
        }
    }
    
    
    if(etaSuffix=="Lt0p4")
    {
        for (int i=0; i<9; ++i)
        {
            cutsW[i] = "(fabs(EtaProbe)<0.4)*"+cutsW[i];
            cutsWSS[i] = "(fabs(EtaProbe)<0.4)*"+cutsWSS[i];
        }
    }
    
    if(etaSuffix=="0p4to0p8")
    {
        for (int i=0; i<9; ++i)
        {
            cutsW[i] = "(fabs(EtaProbe)<0.8 && fabs(EtaProbe)>0.4)*"+cutsW[i];
            cutsWSS[i] = "(fabs(EtaProbe)<0.8 && fabs(EtaProbe)>0.4)*"+cutsWSS[i];
        }
    }
    
    if(etaSuffix=="0p8to1p2")
    {
        for (int i=0; i<9; ++i)
        {
            cutsW[i] = "(fabs(EtaProbe)<1.2 && fabs(EtaProbe)>0.8)*"+cutsW[i];
            cutsWSS[i] = "(fabs(EtaProbe)<1.2 && fabs(EtaProbe)>0.8)*"+cutsWSS[i];
        }
    }

    if(etaSuffix=="1p2to1p7")
    {
        for (int i=0; i<9; ++i)
        {
            cutsW[i] ="(fabs(EtaProbe)<1.7 && fabs(EtaProbe)>1.2)*"+cutsW[i];
            cutsWSS[i] ="(fabs(EtaProbe)<1.7 && fabs(EtaProbe)>1.2)*"+cutsWSS[i];
        }
    }
    if(etaSuffix=="Gt1p7")
    {
        for (int i=0; i<9; ++i)
        {
            cutsW[i] ="(fabs(EtaProbe)>1.7)*"+cutsW[i];
            cutsWSS[i] ="(fabs(EtaProbe)>1.7)*"+cutsWSS[i];
        }
    }
    
    // filling histograms for WJets samples
    for (int i=0; i<9; ++i) { // run over samples
        
        TFile * file = new TFile("./NTuples/"+wSampleNames[i]+".root");
        TTree * tree = (TTree*)file->Get("MuTauFR");
        double normaliza = wNorm[i];
        TString histName = wSampleNames[i] + "_"+varName;
        TString histNameSS = wSampleNames[i] + "_"+varName+"_ss";
        histW[i] = new TH1D(histName,"",nBins,xMin,xMax);
        histW[i]->Sumw2();
        histWSS[i] = new TH1D(histNameSS,"",nBins,xMin,xMax);
        histWSS[i]->Sumw2();
        tree->Draw(varName+">>"+histName,cutsW[i]+extraCuts);
        tree->Draw(varName+">>"+histNameSS,cutsWSS[i]+extraCuts);
        
        for (int iB=1; iB<=nBins; ++iB) {
            
            double x = histW[i]->GetBinContent(iB);
            double e = histW[i]->GetBinError(iB);
            histW[i]->SetBinContent(iB,normaliza*x);
            histW[i]->SetBinError(iB,normaliza*e);
            x = histWSS[i]->GetBinContent(iB);
            e = histWSS[i]->GetBinError(iB);
            histWSS[i]->SetBinContent(iB,normaliza*x);
            histWSS[i]->SetBinError(iB,normaliza*e);
        }
        std::cout << wSampleNames[i] << " -> W = " << histW[i]->GetEntries() << " : " << histW[i]->GetSumOfWeights() << std::endl;
        //    delete file;
    }
    hist[15]   = histW[0];
    histSS[15] = histWSS[0];
    
    for (int iDY=1; iDY<9; ++iDY)
    {
        hist[15]->Add(hist[15],histW[iDY]);
        histSS[15]->Add(histSS[15],histWSS[iDY]);
    }

    // adding up ttbar backgrounds
    for (int iH=2; iH<4; ++iH)
    {
        hist[1]->Add(hist[1],hist[iH]);
    }
    //  adding up VV backgrounds
    for (int iH=5; iH<nSamples; ++iH)
    {
        hist[4]->Add(hist[4],hist[iH]);
    }
    
    // QCD
    for (int iH=1; iH<16; ++iH)
    {
        histSS[0]->Add(histSS[0],histSS[iH],1,-1);
    }
    
    //Getting rid of negative bin in QCD
    for (int iB=1; iB<=nBins; ++iB)
    {
        float ySS = histSS[0]->GetBinContent(iB);
        if (ySS<0)
        {
            histSS[0]->SetBinContent(iB,0.);
            histSS[0]->SetBinError(iB,0.);
        }
        
    }
    
    TFile * file = new TFile("./shapes/MuTauFR"+wp+etaSuffix+suffix+"_"+varName+".root","recreate");
    file->mkdir("MuTauFR"+suffix);
    file->cd("MuTauFR"+suffix);
    TH1D * data_obs = (TH1D*)hist[0]->Clone("data_obs");
    TH1D * QCD = (TH1D*)histSS[0]->Clone("QCD");
    TH1D * ZMM = (TH1D*)hist[11]->Clone("ZMM");
    TH1D * ZJ = (TH1D*)hist[12]->Clone("ZJ");
    TH1D * ZTT_ml = (TH1D*)hist[13]->Clone("ZTT_ml");
    TH1D * ZTT_mt = (TH1D*)hist[14]->Clone("ZTT_mt");
    TH1D * W = (TH1D*)hist[15]->Clone("W");
    TH1D * TT = (TH1D*)hist[1]->Clone("TT");
    TH1D * VV = (TH1D*)hist[4]->Clone("VV");
    
    float totData = data_obs->GetSumOfWeights();
    float totZMM = ZMM->GetSumOfWeights();
    float totZJ = ZJ->GetSumOfWeights();
    float totZTT_ml = ZTT_ml->GetSumOfWeights();
    float totZTT_mt = ZTT_mt->GetSumOfWeights();
    float totQCD = QCD->GetSumOfWeights();
    float totTT = TT->GetSumOfWeights();
    float totVV = VV->GetSumOfWeights();
    float totW = W->GetSumOfWeights();
    float totMC = totZMM+totZJ+totZTT_ml+totZTT_mt+totTT+totW+totVV+totQCD;
    

    std::cout << "data : " << totData << std::endl;
    std::cout << "QCD : " << QCD->GetSumOfWeights() << std::endl;
    std::cout << "VV  : " << VV->GetSumOfWeights() << std::endl;
    std::cout << "W   : " << W->GetSumOfWeights() << std::endl;
    std::cout << "TT  : " << TT->GetSumOfWeights() << std::endl;
    std::cout << "ZMM : " << ZMM->GetSumOfWeights() << std::endl;
    std::cout << "ZJ : " << ZJ->GetSumOfWeights() << std::endl;
    std::cout << "ZTT_ml : " << ZTT_ml->GetSumOfWeights() << std::endl;
    std::cout << "ZTT_mt : " << ZTT_mt->GetSumOfWeights() << std::endl;
    std::cout << "MC : " << totMC << std::endl;
    
    file->Write();
    file->Close();


}

