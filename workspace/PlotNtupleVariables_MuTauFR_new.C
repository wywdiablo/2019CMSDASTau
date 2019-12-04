//----------------------Version of CMSDAS2019Beijing Tau exercise-------------------------//
//Ploting Ntuple variables for mu -> tau_h fake rate studies
//Author: Yiwen Wen
//DESY
//----------------------------------------------------------//
#include "HttStylesNew.cc"
#include "HtoH.h"
#include "TColor.h"
#include "CMS_lumi.C"
void PlotNtupleVariables_MuTauFR_new(
                                 TString varName= "XXX", //--> Input the name of the variable to you need to plot: "m_vis", "mt_1"...
                                 TString xtitle = "XXX", //--> Insert the title on the x-axis: "visible mass [GeV]"
                                 TString ytitle = "Events",
                                 float xLower = 0, // --> Lower bound of the variable you want to plot
                                 float xUpper = 200, // --> Upper bound of the variable
                                 int numberofbins = 40, //--> Bin numbers
                                 bool applyPU = true,
                                 bool passProbe = XXX, //--> Switch to look at "PASS" or "FAIL" region, true for "PASS", false for "FAIL"
                                 bool logY = false,
                                 bool legLeft = false,
                                 bool deepTau = true, //--> Flag to indicate the useage of DNN-based tau discriminator
                                 TString extraCuts = "") //--> You can specify any uniformly extra cuts to all processes
{
	SetStyle();
    
    // Every process has a root file. Here, we speficy all the input root files. For DY and W+jets, we will have the stitiched procedure, save them to later.
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

    //Specifying the luminosity
	float lumi = XXX; //insert the luminosity of the year 2018
	
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
    
    
    //Old tau id discrinimators
    TString whatTauid = "DecayModeFindingProbe>0.5 && taubyTightIsolationMVArun2017v2DBoldDMwLT2017";
    TString whatAntiMu = "tauagainstMuonLoose";
    
    //Initiation of deep tau discriminators
    if (deepTau)
    {
        whatTauid = "XXX"; // insert the cut of tau reconstuction and the VSjet discriminator, formatted as the one above
        whatAntiMu = "XXX"; // insert the VSmu discriminator
    }
    

	//inistiating cuts
	TString cuts[15];
	TString cutsSS[15];

    TString WJetsweight="0.950578*"; // W+jets normalization factor
    TString qcdweight ="1.06*"; // QCD (used same signed shape from Data), normalization factor 1.06 is the OS/SS ratio dervied from control region (relaxed isolation)
    
    
    // "cutsSS" and "histSS" are saved for QCD estimation, the differences are "os < 0.5" and requiring "qcdweight" to be applied
    // the QCD template at the end is the [SS data - (all MC SS)]
	for (int i=0; i<nSamples; ++i)
	{
        cuts[i] = "effweight*mcweight*(mt_1<40&&os>0.5 && iso_1<0.15)";
        cutsSS[i] = qcdweight+"effweight*mcweight*(mt_1<40 && os<0.5 &&iso_1<0.15)";
  	}

    //inititiation of cut on data. No mc weight, iso/id & trigger weight ("effweight") for data, because data doesn't need data-to-MC correction
  	cuts[0] = "(mt_1<40 && os>0.5 &&iso_1<0.15)";
    cutsSS[0] = qcdweight+"(mt_1<40 && os<0.5 && iso_1<0.15)";
    
    
    //defined PASS and FAIL region
    if(passProbe)
    {
        for (int i=0; i<nSamples; ++i)
        {
            cuts[i] ="(XXX&&"+whatTauid+">0.5)*"+cuts[i]; //insert here the definition of PASS and FAIL probe
            cutsSS[i] ="(XXX&&"+whatTauid+">0.5)*"+cutsSS[i]; //insert here the definition of PASS and FAIL probe
        }
        
    }
    else
    {
        for (int i=0; i<nSamples; ++i)
        {
            cuts[i] ="(XXX&&"+whatTauid+">0.5)*"+cuts[i]; //insert here the definition of PASS and FAIL probe
            cutsSS[i] ="(XXX&&"+whatTauid+">0.5)*"+cutsSS[i]; //insert here the definition of PASS and FAIL probe
        }
    }

    if(applyPU)
	{
		for (int i=1; i<nSamples; ++i)
		{
			cuts[i] ="puweight*"+cuts[i];
            cutsSS[i] ="puweight*"+cutsSS[i];
		}
	}
	for (int i=0; i<nSamples; ++i)
	{
		std::cout <<i<<":"+samples[i]<<std::endl;
		TFile * file = new TFile("XXX"+samples[i]+".root"); //reading the root files <-- chane the directory to ready the files
		TH1D * histWeightsH = (TH1D*)file->Get("histWeightsH");
		TTree * tree = (TTree*)file->Get("MuTauFR"); //reading the tree in the root file
        double normaliza = xsec[i]*lumi/histWeightsH->GetSumOfWeights(); //here normaliza the MC with the luminosity, "histWeightsH->GetSumOfWeights()" is the total number of weighted events, for LO with mcweight = 1, one can sure use get entires. But for the case of mc weight != 1, one should use GetSumOfWeights
		TString histName = samples[i] + "_"+varName;
        TString histNameSS = samples[i] + "_"+varName+"_ss";
        hist[i] = new TH1D(histName,"",nBins,xMin,xMax);
		hist[i]->Sumw2();
        histSS[i] = new TH1D(histNameSS,"",nBins,xMin,xMax);
        histSS[i]->Sumw2();
        tree->Draw(varName+">>"+histName,cuts[i] + extraCuts);//Draw variable to the histrogram with cuts applied!!
        tree->Draw(varName+">>"+histNameSS,cutsSS[i] + extraCuts);
		if(i > 0) //here normaliza the MC with the luminosity
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
    
    // ****************************************
    // ***** Stitching of DY+Jets samples *****
    // ****************************************
    
    // stitching means to combine the the inclusive sample with the exclusive (binned in number of jets) sample
    // The purpose is to have more event statistic than using the inclusive sample only
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
        1.137*44.03};//2017 or 2018 values
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
    
    //the stitching is based in the generator information of number of parton in an event.
    //Split the inclusive DY events in five parts based on the GEN number of parton
    //DY jet-binned sample already have the corresponding GEN number of parton applied.
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
        // inintiate the gen matching cut for muon and hadronic tau
        
        // what are the gen matching codes for Z -> mu mu?
        cutsZMM[iDY]   = "zptweight*effweight*mcweight*(mt_1<40&&os>0.5&& XXX && XXX  &&iso_1<0.15"+npartonCuts[iDY]+")";
        cutsZMMSS[iDY] = qcdweight+"zptweight*effweight*mcweight*(mt_1<40&&os<0.5&& XXX && XXX &&iso_1<0.15"+npartonCuts[iDY]+")";
        // The definition of ZJ is to have the hadronic tau matched to a jet
        cutsZJ[iDY]   = "zptweight*effweight*mcweight*(mt_1<40&&os>0.5&& gen_match_2 == 6 &&iso_1<0.15"+npartonCuts[iDY]+")";
        cutsZJSS[iDY] = qcdweight+"zptweight*effweight*mcweight*(mt_1<40&&os<0.5&& gen_match_2 == 6 &&iso_1<0.15"+npartonCuts[iDY]+")";
        // what are the gen matching codes for Z -> TT -> mu mu/ele?
        cutsZTTml[iDY]   = "zptweight*effweight*mcweight*(mt_1<40&&os>0.5&& XXX &&iso_1<0.15"+npartonCuts[iDY]+")";
        cutsZTTmlSS[iDY] = qcdweight+"zptweight*effweight*mcweight*(mt_1<40&&os<0.5&& XXX &&iso_1<0.15"+npartonCuts[iDY]+")";
        // what are the gen matching codes for Z -> TT -> mu tauh?
        cutsZTTmt[iDY]   = "zptweight*effweight*mcweight*(mt_1<40&&os>0.5&& XXX &&iso_1<0.15"+npartonCuts[iDY]+")";
        cutsZTTmtSS[iDY] = qcdweight+"zptweight*effweight*mcweight*(mt_1<40&&os<0.5&& XXX &&iso_1<0.15"+npartonCuts[iDY]+")";
        
    }
    if(passProbe)
    {
        for (int i=0; i<9; ++i)
        {
            cutsZMM[i] ="("+whatAntiMu+">0.5&&"+whatTauid+">0.5)*"+cutsZMM[i];
            cutsZMMSS[i] ="("+whatAntiMu+">0.5&&"+whatTauid+">0.5)*"+cutsZMMSS[i];
            cutsZJ[i] ="("+whatAntiMu+">0.5&&"+whatTauid+">0.5)*"+cutsZJ[i];
            cutsZJSS[i] ="("+whatAntiMu+">0.5&&"+whatTauid+">0.5)*"+cutsZJSS[i];
            cutsZTTml[i] ="("+whatAntiMu+">0.5&&"+whatTauid+">0.5)*"+cutsZTTml[i];
            cutsZTTmlSS[i] ="("+whatAntiMu+">0.5&&"+whatTauid+">0.5)*"+cutsZTTmlSS[i];
            cutsZTTmt[i] ="("+whatAntiMu+">0.5&&"+whatTauid+">0.5)*"+cutsZTTmt[i];
            cutsZTTmtSS[i] ="("+whatAntiMu+">0.5&&"+whatTauid+">0.5)*"+cutsZTTmtSS[i];
        }
        
    }
    else
    {
        for (int i=0; i<9; ++i)
        {
            cutsZMM[i] ="("+whatAntiMu+"<0.5&&"+whatTauid+">0.5)*"+cutsZMM[i];
            cutsZMMSS[i] ="("+whatAntiMu+"<0.5&&"+whatTauid+">0.5)*"+cutsZMMSS[i];
            cutsZJ[i] ="("+whatAntiMu+"<0.5&&"+whatTauid+">0.5)*"+cutsZJ[i];
            cutsZJSS[i] ="("+whatAntiMu+"<0.5&&"+whatTauid+">0.5)*"+cutsZJSS[i];
            cutsZTTml[i] ="("+whatAntiMu+"<0.5&&"+whatTauid+">0.5)*"+cutsZTTml[i];
            cutsZTTmlSS[i] ="("+whatAntiMu+"<0.5&&"+whatTauid+">0.5)*"+cutsZTTmlSS[i];
            cutsZTTmt[i] ="("+whatAntiMu+"<0.5&&"+whatTauid+">0.5)*"+cutsZTTmt[i];
            cutsZTTmtSS[i] ="("+whatAntiMu+"<0.5&&"+whatTauid+">0.5)*"+cutsZTTmtSS[i];
            
        }
    }
    if(applyPU)
    {
        for (int i=0; i<9; ++i)
        {
            cutsZMM[i] ="puweight*"+cutsZMM[i];
            cutsZMMSS[i] ="puweight*"+cutsZMMSS[i];
            cutsZJ[i] ="puweight*"+cutsZJ[i];
            cutsZJSS[i] ="puweight*"+cutsZJSS[i];
            cutsZTTml[i] ="puweight*"+cutsZTTml[i];
            cutsZTTmlSS[i] ="puweight*"+cutsZTTmlSS[i];
            cutsZTTmt[i] ="puweight*"+cutsZTTmt[i];
            cutsZTTmtSS[i] ="puweight*"+cutsZTTmtSS[i];
            
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
    // ****************************************
    // ***** Stitching of W+Jets samples ******
    // ****************************************
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
        cutsW[iDY]   = WJetsweight+"effweight*mcweight*(os>0.5&&mt_1<40&&iso_1<0.15"+npartonCuts[iDY]+")";
        cutsWSS[iDY] = WJetsweight+qcdweight+"effweight*mcweight*(os<0.5&&mt_1<40&&iso_1<0.15"+npartonCuts[iDY]+")";
    }
    
    if(passProbe)
    {
        for (int i=0; i<9; ++i)
        {
            cutsW[i] ="("+whatAntiMu+">0.5&&"+whatTauid+">0.5)*"+cutsW[i];
            cutsWSS[i] ="("+whatAntiMu+">0.5&&"+whatTauid+">0.5)*"+cutsWSS[i];
        }
        
    }
    else
    {
        for (int i=0; i<9; ++i)
        {
            cutsW[i] ="("+whatAntiMu+"<0.5 && "+whatTauid+">0.5)*"+cutsW[i];
            cutsWSS[i] ="("+whatAntiMu+"<0.5 && "+whatTauid+">0.5)*"+cutsWSS[i];
        }
    }
    
    if(applyPU)
    {
        for (int i=0; i<9; ++i)
        {
            cutsW[i] ="puweight*"+cutsW[i];
            cutsWSS[i] ="puweight*"+cutsWSS[i];
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
        tree->Draw(varName+">>"+histName,cutsW[i] + extraCuts);
        tree->Draw(varName+">>"+histNameSS,cutsWSS[i] + extraCuts);
        
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
    
    //QCD
    //the QCD template at the end is the [SS data - (all MC SS)]
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

    
    float sfW = (totData-totZMM-totZJ-totZTT_ml-totZTT_mt-totTT-totQCD-totVV)/totW;
    float sfWerr = TMath::Sqrt(totData)/totW;
    std::cout << "Data: " << data_obs->GetSumOfWeights() << std::endl;
    std::cout << "QCD : " << QCD->GetSumOfWeights() << std::endl;
    std::cout << "VV  : " << VV->GetSumOfWeights() << std::endl;
    std::cout << "W   : " << W->GetSumOfWeights() << std::endl;
    std::cout << "TT  : " << TT->GetSumOfWeights() << std::endl;
    std::cout << "ZMM : " << ZMM->GetSumOfWeights() << std::endl;
    std::cout << "ZJ : " << ZJ->GetSumOfWeights() << std::endl;
    std::cout << "ZTT_ml : " << ZTT_ml->GetSumOfWeights() << std::endl;
    std::cout << "ZTT_mt : " << ZTT_mt->GetSumOfWeights() << std::endl;
    std::cout << "sfW = " << sfW << " +/- " << sfWerr << std::endl;

    //incluce uncertainties in the plot
	TH1D * dummy = (TH1D*)ZMM->Clone("dummy");
    float errLumi = 0.03;
    float errQCD = 0.1;
    float errDY=0.1;
    float errVV = 0.15;
    float errW = 0.15;
    float errTT = 0.15;
    
    for (int iB=1; iB<=nBins; ++iB)
    {
        float eQCD = errQCD*QCD->GetBinContent(iB);
        float eVV = errVV*VV->GetBinContent(iB);
        float eDYMM = errDY*ZMM->GetBinContent(iB);
        float eDYTT_ml = errDY*ZTT_ml->GetBinContent(iB);
        float eDYTT_mt = errDY*ZTT_mt->GetBinContent(iB);
        float eDYJ = errDY*ZJ->GetBinContent(iB);
        float eW = errW*W->GetBinContent(iB);
        float eTT = errTT*TT->GetBinContent(iB);
        float err2 = eQCD*eQCD + eVV*eVV + eW*eW + eTT*eTT+eDYMM*eDYMM+eDYTT_ml*eDYTT_ml+eDYTT_mt*eDYTT_mt+eDYJ*eDYJ;
        float errTot = TMath::Sqrt(err2);
        dummy->SetBinError(iB,errTot);
    }
    
    //Stacking histograms
	W->Add(W,QCD);
	TT->Add(TT,W);
	VV->Add(VV,TT);
	ZTT_mt->Add(ZTT_mt,VV);
    ZTT_ml->Add(ZTT_ml,ZTT_mt);
    ZJ->Add(ZJ,ZTT_ml);
	ZMM->Add(ZMM,ZJ);
	
	TH1D * bkgdErr = (TH1D*)ZMM->Clone("bkgdErr");
  	bkgdErr->SetFillStyle(3013);
  	bkgdErr->SetFillColor(1);
  	bkgdErr->SetMarkerStyle(21);
  	bkgdErr->SetMarkerSize(0);

  	for (int iB=1; iB<=nBins; ++iB)
    {  
        W->SetBinError(iB,0);
        QCD->SetBinError(iB,0);
        VV->SetBinError(iB,0);
        TT->SetBinError(iB,0);
        ZTT_mt->SetBinError(iB,0);
        ZTT_ml->SetBinError(iB,0);
        ZJ->SetBinError(iB,0);
        ZMM->SetBinError(iB,0);
        float eStat =  bkgdErr->GetBinError(iB);
        float X = bkgdErr->GetBinContent(iB);
        float eLumi = errLumi * X;
        float eBkg = dummy->GetBinError(iB);
        float Err = TMath::Sqrt(eStat*eStat+eLumi*eLumi+eBkg*eBkg);
        bkgdErr->SetBinError(iB,Err);
    }  
   	//Colors 
 	Int_t colorZMM = TColor::GetColor("#ffcc66");
    Int_t colorZOther = TColor::GetColor("#4496C8");
	Int_t colorTT = TColor::GetColor("#9999CC");
	Int_t colorVV = TColor::GetColor("#6F2D35");
	Int_t colorW = TColor::GetColor("#DE5A6A");
	Int_t colorQCD = TColor::GetColor("#FFCCFF");

	InitData(data_obs);
	InitHist(QCD,"","",colorQCD,1001);
	InitHist(W,"","",colorW,1001);
	InitHist(TT,"","",colorTT,1001);
	InitHist(VV,"","",colorVV,1001);
	InitHist(ZJ,"","",colorZOther,1001);
	InitHist(ZMM,"","",colorZMM,1001);
	data_obs->GetXaxis()->SetTitle(xtitle);
 	data_obs->GetYaxis()->SetTitle(ytitle);
	data_obs->GetYaxis()->SetTitleOffset(1.5);
	data_obs->GetYaxis()->SetTitleSize(0.06);
	data_obs->GetXaxis()->SetRangeUser(xLower,xUpper);
	float yUpper = data_obs->GetMaximum();
	if (logY)
	data_obs->GetYaxis()->SetRangeUser(0.5,2*yUpper);
	else
	data_obs->GetYaxis()->SetRangeUser(0,1.2*yUpper);
	data_obs->SetMarkerSize(1.5);
	data_obs->GetXaxis()->SetLabelSize(0);
	data_obs->GetYaxis()->SetLabelSize(0.06);
	
	TCanvas * canv1 = MakeCanvas("canv1", "", 700, 800);
	
	TPad *upper = new TPad("upper", "pad",0,0.31,1,1);
	upper->Draw();
	upper->cd();
	upper->SetFillColor(0);
	upper->SetBorderMode(0);
	upper->SetBorderSize(10);
	upper->SetTickx(1);
	upper->SetTicky(1);
	upper->SetLeftMargin(0.17);
	upper->SetRightMargin(0.05);
	upper->SetBottomMargin(0.02);
 	upper->SetFrameFillStyle(0);
	upper->SetFrameLineStyle(0);
  	upper->SetFrameLineWidth(2);
  	upper->SetFrameBorderMode(0);
  	upper->SetFrameBorderSize(10);
  	upper->SetFrameFillStyle(0);
  	upper->SetFrameLineStyle(0);
  	upper->SetFrameLineWidth(2);
  	upper->SetFrameBorderMode(0);
  	upper->SetFrameBorderSize(10);

	//Drawing histogram
  	data_obs->Draw("e1");
	ZMM->Draw("sameh");
  	ZJ->Draw("sameh");//ZJ is ZOther
  	VV->Draw("sameh");
  	TT->Draw("sameh");
	W->Draw("sameh");
	QCD->Draw("sameh");
  	data_obs->Draw("e1same");
	bkgdErr->Draw("e2same");
	//Calculating chi2
	float chi2 = 0;
	for (int iB=1; iB<=nBins; ++iB) 
	{
		float xData = data_obs->GetBinContent(iB);
		float xMC = ZMM->GetBinContent(iB);
		if (xMC>1e-1) 
		{
      			float diff2 = (xData-xMC)*(xData-xMC);
      			chi2 += diff2/xMC;
    		}
  	}
  	std::cout << std::endl;
  	std::cout << "Chi2 = " << chi2 << std::endl;
  	std::cout << std::endl;

  	float x1Leg = 0.60;
  	float x2Leg = 0.85;
  	if (legLeft) 
	{
    		x1Leg = 0.20;
    		x2Leg = 0.45;
  	}
	TLegend * leg = new TLegend(x1Leg,0.6,x2Leg,0.88);
  	SetLegendStyle(leg);
  	leg->SetTextSize(0.05);
  	leg->AddEntry(data_obs,"Data","lp");
  	leg->AddEntry(VV,"Dibosons","f");
  	leg->AddEntry(W,"WJets","f");
 	leg->AddEntry(QCD,"QCD","f");
  	leg->AddEntry(TT,"t#bar{t}","f");
	leg->AddEntry(ZJ,"DY others","f");
  	leg->AddEntry(ZMM,"Z#rightarrow#mu#mu","f");
  	leg->Draw();
  	//plotchannel("e#mu");
    //	if (!applyPU) suffix = "_noPU";

  	TLatex * cms = new TLatex(0.20,0.94,"CMS L = 59.7 fb^{-1} at #sqrt{s} = 13 TeV");

  	cms->SetNDC();
  	cms->SetTextSize(0.05);
  	cms->Draw();
    TLatex * workinprogress = new TLatex(x1Leg,0.55,"Work in progress");
    workinprogress->SetNDC();
    workinprogress->SetTextSize(0.05);
    workinprogress->Draw();
    
    
  	if (logY) upper->SetLogy(true);

 	upper->Draw("SAME");
  	upper->RedrawAxis();
  	upper->Modified();
  	upper->Update();
 	canv1->cd();
	
	TH1D * ratioH = (TH1D*)data_obs->Clone("ratioH");
	TH1D * ratioErrH = (TH1D*)bkgdErr->Clone("ratioErrH");
  	ratioH->SetMarkerColor(1);
 	ratioH->SetMarkerStyle(20);
  	ratioH->SetMarkerSize(1.5);
  	ratioH->SetLineColor(1);
  	ratioH->GetYaxis()->SetRangeUser(0.3,1.8);
 	ratioH->GetYaxis()->SetNdivisions(505);
  	ratioH->GetXaxis()->SetLabelFont(42);
  	ratioH->GetXaxis()->SetLabelOffset(0.04);
  	ratioH->GetXaxis()->SetLabelSize(0.14);
  	ratioH->GetXaxis()->SetTitleSize(0.13);
  	ratioH->GetXaxis()->SetTitleOffset(1.2);
  	ratioH->GetYaxis()->SetTitle("obs/exp");
 	ratioH->GetYaxis()->SetLabelFont(42);
  	ratioH->GetYaxis()->SetLabelOffset(0.015);
  	ratioH->GetYaxis()->SetLabelSize(0.13);
  	ratioH->GetYaxis()->SetTitleSize(0.14);
  	ratioH->GetYaxis()->SetTitleOffset(0.5);
  	ratioH->GetXaxis()->SetTickLength(0.07);
  	ratioH->GetYaxis()->SetTickLength(0.04);
  	ratioH->GetYaxis()->SetLabelOffset(0.01);
	
	for (int iB=1; iB<=nBins; ++iB) 
	{
    		float x1 = data_obs->GetBinContent(iB);
    		float x2 = ZMM->GetBinContent(iB);
    		ratioErrH->SetBinContent(iB,1.0);
    		ratioErrH->SetBinError(iB,0.0);
		float xBkg = bkgdErr->GetBinContent(iB);
   		float errBkg = bkgdErr->GetBinError(iB);
		if (xBkg>0) 
		{
     			float relErr = errBkg/xBkg;
      			ratioErrH->SetBinError(iB,relErr);
    		}
		if (x1>0&&x2>0) 
		{
      			float e1 = data_obs->GetBinError(iB);
      			float ratio = x1/x2;
      			float eratio = e1/x2;
      			ratioH->SetBinContent(iB,ratio);
     			ratioH->SetBinError(iB,eratio);
    		}
   		else 
		{
      			ratioH->SetBinContent(iB,1000);
    		}
  	}


	TPad *lower = new TPad("lower", "pad",0,0,1,0.30);
  	lower->Draw();
  	lower->cd();
  	lower->SetFillColor(0);
  	lower->SetBorderMode(0);
  	lower->SetBorderSize(10);
 	lower->SetGridy();
  	lower->SetTickx(1);
  	lower->SetTicky(1);
  	lower->SetLeftMargin(0.17);
 	lower->SetRightMargin(0.05);
  	lower->SetTopMargin(0.026);
  	lower->SetBottomMargin(0.35);
  	lower->SetFrameFillStyle(0);
  	lower->SetFrameLineStyle(0);
  	lower->SetFrameLineWidth(2);
  	lower->SetFrameBorderMode(0);
  	lower->SetFrameBorderSize(10);
  	lower->SetFrameFillStyle(0);
  	lower->SetFrameLineStyle(0);
  	lower->SetFrameLineWidth(2);
  	lower->SetFrameBorderMode(0);
  	lower->SetFrameBorderSize(10);

  	ratioH->Draw("e1");
	ratioErrH->Draw("e2same");
  	lower->Modified();
  	lower->RedrawAxis();
  	canv1->cd();
  	canv1->Modified();
  	canv1->cd();
  	canv1->SetSelected(canv1);
    TString deepTauSuffix = "";
    if(deepTau)
        deepTauSuffix = "DeepTau";
	canv1->Print(varName+deepTauSuffix+".png");
	canv1->Print(varName+deepTauSuffix+".pdf","Portrait pdf");//--> output the pdf file


}

