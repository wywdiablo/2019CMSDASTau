//----------------------Version 1.0-------------------------//
//Producing shapes input all for Mu -> Tau Fake Rates study
//Author: Yiwen Wen
//DESY
//-----------------------------------------------------------//
#include "HttStylesNew.cc"
#include "TColor.h"
#include "ProduceDatacardInputs_MuTauFR.C"
#include "ProduceDatacardInputs_MuTauFR_m_vis_Scale_ZTT.C"
#include "ProduceDatacardInputs_MuTauFR_m_vis_Scale_ZMM.C"
#include <math.h>

ProduceEveryShapesInput(TString eta="XXX",TString wp = "XXX", bool deepTau = true)//<--- please fill the eta region and anti-mu discriminiator you would like to measure
{
    //nominal
    ProduceDatacardInputs_MuTauFR("m_vis",eta,wp,70,120,10,true,deepTau,"*tauidweight");
    ProduceDatacardInputs_MuTauFR("m_vis",eta,wp,70,120,1,false,deepTau,"*tauidweight");
    
    //Probe Muon energy scale
    // ProduceDatacardInputs_MuTauFR_m_vis_Scale_ZMM("NAME OF THE SYSTEMATIC IN TREE","NAME OF THE SAME TEMPLATE TO BE STORED",eta,wp,LOWER BOUND,UPPER BOUND,BIN NUMBER,PASS OR FAIL,DEEPTAU OR NOT,"EXTRA CUTS");
    
    //<-- Please fill in the following Muon Energy scale shape systematic
    ProduceDatacardInputs_MuTauFR_m_vis_Scale_ZMM("XXX","probemuon_Up",eta,wp,70,120,10,true,deepTau,"*tauidweight");
    ProduceDatacardInputs_MuTauFR_m_vis_Scale_ZMM("XXX","probemuon_Up",eta,wp,70,120,1,false,deepTau,"*tauidweight");
    ProduceDatacardInputs_MuTauFR_m_vis_Scale_ZMM("XXX","probemuon_Down",eta,wp,70,120,10,true,deepTau,"*tauidweight");
    ProduceDatacardInputs_MuTauFR_m_vis_Scale_ZMM("XXX","probemuon_Down",eta,wp,70,120,1,false,deepTau,"*tauidweight");
    
    //Probe Tau energy scale
    //<-- Please fill in the following Tau Energy scale shape systematic

    ProduceDatacardInputs_MuTauFR_m_vis_Scale_ZTT("XXX","probetau_Up",eta,wp,70,120,10,true,deepTau,"*tauidweight");
    ProduceDatacardInputs_MuTauFR_m_vis_Scale_ZTT("XXX","probetau_Up",eta,wp,70,120,1,false,deepTau,"*tauidweight");
    ProduceDatacardInputs_MuTauFR_m_vis_Scale_ZTT("XXX","probetau_Down",eta,wp,70,120,10,true,deepTau,"*tauidweight");
    ProduceDatacardInputs_MuTauFR_m_vis_Scale_ZTT("XXX","probetau_Down",eta,wp,70,120,1,false,deepTau,"*tauidweight");

    
    //visbale mass resolution scale
    //<-- Please fill in the following visible mass resolution shape systematic
    ProduceDatacardInputs_MuTauFR_m_vis_Scale_ZMM("XXX","reso_Up",eta,wp,70,120,10,true,deepTau,"*tauidweight");
    ProduceDatacardInputs_MuTauFR_m_vis_Scale_ZMM("XXX","reso_Up",eta,wp,70,120,1,false,deepTau,"*tauidweight");
    ProduceDatacardInputs_MuTauFR_m_vis_Scale_ZMM("XXX","reso_Down",eta,wp,70,120,10,true,deepTau,"*tauidweight");
    ProduceDatacardInputs_MuTauFR_m_vis_Scale_ZMM("XXX","reso_Down",eta,wp,70,120,1,false,deepTau,"*tauidweight");
    
}

