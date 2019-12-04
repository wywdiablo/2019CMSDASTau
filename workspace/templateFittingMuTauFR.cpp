#include <string>
#include <map>
#include <set>
#include <iostream>
#include <utility>
#include <vector>
#include <cstdlib>
#include "CombineHarvester/CombineTools/interface/CombineHarvester.h"
#include "CombineHarvester/CombineTools/interface/Observation.h"
#include "CombineHarvester/CombineTools/interface/Process.h"
#include "CombineHarvester/CombineTools/interface/Utilities.h"
#include "CombineHarvester/CombineTools/interface/Systematics.h"
#include "CombineHarvester/CombineTools/interface/BinByBin.h"
#include "TRegexp.h"

using namespace std;

int main(int argc, char * argv[]) {
  //! [part1]
  // First define the location of the "auxiliaries" directory where we can
  // source the input files containing the datacard shapes
  string aux_shapes = "XXX"; //<-- HERE, change to the path to your shapes!!!

  // Create an empty CombineHarvester instance that will hold all of the
  // datacard configuration and histograms etc.
  ch::CombineHarvester cb;
  // Uncomment this next line to see a *lot* of debug information
  // cb.SetVerbosity(3);

  // Here we will just define two categories for an 8TeV analysis. Each entry in
  // the vector below specifies a bin name and corresponding bin_id.
  ch::Categories cats = {
      {1, "MuTauFR_pass"},
      {2, "MuTauFR_fail"}
  };
  // ch::Categories is just a typedef of vector<pair<int, string>>
  //! [part1]


  //! [part2]
  //vector<string> masses = ch::MassesFromRange("120-135:5");
  // Or equivalently, specify the mass points explicitly:
    vector<string> masses = {"90"};
  //! [part2]

  //! [part3]
  cb.AddObservations({"*"}, {"MuTauFR"}, {"13TeV"}, {"mt"}, cats);
  //! [part3]

  //! [part4]
  vector<string> bkg_procs = {"XXX","XXX","XXX", "XXX", "XXX", "XXX","XXX"};//<-- Change to your background processes, which are?
  cb.AddProcesses({"*"}, {"MuTauFR"}, {"13TeV"}, {"mt"}, bkg_procs, cats, false);

  vector<string> sig_procs = {"XXX"}; //<-- Change to your signal process, which is?
  cb.AddProcesses(masses, {"MuTauFR"}, {"13TeV"}, {"mt"}, sig_procs, cats, true);
  //! [part4]


  //Some of the code for this is in a nested namespace, so
  // we'll make some using declarations first to simplify things a bit.
  using ch::syst::SystMap;
  using ch::syst::era;
  using ch::syst::bin_id;
  using ch::syst::process;


  //! [part5]
    //most of the systematic has two types: normalization systematic affecting the yield "lnN" and shape systematic "shape"
    //<-- here is the luminosity uncertainty
    cb.cp().process(ch::JoinStr({sig_procs,{"ZJ"},{"ZTT_ml"},{"ZTT_mt"},{"VV"},{"TT"}})).AddSyst(cb, "lumi_$ERA", "lnN", SystMap<era>::init({"13TeV"}, XXX)); //<-- fill in the uncertainty of luminosity of 2018, 1.026
    
  //! [part6]
    cb.cp().process(ch::JoinStr({sig_procs})).AddSyst(cb, "probemuon_", "shape", SystMap<>::init(1));//<-- this is an example
    
    //cb.cp().process(PROCESS).AddSyst(cb, "NAME OF SHAPE IN FILE WITHOUT 'UP'OR'DOWN'", "TYPE OF SYSTEMATIC(NORMALIZATION OR SHAPE?)", SystMap<>::init(VALUE(SHAPE HAS 1 BUT NORMALIZATION YOU CAN SPECIFY YOUR VALUE)));
    
    //<-- here fill in the tau energy scale shape systematic for Z to tau tau to mu tau event
    cb.cp().process(XXX).AddSyst(cb, "XXX", "XXX", SystMap<>::init(1));
    
    //<-- here fill in the visible mass resolution systematic for Z to mu mu event
    cb.cp().process(XXX).AddSyst(cb, "XXX", "XXX", SystMap<>::init(1));
    
    //<-- here fill in the cross section uncertainty for W+jets
    cb.cp().process(XXX).AddSyst(cb, "XXX", "XXX", SystMap<>::init(XXX));
    
    cb.cp().process({"ZJ"}).AddSyst(cb, "jetTauFR", "lnN", SystMap<>::init(1.30));

    cb.cp().process(ch::JoinStr({{"ZMM"},{"ZJ"},{"ZTT_mt"},{"ZTT_ml"}})).AddSyst(cb, "normalizationDY", "lnN", SystMap<>::init(1.03));
    
    //(ignore me!This is a hidden level) this is the scale factor for VSjet on Z->MuMu
    //cb.cp().process(ch::JoinStr({sig_procs})).AddSyst(cb, "normalizationDYMuMu", "rateParam", SystMap<>::init(1));

    cb.cp().process({"VV"}).AddSyst(cb, "normalizationVV", "lnN", SystMap<>::init(1.15));
    cb.cp().process({"TT"}).AddSyst(cb, "normalizationTT", "lnN", SystMap<>::init(1.10));
    
    cb.cp().process(ch::JoinStr({sig_procs,{"ZJ"},{"ZTT_mt"},{"ZTT_ml"},{"VV"},{"TT"}})).AddSyst(cb, "CMS_eff_m", "lnN", SystMap<>::init(1.02));
    cb.cp().process(ch::JoinStr({sig_procs,{"ZTT_mt"},{"VV"},{"TT"}})).AddSyst(cb, "CMS_eff_t", "lnN", SystMap<>::init(1.02));
    cb.cp().process({"QCD"}).AddSyst(cb, "normalizationQCD", "lnN", SystMap<>::init(1.2));
    
    
  //! [part6]

  //! [part7]
  cb.cp().backgrounds().ExtractShapes(
      aux_shapes + argv[1],
      "$BIN/$PROCESS",
      "$BIN/$PROCESS_$SYSTEMATIC");
  cb.cp().signals().ExtractShapes(
      aux_shapes + argv[1],
      "$BIN/$PROCESS",
      "$BIN/$PROCESS_$SYSTEMATIC");
  //! [part7]

  //! [part8]
    auto bbb = ch::BinByBinFactory().SetAddThreshold(0.1).SetFixNorm(true);
    
    bbb.AddBinByBin(cb.cp().backgrounds(),cb);
    
    TString outputdcname = (TString) argv[1];
    TRegexp re(".root");
    outputdcname(re) = ".txt";
  //! [part8]

  //! [part9]
  // First we generate a set of bin names:
  set<string> bins = cb.bin_set();
  // This method will produce a set of unique bin names by considering all
  // Observation, Process and Systematic entries in the CombineHarvester
  // instance.

  // We create the output root file that will contain all the shapes.
  TFile output("htt_mt.input.root", "RECREATE");

  // Finally we iterate through each bin,mass combination and write a
  // datacard.
    cb.cp().WriteDatacard((string)outputdcname, output);
  //! [part9]
    cout << "pre-fit fake rate: "
    
    // please fill in the fake rate for pre fit, HINT: use 'cb.cp().bin({"XXX"}).process({"XXX"}).GetRate()' to get the yield of a process in a region(bin), fake rate is defined as #PASS/#TOTAL
    
    float sigRatePassPre = cb.cp().bin({"MuTauFR_pass"}).process({"ZMM"}).GetRate();
    float sigRateFailPre = cb.cp().bin({"MuTauFR_fail"}).process({"ZMM"}).GetRate();
    float sigErrPassPre = cb.cp().bin({"MuTauFR_pass"}).process({"ZMM"}).GetUncertainty();
    float sigErrFailPre = cb.cp().bin({"MuTauFR_fail"}).process({"ZMM"}).GetUncertainty();
    
    float dfdxPre = sigRateFailPre/((sigRatePassPre+sigRateFailPre)*(sigRatePassPre+sigRateFailPre));
    float dfdyPre = - sigRatePassPre/((sigRatePassPre+sigRateFailPre)*(sigRatePassPre+sigRateFailPre));
    float errfakeratePrefit= sqrt((dfdxPre*sigErrPassPre)*(dfdxPre*sigErrPassPre)+(dfdyPre*sigErrFailPre)*(dfdyPre*sigErrFailPre));
    
    cout << "pre-fit fake rate errors:" << errfakeratePrefit << endl;
    
    
}
