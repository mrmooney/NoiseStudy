#ifndef NoiseAnalyzer_h
#define NoiseAnalyzer_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TH2.h>

#include <vector>
#include <iostream>

using namespace std;

class NoiseAnalyzer {
public :
   TFile          *fFile;
   TH1F           *fHist;
   Int_t           fRunNumber;
   Int_t           fSubrunNumber;
   Int_t           fEventNumber;

   TH1F           *avgMeanHist;
   TH1F           *avgRMSHist;
   TH1F           *avgMaxHist;
   TH1F           *chirpFracHist;
   TH1F           *zigzagFracHist;
   TH1F           *waveFracHist;
   TH1F           *avgFiltRMSHist;
   TH1F           *lowCountFracHist;
   TH1F           *high1CountFracHist;
   TH1F           *high2CountFracHist;

   TH1F           *avgMeanElecHist;
   TH1F           *avgRMSElecHist;
   TH1F           *avgMaxElecHist;
   TH1F           *chirpFracElecHist;
   TH1F           *zigzagFracElecHist;
   TH1F           *waveFracElecHist;
   TH1F           *avgFiltRMSElecHist;
   TH1F           *lowCountFracElecHist;
   TH1F           *high1CountFracElecHist;
   TH1F           *high2CountFracElecHist;

   NoiseAnalyzer(Int_t runNumber = 1532, Int_t subrunNumber = 0, Int_t eventNumber = 0);
   virtual ~NoiseAnalyzer();
   virtual Int_t    LoadHist(Int_t planeNum, Int_t wireNum);
   virtual void     RunAnalyzer();

   virtual pair<Double_t,Double_t> CalcMeanAndRMS();
   virtual pair<Double_t,Double_t> GetExtrema();
   virtual Int_t      CountSigBins(Double_t meanVal);
   virtual Double_t   ChirpStudyAlg();
   virtual Double_t   ZigzagStudyAlg();
   virtual Double_t   WaveStudyAlg();
   virtual Double_t   FiltRMSStudyAlg();
   virtual vector<TH1F*>  GetHists();
};

#endif

#ifdef NoiseAnalyzer_cxx
NoiseAnalyzer::NoiseAnalyzer(Int_t runNumber, Int_t subrunNumber, Int_t eventNumber)
{
   fRunNumber = runNumber;
   fSubrunNumber = subrunNumber;
   fEventNumber = eventNumber;

   fFile = new TFile(Form("data/processed/celltree_Run%d_Subrun%d_Event%d.root",runNumber,subrunNumber,eventNumber));
   fHist = 0;
}

NoiseAnalyzer::~NoiseAnalyzer()
{
   if (!fFile) return;
   delete fFile;
}

Int_t NoiseAnalyzer::LoadHist(Int_t planeNumber, Int_t wireNumber)
{
  if(planeNumber == 0)
    fFile->GetObject(Form("U1_%d",wireNumber),fHist);
  else if(planeNumber == 1)
    fFile->GetObject(Form("V1_%d",wireNumber),fHist);
  else if(planeNumber == 2)
    fFile->GetObject(Form("W1_%d",wireNumber),fHist);

  if(!fHist)
    return -1;
  else
    return 0;
}

#endif // #ifdef NoiseAnalyzer_cxx
