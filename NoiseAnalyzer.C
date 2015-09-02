#define NoiseAnalyzer_cxx
#include "NoiseAnalyzer.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TMath.h>
#include <TVirtualFFT.h>

#include <iostream>
#include <fstream>

using namespace std;

const Char_t* studyMappingFile = "mapping/20150818_uBooneChanMap.root";

const Bool_t runChirpStudyAlg = true;
const Bool_t runZigzagStudyAlg = false;
const Bool_t runWaveStudyAlg = false;
const Bool_t runFiltRMSStudyAlg = false;

const Bool_t saveWaveforms = false;
const Bool_t printHighNoiseChannels = false;
const Bool_t printLowNoiseChannels = false;

const Int_t maxChannels = 8300;
const Int_t maxChannelsElec = 8640;
const Int_t maxTicks = 9594;
const Int_t wireMaxNum[3] = {2400,2400,3456};

const Double_t filtrmsFactorMax = 1.2;
const Double_t filtrmsFactorMin = 0.95;
const Double_t sigAmpCut = 60.0;
const Double_t sigRMSCut = 5.0;
const Double_t sigFiltRMSCut = 2.0;
const Int_t sigBinMaxCut = 30;

void NoiseAnalyzer::RunAnalyzer()
{
  if(fFile == 0) return;

  avgMeanHist = new TH1F("avgMeanHist",";Channel;ADC Value Mean (Averaged Over Events)", maxChannels,-0.5,maxChannels-0.5);
  avgRMSHist = new TH1F("avgRMSHist", ";Channel;ADC Value RMS (Averaged Over Events)", maxChannels,-0.5,maxChannels-0.5);
  avgMaxHist = new TH1F("avgMaxHist", ";Channel;ADC Value Max-Mean (Averaged Over Events)", maxChannels,-0.5,maxChannels-0.5);
  chirpFracHist = new TH1F("chirpFracHist", ";Channel;Chirpiness (Averaged Over Events)", maxChannels,-0.5,maxChannels-0.5);
  zigzagFracHist = new TH1F("zigzagFracHist", ";Channel;Zigzagginess (Averaged Over Events)", maxChannels,-0.5,maxChannels-0.5);
  waveFracHist = new TH1F("waveFracHist", ";Channel;Waviness (Averaged Over Events)", maxChannels,-0.5,maxChannels-0.5);
  avgFiltRMSHist = new TH1F("avgFiltRMSHist", ";Channel;Post-filter ADC Value RMS (Averaged Over Events)", maxChannels,-0.5,maxChannels-0.5);
  lowCountFracHist = new TH1F("lowCountFracHist", ";Channel;Low RMS Channels (Averaged Over Events)", maxChannels,-0.5,maxChannels-0.5);
  high1CountFracHist = new TH1F("high1CountFracHist", ";Channel;High RMS Channels (Averaged Over Events)", maxChannels,-0.5,maxChannels-0.5);
  high2CountFracHist = new TH1F("high2CountFracHist", ";Channel;Signal Channels (Averaged Over Events)", maxChannels,-0.5,maxChannels-0.5);

  avgMeanElecHist = new TH1F("avgMeanElecHist",";ElecChannel;ADC Value Mean (Averaged Over Events)", maxChannelsElec,-0.5,maxChannelsElec-0.5);
  avgRMSElecHist = new TH1F("avgRMSElecHist", ";ElecChannel;ADC Value RMS (Averaged Over Events)", maxChannelsElec,-0.5,maxChannelsElec-0.5);
  avgMaxElecHist = new TH1F("avgMaxElecHist", ";ElecChannel;ADC Value Max-Mean (Averaged Over Events)", maxChannelsElec,-0.5,maxChannelsElec-0.5);
  chirpFracElecHist = new TH1F("chirpFracElecHist", ";ElecChannel;Chirpiness (Averaged Over Events)", maxChannelsElec,-0.5,maxChannelsElec-0.5);
  zigzagFracElecHist = new TH1F("zigzagFracElecHist", ";ElecChannel;Zigzagginess (Averaged Over Events)", maxChannelsElec,-0.5,maxChannelsElec-0.5);
  waveFracElecHist = new TH1F("waveFracElecHist", ";ElecChannel;Waviness (Averaged Over Events)", maxChannelsElec,-0.5,maxChannelsElec-0.5);
  avgFiltRMSElecHist = new TH1F("avgFiltRMSElecHist", ";ElecChannel;Post-filter ADC Value RMS (Averaged Over Events)", maxChannelsElec,-0.5,maxChannelsElec-0.5);
  lowCountFracElecHist = new TH1F("lowCountFracElecHist", ";ElecChannel;Low RMS Channels (Averaged Over Events)", maxChannelsElec,-0.5,maxChannelsElec-0.5);
  high1CountFracElecHist = new TH1F("high1CountFracElecHist", ";ElecChannel;High RMS Channels (Averaged Over Events)", maxChannelsElec,-0.5,maxChannelsElec-0.5);
  high2CountFracElecHist = new TH1F("high2CountFracElecHist", ";ElecChannel;Signal Channels (Averaged Over Events)", maxChannelsElec,-0.5,maxChannelsElec-0.5);

  TFile* mapFile = new TFile(studyMappingFile,"READ");
  TTree* mapTree = (TTree*) mapFile->Get("ubooneChanMap");

  Short_t varLarch[10000];
  Short_t varCrate[10000];
  Short_t varSlot[10000];
  Short_t varFemch[10000];
  Short_t varFeed[10000];
  Short_t varPlane[10000];
  Short_t varMboard[10000];
  Short_t varAsic[10000];
  Short_t varWirenum[10000];
  Short_t varLarwire[10000];

  Short_t larch;
  Short_t crate;
  Short_t slot;
  Short_t femch;
  Short_t feed;
  Short_t plane;
  Short_t mboard;
  Short_t asic;
  Short_t wirenum;
  Short_t larwire;

  mapTree->SetBranchAddress("larch", &larch);
  mapTree->SetBranchAddress("crate", &crate);
  mapTree->SetBranchAddress("slot", &slot);
  mapTree->SetBranchAddress("femch", &femch);
  mapTree->SetBranchAddress("feed", &feed);
  mapTree->SetBranchAddress("plane", &plane);
  mapTree->SetBranchAddress("mboard", &mboard);
  mapTree->SetBranchAddress("asic", &asic);
  mapTree->SetBranchAddress("wirenum", &wirenum);
  mapTree->SetBranchAddress("larwire", &larwire);

  Long64_t nentries = mapTree->GetEntries();
  Long64_t nbytes = 0, nb = 0;
  for(Long64_t jentry = 0; jentry < nentries; jentry++)
  {
    nb = mapTree->GetEntry(jentry);
    nbytes += nb;

    varLarch[larch] = larch;
    varCrate[larch] = crate;
    varSlot[larch] = slot;
    varFemch[larch] = femch;
    varFeed[larch] = feed;
    varPlane[larch] = plane;
    varMboard[larch] = mboard;
    varAsic[larch] = asic;
    varWirenum[larch] = wirenum;
    varLarwire[larch] = larwire;
  }
  mapFile->Close();

  TFile* outputFile = new TFile(Form("data/study/study_Run%d_Subrun%d_Event%d.root",fRunNumber,fSubrunNumber,fEventNumber),"RECREATE");

  ofstream outputText;
  if((printLowNoiseChannels == true) || (printHighNoiseChannels == true))
    outputText.open(Form("noiseSummary_Run%d_Subrun%d_Event%d.txt",fRunNumber,fSubrunNumber,fEventNumber));

  TH1F* meanHist = new TH1F("meanHist",";Channel;ADC Value Mean", maxChannels,-0.5,maxChannels-0.5);
  TH1F* rmsHist = new TH1F("rmsHist", ";Channel;ADC Value RMS", maxChannels,-0.5,maxChannels-0.5);
  TH1F* maxHist = new TH1F("maxHist", ";Channel;ADC Value Max-Mean", maxChannels,-0.5,maxChannels-0.5);
  TH1F* chirpHist = new TH1F("chirpHist", ";Channel;Chirpiness", maxChannels,-0.5,maxChannels-0.5);
  TH1F* zigzagHist = new TH1F("zigzagHist", ";Channel;Zigzagginess", maxChannels,-0.5,maxChannels-0.5);
  TH1F* waveHist = new TH1F("waveHist", ";Channel;Waviness", maxChannels,-0.5,maxChannels-0.5);
  TH1F* filtrmsHist = new TH1F("filtrmsHist", ";Channel;Post-filter ADC Value RMS", maxChannels,-0.5,maxChannels-0.5);
  TH1F* lowcountHist = new TH1F("lowcountHist", ";Channel;Low RMS Channels", maxChannels,-0.5,maxChannels-0.5);
  TH1F* high1countHist = new TH1F("high1countHist", ";Channel;High RMS Channels", maxChannels,-0.5,maxChannels-0.5);
  TH1F* high2countHist = new TH1F("high2countHist", ";Channel;Signal Channels", maxChannels,-0.5,maxChannels-0.5);

  TH1F* meanElecHist = new TH1F("meanElecHist",";ElecChannel;ADC Value Mean", maxChannelsElec,-0.5,maxChannelsElec-0.5);
  TH1F* rmsElecHist = new TH1F("rmsElecHist", ";ElecChannel;ADC Value RMS", maxChannelsElec,-0.5,maxChannelsElec-0.5);
  TH1F* maxElecHist = new TH1F("maxElecHist", ";ElecChannel;ADC Value Max-Mean", maxChannelsElec,-0.5,maxChannelsElec-0.5);
  TH1F* chirpElecHist = new TH1F("chirpElecHist", ";ElecChannel;Chirpiness", maxChannelsElec,-0.5,maxChannelsElec-0.5);
  TH1F* zigzagElecHist = new TH1F("zigzagElecHist", ";ElecChannel;Zigzagginess", maxChannelsElec,-0.5,maxChannelsElec-0.5);
  TH1F* waveElecHist = new TH1F("waveElecHist", ";ElecChannel;Waviness", maxChannelsElec,-0.5,maxChannelsElec-0.5);
  TH1F* filtrmsElecHist = new TH1F("filtrmsElecHist", ";ElecChannel;Post-filter ADC Value RMS", maxChannelsElec,-0.5,maxChannelsElec-0.5);
  TH1F* lowcountElecHist = new TH1F("lowcountElecHist", ";ElecChannel;Low RMS Channels", maxChannelsElec,-0.5,maxChannelsElec-0.5);
  TH1F* high1countElecHist = new TH1F("high1countElecHist", ";ElecChannel;High RMS Channels", maxChannelsElec,-0.5,maxChannelsElec-0.5);
  TH1F* high2countElecHist = new TH1F("high2countElecHist", ";ElecChannel;Signal Channels", maxChannelsElec,-0.5,maxChannelsElec-0.5);

  TH1F** noiseHist = new TH1F * [maxChannels];
  for(Int_t i = 0; i < maxChannels; i++)
  {
    noiseHist[i] = new TH1F(Form("noiseHist_channel%d",i),";Tick;ADC Value", maxTicks,-0.5,maxTicks-0.5);
  }

  Int_t channelNum;
  Int_t channelNumElec;
  Int_t loadVal = 0;
  Int_t chanNum = -1;
  for(Int_t planeNum = 0; planeNum <= 2; planeNum++)
  {
    for(Int_t wireNum = 0; wireNum < wireMaxNum[planeNum]; wireNum++)
    {
      chanNum++;
      loadVal = LoadHist(planeNum,wireNum);
      if(loadVal < 0) break;

      channelNum = chanNum;
      channelNumElec = 15*64*(varCrate[chanNum]-1)+64*(18-varSlot[chanNum])+varFemch[chanNum];

      pair<Double_t,Double_t> meanAndRMS = CalcMeanAndRMS();      
      Double_t meanVal = meanAndRMS.first;
      Double_t rmsVal = meanAndRMS.second;

      pair<Double_t,Double_t> extremeVals = GetExtrema();
      Double_t minVal = extremeVals.first;
      Double_t maxVal = extremeVals.second;
      
      meanHist->Fill(channelNum,meanVal);
      rmsHist->Fill(channelNum,rmsVal);
      maxHist->Fill(channelNum,maxVal-meanVal);
      meanElecHist->Fill(channelNumElec,meanVal);
      rmsElecHist->Fill(channelNumElec,rmsVal);
      maxElecHist->Fill(channelNumElec,maxVal-meanVal);
      
      avgMeanHist->Fill(channelNum,meanVal);
      avgRMSHist->Fill(channelNum,rmsVal);
      avgMaxHist->Fill(channelNum,maxVal-meanVal);
      avgMeanElecHist->Fill(channelNumElec,meanVal);
      avgRMSElecHist->Fill(channelNumElec,rmsVal);
      avgMaxElecHist->Fill(channelNumElec,maxVal-meanVal);
      
      if(saveWaveforms == true)
      {
        for(Int_t i = 0; i < fHist->GetNbinsX(); i++)
        {
          noiseHist[channelNum]->SetBinContent(i+1,fHist->GetBinContent(i+1));
        }
      }
      
      if(runChirpStudyAlg == true)
      {
        Double_t chirpVal = ChirpStudyAlg();
      
        if(chirpVal > 0.0)
        {
          chirpHist->Fill(channelNum,chirpVal);
          chirpElecHist->Fill(channelNumElec,chirpVal);
      
          chirpFracHist->Fill(channelNum,chirpVal);
          chirpFracElecHist->Fill(channelNumElec,chirpVal);
        }
      }
      
      if(runZigzagStudyAlg == true)
      {
        Double_t zigzagVal = ZigzagStudyAlg();
      
        if(zigzagVal > 0.0)
        {
          zigzagHist->Fill(channelNum,zigzagVal);
          zigzagElecHist->Fill(channelNumElec,zigzagVal);
      
          zigzagFracHist->Fill(channelNum,zigzagVal);
          zigzagFracElecHist->Fill(channelNumElec,zigzagVal);
        }
      }
      
      if(runWaveStudyAlg == true)
      {
        Double_t waveVal = WaveStudyAlg();
      
        if(waveVal > 0.0)
        {
          waveHist->Fill(channelNum,waveVal);
          waveElecHist->Fill(channelNumElec,waveVal);
      
          waveFracHist->Fill(channelNum,waveVal);
          waveFracElecHist->Fill(channelNumElec,waveVal);
        }
      }
      
      if(runFiltRMSStudyAlg == true)
      {
        Double_t filtrmsVal = FiltRMSStudyAlg(); // For noise runs
        //Double_t filtrmsVal = 0.0; // For pulser runs
      
        filtrmsHist->Fill(channelNum,filtrmsVal);
        filtrmsElecHist->Fill(channelNumElec,filtrmsVal);
      
        avgFiltRMSHist->Fill(channelNum,filtrmsVal);
        avgFiltRMSElecHist->Fill(channelNumElec,filtrmsVal);
      
        Double_t expectedFiltRMS;
        if(channelNum < 670.0)
          expectedFiltRMS = 0.63 + 0.64*((channelNum-0.0)/670.0);
        else if(channelNum < 1730.0)
          expectedFiltRMS = 1.27;
        else if(channelNum < 2400.0)
          expectedFiltRMS = 1.27 - 0.66*((channelNum-1730.0)/670.0);
        else if(channelNum < 3070.0)
          expectedFiltRMS = 0.61 + 0.73*((channelNum-2400.0)/670.0);
        else if(channelNum < 4130.0)
          expectedFiltRMS = 1.34;
        else if(channelNum < 4800.0)
          expectedFiltRMS = 1.34 - 0.72*((channelNum-4130.0)/670.0);
        else
          expectedFiltRMS = 0.93;
      
        Double_t expectedFiltRMSplaneMin;
        if(channelNum < 2400.0)
          expectedFiltRMSplaneMin = 0.63;
        else if(channelNum < 4800.0)
          expectedFiltRMSplaneMin = 0.61;
        else
          expectedFiltRMSplaneMin = 0.85; // Edited to prevent false positives on Y plane
      
        Double_t expectedFiltRMSplaneMax;
        if(channelNum < 2400.0)
          expectedFiltRMSplaneMax = 1.27;
        else if(channelNum < 4800.0)
          expectedFiltRMSplaneMax = 1.34;
        else
          expectedFiltRMSplaneMax = 0.93;
      
        Int_t sigBins = CountSigBins(meanVal); // For noise runs
        //Int_t sigBins = 0; // For pulser runs
      
        // ONCE CHANNEL MAPPING IS FIXED:  expectedFiltRMSplaneMin/Max --> expectedFiltRMS
        
        if((filtrmsVal > filtrmsFactorMax*expectedFiltRMSplaneMax) && (maxVal-minVal > sigAmpCut) && (rmsVal < sigRMSCut) && (filtrmsVal < sigFiltRMSCut) && (sigBins > 0) && (sigBins < sigBinMaxCut))
        {
          high2countHist->Fill(channelNum,1.0);
          high2countElecHist->Fill(channelNumElec,1.0);
      
          high2CountFracHist->Fill(channelNum,1.0);
          high2CountFracElecHist->Fill(channelNumElec,1.0);
      
          if(printHighNoiseChannels == true)
            outputText << "HIGH2:  " << channelNum << " " << filtrmsVal << " " << rmsVal << " " << maxVal-meanVal << endl;
        }
        else if(filtrmsVal > filtrmsFactorMax*expectedFiltRMSplaneMax)
        {
          high1countHist->Fill(channelNum,1.0);
          high1countElecHist->Fill(channelNumElec,1.0);
      
          high1CountFracHist->Fill(channelNum,1.0);
          high1CountFracElecHist->Fill(channelNumElec,1.0);
      
          if(printHighNoiseChannels == true)
            outputText << "HIGH1:  " << channelNum << " " << filtrmsVal << " " << rmsVal << " " << maxVal-meanVal << endl;
        }
      
        if(filtrmsVal < filtrmsFactorMin*expectedFiltRMSplaneMin)
        //if((filtrmsVal < filtrmsFactorMin*expectedFiltRMSplaneMin) && ((channelNum < 6528) || (channelNum > 6623))) // Hack to remove string of low RMS Y channels that still see pulses
        //if(maxVal-meanVal < 450.0) // Hack for low RMS channels in pulser runs
        {
          lowcountHist->Fill(channelNum,1.0);
          lowcountElecHist->Fill(channelNumElec,1.0);
      
          lowCountFracHist->Fill(channelNum,1.0);
          lowCountFracElecHist->Fill(channelNumElec,1.0);
      
          if(printLowNoiseChannels == true)
            outputText << "LOW:    " << channelNum << " " << filtrmsVal << " " << rmsVal << " " << maxVal-meanVal << endl;
        }
      }
    }

    if(loadVal < 0) break;
  }

  gStyle->SetOptStat(0);

  meanHist->SetOption("HIST");
  rmsHist->SetOption("HIST");
  maxHist->SetOption("HIST");
  meanHist->Write();
  rmsHist->Write();
  maxHist->Write();

  meanElecHist->SetOption("HIST");
  rmsElecHist->SetOption("HIST");
  maxElecHist->SetOption("HIST");
  meanElecHist->Write();
  rmsElecHist->Write();
  maxElecHist->Write();
  
  if(runChirpStudyAlg == true)
  {
    chirpHist->SetOption("HIST");
    chirpHist->Write();

    chirpElecHist->SetOption("HIST");
    chirpElecHist->Write();
  }

  if(runZigzagStudyAlg == true)
  {
    zigzagHist->SetOption("HIST");
    zigzagHist->Write();

    zigzagElecHist->SetOption("HIST");
    zigzagElecHist->Write();
  }

  if(runWaveStudyAlg == true)
  {
    waveHist->SetOption("HIST");
    waveHist->Write();

    waveElecHist->SetOption("HIST");
    waveElecHist->Write();
  }

  if(runFiltRMSStudyAlg == true)
  {
    filtrmsHist->SetOption("HIST");
    lowcountHist->SetOption("HIST");
    high1countHist->SetOption("HIST");
    high2countHist->SetOption("HIST");
    filtrmsHist->Write();
    lowcountHist->Write();
    high1countHist->Write();
    high2countHist->Write();

    filtrmsElecHist->SetOption("HIST");
    lowcountElecHist->SetOption("HIST");
    high1countElecHist->SetOption("HIST");
    high2countElecHist->SetOption("HIST");
    filtrmsElecHist->Write();
    lowcountElecHist->Write();
    high1countElecHist->Write();
    high2countElecHist->Write();
  }

  if(saveWaveforms == true)
  {
    for(Int_t i = 0; i < maxChannels; i++)
    {
      noiseHist[i]->SetOption("HIST");
      noiseHist[i]->Write();
    }
  }

  for(Int_t i = 0; i < maxChannels; i++)
  {
    delete noiseHist[i];
  }
  delete[] noiseHist;

  outputFile->Close();

  if((printLowNoiseChannels == true) || (printHighNoiseChannels == true))
    outputText.close();
}

pair<Double_t,Double_t> NoiseAnalyzer::CalcMeanAndRMS()
{
  pair<Double_t,Double_t> meanAndRMS;

  Double_t ADCval;
  Double_t theMean = 0.0;
  Double_t theRMS = 0.0;
  Int_t waveformSize = fHist->GetNbinsX();
  for(Int_t i = 0; i < waveformSize; i++)
  {
    ADCval = fHist->GetBinContent(i+1);
    theMean += ADCval;
    theRMS += TMath::Power(ADCval,2.0);
  }
  theMean /= (Double_t)waveformSize;
  theRMS /= (Double_t)waveformSize;
  theRMS = TMath::Sqrt(theRMS-TMath::Power(theMean,2.0));

  meanAndRMS.first = theMean;
  meanAndRMS.second = theRMS;

  return meanAndRMS;
}

pair<Double_t,Double_t> NoiseAnalyzer::GetExtrema()
{
  Double_t ADCval;
  Double_t minADC = 99999999.0;
  Double_t maxADC = -1.0;
  for(Int_t i = 0; i < fHist->GetNbinsX(); i++)
  {
    ADCval = fHist->GetBinContent(i+1);

    if(ADCval < minADC)
      minADC = ADCval;

    if(ADCval > maxADC)
      maxADC = ADCval;
  }

  pair<Double_t,Double_t> extremaADC;
  extremaADC.first = minADC;
  extremaADC.second = maxADC;

  return extremaADC;
}

Int_t NoiseAnalyzer::CountSigBins(Double_t meanVal)
{
  Double_t ADCval;
  Int_t sigBinCount = 0;

  for(Int_t i = 0; i < fHist->GetNbinsX(); i++)
  {
    ADCval = fHist->GetBinContent(i+1);

    if((ADCval >= meanVal+sigAmpCut/2.0) || (ADCval <= meanVal-sigAmpCut/2.0))
      sigBinCount++;
  }

  return sigBinCount;
}

Double_t NoiseAnalyzer::ChirpStudyAlg()
{
  const Int_t windowSize = 20;
  const Double_t chirpAmp = 0.0;
  //const Double_t chirpMinRMS = 0.55;
  const Double_t chirpMinRMS = 10000000.0;
  const Double_t chirpMaxRMS = 0.0;
  const Double_t chirpRatioMaxMin = 0.0;
  const Double_t chirpMinRMS_forNumLowCalc = 0.9;
  const Double_t maxNormalNeighborFrac = 0.20;

  Int_t counter = 0;
  Double_t ADCval;
  Double_t meanAmp = 0.0;
  Double_t runningAmpMean = 0.0;
  Double_t runningAmpRMS = 0.0;
  Double_t maxAmp = -1.0;
  Double_t maxRMS = -1.0;
  Double_t minRMS = 9999999.0;
  Int_t minRMSBin = -1;
  Int_t maxRMSBin = -1;
  Int_t numLowRMS = 0;
  Int_t firstLowRMSBin = -1;
  Int_t lastLowRMSBin = -1;
  Bool_t lowRMSFlag = false;
  Double_t RMSfirst = 0.0;
  Double_t RMSsecond = 0.0;
  Double_t RMSthird = 0.0;
  Int_t numNormalNeighbors = 0;
  for(Int_t i = 0; i < fHist->GetNbinsX(); i++)
  {
    ADCval = fHist->GetBinContent(i+1);
    meanAmp += ADCval;
    runningAmpMean += ADCval;
    runningAmpRMS += TMath::Power(ADCval,2.0);

    if(ADCval > maxAmp)
      maxAmp = ADCval;

    counter++;
    if(counter == windowSize)
    {
      runningAmpMean /= (Double_t)windowSize;
      runningAmpRMS /= (Double_t)windowSize;
      runningAmpRMS = TMath::Sqrt(runningAmpRMS-TMath::Power(runningAmpMean,2.0));

      RMSfirst = RMSsecond;
      RMSsecond = RMSthird;
      RMSthird = runningAmpRMS;

      if(runningAmpRMS > maxRMS)
      {
        maxRMS = runningAmpRMS;
        maxRMSBin = i-windowSize+1;
      }

      if(runningAmpRMS < minRMS)
      {
        minRMS = runningAmpRMS;
        minRMSBin = i-windowSize+1;
      }

      if(runningAmpRMS < chirpMinRMS_forNumLowCalc)
      {
        numLowRMS++;
        if(lowRMSFlag == false)
 	{
          lowRMSFlag = true;
          firstLowRMSBin = i-windowSize+1;
          lastLowRMSBin = i-windowSize+1;
	}
	else
	{
          lastLowRMSBin = i-windowSize+1;
	}
      }

      if(i >= 3*windowSize-1)
      {
        if((RMSsecond < chirpMinRMS_forNumLowCalc) && ((RMSfirst > chirpMinRMS_forNumLowCalc) || (RMSthird > chirpMinRMS_forNumLowCalc)))
          numNormalNeighbors++;
      }

      counter = 0;
      runningAmpMean = 0.0;
      runningAmpRMS = 0.0;
    }
  }
  meanAmp /= (Double_t)fHist->GetNbinsX();

  Double_t ratioMaxMin = maxRMS/minRMS;
  Double_t chirpFrac = ((Double_t) numLowRMS)/(((Double_t) maxTicks)/((Double_t) windowSize));
  Double_t normalNeighborFrac = ((Double_t) numNormalNeighbors)/((Double_t) numLowRMS);

  if((maxAmp-meanAmp > chirpAmp) && (minRMS < chirpMinRMS) && (maxRMS > chirpMaxRMS) && (ratioMaxMin > chirpRatioMaxMin) && ((normalNeighborFrac < maxNormalNeighborFrac) || ((numLowRMS < 2.0/maxNormalNeighborFrac) && (lastLowRMSBin-firstLowRMSBin == numLowRMS*windowSize))) && (numLowRMS > 4))
  {
    return chirpFrac;
  }
  else
  {
    return 0.0;
  }
}

Double_t NoiseAnalyzer::ZigzagStudyAlg()
{
  const Int_t startSidebandBin = 3700;
  const Int_t endSidebandBin = 4000;
  const Int_t startPeakBin = 4000;
  const Int_t endPeakBin = 4744;

  Double_t sidebandAvg;
  Double_t peakAvg;

  TH1F *currentHist = new TH1F("","",maxTicks,-0.5,maxTicks-0.5);
  for(Int_t i = 0; i < fHist->GetNbinsX(); i++)
  {
    currentHist->SetBinContent(i+1,fHist->GetBinContent(i+1));
  }

  TH1F *currentFFTHist = new TH1F("","",maxTicks,-0.5,maxTicks-0.5);
  currentFFTHist = (TH1F*)currentHist->FFT(currentFFTHist,"MAG");

  Double_t integral = currentFFTHist->Integral(startPeakBin,startPeakBin+50);
  Double_t maxIntegral = integral;
  for(int i = startPeakBin+1; i <= endPeakBin; i++)
  {
    integral -= currentFFTHist->GetBinContent(i-1);
    integral += currentFFTHist->GetBinContent(i+50);

    if(integral > maxIntegral)
      maxIntegral = integral;
  }
  peakAvg = maxIntegral/50.0;
  sidebandAvg = TMath::Max(0.0001,currentFFTHist->Integral(startSidebandBin,endSidebandBin)/((Double_t) endSidebandBin-startSidebandBin+1));

  delete currentHist;
  delete currentFFTHist;

  if(peakAvg/sidebandAvg > 0.0)
    return peakAvg/sidebandAvg;
  else
    return 0.0;
}

Double_t NoiseAnalyzer::WaveStudyAlg()
{
  const Int_t startSidebandBin = 3700;
  const Int_t endSidebandBin = 4000;
  const Int_t startPeakBin = 2;
  const Int_t endPeakBin = 248;

  Double_t sidebandAvg;
  Double_t peakAvg;

  TH1F *currentHist = new TH1F("","",maxTicks,-0.5,maxTicks-0.5);
  for(Int_t i = 0; i < fHist->GetNbinsX(); i++)
  {
    currentHist->SetBinContent(i+1,fHist->GetBinContent(i+1));
  }

  TH1F *currentFFTHist = new TH1F("","",maxTicks,-0.5,maxTicks-0.5);
  currentFFTHist = (TH1F*)currentHist->FFT(currentFFTHist,"MAG");

  Double_t integral = currentFFTHist->Integral(startPeakBin,startPeakBin+50);
  Double_t maxIntegral = integral;
  for(int i = startPeakBin+1; i <= endPeakBin; i++)
  {
    integral -= currentFFTHist->GetBinContent(i-1);
    integral += currentFFTHist->GetBinContent(i+50);

    if(integral > maxIntegral)
      maxIntegral = integral;
  }
  peakAvg = maxIntegral/50.0;
  sidebandAvg = TMath::Max(0.0001,currentFFTHist->Integral(startSidebandBin,endSidebandBin)/((Double_t) endSidebandBin-startSidebandBin+1));

  delete currentHist;
  delete currentFFTHist;

  if(peakAvg/sidebandAvg > 0.0)
    return peakAvg/sidebandAvg;
  else
    return 0.0;
}

Double_t NoiseAnalyzer::FiltRMSStudyAlg()
{
  const Int_t startFiltBin = 500;
  const Int_t endFiltBin = 3700;

  TH1F *currentHist = new TH1F("","",maxTicks,-0.5,maxTicks-0.5);
  for(Int_t i = 0; i < fHist->GetNbinsX(); i++)
  {
    currentHist->SetBinContent(i+1,fHist->GetBinContent(i+1));
  }

  Double_t *re = new Double_t[maxTicks];
  Double_t *im = new Double_t[maxTicks];
  Double_t *reFilt = new Double_t[maxTicks];
  Double_t *imFilt = new Double_t[maxTicks];

  TH1F *currentFFTHist = new TH1F("","",maxTicks,-0.5,maxTicks-0.5);
  currentFFTHist = (TH1F*)currentHist->FFT(currentFFTHist,"MAG");
  TVirtualFFT *currentFFTObject = TVirtualFFT::GetCurrentTransform();
  currentFFTObject->GetPointsComplex(re,im);

  Int_t waveformSize = fHist->GetNbinsX();
  for(Int_t i = 0; i < maxTicks; i++)
  {
    if(i+1 > waveformSize)
    {
      reFilt[i] = 0.0;
      imFilt[i] = 0.0;
    }
    else if((i+1 > startFiltBin) && (i+1 < endFiltBin))
    {
      reFilt[i] = re[i]/maxTicks;
      imFilt[i] = im[i]/maxTicks;
    }
    else
    {
      reFilt[i] = 0.0;
      imFilt[i] = 0.0;
    }
  }

  Int_t nFreqBins = maxTicks;
  TVirtualFFT *invCurrentFFTObject = TVirtualFFT::FFT(1,&nFreqBins,"C2R M K");
  invCurrentFFTObject->SetPointsComplex(reFilt,imFilt);
  invCurrentFFTObject->Transform();
  TH1F *newHist = new TH1F("","",maxTicks,-0.5,maxTicks-0.5);
  newHist = (TH1F*)TH1::TransformHisto(invCurrentFFTObject,newHist,"MAG");

  Double_t ADCval;
  Double_t filtMean = 0.0;
  Double_t filtRMS = 0.0;
  for(Int_t i = 0; i < waveformSize; i++)
  {
    ADCval = newHist->GetBinContent(i+1);
    filtMean += ADCval;
    filtRMS += TMath::Power(ADCval,2.0);
  }
  filtMean /= (Double_t)waveformSize;
  filtRMS /= (Double_t)waveformSize;
  filtRMS = TMath::Sqrt(filtRMS-TMath::Power(filtMean,2.0));

  delete currentHist;
  delete currentFFTHist;
  delete newHist;
  delete currentFFTObject;
  delete invCurrentFFTObject;
  delete[] re;
  delete[] im;
  delete[] reFilt;
  delete[] imFilt;

  return filtRMS;
}

vector<TH1F*> NoiseAnalyzer::GetHists()
{
  vector<TH1F*> histVec;

  histVec.push_back((TH1F*)avgMeanHist->Clone());
  histVec.push_back((TH1F*)avgRMSHist->Clone());
  histVec.push_back((TH1F*)avgMaxHist->Clone());
  histVec.push_back((TH1F*)chirpFracHist->Clone());
  histVec.push_back((TH1F*)zigzagFracHist->Clone());
  histVec.push_back((TH1F*)waveFracHist->Clone());
  histVec.push_back((TH1F*)avgFiltRMSHist->Clone());
  histVec.push_back((TH1F*)lowCountFracHist->Clone());
  histVec.push_back((TH1F*)high1CountFracHist->Clone());
  histVec.push_back((TH1F*)high2CountFracHist->Clone());

  histVec.push_back((TH1F*)avgMeanElecHist->Clone());
  histVec.push_back((TH1F*)avgRMSElecHist->Clone());
  histVec.push_back((TH1F*)avgMaxElecHist->Clone());
  histVec.push_back((TH1F*)chirpFracElecHist->Clone());
  histVec.push_back((TH1F*)zigzagFracElecHist->Clone());
  histVec.push_back((TH1F*)waveFracElecHist->Clone());
  histVec.push_back((TH1F*)avgFiltRMSElecHist->Clone());
  histVec.push_back((TH1F*)lowCountFracElecHist->Clone());
  histVec.push_back((TH1F*)high1CountFracElecHist->Clone());
  histVec.push_back((TH1F*)high2CountFracElecHist->Clone());

  return histVec;
}
