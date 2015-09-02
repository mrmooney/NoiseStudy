#include "NoiseAnalyzer.h"
#include "TH2.h"
#include "TFile.h"
#include "TStyle.h"

vector<TH1F*> studyEvent(Int_t runNum, Int_t subrunNum, Int_t eventNum);

int main(int argc, char **argv)
{
  const Int_t runMode = 2;

  const Int_t maxChannels = 8300;
  const Int_t maxChannelsElec = 8640;

  const Int_t minRunNum = 1532;
  const Int_t maxRunNum = 1532;

  const Int_t minEventNum = 0;
  const Int_t maxEventNum = 9;

  Char_t *histString = (Char_t*) "";
  Int_t minVal;
  Int_t maxVal;
  if(runMode == 1)
  {
    histString = (Char_t*) "Run";
    minVal = minRunNum;
    maxVal = maxRunNum;
  }
  else if(runMode == 2)
  {
    histString = (Char_t*) "Event";
    minVal = minEventNum+1;
    maxVal = maxEventNum+1;
  }

  TFile* outputFile2D = new TFile("StudyNoiseHists.root","RECREATE");

  TH2F *meanHist2D = new TH2F("meanHist2D", Form("ADC Value Mean;Channel;%s #",histString), maxChannels,-0.5,maxChannels-0.5,maxVal-minVal+1,minVal-0.5,maxVal+0.5);
  TH2F *rmsHist2D = new TH2F("rmsHist2D", Form("ADC Value RMS;Channel;%s #",histString), maxChannels,-0.5,maxChannels-0.5,maxVal-minVal+1,minVal-0.5,maxVal+0.5);
  TH2F *maxHist2D = new TH2F("maxHist2D", Form("ADC Value Max-Mean;Channel;%s #",histString), maxChannels,-0.5,maxChannels-0.5,maxVal-minVal+1,minVal-0.5,maxVal+0.5);
  TH2F *chirpHist2D = new TH2F("chirpHist2D", Form("Chirpiness;Channel;%s #",histString), maxChannels,-0.5,maxChannels-0.5,maxVal-minVal+1,minVal-0.5,maxVal+0.5);
  TH2F *zigzagHist2D = new TH2F("zigzagHist2D", Form("Zigzagginess;Channel;%s #",histString), maxChannels,-0.5,maxChannels-0.5,maxVal-minVal+1,minVal-0.5,maxVal+0.5);
  TH2F *waveHist2D = new TH2F("waveHist2D", Form("Waviness;Channel;%s #",histString), maxChannels,-0.5,maxChannels-0.5,maxVal-minVal+1,minVal-0.5,maxVal+0.5);
  TH2F *filtrmsHist2D = new TH2F("filtrmsHist2D", Form("Post-filter ADC Value RMS;Channel;%s #",histString), maxChannels,-0.5,maxChannels-0.5,maxVal-minVal+1,minVal-0.5,maxVal+0.5);
  TH2F *lowcountHist2D = new TH2F("lowcountHist2D", Form("Low RMS Channels;Channel;%s #",histString), maxChannels,-0.5,maxChannels-0.5,maxVal-minVal+1,minVal-0.5,maxVal+0.5);
  TH2F *high1countHist2D = new TH2F("high1countHist2D", Form("High RMS Channels;Channel;%s #",histString), maxChannels,-0.5,maxChannels-0.5,maxVal-minVal+1,minVal-0.5,maxVal+0.5);
  TH2F *high2countHist2D = new TH2F("high2countHist2D", Form("Signal Channels;Channel;%s #",histString), maxChannels,-0.5,maxChannels-0.5,maxVal-minVal+1,minVal-0.5,maxVal+0.5);

  TH2F *meanElecHist2D = new TH2F("meanElecHist2D", Form("ADC Value Mean;ElecChannel;%s #",histString), maxChannelsElec,-0.5,maxChannelsElec-0.5,maxVal-minVal+1,minVal-0.5,maxVal+0.5);
  TH2F *rmsElecHist2D = new TH2F("rmsElecHist2D", Form("ADC Value RMS;ElecChannel;%s #",histString), maxChannelsElec,-0.5,maxChannelsElec-0.5,maxVal-minVal+1,minVal-0.5,maxVal+0.5);
  TH2F *maxElecHist2D = new TH2F("maxElecHist2D", Form("ADC Value Max-Mean;ElecChannel;%s #",histString), maxChannelsElec,-0.5,maxChannelsElec-0.5,maxVal-minVal+1,minVal-0.5,maxVal+0.5);
  TH2F *chirpElecHist2D = new TH2F("chirpElecHist2D", Form("Chirpiness;ElecChannel;%s #",histString), maxChannelsElec,-0.5,maxChannelsElec-0.5,maxVal-minVal+1,minVal-0.5,maxVal+0.5);
  TH2F *zigzagElecHist2D = new TH2F("zigzagElecHist2D", Form("Zigzagginess;ElecChannel;%s #",histString), maxChannelsElec,-0.5,maxChannelsElec-0.5,maxVal-minVal+1,minVal-0.5,maxVal+0.5);
  TH2F *waveElecHist2D = new TH2F("waveElecHist2D", Form("Waviness;ElecChannel;%s #",histString), maxChannelsElec,-0.5,maxChannelsElec-0.5,maxVal-minVal+1,minVal-0.5,maxVal+0.5);
  TH2F *filtrmsElecHist2D = new TH2F("filtrmsElecHist2D", Form("Post-filter ADC Value RMS;ElecChannel;%s #",histString), maxChannelsElec,-0.5,maxChannelsElec-0.5,maxVal-minVal+1,minVal-0.5,maxVal+0.5);
  TH2F *lowcountElecHist2D = new TH2F("lowcountElecHist2D", Form("Low RMS Channels;ElecChannel;%s #",histString), maxChannelsElec,-0.5,maxChannelsElec-0.5,maxVal-minVal+1,minVal-0.5,maxVal+0.5);
  TH2F *high1countElecHist2D = new TH2F("high1countElecHist2D", Form("High RMS Channels;ElecChannel;%s #",histString), maxChannelsElec,-0.5,maxChannelsElec-0.5,maxVal-minVal+1,minVal-0.5,maxVal+0.5);
  TH2F *high2countElecHist2D = new TH2F("high2countElecHist2D", Form("Signal Channels;ElecChannel;%s #",histString), maxChannelsElec,-0.5,maxChannelsElec-0.5,maxVal-minVal+1,minVal-0.5,maxVal+0.5);

  vector<TH1F*> histSet;
  for(Int_t i = minVal; i <= maxVal; i++)
  {
    histSet.clear();

    if(runMode == 1)
    {
      cout << "Analyzing Run #" << i << endl;
      histSet = studyEvent(i,0,0);
    }
    else if(runMode == 2)
    {
      cout << "Analyzing Event #" << i << endl;
      histSet = studyEvent(minRunNum,0,i-1);
    }

    for(Int_t j = 0; j < maxChannels; j++)
    {
      meanHist2D->Fill(j,i,histSet[0]->GetBinContent(j+1));
      rmsHist2D->Fill(j,i,histSet[1]->GetBinContent(j+1));
      maxHist2D->Fill(j,i,histSet[2]->GetBinContent(j+1));
      chirpHist2D->Fill(j,i,histSet[3]->GetBinContent(j+1));
      zigzagHist2D->Fill(j,i,histSet[4]->GetBinContent(j+1));
      waveHist2D->Fill(j,i,histSet[5]->GetBinContent(j+1));
      filtrmsHist2D->Fill(j,i,histSet[6]->GetBinContent(j+1));
      lowcountHist2D->Fill(j,i,histSet[7]->GetBinContent(j+1));
      high1countHist2D->Fill(j,i,histSet[8]->GetBinContent(j+1));
      high2countHist2D->Fill(j,i,histSet[9]->GetBinContent(j+1));

      meanElecHist2D->Fill(j,i,histSet[10]->GetBinContent(j+1));
      rmsElecHist2D->Fill(j,i,histSet[11]->GetBinContent(j+1));
      maxElecHist2D->Fill(j,i,histSet[12]->GetBinContent(j+1));
      chirpElecHist2D->Fill(j,i,histSet[13]->GetBinContent(j+1));
      zigzagElecHist2D->Fill(j,i,histSet[14]->GetBinContent(j+1));
      waveElecHist2D->Fill(j,i,histSet[15]->GetBinContent(j+1));
      filtrmsElecHist2D->Fill(j,i,histSet[16]->GetBinContent(j+1));
      lowcountElecHist2D->Fill(j,i,histSet[17]->GetBinContent(j+1));
      high1countElecHist2D->Fill(j,i,histSet[18]->GetBinContent(j+1));
      high2countElecHist2D->Fill(j,i,histSet[19]->GetBinContent(j+1));
    }
  }

  outputFile2D->cd();

  meanHist2D->Write();
  rmsHist2D->Write();
  maxHist2D->Write();
  chirpHist2D->Write();
  zigzagHist2D->Write();
  waveHist2D->Write();
  filtrmsHist2D->Write();
  lowcountHist2D->Write();
  high1countHist2D->Write();
  high2countHist2D->Write();

  meanElecHist2D->Write();
  rmsElecHist2D->Write();
  maxElecHist2D->Write();
  chirpElecHist2D->Write();
  zigzagElecHist2D->Write();
  waveElecHist2D->Write();
  filtrmsElecHist2D->Write();
  lowcountElecHist2D->Write();
  high1countElecHist2D->Write();
  high2countElecHist2D->Write();

  outputFile2D->Close();

  return 0;
}

vector<TH1F*> studyEvent(Int_t runNum, Int_t subrunNum, Int_t eventNum)
{
  NoiseAnalyzer *noiseAna = new NoiseAnalyzer(runNum,subrunNum,eventNum);
  noiseAna->RunAnalyzer(); 

  vector<TH1F*> hists = noiseAna->GetHists();
  
  delete noiseAna;
  return hists;
}
