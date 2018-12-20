#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TMath.h"
#include "TCanvas.h"
#include "Function.h"

////////////////////////////////////////////////////////////////////////////////////////////////////

void Energy_MultiCal_B()
{
//for(int numb = 0 ; numb < 1 ; numb++)
// {
  const int numb = 11;
  gStyle -> SetOptFit();

  TChain* chain = new TChain("EventTree");
  TString filename;
  ifstream filelist("CosmicRay.txt");

  while(filelist >> filename)
	chain -> Add(filename);


  //TFile* filename = new TFile("RUN_772_event_1.root","read");
  //TTree* tree = (TTree*)filename -> Get("EventTree");

  TFile* file = new TFile("Position_MultiOffset_B.root","read");
  TTree* treeoff = (TTree*)file -> Get("offsetVal");

  double offsetvalue[24];

  treeoff -> SetBranchAddress("offsetvalue",offsetvalue);

  //Wave Form Histogram
  TH1D* histR[24];
  TH1D* histL[24];
  //TH2D* ADC_Pos[24];
  //TH2D* QDC_Pos[24];

  for(int i = 0 ; i < 24 ; i++)
  { 
	histR[i] = new TH1D(Form("histR%d",i),Form("Waveform_Right_%d_Bar",i),240,0,240);
	histL[i] = new TH1D(Form("histL%d",i),Form("Waveform_Left_%d_Bar",i),240,0,240);

	//ADC_Pos[i] = new TH2D(Form("ADC_Pos%d",i),Form("ADC vs Position %d Bar",i),150,-150,150,400,0,800);
	//QDC_Pos[i] = new TH2D(Form("QDC_Pos%d",i),Form("QDC vs Position %d Bar",i),200,-200,200,4000,0,8000);
  }  
  ///////////////////////////////////////////////////////////////// 
	TFile* OutPut = new TFile(Form("Bar%iposition_QDC_B.root",numb),"recreate");//change
  /* TCanvas*  can1[6];
	 for(int i = 0 ; i < 6 ; i++)
	 { 
	 can1[i] = new TCanvas(Form(" can1%d",i),"",1000,1000);
	 can1[i] ->Divide(2,2);
	 } */

  TCanvas*  histqdc[7];
  for(int i = 0 ; i < 7 ; i++)
  { 
	histqdc[i] = new TCanvas(Form(" histqdc%d",i),"",1000,1000);
	histqdc[i] ->Divide(2,2);
  } 


  TH1D* qdcPos[25];
  for(int i = 0 ; i < 25 ; i++)
  {
	qdcPos[i] = new TH1D(Form("qdcPos%d",i),Form("qdc_position From %d to %d [cm]",-100+(i*8),-100+((i+1)*8) ),100,0,8000);
  }  


  //data from root file
  const int Mbin = 240; 
  const int MWAVE  = 93;

  int ient;

  int nWAVE; 
  UShort_t WAVE[MWAVE][Mbin];
  UShort_t Channel[MWAVE];

  chain -> SetBranchAddress("nWAVE",&nWAVE);
  chain -> SetBranchAddress("WAVE",WAVE);
  chain -> SetBranchAddress("Channel",Channel);


  //for using function

  double Cross_PR[24];
  double Cross_PL[24]; 
  //double bar_ADC_R[24];
  //double bar_ADC_L[24];
  double bar_QDC_R[25];
  double bar_QDC_L[25];

  double QDC_Geo[24];
  //double QDC_Geo[25];
  double Cross_diff[24];
  double Modi_Cross_diff[24];

  int entries = chain -> GetEntries();

  //cout<<entries<<endl;
  double checkP[26];
  for(int i = 0 ; i <26 ; i++)
  { 
	checkP[i] = -100. + (8*i);
  }  

  for(Int_t iEntry = 0 ; iEntry < entries ; iEntry++)
  {	
	if(iEntry == 215703) continue;
	if(iEntry == 259896) continue;
	if(iEntry == 2516573) continue;
	if(iEntry == 3403010) continue;
	if(iEntry == 7004675) continue;

	chain -> GetEntry(iEntry);

	if(nWAVE >1)
	{
	  for(int nwave = 0 ; nwave < nWAVE ; nwave++)
	  {
		// WAVE
		if( 56 <= Channel[nwave] &&  Channel[nwave] <= 103)
		{
		  for(int i = 0 ; i < nwave ; i++)
		  {
			for(int j = 0 ; j < nwave ; j++)
			{
			  if( 56 <= Channel[i] &&  Channel[i] <= 79 && 80 <= Channel[j] &&  Channel[j] <= 103)
			  {
				if(Channel[i] +24 == Channel[j])
				{

				  for(int h = 0 ; h < 24 ; h++)
				  {
					if(Channel[i] == 79 - h )
					{
					  for(int bin = 0 ; bin < Mbin ; bin++)
					  {
						histR[h] -> SetBinContent(bin + 1, WAVE[i][bin]);
					  }
					  Cross_PR[h] = TDC(Mbin,histR[h]);
					  bar_QDC_R[h] = QDC(Mbin,histR[h]);
					}

					if(Channel[j] == 103 - h )
					{
					  for(int bin = 0 ; bin < Mbin ; bin++)
					  {
						histL[h] -> SetBinContent(bin + 1, WAVE[j][bin]);
					  }
					  Cross_PL[h] = TDC(Mbin,histL[h]);
					  bar_QDC_L[h] = QDC(Mbin,histL[h]);
					}

					treeoff -> GetEntry(h);
					Cross_diff[h] =  (Cross_PL[h] - Cross_PR[h])*7.7;
					Modi_Cross_diff[h] = Cross_diff[h]-offsetvalue[h];

					QDC_Geo[h] = TMath::Sqrt(bar_QDC_R[h] * bar_QDC_L[h]);
					//ADC_Pos[h] -> Fill(Modi_Cross_diff[h],ADC_Geo[h]) ;
				  }

				  for(int k = 0 ; k < 25; k++)
				  {
					if( checkP[k] <= Modi_Cross_diff[numb] && Modi_Cross_diff[numb] <= checkP[k+1])//change
					{ 
					  qdcPos[k] -> Fill(QDC_Geo[numb]);//change 
					}
				  }
				}
			  }
			}
		  }
		}
	  }
	}
  }
  
  //////////////////////////////////////////
  for(int k = 0 ; k < 25; k++)
  {
	int cannumber = k/4;

	histqdc[cannumber] -> cd(k%4+1);
	//adcPos[k] -> Fit("landau","","",130,800);
	qdcPos[k] -> Draw();
	qdcPos[k] -> Write();
  }
  /////////////////////////////////////
  /*for(int k = 0 ; k < 24 ; k++)
	{
	int canumber = k/4;

	can1[canumber] -> cd(k%4+1);
	can1[canumber] -> SetLogz();
	ADC_Pos[k] -> Draw("colz");

  // pan1[canumber] -> cd(k%4+1);
  // QDC_GeoMean[k] -> Draw();
  }*/

}
//}
