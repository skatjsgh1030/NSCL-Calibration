#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TF1.h"
#include "TMath.h"
#include "TCanvas.h"
#include "Function.h"

void TDC_Multi_B()
{
    TChain* chain = new TChain("EventTree");
	  TString filename;
	  ifstream filelist("CosmicRay.txt");

	  while(filelist >> filename)
	  chain -> Add(filename);
  //TFile* file = new TFile("RUN_772_event_1.root","read");
  //TTree* tree = (TTree*)file -> Get("EventTree");

  //Wave Form Histogram
  TH1D* histR[24];
  TH1D* histL[24];
  TH1D* tdcLR[24];
  TH1D* tdcR[24];
  TH1D* tdcL[24];
  TH1D* TDC_diff[24];
  for(int k = 0 ; k < 24 ; k++)
  { 
	histR[k] = new TH1D(Form("histR%d",k),Form("Waveform_Right_%d_Bar",k),240,0,240);
	histL[k] = new TH1D(Form("histL%d",k),Form("Waveform_Left_%d_Bar",k),240,0,240);
	tdcLR[k] = new TH1D(Form("tdcLR%d",k),Form("%d Bar Position",k),300,-300,300);
	tdcR[k] = new TH1D(Form("tdcR%d",k),Form("TDC_diff_Right_%d_Bar",k),150,50,200);
	tdcL[k] = new TH1D(Form("tdcL%d",k),Form("TDC_diff_Left_%d_Bar",k),150,50,200);
	TDC_diff[k] = new TH1D(Form("TDC_diff%d",k),Form("TDC_diff_%d_Bar",k),300,-300,300);
  }  

  TFile* Output = new TFile("TDC_diff_Multi_B.root","recreate");
  
  // 2D Histogram for arrangment
  //TH2D* track =  new TH2D("track","Track of particle",300,-300,300,25,0,25);

  //for Out TDC difference Histogram

  TCanvas* pan[6];
  for(int i = 0 ; i<6 ; i++)
  {	
	pan[i] = new TCanvas(Form("pan%d",i),"",1000,1000);
	pan[i] -> Divide(2,2);
  }

  //data from root file
  //////////////////////////////////////////////////////////////////////////////////////

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

  double Cross_diff[24];
  //double TDC_diff[24];

  int entries = chain -> GetEntries();

  for(Int_t iEntry = 0 ; iEntry < entries ; iEntry++)
  {	
	if(iEntry == 215703) continue;
	if(iEntry == 259896) continue;
	if(iEntry == 2516573) continue;
	if(iEntry == 3403010) continue;
	if(iEntry == 7004675) continue;

	chain -> GetEntry(iEntry);

	//cout <<"Entry :"<< iEntry << endl;//check entry
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
				  Cross_PR[h] =  -9999;
				  Cross_PL[h] =  -9999;
					if(Channel[i] == 79 - h )
					{
					  //cout << Channel[i] << endl; 
					  for(int bin = 0 ; bin < Mbin ; bin++)
					  {	
						histR[h] -> SetBinContent(bin + 1, WAVE[i][bin]);
					  }	
					  Cross_PR[h] = TDC(Mbin,histR[h]);
					  tdcR[h] -> Fill(Cross_PR[h]);
					  histR[h] -> Reset();
					}

					if(Channel[j] == 103 - h )
					{
					  //cout << Channel[j] << endl; 
					  for(int bin = 0 ; bin < Mbin ; bin++)
					  {	
						histL[h] -> SetBinContent(bin + 1, WAVE[j][bin]);
					  }	
					  Cross_PL[h] = TDC(Mbin,histL[h]);
					  tdcL[h] -> Fill(Cross_PL[h]);
					  histL[h] -> Reset();
					}
					Cross_diff[h] =  (Cross_PL[h] - Cross_PR[h])*7.7;
					if(Cross_PR[h] == -9999 ||Cross_PL[h] == -9999) continue;
					tdcLR[h] -> Fill(Cross_diff[h]);
					TDC_diff[h] -> Fill(Cross_diff[h]);
				  }
				}
			  }
			}
		  }
		}
	  }
	}
  }

  for(int i = 0 ; i < 24 ; i++)
  {
	int canumber = i/4;
	pan[canumber] -> cd(i%4+1);
	tdcLR[i] -> Draw();
	
	TDC_diff[i] -> Write();
  }
}

