#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TMath.h"

void PositionQDC_fit_B()
{
  gStyle -> SetOptFit();
	const int num = 11;
  TFile* file = new TFile(Form("Bar%iposition_QDC_B.root",num),"read");

  TFile* fitPeak = new TFile(Form("Bar%iQDC_value_B.root",num),"recreate");
  TTree* tree = new TTree("fitPeak","PeakofQDC");


  //double peakval1[25];
  double peakval[25];
  double Chi2val[25];
  double ndfval[25];
  double fitPrecise[25];
  //double peakval3[25];
  double peakvalue[25];
  double Chi2value[25];
  double NDFvalue[25];
  
  tree -> Branch("peakvalue",&peakvalue,"peakvalue[25]/D");
  tree -> Branch("Chi2value",&Chi2value,"Chi2value[25]/D");
  tree -> Branch("NDFvalue",&NDFvalue,"NDFvalue[25]/D");

 /* TCanvas* c1[7];
  for(int i =0 ; i<7 ;i++ )
  {	
	c1[i] = new TCanvas(Form("c1%d",i),"",1000,1000);
	c1[i]->Divide(2,2);
  }*/

  TCanvas* c2[25];
  for(int i = 0 ; i<25 ;i++ )
  {	
	c2[i] = new TCanvas(Form("c2%d",i),"",800,800);
	c2[i]->Divide(1,1);
  }	

  TH1D* qdcPos_R[25];
  for(int j = 1 ; j < 24 ; j++)
  {	
	qdcPos_R[j] = (TH1D*) file -> Get(Form("qdcPos%d",j));

	int cnu = j/4;
	int ccnu = j/1;
	

	/*c1[cnu] -> cd(j%4+1);
	qdcPos_R[j] -> Fit("landau","q","",(0.00965j*j)-(0.01751*j)+1100,(0.00965j*j)-(0.01751*j)+4900);
	qdcPos_R[j] -> Draw();*/
	
	c2[ccnu] ->cd(j%1+1);
	qdcPos_R[j] -> Fit("landau","q","",(0.05266j*j)-(0.14492*j)+1100,(0.05266j*j)-(0.14492*j)+4900);
	qdcPos_R[j] -> Draw();
	//c2[ccnu] -> SaveAs(Form("positionQDC%i.pdf",ccnu));
	
	auto landau = (TF1*) gROOT -> GetFunction("landau");
	peakval[j] = landau -> GetParameter(1);
	Chi2val[j] = landau -> GetChisquare();
	ndfval[j] = landau -> GetNDF();
	fitPrecise[j] = Chi2val[j]/ndfval[j];
	
	const int i = -92;
	const int k = 92;
	peakvalue[0] = (0.05266*i*i)-(0.14492*i)+1373;
	peakvalue[24] = (0.05266*k*k)-(0.14492*k)+1373;
	
	peakvalue[j] = peakval[j];
	Chi2value[j] = Chi2val[j];
	NDFvalue[j] = ndfval[j];
	tree -> Fill();
	
	//cout << Form( "peak value %d : ",j)<< peakval[j] << endl;
	//cout << Form( "Chi2 %d : ",j)<< Chi2val[j] << endl;
	//cout << Form( "NDF %d : ",j)<< ndfval[j] << endl;
	//cout << Form( "***Chi^2 / ndf %d : ",j)<< fitPrecise[j]<<"***" << endl;

  }	
  	
  for(int j = 0 ; j < 25 ; j++)
  {	
	cout << Form( "peak value %d : ",j)<< peakvalue[j] << endl;
  }
	tree -> Write();
	fitPeak ->Close();
}
