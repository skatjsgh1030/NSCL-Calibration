#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "TTree.h"
#include "TF1.h"
#define upper 0.440//0.440
#define lower 0.0275//0.0275

/////////////////////////////////////////////////////////////////////////////////////////////
//cross point of Top side of TDC differ Histogram
/////////////////////////////////////////////////////////////////////////////////////

double U_CrossPoint (int range, TH1D* hist)
{
        double U_Level = (hist -> GetMaximum())*upper; 
        double U_position = 0;
        int n = 0 ; 
        for(n = -300 ; n < range; n++)
        {   
                if(U_Level > hist -> GetBinContent(n+1) && U_Level < hist -> GetBinContent(n+2))
                {   
                        U_position =(U_Level - (hist -> GetBinContent(n+1) - (hist -> GetBinContent(n+2) - hist -> GetBinContent(n+1))*(n+1))) / ( hist -> GetBinContent(n+2) - hist -> GetBinContent(n+1) );
                        break;
                }   
                else if(U_Level == hist -> GetBinContent(n+1))
                {   
                        U_position = (n+1);
                }   
        }   
        return U_position;
}

double D_CrossPoint (int range, TH1D* hist)
{
        double D_Level = (hist -> GetMaximum())*upper; 
        double D_position = 0;
        int n = 0 ; 
        for(n = 300 ; n > range; n--)
        {   
                if(D_Level < hist -> GetBinContent(n+1) && D_Level > hist -> GetBinContent(n+2))
                {   
                        D_position =(D_Level - (hist -> GetBinContent(n+1) - (hist -> GetBinContent(n+2) - hist -> GetBinContent(n+1))*(n+1))) / ( hist -> GetBinContent(n+2) - hist -> GetBinContent(n+1) );
                        break;
                }   
                else if(D_Level == hist -> GetBinContent(n+1))
                {   
                        D_position = (n+1);
                }   
        }   
        return D_position;
}


///////////////////////////////////////////////////////////////////////////////////////////////////
//cross point of Bottom side of TDC differ Histogram
////////////////////////////////////////////////////////////////////////////////////////////
double UL_CrossPoint (int range, TH1D* hist)
{
        double D_Level = (hist -> GetMaximum())*lower; 
        double D_position = 0;
        int n = 0 ;
        for(n = -300 ; n < range; n++)
        {
                if(D_Level > hist -> GetBinContent(n+1) && D_Level < hist -> GetBinContent(n+2))
                {
                        D_position =(D_Level - (hist -> GetBinContent(n+1) - (hist -> GetBinContent(n+2) - hist -> GetBinContent(n+1))*(n+1))) / ( hist -> GetBinContent(n+2) - hist -> GetBinContent(n+1) );
                        break;
                }
                else if(D_Level == hist -> GetBinContent(n+1))
                {
                        D_position = (n+1);
                }
        }
        return D_position;
}

double DL_CrossPoint (int range, TH1D* hist)
{
        double DL_Level = (hist -> GetMaximum())*lower;
        double DL_position = 0;
        int n = 0 ;
        for(n = 300 ; n > range; n--)
        {
                if(DL_Level < hist -> GetBinContent(n+1) && DL_Level > hist -> GetBinContent(n+2))
                {
                        DL_position =(DL_Level - (hist -> GetBinContent(n+1) - (hist -> GetBinContent(n+2) - hist -> GetBinContent(n+1))*(n+1))) / ( hist -> GetBinContent(n+2) - hist -> GetBinContent(n+1) );
                        break;
                }
                else if(DL_Level == hist -> GetBinContent(n+1))
                {
                        DL_position = (n+1);
                }
        }
        return DL_position;
}
//////////////////////////////////////////////////////////////////////


void TDC_Multi_fit_B ()
{ 
        TFile* file = new TFile("TDC_diff_Multi_B.root","read");  


        //TFile* filename = new TFile("RUN_772_event_1.root","read");
        //TTree* tree = (TTree*)filename -> Get("EventTree");  

        TF1* UpH_fit[24];
        TF1* DpH_fit[24];

        ////////////////////////////////////////////////////////////////////////////////
        TCanvas* can[24];
        for(int j =0 ; j < 24 ; j++ )
        { 
                can[j] = new TCanvas(Form("can%d",j),"",1000,1000);
                can[j] -> Divide(1,1);
        }

        TH1D* TDC_diff[24];
        for(int i = 0 ; i < 24 ; i++)
        { 
                TDC_diff[i]  = (TH1D*) file -> Get(Form("TDC_diff%d",i));
        }
        
        /////////////////////////////////////////////////////////////
        TFile* offset = new TFile("Position_MultiOffset_B.root","recreate");
        TTree* treeOff = new TTree("offsetVal","Offset_value");
        /////////////////////////////////////////////////////////////
        
        //make offset value tree//
        double offsetvalue[24];
        treeOff -> Branch("offsetvalue",&offsetvalue,"offsetvalue[24]/D");

        //for making Fitting Function
        ///////////////////////////////////////////////////////////////////////////////////

        double UCPoT[24];
        double DCPoT[24];
        double ULCPoT[24];
        double DLCPoT[24];
        double Utilt[24];
        double Dtilt[24];

        double Leftvalue[24];
        double Rightvalue[24];
        double Offset_val[24];

        for(int i = 0 ; i < 24 ; i++)
        { 
                UCPoT[i] = U_CrossPoint(300,TDC_diff[i]); 
                DCPoT[i] = D_CrossPoint(-300,TDC_diff[i]);  
                ULCPoT[i] = UL_CrossPoint(300,TDC_diff[i]); 
                DLCPoT[i] = DL_CrossPoint(-300,TDC_diff[i]); 
                Utilt[i] = 1000/(UCPoT[i]-ULCPoT[i]);
                Dtilt[i] = 1000/(DLCPoT[i]-DCPoT[i]);
        }

        for(int k = 0 ; k < 24 ;k++)
        {   
                UpH_fit[k] = new TF1(Form("UpH_fit%d",k),"[0]*exp(-1*(min(x,[1])-[1])^2/(2*[2]^2))",-300,UCPoT[k]-20);//
                UpH_fit[k] -> SetParName(0,Form("UC0_%d",k));
                UpH_fit[k] -> SetParName(1,Form("UC1_%d",k));
                UpH_fit[k] -> SetParName(2,Form("UC2_%d",k));

                UpH_fit[k] -> SetParameter(0,1200);
                UpH_fit[k] -> SetParameter(1,UCPoT[k]);
                UpH_fit[k] -> SetParameter(2,Utilt[k]);



                DpH_fit[k] = new TF1(Form("DpH_fit%d",k),"[0]*exp(-1*([1]-max(x,[1]))^2/(2*[2]^2))",DCPoT[k]-100,300);
                DpH_fit[k] -> SetParName(0,Form("DC0_%d",k));
                DpH_fit[k] -> SetParName(1,Form("DC1_%d",k));
                DpH_fit[k] -> SetParName(2,Form("DC2_%d",k));

                DpH_fit[k] -> SetParameter(0,1200);
                DpH_fit[k] -> SetParameter(1,DCPoT[k]);
                DpH_fit[k] -> SetParameter(2,Dtilt[k]);


                // cout << Form("%d bar  ",k) << endl;

                int cannumber=k/1;
                can[cannumber] -> cd(k%1+1);
                // std::cout<<i<<"  number : "<<cannumber<<i%4 <<std::endl;
                TDC_diff[k] -> Fit(Form("UpH_fit%d",k),"","",TDC_diff[k]->GetMean() -200 , TDC_diff[k]->GetMean());
                TDC_diff[k] -> Fit(Form("DpH_fit%d",k),"","",TDC_diff[k]->GetMean(),TDC_diff[k]->GetMean()+200);//  
                TDC_diff[k] -> Draw();
                UpH_fit[k] -> Draw("same");
                //can[cannumber] -> SaveAs(Form("TDC_fitted%i.pdf",cannumber));

        ///////////////////////////////////////////////////////////
                Leftvalue[k] = UpH_fit[k] -> GetParameter(1);
                Rightvalue[k] = DpH_fit[k] -> GetParameter(1);
                Offset_val[k] = (Rightvalue[k] + Leftvalue[k])/2;
                
                offsetvalue[k] = Offset_val[k];
                treeOff -> Fill();
                
                cout << Form("%d bar Left position :  ",k) << Leftvalue[k] << endl;
                cout << Form("%d bar Right position :  ",k) << Rightvalue[k] << endl;
                cout << Form("%d offset value : ",k) << Offset_val[k] << endl;
                cout <<"   " << endl;

                //treeOff -> Fill();
        ///////////////////////////////////////////////////////////
        }
        treeOff -> Write();
        offset -> Close();
        
        /* for(int nb = 0 ; nb < NB ; nb++)
           {  
           array -> Fill(TDC_diff[BID[nb]],BID[nb]);
           }

           pan -> cd();
           pan -> SetLogz();
           array -> Draw("colz");
           */        

}
