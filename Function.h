#ifndef FUNCTION_H
#define FUNCTION_H

#include "TMath.h"

/*class FUNCTION 
{ 
public :
	
	double pedAvg(TH1D* hist);
	double ADC (TH1D* hist);
	double TDC(int range, TH1D* hist);
};*/

 double pedAvg(TH1D* hist)
  {
    double sum = 0;
    int n = 10;
    double avg = 0;
    for(int k = 0 ; k < n ; k++)
    {
      sum += hist -> GetBinContent(k+1);
    }
    avg = sum/n;
    return avg;
  }

  double ADC (TH1D* hist)
  {
    double hight = 0;
    double avg;

    avg = pedAvg(hist);    
    hight =  hist->GetMaximum() - avg;

    return hight;
  }

double QDC(int range,TH1D* hist)
  {
  
  double avg;
  avg = pedAvg(hist);
  
  double k = 0;
  double sumHigh = 0;  
  for(k = 0 ; k < range ; k++)
  {   
  sumHigh += abs(hist -> GetBinContent(k+1))-avg;
  }   
  return sumHigh;
  }

double TDC(int range, TH1D* hist)
  {
    double avg;
    avg = pedAvg(hist);    
    double thersH =  ADC(hist)/5 + avg ; 
    double position = 0;
    int n = 0 ; 

    for(n = 0; n < range; n++)
    {   

      if(thersH > hist -> GetBinContent(n+1) && thersH < hist -> GetBinContent(n+2))
      {
        position =2* (thersH - (hist -> GetBinContent(n+1) - (hist -> GetBinContent(n+2) - hist -> GetBinContent(n+1))*(n+1))) / ( hist -> GetBinContent(n+2) - hist -> GetBinContent(n+1) );
        break;
      }

      else if(thersH == hist -> GetBinContent(n+1))
      {
        position = 2* (n+1);
      }

    }
      return position;
  }


#endif
