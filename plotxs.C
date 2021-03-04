#include "math.h"
#include "TROOT.h"
#include "TF1.h"
#include "TCanvas.h"

double xs_vs_theta(double *x, double *par)
{
  const double deg=atan(1.0)/45.;
  double Theta=x[0]*deg;
  int Z=(int)(par[0]);
  int N=(int)(par[1]);
  double E0=par[2];
  double Ef=par[3];
  double xs=0.0;
  if(E0>Ef) xs = XEM::GetXS(Z,N,E0,Ef,Theta);
  //if(isnan(xs)) xs=0.0;
  
  return xs;
}

double xs_vs_p(double *x, double *par)
{
  const double deg=atan(1.0)/45.;
  double Ef=x[0];
  
  int Z=int(par[0]);
  int N=int(par[1]);
  double E0=par[2];
  double Theta=par[3]*deg;
  double xs=0.0;
  if(E0>Ef) xs = XEM::GetXS(Z,N,E0,Ef,Theta);
  //if(isnan(xs)) xs=0.0;
  
  return xs;
}


void plotxs1()
{
  TF1 *f1 = new TF1("f1",xs_vs_theta,5.,25.,4);
  f1->SetParameters(1,0,2.100,1.8);
  f1->Draw();
}

void plotxs()
{
  TCanvas *c1 = new TCanvas("c1","",800,600);
  
  TF1 *f1 = new TF1("f1",xs_vs_p,0.8,2.2,4);
  f1->SetParameters(6,6,2.100,11.7);
  f1->SetLineColor(1);
  f1->Draw();
    
  TF1 *f2 = new TF1("f2",xs_vs_p,0.8,2.2,4);
  f2->SetParameters(2,1,2.100,11.7);
  f2->SetLineColor(2);
  f2->Draw("same");
  
  TF1 *f3 = new TF1("f3",xs_vs_p,0.8,2.2,4);
  f3->SetParameters(1,0,2.100,11.7);
  f3->SetLineColor(3);
  f3->Draw("same");
}

