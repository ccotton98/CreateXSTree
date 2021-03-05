#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <math.h>
#include <time.h>

#include <TROOT.h>
#include <TTree.h>
#include <TFile.h>
#include <TMath.h>
#include <TVector3.h>
#include <TRandom.h>

#include "XSModel.hh"          //XSModel return DiffXS in ub/MeV-Sr for inelas or ub for elas.
#include "G2PRand.hh"
#include "HRSTransform_TCSNHCS.hh"
#include "ACCInc.h"

using namespace std;

//uncomment the next line if you want to compare these 2 function: the old one is binned in Theta, the new one is binnedin cosTheta 
//#define CMP_GETINTEXS

//DEBUG level will be used to control what message to print:
// >=1: will print pressure curve rates for H2, N2 and 3He target
// >=2: will print rates for each x bin for 3He target
// >=3: add VZ bin rates
//      if CMP_GETINTEXS is defined, will also print the comparison between old and new GetInteXS()
// >=4: print a lot of debug information in GetInteXS()     
// >=5: print more debug information in GetInteXS()      
// >=6: print even more debug information in GetInteXS() and GetXS()     
#define DEBUG 2

static bool  gSkipC12 = false;
static bool  gSkipPresureCurve = false;
static const double kC12FoilThick = 0.010*2.54;    //10 mil C12 foil thickness

static const double deg = acos(0.0)/90.0;
double GetXS(float Ei, float Ef, float Theta, float Tb, float Ta, int ElasOnly=0);
double GetXS_GE180(float Ei, float Ef, float Theta, float Tb, float Ta, int ElasOnly=0);

////////////////////////////////////////////////////////////////////////////
//return in unit of 10^30, corresponsed to microbarn, number of nuclei per second per cm2
double GetLumi10pow30(double current_na, int nuclei_per_molecule, double mol_mass_gpermol, double thickXdens_gpercm2)
{
  double I=current_na * 1.0e-9;
  double T=thickXdens_gpercm2 ;

  const double N_A = 6.022e+23; //Avagadro's number
  const double q_e = 1.602e-19; //electron charge

  double Lumi=0;
  Lumi = nuclei_per_molecule * ( I/q_e *N_A* T / mol_mass_gpermol ) /1.0e30;

  //L=3.759*1.0e33 * current_na*thickXdens_gpercm2/mol_mass_gpermol;
  if(DEBUG>=6) cout<<"Luminosity="<<Lumi<<" X 1.0E+30 #/s/cm2"<<endl;
  return Lumi;
}

////////////////////////////////////////////////////////////////////////////
//return elastic kinimatics for given beam energy and theta
double GetElasKin(double beam_gev, double theta_e_rad, double M_gev, 
                  double &Ef, double &ptot_p, double & theta_p_rad)
{
  double Ei=beam_gev, t=theta_e_rad, M=M_gev;
  Ef = Ei / (1+2*Ei/M*pow(sin(t/2),2.0));

  double Pperp_p,Pz_p;
  Pperp_p = -Ef * sin(t);
  Pz_p = Ei - Ef * cos(t);
  ptot_p = sqrt(Pperp_p*Pperp_p+Pz_p*Pz_p);
  theta_p_rad = asin(Pperp_p/ptot_p);

  if(DEBUG>=6) {
    cout<< "Ei="<<Ei<< " Theta_e="<<t/deg<<" M="<<M<<"  ==>  Ef=" <<Ef
        <<"  P_p=" <<ptot_p<<"   Theta_p=" <<theta_p_rad/deg<<endl;
  }

  return Ef; 
}

double GetElasEprime(double beam_gev, double theta_e_rad, double M_gev)
{
  double Ei=beam_gev, t=theta_e_rad, M=M_gev;
  double Ef = Ei / (1+2*Ei/M*pow(sin(t/2),2.0));
  return Ef; 
}

//xs wrapper
double GetXS(int Z, int N, float Ei, float Ef, float Theta, float Tb, float Ta, int ElasOnly)
{
#ifdef DEBUG 
  if(DEBUG>=6) cout<<" Executing GetXS(Z="<<Z<<", N="<<N<<", Ei="<<Ei<<", Ef="<<Ef<<", Theta="<<Theta<<", Tb="<<Tb<<", Ta="<<Ta<<")\n";
#endif

  if(Z<0) return GetXS_GE180(Ei,Ef,Theta,Tb,Ta,ElasOnly);
  
  double pXs=0.0;
  if(ElasOnly) {
    pXs = ElasModel::GetXS(Z, N, (double)Ei, (double)Theta);
    if(isnanf(pXs)) {
      pXs=0.0;
#ifdef DEBUG 
      if(DEBUG>=4) cout<<" isnan  from ElasModel::GetXS(Z="<<Z<<", N="<<N<<", Ei="<<Ei<<", Theta="<<Theta<<")\n";
#endif
    }
    pXs *= 1.0;   //already in ub/Sr
  } else {
    //note that XEM::GetXS() will not work for Q2>11, I give up these points
    pXs = XEM::GetXS(Z, N, (double)Ei, (double)Ef, (double)Theta, (double)Tb, (double)Ta);
    if(isnanf(pXs)) {
      pXs=0.0;
#ifdef DEBUG 
      if(DEBUG>=4) cout<<" isnan  from XEM::GetXS(Z="<<Z<<", N="<<N<<", Ei="<<Ei<<", Ef="<<Ef<<", Theta="<<Theta<<", Tb="<<Tb<<", Ta="<<Ta<<")\n";
      //char cc[100];cout<<"\nPress any key to continue ...";cin>>cc;
#endif
    }
    pXs *= 1000.0;   //turn from ub/MeV/Sr into ub/GeV/Sr
  }
  return pXs;
}


//XS for GE180 compound
double GetXS_GE180(float Ei, float Ef, float Theta, float Tb, float Ta, int ElasOnly)
{
  //http://galileo.phys.virginia.edu/research/groups/spinphysics/glass_properties.html
  //GE180 glass composition:
  //SiO2  60.3%
  //BaO   18.2%
  //Al2O3 14.3%
  //CaO   6.5%
  //SrO   0.25%
  // Z/A = 0.4829
  // The total of the above is 99.55%, what is the rest? I treat it as H
  //H    0.45%
  
  double MolMass_SiO2=(28.0856*1+15.9994*2);  //in g/mol   
  double MolMass_BaO =(137.327*1+15.9994*1);  //in g/mol
  double MolMass_Al2O3=(26.982*2+15.9994*3);  //in g/mol
  double MolMass_CaO =(40.078*1+15.9994*1);   //in g/mol
  double MolMass_SrO =(87.62*1+15.9994*1);    //in g/mol
  double MolMass_pr=1.00794;                  //in g/mol
    
  //const char*  kGE180_MName[]={"SiO2","BaO","Al2O3","CaO","SrO","H"};
  const double kGE180_MFraction[]={0.603,0.182,0.143,0.065,0.0025,0.0045};
  const double kGE180_MMolMass[]={MolMass_SiO2,MolMass_BaO,MolMass_Al2O3,MolMass_CaO,MolMass_SrO,MolMass_pr};
        
  double pMolMass_aver_ge180_rev = 0;
  for(int i=0;i<6;i++) {
    pMolMass_aver_ge180_rev += kGE180_MFraction[i]/kGE180_MMolMass[i];
  }
  double pMolMass_aver_ge180 = 1.0 / pMolMass_aver_ge180_rev;
  int pZ_aver_ge180 = lround(pMolMass_aver_ge180 * 0.4829);
  int pN_aver_ge180 = lround(pMolMass_aver_ge180-pZ_aver_ge180);  //round up to nearest int
    
  //comput xs for each element
  double pXs_O  = GetXS( 8,  8, Ei, Ef, Theta, Tb, Ta, ElasOnly);
  double pXs_Si = GetXS(14, 14, Ei, Ef, Theta, Tb, Ta, ElasOnly);
  double pXs_Ba = GetXS(56, 82, Ei, Ef, Theta, Tb, Ta, ElasOnly);
  double pXs_Al = GetXS(13, 14, Ei, Ef, Theta, Tb, Ta, ElasOnly);
  double pXs_Ca = GetXS(20, 20, Ei, Ef, Theta, Tb, Ta, ElasOnly);
  double pXs_Sr = GetXS(38, 50, Ei, Ef, Theta, Tb, Ta, ElasOnly);
  double pXs_H  = GetXS( 1,  0, Ei, Ef, Theta, Tb, Ta, ElasOnly);

  const double kGE180_MXs[] = {pXs_Si+2*pXs_O, pXs_Ba+pXs_O, 2*pXs_Al+3*pXs_O,pXs_Ca+pXs_O, pXs_Sr+pXs_O, pXs_H};
    
  double pXs_ge180=0.0;
  for(int i=0;i<6;i++) {
    pXs_ge180 += pMolMass_aver_ge180*kGE180_MFraction[i]/kGE180_MMolMass[i]*kGE180_MXs[i];
  }

#ifdef DEBUG 
  if(DEBUG>=5) {
    //comput xs using average Z and N
    double pXs=0.0;
    pXs = GetXS(pZ_aver_ge180, pN_aver_ge180, Ei, Ef, Theta, Tb, Ta, ElasOnly);
    
    cout<<" using average Z and N, XS_GE180 = "<<pXs
        <<"  sum up all elements, XS_GE180 = "<<pXs_ge180<<endl;
  }
#endif

  return pXs_ge180;
}


void GetA1NTaTb(double Z, double N, double VZ, double Theta, double& Tb, double& Ta)
{
  Ta=Tb=0.0;
  //uncomment the next line if you want to turn on rad. corr.
  //By Jixie: the xs model has bug in doing R.C., do not use this feature until this bug is fixed 
  return;
  
  double CosTh = cos(Theta);
  double SinTh = sin(Theta);
  //rad. corr., gas contribution is negligible, therefore I use 12 AMG of 3He for 3He at 10 ATM  
  const double kX0_GE180 = 7.040;   //cm
  const double kX0_3He = 41959.174; //cm, 3He at 12 amg, @ 10 ATM
  const double kX0_H2 = 75510.241;  //cm, @ 10 ATM
  const double kX0_N2 = 3274.254;   //cm, @ 10 ATM
  const double kX0_12C = 18.834;    //cm
  const double kCellWallThick = 0.10;    //1 mm
  const double kCellRadiusIn = 0.87*2.54/2.0-kCellWallThick;    //0.87 inch outer diameter

  //get the theta limit that a particle can hit the wall, assuming the cell is a perfect cylinder with 2 end-plate, not sphere shape end
  double Theta_W = atan2(kCellRadiusIn,20.0-VZ);
  if(Theta_W<0.0) Theta_W += 3.1415926535;
  
  //get the path length thrught which a paritle going out will go
  double PathLength_wall = kCellWallThick/SinTh / kX0_GE180; // in rad.len. unit
  double PathLength_end = 0.014/CosTh / kX0_GE180;  // in rad.len. unit
  double Pathlength_out = (Theta>Theta_W) ? PathLength_wall : PathLength_end;

  //add the thickness that outside the cell 
  double Pathlength_exbefore=0.0004;// (he4 bag and exit window of beam pipe)
  double Pathlength_extra = 0.01;//(he4 bag and entrance|exit window of the spectrometer)
  Pathlength_out += Pathlength_extra;
  
  double kWindow_thick=0.014;//GE180 window thickness
  double kTarget_length=40.0;//3He target cell length
  if(Z<0.0) {
    //GE180
    if(VZ<-15.0) {
      //upstream window: //0.007cm GE180
      Tb = 0.5*kWindow_thick/kX0_GE180+Pathlength_exbefore;
      if (Theta>Theta_W) {
        Ta = 0.5*kWindow_thick/CosTh/kX0_GE180 +kCellRadiusIn/SinTh/kX0_3He+ Pathlength_out;
      } else {
        Ta = 0.5*kWindow_thick/CosTh/kX0_GE180+(kTarget_length)/CosTh/kX0_3He + Pathlength_out;
      } 
    }
    if(VZ>15.0) {
      //downstream window: //0.014cm GE180 + 40cm 3He + 0.007cm GE180
      Tb = 1.5*kWindow_thick/kX0_GE180 + kTarget_length/kX0_3He+Pathlength_exbefore;
      Ta = 0.5*kWindow_thick/CosTh/kX0_GE180 + Pathlength_extra;
    }
  } else if(fabs(Z-1)<0.1 && fabs(N-0)<0.1) {
    //H2:   //0.014cm GE180 + (20+Z)cm H2 
    Tb = 1.0*kWindow_thick/kX0_GE180 + (0.5*kTarget_length+VZ)/kX0_H2+Pathlength_exbefore;
    if (Theta>Theta_W) {
      Ta = kCellRadiusIn/SinTh/kX0_H2 + Pathlength_out;
    } else {
      Ta = (0.5*kTarget_length-VZ)/CosTh/kX0_H2 + Pathlength_out;
    }
  } else if(fabs(Z-2)<0.1 && fabs(N-1)<0.1) {
    //He3:   //0.014cm GE180 + (20+Z)cm 3He 
    Tb =  1.0*kWindow_thick/kX0_GE180 + (0.5*kTarget_length+VZ)/kX0_3He+Pathlength_exbefore;
    if (Theta>Theta_W) {
      Ta = kCellRadiusIn/SinTh/kX0_3He + Pathlength_out;
    } else {
      Ta = (0.5*kTarget_length-VZ)/CosTh/kX0_3He + Pathlength_out;
    }
  } else if(fabs(Z-7)<0.1 && fabs(N-7)<0.1) {
    //N2:   //0.014cm GE180 + (20+Z)cm N2 
    Tb = 1.0*kWindow_thick/kX0_GE180 + (0.5*kTarget_length+VZ)/kX0_N2+Pathlength_exbefore;
    if (Theta>Theta_W) {
      Ta = kCellRadiusIn/SinTh/kX0_N2 + Pathlength_out;
    } else {
      Ta = (0.5*kTarget_length-VZ)/CosTh/kX0_N2 + Pathlength_out;
    }  
  } else if(fabs(Z-6)<0.1 && fabs(N-6)<0.1) {
    //C12:   //5 mil 12C
    Tb = 0.5*kC12FoilThick/kX0_12C+Pathlength_exbefore;
    Ta = 0.5*kC12FoilThick/CosTh/kX0_12C+Pathlength_extra;  
  } else {
    Tb = 0.0; 
    Ta = 0.0;
  }

  if(DEBUG>=6) {
    cout<<" Z="<<Z<<" N="<<N<<" VZ="<<VZ<<" Theta_W="<<Theta_W/deg<<" Theta="<<Theta/deg<<" (deg) ==> Ta="<<Ta<<",  Tb="<<Tb<<endl;
  }
  return;
}
  
////////////////////////////////////////////////////////////////////////////
//old version is binned in Theta-Phi-P, new version is binned in CosTheta-Phi-P
//get integrated xs for given element at given vertex z,
//applying x, q2 and w cuts if their upper limits are larger than zero
//ElasOnly=-1: pure inelastic for full acceptance
//ElasOnly=0:  inelastic + elastic for full acceptance
//ElasOnly=1:  pure elastic for full acceptance
//ElasOnly=-30: pure inelastic for 2-SC-Bar acceptance
//ElasOnly=30: inelastic + elastic for 2-SC-Bar acceptance
//ElasOnly=31: pure elastic for 2-SC-Bar acceptance
double GetInteXS(double pBeamE, double pAngle, double pMomentum, double Z, int N, double VZ, string Name, int ElasOnly=0,
                 double pXmin=-1.0, double pXmax=-1.0, double pQ2min=-1.0, double pQ2max=-1.0,double pWmin=-1.0, double pWmax=-1.0)
{
  int pre = cout.precision();   //get the original precision of cout
  cout.precision(4);
  //for HMS
  int run = 1;
  double pTheta_tg_max=0.058;   //to better match method2
  double pTheta_tg_min=-0.054;
  double pPhi_tg_max=0.031;    
  double pPhi_tg_min=-0.031; 
  double pEprime_max=pMomentum*1.13;
  double pEprime_min=pMomentum*0.89;
  if(Name=="SHMS") {
    run = 2;
    pTheta_tg_max=0.046;
    pTheta_tg_min=-0.046;    //to better match method2
    pPhi_tg_max=0.025;
    pPhi_tg_min=-0.028;
    pEprime_max=pMomentum*1.27;
    pEprime_min=pMomentum*0.82;
  }
  if(pEprime_max>pBeamE) pEprime_max=pBeamE;

  double Tb=0.0,Ta=0.0;
  
  const double AMU = 0.9314941;
  double pMtg=(fabs(Z)+fabs(N))*AMU;
  double pP_elas=0.0;

  double pTheta_tg,pPhi_tg=0;
  double pTheta=-999.,pCosTh=-999.,pPhi=-999.,pEprime=-999.;
  double pXs_elas=0.0,pXs=0.0,pInteXs=0.0;

  double pTheta_min = pAngle - 20*deg;
  double pTheta_max = pAngle + 20*deg;
  if(pTheta_min < 5*deg) pTheta_min = 5*deg;
  double pCosTh_min = cos(pTheta_max);
  double pCosTh_max = cos(pTheta_min);  
  
  double pPhi_min = -80 *deg ;
  double pPhi_max =  80 *deg ;

  //double dCosTh = 0.0005;  //xs changes rapidly in theta, therefore put tiny bin width here
  //double dPhi = 0.002;
  //double dEprime = 0.002;
  //If Change dPhi from 0.002 to 0.01, geometry acceptance change is <0.3% for all targets
  //If change dCosTh from 0.0005 to 0.001, 3He geometry acceptance will change 4%
  //but geometry acceptance  for GE180(glass) at vz=-20cm will change ~7%
  //If change dEprime from 0.002 to 0.005, geometry acceptance change is <0.3% for all targets
  //except geometry acceptance for GE180(glass) at vz=-20cm will change ~7%
  double dCosTh = 0.0005;  //xs changes rapidly in theta, therefore put tiny bin width here
  double dPhi = 0.01;
  double dEprime = 0.002;
  
  double x_tg,y_tg,z_tg=0;

  pCosTh = pCosTh_min - 0.5*dCosTh;
  while (pCosTh < pCosTh_max) {
    
    pCosTh += dCosTh;
    pTheta = acos(pCosTh);
    //get elastic Eprime 
    pP_elas = GetElasEprime(pBeamE, pTheta, pMtg);
    double dOmega = dCosTh*dPhi;

    //reset pEprime_max for each theta since pP_elas is the maximum momentum of a physics event
    double thisEprime_max = pEprime_max;
    if(thisEprime_max>pP_elas) thisEprime_max=pP_elas;
    if(DEBUG>=4) {
      cout<<" pTheta="<<pTheta/deg<<"(deg),  pEprime_min="<<pEprime_min<<"  pEprime_max="<<pEprime_max<<"  thisEprime_max="<<thisEprime_max<<endl;
    }

    //get the radation length of the path length
    GetA1NTaTb(Z,N,VZ,pTheta,Tb,Ta);
  
    pPhi = pPhi_min - 0.5*dPhi;
    while (pPhi < pPhi_max) {
      pPhi += dPhi;

      //convert Lab frame to Tranportation frame
      //void  P_HCS2TCS(double Theta_hall, double Phi_hall, double EndPlaneTheta_hall, double &Theta_tr, double &Phi_tr)
      Transform::P_HCS2TCS(pTheta, pPhi, pAngle, pTheta_tg, pPhi_tg);
      //Jixie: do not need this cut any more, but I found that method 1 is 10% more than method 2
      //I am using these 2 cuts to make them match ( I could be wrong, since method 2 might underestimate rates by itself ...)
      //conclusion:
      //1) Theta_tg cut will reduce up to 15% geometry acceptance for HMS, up to 7% geometry acceptance for SHMS
      //2) Phi_tg cut will reduce ~3.5% geometry acceptance.
      //3) I would like to cut only onto Theta_tg, since Phi_tg has strong VZ dependence.  If I cut Phi_tg, GE180 calculation
      //   will not be reliable.
      if(pTheta_tg>pTheta_tg_max || pTheta_tg<pTheta_tg_min)  continue;
      //if(pPhi_tg>pPhi_tg_max || pPhi_tg<pPhi_tg_min)  continue;  
      
      //project to target plane     
      //void Project(double x,double y,double z,double z_drift,double theta,double phi,double &x_out, double &y_out, double &z_out)
      x_tg=y_tg=z_tg=0;
      
      //void  X_HCS2TCS(double x, double y, double z,double EndPlaneTheta_hall, double &x_tr, double &y_tr, double &z_tr)
      double x_tr=0, y_tr=0, z_tr=0;
      Transform::X_HCS2TCS(0,0,VZ,pAngle,x_tr,y_tr,z_tr);
      Transform::Project(x_tr, y_tr, z_tr, -z_tr,pTheta_tg,pPhi_tg,x_tg,y_tg,z_tg);
      if(DEBUG>=4) {
        cout<<" x_tr="<<x_tr<<"  y_tr="<<y_tr<<"  z_tr="<<z_tr<<"  theta_tr="<<pTheta_tg<<"  phi_tr="<<pPhi_tg
            <<" ==> x_tg="<<x_tg<<"  y_tg="<<y_tg<<"  z_tg="<<z_tg<<endl;
      }

      pEprime = pEprime_min - 0.5*dEprime;
      double pEprime_up = pEprime + 0.5*dEprime;
      while (pEprime < thisEprime_max) {
        //check if this is the last bin
        double deltaEprime = 0.0;
        double pEprimeLeft = pEprime_max - pEprime_up;   //how much Eprime has not been integrated yet
        if(pEprimeLeft <= 0.0) break;
        if(pEprimeLeft > dEprime) {
          pEprime += dEprime;
          pEprime_up += dEprime;
          deltaEprime = dEprime;
        } else {
          pEprime = pEprime_up + pEprimeLeft/2.0;
          pEprime_up += pEprimeLeft;
          deltaEprime = pEprimeLeft;
        }
        
        ////////////////////////////////////////////////////////////////////////
        //now apply acc cut
        //double ACCEPTANCE::GetAcceptance(int run,double pYtar,double pDelta,double pTheta, double pPhi,int type);
        int pType = (ElasOnly==30 || ElasOnly==31)?2:1;   //if (ElasOnly==3), only turn on 2 SC bars
        double pDelta = (pEprime - pMomentum) / pMomentum * 100.;
        double pAcc = ACCEPTANCE::GetAcceptance(run,y_tg,pDelta,pTheta_tg,pPhi_tg,pType);
        if(DEBUG>=4) cout<<" y_tg="<<y_tg<<" theta_tg="<<pTheta_tg<<"  phi_tg="<<pPhi_tg<<"  delta="<<pDelta<<"  ==> pAcc="<<pAcc<<endl;
        if(pAcc<0.1) continue;
        
        ////////////////////////////////////////////////////////////////////////
        double pQ2=2.0*pBeamE*pEprime*(1.0-cos(pTheta));
        double pXbj=pQ2/(2.0*0.9383*(pBeamE-pEprime));
        double M_p=0.9383;
        double nu=pBeamE-pEprime;
        double pW=sqrt(pow(M_p,2.0)+2*M_p*nu-pQ2);
        //applying X cuts
        if(pXmax>0.0 && pXmax>pXmin) {
          if(pXbj<pXmin || pXbj>pXmax) continue;
        }
        //applying Q2 cuts
        if(pQ2max>0.0 && pQ2max>pQ2min) {
          if(pQ2<pQ2min || pQ2>pQ2max) continue;
        }
        //applying W cuts: for Delta transverse Asymmetry, need to apply W cut for 3He Delta resonance peak        
        if(pWmax>0.0 && pWmax>pWmin) {
          if(pW<pWmin || pW>pWmax) continue;//continue for jump the rest of the commands
        }
        
        if(DEBUG>=4) cout<<" pTheta="<<pTheta/deg<<"  pPhi="<<pPhi/deg<<"  pEprime="<<pEprime<<"  deltaEprime="<<deltaEprime<<"  pP_elas="<<pP_elas<<endl;
        
        if(ElasOnly!=1 && ElasOnly!=31) {
          pXs = 0.0;
          if(pQ2 < 11.0) pXs = GetXS(Z, N, pBeamE, pEprime, pTheta, Tb, Ta, 0);
          if(DEBUG>=4) cout<<" Xbj = "<< pXbj<<"  Q2 = "<<pQ2<<"  XEM::GetXS() = "<<pXs*1.0E3<<" (nb/GeV/Sr)"<<endl;
          pInteXs += pXs*deltaEprime*dOmega*pAcc;
        }
        
        //check if elastic events are accepted or not, add only once per dOmega
        //please note that elas xs is for per nucleus
        if(fabs(pP_elas-pEprime)<0.51*deltaEprime && ElasOnly>=0) {
          pXs_elas = GetXS(Z, N, pBeamE, pEprime, pTheta, Tb, Ta, 1);
          if(DEBUG>=4) cout<<" pP_elas="<<pP_elas<<" deltaEprime="<<deltaEprime<<",  ElasModel::GetXS() = "<<pXs_elas*1.0E3<<" (nb/Sr)"<<endl;
          pInteXs += pXs_elas*dOmega*pAcc;
          if(DEBUG>=7) {char cc[100];cout<<"\nPress any key to continue ...";cin>>cc;}
        } 
        
      }
    }    
  }
  
  cout.precision(pre);    //recover cout default precision
  return pInteXs;
}

////////////////////////////////////////////////////////////////////////////
//old version is binned in Theta-Phi-P, new version is binned in CosTheta-Phi-P
//get integrated xs for given element at given vertex z,
//applying x and q2 cuts if their upper limits are larger than zero
//ElasOnly=-1: pure inelastic for full acceptance
//ElasOnly=0:  inelastic + elastic for full acceptance
//ElasOnly=1:  pure elastic for full acceptance
//ElasOnly=-30: pure inelastic for 2-SC-Bar acceptance
//ElasOnly=30: inelastic + elastic for 2-SC-Bar acceptance
//ElasOnly=31: pure elastic for 2-SC-Bar acceptance
double GetInteXS_old(double pBeamE, double pAngle, double pMomentum, double Z, int N, double VZ, string Name, int ElasOnly=0,
                     double pXmin=-1.0, double pXmax=-1.0, double pQ2min=-1.0, double pQ2max=-1.0,double pWmin=-1.0, double pWmax=-1.0)
{
  int pre = cout.precision();   //get the original precision of cout
  cout.precision(4);
  //for HMS
  int run = 1;
  double pTheta_tg_max=0.058;    //to better match method2
  double pTheta_tg_min=-0.054;
  double pPhi_tg_max=0.031;
  double pPhi_tg_min=-0.031; 
  double pEprime_max=pMomentum*1.13;
  double pEprime_min=pMomentum*0.89;
  if(Name=="SHMS") {
    run = 2;
    pTheta_tg_max=0.046;
    pTheta_tg_min=-0.046;
    pPhi_tg_max=0.025;
    pPhi_tg_min=-0.028;
    pEprime_max=pMomentum*1.27;
    pEprime_min=pMomentum*0.82;
  }
  if(pEprime_max>pBeamE) pEprime_max=pBeamE;
    
  double Tb=0.0,Ta=0.0;
  
  const double AMU = 0.9314941;
  double pMtg=(fabs(Z)+fabs(N))*AMU;
  double pP_elas=0.0;

  double pTheta_tg,pPhi_tg=0;
  double pTheta=-999.,pPhi=-999.,pEprime=-999.;
  double pXs_elas=0.0,pXs=0.0,pInteXs=0.0;

  double pTheta_min = pAngle - 20*deg;
  double pTheta_max = pAngle + 20*deg;
  if(pTheta_min < 5*deg) pTheta_min = 5*deg;
  
  double pPhi_min = -80 *deg ;
  double pPhi_max =  80 *deg ;

  double dTheta = 0.2*deg;  //xs changes rapidly in theta, therefore put tiny bin width here
  double dPhi = 0.5*deg;
  double dEprime = 0.005;
  //double dTheta = 0.5*deg;
  //double dPhi = 1.0*deg;
  //double dEprime = 0.01;
  
  double x_tg,y_tg,z_tg=0;

  pTheta = pTheta_min - 0.5*dTheta;
  while (pTheta < pTheta_max) {
    
    pTheta += dTheta;
    //get elastic Eprime 
    pP_elas = GetElasEprime(pBeamE, pTheta, pMtg);
    double dOmega = (cos(pTheta-0.5*dTheta)-cos(pTheta+0.5*dTheta))*dPhi;

    //reset pEprime_max for each theta since pP_elas is the maximum momentum of a physics event
    double thisEprime_max = pEprime_max;
    if(thisEprime_max>pP_elas) thisEprime_max=pP_elas;
    if(DEBUG>=4) {
      cout<<" pTheta="<<pTheta/deg<<"(deg),  pEprime_min="<<pEprime_min<<"  pEprime_max="<<pEprime_max<<"  thisEprime_max="<<thisEprime_max<<endl;
    }
    
    //get the radation length of the path length
    GetA1NTaTb(Z,N,VZ,pTheta,Tb,Ta);

    pPhi = pPhi_min - 0.5*dPhi;
    while (pPhi < pPhi_max) {
      pPhi += dPhi;

      //convert Lab frame to Tranportation frame
      //void  P_HCS2TCS(double Theta_hall, double Phi_hall, double EndPlaneTheta_hall, double &Theta_tr, double &Phi_tr)
      Transform::P_HCS2TCS(pTheta, pPhi, pAngle, pTheta_tg, pPhi_tg);
      //Jixie: do not need this cut any more
      if(pTheta_tg>pTheta_tg_max || pTheta_tg<pTheta_tg_min)  continue;
      //if(pPhi_tg>pPhi_tg_max || pPhi_tg<pPhi_tg_min)  continue;
      
      //project to target plane     
      //void Project(double x,double y,double z,double z_drift,double theta,double phi,double &x_out, double &y_out, double &z_out)
      x_tg=y_tg=z_tg=0;
      
      //void  X_HCS2TCS(double x, double y, double z,double EndPlaneTheta_hall, double &x_tr, double &y_tr, double &z_tr)
      double x_tr=0, y_tr=0, z_tr=0;
      Transform::X_HCS2TCS(0,0,VZ,pAngle,x_tr,y_tr,z_tr);
      Transform::Project(x_tr, y_tr, z_tr, -z_tr,pTheta_tg,pPhi_tg,x_tg,y_tg,z_tg);
      if(DEBUG>=4) {
        cout<<" x_tr="<<x_tr<<"  y_tr="<<y_tr<<"  z_tr="<<z_tr<<"  theta_tr="<<pTheta_tg<<"  phi_tr="<<pPhi_tg
            <<" ==> x_tg="<<x_tg<<"  y_tg="<<y_tg<<"  z_tg="<<z_tg<<endl;
      }

      pEprime = pEprime_min - 0.5*dEprime;
      double pEprime_up = pEprime + 0.5*dEprime;
      while (pEprime < thisEprime_max) {
        //check if this is the last bin
        double deltaEprime = 0.0;
        double pEprimeLeft = pEprime_max - pEprime_up;   //how much Eprime has not been integrated yet
        if(pEprimeLeft <= 0.0) break;
        if(pEprimeLeft > dEprime) {
          pEprime += dEprime;
          pEprime_up += dEprime;
          deltaEprime = dEprime;
        } else {
          pEprime = pEprime_up + pEprimeLeft/2.0;
          pEprime_up += pEprimeLeft;
          deltaEprime = pEprimeLeft;
        }
        
        ////////////////////////////////////////////////////////////////////////
        //now apply acc cut
        //double ACCEPTANCE::GetAcceptance(int run,double pYtar,double pDelta,double pTheta, double pPhi,int type);
        int pType = (ElasOnly==30 || ElasOnly==31)?2:1;   //if (ElasOnly==3), only turn on 2 SC bars
        double pDelta = (pEprime - pMomentum) / pMomentum * 100.;
        double pAcc = ACCEPTANCE::GetAcceptance(run,y_tg,pDelta,pTheta_tg,pPhi_tg,pType);
        if(DEBUG>=4) cout<<" y_tg="<<y_tg<<" theta_tg="<<pTheta_tg<<"  phi_tg="<<pPhi_tg<<"  delta="<<pDelta<<"  ==> pAcc="<<pAcc<<endl;
        if(pAcc<0.1) continue;
        
        ////////////////////////////////////////////////////////////////////////
        double pQ2=2.0*pBeamE*pEprime*(1.0-cos(pTheta));
        double pXbj=pQ2/(2.0*0.9383*(pBeamE-pEprime));
        double M_p=0.9383;
        double nu=pBeamE-pEprime;
        double pW=sqrt(pow(M_p,2.0)+2*M_p*nu-pQ2);
        //applying X cuts
        if(pXmax>0.0 && pXmax>pXmin) {
          if(pXbj<pXmin || pXbj>pXmax) continue;
        }
        //applying Q2 cuts
        if(pQ2max>0.0 && pQ2max>pQ2min) {
          if(pQ2<pQ2min || pQ2>pQ2max) continue;
        }
        //applying W cuts: for Delta transverse Asymmetry, need to apply W cut for 3He Delta resonance peak        
        if(pWmax>0.0 && pWmax>pWmin) {
          if(pW<pWmin || pW>pWmax) continue;//continue for jump the rest of the commands
        }
        
        if(DEBUG>=4) cout<<" pTheta="<<pTheta/deg<<"  pPhi="<<pPhi/deg<<"  pEprime="<<pEprime<<"  deltaEprime="<<deltaEprime<<"  pP_elas="<<pP_elas<<endl;
        
        if(ElasOnly!=1 && ElasOnly!=31) {
          pXs = 0.0;
          if(pQ2 < 11.0) pXs = GetXS(Z, N, pBeamE, pEprime, pTheta, Tb, Ta, 0);
          if(DEBUG>=4) cout<<" Xbj = "<< pXbj<<"  Q2 = "<<pQ2<<"  XEM::GetXS() = "<<pXs*1.0E3<<" (nb/GeV/Sr)"<<endl;
          pInteXs += pXs*deltaEprime*dOmega*pAcc;
        }
        
        //check if elastic events are accepted or not, add only once per dOmega
        //please note that elas xs is for per nucleus
        if(fabs(pP_elas-pEprime)<0.51*deltaEprime && ElasOnly>=0) {
          pXs_elas = GetXS(Z, N, pBeamE, pEprime, pTheta, Tb, Ta, 1);
          if(DEBUG>=4) cout<<" pP_elas="<<pP_elas<<" deltaEprime="<<deltaEprime<<",  ElasModel::GetXS() = "<<pXs_elas*1.0E3<<" (nb/Sr)"<<endl;
          pInteXs += pXs_elas*dOmega*pAcc;
          if(DEBUG>=7) {char cc[100];cout<<"\nPress any key to continue ...";cin>>cc;}
        } 
        
      }
    }    
  }
  
  cout.precision(pre);    //recover cout default precision
  return pInteXs;
}

////////////////////////////////////////////////////////////////////////////
//Beam Current in uA, All energies are in GeV unit. All angles are in radian unit.
//return array contains rates for the following 7 materials: {"C12","He3","GE-180","GE-180","H2","N2","He3"};
double* GetRate(double pBeamCurrent, double pBeamE, double pDetectorAngle, double pDetectorMomentum, string pDetectorName, int pElasOnly=0)
{
  if(pElasOnly==-1)  cout<<"\n===================Full acceptance:  Pure Inelastic =======================\n";
  if(pElasOnly==0)   cout<<"\n===================Full acceptance:  Inelastic + Elastic ==================\n";
  if(pElasOnly==1)   cout<<"\n===================Full acceptance:  Pure elastic =========================\n";
  if(pElasOnly==2)   cout<<"\n===================Full acceptance:  With Delta (1.1<W<1.35) Cut ==========\n";
  if(pElasOnly==4)   cout<<"\n===================Full acceptance:  With DIS (2.0<W<100.0) Cut ===========\n";
  if(pElasOnly==5)   cout<<"\n===================Full acceptance:  With Resonance (0.0<W<2.0) Cut =======\n";
  if(pElasOnly==-30) cout<<"\n===================Only 2 SC bars:  Pure Inelastic ========================\n";
  if(pElasOnly==30)  cout<<"\n===================Only 2 SC bars:  Inelastic + Elastic ===================\n";
  if(pElasOnly==31)  cout<<"\n===================Only 2 SC bars:  Pure elastic ==========================\n";
  int pre = cout.precision();   //get the original precision of cout
  //define material 

  //const double amu = 0.93149403;//molar mass constant in GeV/c2

  double MolMass_pr=1.00794;    //in g/mol
  double MolMass_he3=3.01603;   //in g/mol
  double MolMass_c12=12.01078;  //in g/mol
  double MolMass_n14=14.0067;   //in g/mol
  double MolMass_ge180 = 54.7251;
  
  double Z_ge180 = MolMass_ge180 * 0.4829;
  int N_ge180 = lround(MolMass_ge180-Z_ge180);  //round up to nearest int
  Z_ge180 = MolMass_ge180 - N_ge180;
  
  if(DEBUG>=3) cout<<"  Z_ge180="<<Z_ge180<<"  N_ge180="<<N_ge180<<"    MolMass_ge180="<<MolMass_ge180<<endl;
  
  //double DH2=8.349E-4; // g/cm3; H2 gas density at 21.1°C (70°F) 10atm; this is correct
  //double DN2=0.0116;   // g/cm3; N2 gas density at 21.1°C (70°F) 10atm;
  double Dc12=2.267;     // g/cm3;
  double Dhe3=1.6147E-3; // g/cm3;  //Density at 12 amg, 
  double Dge180 = 2.77;  // g/cm3;

  //comput density using ideal gas law
  double R_const=82.057; //calculate gas density from ideal gas law, R=82.057 cm^3*atm*k^(-1)*mol^(-1)
  double Temp=21.1+273.15;//Density at 21.1°C   
  
  double DH2=10*MolMass_pr*2.0/(Temp*R_const);       // g/cm3; H2 gas density at 21.1°C (70°F) 10atm;
  double DN2=10*MolMass_n14*2.0/(Temp*R_const);      // g/cm3; N2 gas density at 21.1°C (70°F) 10atm; 
  double D3He=10*MolMass_he3/(Temp*R_const);         // g/cm3; 3He gas density at 21.1°C (70°F) 10atm;
  
  if(DEBUG>=3) cout<<"  density_he3_@_12amg="<<Dhe3<<" g/cm^3"<<", (according to ideal gas law, density_@_12_atm_294.25k="<<D3He*1.2<<")"<<endl;
  
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
   
  double pI_na=pBeamCurrent*1000.0;
  double pThickXDens;// in g/cm2
  double pLumi=0, pInteXs, pInteRate;

  cout <<"\n Detector= "<<pDetectorName<<",  Beam_current= "<<pBeamCurrent<<" uA,  Angle= "<<pDetectorAngle/deg<<" deg"<<endl;

  //this array will be returned to the caller, it has to be static, otherwise this array will be erased once this function finish 
  static double myRates[7];  //use this to keep rates for each target material in L1 array
  for(int i=0;i<7;i++) {myRates[i]=0.0;}
   
  //Loop1, 
  const char *L1Name[]={"C12","He3","GE-180","GE-180","H2","N2","He3"};
  double L1Z[]={6., 2., double(-Z_ge180), double(-Z_ge180), 1.0, 7.0, 2.0};  //number of protons in one atom, for GE180, use negative z
  double L1N[]={6., 1., double( N_ge180), double( N_ge180), 0.0, 7.0, 1.0};  //number of neutrons in one atom
  int L1Nuclei_per_molecule[]={1, 1, 1, 1, 2, 2, 1};  //number of atoms per molecule, use 1 for all compounds, n for Hn, Dn
  double L1Dens[]={Dc12, Dhe3, Dge180, Dge180, DH2, DN2, D3He};  //g/cm3
  double L1MolMass[]={MolMass_c12, MolMass_he3, MolMass_ge180, MolMass_ge180, MolMass_pr*2, MolMass_n14*2, MolMass_he3};  //g
  double L1Thick[]={kC12FoilThick, 40.0, 0.014, 0.014, 40.0, 40.0, 40.0}; //cm
  double L1VZ[]={0.0, 0.0, -20.0, 20.0, 0.0, 0.0, 0.0};  // vertex location in cm
  
  //x an Q2 binning
  //v/cr vxcen(15) R 0.250 0.300 0.350 0.400 0.450 0.500 0.550 0.600 0.650 0.705 0.765 0.825 0.885 0.945 1.005
  //v/cr vxbin(15) R 0.025 0.025 0.025 0.025 0.025 0.025 0.025 0.025 0.025 0.030 0.030 0.030 0.030 0.030 0.030
  double L1Xbj[]={0.025, 0.225, 0.275, 0.325, 0.375, 0.425, 0.475, 0.525, 0.575, 0.625, 0.675, 0.735, 0.795, 0.855, 0.915, 0.975};  // x boundary
  double L1Q2[]={2.0, 2.5, 3.0, 3.5, 4.0, 5.0, 6.0, 7.0};  //Q2 boundary


  cout.setf(ios_base::fixed); 
  cout.precision(4);
  cout<<setw(4)<<"Det"<<setw(8)<<"Target"
      <<setw(10)<<"BeamE"<<setw(10)<<"Current"<<setw(10)<<"DetP0"
      <<setw(10)<<"DetAngle"<<setw(10)<<"VZ(cm)"<<setw(10)<<" Thick(cm)"
      <<setw(12)<<"Lumi(10^33)"<<setw(12)<<"InteXS(pb)"<<setw(12)<<"Rate(Hz)"
      <<endl;
  //just for compare to method2, in this case: pElasOnly<0, inelastic only, only calculate rates for i>=4
  int istart = (pElasOnly<0)?4:0;
  for(int i=istart;i<7;i++) {
    cout.precision(3);
    pThickXDens=L1Dens[i]*L1Thick[i]; // (g/cm2);
    pLumi=GetLumi10pow30(pI_na,L1Nuclei_per_molecule[i],L1MolMass[i],pThickXDens);
    
    pInteXs = 0.0;

    //I am use these 2 lines to control if to run C12 and Pressure curve targets
    if(gSkipC12 && i==0) continue;
    if(gSkipPresureCurve && i>=4 && i<=6) continue;
    
    //do C12 or pressure curve only for 1-pass elastic kinematic points, Perform DIS runs for reference cell 3He, N2, Hz at 10 atm.
    //this part is now controlled in A1NRates()
    //if((i==0)  && !(fabs(pBeamE-2.1)<0.1 && pDetectorMomentum>2.0)) continue;
    
    //I will also do x binning for He3 DIS run
    if(i==1 && pBeamE>8.0) {
      //He3, calculate xs for each x-z bin
      for(int ix=0;ix<15;ix++) {
      
        int nTry = 20;
        double pXS_he3 = 0, pXS_tot_he3 = 0;
        for(int iz=0;iz<nTry;iz++) {
          double pVZ = -0.5*L1Thick[i] + (iz+0.5)* (L1Thick[i]/nTry);
          
          if(pElasOnly==2)
            pXS_he3=GetInteXS(pBeamE, pDetectorAngle, pDetectorMomentum, L1Z[i], L1N[i], pVZ, pDetectorName, pElasOnly,L1Xbj[ix],L1Xbj[ix+1],-1.,-1.,1.10,1.35);
          else if(pElasOnly==4)
            pXS_he3=GetInteXS(pBeamE, pDetectorAngle, pDetectorMomentum, L1Z[i], L1N[i], pVZ, pDetectorName, pElasOnly,L1Xbj[ix],L1Xbj[ix+1],-1.,-1.,2.00,100.0);
          else if(pElasOnly==5)
            pXS_he3=GetInteXS(pBeamE, pDetectorAngle, pDetectorMomentum, L1Z[i], L1N[i], pVZ, pDetectorName, pElasOnly,L1Xbj[ix],L1Xbj[ix+1],-1.,-1.,0.00,2.0);          
          else
            pXS_he3=GetInteXS(pBeamE, pDetectorAngle, pDetectorMomentum, L1Z[i], L1N[i], pVZ, pDetectorName, pElasOnly, L1Xbj[ix],L1Xbj[ix+1]);

          if(DEBUG>=3) cout<<"  x = "<<(L1Xbj[ix]+L1Xbj[ix+1])/2.0<<" +/- "<<(L1Xbj[ix+1]-L1Xbj[ix])/2.0<<"  VZ = "<<setw(6)<<pVZ<<"  InteXS_he3 = "<<setw(8)<<pXS_he3*1.0E6<<" (pb)"<<endl;
          pXS_tot_he3 += pXS_he3;
        }
        if(DEBUG>=3) cout<<"  Averaged  InteXS_he3 = "<<pXS_tot_he3/nTry*1.0E6<<" (pb)"<<endl;
      
        double pInteXs_xq=pXS_tot_he3/nTry;  //Need to get average, or change the luminosity during integration
        double pRate_xq=pLumi*pInteXs_xq;
        if(DEBUG>=2) cout<<"  x="<<(L1Xbj[ix]+L1Xbj[ix+1])/2.0<<" +/- "<<(L1Xbj[ix+1]-L1Xbj[ix])/2.0<<"  pInteXs_xq = "<<setw(8)<<pInteXs_xq*1.0E6<<" (pb)  Rate="<<pRate_xq<<" (Hz)"<<endl;

        pInteXs += pInteXs_xq;  
      }
    }
    else if (i==1 || (i>=4 && i<=6)) {
      //long target, calculate xs for each z bin     
      int nTry = 20;
      double pXS_long = 0, pXS_tot_long = 0;
#ifdef CMP_GETINTEXS        
      double pXS_long_old, pXS_tot_long_old=0;
#endif      
      //for He3 target, calculate rate for each x bin
      for(int iz=0;iz<nTry;iz++) {
        double pVZ = -0.5*L1Thick[i] + (iz+0.5)* (L1Thick[i]/nTry);
        if(pElasOnly==2)
          pXS_long=GetInteXS(pBeamE, pDetectorAngle, pDetectorMomentum, L1Z[i], L1N[i], pVZ, pDetectorName, pElasOnly,-1.,-1.,-1.,-1.,1.10,1.35);
        else if(pElasOnly==4)
          pXS_long=GetInteXS(pBeamE, pDetectorAngle, pDetectorMomentum, L1Z[i], L1N[i], pVZ, pDetectorName, pElasOnly,-1.,-1.,-1.,-1.,2.00,100.0);
        else if(pElasOnly==5)
          pXS_long=GetInteXS(pBeamE, pDetectorAngle, pDetectorMomentum, L1Z[i], L1N[i], pVZ, pDetectorName, pElasOnly,-1.,-1.,-1.,-1.,0.00,2.0);
        else
          pXS_long=GetInteXS(pBeamE, pDetectorAngle, pDetectorMomentum, L1Z[i], L1N[i], pVZ, pDetectorName, pElasOnly);
        if(DEBUG>=3) cout<<"  VZ = "<<setw(6)<<pVZ<<"  InteXS_"<<L1Name[i]<<" = "<<setw(8)<<pXS_long*1.0E6<<" (pb)"<<endl;
        pXS_tot_long += pXS_long;

        //,,,,,,,,,,,,,,,,compare old and new GetInteXS() start ,,,,,,,,,,,,,,,,,
#ifdef CMP_GETINTEXS        
        if(DEBUG>=3) {
          if(pElasOnly==2)
            pXS_long_old=GetInteXS_old(pBeamE, pDetectorAngle, pDetectorMomentum, L1Z[i], L1N[i], pVZ, pDetectorName, pElasOnly,-1.,-1.,-1.,-1.,1.10,1.35);
          else if(pElasOnly==4)
            pXS_long_old=GetInteXS_old(pBeamE, pDetectorAngle, pDetectorMomentum, L1Z[i], L1N[i], pVZ, pDetectorName, pElasOnly,-1.,-1.,-1.,-1.,2.00,100.0);
          else if(pElasOnly==5)
            pXS_long_old=GetInteXS_old(pBeamE, pDetectorAngle, pDetectorMomentum, L1Z[i], L1N[i], pVZ, pDetectorName, pElasOnly,-1.,-1.,-1.,-1.,0.00,2.0);
          else
            pXS_long_old=GetInteXS_old(pBeamE, pDetectorAngle, pDetectorMomentum, L1Z[i], L1N[i], pVZ, pDetectorName, pElasOnly);
          pXS_tot_long_old += pXS_long_old;
          cout<<" New - OLD:  VZ = "<<setw(6)<<pVZ<<"  InteXS_"<<L1Name[i]<<" = "<<setw(8)<<(pXS_long-pXS_long_old)*1.0E6<<" (pb)"<<endl;
        }
#endif        
        //,,,,,,,,,,,,,,,,,,compare old and new GetInteXS() end,,,,,,,,,,,,,,,,,,
      }
      if(DEBUG>=3) cout<<"  Averaged  InteXS_"<<L1Name[i]<<" = "<<pXS_tot_long/nTry*1.0E6<<" (pb)"<<endl;
#ifdef CMP_GETINTEXS        
      if(DEBUG>=3) cout<<"  OLD  Averaged  InteXS_"<<L1Name[i]<<" = "<<pXS_tot_long_old/nTry*1.0E6<<" (pb)"<<endl;
#endif      
    
      double pInteXs_long=pXS_tot_long/nTry;  //Need to get average, or change the luminosity during integration
      double pRate_long=pLumi*pInteXs_long;
      if(DEBUG>=3) cout<<"  pInteXs_"<<L1Name[i]<<" = "<<setw(8)<<pInteXs_long*1.0E6<<" (pb)  Rate="<<pRate_long<<" (Hz)"<<endl;

      pInteXs += pInteXs_long;
    }
    else {
      if(pElasOnly==2)
        pInteXs=GetInteXS(pBeamE, pDetectorAngle, pDetectorMomentum, L1Z[i], L1N[i], L1VZ[i], pDetectorName, pElasOnly,-1.,-1.,-1.,-1.,1.10,1.35);
      else if(pElasOnly==4)
        pInteXs=GetInteXS(pBeamE, pDetectorAngle, pDetectorMomentum, L1Z[i], L1N[i], L1VZ[i], pDetectorName, pElasOnly,-1.,-1.,-1.,-1.,2.00,100.0);
      else if(pElasOnly==5)
        pInteXs=GetInteXS(pBeamE, pDetectorAngle, pDetectorMomentum, L1Z[i], L1N[i], L1VZ[i], pDetectorName, pElasOnly,-1.,-1.,-1.,-1.,0.00,2.0);
      else
        pInteXs=GetInteXS(pBeamE, pDetectorAngle, pDetectorMomentum, L1Z[i], L1N[i], L1VZ[i], pDetectorName, pElasOnly);
    }
    pInteRate=pLumi*pInteXs;
    
    myRates[i]=pInteRate;  //store the result into array for return
  
    cout<<setw(4)<<pDetectorName<<setw(8)<<L1Name[i]
        <<setprecision(4)<<setw(10)<<pBeamE
        <<setw(10)<<setprecision(2)<<pBeamCurrent<<setw(10)<<setprecision(3)<<pDetectorMomentum
        <<setprecision(2)<<setw(10)<<pDetectorAngle/deg<<setw(10)<<L1VZ[i]
        <<setw(10)<<setprecision(4)<<L1Thick[i]
        <<setprecision(2)<<setw(12)<<pLumi/1000.<<" "<<setw(11)<<pInteXs*1.0E6<<" "<<setw(11)<<pInteRate
        <<endl;

    //print out the pressure curve rates, just scale them
    if(DEBUG>=1) {
      if(i>=4 && pDetectorMomentum>2.0 && pDetectorMomentum<2.2) {
        // Elastic (P_bP_t) pressure curve
        for(int iden=1;iden<=5;iden++) {
          cout.precision(3);
          cout<<"  ["<<L1Name[i]<<"] Pressure ="<<setw(6)<<2.0*iden<<" (atm)"<<" Lumi="<<setw(8)<<pLumi/1000.*0.2*iden<<" (10^33)"
              <<" pInteXs="<<setw(8)<<pInteXs*1.0E6<<" (pb)"<<" Rate="<<setw(8)<<pInteXs*pLumi*0.2*iden<<" (Hz)"<<endl;
        }
      }
    } 
  }   
  cout<<endl;
  cout.precision(pre);    //recover cout default precision
  cout.unsetf(std::ios::floatfield);  //recover cout floatfield properties
  
  return myRates;
}

////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////
//Beam Current in uA, All energies are in GeV unit. All angles are in radian unit.
double* GetOpticsRate(double pBeamCurrent, double pBeamE, double pDetectorAngle, double pDetectorMomentum, string pDetectorName, int pElasOnly=0)
{
  if(pElasOnly==-1)  cout<<"\n===================Full acceptance:  Pure Inelastic =======================\n";
  if(pElasOnly==0)   cout<<"\n===================Full acceptance:  Inelastic + Elastic ==================\n";
  if(pElasOnly==1)   cout<<"\n===================Full acceptance:  Pure elastic =========================\n";
  if(pElasOnly==2)   cout<<"\n===================Full acceptance:  With Delta (1.1<W<1.35) Cut ==========\n";
  if(pElasOnly==-30) cout<<"\n===================Only 2 SC bars:  Pure Inelastic ========================\n";
  if(pElasOnly==30)  cout<<"\n===================Only 2 SC bars:  Inelastic + Elastic ===================\n";
  if(pElasOnly==31)  cout<<"\n===================Only 2 SC bars:  Pure elastic ==========================\n";

  int pre = cout.precision();   //get the original precision of cout

  ////////////////////////////////////////////////////////////////////////////////
  //define material 
  double MolMass_c12=12.01078;  //in g/mol
  double Dc12=2.267;     // g/cm3;

   
  double pI_na=pBeamCurrent*1000.0;
  double pThickXDens=0.0; // (g/cm2);
  double pLumi[9], pInteXs[9], pInteRate[9];
  
  //this array will be returned to the caller, it has to be static, otherwise this array will be erased once this function finish 
  static double myRate[3];

  cout <<"\n Detector= "<<pDetectorName<<" Optics,  Beam_current= "<<pBeamCurrent<<" uA,  Angle= "<<pDetectorAngle/deg<<" deg"<<endl;
   
  //Loop1, 
  const char *L1Name[]={"C12_0","C12_1","C12_2","C12_3","C12_4","C12_5","C12_6","C12_3N5","C12_All"};
  int L1Nuclei_per_molecule[]={1, 1, 1, 1, 1, 1, 1, 1, 1};  //number of atoms per molecule, use 1 for all compounds, n for Hn, Dn
  double L1Z[]={6., 6., 6., 6., 6., 6., 6., 6., 6.};  //number of protons in one atom, for GE180, use negative z
  double L1N[]={6., 6., 6., 6., 6., 6., 6., 6., 6.};  //number of neutrons in one atom
  double L1Dens[]={Dc12, Dc12, Dc12, Dc12, Dc12, Dc12, Dc12, Dc12, Dc12};  //g/cm3
  double L1MolMass[]={MolMass_c12, MolMass_c12, MolMass_c12, MolMass_c12, MolMass_c12, MolMass_c12, MolMass_c12, MolMass_c12, MolMass_c12};  //g  
  double L1Thick[]={kC12FoilThick, kC12FoilThick, kC12FoilThick, kC12FoilThick, kC12FoilThick, kC12FoilThick, kC12FoilThick, kC12FoilThick*2, kC12FoilThick*7}; //cm
  double L1VZ[]={-20.0, -13.33, -6.67, 0.0, 6.67, 13.33, 20.0, 0.0, 0.0};  // vertex location in cm
  
    
  cout.setf(ios_base::fixed); 
  cout.precision(4);
  cout<<setw(4)<<"Det"<<setw(8)<<"Target"
      <<setw(10)<<"BeamE"<<setw(10)<<"Current"<<setw(10)<<"DetP0"
      <<setw(10)<<"DetAngle"<<setw(10)<<"VZ(cm)"<<setw(10)<<" Thick(cm)"
      <<setw(12)<<"Lumi(10^33)"<<setw(12)<<"InteXS(pb)"<<setw(12)<<"Rate(Hz)"
      <<endl;

      
  for(int i=0;i<9;i++) {
    pInteRate[i] =0.0;  //reset
    
    cout.precision(3);
    pThickXDens = L1Dens[i]*L1Thick[i]; // (g/cm2);
    pLumi[i] = GetLumi10pow30(pI_na,L1Nuclei_per_molecule[i],L1MolMass[i],pThickXDens);

    //only calculate each single foil, 2-foil and 7-foil cases will be caculated from the previous 7 results
    if(i<7) {
      pInteXs[i] = GetInteXS(pBeamE, pDetectorAngle, pDetectorMomentum, L1Z[i], L1N[i], L1VZ[i], pDetectorName, pElasOnly);
      pInteRate[i] = pLumi[i]*pInteXs[i];
    }
    else if(i==7)
    {
      pInteXs[7] = (pInteXs[3] + pInteXs[5])/2;
      pInteRate[7] = pInteRate[3] + pInteRate[5];
    }      
    else if(i==8)
    {
      pInteXs[8]=pInteRate[8]=0.0;
      for(int j=0;j<7;j++) {
        pInteXs[8] += pInteXs[j]/7;
        pInteRate[8] += pInteRate[j];
      }
    }
  
    cout<<setw(4)<<pDetectorName<<setw(8)<<L1Name[i]
        <<setprecision(4)<<setw(10)<<pBeamE
        <<setw(10)<<setprecision(2)<<pBeamCurrent<<setw(10)<<setprecision(3)<<pDetectorMomentum
        <<setprecision(2)<<setw(10)<<pDetectorAngle/deg<<setw(10)<<L1VZ[i]
        <<setw(10)<<setprecision(4)<<L1Thick[i]
        <<setprecision(2)<<setw(12)<<pLumi[i]/1000.<<" "<<setw(11)<<pInteXs[i]*1.0E6<<" "<<setw(11)<<pInteRate[i]
        <<endl;
    
  }   
  cout<<endl;
  cout.precision(pre);    //recover cout default precision
  cout.unsetf(std::ios::floatfield);  //recover cout floatfield properties

  //store for outut
  myRate[0]=pInteRate[3];
  myRate[1]=pInteRate[7];
  myRate[2]=pInteRate[8];
  return myRate;
}

////////////////////////////////////////////////////////////////////////////
//calcualte suggested beam current and prescale factor
void GetBeamPara(double pBeamI_uA, double pRate_Hz, double fixedBeamI_uA,
                 double& SuggestedBeamI, int& SuggestedPS, double& SuggestedRate)
{
  const double kDAQLimit = 4500.0;  //Hz
  const double kMaxBeamI = 30.0;    //uA
  const double kMinBeamI = 0.10;    //uA
  
  double theRate=pRate_Hz;
  double theBeamI=pBeamI_uA;
  
  //Get the Beam Current which will give the rate match to DAQLimit
  if(fixedBeamI_uA > 0) {
    SuggestedBeamI = fixedBeamI_uA;
  } else {
    SuggestedBeamI = theBeamI/theRate * kDAQLimit;
    //Keep it as accurate as 0.01 uA
    SuggestedBeamI = int(SuggestedBeamI*100.0)/100.;
  }

  //In case the suggested BeamI exceed maximum allowed beam current
  if( SuggestedBeamI > kMaxBeamI ) SuggestedBeamI = kMaxBeamI;
  if( SuggestedBeamI < kMinBeamI ) SuggestedBeamI = kMinBeamI;
  
  double tmpRate = SuggestedBeamI/theBeamI * theRate;

  //get the prescale factor, if it is not exceed N.15, still round it to N
  SuggestedPS = int(ceil(tmpRate/kDAQLimit-0.15));
  if(SuggestedPS==0) SuggestedPS=1;
  
  //apply the prescale factor to the rate
  SuggestedRate = tmpRate/SuggestedPS;
}


//////////////////////////////////////////////////////////////////////////////////////
//The input rate array contains rates for the following materials:{"C12","He3","GE-180","GE-180","H2","N2","He3"};
void GetBeamPara(double pBeamI_uA, double* pRate_HMS, double* pRate_SHMS, double* pRate_HMS_good, double* pRate_SHMS_good,
                 double pBeam, double pHMSAngle, double pHMSP0, double pSHMSAngle, double pSHMSP0,
                 const char *pStrHMS="", const char *pStrSHMS="")
{
  //the given pRate contains rates for the following target material
  const char *TgName[]={"C12_10Mil","He3_12AMG","GE-180_Up","GE-180-Dn","H2_10ATM","N2_10ATM","He3_10ATM"};
  double theRate;
  double fixedBeamI_uA=-1.0;
  double SuggestedBeamI=-1.0;
  int    SuggestedPS;
  double pRate, pRate_good;

  cout<<"\n Det     Target   Beam  Angle Momentum Current Prescale  DAQRate GoodRate  Comment"<<endl;
  for(int i=0;i<7;i++) {
    if(i==2 || i==3) continue;   //skip glass window
    fixedBeamI_uA = SuggestedBeamI = -1.0;
    //use HMS to determine the beam current, if HMS rate is zero, will use SHMS to determine beam current
    if(pRate_HMS[i] > 0) {
      if(i==0) {theRate = pRate_HMS[i];}
      else {theRate = pRate_HMS[i] + pRate_HMS[2] + pRate_HMS[3];}
      GetBeamPara(pBeamI_uA,theRate,fixedBeamI_uA, SuggestedBeamI, SuggestedPS, pRate);
      pRate_good = pRate_HMS_good[i]*SuggestedBeamI/pBeamI_uA/SuggestedPS;
      printf(" HMS %10s %6.3f %6.1f %8.3f %7.2f %8d %8.2f %8.2f  %s\n",
             TgName[i], pBeam, pHMSAngle, pHMSP0, SuggestedBeamI, SuggestedPS, pRate, pRate_good, pStrHMS);
    }
    if(pRate_SHMS[i] > 0) {
      fixedBeamI_uA = SuggestedBeamI;
      if(i==0) {theRate = pRate_SHMS[i];}
      else {theRate = pRate_SHMS[i] + pRate_SHMS[2] + pRate_SHMS[3];}
      GetBeamPara(pBeamI_uA,theRate,fixedBeamI_uA, SuggestedBeamI, SuggestedPS, pRate);
      pRate_good = pRate_SHMS_good[i]*SuggestedBeamI/pBeamI_uA/SuggestedPS;
      printf("SHMS %10s %6.3f %6.1f %8.3f %7.2f %8d %8.2f %8.2f  %s\n",
             TgName[i], pBeam, pSHMSAngle, pSHMSP0, SuggestedBeamI, SuggestedPS, pRate, pRate_good, pStrSHMS);
    }
  }
}


  
////////////////////////////////////////////////////////////////////////////
//print out beam time table, always use sieve
void GetOpticsBeamPara(double pBeamI_uA, double* pRate_HMS, double* pRate_SHMS, double* pRate_HMS_Elas, double* pRate_SHMS_Elas,
                       double pBeam, double pHMSAngle, double pHMSP0, double pSHMSAngle, double pSHMSP0,
                       const char *pStrHMS="", const char *pStrSHMS="")
{
  //the given pRate contains rates for the following target material
  const char *TgName[]={"C12_1Foil","C12_2Foil","C12_7Foil"};
  //to save time, require center hole to have 25 event for single foil, 50 events for 2-foil, 700 event for 7 foil 
  int kCenterHoleEvents_HMS[]={25,25*2,25*7};
  int kCenterHoleEvents_SHMS[]={25,25*2,25*7};
  double theRate;
  double fixedBeamI_uA=-1.0;
  double SuggestedBeamI=-1.0;
  int    SuggestedPS;
  double pRate, pRate_Elas, pRate_Elas_Sieve, pRate_Elas_CenterHole, pMinute;
  double kSieve2NoSieve_HMS=0.046;
  double kCenterHoleFrac_HMS=0.0072;
  double kSieve2NoSieve_SHMS=0.062;
  double kCenterHoleFrac_SHMS=0.0015;

  if(pHMSAngle>15) {
    kSieve2NoSieve_HMS=0.0477;
    kCenterHoleFrac_HMS=0.0065;
    kCenterHoleEvents_HMS[0]=100;
    kCenterHoleEvents_HMS[1]=200;
    kCenterHoleEvents_HMS[2]=700;
  }
  if(pSHMSAngle>15) {
    kSieve2NoSieve_SHMS=0.0661;
    kCenterHoleFrac_SHMS=0.0048;
    kCenterHoleEvents_SHMS[0]=100;
    kCenterHoleEvents_SHMS[1]=200;
    kCenterHoleEvents_SHMS[2]=700;
  }

  cout<<endl;
  cout<<" === HMS: Suggested beam Current and time to get "<<kCenterHoleEvents_HMS[0]<<"(1-Foil)|"
      <<kCenterHoleEvents_HMS[1]<<"(2-Foil)|"<<kCenterHoleEvents_HMS[2]<<"(7-Foil) events in center sieve hole==="<<endl;
  cout<<" ===SHMS: Suggested beam Current and time to get "<<kCenterHoleEvents_SHMS[0]<<"(1-Foil)|"
      <<kCenterHoleEvents_SHMS[1]<<"(2-Foil)|"<<kCenterHoleEvents_SHMS[2]<<"(7-Foil) events in center sieve hole===\n"<<endl;
  cout<<" Det     Target   Beam  Angle Momentum Current Prescale  DAQRate GoodRate SieveRate CenterHoleRate BeamTime(m)  Comment"<<endl;
  for(int i=0;i<3;i++) {
    fixedBeamI_uA = SuggestedBeamI = -1.0; //reset the trigger
    if(pRate_HMS[i] > 0) {
      theRate = pRate_HMS[i]*kSieve2NoSieve_HMS;
      GetBeamPara(pBeamI_uA,theRate,fixedBeamI_uA, SuggestedBeamI, SuggestedPS, pRate);
      pRate /= kSieve2NoSieve_HMS;
      pRate_Elas = pRate_HMS_Elas[i]*SuggestedBeamI/pBeamI_uA/SuggestedPS;
      pRate_Elas_Sieve = pRate_Elas*kSieve2NoSieve_HMS;
      pRate_Elas_CenterHole = pRate_Elas_Sieve*kCenterHoleFrac_HMS;
      pMinute = kCenterHoleEvents_HMS[i]/(pRate_Elas_CenterHole)/60.0;
      printf(" HMS %10s %6.3f %6.1f %8.3f %7.2f %8d %8.2f %8.2f %9.4f %14.4f %11.1f  %s\n",
             TgName[i], pBeam, pHMSAngle, pHMSP0, SuggestedBeamI, SuggestedPS,
             pRate, pRate_Elas, pRate_Elas_Sieve, pRate_Elas_CenterHole, pMinute, pStrHMS);
    }
    if(pRate_SHMS[i] > 0) {
      fixedBeamI_uA = SuggestedBeamI;
      theRate = pRate_SHMS[i]*kSieve2NoSieve_SHMS; 
      GetBeamPara(pBeamI_uA,theRate,fixedBeamI_uA, SuggestedBeamI, SuggestedPS, pRate);
      pRate /= kSieve2NoSieve_SHMS;
      pRate_Elas = pRate_SHMS_Elas[i]*SuggestedBeamI/pBeamI_uA/SuggestedPS;
      pRate_Elas_Sieve = pRate_Elas*kSieve2NoSieve_SHMS;
      pRate_Elas_CenterHole = pRate_Elas_Sieve*kCenterHoleFrac_SHMS;
      pMinute = kCenterHoleEvents_SHMS[i]/(pRate_Elas_CenterHole)/60.0;
      printf("SHMS %10s %6.3f %6.1f %8.3f %7.2f %8d %8.2f %8.2f %9.4f %14.4f %11.1f  %s\n",
             TgName[i], pBeam, pSHMSAngle, pSHMSP0, SuggestedBeamI, SuggestedPS,
             pRate, pRate_Elas, pRate_Elas_Sieve, pRate_Elas_CenterHole, pMinute, pStrSHMS);
    }
  }
}

////////////////////////////////////////////////////////////////////////////
void A1NRates()
{
  cout<<"\n====================Calculate rates for A1N using method 1 =================\n";
  const double kBeamI[]={1.0,1.0,30.0,30.0,30.0,30.0};  //in uA
  const double kBeamE[]={2.1,2.1,10.5,10.5,10.5,10.5};
  const double kHMSAngle[]={11.7,11.7,12.5,12.5,30.0,30.0};
  const double kHMSP0[]={2.068,1.682,5.7,6.8,2.9,3.5};
  const double kSHMSAngle[]={8.5,8.5,12.5,12.5,30.0,30.0};
  const double kSHMSP0[]={2.083,1.718,5.8,7.5,2.4,3.4};

  double *pRate;
  double pRate_HMS[7], pRate_HMS_good[7], pRate_SHMS[7], pRate_SHMS_good[7];
  for(int i=0;i<7;i++) pRate_HMS[i] = pRate_SHMS[i] = pRate_HMS_good[i] = pRate_SHMS_good[i] = 0.0;
  
  for(int j=0;j<6;j++) {    
    //only run C12 an dpressure curve for j==0
    if(j>0) {
      gSkipC12 = true;
      gSkipPresureCurve = true;
    }
    
    if(j==0) {
      //elas and PbPt
      
      //inelas + elas
      for(int i=0;i<7;i++) pRate_HMS[i] = pRate_SHMS[i] = pRate_HMS_good[i] = pRate_SHMS_good[i] = 0.0;  //reset
      pRate = GetRate(kBeamI[j],kBeamE[j],kHMSAngle[j]*deg,kHMSP0[j],"HMS",0);
      for(int i=0;i<7;i++) pRate_HMS[i] = pRate[i]; 
      pRate = GetRate(kBeamI[j],kBeamE[j],kSHMSAngle[j]*deg,kSHMSP0[j],"SHMS",0);
      for(int i=0;i<7;i++) pRate_SHMS[i] = pRate[i];
      
      //pure elas
      pRate = GetRate(kBeamI[j],kBeamE[j],kHMSAngle[j]*deg,kHMSP0[j],"HMS",1);
      for(int i=0;i<7;i++) pRate_HMS_good[i] = pRate[i]; 
      pRate = GetRate(kBeamI[j],kBeamE[j],kSHMSAngle[j]*deg,kSHMSP0[j],"SHMS",1);
      for(int i=0;i<7;i++) pRate_SHMS_good[i] = pRate[i];
      GetBeamPara(kBeamI[j],pRate_HMS,pRate_SHMS,pRate_HMS_good,pRate_SHMS_good,kBeamE[j],
                  kHMSAngle[j],kHMSP0[j],kSHMSAngle[j],kSHMSP0[j],"Full-Acc","Full-Acc");
      
      //pure inelastic
      //for(int i=0;i<7;i++) pRate_HMS[i] = pRate_SHMS[i] = 0.0;  //reset
      //pRate = GetRate(kBeamI[j],kBeamE[j],kHMSAngle[j]*deg,kHMSP0[j],"HMS",-1);  
      //for(int i=0;i<7;i++) pRate_HMS[i] = pRate[i]; 
      //pRate = GetRate(kBeamI[j],kBeamE[j],kSHMSAngle[j]*deg,kSHMSP0[j],"SHMS",-1);
      //for(int i=0;i<7;i++) pRate_SHMS[i] = pRate[i]; 
      
      //2-sc-bar only
      for(int i=0;i<7;i++) pRate_HMS[i] = pRate_SHMS[i] = pRate_HMS_good[i] = pRate_SHMS_good[i] = 0.0;  //reset
      pRate = GetRate(kBeamI[j],kBeamE[j],kHMSAngle[j]*deg,kHMSP0[j],"HMS",30);
      for(int i=0;i<7;i++) pRate_HMS[i] = pRate[i]; 
      pRate = GetRate(kBeamI[j],kBeamE[j],kSHMSAngle[j]*deg,kSHMSP0[j],"SHMS",30);
      for(int i=0;i<7;i++) pRate_SHMS[i] = pRate[i];
       
      //pure elas and 2-sc-bar only
      pRate = GetRate(kBeamI[j],kBeamE[j],kHMSAngle[j]*deg,kHMSP0[j],"HMS",31);
      for(int i=0;i<7;i++) pRate_HMS_good[i] = pRate[i]; 
      pRate = GetRate(kBeamI[j],kBeamE[j],kSHMSAngle[j]*deg,kSHMSP0[j],"SHMS",31);
      for(int i=0;i<7;i++) pRate_SHMS_good[i] = pRate[i]; 
      GetBeamPara(kBeamI[j],pRate_HMS,pRate_SHMS,pRate_HMS_good,pRate_SHMS_good,kBeamE[j],
                  kHMSAngle[j],kHMSP0[j],kSHMSAngle[j],kSHMSP0[j],"2-Bar-Acc","2-Bar-Acc");
    }
    else if(j==1) {
      //Delta peak
      for(int i=0;i<7;i++) pRate_HMS[i] = pRate_SHMS[i] = pRate_HMS_good[i] = pRate_SHMS_good[i] = 0.0;  //reset
      pRate = GetRate(kBeamI[j],kBeamE[j],kHMSAngle[j]*deg,kHMSP0[j],"HMS",0);
      for(int i=0;i<7;i++) pRate_HMS[i] = pRate[i]; 
      pRate = GetRate(kBeamI[j],kBeamE[j],kSHMSAngle[j]*deg,kSHMSP0[j],"SHMS",0);
      for(int i=0;i<7;i++) pRate_SHMS[i] = pRate[i];
      
      //apply W cuts for Delta peak
      pRate = GetRate(kBeamI[j],kBeamE[j],kHMSAngle[j]*deg,kHMSP0[j],"HMS",2);
      for(int i=0;i<7;i++) pRate_HMS_good[i] = pRate[i]; 
      pRate = GetRate(kBeamI[j],kBeamE[j],kSHMSAngle[j]*deg,kSHMSP0[j],"SHMS",2);
      for(int i=0;i<7;i++) pRate_SHMS_good[i] = pRate[i]; 
      GetBeamPara(kBeamI[j],pRate_HMS,pRate_SHMS,pRate_HMS_good,pRate_SHMS_good,kBeamE[j],
                  kHMSAngle[j],kHMSP0[j],kSHMSAngle[j],kSHMSP0[j],"1.1<W<1.35","1.1<W<1.35");
    }
    else if (j==3) {
      //DIS
      for(int i=0;i<7;i++) pRate_HMS[i] = pRate_SHMS[i] = pRate_HMS_good[i] = pRate_SHMS_good[i] = 0.0;  //reset
      pRate = GetRate(kBeamI[j],kBeamE[j],kHMSAngle[j]*deg,kHMSP0[j],"HMS",0);
      for(int i=0;i<7;i++) pRate_HMS[i] = pRate[i];    
      pRate = GetRate(kBeamI[j],kBeamE[j],kSHMSAngle[j]*deg,kSHMSP0[j],"SHMS",0);
      for(int i=0;i<7;i++) pRate_SHMS[i] = pRate[i];
      
      //apply W cuts for DIS and Resonance for SHMS kine D
      pRate = GetRate(kBeamI[j],kBeamE[j],kHMSAngle[j]*deg,kHMSP0[j],"HMS",4);
      for(int i=0;i<7;i++) pRate_HMS_good[i] = pRate[i]; 
      pRate = GetRate(kBeamI[j],kBeamE[j],kSHMSAngle[j]*deg,kSHMSP0[j],"SHMS",5);
      for(int i=0;i<7;i++) pRate_SHMS_good[i] = pRate[i];
      GetBeamPara(kBeamI[j],pRate_HMS,pRate_SHMS,pRate_HMS_good,pRate_SHMS_good,kBeamE[j],
                  kHMSAngle[j],kHMSP0[j],kSHMSAngle[j],kSHMSP0[j],"W>2.0","W<2.0");
   }     
    else {
      //DIS
      for(int i=0;i<7;i++) pRate_HMS[i] = pRate_SHMS[i] = pRate_HMS_good[i] = pRate_SHMS_good[i] = 0.0;  //reset
      pRate = GetRate(kBeamI[j],kBeamE[j],kHMSAngle[j]*deg,kHMSP0[j],"HMS",0);
      for(int i=0;i<7;i++) pRate_HMS[i] = pRate[i];    
      pRate = GetRate(kBeamI[j],kBeamE[j],kSHMSAngle[j]*deg,kSHMSP0[j],"SHMS",0);
      for(int i=0;i<7;i++) pRate_SHMS[i] = pRate[i];
      
      //apply W cuts for DIS
      pRate = GetRate(kBeamI[j],kBeamE[j],kHMSAngle[j]*deg,kHMSP0[j],"HMS",4);
      for(int i=0;i<7;i++) pRate_HMS_good[i] = pRate[i]; 
      pRate = GetRate(kBeamI[j],kBeamE[j],kSHMSAngle[j]*deg,kSHMSP0[j],"SHMS",4);
      for(int i=0;i<7;i++) pRate_SHMS_good[i] = pRate[i]; 
      GetBeamPara(kBeamI[j],pRate_HMS,pRate_SHMS,pRate_HMS_good,pRate_SHMS_good,kBeamE[j],
                  kHMSAngle[j],kHMSP0[j],kSHMSAngle[j],kSHMSP0[j],"W>2.0","W>2.0");
    }
  }
}


void D2NRates()
{
  cout<<"\n===================Calculate rates for d2n using method 1 =================\n";
  const double kBeamI[]    ={30.0, 30.0, 30.0, 30.0};  //in uA
  const double kBeamE[]    ={10.5, 10.5, 10.5, 10.5};
  const double kHMSAngle[] ={13.5, 16.4, 20.0, 25.0};
  const double kHMSP0[]    ={ 4.3,  5.1,  4.0,  2.5};
  const double kSHMSAngle[]={11.0, 13.3, 15.5, 18.0};
  const double kSHMSP0[]   ={ 7.5,  7.0,  6.3,  5.6};

  double *pRate;
  double pRate_HMS[7], pRate_SHMS[7];
  for(int i=0;i<7;i++) pRate_HMS[i] = pRate_SHMS[i] = 0.0;
  
  gSkipC12 = true;
  gSkipPresureCurve = true;
  
  for(int j=0;j<4;j++) {
    //DIS
    for(int i=0;i<7;i++) pRate_HMS[i] = pRate_SHMS[i] = 0.0;  //reset
    pRate = GetRate(kBeamI[j],kBeamE[j],kHMSAngle[j]*deg,kHMSP0[j],"HMS",0);
    for(int i=0;i<7;i++) pRate_HMS[i] = pRate[i];    
    pRate = GetRate(kBeamI[j],kBeamE[j],kSHMSAngle[j]*deg,kSHMSP0[j],"SHMS",0);
    for(int i=0;i<7;i++) pRate_SHMS[i] = pRate[i];
    GetBeamPara(kBeamI[j],pRate_HMS,pRate_SHMS,pRate_HMS,pRate_SHMS,kBeamE[j],
                kHMSAngle[j],kHMSP0[j],kSHMSAngle[j],kSHMSP0[j],"Full-Acc","Full-Acc");
    /*
    //apply W cuts for DIS
    for(int i=0;i<7;i++) pRate_HMS[i] = pRate_SHMS[i] = 0.0;  //reset
    pRate = GetRate(kBeamI[j],kBeamE[j],kHMSAngle[j]*deg,kHMSP0[j],"HMS",4);
    for(int i=0;i<7;i++) pRate_HMS[i] = pRate[i]; 
    pRate = GetRate(kBeamI[j],kBeamE[j],kSHMSAngle[j]*deg,kSHMSP0[j],"SHMS",4);
    for(int i=0;i<7;i++) pRate_SHMS[i] = pRate[i]; 
    */
  }
}

void He3FormFactorRates()
{
  cout<<"\n===================Calculate rates for 3He Form Factor=================\n";
  //if P0 less than 0.5, this point will not run
  const double kBeamI[]    ={ 1.0, 1.0, 1.0, 1.0, 1.0};  //in uA
  const double kBeamE[]    ={ 2.1, 2.1, 2.1, 2.1, 2.1};
  const double kHMSAngle[] ={11.0,13.0,15.0,17.0,19.0};
  const double kSHMSAngle[]={21.0,21.0,21.0,21.0,21.0};
  const double kHMSP0[]    ={2.0714, 2.0603, 2.0476, 2.0333, 2.0174};
  const double kSHMSP0[]   ={2.0002, 2.0002, 2.0002, 2.0002, 2.0002};

  double *pRate;
  double pRate_HMS[7], pRate_HMS_Elas[7], pRate_SHMS[7], pRate_SHMS_Elas[7];

  //skip C12 and pressure curve
  gSkipC12 = true;
  gSkipPresureCurve = true;
  
  bool bSkipHMS = false;
  bool bSkipSHMS = false;
  for(int j=10;j<5;j++) {
    //Full-Acc
    //for(int i=0;i<7;i++) pRate_HMS[i]=pRate_SHMS[i]=pRate_HMS_Elas[i]=pRate_SHMS_Elas[i]=0.0;  //reset
    if(j>0) bSkipSHMS=true;
    for(int i=0;i<7;i++) pRate_HMS[i]=pRate_HMS_Elas[i]=0.0;  //do not reset SHMS
    if(kHMSP0[j]>0.5 && !bSkipHMS) {
      pRate = GetRate(kBeamI[j],kBeamE[j],kHMSAngle[j]*deg,kHMSP0[j],"HMS",0);
      for(int i=0;i<7;i++) pRate_HMS[i] = pRate[i];
    }
    if(kSHMSP0[j]>0.5 && !bSkipSHMS) {
      pRate = GetRate(kBeamI[j],kBeamE[j],kSHMSAngle[j]*deg,kSHMSP0[j],"SHMS",0);
      for(int i=0;i<7;i++) pRate_SHMS[i] = pRate[i];
    }
    if(kHMSP0[j]>0.5 && !bSkipHMS) {
      pRate = GetRate(kBeamI[j],kBeamE[j],kHMSAngle[j]*deg,kHMSP0[j],"HMS",1);
      for(int i=0;i<7;i++) pRate_HMS_Elas[i] = pRate[i];
    }
    if(kSHMSP0[j]>0.5 && !bSkipSHMS) {
      pRate = GetRate(kBeamI[j],kBeamE[j],kSHMSAngle[j]*deg,kSHMSP0[j],"SHMS",1);
      for(int i=0;i<7;i++) pRate_SHMS_Elas[i] = pRate[i];
    }
    GetBeamPara(kBeamI[j],pRate_HMS,pRate_SHMS,pRate_HMS_Elas,pRate_SHMS_Elas,kBeamE[j],
                kHMSAngle[j],kHMSP0[j],kSHMSAngle[j],kSHMSP0[j],"Full-Acc","Full-Acc");
  }
  
  bSkipSHMS=false;
  for(int j=0;j<5;j++) {
    //2-Bar-Acc
    //for(int i=0;i<7;i++) pRate_HMS[i]=pRate_SHMS[i]=pRate_HMS_Elas[i]=pRate_SHMS_Elas[i]=0.0;  //reset
    if(j>0) bSkipSHMS=true;
    for(int i=0;i<7;i++) pRate_HMS[i]=pRate_HMS_Elas[i]=0.0;  //do not reset SHMS
    if(kHMSP0[j]>0.5 && !bSkipHMS) {
      pRate = GetRate(kBeamI[j],kBeamE[j],kHMSAngle[j]*deg,kHMSP0[j],"HMS",30);
      for(int i=0;i<7;i++) pRate_HMS[i] = pRate[i];
    }
    if(kSHMSP0[j]>0.5 && !bSkipSHMS) {
      pRate = GetRate(kBeamI[j],kBeamE[j],kSHMSAngle[j]*deg,kSHMSP0[j],"SHMS",30);
      for(int i=0;i<7;i++) pRate_SHMS[i] = pRate[i];
    }
    if(kHMSP0[j]>0.5 && !bSkipHMS) {
      pRate = GetRate(kBeamI[j],kBeamE[j],kHMSAngle[j]*deg,kHMSP0[j],"HMS",31);
      for(int i=0;i<7;i++) pRate_HMS_Elas[i] = pRate[i];
    }
    if(kSHMSP0[j]>0.5 && !bSkipSHMS) {
      pRate = GetRate(kBeamI[j],kBeamE[j],kSHMSAngle[j]*deg,kSHMSP0[j],"SHMS",31);
      for(int i=0;i<7;i++) pRate_SHMS_Elas[i] = pRate[i];
    }
    GetBeamPara(kBeamI[j],pRate_HMS,pRate_SHMS,pRate_HMS_Elas,pRate_SHMS_Elas,kBeamE[j],
                kHMSAngle[j],kHMSP0[j],kSHMSAngle[j],kSHMSP0[j],"2-Bar-Acc","2-Bar-Acc");
  }
}


void A1NOptics()
{
  cout<<"\n===================Calculate rates for A1N|d2n Optics =================\n";
  //if P0 less than 0.5, this point will not run
  const double kBeamI[]    ={ 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 30.0, 30.0, 30.0};  //in uA
  const double kBeamE[]    ={ 2.1, 2.1, 2.1, 2.1, 2.1, 2.1, 2.1, 10.5, 10.5, 4.20};
  const double kHMSAngle[] ={11.7,11.7,11.7,11.7,11.7,11.7,11.7, 30.0, 35.0, 35.0};
  const double kSHMSAngle[]={ 8.5, 8.5, 8.5, 8.5, 8.5, 8.5, 8.5, 30.0, 35.0, 35.0};
  const double kHMSP0[]    ={2.0918, 2.0918/1.05, 2.0918/1.08,         0.0,         0.0, 2.0918/0.92, 2.0918/0.95, 1.5, 1.5, 1.5};
  const double kSHMSP0[]   ={2.0957, 2.0957/1.05, 2.0957/1.10, 2.0957/1.15, 2.0957/0.85, 2.0957/0.90, 2.0957/0.95, 1.8, 1.8, 1.8};
  const char* kCommentHMS[] = {"dP=0%","dP=+5%","dP=+8%",        "",       "", "dP=-8%","dP=-5%","","",""};
  const char* kCommentSHMS[]= {"dP=0%","dP=+5%","dP=+10%","dP=+15%","dP=-15%","dP=-10%","dP=-5%","","",""};
  const char* kCommentHMS_2Bar[] = {"dP=0%, 2-Bar-Acc","dP=+5%, 2-Bar-Acc", "dP=+8%, 2-Bar-Acc",                  "",                  "", "dP=-8%, 2-Bar-Acc","dP=-5%, 2-Bar-Acc","","",""};
  const char* kCommentSHMS_2Bar[]= {"dP=0%, 2-Bar-Acc","dP=+5%, 2-Bar-Acc","dP=+10%, 2-Bar-Acc","dP=+15%, 2-Bar-Acc","dP=-15%, 2-Bar-Acc","dP=-10%, 2-Bar-Acc","dP=-5%, 2-Bar-Acc","","",""};


  double *pRate;
  double pRate_HMS[3], pRate_HMS_Elas[3], pRate_SHMS[3], pRate_SHMS_Elas[3];
  
  for(int j=0;j<7;j++) {
    //C12 delta scan
    
    for(int i=0;i<3;i++) pRate_HMS[i]=pRate_SHMS[i]=pRate_HMS_Elas[i]=pRate_SHMS_Elas[i]=0.0;  //reset
    if(kHMSP0[j]>0.5) {
      pRate = GetOpticsRate(kBeamI[j],kBeamE[j],kHMSAngle[j]*deg,kHMSP0[j],"HMS",0);
      for(int i=0;i<3;i++) pRate_HMS[i] = pRate[i];
    }
    if(kSHMSP0[j]>0.5) {
      pRate = GetOpticsRate(kBeamI[j],kBeamE[j],kSHMSAngle[j]*deg,kSHMSP0[j],"SHMS",0);
      for(int i=0;i<3;i++) pRate_SHMS[i] = pRate[i];
    }
    if(kHMSP0[j]>0.5) {
      pRate = GetOpticsRate(kBeamI[j],kBeamE[j],kHMSAngle[j]*deg,kHMSP0[j],"HMS",1);
      for(int i=0;i<3;i++) pRate_HMS_Elas[i] = pRate[i];
    }
    if(kSHMSP0[j]>0.5) {
      pRate = GetOpticsRate(kBeamI[j],kBeamE[j],kSHMSAngle[j]*deg,kSHMSP0[j],"SHMS",1);
      for(int i=0;i<3;i++) pRate_SHMS_Elas[i] = pRate[i];
    }
    GetOpticsBeamPara(kBeamI[j],pRate_HMS,pRate_SHMS,pRate_HMS_Elas,pRate_SHMS_Elas,kBeamE[j],
                      kHMSAngle[j],kHMSP0[j],kSHMSAngle[j],kSHMSP0[j],kCommentHMS[j],kCommentSHMS[j]);
    
    if(j==0)
    {
      for(int i=0;i<3;i++) pRate_HMS[i]=pRate_SHMS[i]=pRate_HMS_Elas[i]=pRate_SHMS_Elas[i]=0.0;  //reset
      pRate = GetOpticsRate(kBeamI[j],kBeamE[j],kHMSAngle[j]*deg,kHMSP0[j],"HMS",30);
      for(int i=0;i<3;i++) pRate_HMS[i] = pRate[i];    
      pRate = GetOpticsRate(kBeamI[j],kBeamE[j],kSHMSAngle[j]*deg,kSHMSP0[j],"SHMS",30);
      for(int i=0;i<3;i++) pRate_SHMS[i] = pRate[i];    
      pRate = GetOpticsRate(kBeamI[j],kBeamE[j],kHMSAngle[j]*deg,kHMSP0[j],"HMS",31);
      for(int i=0;i<3;i++) pRate_HMS_Elas[i] = pRate[i];    
      pRate = GetOpticsRate(kBeamI[j],kBeamE[j],kSHMSAngle[j]*deg,kSHMSP0[j],"SHMS",31);
      for(int i=0;i<3;i++) pRate_SHMS_Elas[i] = pRate[i];
      GetOpticsBeamPara(kBeamI[j],pRate_HMS,pRate_SHMS,pRate_HMS_Elas,pRate_SHMS_Elas,kBeamE[j],
                        kHMSAngle[j],kHMSP0[j],kSHMSAngle[j],kSHMSP0[j],kCommentHMS_2Bar[j],kCommentSHMS_2Bar[j]);
    } 
  }
  for(int j=7;j<10;j++) {
    //angle and vz calibration
    for(int i=0;i<3;i++) pRate_HMS[i]=pRate_SHMS[i]=pRate_HMS_Elas[i]=pRate_SHMS_Elas[i]=0.0;  //reset
    if(kHMSP0[j]>0.5) {
      pRate = GetOpticsRate(kBeamI[j],kBeamE[j],kHMSAngle[j]*deg,kHMSP0[j],"HMS",0);
      for(int i=0;i<3;i++) pRate_HMS[i] = pRate[i];
    }
    if(kSHMSP0[j]>0.5) {
      pRate = GetOpticsRate(kBeamI[j],kBeamE[j],kSHMSAngle[j]*deg,kSHMSP0[j],"SHMS",0);
      for(int i=0;i<3;i++) pRate_SHMS[i] = pRate[i];
    }
    
    for(int i=0;i<3;i++) {
      //for this setting, no elas event, just physicsal event
      pRate_HMS_Elas[i]=pRate_HMS[i];
      pRate_SHMS_Elas[i]=pRate_SHMS[i];
    }
    GetOpticsBeamPara(kBeamI[j],pRate_HMS,pRate_SHMS,pRate_HMS_Elas,pRate_SHMS_Elas,kBeamE[j],
                      kHMSAngle[j],kHMSP0[j],kSHMSAngle[j],kSHMSP0[j]);
  }
}


/*
//check E' @ W=1232 peak
double W=1.232;
double M=0.9383;
double sin2th = pow(sin(11.7/57.3/2.),2.);
double E0=2.10;
double nu = (W*W-M*M+4*E0*E0*sin2th)/(2*M+4*E0*sin2th);
cout<<" E0="<<E0<<",  nu="<<nu<<"  ==> E'="<<E0-nu<<"\n";

double W=1.232;
double M=0.9383;
double sin2th = pow(sin(8.5/57.3/2.),2.);
double E0=2.10;
double nu = (W*W-M*M+4*E0*E0*sin2th)/(2*M+4*E0*sin2th);
cout<<" E0="<<E0<<",  nu="<<nu<<"  ==> E'="<<E0-nu<<"\n";


double Ep = 1.682;
double Q2 = 4*E0*Ep*sin2th;
double W2 =-Q2+M*M+2*M*(E0-Ep);
cout<<" E0="<<E0<<",  E'="<<Ep<<"  ==> W="<<sqrt(W2)<<"\n";
*/
