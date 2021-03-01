//SHMSXSTree.h 
// create my own tree from mc-single-arm
// ---------------------------------------------

#ifndef _SHMSXSTree_
#define _SHMSXSTree_

#include <stdio.h>
#include <vector>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <exception>
#include <TROOT.h>
#include <TStyle.h>
#include <TF1.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TFile.h>
#include <TDirectory.h>
#include <TObject.h>
#include <TVector3.h>
#include <TLorentzVector.h>
#include <TCanvas.h>
#include <TPad.h>
#include <TProfile.h>
#include <TString.h>
#include "ReadSHMS.h"

using namespace std;


#define SHMSXSTree_Debug 4
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
class SHMSXSTree:ReadSHMS
{
public:
    SHMSXSTree(const char* filename);
    virtual ~SHMSXSTree();

    void  BeginOfRun();
    void  Run(int nevents_to_process=-1);
    void  EndOfRun();
    ////////////////////////////////////////////////////////////////////
    //P. Bosted XS model wrapper
    double GetXS(int Z, int N, float Ei, float Ef, float Theta, float Tb=-0.000, float Ta=-0.000);
    double GetXS_GE180(float Ei, float Ef, float Theta, float Tb=-0.000, float Ta=-0.000);
    
public:
    void SetTreeLevel(int pLevel) {mTreeLevel=pLevel;cout<<"Update SHMSXSTree::mTreeLevel to "<<mTreeLevel<<endl; };
    int  GetTreeLevel() {return mTreeLevel;};
    void Reset();  //reset the tree variables    
   
    void SetElas(int v) {mElasOnly=v;};
    void SetPara(float pBeam, float pTheta, float pMom) {mBeam=pBeam;mDetAngle=pTheta;mDetMom=pMom;};
    //================================================================
private:
    //mTreeLevel controls the tree output level
    //0 fill all into the tree
    //1 skip 2-particle events, keep only D(e,e'pi-)X or D(e,e'Ps)X events
    //2 keep only D(e,e'pi-)X events and 0.7<Missing Mass<1.2 only
    int mTreeLevel;
    int mElasOnly;

private:
    TFile*   mFile;
    TTree*   mTree;
    TString  mOutFileName;
    
    //-----------------------------block3 start  -----------------------
    //ini
    float vx0,vy0,vz0,p0,theta0,phi0;
    float xtg0,ytg0,ztg0,xptg0,yptg0,delta0;

    //intermedia
    int    ixs,iys;  //sieve hole index
    float  xsieve,ysieve;
    float  xfp,yfp,xpfp,ypfp;
    int    istop;

    //recon
    float  vx,vy,vz,p,theta,phi;
    float  xtg,ytg,ztg,xptg,yptg,delta;

    //kin
    float  nu,Q2,W,xbj;
    //rate
    float p_rate_4he, p_rate_2h, p_yield_4he, p_yield_2h;
    float p_rate_4he_pi, p_yield_4he_pi;
    float tar_len_4he, tar_len_2h, dens_4he, dens_2h;
    float p_spec, p_accept, th_accept, ph_accept, n_trials;

    //-----------------------------block3 end  -------------------------

    //-----------------------------block2 start  -----------------------
    //this block will not be filled if mTreeLevel>=2
    float  xs_1h,xs_3he,xs_4he,xs_12c,xs_14n,xs_27al,xs_ge180, xs_4he_pi; 
    //-----------------------------block2 end  -------------------------
    
    //-----------------------------block1 start  -----------------------
    //this block will not be filled if mTreeLevel>=1
    int    mDet;
    float  mBeam,mDetAngle,mDetMom;
    //-----------------------------block1 end  -------------------------

};
#endif //_SHMSXSTree_
