//implement of  SHMSXSTree class
#include "SHMSXSTree.h"
#include "HRSTransform_TCSNHCS.hh"
#include "XSModel.hh"   //XSModel return DiffXS in ub/MeV-Sr for inelas or ub for elas.
#include "G2PRand.hh"
#include <math.h>

////////////////////////////////////////////////////////////////////////
SHMSXSTree::SHMSXSTree(const char* filename):
  ReadSHMS(filename)
{
#ifdef SHMSXSTree_Debug
  if (SHMSXSTree_Debug>=1) cout<<" >>>>>Init SHMSXSTree ... <<<<<"<<endl;
#endif
  mElasOnly=0;
  mTreeLevel=1;
  mOutFileName = filename;
  mOutFileName.ReplaceAll(".root","_xstree.root");

  mDet=2;  //SHMS
}

////////////////////////////////////////////////////////////////////////
SHMSXSTree::~SHMSXSTree()
{
}

////////////////////////////////////////////////////////////////////////
void SHMSXSTree::BeginOfRun()
{
#ifdef SHMSXSTree_Debug
  if (SHMSXSTree_Debug>=4) cout<<" SHMSXSTree::BeginOfRun()"<<endl;
#endif
  mFile = new TFile(mOutFileName.Data(),"recreate","Jixie's mc-single-arm tree");
  /////////////////////////////////////////////////////////////////////
  mTree=new TTree("T","Jixie's tree based on mc-single-arm hms or shms");

  //ini
  mTree->Branch("vx0",&vx0,"vx0/F");
  mTree->Branch("vy0",&vy0,"vy0/F");
  mTree->Branch("vz0",&vz0,"vz0/F");
  mTree->Branch("p0",&p0,"p0/F");
  mTree->Branch("theta0",&theta0,"theta0/F");
  mTree->Branch("phi0",&phi0,"phi0/F");
    
  mTree->Branch("xtg0",&xtg0,"xtg0/F");
  mTree->Branch("ytg0",&ytg0,"ytg0/F");
  mTree->Branch("ztg0",&ztg0,"ztg0/F");
  mTree->Branch("xptg0",&xptg0,"xptg0/F");
  mTree->Branch("yptg0",&yptg0,"yptg0/F");
  mTree->Branch("delta0",&delta0,"delta0/F");
    
  //intermedia
  mTree->Branch("ixs",&ixs,"ixs/I");
  mTree->Branch("iys",&iys,"iys/I");
        
  mTree->Branch("xsieve",&xsieve,"xsieve/F");
  mTree->Branch("ysieve",&ysieve,"ysieve/F");
    
  mTree->Branch("xfp",&xfp,"xfp/F");
  mTree->Branch("yfp",&yfp,"yfp/F");
  mTree->Branch("xpfp",&xpfp,"xpfp/F");
  mTree->Branch("ypfp",&ypfp,"ypfp/F");
    
  mTree->Branch("istop",&istop,"istop/I");

  //recon
  mTree->Branch("vx",&vx,"vx/F");
  mTree->Branch("vy",&vy,"vy/F");
  mTree->Branch("vz",&vz,"vz/F");
  mTree->Branch("p",&p,"p/F");
  mTree->Branch("theta",&theta,"theta/F");
  mTree->Branch("phi",&phi,"phi/F");
    
  mTree->Branch("xtg",&xtg,"xtg/F");
  mTree->Branch("ytg",&ytg,"ytg/F");
  mTree->Branch("ztg",&ztg,"ztg/F");
  mTree->Branch("xptg",&xptg,"xptg/F");
  mTree->Branch("yptg",&yptg,"yptg/F");
  mTree->Branch("delta",&delta,"delta/F");

  //kin
  mTree->Branch("nu",&nu,"nu/F");
  mTree->Branch("Q2",&Q2,"Q2/F");
  mTree->Branch("W",&W,"W/F");
  mTree->Branch("xbj",&xbj,"xbj/F");
    

  //the following will not be filled if mTreeLevel>=2/////////////////
  if(mTreeLevel<2)
    {
      mTree->Branch("xs_1h",&xs_1h,"xs_1h/F");
      mTree->Branch("xs_3he",&xs_3he,"xs_3he/F");
      mTree->Branch("xs_4he",&xs_4he,"xs_4he/F");
      mTree->Branch("xs_12c",&xs_12c,"xs_12c/F");
      mTree->Branch("xs_14n",&xs_14n,"xs_14n/F");
      mTree->Branch("xs_27al",&xs_27al,"xs_27al/F");
      mTree->Branch("xs_ge180",&xs_ge180,"xs_ge180/F");
      mTree->Branch("xs_4he_pi",&xs_4he_pi,"xs_4he_pi/F");
      mTree->Branch("p_rate_4he_pi",&p_rate_4he_pi,"p_rate_4he_pi/F");
      mTree->Branch("p_yield_4he_pi",&p_yield_4he_pi,"p_yield_4he_pi/F");
      mTree->Branch("p_rate_4he",&p_rate_4he,"p_rate_4he/F");
      mTree->Branch("p_rate_2h",&p_rate_2h,"p_rate_2h/F");
      mTree->Branch("p_yield_4he",&p_yield_4he,"p_yield_4he/F");
      mTree->Branch("p_yield_2h",&p_yield_2h,"p_yield_2h/F");
    }
  //the above will not be filled if mTreeLevel>=2 /////////////////////


  //the following will not be filled if mTreeLevel>=1/////////////////
  if(mTreeLevel<1)
    {
      mTree->Branch("mDet",&mDet,"mDet/I");
      mTree->Branch("mBeam",&mBeam,"mBeam/F");
      mTree->Branch("mDetAngle",&mDetAngle,"mDetAngle/F");
      mTree->Branch("mDetMom",&mDetMom,"mDetMom/F");
    }
  //the above will not be filled if mTreeLevel>=1 /////////////////////

  //set the maximum tree size, without this line it can not exceed 2G
  mTree->SetMaxTreeSize((Long64_t)(20000000000.0)); //20G
}

/////////////////////////////////////////////////////////////////////////////////////
void SHMSXSTree::Run(int nevents_to_process)
{
  BeginOfRun();
  const double degree = asin(1.0)/90.0;
#ifdef SHMSXSTree_Debug 
  if(SHMSXSTree_Debug>=4)
    cout<<" SHMSXSTree::Run() "<<endl;
#endif
  const double kMp = 0.9383;
  Long64_t nentries_save = 0;
    
  Long64_t nentries = fChain->GetEntriesFast();
  if(nevents_to_process>0 && nevents_to_process<nentries) nentries=nevents_to_process;
  
#ifdef SHMSXSTree_Debug 
  if(SHMSXSTree_Debug>=5) nentries=100;
#endif

  for (Long64_t jentry=0; jentry<nentries;jentry++) {
#ifdef SHMSXSTree_Debug 
    if(SHMSXSTree_Debug>=5)
      cout<<" SHMSXSTree::Run() is processing event "<<std::setw(6)<<jentry<<"\n";
#endif

#ifdef SHMSXSTree_Debug 
    if( ((jentry+1)%10000) == 0)
      cout<<" SHMSXSTree::Run() is processing event "<<std::setw(6)<<jentry+1<<" ... \n";
#endif
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    fChain->GetEntry(jentry);
      
    //apply cuts
    //if(mTreeLevel>=3 && stop_id!=0) continue;

    //set variables

    //ini
    xtg0 = psxtari;
    ytg0 = psytari;
    ztg0 = psztari;
    xptg0 = psxptari;
    yptg0 = psyptari;
    delta0 = psdeltai;
    
    vx0 = vxi;   //include spectrometer offset
    vy0 = vyi;   //include spectrometer offset
    vz0 = vzi;   //include spectrometer offset
    p0 = (1.0 + delta0/100.) * mDetMom;
    double tmpTh=0.0,tmpPh=0.0;
    Transform::P_TCS2HCS(xptg0,yptg0,mDetAngle,tmpTh,tmpPh);
    theta0=tmpTh;
    phi0=tmpPh;
    if(p0>=mBeam) continue;

    //intermedia
    ixs = lroundf(xsnum);
    iys = lroundf(ysnum);  //sieve hole index
    xsieve = ReadSHMS::xsieve;
    ysieve = ReadSHMS::ysieve;
        
    xfp = psxfp;
    yfp = psyfp;
    xpfp = psxpfp;
    ypfp = psypfp;
    istop = lroundf(stop_id);

    //recon
    xtg = 0; //not available from original tree
    ytg = psytar;
    ztg = psztar;
    xptg = psxptar;
    yptg = psyptar;
    delta = psdelta;
    
    vx = vx0 + G2PRand::Gaus(0.0,0.02);   //add bpm resolution 0.02cm
    vy = vy0 + G2PRand::Gaus(0.0,0.02);   //add bpm resolution 0.02cm
    //vz could not be determined if vx is not known,try to use whatever vx above
    double tan2phitg = yptg*yptg;
    double sinphitg = sqrt(tan2phitg/(1.0+tan2phitg));
    double sintheangle = sin(mDetAngle - atan(yptg));
    double theside = (ytg-vx/cos(mDetAngle))/tan(mDetAngle);
        
    double vz_part0 = vy*tan(mDetAngle);
    double vz_part1 = (ytg-vx/cos(mDetAngle))/sin(mDetAngle);
    double vz_part2 = (theside/sintheangle) * sinphitg;
    vz = -(vz_part0 + vz_part1 + vz_part2);
        
    p = (1.0 + delta/100.) * mDetMom;
    Transform::P_TCS2HCS(xptg,yptg,mDetAngle,tmpTh,tmpPh);
    theta=tmpTh;
    phi=tmpPh;
        
    //kin
    nu = mBeam - p;
    Q2 = 2*mBeam*p*(1-cos(theta));
    W = sqrt(kMp*kMp - Q2 + 2*kMp*nu);
    xbj = Q2/(2*kMp*nu);

#ifdef SHMSXSTree_Debug 
    if(SHMSXSTree_Debug>=5) {
      cout<<" p="<<p<<"  xptg="<<xptg<<"  yptg="<<yptg<<"  theta_deg="<<theta/degree
          <<"  mDetAngle="<<mDetAngle/degree<<"  phi_deg="<<phi/degree<<endl;
      cout<<" nu="<<nu<<" Q2="<<Q2<<"  W="<<W<<"  xbj="<<xbj<<endl;
    }
#endif
    ////////////////////////////////////////////////////////////////
    //This part will be filled if mTreeLevel<2
    //I do not use p and theta because the resolution will sometimes make p>Beam
    //xs_1h,xs_3he,xs_4he,xs_12c,xs_14n,xs_27al,xs_ge180;
    if(mTreeLevel<2 && istop==0) {
      xs_1h  = GetXS(1,0,mBeam,p0,theta0,0.0,0.0);
      xs_3he = GetXS(2,1,mBeam,p0,theta0,0.0,0.0);
      xs_4he = GetXS(2,2,mBeam,p0,theta0,0.0,0.0);
      xs_12c = GetXS(6,6,mBeam,p0,theta0,0.0,0.0);
      xs_14n = GetXS(7,7,mBeam,p0,theta0,0.0,0.0);
      xs_27al = GetXS(13,14,mBeam,p0,theta0,0.0,0.0);
      xs_ge180 = GetXS_GE180(mBeam,p0,theta0,0.0,0.0);
      xs_4he_pi = WISER::GetXS(-211,2,2,mBeam,p0,theta0,0.067); //first have cross section per nucleon in mb/GeV/Sr
      xs_4he_pi *= (1000.0); //convert to ub/GeV/Sr
      xs_4he_pi *= 3.03; //convert to cross section per nuclei (A**0.8), A = 4
      //rate
      p_accept = (30.0+abs(-15.0))/100.0; 
      th_accept = (55.0+abs(-55.0))/1000.0;
      ph_accept = (50.0+abs(-50.0))/1000.0;
      n_trials=100000.0;
      tar_len_4he=10.0; // in cm
      tar_len_2h=10.0;
      p_spec=mDetMom; // in GeV
      dens_4he=2.41e28; // in number per m^3, using a density of 0.16 g/cm^3
      dens_2h=5.03e28;
      p_rate_2h=GetXS(1,1,mBeam,p0,theta0,0.0,0.0)*p_spec*p_accept*th_accept*ph_accept*dens_2h*1e-34*40e-6*tar_len_2h/100/(1.6*1e-19*n_trials);
      p_rate_4he=GetXS(2,2,mBeam,p0,theta0,0.0,0.0)*p_spec*p_accept*th_accept*ph_accept*dens_4he*1e-34*60e-6*tar_len_4he/100.0/(1.6*1e-19*n_trials);
      p_rate_4he_pi=xs_4he_pi*p_spec*p_accept*th_accept*ph_accept*dens_4he*1e-34*60e-6*tar_len_4he/100.0/(1.6*1e-19*n_trials);
      p_yield_2h=p_rate_2h/(40e-6);
      p_yield_4he=p_rate_4he/(60e-6);
      p_yield_4he_pi=p_rate_4he_pi/(60e-6);
    }

    //fill the tree
    mTree->Fill(); nentries_save++;
  }
  cout<<" number of events saved: "<<nentries_save<<endl;
  EndOfRun();
  return;
}

/////////////////////////////////////////////////////////////////////////////////////
void SHMSXSTree::EndOfRun()
{
#ifdef SHMSXSTree_Debug 
  if(SHMSXSTree_Debug>=4) cout<<"SHMSXSTree::EndOfStep() "<<endl; 
#endif

  mTree->Write();
  mTree->Delete();
  mFile->Write("",TObject::kOverwrite);
  mFile->Close();
  mFile->Delete();
}

/////////////////////////////////////////////////////////////////////////////////////
void SHMSXSTree::Reset()
{
  //ini
  vx0=vy0=vz0=p0=theta0=phi0=-9999.0;
  xtg0=ytg0=ztg0=xptg0=yptg0=delta0=-9999.0;

  //intermedia
  ixs=iys=-9999; 
  xsieve=ysieve=-9999.0;
  xfp=yfp=xpfp=ypfp=-9999.0;
  istop=-9999;

  //recon
  vx=vy=vz=p=theta=phi=-9999.0;
  xtg=ytg=ztg=xptg=yptg=delta=-9999.0;

  //kin
  nu=Q2=W=xbj=-9999.0;

  //this block will not be filled if mTreeLevel>=2
  xs_1h=xs_3he=xs_4he=xs_12c=xs_14n=xs_27al=xs_ge180=xs_4he_pi=p_rate_4he=p_yield_4he=p_rate_2h=p_yield_2h=-9999.0; 
  
}

//I put this wrapper here to avoid checking Q2 and isnan each call
double SHMSXSTree::GetXS(int Z, int N, float Ei, float Ef, float Theta, float Tb, float Ta)
{
  double pXs=0.0;
  if(mElasOnly) {
    pXs = ElasModel::GetXS(Z, N, (double)Ei, (double)Theta);
    if(isnanf(pXs)) {
      pXs=0.0;
#ifdef SHMSXSTree_Debug
      cout<<" isnan  from ElasModel::GetXS(Z="<<Z<<", N="<<N<<", Ei="<<Ei<<", Theta="<<Theta<<")\n";
#endif
    }
    pXs *= 1.0;   //already in ub/Sr
  } else {
    //note that XEM::GetXS() will not work for Q2>11, I give up these points
    if(Q2 < 11.0) pXs = XEM::GetXS(Z, N, (double)Ei, (double)Ef, (double)Theta, (double)Tb, (double)Ta);
    if(isnanf(pXs)) {
      pXs=0.0;
#ifdef SHMSXSTree_Debug
      cout<<" isnan  from XEM::GetXS(Z="<<Z<<", N="<<N<<", Ei="<<Ei<<", Ef="<<Ef<<", Theta="<<Theta<<", Tb="<<Tb<<", Ta="<<Ta<<")\n";
      //char cc[100];cout<<"\nPress any key to continue ...";cin>>cc;
#endif
    }
    pXs *= 1000.0;   //turn from ub/MeV/Sr into ub/GeV/Sr
  }
  return pXs;
}

//XS for GE180 compound
double SHMSXSTree::GetXS_GE180(float Ei, float Ef, float Theta, float Tb, float Ta)
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
    

  //comput xs using average Z and N
  double pXs=0.0;
  pXs = GetXS(pZ_aver_ge180, pN_aver_ge180, Ei, Ef, Theta, Tb, Ta);
    
  //comput xs for each element
  double pXs_O  = GetXS( 8,  8, Ei, Ef, Theta, Tb, Ta);
  double pXs_Si = GetXS(14, 14, Ei, Ef, Theta, Tb, Ta);
  double pXs_Ba = GetXS(56, 82, Ei, Ef, Theta, Tb, Ta);
  double pXs_Al = GetXS(13, 14, Ei, Ef, Theta, Tb, Ta);
  double pXs_Ca = GetXS(20, 20, Ei, Ef, Theta, Tb, Ta);
  double pXs_Sr = GetXS(38, 50, Ei, Ef, Theta, Tb, Ta);
  double pXs_H  = GetXS( 1,  0, Ei, Ef, Theta, Tb, Ta);

  const double kGE180_MXs[] = {pXs_Si+2*pXs_O, pXs_Ba+pXs_O, 2*pXs_Al+3*pXs_O,pXs_Ca+pXs_O, pXs_Sr+pXs_O, pXs_H};
    
  double pXs_ge180=0.0;
  for(int i=0;i<6;i++) {
    pXs_ge180 += pMolMass_aver_ge180*kGE180_MFraction[i]/kGE180_MMolMass[i]*kGE180_MXs[i];
  }

#ifdef SHMSXSTree_Debug 
  if(SHMSXSTree_Debug>=5) {
    cout<<" using average Z and N, XS_GE180 = "<<pXs
        <<"  sum up all elements, XS_GE180 = "<<pXs_ge180<<endl;
  }
#endif
  return pXs_ge180;
}

