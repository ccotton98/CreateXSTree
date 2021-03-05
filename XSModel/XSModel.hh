// This is all available XSModel
//I simply copy all headers here, one can include each individual header
//or include this one for all

#ifndef _ElasModel_H
#define _ElasModel_H

//return elas xs in micro barn
namespace  ElasModel
{
	//input
	//		Z,N: proton and neutron number of the nucleus.	;
	//		Ei: incoming electron energy in GeV;
	//		Theta: scattering angle for outgoing particle in radian;
	//		Mtg_GeV: target mass if it is given
	//		iUseLarry: in case for C12, use larry's model
	double GetXS(int Z, int N, double Ei, double Theta, double Mtg_GeV=0.0, int iUseLarry=0);

	double XS_H(double Ei, double theta);
	double XS_He3(double Ei, double theta);
	double XS_He4(double Ei, double theta);
	double XS_Ta184(double Ei, double theta);
	//wrapper for Stansfield's elas C12 XS, PRC3.1448 for details
	double XS_C12(double Ei, double theta);
	//wrapper for L. Cardman's elas C12 XS, very slow, use with caution
	//read Physics Letters V91B number 2, pp203-206 for details
	double Carbon(double Ei, double theta);
	// Stansfield's parameterized model for all nuclei, not good for high z
	//read PRC3.1448 for details
	double GetXS_All(double Ei, double theta, int iZ, double Mtg_GeV);
}
#endif


#ifndef _PBosted_H
#define _PBosted_H

namespace PBosted 
{
	//input
	//       Z,N: proton and neutron number of the nucleus.	;
	//Ei, Ef: incoming and outgoing electron energy in GeV;
	//Theta: scattering angle for outgoing particle in radian;
	//Tb and Ta will be used for radiated XS only if they are both positive
	//Tb: material thickness in unit of radiation length before scattering;
	//Ta: material thickness in unit of radiation length after scattering;
	double GetXS(int Z, int N, double Ei, double Ef, double theta, double Tb=-0.001, double Ta=-0.001);
}
#endif

#ifndef _XEM_H
#define _XEM_H

namespace XEM 
{
	//input
	//       Z,N: proton and neutron number of the nucleus.	;
	//Ei, Ef: incoming and outgoing electron energy in GeV;
	//Theta: scattering angle for outgoing particle in radian;
	//Tb and Ta will be used for radiated XS only if they are both positive
	//Tb: material thickness in unit of radiation length before scattering;
	//Ta: material thickness in unit of radiation length after scattering;
	double GetXS(int Z, int N, double Ei, double Ef, double theta, double Tb=-0.001, double Ta=-0.001);
}
#endif

#ifndef _WISER_H
#define _WISER_H

namespace WISER 
{
	//input:
	//       Z,N: proton and neutron number of the nucleus.
	//       PART: praticle type, 1(pi+), 2(pi-), 3(k+), 4(k-),5(pr)		;
	//Ei, Ef: incoming and outgoing electron energy in GeV;
	//Theta: scattering angle for outgoing particle in radian;
	//RadLen: material thick in unit of radiation length;
	double WISER(int PART, int Z, int N, double Ei, double Ef, double theta, double radlen);


	//input
	//		Z,N: proton and neutron number of the nucleus.
	//		Pid: praticle ID, following the PDG definition
	//       	2212	for p;	+/-211	for pi+/-;  +/-321	for k+/-		;
	//		Ei, Ef: incoming and outgoing electron energy in GeV;
	//		Theta: scattering angle for outgoing particle in radian;
	//		RadLen: material thick in unit of radiation length;
	double GetXS(int Pid, int Z, int N,  double Ei, double Ef, double Theta, double RadLen);
} 

#endif


#ifndef _QFS_N_EPC_H
#define _QFS_N_EPC_H

//This function is provided to calculate quasi-elastic and resonance cross section for different kinds of particles
//Use EPC model for proton, neutron, pions
//Use QFS model for electron
//Do not support positron any more
//Use a subroutine to calculate photons
//The photon cross section is calculated from the pi0 decay so it could not be used with pi0 simultaneously
//
//Usage: There are several forms of this getQElXS function
//getQElXS(PID,Z,N,Eb,theta,pf);
//getQElXS(PID,Z,N,Eb,theta,pf,EPS,EPSD,FP);
//getQElXS(PID,Z,N,Eb,theta,pf,Tb,Ta);
//getQElXS(PID,Z,N,Eb,theta,pf,EPS,EPSD,FP,Tb,Ta);
//
//Meaning of parameters:
//PID: praticle ID, following the PDG definition
//=	2212	for p		;	2112	for n	;	211	for pi+	;
//-211	for pi-		;	111		for pi0	;	11	for e-	;
//22		for photon	;
//Z,N: proton and neutron number of the nucleus.
//Eb: incoming electron energy in GeV;
//theta: scattering angle for outgoing particle in radian;
//pf: outgoing particle momentum in GeV/c;
//
//EPS,EPSD,FP: nucleus parameters, which used in the QFS code. EPS - seperation energy in MeV, 
//EPSD - delta seperation energy in MeV, FP - Fermi momentum in MeV/c; their values can be changed 
//to fit experimetal results, they will be set with the recommended value if not provided (EPS=10, EPSD=-10, FP=220);
//Tb,Ta: target parameters, Tb - total radiative length before scattering, Ta - total radiative 
//length after scattering, both in the unit of Radiation Length, if they are not provided, they will 
//be set with 0, which means not taking target radiative length into account;

namespace QFS_N_EPC
{
	double getQElXS(int PID, int Z, int N, double Eb, double theta, double pf);
	double getQElXS(int PID, int Z, int N, double Eb, double theta, double pf, 
		double EPS, double EPSD, double FP);
	double getQElXS(int PID, int Z, int N, double Eb, double theta, double pf, 
		double Tb, double Ta);
	double getQElXS(int PID, int Z, int N, double Eb, double theta, double pf, 
		double EPS, double EPSD, double FP, double Tb, double Ta);
}
#endif //_QFS_N_EPC_H



#ifndef _COMPTON_H
#define _COMPTON_H

namespace Compton
{
	//input: photon incident energy in GeV and 
	//       photon scattering angle in CM frame in radian 
	//return RCS XS in nb/GeV^2/sr
	double GetRCSXS(double Ei, double Theta_cm_g);
	double GetRCSKLL(double Ei, double Theta_cm_g);
	double GetRCSALL(double Ei, double Theta_cm_g);
	  
	double GetTheta_cm(double Ei, double Theta_lab_g);

	double GetRCSXS_lab(double Ei, double Theta_lab_g);
	double GetRCSKLL_lab(double Ei, double Theta_lab_g);
	double GetRCSALL_lab(double Ei, double Theta_lab_g);

  //root fcn
  double XSVsE_FixedThetalab_deg(double *x, double *par);
  double XSVsE_FixedThetacm_deg(double *x, double *par);

  //x is cos(Theta_cm)
  double RCSXSVsCosTh_cm(double *x, double *par);
  double RCSKLLVsCosTh_cm(double *x, double *par);
  double RCSALLVsCosTh_cm(double *x, double *par);

  //x is Theta_cm
  double RCSXSVsTh_cm(double *x, double *par);
  double RCSKLLVsTh_cm(double *x, double *par);
  double RCSALLVsTh_cm(double *x, double *par);

}


#endif  //_COMPTON_H


