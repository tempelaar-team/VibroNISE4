#ifndef ParmsMod
#define ParmsMod

#include "GlobalsMod.h"

#define ClubTotal 8

static const RealType SqrtHalfOfThree = (RealType) 1.22474487;
static const RealType Pi = (RealType) 3.14159265;
static const RealType LightSpeed = (RealType) 2.99792458e-5;							// [cm/fs]
static const RealType InvCmToFreq = (RealType) (2 * 3.14159265 * 2.99792458e-5);		// 2Pi c [cm/fs]
static const RealType MeVToFreq = (RealType) (1 / 6.58211951e+2);						// [/meV /fs]

enum ModeList { ModeNothing, ModeDebug, ModeDirect1D, ModeFTMapAnalysis, ModeResponse1D, ModeResponse2D, ModeDynamics };
enum AbsModeList { AbsModeDefault, AbsModeTim, AbsModeNick };
enum DDCouplingsFormatList { DDCouplingsFormatTim, DDCouplingsFormatNick, DDCouplingsFormatOld };
enum EnergyUnitList { WaveNumbers, MilliElectronVolt, Frequency, UnitLess };

typedef struct {
	
	enum ModeList Mode;
	enum AbsModeList AbsMode;
	enum DDCouplingsFormatList DDCouplingsFormat;
	enum EnergyUnitList EnergyUnit;
	
	int DAndC, K0Approach, Adiab, NoSE;
	int NoCTCTCoupling, Gaussian, Normalize, PrepareAdiab, TT_nOnly;
	int RandomSeed[2];
	int SampleTotal;
	int TimeTotal[2];
	int GBasisBound[ClubTotal], SBasisBound[ClubTotal], DBasisBound[ClubTotal], K0BasisBound[ClubTotal];
	int ACellTotal, BCellTotal, CellTotal, SubTotal, SiteTotal;
	int VTotal;
	int Periodic;
	int PolResolved, DoubleCrossPol;
	int SpectrumTotal;
	int ModeParm;

	RealType EnergyToFreq;
	RealType TimeStep;
	RealType AConstant, BConstant, ABAngle;
	RealType SEnergyMean, TTEnergyMean, TT_nEnergyMean;
	RealType DDCouplingNN, TripletDDCouplingNN, LLCouplingNN, HHCouplingNN, HLCouplingNN, LLHLCouplingNN, LLHHCouplingNN, LnLnHHCouplingNN;
	RealType TT_nDipole, RescaleEA;
	RealType HuangRhys[6];
	RealType VSlope, VRadius, EHRadius, TTRadius;
	RealType EHSeparated, EHScaling;
	RealType VEnergy;
	RealType BathReorganization;
    RealType LLSigma, HHSigma, HLSigma;
	RealType HBroad[3];
    RealType StokesShift;
	RealType LifeTime[3];
    RealType Temperature;
	RealType OptDielectric;
	RealType EnergyBounds[2], LaserBounds[2];
	RealType Polarizations[8];

	char DipolesFile[100], BathFile[100];
	char DDCouplingsFile[100], TripletDDCouplingsFile[100];
	char LLCouplingsFile[100], HHCouplingsFile[100], HLCouplingsFile[100], LLHLCouplingsFile[100], LLHHCouplingsFile[100], LnLnHHCouplingsFile[100];
	char EHEnergiesFile[100];

} ParmStruct;

#endif