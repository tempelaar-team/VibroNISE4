#ifndef BasisMod
#define BasisMod


#define PTotal 4
#define CTotal 6


enum Character { Ground, Singlet, Electron, Hole, Triplet, Triplet_n
};


typedef struct {
	
	int S;				// Site
	enum Character C;	// Character
	int V[CTotal];		// Vibrational states
	
} ParticleStruct;


typedef struct {
	
	ParticleStruct P[PTotal];
	
} BasisStruct;


void BuildBases(BasisStruct **GBasis, BasisStruct **SBasis, BasisStruct **DBasis, BasisStruct **K0Basis, int *GBasisBound, int *SBasisBound, int *DBasisBound, int *K0BasisBound, const RealType *Distances, const RealType VSlope, const RealType VRadius, const RealType EHRadius, const RealType TTRadius, const int SiteTotal, const int VTotal, const int CellSiteTotal, const int TT_nOnly);
void SwitchFirstTwoParticles(BasisStruct *BS);
void PrintBasisState(const BasisStruct BS);
void PrintBasisStateAdvanced(const BasisStruct BS);
void SaveOperator(const char *FileName, const RealType *Operator, const BasisStruct *PriBasis, const BasisStruct *SecBasis, const int PriDimension, const int SecDimension);

int Present(const BasisStruct *BS, const int C);
int CreateS(BasisStruct *BS, const int UnitId);
int AnnihilateS(BasisStruct *BS, const int UnitId);
int CreateE(BasisStruct *BS, const int UnitId);
int CreateH(BasisStruct *BS, const int UnitId);
int CreateT(BasisStruct *BS, const int UnitId);
int DipoleCreateT_n(BasisStruct *BS, const int UnitId);
int CreateT_n(BasisStruct *BS, const int UnitId);
int AnnihilateT_n(BasisStruct *BS, const int UnitId);
int CountParticles(const BasisStruct *BS);
int CountVibrations(const BasisStruct *BS);
int IdentifyBasisState(RealType *EffectiveVOverlap, const BasisStruct *BS, const BasisStruct *RefBS, RealType * const *VOverlaps, const int VTotal);
int IdentifyFirstTwoVibrations(const BasisStruct *BS, const BasisStruct *RefBS);
int IdentifyPurelyElectronic(const BasisStruct *BS, const BasisStruct *RefBS);


#endif