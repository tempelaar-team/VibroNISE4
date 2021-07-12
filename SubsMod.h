typedef struct {

	int RefCount, RefTotal, Count, Total;
	long int Timer;

} ProgressStruct;


typedef struct {
	
	int TwoTimesPosA, TwoTimesPosB;
	int CellA, CellB;
	int Sub;
	
} CrystalStruct;


typedef struct {

	int TwoTimesDPosA, TwoTimesDPosB; // TwoTimesDPosA/B describe the distance in lattice spacings from row entry to column entry (column - row).
	int DCellA, DCellB;
	
} DCrystalStruct;


typedef struct {
	
	RealType *EHEnergies;
	MaskRArrayStruct DD, TripletDD, LL, HH, HL, LLHL, LLHH, LnLnHH; // All couplings describe transferring a particle from row entry to column entry.

} CouplingsStruct;


typedef struct {

	int Pulses[4];
	int Active[3];
	int Id;
	RealType Weight;

} PolStruct;


void InitProgress(ProgressStruct *Progress, const PolStruct *Pols, const ParmStruct *Parms);
void UpdateProgress(ProgressStruct *Progress);
void BuildCrystal(CrystalStruct **Crystal, const ParmStruct *Parms);
void BuildDCrystal(DCrystalStruct **DCrystal, RealType **Distances, const CrystalStruct *Crystal, const ParmStruct *Parms);
void BuildTransDistances(int **TransDistances, const CrystalStruct *Crystal, const ParmStruct *Parms);
void BuildPermDipoles(RealType **PermDipoles, PolStruct *Pols, const CrystalStruct *Crystal, const ParmStruct *Parms);
void BuildDDCouplings(MaskRArrayStruct *DDCouplings, const CrystalStruct *Crystal, const DCrystalStruct *DCrystal, const ParmStruct *Parms);
void BuildCouplings(MaskRArrayStruct *Couplings, const char *FileName, const RealType CouplingNN, const CrystalStruct *Crystal, const DCrystalStruct *DCrystal, const ParmStruct *Parms);
void BuildEHEnergies(RealType **EHEnergies, const DCrystalStruct *DCrystal, const RealType *Distances, ParmStruct *Parms);
void BuildK0Vectors(int **K0Vectors, const int *TransDistances, const CrystalStruct *Crystal, const BasisStruct *K0Basis, const BasisStruct *SBasis, const ParmStruct *Parms);
void BuildEnergyGrid(RealType **EnergyGrid, const ParmStruct *Parms);
void BuildDipoles(const char *FileName, MaskRArrayStruct *Dipoles, const BasisStruct *LoBasis, const BasisStruct *HiBasis, const int LoBasisTotal, const int HiBasisTotal, RealType * const *VOverlaps, const ParmStruct *Parms);
void BuildGEigenValues(RealType **GEigenValues, const BasisStruct *GBasis, const ParmStruct *Parms);
void BuildGPropagator(ComplexType **GPropagator, const BasisStruct *GBasis, const ParmStruct *Parms);
void BuildPropagator(ComplexType **Propagator, const RealType *EigenValues, const RealType TimeStepAsInvEnergy, const int BasisTotal);
void BuildHamiltonian(const char *FileName, MaskRArrayStruct *Hamiltonian, RealType **ToDiag, RealType **EigenValues, int *ToDiagTotal, const CouplingsStruct *Couplings, RealType *const *VOverlaps, const BasisStruct *Basis, const int *BasisBound, const ParmStruct *Parms);
void DrawOffDiagDisorder(MaskRArrayStruct *OffDiagDisorder, const MaskRArrayStruct Couplings, const RealType Sigma, const ParmStruct *Parms, const int Symmetric);
void AddOffDiagDisorder(MaskRArrayStruct Hamiltonian, const CouplingsStruct *OffDiagDisorder, RealType *const *VOverlaps, const BasisStruct *Basis, const int *BasisBound, const ParmStruct *Parms);
void BuildK0Hamiltonian(RealType **K0Hamiltonian, RealType **ToDiag, RealType **SEigenValues, int *ToDiagTotal, MaskRArrayStruct *SHamiltonian, const int *K0Vectors, const ParmStruct *Parms);
void BuildBathFT(ComplexType **BathFT, const RealType *EigenValues, const BasisStruct *Basis, const int BasisTotal, const int CutOff, const ParmStruct *Parms);
void BuildReperator(RealType **Reperator, const RealType *EigenVectors, const BasisStruct *Basis, const int BasisTotal, const int CutOff, const ParmStruct *Parms);
void BuildRedfieldTensor(MaskCArrayStruct *RedfieldTensor, const ComplexType *KetBathFT, const ComplexType *BraBathFT, const RealType *KetReperator, const RealType *BraReperator, const RealType *KetEigenValues, const RealType *BraEigenValues, const int KetBasisTotal, const int BraBasisTotal, const int KetCutOff, const int BraCutOff, const ParmStruct *Parms);
void BuildAdiabaticDipoles(RealType **ADipoles, const MaskRArrayStruct Dipoles, const RealType *PermDipoles, const RealType *LoEigenVectors, const RealType *HiEigenVectors, const RealType *LoEigenValues, const RealType *HiEigenValues, const int LoBasisTotal, const int HiBasisTotal, const PolStruct *Pols, const ParmStruct *Parms);
void K0ToSEigenVectors(RealType **SEigenVectors, const RealType *K0EigenVectors, const int *K0Vectors, const ParmStruct *Parms);
void SetPolarizations(PolStruct *Pols, const int DoubleCrossPol);
void CreatePermDipolesSeq(RealType *PermDipolesSeq, PolStruct *Pols, const RealType *PermDipoles, const ParmStruct *Parms);
void DirectAbsorption(RealType **Absorption, RealType **LineStrengths, const RealType *EnergyGrid, const RealType *PermDipoles, const RealType *SEigenVectors, const RealType *SEigenValues, const MaskRArrayStruct *GSDipoles, const PolStruct *Pols, const ParmStruct *Parms, const int StateTotal);
void DirectFluorescence(RealType **Fluorescence, RealType **FluStrengths, const RealType *EnergyGrid, const RealType *PermDipoles, const RealType *SEigenVectors, const RealType *SEigenValues, const MaskRArrayStruct *GSDipoles, const BasisStruct *GBasis, const PolStruct *Pols, const ParmStruct *Parms, const int StateTotal);
void DirectCPopulations(RealType **CPopulations,  const RealType *SEigenVectors, const BasisStruct *SBasis, const ParmStruct *Parms, const int StateTotal);
void DirectCPopulationsExtended(const RealType *SEigenVectors, const RealType *LineStrengths, const RealType *SEigenValues, const BasisStruct *SBasis, const ParmStruct *Parms, const int StateTotal);
void DirectTDensity(const RealType *SEigenVectors, const RealType *LineStrengths, const int *TransDistances, const DCrystalStruct *DCrystal, const BasisStruct *SBasis, const ParmStruct *Parms, const int StateTotal);
void DirectCTDensity(const RealType *SEigenVectors, const RealType *LineStrengths, const int *TransDistances, const DCrystalStruct *DCrystal, const BasisStruct *SBasis, const ParmStruct *Parms, const int StateTotal);
void DirectLineStrengths(RealType **LineStrengths, const RealType *SEigenValues, const RealType *GSAdiabDipoles, const ParmStruct *Parms, const int StateTotal);
void DirectFTMapAnalysis(PolStruct *Pols, const RealType *LineStrengths, const RealType *CPopulations, const RealType *GSAdiabDipoles, const RealType *SDAdiabDipoles, const RealType *SEigenValues, const RealType *DEigenValues, const BasisStruct *GBasis, const ParmStruct *Parms);
void Excite(ComplexType *HiCoef, const ComplexType *LoCoef, const RealType *LoHiDipoles, const int HiTotal, const int LoTotal);
void DeExcite(ComplexType *LoCoef, const ComplexType *HiCoef, const RealType *LoHiDipoles, const int LoTotal, const int HiTotal);
void ExciteVacuumDiabatic(ComplexType *Coef, const RealType *PermDipoles, const MaskRArrayStruct *LoHiDipoles, const int Total);
void ExciteDiabatic(ComplexType *HiCoef, const ComplexType *LoCoef, const RealType *PermDipoles, const MaskRArrayStruct *LoHiDipoles, const int HiTotal, const int LoTotal);
void DeExciteDiabatic(ComplexType *LoCoef, const ComplexType *HiCoef, const RealType *PermDipoles, const MaskRArrayStruct *LoHiDipoles, const int LoTotal, const int HiTotal);
void Propagate(ComplexType *Coef, ComplexType *CoefWork, const MaskRArrayStruct Hamiltonian, const int BasisTotal);
void Propagate_Redfield(ComplexType *Coef, ComplexType *CoefWork, const RealType *KetEigenValues, const RealType *BraEigenValues, const MaskCArrayStruct RedfieldTensor, const int KetBasisTotal, const int BraBasisTotal, const int KetCutOff, const int BraCutOff, const ParmStruct *Parms);
void ExtractSPopulations(RealType  *SPopulations, const RealType *Populations, const BasisStruct *SBasis, const ParmStruct *Parms);
void ExtractCPopulations(RealType  *CPopulations, const RealType *Populations, const BasisStruct *SBasis, const ParmStruct *Parms);
void SaveLineStrengths(const RealType *LineStrengths, const RealType *FluStrengths, const RealType *SEigenValues, const int StateTotal);
void SaveCPopulations(const RealType *CPopulations, const RealType *LineStrengths, const RealType *SEigenValues, const int StateTotal);
void SaveResponse2D(const char *FileName, const ComplexType *Response, const ParmStruct *Parms, const int MemSkip, const int Minus);

int DetermineStateCutOff(const RealType *EigenValues, const int BasisTotal, const ParmStruct *Parms);
ComplexType DeExciteVacuum(const ComplexType *Coef, const RealType *LoHiDipoles, const int Total);
ComplexType DeExciteVacuumDiabatic(const ComplexType *Coef, const RealType *PermDipoles, const MaskRArrayStruct *LoHiDipoles, const int Total);