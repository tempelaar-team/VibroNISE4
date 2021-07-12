#include "GlobalsMod.h"
#include "ParmsMod.h"
#include "RandomMod.h"
#include "ToolsMod.h"
#include "BasisMod.h"
#include "VOverlapsMod.h"
#include "ParmsHandleMod.h"
#include "SubsMod.h"

#include <string.h>
#include <math.h>


static void Debug(void);
static void Direct1D(ParmStruct *Parms, const PolStruct *Pols, const RealType *PermDipoles, const CouplingsStruct *Couplings, RealType * const *VOverlaps, const CrystalStruct *Crystal, const DCrystalStruct *DCrystal, const RealType *Distances);
static void FTMapAnalysis(ParmStruct *Parms, PolStruct *Pols, const RealType *PermDipoles, const CouplingsStruct *Couplings, RealType *const *VOverlaps, const RealType *Distances);
static void Response1D(ParmStruct *Parms, PolStruct *Pols, const RealType *PermDipoles, const CouplingsStruct *Couplings, RealType * const *VOverlaps, const RealType *Distances);
static void Response1D_Adiab(ParmStruct *Parms, PolStruct *Pols, const RealType *PermDipoles, const CouplingsStruct *Couplings, RealType * const *VOverlaps, const RealType *Distances);
static void Response2D(ParmStruct *Parms, PolStruct *Pols, const RealType *PermDipoles, const CouplingsStruct *Couplings, RealType *const *VOverlaps, const RealType *Distances);
static void Response2D_Adiab(ParmStruct *Parms, PolStruct *Pols, const RealType *PermDipoles, const CouplingsStruct *Couplings, RealType *const *VOverlaps, const RealType *Distances);
static void Dynamics(ParmStruct *Parms, const CouplingsStruct *Couplings, RealType * const *VOverlaps, const RealType *Distances);


int main(int ArgsTotal, char *Args[])
{
	CouplingsStruct Couplings;

	RealType **VOverlaps;

	ParmStruct *Parms = calloc(1, sizeof(ParmStruct));
	PolStruct *Pols = calloc(1, sizeof(PolStruct));
	CrystalStruct *Crystal;
	DCrystalStruct *DCrystal;
	RealType *Distances, *PermDipoles;

	long int InitialTimer, FinalTimer;
	char Timing[10], DisplayNumber[10];
	
	
	if (ArgsTotal < 3) Embarrassing = 0;
	else sscanf(Args[2], "%d", &Embarrassing);
	
	InitializeLog(); PrintTitle();
	
	if (Embarrassing)
		IntToString(DisplayNumber, Embarrassing),
		PrintMessage("Running embarrassingly parallel session ", DisplayNumber, ".");

	InitializeParms(Parms);
	ReadParms(Args, Parms);
	DerivedParms(Parms);

	InitialTimer = SetTimer();
	PrintTime(); PrintMessage("Starting calculations!", "", "");
	
	
	PrintMessage("Constructing crystal parameters.", "", "");
	BuildCrystal(&Crystal, Parms);
	BuildDCrystal(&DCrystal, &Distances, Crystal, Parms);

	PrintMessage("Constructing overlap factors, couplings and permanent dipoles.", "", "");
	BuildVOverlaps(&VOverlaps, Parms->HuangRhys, CTotal, Parms->VTotal, Embarrassing);
	BuildPermDipoles(&PermDipoles, Pols, Crystal, Parms);
	BuildDDCouplings(&Couplings.DD, Crystal, DCrystal, Parms);
	BuildCouplings(&Couplings.TripletDD, Parms->TripletDDCouplingsFile, Parms->TripletDDCouplingNN, Crystal, DCrystal, Parms);
	BuildCouplings(&Couplings.LL, Parms->LLCouplingsFile, Parms->LLCouplingNN, Crystal, DCrystal, Parms);
	BuildCouplings(&Couplings.HH, Parms->HHCouplingsFile, Parms->HHCouplingNN, Crystal, DCrystal, Parms);
	BuildCouplings(&Couplings.HL, Parms->HLCouplingsFile, Parms->HLCouplingNN, Crystal, DCrystal, Parms);
	BuildCouplings(&Couplings.LLHL, Parms->LLHLCouplingsFile, Parms->LLHLCouplingNN, Crystal, DCrystal, Parms);
	BuildCouplings(&Couplings.LLHH, Parms->LLHHCouplingsFile, Parms->LLHHCouplingNN, Crystal, DCrystal, Parms);
	BuildCouplings(&Couplings.LnLnHH, Parms->LnLnHHCouplingsFile, Parms->LnLnHHCouplingNN, Crystal, DCrystal, Parms);
	BuildEHEnergies(&Couplings.EHEnergies, DCrystal, Distances, Parms);


	if (Parms->Mode == ModeDebug)				Debug();
	else if (Parms->Mode == ModeDirect1D)		Direct1D(Parms, Pols, PermDipoles, &Couplings, VOverlaps, Crystal, DCrystal, Distances);
	else if (Parms->Mode == ModeFTMapAnalysis)	FTMapAnalysis(Parms, Pols, PermDipoles, &Couplings, VOverlaps, Distances);
	else if (Parms->Mode == ModeResponse1D) {
		if (Parms->Adiab)						Response1D_Adiab(Parms, Pols, PermDipoles, &Couplings, VOverlaps, Distances);
		else									Response1D(Parms, Pols, PermDipoles, &Couplings, VOverlaps, Distances);
	}
	else if (Parms->Mode == ModeResponse2D) {
		if (Parms->Adiab)						Response2D_Adiab(Parms, Pols, PermDipoles, &Couplings, VOverlaps, Distances);
		else									Response2D(Parms, Pols, PermDipoles, &Couplings, VOverlaps, Distances);
	}
	else if (Parms->Mode == ModeDynamics)		Dynamics(Parms, &Couplings, VOverlaps, Distances);
	
	
	FinalTimer = SetTimer();
	
	IntToString(Timing, (int) (FinalTimer - InitialTimer));
	PrintTime(); PrintMessage("Finished in ", Timing, " seconds!");


	FreeMaskArray(&Couplings.DD); FreeMaskArray(&Couplings.TripletDD);
	FreeMaskArray(&Couplings.LL); FreeMaskArray(&Couplings.HH); FreeMaskArray(&Couplings.HL); FreeMaskArray(&Couplings.LLHL);
	FreeMaskArray(&Couplings.LLHH); FreeMaskArray(&Couplings.LnLnHH);
	free(Couplings.EHEnergies);

	FreeVOverlaps(&VOverlaps, CTotal);

	free(Parms); free(Pols); free(Crystal); free(DCrystal);
	free(Distances); free(PermDipoles);
	
	return 0;
}


static void Debug()
{

}


static void Direct1D(ParmStruct *Parms, const PolStruct *Pols, const RealType *PermDipoles, const CouplingsStruct *Couplings, RealType * const *VOverlaps, const CrystalStruct *Crystal, const DCrystalStruct *DCrystal, const RealType *Distances)
{
	MaskRArrayStruct SMatrix, GSDipoles;

	BasisStruct *GBasis, *SBasis, *K0Basis;

	RealType *EnergyGrid, *K0Matrix, *ToDiag, *SEigenValues, *Absorption, *Fluorescence, *LineStrengths, *FluStrengths, *CPopulations;

	int *TransDistances, *K0Vectors;

	int ToDiagTotal;

	BuildBases(&GBasis, &SBasis, NULL, &K0Basis, Parms->GBasisBound, Parms->SBasisBound, NULL, Parms->K0BasisBound, Distances, Parms->VSlope, Parms->VRadius, Parms->EHRadius, Parms->TTRadius, Parms->SiteTotal, Parms->VTotal, Parms->SubTotal, 0);

	BuildEnergyGrid(&EnergyGrid, Parms);
	BuildDipoles("GSDipoles", &GSDipoles, GBasis, SBasis, *Parms->GBasisBound, *Parms->SBasisBound, VOverlaps, Parms);
	BuildHamiltonian("SHamiltonian", &SMatrix, &ToDiag, &SEigenValues, &ToDiagTotal, Couplings, VOverlaps, SBasis, Parms->SBasisBound, Parms);

	BuildTransDistances(&TransDistances, Crystal, Parms);
	if (Parms->K0Approach) {
		BuildK0Vectors(&K0Vectors, TransDistances, Crystal, K0Basis, SBasis, Parms);
		BuildK0Hamiltonian(&K0Matrix, &ToDiag, &SEigenValues, &ToDiagTotal, &SMatrix, K0Vectors, Parms);
	} else K0Matrix = NULL, K0Vectors = NULL;

	InitializeDiag(ToDiag, SEigenValues, ToDiagTotal, Parms->DAndC);

	DiagonalizeMatrix(ToDiag, SEigenValues, ToDiagTotal, Parms->DAndC);

	if (Parms->K0Approach) K0ToSEigenVectors(&SMatrix.Array, K0Matrix, K0Vectors, Parms);

	DirectAbsorption(&Absorption, &LineStrengths, EnergyGrid, PermDipoles, SMatrix.Array, SEigenValues, &GSDipoles, Pols, Parms, ToDiagTotal);
    DirectFluorescence(&Fluorescence, &FluStrengths, EnergyGrid, PermDipoles, SMatrix.Array, SEigenValues, &GSDipoles, GBasis, Pols, Parms, ToDiagTotal);

	RealPlusRArray(EnergyGrid, Parms->SEnergyMean, Parms->SpectrumTotal);
	if (Parms->Normalize) {
        NormalizeRArray(Absorption, 4 * Parms->SpectrumTotal);
        NormalizeRArray(Fluorescence, 4 * Parms->SpectrumTotal);
    }
	SaveRArray("DirectAbs", Absorption, Parms->SpectrumTotal, 4, false, false, EnergyGrid, false);
    SaveRArray("DirectFlu", Fluorescence, Parms->SpectrumTotal, 4, false, false, EnergyGrid, false);

	RealPlusRArray(SEigenValues, Parms->SEnergyMean, ToDiagTotal);
	SaveLineStrengths(LineStrengths, FluStrengths, SEigenValues, ToDiagTotal);
	DirectCPopulations(&CPopulations, SMatrix.Array, SBasis, Parms, ToDiagTotal);
	SaveCPopulations(CPopulations, LineStrengths, SEigenValues, ToDiagTotal);
	DirectTDensity(SMatrix.Array, LineStrengths, TransDistances, DCrystal, SBasis, Parms, ToDiagTotal);
    DirectCTDensity(SMatrix.Array, LineStrengths, TransDistances, DCrystal, SBasis, Parms, ToDiagTotal);

	CloseDiag(Parms->DAndC);

	if (Parms->K0Approach) free(K0Matrix), free(K0Vectors);

	FreeMaskArray(&SMatrix); FreeMaskArray(&GSDipoles);

	free(GBasis); free(SBasis); free(K0Basis);
	free(EnergyGrid); free(SEigenValues); free(Absorption); free(Fluorescence); free(LineStrengths); free(FluStrengths); free(CPopulations);
	free(TransDistances);
}


static void FTMapAnalysis(ParmStruct *Parms, PolStruct *Pols, const RealType *PermDipoles, const CouplingsStruct *Couplings, RealType *const *VOverlaps, const RealType *Distances)
{
	MaskRArrayStruct SMatrix, DMatrix, GSDipoles, SDDipoles;

	BasisStruct *GBasis, *SBasis, *DBasis;

	RealType *GEigenValues, *SEigenValues, *DEigenValues, *GSAdiabDipoles, *SDAdiabDipoles, *LineStrengths, *CPopulations;


	BuildBases(&GBasis, &SBasis, &DBasis, NULL, Parms->GBasisBound, Parms->SBasisBound, Parms->DBasisBound, NULL, Distances, Parms->VSlope, Parms->VRadius, Parms->EHRadius, Parms->TTRadius, Parms->SiteTotal, Parms->VTotal, Parms->SubTotal, Parms->TT_nOnly);

	BuildDipoles("GSDipoles", &GSDipoles, GBasis, SBasis, *Parms->GBasisBound, *Parms->SBasisBound, VOverlaps, Parms);
	BuildDipoles("SDDipoles", &SDDipoles, SBasis, DBasis, *Parms->SBasisBound, *Parms->DBasisBound, VOverlaps, Parms);
	BuildHamiltonian("SHamiltonian", &SMatrix, NULL, &SEigenValues, NULL, Couplings, VOverlaps, SBasis, Parms->SBasisBound, Parms);
	BuildHamiltonian("DHamiltonian", &DMatrix, NULL, &DEigenValues, NULL, Couplings, VOverlaps, DBasis, Parms->DBasisBound, Parms);

	if (*Parms->DBasisBound > *Parms->SBasisBound)	InitializeDiag(DMatrix.Array, DEigenValues, *Parms->DBasisBound, Parms->DAndC);
	else											InitializeDiag(SMatrix.Array, SEigenValues, *Parms->SBasisBound, Parms->DAndC);

	DiagonalizeMatrix(SMatrix.Array, SEigenValues, *Parms->SBasisBound, Parms->DAndC);
	DiagonalizeMatrix(DMatrix.Array, DEigenValues, *Parms->DBasisBound, Parms->DAndC);

	BuildGEigenValues(&GEigenValues, GBasis, Parms);

	PrintTime(); PrintMessage("Building dipoles.", "", "");

	BuildAdiabaticDipoles(&GSAdiabDipoles, GSDipoles, PermDipoles, NULL, SMatrix.Array, GEigenValues, SEigenValues, *Parms->GBasisBound, *Parms->SBasisBound, Pols, Parms);
	BuildAdiabaticDipoles(&SDAdiabDipoles, SDDipoles, PermDipoles, SMatrix.Array, DMatrix.Array, SEigenValues, DEigenValues, *Parms->SBasisBound, *Parms->DBasisBound, Pols, Parms);

	RealPlusRArray(SEigenValues, Parms->SEnergyMean, *Parms->SBasisBound);
	RealPlusRArray(DEigenValues, 2 * Parms->SEnergyMean, *Parms->DBasisBound);
	DirectCPopulations(&CPopulations, SMatrix.Array, SBasis, Parms, *Parms->SBasisBound);
	DirectLineStrengths(&LineStrengths, SEigenValues, GSAdiabDipoles, Parms, *Parms->SBasisBound);
	DirectCPopulationsExtended(SMatrix.Array, LineStrengths, SEigenValues, SBasis, Parms, *Parms->SBasisBound);
	SaveLineStrengths(LineStrengths, NULL, SEigenValues, *Parms->SBasisBound);
	DirectFTMapAnalysis(Pols, LineStrengths, CPopulations, GSAdiabDipoles, SDAdiabDipoles, SEigenValues, DEigenValues, GBasis, Parms);


	CloseDiag(Parms->DAndC);

	FreeMaskArray(&SMatrix); FreeMaskArray(&DMatrix); FreeMaskArray(&GSDipoles); FreeMaskArray(&SDDipoles);

	free(GBasis); free(SBasis); free(DBasis);
	free(GEigenValues); free(SEigenValues); free(DEigenValues); free(GSAdiabDipoles); free(SDAdiabDipoles); free(LineStrengths); free(CPopulations);
}


static void Response1D(ParmStruct *Parms, PolStruct *Pols, const RealType *PermDipoles, const CouplingsStruct *Couplings, RealType * const *VOverlaps, const RealType *Distances)
{
	MaskRArrayStruct SMatrix, GSDipoles;

	BasisStruct *GBasis, *SBasis;

	RealType *TimeGrid = malloc(*Parms->TimeTotal * sizeof(RealType));

	ComplexType *SCoef, *SCoefWork;

	ComplexType *Response = calloc((size_t) 3 * *Parms->TimeTotal, sizeof(ComplexType));

	int TimeCount;

	ProgressStruct Progress;


	BuildBases(&GBasis, &SBasis, NULL, NULL, Parms->GBasisBound, Parms->SBasisBound, NULL, NULL, Distances, Parms->VSlope, Parms->VRadius, Parms->EHRadius, Parms->TTRadius, Parms->SiteTotal, Parms->VTotal, Parms->SubTotal, 0);

	BuildDipoles("GSDipoles", &GSDipoles, GBasis, SBasis, *Parms->GBasisBound, *Parms->SBasisBound, VOverlaps, Parms);
	BuildHamiltonian("SHamiltonian", &SMatrix, NULL, NULL, NULL, Couplings, VOverlaps, SBasis, Parms->SBasisBound, Parms);
	MirrorRMatrix(SMatrix.Array, *Parms->SBasisBound); MirrorIMatrix(SMatrix.Mask, *Parms->SBasisBound);
	RealTimesRArray(SMatrix.Array, Parms->EnergyToFreq * Parms->TimeStep, IntPow2(*Parms->SBasisBound));

	SCoef = calloc((size_t) *Parms->SBasisBound, sizeof(ComplexType)); SCoefWork = malloc(7 * *Parms->SBasisBound * sizeof(ComplexType));

	InitProgress(&Progress, Pols, Parms);
	Pols->Id = -1;
	while (true) {
		SetPolarizations(Pols, false);

		if (Pols->Id > 2) break;

		ExciteVacuumDiabatic(SCoef, &PermDipoles[Pols->Pulses[0] * 2 * Parms->SiteTotal], &GSDipoles, *Parms->SBasisBound);

		for (TimeCount = 0; TimeCount < *Parms->TimeTotal; TimeCount++) {
			UpdateProgress(&Progress);
			Response[3 * TimeCount + Pols->Id] = DeExciteVacuumDiabatic(SCoef, &PermDipoles[Pols->Pulses[1] * 2 * Parms->SiteTotal], &GSDipoles, *Parms->SBasisBound);

			Propagate(SCoef, SCoefWork, SMatrix, *Parms->SBasisBound);
		}
	}

	CreateGrid(TimeGrid, 0, (*Parms->TimeTotal - 1) * Parms->TimeStep, *Parms->TimeTotal);
	SaveCArray("Response1D", Response, *Parms->TimeTotal, 3, TimeGrid);


	FreeMaskArray(&SMatrix); FreeMaskArray(&GSDipoles);

	free(GBasis); free(SBasis);
	free(TimeGrid);
	free(SCoef); free(SCoefWork); free(Response);
}


static void Response1D_Adiab(ParmStruct *Parms, PolStruct *Pols, const RealType *PermDipoles, const CouplingsStruct *Couplings, RealType * const *VOverlaps, const RealType *Distances)
{
	MaskRArrayStruct SMatrix, GSDipoles;

	BasisStruct *GBasis, *SBasis;

	RealType *SEigenValues, *GSAdiabDipoles;

	RealType *TimeGrid = malloc(*Parms->TimeTotal * sizeof(RealType));

	ComplexType *SCoef, *SPropagator;

	ComplexType *Response = calloc((size_t) 3 * *Parms->TimeTotal, sizeof(ComplexType));

	int TimeCount;

	ProgressStruct Progress;


	BuildBases(&GBasis, &SBasis, NULL, NULL, Parms->GBasisBound, Parms->SBasisBound, NULL, NULL, Distances, Parms->VSlope, Parms->VRadius, Parms->EHRadius, Parms->TTRadius, Parms->SiteTotal, Parms->VTotal, Parms->SubTotal, 0);

	BuildDipoles("GSDipoles", &GSDipoles, GBasis, SBasis, *Parms->GBasisBound, *Parms->SBasisBound, VOverlaps, Parms);
	BuildHamiltonian("SHamiltonian", &SMatrix, NULL, &SEigenValues, NULL, Couplings, VOverlaps, SBasis, Parms->SBasisBound, Parms);
	InitializeDiag(SMatrix.Array, SEigenValues, *Parms->SBasisBound, Parms->DAndC);
	DiagonalizeMatrix(SMatrix.Array, SEigenValues, *Parms->SBasisBound, Parms->DAndC);
	BuildAdiabaticDipoles(&GSAdiabDipoles, GSDipoles, PermDipoles, NULL, SMatrix.Array, NULL, NULL, *Parms->GBasisBound, *Parms->SBasisBound, Pols, Parms);
	BuildPropagator(&SPropagator, SEigenValues, Parms->EnergyToFreq * Parms->TimeStep, *Parms->SBasisBound);

	SCoef = calloc((size_t) *Parms->SBasisBound, sizeof(ComplexType));

	InitProgress(&Progress, Pols, Parms);
	Pols->Id = -1;
	while (true) {
		SetPolarizations(Pols, false);

		if (Pols->Id > 2) break;

		RArrayToCArray(SCoef, &GSAdiabDipoles[Pols->Pulses[0] * *Parms->GBasisBound * *Parms->SBasisBound], *Parms->SBasisBound);

		for (TimeCount = 0; TimeCount < *Parms->TimeTotal; TimeCount++) {
			UpdateProgress(&Progress);
			Response[3 * TimeCount + Pols->Id] = DeExciteVacuum(SCoef, &GSAdiabDipoles[Pols->Pulses[1] * *Parms->GBasisBound * *Parms->SBasisBound], *Parms->SBasisBound);

			CArrayTimesCArray(SCoef, SPropagator, *Parms->SBasisBound);
		}
	}

	CreateGrid(TimeGrid, 0, (*Parms->TimeTotal - 1) * Parms->TimeStep, *Parms->TimeTotal);
	SaveCArray("Response1D", Response, *Parms->TimeTotal, 3, TimeGrid);


	CloseDiag(Parms->DAndC);

	FreeMaskArray(&SMatrix); FreeMaskArray(&GSDipoles);

	free(GBasis); free(SBasis);
	free(SEigenValues); free(GSAdiabDipoles); free(TimeGrid);
	free(SCoef); free(SPropagator); free(Response);
}


static void Response2D(ParmStruct *Parms, PolStruct *Pols, const RealType *PermDipoles, const CouplingsStruct *Couplings, RealType *const *VOverlaps, const RealType *Distances)
{
	typedef struct { ComplexType GBRR, GBNR, SERR, SENR, EARR, EANR; } ResponseStruct;

	ResponseStruct *Response = calloc((size_t) IntPow2(*Parms->TimeTotal), sizeof(ResponseStruct));

	MaskRArrayStruct SMatrix, DMatrix, GSDipoles, SDDipoles;

	BasisStruct *GBasis, *SBasis, *DBasis;

	RealType *PermDipolesSeq = malloc(8 * Parms->SiteTotal * sizeof(RealType));

	ComplexType *GCoefKet, *GCoefBra, *GCoefDeX, *SCoefKet, *SCoefBra, *SCoefDeX, *PriSCoefKet, *WaitSCoefKet, *WaitSCoefBra, *DCoefKet, *DCoefBra, *CoefWork, *GPropagator;

	int TimeCount[3];

	ProgressStruct Progress;


	BuildBases(&GBasis, &SBasis, &DBasis, NULL, Parms->GBasisBound, Parms->SBasisBound, Parms->DBasisBound, NULL, Distances, Parms->VSlope, Parms->VRadius, Parms->EHRadius, Parms->TTRadius, Parms->SiteTotal, Parms->VTotal, Parms->SubTotal, Parms->TT_nOnly);

	BuildDipoles("GSDipoles", &GSDipoles, GBasis, SBasis, *Parms->GBasisBound, *Parms->SBasisBound, VOverlaps, Parms);
	BuildDipoles("SDDipoles", &SDDipoles, SBasis, DBasis, *Parms->SBasisBound, *Parms->DBasisBound, VOverlaps, Parms);
	BuildHamiltonian("SHamiltonian", &SMatrix, NULL, NULL, NULL, Couplings, VOverlaps, SBasis, Parms->SBasisBound, Parms);
	MirrorRMatrix(SMatrix.Array, *Parms->SBasisBound); MirrorIMatrix(SMatrix.Mask, *Parms->SBasisBound);
	RealTimesRArray(SMatrix.Array, Parms->EnergyToFreq * Parms->TimeStep, IntPow2(*Parms->SBasisBound));
	BuildHamiltonian("DHamiltonian", &DMatrix, NULL, NULL, NULL, Couplings, VOverlaps, DBasis, Parms->DBasisBound, Parms);
	MirrorRMatrix(DMatrix.Array, *Parms->DBasisBound); MirrorIMatrix(DMatrix.Mask, *Parms->DBasisBound);
	RealTimesRArray(DMatrix.Array, Parms->EnergyToFreq * Parms->TimeStep, IntPow2(*Parms->DBasisBound));
	BuildGPropagator(&GPropagator, GBasis, Parms);

	GCoefKet = malloc(*Parms->GBasisBound * sizeof(ComplexType)); GCoefBra = malloc(*Parms->GBasisBound * sizeof(ComplexType));
	GCoefDeX = malloc(*Parms->GBasisBound * sizeof(ComplexType));
	SCoefKet = malloc(*Parms->SBasisBound * sizeof(ComplexType)); SCoefBra = malloc(*Parms->SBasisBound * sizeof(ComplexType));
	SCoefDeX = malloc(*Parms->SBasisBound * sizeof(ComplexType));
	PriSCoefKet = malloc(*Parms->SBasisBound * sizeof(ComplexType));
	WaitSCoefKet = malloc(*Parms->SBasisBound * sizeof(ComplexType)); WaitSCoefBra = malloc(*Parms->SBasisBound * sizeof(ComplexType));
	DCoefKet = (ComplexType *) malloc(*Parms->DBasisBound * sizeof(ComplexType)); DCoefBra = (ComplexType *) malloc(*Parms->DBasisBound * sizeof(ComplexType));
	CoefWork = malloc(7 * IntMax(*Parms->SBasisBound, *Parms->DBasisBound) * sizeof(ComplexType));


	PrintTime(); PrintMessage("Starting time evolution.", "", "");

	InitProgress(&Progress, Pols, Parms);
	Pols->Id = -1;
	while (true) {
		SetPolarizations(Pols, Parms->DoubleCrossPol);

		if (Pols->Id > 20) break;

		CreatePermDipolesSeq(PermDipolesSeq, Pols, PermDipoles, Parms);

		ExciteVacuumDiabatic(SCoefKet, PermDipolesSeq, &GSDipoles, *Parms->SBasisBound);

		for (TimeCount[0] = 0; TimeCount[0] < Parms->TimeTotal[0]; TimeCount[0]++) {
			UpdateProgress(&Progress);
			memcpy(PriSCoefKet, SCoefKet, *Parms->SBasisBound * sizeof(ComplexType));

			ExciteVacuumDiabatic(SCoefBra, &PermDipolesSeq[2 * Parms->SiteTotal], &GSDipoles, *Parms->SBasisBound);
			DeExciteDiabatic(GCoefKet, SCoefKet, &PermDipolesSeq[2 * Parms->SiteTotal], &GSDipoles, *Parms->GBasisBound, *Parms->SBasisBound);

			for (TimeCount[1] = 0; TimeCount[1] < Parms->TimeTotal[1]; TimeCount[1]++) {
				Propagate(SCoefKet, CoefWork, SMatrix, *Parms->SBasisBound);
				Propagate(SCoefBra, CoefWork, SMatrix, *Parms->SBasisBound);
				CArrayTimesCArray(GCoefKet, GPropagator, *Parms->GBasisBound);
			}

			memcpy(WaitSCoefKet, SCoefKet, *Parms->SBasisBound * sizeof(ComplexType)); memcpy(WaitSCoefBra, SCoefBra, *Parms->SBasisBound * sizeof(ComplexType));

			ExciteDiabatic(SCoefKet, GCoefKet, &PermDipolesSeq[4 * Parms->SiteTotal], &GSDipoles, *Parms->SBasisBound, *Parms->GBasisBound);
			ExciteVacuumDiabatic(SCoefBra, &PermDipolesSeq[4 * Parms->SiteTotal], &GSDipoles, *Parms->SBasisBound);

			for (TimeCount[2] = 0; TimeCount[2] < *Parms->TimeTotal; TimeCount[2]++) {
				DeExciteDiabatic(GCoefBra, SCoefBra, &PermDipolesSeq[6 * Parms->SiteTotal], &GSDipoles, *Parms->GBasisBound, *Parms->SBasisBound);
				Response[TimeCount[0] * *Parms->TimeTotal + TimeCount[2]].GBRR += Pols->Weight * CInnerProduct(GCoefKet, GCoefBra, *Parms->GBasisBound);

				Response[TimeCount[0] * *Parms->TimeTotal + TimeCount[2]].GBNR += Pols->Weight * DeExciteVacuumDiabatic(SCoefKet, &PermDipolesSeq[6 * Parms->SiteTotal], &GSDipoles, *Parms->SBasisBound);

				Propagate(SCoefKet, CoefWork, SMatrix, *Parms->SBasisBound);
				Propagate(SCoefBra, CoefWork, SMatrix, *Parms->SBasisBound);
				CArrayTimesCArray(GCoefKet, GPropagator, *Parms->GBasisBound);
			}

			memcpy(SCoefKet, WaitSCoefKet, *Parms->SBasisBound * sizeof(ComplexType)); memcpy(SCoefBra, WaitSCoefBra, *Parms->SBasisBound * sizeof(ComplexType));

			DeExciteDiabatic(GCoefKet, SCoefKet, &PermDipolesSeq[4 * Parms->SiteTotal], &GSDipoles, *Parms->GBasisBound, *Parms->SBasisBound);
			DeExciteDiabatic(GCoefBra, SCoefBra, &PermDipolesSeq[4 * Parms->SiteTotal], &GSDipoles, *Parms->GBasisBound, *Parms->SBasisBound);
			ExciteDiabatic(DCoefKet, SCoefKet, &PermDipolesSeq[4 * Parms->SiteTotal], &SDDipoles, *Parms->DBasisBound, *Parms->SBasisBound);
			ExciteDiabatic(DCoefBra, SCoefBra, &PermDipolesSeq[4 * Parms->SiteTotal], &SDDipoles, *Parms->DBasisBound, *Parms->SBasisBound);

			for (TimeCount[2] = 0; TimeCount[2] < *Parms->TimeTotal; TimeCount[2]++) {
				DeExciteDiabatic(GCoefDeX, SCoefBra, &PermDipolesSeq[6 * Parms->SiteTotal], &GSDipoles, *Parms->GBasisBound, *Parms->SBasisBound);
				Response[TimeCount[0] * *Parms->TimeTotal + TimeCount[2]].SERR += Pols->Weight * CInnerProduct(GCoefKet, GCoefDeX, *Parms->GBasisBound);

				DeExciteDiabatic(GCoefDeX, SCoefKet, &PermDipolesSeq[6 * Parms->SiteTotal], &GSDipoles, *Parms->GBasisBound, *Parms->SBasisBound);
				Response[TimeCount[0] * *Parms->TimeTotal + TimeCount[2]].SENR += Pols->Weight * CInnerProduct(GCoefBra, GCoefDeX, *Parms->GBasisBound);

				DeExciteDiabatic(SCoefDeX, DCoefBra, &PermDipolesSeq[6 * Parms->SiteTotal], &SDDipoles, *Parms->SBasisBound, *Parms->DBasisBound);
				Response[TimeCount[0] * *Parms->TimeTotal + TimeCount[2]].EARR += Pols->Weight * CInnerProduct(SCoefKet, SCoefDeX, *Parms->SBasisBound);

				DeExciteDiabatic(SCoefDeX, DCoefKet, &PermDipolesSeq[6 * Parms->SiteTotal], &SDDipoles, *Parms->SBasisBound, *Parms->DBasisBound);
				Response[TimeCount[0] * *Parms->TimeTotal + TimeCount[2]].EANR += Pols->Weight * CInnerProduct(SCoefBra, SCoefDeX, *Parms->SBasisBound);

				Propagate(DCoefKet, CoefWork, DMatrix, *Parms->DBasisBound);
				Propagate(DCoefBra, CoefWork, DMatrix, *Parms->DBasisBound);
				Propagate(SCoefKet, CoefWork, SMatrix, *Parms->SBasisBound);
				Propagate(SCoefBra, CoefWork, SMatrix, *Parms->SBasisBound);
				CArrayTimesCArray(GCoefKet, GPropagator, *Parms->GBasisBound);
				CArrayTimesCArray(GCoefBra, GPropagator, *Parms->GBasisBound);
			}

			memcpy(SCoefKet, PriSCoefKet, *Parms->SBasisBound * sizeof(ComplexType));
			Propagate(SCoefKet, CoefWork, SMatrix, *Parms->SBasisBound);
		}
	}

	PrintTime(); PrintMessage("Saving data.", "", "");

	SaveResponse2D("ResponseGBRR", &(*Response).GBRR, Parms, sizeof(ResponseStruct) / sizeof(Response[0].GBRR), false);
	SaveResponse2D("ResponseGBNR", &(*Response).GBNR, Parms, sizeof(ResponseStruct) / sizeof(Response[0].GBRR), false);
	SaveResponse2D("ResponseSERR", &(*Response).SERR, Parms, sizeof(ResponseStruct) / sizeof(Response[0].GBRR), false);
	SaveResponse2D("ResponseSENR", &(*Response).SENR, Parms, sizeof(ResponseStruct) / sizeof(Response[0].GBRR), false);
	SaveResponse2D("ResponseEARR", &(*Response).EARR, Parms, sizeof(ResponseStruct) / sizeof(Response[0].GBRR), true);
	SaveResponse2D("ResponseEANR", &(*Response).EANR, Parms, sizeof(ResponseStruct) / sizeof(Response[0].GBRR), true);


	FreeMaskArray(&SMatrix); FreeMaskArray(&DMatrix); FreeMaskArray(&GSDipoles); FreeMaskArray(&SDDipoles);

	free(Response);
	free(GBasis); free(SBasis); free(DBasis);
	free(PermDipolesSeq);
	free(GCoefKet); free(GCoefBra); free(GCoefDeX); free(SCoefKet); free(SCoefBra); free(SCoefDeX); free(PriSCoefKet); free(WaitSCoefKet); free(WaitSCoefBra);
	free(DCoefKet); free(DCoefBra); free(CoefWork); free(GPropagator);
}


static void Response2D_Adiab(ParmStruct *Parms, PolStruct *Pols, const RealType *PermDipoles, const CouplingsStruct *Couplings, RealType *const *VOverlaps, const RealType *Distances)
{
	typedef struct { ComplexType GBRR, GBNR, SERR, SENR, EARR, EANR; } ResponseStruct;

	ResponseStruct *Response = calloc((size_t) IntPow2(*Parms->TimeTotal), sizeof(ResponseStruct));

	MaskRArrayStruct SMatrix, DMatrix, GSDipoles, SDDipoles;

	BasisStruct *GBasis, *SBasis, *DBasis;

	RealType *GEigenValues, *SEigenValues, *DEigenValues, *GSAdiabDipoles, *SDAdiabDipoles;

	ComplexType *GCoefKet, *GCoefBra, *GCoefDeX, *SCoefKet, *SCoefBra, *SCoefDeX, *PriSCoefKet, *WaitSCoefKet, *WaitSCoefBra, *DCoefKet, *DCoefBra;
	ComplexType *GPropagator, *SPropagator, *DPropagator;

	int TimeCount[3];

	ProgressStruct Progress;


	BuildBases(&GBasis, &SBasis, &DBasis, NULL, Parms->GBasisBound, Parms->SBasisBound, Parms->DBasisBound, NULL, Distances, Parms->VSlope, Parms->VRadius, Parms->EHRadius, Parms->TTRadius, Parms->SiteTotal, Parms->VTotal, Parms->SubTotal, Parms->TT_nOnly);

	BuildDipoles("GSDipoles", &GSDipoles, GBasis, SBasis, *Parms->GBasisBound, *Parms->SBasisBound, VOverlaps, Parms);
	BuildDipoles("SDDipoles", &SDDipoles, SBasis, DBasis, *Parms->SBasisBound, *Parms->DBasisBound, VOverlaps, Parms);
	BuildHamiltonian("SHamiltonian", &SMatrix, NULL, &SEigenValues, NULL, Couplings, VOverlaps, SBasis, Parms->SBasisBound, Parms);
	BuildHamiltonian("DHamiltonian", &DMatrix, NULL, &DEigenValues, NULL, Couplings, VOverlaps, DBasis, Parms->DBasisBound, Parms);

	if (*Parms->DBasisBound > *Parms->SBasisBound)	InitializeDiag(DMatrix.Array, DEigenValues, *Parms->DBasisBound, Parms->DAndC);
	else											InitializeDiag(SMatrix.Array, SEigenValues, *Parms->SBasisBound, Parms->DAndC);

	DiagonalizeMatrix(SMatrix.Array, SEigenValues, *Parms->SBasisBound, Parms->DAndC);
	DiagonalizeMatrix(DMatrix.Array, DEigenValues, *Parms->DBasisBound, Parms->DAndC);

	PrintTime(); PrintMessage("Building propagators and dipoles.", "", "");

	BuildPropagator(&SPropagator, SEigenValues, Parms->EnergyToFreq * Parms->TimeStep, *Parms->SBasisBound);
	BuildPropagator(&DPropagator, DEigenValues, Parms->EnergyToFreq * Parms->TimeStep, *Parms->DBasisBound);

	BuildGEigenValues(&GEigenValues, GBasis, Parms);
	BuildPropagator(&GPropagator, GEigenValues, Parms->EnergyToFreq * Parms->TimeStep, *Parms->GBasisBound);

	BuildAdiabaticDipoles(&GSAdiabDipoles, GSDipoles, PermDipoles, NULL, SMatrix.Array, GEigenValues, SEigenValues, *Parms->GBasisBound, *Parms->SBasisBound, Pols, Parms);
	BuildAdiabaticDipoles(&SDAdiabDipoles, SDDipoles, PermDipoles, SMatrix.Array, DMatrix.Array, SEigenValues, DEigenValues, *Parms->SBasisBound, *Parms->DBasisBound, Pols, Parms);

	GCoefKet = malloc(*Parms->GBasisBound * sizeof(ComplexType)); GCoefBra = malloc(*Parms->GBasisBound * sizeof(ComplexType));
	GCoefDeX = malloc(*Parms->GBasisBound * sizeof(ComplexType));
	SCoefKet = malloc(*Parms->SBasisBound * sizeof(ComplexType)); SCoefBra = malloc(*Parms->SBasisBound * sizeof(ComplexType));
	SCoefDeX = malloc(*Parms->SBasisBound * sizeof(ComplexType));
	PriSCoefKet = malloc(*Parms->SBasisBound * sizeof(ComplexType));
	WaitSCoefKet = malloc(*Parms->SBasisBound * sizeof(ComplexType)); WaitSCoefBra = malloc(*Parms->SBasisBound * sizeof(ComplexType));
	DCoefKet = (ComplexType *) malloc(*Parms->DBasisBound * sizeof(ComplexType)); DCoefBra = (ComplexType *) malloc(*Parms->DBasisBound * sizeof(ComplexType));


	PrintTime(); PrintMessage("Starting time evolution.", "", "");

	InitProgress(&Progress, Pols, Parms);
	Pols->Id = -1;
	while (true) {
		SetPolarizations(Pols, Parms->DoubleCrossPol);

		if (Pols->Id > 20) break;

		RArrayToCArray(SCoefKet, &GSAdiabDipoles[Pols->Pulses[0] * *Parms->GBasisBound * *Parms->SBasisBound], *Parms->SBasisBound);

		for (TimeCount[0] = 0; TimeCount[0] < Parms->TimeTotal[0]; TimeCount[0]++) {
			UpdateProgress(&Progress);
			memcpy(PriSCoefKet, SCoefKet, *Parms->SBasisBound * sizeof(ComplexType));

			RArrayToCArray(SCoefBra, &GSAdiabDipoles[Pols->Pulses[1] * *Parms->GBasisBound * *Parms->SBasisBound], *Parms->SBasisBound);
			DeExcite(GCoefKet, SCoefKet, &GSAdiabDipoles[Pols->Pulses[1] * *Parms->GBasisBound * *Parms->SBasisBound], *Parms->GBasisBound, *Parms->SBasisBound);

			for (TimeCount[1] = 0; TimeCount[1] < Parms->TimeTotal[1]; TimeCount[1]++) {
				CArrayTimesCArray(SCoefKet, SPropagator, *Parms->SBasisBound);
				CArrayTimesCArray(SCoefBra, SPropagator, *Parms->SBasisBound);
				CArrayTimesCArray(GCoefKet, GPropagator, *Parms->GBasisBound);
			}

			memcpy(WaitSCoefKet, SCoefKet, *Parms->SBasisBound * sizeof(ComplexType)); memcpy(WaitSCoefBra, SCoefBra, *Parms->SBasisBound * sizeof(ComplexType));

			Excite(SCoefKet, GCoefKet, &GSAdiabDipoles[Pols->Pulses[2] * *Parms->GBasisBound * *Parms->SBasisBound], *Parms->SBasisBound, *Parms->GBasisBound);
			RArrayToCArray(SCoefBra, &GSAdiabDipoles[Pols->Pulses[2] * *Parms->GBasisBound * *Parms->SBasisBound], *Parms->SBasisBound);

			for (TimeCount[2] = 0; TimeCount[2] < *Parms->TimeTotal; TimeCount[2]++) {
				DeExcite(GCoefBra, SCoefBra, &GSAdiabDipoles[Pols->Pulses[3] * *Parms->GBasisBound * *Parms->SBasisBound], *Parms->GBasisBound, *Parms->SBasisBound);
				Response[TimeCount[0] * *Parms->TimeTotal + TimeCount[2]].GBRR += Pols->Weight * CInnerProduct(GCoefKet, GCoefBra, *Parms->GBasisBound);

				Response[TimeCount[0] * *Parms->TimeTotal + TimeCount[2]].GBNR += Pols->Weight * DeExciteVacuum(SCoefKet, &GSAdiabDipoles[Pols->Pulses[3] * *Parms->GBasisBound * *Parms->SBasisBound], *Parms->SBasisBound);

				CArrayTimesCArray(SCoefKet, SPropagator, *Parms->SBasisBound);
				CArrayTimesCArray(SCoefBra, SPropagator, *Parms->SBasisBound);
				CArrayTimesCArray(GCoefKet, GPropagator, *Parms->GBasisBound);
			}

			memcpy(SCoefKet, WaitSCoefKet, *Parms->SBasisBound * sizeof(ComplexType)); memcpy(SCoefBra, WaitSCoefBra, *Parms->SBasisBound * sizeof(ComplexType));

			DeExcite(GCoefKet, SCoefKet, &GSAdiabDipoles[Pols->Pulses[2] * *Parms->GBasisBound * *Parms->SBasisBound], *Parms->GBasisBound, *Parms->SBasisBound);
			DeExcite(GCoefBra, SCoefBra, &GSAdiabDipoles[Pols->Pulses[2] * *Parms->GBasisBound * *Parms->SBasisBound], *Parms->GBasisBound, *Parms->SBasisBound);
			Excite(DCoefKet, SCoefKet, &SDAdiabDipoles[Pols->Pulses[2] * *Parms->SBasisBound * *Parms->DBasisBound], *Parms->DBasisBound, *Parms->SBasisBound);
			Excite(DCoefBra, SCoefBra, &SDAdiabDipoles[Pols->Pulses[2] * *Parms->SBasisBound * *Parms->DBasisBound], *Parms->DBasisBound, *Parms->SBasisBound);

			for (TimeCount[2] = 0; TimeCount[2] < *Parms->TimeTotal; TimeCount[2]++) {
				DeExcite(GCoefDeX, SCoefBra, &GSAdiabDipoles[Pols->Pulses[3] * *Parms->GBasisBound * *Parms->SBasisBound], *Parms->GBasisBound, *Parms->SBasisBound);
				Response[TimeCount[0] * *Parms->TimeTotal + TimeCount[2]].SERR += Pols->Weight * CInnerProduct(GCoefKet, GCoefDeX, *Parms->GBasisBound);

				DeExcite(GCoefDeX, SCoefKet, &GSAdiabDipoles[Pols->Pulses[3] * *Parms->GBasisBound * *Parms->SBasisBound], *Parms->GBasisBound, *Parms->SBasisBound);
				Response[TimeCount[0] * *Parms->TimeTotal + TimeCount[2]].SENR += Pols->Weight * CInnerProduct(GCoefBra, GCoefDeX, *Parms->GBasisBound);

				DeExcite(SCoefDeX, DCoefBra, &SDAdiabDipoles[Pols->Pulses[3] * *Parms->SBasisBound * *Parms->DBasisBound], *Parms->SBasisBound, *Parms->DBasisBound);
				Response[TimeCount[0] * *Parms->TimeTotal + TimeCount[2]].EARR += Pols->Weight * CInnerProduct(SCoefKet, SCoefDeX, *Parms->SBasisBound);

				DeExcite(SCoefDeX, DCoefKet, &SDAdiabDipoles[Pols->Pulses[3] * *Parms->SBasisBound * *Parms->DBasisBound], *Parms->SBasisBound, *Parms->DBasisBound);
				Response[TimeCount[0] * *Parms->TimeTotal + TimeCount[2]].EANR += Pols->Weight * CInnerProduct(SCoefBra, SCoefDeX, *Parms->SBasisBound);

				CArrayTimesCArray(DCoefKet, DPropagator, *Parms->DBasisBound);
				CArrayTimesCArray(DCoefBra, DPropagator, *Parms->DBasisBound);
				CArrayTimesCArray(SCoefKet, SPropagator, *Parms->SBasisBound);
				CArrayTimesCArray(SCoefBra, SPropagator, *Parms->SBasisBound);
				CArrayTimesCArray(GCoefKet, GPropagator, *Parms->GBasisBound);
				CArrayTimesCArray(GCoefBra, GPropagator, *Parms->GBasisBound);
			}

			memcpy(SCoefKet, PriSCoefKet, *Parms->SBasisBound * sizeof(ComplexType));
			CArrayTimesCArray(SCoefKet, SPropagator, *Parms->SBasisBound);
		}
	}

	PrintTime(); PrintMessage("Saving data.", "", "");

	SaveResponse2D("ResponseGBRR", &(*Response).GBRR, Parms, sizeof(ResponseStruct) / sizeof(Response[0].GBRR), false);
	SaveResponse2D("ResponseGBNR", &(*Response).GBNR, Parms, sizeof(ResponseStruct) / sizeof(Response[0].GBRR), false);
	SaveResponse2D("ResponseSERR", &(*Response).SERR, Parms, sizeof(ResponseStruct) / sizeof(Response[0].GBRR), false);
	SaveResponse2D("ResponseSENR", &(*Response).SENR, Parms, sizeof(ResponseStruct) / sizeof(Response[0].GBRR), false);
	SaveResponse2D("ResponseEARR", &(*Response).EARR, Parms, sizeof(ResponseStruct) / sizeof(Response[0].GBRR), true);
	SaveResponse2D("ResponseEANR", &(*Response).EANR, Parms, sizeof(ResponseStruct) / sizeof(Response[0].GBRR), true);

	CloseDiag(Parms->DAndC);

	FreeMaskArray(&SMatrix); FreeMaskArray(&DMatrix); FreeMaskArray(&GSDipoles); FreeMaskArray(&SDDipoles);

	free(Response);
	free(GBasis); free(SBasis); free(DBasis);
	free(GEigenValues); free(SEigenValues); free(DEigenValues); free(GSAdiabDipoles); free(SDAdiabDipoles);
	free(GCoefKet); free(GCoefBra); free(GCoefDeX); free(SCoefKet); free(SCoefBra); free(SCoefDeX); free(PriSCoefKet); free(WaitSCoefKet); free(WaitSCoefBra);
	free(DCoefKet); free(DCoefBra); free(GPropagator); free(SPropagator); free(DPropagator);
}


static void Dynamics(ParmStruct *Parms, const CouplingsStruct *Couplings, RealType * const *VOverlaps, const RealType *Distances)
{
	MaskRArrayStruct StaticHamiltonian, SMatrix;
	MaskCArrayStruct SSRedfieldTensor;
    CouplingsStruct OffDiagDisorder;

	BasisStruct *SBasis;

	RealType *SReperator, *SEigenVectorsT, *SEigenValues, *Populations, *PopulationsAdiab, *SPopulations, *PopulationsElement, *PopulationsAdiabElement;

	RealType *TimeGrid = malloc(*Parms->TimeTotal * sizeof(RealType)), *CPopulations = calloc((size_t) 3 * *Parms->TimeTotal, sizeof(RealType));

	ComplexType *SBathFT, *SSCoefDia, *SSCoef, *SSCoefWork;

    RealType MaxPop;

	int SampleCount, TimeCount, PriSCount, SecSCount, CutOff, MaxPopId;

	ProgressStruct Progress;


	BuildBases(NULL, &SBasis, NULL, NULL, NULL, Parms->SBasisBound, NULL, NULL, Distances, Parms->VSlope, Parms->VRadius, Parms->EHRadius, Parms->TTRadius, Parms->SiteTotal, Parms->VTotal, Parms->SubTotal, 0);
	BuildHamiltonian("SHamiltonian", &StaticHamiltonian, NULL, &SEigenValues, NULL, Couplings, VOverlaps, SBasis, Parms->SBasisBound, Parms);

	SEigenVectorsT = calloc((size_t) IntPow2(*Parms->SBasisBound), sizeof(RealType));
	SSCoefDia = calloc((size_t) IntPow2(*Parms->SBasisBound), sizeof(ComplexType));
	SSCoef = malloc(IntPow2(*Parms->SBasisBound) * sizeof(ComplexType));
	SSCoefWork = malloc(7 * IntPow2(*Parms->SBasisBound) * sizeof(ComplexType));
	Populations = calloc((size_t) *Parms->TimeTotal * *Parms->SBasisBound, sizeof(RealType));
    PopulationsAdiab = calloc((size_t) *Parms->TimeTotal * *Parms->SBasisBound, sizeof(RealType));
    SPopulations = calloc((size_t) *Parms->TimeTotal * Parms->SiteTotal, sizeof(RealType));

    AllocateMaskArray(&OffDiagDisorder.LL, IntPow2(Parms->SiteTotal), true, false);
    AllocateMaskArray(&OffDiagDisorder.HH, IntPow2(Parms->SiteTotal), true, false);
    AllocateMaskArray(&OffDiagDisorder.HL, IntPow2(Parms->SiteTotal), true, false);
    memcpy(OffDiagDisorder.LL.Mask, (*Couplings).LL.Mask, IntPow2(Parms->SiteTotal) * sizeof(int));
    memcpy(OffDiagDisorder.HH.Mask, (*Couplings).HH.Mask, IntPow2(Parms->SiteTotal) * sizeof(int));
    memcpy(OffDiagDisorder.HL.Mask, (*Couplings).HL.Mask, IntPow2(Parms->SiteTotal) * sizeof(int));


    AllocateMaskArray(&SMatrix, IntPow2(*Parms->SBasisBound), true, true);

    RandomInitialise(Parms->RandomSeed[0], Parms->RandomSeed[1] + Embarrassing);
	InitializeDiag(SMatrix.Array, SEigenValues, *Parms->SBasisBound, Parms->DAndC);

    PrintTime(); PrintMessage("Starting time evolution.", "", "");
    InitProgress(&Progress, NULL, Parms);

    for (SampleCount = 0; SampleCount < Parms->SampleTotal; SampleCount++) {

        DrawOffDiagDisorder(&OffDiagDisorder.LL, (*Couplings).LL, Parms->LLSigma, Parms, false);
        DrawOffDiagDisorder(&OffDiagDisorder.HH, (*Couplings).HH, Parms->HHSigma, Parms, false);
        DrawOffDiagDisorder(&OffDiagDisorder.HL, (*Couplings).HL, Parms->HLSigma, Parms, true);

        memcpy(SMatrix.Array, StaticHamiltonian.Array, IntPow2(*Parms->SBasisBound) * sizeof(RealType));
        AddOffDiagDisorder(SMatrix, &OffDiagDisorder, VOverlaps, SBasis, Parms->SBasisBound, Parms);
        DiagonalizeMatrix(SMatrix.Array, SEigenValues, *Parms->SBasisBound, Parms->DAndC);
        TransposeRMatrix(SEigenVectorsT, SMatrix.Array, *Parms->SBasisBound);

        CutOff = DetermineStateCutOff(SEigenValues, *Parms->SBasisBound, Parms);

        BuildBathFT(&SBathFT, SEigenValues, SBasis, *Parms->SBasisBound, CutOff, Parms);
        BuildReperator(&SReperator, SMatrix.Array, SBasis, *Parms->SBasisBound, CutOff, Parms);
        BuildRedfieldTensor(&SSRedfieldTensor, SBathFT, SBathFT, SReperator, SReperator, SEigenValues, SEigenValues, *Parms->SBasisBound, *Parms->SBasisBound, CutOff, CutOff, Parms);

        if (Parms->ModeParm == -1) {
            FillCArray(SSCoefDia, 0, IntPow2(*Parms->SBasisBound));
            for (PriSCount = 0; PriSCount < Parms->SBasisBound[2]; PriSCount++)
                if (SBasis[PriSCount].P[0].V[Singlet] == 0) {
                    SSCoefDia[PriSCount * (*Parms->SBasisBound + 1)] = (ComplexType) 1. / Parms->SiteTotal;
                    for (SecSCount = 0; SecSCount < PriSCount; SecSCount++)
                        if (SBasis[SecSCount].P[0].V[Singlet] == 0) {
                            SSCoefDia[PriSCount * *Parms->SBasisBound + SecSCount] = (ComplexType) 1. / Parms->SiteTotal;
                            SSCoefDia[SecSCount * *Parms->SBasisBound + PriSCount] = (ComplexType) 1. / Parms->SiteTotal;
                        }
                }
            CMatrixTimesRMatrix(SSCoefWork, SSCoefDia, SEigenVectorsT, *Parms->SBasisBound);
            RMatrixTimesCMatrix(SSCoef, SMatrix.Array, SSCoefWork, *Parms->SBasisBound);

            MaxPop = (RealType) creal(SSCoef[0]); MaxPopId = 0;
            for (PriSCount = 1; PriSCount < CutOff; PriSCount++)
                if (creal(SSCoef[PriSCount * (*Parms->SBasisBound + 1)]) > MaxPop) {
                    MaxPopId = PriSCount;
                    MaxPop = (RealType) SSCoef[PriSCount * (*Parms->SBasisBound + 1)];
                }

            FillCArray(SSCoef, 0, IntPow2(*Parms->SBasisBound));
            SSCoef[MaxPopId * (*Parms->SBasisBound + 1)] = 1;
            CMatrixTimesRMatrix(SSCoefWork, SSCoef, SMatrix.Array, *Parms->SBasisBound);
            RMatrixTimesCMatrix(SSCoefDia, SEigenVectorsT, SSCoefWork, *Parms->SBasisBound);
        }
        else {
            if (Parms->PrepareAdiab) {
                FillCArray(SSCoef, 0, IntPow2(*Parms->SBasisBound));
                SSCoef[Parms->ModeParm * (*Parms->SBasisBound + 1)] = 1;
                CMatrixTimesRMatrix(SSCoefWork, SSCoef, SMatrix.Array, *Parms->SBasisBound);
                RMatrixTimesCMatrix(SSCoefDia, SEigenVectorsT, SSCoefWork, *Parms->SBasisBound);
            }
            else {
                FillCArray(SSCoefDia, 0, IntPow2(*Parms->SBasisBound));
                SSCoefDia[Parms->ModeParm * (*Parms->SBasisBound + 1)] = 1;
                CMatrixTimesRMatrix(SSCoefWork, SSCoefDia, SEigenVectorsT, *Parms->SBasisBound);
                RMatrixTimesCMatrix(SSCoef, SMatrix.Array, SSCoefWork, *Parms->SBasisBound);
            }
        }

        PopulationsElement = Populations; PopulationsAdiabElement = PopulationsAdiab;
        for (TimeCount = 0; TimeCount < *Parms->TimeTotal; TimeCount++) {
            UpdateProgress(&Progress);

            for (PriSCount = 0; PriSCount < *Parms->SBasisBound; PriSCount++) {
                (*PopulationsElement++) += (RealType) creal(SSCoefDia[(*Parms->SBasisBound + 1) * PriSCount]);
                (*PopulationsAdiabElement++) += (RealType) creal(SSCoef[(*Parms->SBasisBound + 1) * PriSCount]);
            }
            Propagate_Redfield(SSCoef, SSCoefWork, SEigenValues, SEigenValues, SSRedfieldTensor, *Parms->SBasisBound, *Parms->SBasisBound, CutOff, CutOff, Parms);

            CMatrixTimesRMatrixCutOff(SSCoefWork, SSCoef, SMatrix.Array, *Parms->SBasisBound, CutOff);
            RMatrixTimesCMatrixCutOff(SSCoefDia, SEigenVectorsT, SSCoefWork, *Parms->SBasisBound, CutOff);
        }

        FreeMaskCArray(&SSRedfieldTensor);
        free(SBathFT); free(SReperator);
    }

    RealTimesRArray(Populations, (RealType) 1 / Parms->SampleTotal, *Parms->SBasisBound * *Parms->TimeTotal);
    RealTimesRArray(PopulationsAdiab, (RealType) 1 / Parms->SampleTotal, *Parms->SBasisBound * *Parms->TimeTotal);
    ExtractSPopulations(SPopulations, Populations, SBasis, Parms);
	ExtractCPopulations(CPopulations, Populations, SBasis, Parms);

	CreateGrid(TimeGrid, 0, (*Parms->TimeTotal - 1) * Parms->TimeStep, *Parms->TimeTotal);
	SaveRArray("Populations", Populations, *Parms->TimeTotal, *Parms->SBasisBound, false, 10, TimeGrid, false);
	SaveRArray("APopulations", PopulationsAdiab, *Parms->TimeTotal, *Parms->SBasisBound, false, 10, TimeGrid, false);
    SaveRArray("SPopulations", SPopulations, *Parms->TimeTotal, Parms->SiteTotal, false, false, TimeGrid, false);
	SaveRArray("CPopulations", CPopulations, *Parms->TimeTotal, 3, false, false, TimeGrid, false);

	CloseDiag(Parms->DAndC);

    FreeMaskArray(&OffDiagDisorder.LL); FreeMaskArray(&OffDiagDisorder.HH); FreeMaskArray(&OffDiagDisorder.HL);
    FreeMaskArray(&StaticHamiltonian); FreeMaskArray(&SMatrix);

	free(SBasis);
	free(SEigenVectorsT); free(SEigenValues); free(Populations); free(PopulationsAdiab); free(SPopulations);
	free(TimeGrid); free(CPopulations);
	free(SSCoefDia); free(SSCoef); free(SSCoefWork);
}