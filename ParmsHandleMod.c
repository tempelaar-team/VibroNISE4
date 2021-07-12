#include "GlobalsMod.h"
#include "ParmsMod.h"
#include "ToolsMod.h"
#include "ParmsHandleMod.h"

#include <string.h>


static void ReadInputFile(ParmStruct *Parms, char *EnergyChar, char *EnergyDerivedChar, char *TimeChar, const char *FileName);
static int ReadInt(int *Variable, char *ReadLine, char *DisplayNumber);
static int ReadDouble(RealType *Variable, char *ReadLine, char *DisplayNumber);

#include <math.h>


/////////////
// Globals //
/////////////


void InitializeParms(ParmStruct *Parms)
{
	PrintMessage("Applying default settings.", "", "");

	Parms->DAndC = true;
    Parms->SampleTotal = 1;
	Parms->TimeStep = 1;
	Parms->ACellTotal = 1;
	Parms->BCellTotal = 1;
	Parms->SubTotal = 1;
	Parms->Periodic = true;
	Parms->AConstant = 1; Parms->BConstant = 1;
	Parms->ABAngle = 90;
	Parms->VTotal = 1;
	Parms->SpectrumTotal = 500;
	Parms->OptDielectric = 1;
	Parms->RescaleEA = 1;
	Parms->Polarizations[0] = -1;

	sprintf(Parms->DipolesFile, "|");
	sprintf(Parms->DDCouplingsFile, "|");
	sprintf(Parms->TripletDDCouplingsFile, "|");
	sprintf(Parms->LLCouplingsFile, "|");
	sprintf(Parms->HHCouplingsFile, "|");
	sprintf(Parms->HLCouplingsFile, "|");
	sprintf(Parms->LLHLCouplingsFile, "|");
	sprintf(Parms->LLHHCouplingsFile, "|");
	sprintf(Parms->LnLnHHCouplingsFile, "|");
	sprintf(Parms->EHEnergiesFile, "|");
	sprintf(Parms->BathFile, "|");
}


void ReadParms(char *Args[], ParmStruct *Parms)
{
	char *EnergyChar = malloc(6), *EnergyDerivedChar = malloc(8), *TimeChar = malloc(5);
	char InputFile[100];

	strcpy(EnergyChar, " /cm."); strcpy(EnergyDerivedChar, " A/cm."); strcpy(TimeChar, " fs.");

	sprintf(InputFile, "%s", Args[1]);
	ReadInputFile(Parms, EnergyChar, EnergyDerivedChar, TimeChar, InputFile);

	free(EnergyChar); free(EnergyDerivedChar); free(TimeChar);
}


void DerivedParms(ParmStruct *Parms)
{
	char DisplayNumber[100];
	int Count;

	PrintMessage("Derived parameters are determined.", "", "");

	if (Parms->EnergyUnit == WaveNumbers)				Parms->EnergyToFreq = InvCmToFreq;
	else if (Parms->EnergyUnit == MilliElectronVolt)	Parms->EnergyToFreq = MeVToFreq;
	else												Parms->EnergyToFreq = 1;

	for (Count = 0; Count < 3; Count++) {
		if (!Parms->HBroad[Count] && !Parms->LifeTime[Count]) {
			if (Parms->EnergyUnit == WaveNumbers)				Parms->HBroad[Count] = 40;
			else if (Parms->EnergyUnit == MilliElectronVolt)	Parms->HBroad[Count] = 5;
			else												Parms->HBroad[Count] = .05;
		}

		if (!Parms->LifeTime[Count]) {
			Parms->LifeTime[Count] = 1 / (2 * Parms->HBroad[Count] * Parms->EnergyToFreq);
			if (Count == 0) {
				RealToString(DisplayNumber, Parms->LifeTime[0]);
				PrintMessage(" - A-polarized excited state lifetime equals ", DisplayNumber, " fs.");
			}
		}
		else if (!Parms->HBroad[Count]) Parms->HBroad[Count] = Parms->EnergyToFreq / (2 * Parms->LifeTime[Count]);
	}

	Parms->TTEnergyMean -= Parms->SEnergyMean;
	Parms->TT_nEnergyMean -= 2 * Parms->SEnergyMean;

	Parms->CellTotal = Parms->ACellTotal * Parms->BCellTotal;
	Parms->SiteTotal = Parms->CellTotal *  Parms->SubTotal;
	Parms->ABAngle *= Pi / 180;

	if (*Parms->Polarizations != -1)
		for (Count = 0; Count < 4; Count++) {
			Parms->Polarizations[4 + Count] = (RealType) sin((double) Parms->Polarizations[Count] * Pi / 180);
			Parms->Polarizations[Count]		= (RealType) cos((double) Parms->Polarizations[Count] * Pi / 180);
		}
}


////////////
// Locals //
////////////


static void ReadInputFile(ParmStruct *Parms, char *EnergyChar, char *EnergyDerivedChar, char *TimeChar, const char *FileName)
{
	char ReadLine[100], DisplayNumber[100];
	FILE *FileID;


	PrintMessage("Reading input file '", FileName, "'.");

	FileID = fopen(FileName, "r");

	while (!feof(FileID)) {
		fscanf(FileID, "%s\n", ReadLine);

		if (!strcmp(ReadLine, "InputFile")) {
			fscanf(FileID, "%s\n", ReadLine);
			PrintMessage(" - Additional input will be read from file '", ReadLine, "'.");
			ReadInputFile(Parms, EnergyChar, EnergyDerivedChar, TimeChar, ReadLine);
		}
		if (!strcmp(ReadLine, "Mode")) {
			fscanf(FileID, "%s\n", ReadLine);
			if (!strcmp(ReadLine, "Debug"))
				Parms->Mode = ModeDebug,			PrintMessage(" - Running in debugging mode.", "", "");
			else if (!strcmp(ReadLine, "Direct1D"))
				Parms->Mode = ModeDirect1D,			PrintMessage(" - Calculating linear spectra through Fermi's Golden Rule.", "", "");
			else if (!strcmp(ReadLine, "FTMapAnalysis"))
				Parms->Mode = ModeFTMapAnalysis,		PrintMessage(" - Performing FT amplitude map analysis.", "", "");
			else if (!strcmp(ReadLine, "Response1D"))
				Parms->Mode = ModeResponse1D,		PrintMessage(" - Calculating linear response.", "", "");
			else if (!strcmp(ReadLine, "Response2D"))
				Parms->Mode = ModeResponse2D,		PrintMessage(" - Calculating 2D response.", "", "");
			else if (!strcmp(ReadLine, "Dynamics")) {
				Parms->Mode = ModeDynamics,			PrintMessage(" - Calculating dynamics.", "", "");
				fgets(ReadLine, sizeof(ReadLine), FileID);
				ReadInt(&Parms->ModeParm, ReadLine, DisplayNumber);
				PrintMessage(" - Exciting density matrix element ", DisplayNumber, ".");
			}
		}
		if (!strcmp(ReadLine, "NoDAndC")) {
			PrintMessage(" - Divide & Conquer is not used.", "", "");
			Parms->DAndC = false;
		}
		else if (!strcmp(ReadLine, "WaveNumbers")) {
			PrintMessage(" - Energies taken in /cm.", "", "");
			strcpy(EnergyChar, " /cm."); strcpy(EnergyDerivedChar, " A/cm."); strcpy(TimeChar, " fs.");
			Parms->EnergyUnit = WaveNumbers;
		}
		else if (!strcmp(ReadLine, "meV")) {
			PrintMessage(" - Energies taken in millielectronvolt.", "", "");
			strcpy(EnergyChar, " meV."); strcpy(EnergyDerivedChar, " A meV."); strcpy(TimeChar, " fs.");
			Parms->EnergyUnit = MilliElectronVolt;
		}
		else if (!strcmp(ReadLine, "Frequency")) {
			PrintMessage(" - Energies taken in /fs.", "", "");
			strcpy(EnergyChar, " /fs."); strcpy(EnergyDerivedChar, " A/fs."); strcpy(TimeChar, " fs.");
			Parms->EnergyUnit = Frequency;
		}
		else if (!strcmp(ReadLine, "UnitLess")) {
			PrintMessage(" - Working unitless.", "", "");
			strcpy(EnergyChar, "."); strcpy(EnergyDerivedChar, "."); strcpy(TimeChar, ".");
			Parms->EnergyUnit = UnitLess;
		}
		else if (!strcmp(ReadLine, "NoSE")) {
			Parms->NoSE = true;						PrintMessage(" - Stimulated emission not included in 2D spectra.", "", "");
		}
		else if (!strcmp(ReadLine, "K0Approach")) {
			Parms->K0Approach = true;				PrintMessage(" - K = 0 approach enabled.", "", "");
		}
		else if (!strcmp(ReadLine, "PolResolved")) {
			Parms->PolResolved = true;				PrintMessage(" - Resolving absorption polarizations.", "", "");
		}
		else if (!strcmp(ReadLine, "DoubleCrossPol")) {
			Parms->DoubleCrossPol = true;			PrintMessage(" - Applying double-cross polarization sequence.", "", "");
		}
		else if (!strcmp(ReadLine, "DDCouplingsFormatNick")) {
			Parms->DDCouplingsFormat = DDCouplingsFormatNick;
			PrintMessage(" - Adopting Nick's format for reading dipole-dipole couplings.", "", "");
		}
		else if (!strcmp(ReadLine, "DDCouplingsFormatOld")) {
			Parms->DDCouplingsFormat = DDCouplingsFormatOld;
			PrintMessage(" - Adopting old format for reading dipole-dipole couplings.", "", "");
		}
		else if (!strcmp(ReadLine, "AbsModeNick")) {
			Parms->AbsMode = AbsModeNick;			PrintMessage(" - Adopting Nick's formula for absorption.", "", "");
		}
		else if (!strcmp(ReadLine, "AbsModeTim")) {
			Parms->AbsMode = AbsModeTim;			PrintMessage(" - Adopting Tim's formula for absorption.", "", "");
		}
		else if (!strcmp(ReadLine, "Normalize")) {
			Parms->Normalize = true;				PrintMessage(" - Linear spectra are normalized.", "", "");
		}
		else if (!strcmp(ReadLine, "Adiabatic")) {
			Parms->Adiab = true;					PrintMessage(" - Calculations are performed in the adiabatic basis (if applicable).", "", "");
		}
		else if (!strcmp(ReadLine, "PrepareAdiab")) {
			Parms->PrepareAdiab = true;				PrintMessage(" - Adiabatic state is prepared.", "", "");
		}
		else if (!strcmp(ReadLine, "NoPeriodic")) {
			Parms->Periodic = false;				PrintMessage(" - Periodic boundaries disabled.", "", "");
		}
		else if (!strcmp(ReadLine, "TwoSubs")) {
			Parms->SubTotal = 2;					PrintMessage(" - Enabling two sublattices.", "", "");
		}
		else if (!strcmp(ReadLine, "NoCTCTCoupling")) {
			Parms->NoCTCTCoupling = true;			PrintMessage(" - CT-CT coupling disabled.", "", "");
		}
		else if (!strcmp(ReadLine, "TTnOnly")) {
			Parms->TT_nOnly = true;			PrintMessage(" - Double singlet states excluded.", "", "");
		}
		else if (!strcmp(ReadLine, "Gaussian")) {
			Parms->Gaussian = true;					PrintMessage(" - Applying Gaussian lineshapes.", "", "");
		}
		else if (!strcmp(ReadLine, "VRadiusInf")) {
			Parms->VRadius = -1;					PrintMessage(" - No maximum vibrational separation applied.", "", "");
		}
		else if (!strcmp(ReadLine, "EHRadiusInf")) {
			Parms->EHRadius = -1;					PrintMessage(" - No maximum electron-hole separation applied.", "", "");
		}
		else if (!strcmp(ReadLine, "TTRadiusInf")) {
			Parms->TTRadius = -1;					PrintMessage(" - No maximum triplet-triplet separation applied.", "", "");
		}
		else if (!strcmp(ReadLine, "Dipoles")) {
			fscanf(FileID, "%s\n", Parms->DipolesFile);
			PrintMessage(" - Associating the permanent dipoles with file '", Parms->DipolesFile, "'.");
		}
        else if (!strcmp(ReadLine, "SampleTotal")) {
            fgets(ReadLine, sizeof(ReadLine), FileID);
            ReadInt(&Parms->SampleTotal, ReadLine, DisplayNumber);
            PrintMessage(" - Setting the number of samples to ", DisplayNumber, ".");
        }
		else if (!strcmp(ReadLine, "TimeStep")) {
			fgets(ReadLine, sizeof(ReadLine), FileID);
			ReadDouble(&Parms->TimeStep, ReadLine, DisplayNumber);
			PrintMessage(" - Setting the time step size to ", DisplayNumber, TimeChar);
		}
		else if (!strcmp(ReadLine, "TimeTotal")) {
			fgets(ReadLine, sizeof(ReadLine), FileID);
			ReadInt(&Parms->TimeTotal[0], ReadLine, DisplayNumber);
			PrintMessage(" - Setting first/third time interval to ", DisplayNumber, " time steps.");
			if (ReadInt(&Parms->TimeTotal[1], ReadLine, DisplayNumber))
				PrintMessage(" - Setting second time interval to ", DisplayNumber, " time steps.");
		}
		else if (!strcmp(ReadLine, "CellTotal")) {
			fgets(ReadLine, sizeof(ReadLine), FileID);
			ReadInt(&Parms->ACellTotal, ReadLine, DisplayNumber);
			PrintMessage(" - Setting the number of cells in the a-direction to ", DisplayNumber, ".");
			if (ReadInt(&Parms->BCellTotal, ReadLine, DisplayNumber))
				PrintMessage(" - Setting the number of cells in the b-direction to ", DisplayNumber, ".");
		}
		else if (!strcmp(ReadLine, "LatticeConstant")) {
			fgets(ReadLine, sizeof(ReadLine), FileID);
			ReadDouble(&Parms->AConstant, ReadLine, DisplayNumber);
			PrintMessage(" - Setting the lattice constant in the a-direction to ", DisplayNumber, " A.");
			if (ReadDouble(&Parms->BConstant, ReadLine, DisplayNumber))
				PrintMessage(" - Setting the lattice constant in the b-direction to ", DisplayNumber, " A.");
		}
		else if (!strcmp(ReadLine, "VSlope")) {
			fgets(ReadLine, sizeof(ReadLine), FileID);
			ReadDouble(&Parms->VSlope, ReadLine, DisplayNumber);
			PrintMessage(" - Setting the vibrational slope to ", DisplayNumber, " A.");
		}
		else if (!strcmp(ReadLine, "VRadius")) {
			fgets(ReadLine, sizeof(ReadLine), FileID);
			ReadDouble(&Parms->VRadius, ReadLine, DisplayNumber);
			PrintMessage(" - Setting the maximum vibrational separation to ", DisplayNumber, " A.");
		}
		else if (!strcmp(ReadLine, "EHRadius")) {
			fgets(ReadLine, sizeof(ReadLine), FileID);
			ReadDouble(&Parms->EHRadius, ReadLine, DisplayNumber);
			PrintMessage(" - Setting the maximum electron-hole separation to ", DisplayNumber, " A.");
		}
		else if (!strcmp(ReadLine, "TTRadius")) {
			fgets(ReadLine, sizeof(ReadLine), FileID);
			ReadDouble(&Parms->TTRadius, ReadLine, DisplayNumber);
			PrintMessage(" - Setting the maximum triplet-triplet separation to ", DisplayNumber, " A.");
		}
		else if (!strcmp(ReadLine, "ABAngle")) {
			fgets(ReadLine, sizeof(ReadLine), FileID);
			ReadDouble(&Parms->ABAngle, ReadLine, DisplayNumber);
			PrintMessage(" - Setting the angle between a and b to ", DisplayNumber, " degrees.");
		}
		else if (!strcmp(ReadLine, "SEnergyMean")) {
			fgets(ReadLine, sizeof(ReadLine), FileID);
			ReadDouble(&Parms->SEnergyMean, ReadLine, DisplayNumber);
			PrintMessage(" - Setting the average singlet transition energy to ", DisplayNumber, EnergyChar);
		}
		else if (!strcmp(ReadLine, "CouplingNN")) {
			fgets(ReadLine, sizeof(ReadLine), FileID);
			ReadDouble(&Parms->DDCouplingNN, ReadLine, DisplayNumber);
			PrintMessage(" - Setting the nearest-neighbor dipole-dipole coupling to ", DisplayNumber, EnergyChar);
			ReadDouble(&Parms->TripletDDCouplingNN, ReadLine, DisplayNumber);
			PrintMessage(" - Setting the nearest-neighbor triplet dipole-dipole coupling to ", DisplayNumber, EnergyChar);
			if (ReadDouble(&Parms->LLCouplingNN, ReadLine, DisplayNumber))
			if (Parms->LLCouplingNN)
				PrintMessage(" - Setting the nearest-neighbor LUMO-LUMO coupling to ", DisplayNumber, EnergyChar);
			if (ReadDouble(&Parms->HHCouplingNN, ReadLine, DisplayNumber))
			if (Parms->HHCouplingNN)
				PrintMessage(" - Setting the nearest-neighbor HOMO-HOMO coupling to ", DisplayNumber, EnergyChar);
			if (ReadDouble(&Parms->HLCouplingNN, ReadLine, DisplayNumber))
			if (Parms->HLCouplingNN)
				PrintMessage(" - Setting the nearest-neighbor HOMO-LUMO coupling to ", DisplayNumber, EnergyChar);
			if (ReadDouble(&Parms->LLHLCouplingNN, ReadLine, DisplayNumber))
			if (Parms->LLHLCouplingNN)
				PrintMessage(" - Setting the nearest-neighbor double electron transfer (LUMO-LUMO & HOMO-LUMO) coupling to ", DisplayNumber, EnergyChar);
			if (ReadDouble(&Parms->LLHHCouplingNN, ReadLine, DisplayNumber))
				PrintMessage(" - Setting the nearest-neighbor double electron transfer (LUMO-LUMO & HOMO-HOMO) coupling to ", DisplayNumber, EnergyChar);
		}
		else if (!strcmp(ReadLine, "DDCouplings")) {
			fscanf(FileID, "%s\n", Parms->DDCouplingsFile);
			PrintMessage(" - Associating dipole-dipole couplings with file '", Parms->DDCouplingsFile, "'.");
		}
		else if (!strcmp(ReadLine, "TripletDDCouplings")) {
			fscanf(FileID, "%s\n", Parms->TripletDDCouplingsFile);
			PrintMessage(" - Associating triplet dipole-dipole couplings with file '", Parms->TripletDDCouplingsFile, "'.");
		}
		else if (!strcmp(ReadLine, "LLCouplings")) {
			fscanf(FileID, "%s\n", Parms->LLCouplingsFile);
			PrintMessage(" - Associating LUMO-LUMO couplings with file '", Parms->LLCouplingsFile, "'.");
		}
		else if (!strcmp(ReadLine, "HHCouplings")) {
			fscanf(FileID, "%s\n", Parms->HHCouplingsFile);
			PrintMessage(" - Associating HOMO-HOMO couplings with file '", Parms->HHCouplingsFile, "'.");
		}
		else if (!strcmp(ReadLine, "HLCouplings")) {
			fscanf(FileID, "%s\n", Parms->HLCouplingsFile);
			PrintMessage(" - Associating HOMO-LUMO couplings with file '", Parms->HLCouplingsFile, "'.");
		}
		else if (!strcmp(ReadLine, "LLHLCouplings")) {
			fscanf(FileID, "%s\n", Parms->LLHLCouplingsFile);
			PrintMessage(" - Associating LUMO-LUMO & HOMO-LUMO couplings with file '", Parms->LLHLCouplingsFile, "'.");
		}
		else if (!strcmp(ReadLine, "LLHHCouplings")) {
			fscanf(FileID, "%s\n", Parms->LLHHCouplingsFile);
			PrintMessage(" - Associating LUMO-LUMO & HOMO-HUMO couplings with file '", Parms->LLHHCouplingsFile, "'.");
		}
		else if (!strcmp(ReadLine, "LnLnHHCouplings")) {
			fscanf(FileID, "%s\n", Parms->LnLnHHCouplingsFile);
			PrintMessage(" - Associating LUMO_n-LUMO_n & HOMO-HUMO couplings with file '", Parms->LnLnHHCouplingsFile, "'.");
		}
		else if (!strcmp(ReadLine, "EHEnergies")) {
			fscanf(FileID, "%s\n", Parms->EHEnergiesFile);
			PrintMessage(" - Associating Coulomb energies with file '", Parms->EHEnergiesFile, "'.");
		}
		else if (!strcmp(ReadLine, "EHSeparated")) {
			fgets(ReadLine, sizeof(ReadLine), FileID);
			ReadDouble(&Parms->EHSeparated, ReadLine, DisplayNumber);
			PrintMessage(" - Setting the energy of a separated electron-hole pair to ", DisplayNumber, EnergyChar);
		}
		else if (!strcmp(ReadLine, "EHScaling")) {
			fgets(ReadLine, sizeof(ReadLine), FileID);
			ReadDouble(&Parms->EHScaling, ReadLine, DisplayNumber);
			PrintMessage(" - Setting the electron-hole 1/separation energy constant to ", DisplayNumber, EnergyDerivedChar);
		}
		else if (!strcmp(ReadLine, "TTEnergyMean")) {
			fgets(ReadLine, sizeof(ReadLine), FileID);
			ReadDouble(&Parms->TTEnergyMean, ReadLine, DisplayNumber);
			PrintMessage(" - Setting the average triplet-triplet energy to ", DisplayNumber, EnergyChar);
		}
		else if (!strcmp(ReadLine, "TTnEnergyMean")) {
			fgets(ReadLine, sizeof(ReadLine), FileID);
			ReadDouble(&Parms->TT_nEnergyMean, ReadLine, DisplayNumber);
			PrintMessage(" - Setting the average triplet-triplet_n energy to ", DisplayNumber, EnergyChar);
		}
		else if (!strcmp(ReadLine, "TTnDipole")) {
			fgets(ReadLine, sizeof(ReadLine), FileID);
			ReadDouble(&Parms->TT_nDipole, ReadLine, DisplayNumber);
			PrintMessage(" - Setting the (relative) triplet-triplet_n dipole moment to ", DisplayNumber, ".");
		}
		else if (!strcmp(ReadLine, "RescaleEA")) {
			fgets(ReadLine, sizeof(ReadLine), FileID);
			ReadDouble(&Parms->RescaleEA, ReadLine, DisplayNumber);
			PrintMessage(" - Setting the relative excited state absorption rescaling to ", DisplayNumber, ".");
		}
		else if (!strcmp(ReadLine, "VMax")) {
			fgets(ReadLine, sizeof(ReadLine), FileID);
			ReadInt(&Parms->VTotal, ReadLine, DisplayNumber);
			PrintMessage(" - Setting the maximum number of vibrations to ", DisplayNumber, ".");
			Parms->VTotal++;
		}
		else if (!strcmp(ReadLine, "VEnergy")) {
			fgets(ReadLine, sizeof(ReadLine), FileID);
			ReadDouble(&Parms->VEnergy, ReadLine, DisplayNumber);
			PrintMessage(" - Setting the vibrational energy to ", DisplayNumber, EnergyChar);
		}
		else if (!strcmp(ReadLine, "HuangRhys")) {
			fgets(ReadLine, sizeof(ReadLine), FileID);
			ReadDouble(&Parms->HuangRhys[1], ReadLine, DisplayNumber);
			PrintMessage(" - Setting singlet-ground Huang-Rhys factor to ", DisplayNumber, ".");
			if (ReadDouble(&Parms->HuangRhys[2], ReadLine, DisplayNumber))
				PrintMessage(" - Setting electron-ground Huang-Rhys factor to ", DisplayNumber, ".");
			if (ReadDouble(&Parms->HuangRhys[3], ReadLine, DisplayNumber))
				PrintMessage(" - Setting hole-ground Huang Rhys-factor to ", DisplayNumber, ".");
			if (ReadDouble(&Parms->HuangRhys[4], ReadLine, DisplayNumber))
				PrintMessage(" - Setting triplet-ground Huang-Rhys factor to ", DisplayNumber, ".");
			if (ReadDouble(&Parms->HuangRhys[5], ReadLine, DisplayNumber))
				PrintMessage(" - Setting triplet_n-ground Huang-Rhys factor to ", DisplayNumber, ".");
		}
		else if (!strcmp(ReadLine, "BathCorrelation")) {
			fscanf(FileID, "%s\n", Parms->BathFile);
			PrintMessage(" - Associating bath spectral density with file '", Parms->BathFile, "'.");
		}
		else if (!strcmp(ReadLine, "BathReorganization")) {
			fgets(ReadLine, sizeof(ReadLine), FileID);
			ReadDouble(&Parms->BathReorganization, ReadLine, DisplayNumber);
			PrintMessage(" - Setting the bath reorganization energy to ", DisplayNumber, EnergyChar);
		}
        else if (!strcmp(ReadLine, "LLSigma")) {
            fgets(ReadLine, sizeof(ReadLine), FileID);
            ReadDouble(&Parms->LLSigma, ReadLine, DisplayNumber);
            PrintMessage(" - Setting the relative LUMO-LUMO couplings disorder width to ", DisplayNumber, ".");
        }
        else if (!strcmp(ReadLine, "HHSigma")) {
            fgets(ReadLine, sizeof(ReadLine), FileID);
            ReadDouble(&Parms->HHSigma, ReadLine, DisplayNumber);
            PrintMessage(" - Setting the relative HOMO-HOMO couplings disorder width to ", DisplayNumber, ".");
        }
        else if (!strcmp(ReadLine, "HLSigma")) {
            fgets(ReadLine, sizeof(ReadLine), FileID);
            ReadDouble(&Parms->HLSigma, ReadLine, DisplayNumber);
            PrintMessage(" - Setting the relative HOMO-LUMO couplings disorder width to ", DisplayNumber, ".");
        }
		else if (!strcmp(ReadLine, "HBroad")) {
			fgets(ReadLine, sizeof(ReadLine), FileID);
			ReadDouble(&Parms->HBroad[0], ReadLine, DisplayNumber);
			PrintMessage(" - Setting a-polarized (and b/c unless mentioned) homogeneous broadening to ", DisplayNumber, EnergyChar);
			if (ReadDouble(&Parms->HBroad[1], ReadLine, DisplayNumber))
				PrintMessage(" - Setting b-polarized (and c unless mentioned) homogeneous broadening to ", DisplayNumber, EnergyChar);
			else Parms->HBroad[1] = Parms->HBroad[0];
			if (ReadDouble(&Parms->HBroad[2], ReadLine, DisplayNumber))
				PrintMessage(" - Setting c-polarized homogeneous broadening to ", DisplayNumber, EnergyChar);
			else Parms->HBroad[2] = Parms->HBroad[1];
		}
		else if (!strcmp(ReadLine, "LifeTime")) { // Do we want to keep this, or solely have HBroad?
			fgets(ReadLine, sizeof(ReadLine), FileID);
			ReadDouble(&Parms->LifeTime[0], ReadLine, DisplayNumber);
			PrintMessage(" - Setting a-polarized excited state lifetime to ", DisplayNumber, TimeChar);
			if (ReadDouble(&Parms->LifeTime[1], ReadLine, DisplayNumber))
				PrintMessage(" - Setting b-polarized excited state lifetime to ", DisplayNumber, TimeChar);
			else {
				Parms->LifeTime[1] = Parms->LifeTime[0];
				PrintMessage(" - Setting b-polarized excited state lifetime accordingly.", "", "");
			}
			if (ReadDouble(&Parms->LifeTime[2], ReadLine, DisplayNumber))
				PrintMessage(" - Setting c-polarized excited state lifetime to ", DisplayNumber, TimeChar);
			else {
				Parms->LifeTime[2] = Parms->LifeTime[1];
				PrintMessage(" - Setting c-polarized excited state lifetime accordingly.", "", "");
			}
		}
        else if (!strcmp(ReadLine, "Temperature")) {
            fgets(ReadLine, sizeof(ReadLine), FileID);
            ReadDouble(&Parms->Temperature, ReadLine, DisplayNumber);
            PrintMessage(" - Setting thermal quantum to ", DisplayNumber, EnergyChar);
        }
		else if (!strcmp(ReadLine, "OptDielectric")) {
			fgets(ReadLine, sizeof(ReadLine), FileID);
			ReadDouble(&Parms->OptDielectric, ReadLine, DisplayNumber);
			PrintMessage(" - Setting optical dielectric constant to ", DisplayNumber, ".");
		}
		else if (!strcmp(ReadLine, "EnergyBounds")) {
			fgets(ReadLine, sizeof(ReadLine), FileID);
			ReadDouble(&Parms->EnergyBounds[0], ReadLine, DisplayNumber);
			PrintMessage(" - Setting the energy lower bound to ", DisplayNumber, EnergyChar);
			ReadDouble(&Parms->EnergyBounds[1], ReadLine, DisplayNumber);
			PrintMessage(" - Setting the energy upper bound to ", DisplayNumber, EnergyChar);
		}
		else if (!strcmp(ReadLine, "LaserBounds")) {
			fgets(ReadLine, sizeof(ReadLine), FileID);
			ReadDouble(&Parms->LaserBounds[0], ReadLine, DisplayNumber);
			PrintMessage(" - Setting the laser lower bound to ", DisplayNumber, EnergyChar);
			ReadDouble(&Parms->LaserBounds[1], ReadLine, DisplayNumber);
			PrintMessage(" - Setting the laser upper bound to ", DisplayNumber, EnergyChar);
		}
		else if (!strcmp(ReadLine, "SpectrumTotal")) {
			fgets(ReadLine, sizeof(ReadLine), FileID);
			ReadInt(&Parms->SpectrumTotal, ReadLine, DisplayNumber);
			PrintMessage(" - Setting the spectral resolution to ", DisplayNumber, ".");
			Parms->SpectrumTotal++;
		}
		else if (!strcmp(ReadLine, "Polarizations")) {
			fgets(ReadLine, sizeof(ReadLine), FileID);
			ReadDouble(&Parms->Polarizations[0], ReadLine, DisplayNumber);
			PrintMessage(" - First pulse polarized as ", DisplayNumber, " degrees.");
			ReadDouble(&Parms->Polarizations[1], ReadLine, DisplayNumber);
			PrintMessage(" - Second pulse polarized as ", DisplayNumber, " degrees.");
			ReadDouble(&Parms->Polarizations[2], ReadLine, DisplayNumber);
			PrintMessage(" - Third pulse polarized as ", DisplayNumber, " degrees.");
			ReadDouble(&Parms->Polarizations[3], ReadLine, DisplayNumber);
			PrintMessage(" - Fourth pulse polarized as ", DisplayNumber, " degrees.");
		}
	}
	fclose(FileID);
	PrintMessage("File '", FileName, "' closed.");
}


static int ReadInt(int *Variable, char *ReadLine, char *DisplayNumber)
{
	long int ReadValue;
	char *ReadLineEntry = ReadLine;
	
	if (*ReadLineEntry == '\n') return false;
	ReadValue = strtol(ReadLineEntry, &ReadLineEntry, 10);
	*Variable = (int) ReadValue;
	IntToString(DisplayNumber, *Variable);
	sprintf(ReadLine, "%s\n", ReadLineEntry);
	
	return true;
}


static int ReadDouble(RealType *Variable, char *ReadLine, char *DisplayNumber)
{
	double ReadValue;
	char *ReadLineEntry = ReadLine;
	
	if (*ReadLineEntry == '\n') return false;
	ReadValue = strtod(ReadLineEntry, &ReadLineEntry);
	*Variable = (RealType) ReadValue;
	RealToString(DisplayNumber, *Variable);
	sprintf(ReadLine, "%s\n", ReadLineEntry);
	
	return true;
}
