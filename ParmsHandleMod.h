#ifndef ParmsHandleMod
#define ParmsHandleMod

void InitializeParms(ParmStruct *Parms);
void ReadParms(char *Args[], ParmStruct *Parms);
void DerivedParms(ParmStruct *Parms);

#endif