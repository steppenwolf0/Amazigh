#ifndef CELL_H
#define CELL_H
void
initConstsFK(double* CONSTANTS, double* RATES, double *STATES);
void
computeRatesFK(double VOI, double* CONSTANTS, double* RATES, double* STATES, double* ALGEBRAIC);
void
computeVariablesFK(double VOI, double* CONSTANTS, double* RATES, double* STATES, double* ALGEBRAIC);
void
run_FemtomKarma();
void
run_FemtomKarmaRK2();
void
run_LuoRudyRK2();
void
run_Hodgkin();
void 
run_LuoRudyRK4();
void
run_FemtomKarmaMod();

//cardio
void initConsts_FemtomKarmaMod(double* CONSTANTS, double* RATES, double *STATES);
void computeRates_FemtomKarmaMod(double VOI, double* CONSTANTS, double* RATES, double* STATES, double* ALGEBRAIC, double Istim);
void computeVariables_FemtomKarmaMod(double VOI, double* CONSTANTS, double* RATES, double* STATES, double* ALGEBRAIC);
#endif