#ifndef FEMTESTS_H
#define FEMTESTS_H
void FEMTest1();
void FEMTest2();
double stiffFunc(double x, double y, double z);
double ux(double x, double y, double z, double t);
double uy(double x, double y, double z, double t);
double uz(double x, double y, double z, double t);
double fFunc(double x, double y, double z, double t);
double uFunc(double x, double y, double z, double t);

void FEMTestTime1();
void FEMTestTime2();
void FEMTestTime3();
void FEMTestTime4();
#endif