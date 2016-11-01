#ifndef FUNCTION_H
#define FUNCTION_H
#define _USE_MATH_DEFINES

#include <cmath>
#include <vector>
#include <stdio.h>
#include <iostream>
#include "iodata.h"
#include "simdata.h"
#include "coefficient_thomas.h"
#include <string>
#include <fstream>
#include <omp.h>
int getRows(char* filename);
void writeDataPlot(char* filename, iodata::odata odata);
void dataRW(char* filename, const char* o, iodata *io);
void dataRW(char* filename, const char* o, iodata::odata& odata);
double fitting(iodata& io, simdata &sd);



void Calculate_ThomasCoefficient(std::vector<double> X,coefficient_Thomas& coe,double deltaT);
void Calculate_ThomasCoefficient(std::vector<double> X, coefficient_Thomas& coe, double *factor, const simdata& sd);
void getConcentration(coefficient_Thomas& coe,std::vector<double> &C );
simdata gridScale(iodata::idata id);
simdata dimensionTrans(iodata::idata id);
simdata dimensionTrans(simdata sd);
void simulation(simdata& sd);

#endif // FUNCTION_H
