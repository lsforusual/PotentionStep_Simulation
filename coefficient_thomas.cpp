#include "coefficient_thomas.h"


bool coefficient_Thomas::isempty(std::vector<double> &X){return X.empty();}
unsigned int coefficient_Thomas::size(std::vector<double> &X){return X.size();}

void coefficient_Thomas::setAlpha(double value){Alpha.push_back(value); }
void coefficient_Thomas::setBeta (double value){Beta.push_back(value);  }
void coefficient_Thomas::setGamma(double value){Gamma.push_back(value); }
void coefficient_Thomas::setDelta(double value){Delta.push_back(value); }

double coefficient_Thomas::getAlpha(int index){return Alpha[index];}
double coefficient_Thomas::getBeta(int index) {return Beta[index]; }
double coefficient_Thomas::getGamma(int index){return Gamma[index];}
double coefficient_Thomas::getG_mod(int index){return G_mod[index];}
double coefficient_Thomas::getDelta(int index){return Delta[index];}
double coefficient_Thomas::getD_mod(int index){return Delta_mod[index];}




void coefficient_Thomas::setG_mod()
{
  if(G_mod.empty())
    G_mod.push_back(Gamma[0]/Beta[0]);
  else
    G_mod[0]= Gamma[0]/Beta[0];

  for(unsigned int i=1; i<Alpha.size(); i++)
    G_mod.push_back(Gamma[i]/(Beta[i] - G_mod[i-1]*Alpha[i] ));
}


void coefficient_Thomas::setDelta_mod()
{
  if(Delta_mod.empty())
    Delta_mod.push_back(Delta[0]/Beta[0]);
  else
    Delta_mod[0]=Delta[0]/Beta[0];

  for(unsigned int i=1; i<Alpha.size(); i++)
    Delta_mod.push_back(( Delta[i] - Delta_mod[i-1] * Alpha[i] ) / ( Beta[i] - G_mod[i-1] * Alpha[i] ));
}

void coefficient_Thomas::setValue(std::vector<double>& X, double value,int index){X[index]=value;}
void coefficient_Thomas::setFirst(std::vector<double>& X, double value){X[0]=value;}
void coefficient_Thomas::setLast(std::vector<double>& X, double value){X[X.size()-1]=value;}

void coefficient_Thomas::asignValue(std::vector<double>& X, double value, int length){X.resize(length);X.assign(length,value);}
void coefficient_Thomas::asignValue(std::vector<double>& X, double value){X.assign(X.size(),value);}
void coefficient_Thomas::asignValue(const char['ALL'] ,double value,int length)
{
  /*
   * reset all parameter except DELTA
   */
      asignValue(alpha,value,length);
      asignValue(beta,value,length);
      asignValue(gamma,value,length);
      asignValue(g_mod,value,length);
      asignValue(d_mod,value,length);
}

