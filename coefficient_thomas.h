#ifndef COEFFICIENT_THOMAS_H
#define COEFFICIENT_THOMAS_H

#include <vector>


class coefficient_Thomas{
	
public:
  std::vector<double>& alpha= Alpha;
  std::vector<double>& beta= Beta;
  std::vector<double>& gamma= Gamma;
  std::vector<double>& delta= Delta;
  std::vector<double>& g_mod= G_mod;
  std::vector<double>& d_mod= Delta_mod;

  bool  isempty(std::vector<double> &X);
  unsigned int  size(std::vector<double> &X);

  void  setAlpha(double value);
  void  setBeta (double value);
  void  setGamma(double value);
  void  setDelta(double value);

  double  getAlpha(int index);
  double  getBeta(int index) ;
  double  getGamma(int index);
  double  getG_mod(int index);
  double  getDelta(int index);
  double  getD_mod(int index);
  void  setG_mod();
  void  setDelta_mod();

  void  setValue(std::vector<double>& X, double value,int index);
  void  setFirst(std::vector<double>& X, double value);
  void  setLast(std::vector<double>& X, double value);
  void  asignValue(std::vector<double>& X, double value, int length);
  void  asignValue(std::vector<double>& X, double value);
  void  asignValue(const char['ALL'] ,double value,int length);


private:
  std::vector<double> Alpha;
  std::vector<double> Beta;
  std::vector<double> Gamma;
  std::vector<double> Delta;

  std::vector<double> G_mod;
  std::vector<double> Delta_mod;
};

#endif // COEFFICIENT_THOMAS_H
