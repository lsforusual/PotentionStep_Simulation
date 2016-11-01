#ifndef SIMDATA_H
#define SIMDATA_H
#include <iodata.h>
#include <vector>
class simdata: public iodata::idata
{
public:
  simdata(iodata::idata id);
  simdata();
  void setOdata(iodata &io);
  std::vector<double> r,z,t,A,B;
  std::vector<double> R,Z,dtau;
  std::vector<std::vector<std::vector<double> > >& a (std::vector<std::vector<double> >(std::vector<double> )) ;
  std::vector<std::vector<std::vector<double> > >& b (std::vector<std::vector<double> >(std::vector<double> )) ;

  double d;
  double tau0,tau1;
//  simdata dimensionTrans(simdata io, char[]);
};

#endif // SIMDATA_H
