#ifndef IODATA_H
#define IODATA_H
#include <vector>
#include <cmath>
#include <string>
using namespace std;
class iodata
{
public:
  iodata();

  double ssrFunc(std::vector<double> id, std::vector<double> od);
  double ssrFunc();

  double t0,t1;
  double DA,DB;
  double rd,A0,B0;

  class idata
  {
  public:
    double t0,t1;
    double DA,DB;
    double rd,A0,B0;
    double theta0, theta1, theta2;
    double T;

    idata();
    void setData(string str);

    std::vector<double> seriesTime();

    std::vector<double> time;
    std::vector<double> current;
    std::vector<double> potention;

    void readData(char* filename);


  };

  class odata
  {
  public:
    odata();
    string getPlotData();
    string getTitle();
    string getData(const char* o);
    string getData();

    std::vector<double> time;
    std::vector<double> current;
    std::vector<double> potention;
    double ssr;
    double DA,DB;
    friend class idata;

  };

  idata id;
  odata od;
};

#endif // IODATA_H
