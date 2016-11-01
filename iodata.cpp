#include "iodata.h"
#include "functions.h"
#include<iostream>
#include<sstream>
#include<fstream>

iodata::iodata()
{
  iodata::idata id = *(new idata);
  iodata::odata od = *(new odata);
}

double iodata::ssrFunc(std::vector<double> id, std::vector<double> od)
{
  double ssr=0.0;
  int N;
  id.size()<od.size() ? N=id.size():N=od.size();
  for(int i=0; i<N-1; i++)
    ssr += std::abs((id[i]-od[i])*(id[i]-od[i])/(od[i]+1));
//  ssr /= N;
  return ssr;
}
double iodata::ssrFunc()
{
  double ssr=ssrFunc(this->id.current,this->od.current);
  return ssr;
}

iodata::idata::idata()
  :t0(1.0),t1(3.0),DA(1E-7),DB(DA),
    rd(1E-5),A0(1),B0(0),T(298.15){}

void iodata::idata::setData(string str="0 0 0 0 0 0 0 0 0 0 0 0 0")
{
  std::stringstream ss(str);
  int i=0;
  while(i<11)
    {
          switch (i % 11)
            {
            case 0:
              ss >> this->t0;
              std::cout << "t0 =" << this->t0 <<";";
              break;
            case 1:
              ss >> this->t1;
              std::cout << "t1 = " << this->t1 << ";" ;
              std::cout << "\n";
              break;
            case 2:
              ss >> this->DA;
              std::cout << "DA = " << this->DA << ";" ;
              break;
            case 3:
              ss >> this->DB;
              std::cout << "DB = "  << this->DB  << ";" ;
              break;
              std::cout << "\n";
            case 4:
              ss >> this->rd;
              std::cout << "rd = " << this->rd << ";" ;
              break;
            case 5:
              ss >> this->A0;
              std::cout << "A0 = "  << this->A0  << ";" ;
              break;
            case 6:
              ss >> this->B0;
              std::cout << "B0 = "  << this->B0  << ";" ;
              break;
              std::cout << "\n";
            case 7:
              ss >> this->theta0;
              std::cout << "theta0 = " << this->theta0 << ";" ;
              break;
            case 8:
              ss >> this->theta1;
              std::cout << "theta1 = " << this->theta1 << ";" ;
              break;
            case 9:
              ss >> this->theta2;
              std::cout << "theta1 = "  << this->theta1  << ";" ;
              std::cout << std::endl;
              break;
            case 10:
              ss >> this->T;
              std::cout << "T = "  << this->T  << ";" ;
              std::cout << std::endl;
              break;
            default:
              std::cout << "Can not split the data correctly, please check the dataformat!!" <<std::endl;
              exit(0);
            }

          i++;
        }

}

std::vector<double> iodata::idata::seriesTime()
{
  double t,dt,w;
  int N=200;
  std::vector<double> time(0);
  for(t=0.0,dt=t0/N; t<t0; t+=dt)
    time.push_back(t);
  for(t=t0,dt=t0/N,w=0.1; t<t1;t+=(dt*=(1+w)))
    time.push_back(t);
  for(t=t0,dt=t0/N,w=0.1; t<t1;t+=(dt*=(1+w)))
    time.push_back(t);
  return time;
}

void iodata::idata::readData(char* filename)
{
  int rows,i=0;
  rows = getRows(filename);
  fstream file(filename);
  this->time.resize(rows);
  this->current.resize(rows);
  this->potention.resize(rows);

  if(!file.is_open())
    {cout << "Can not opening file !!!\n";system("pause");exit(1);}
  while(file.peek() != EOF)
    {
      char str[256];

      file.getline(str,256);

      stringstream ss(str);
      ss>>this->time[i];
      ss>>this->current[i];
      ss>>this->potention[i];
      //return idata;
      i++;
    }
  file.close();

}


iodata::odata::odata()
  :ssr(0),DA(1E-7),DB(DA){}

string iodata::odata::getPlotData()
{
  std::stringstream str;
  str.setf(ios::scientific);
  str.precision(9);
  for(unsigned int i=0 ;i<this->current.size();i++)
    if((this->time.size() == this->current.size()) && this->current.size() !=0 )
      str << this->time[i] << "," <<this->current[i] <<"\n";
  return str.str();
}
string iodata::odata::getTitle()
{
  std::stringstream str;
  str << "DA=" << this->DA << ",\t DB="<< this->DB << ",\t SSR="<<this->ssr <<";\n";
  return str.str();
}

string iodata::odata::getData()
{
  std::stringstream str;
  str.setf(ios::scientific);
  str.precision(9);
  str << this->DA;
  str << this->DB;
  str << this->ssr;

  return str.str();
}

string iodata::odata::getData(const char* o)
{
  if(*o)
    {
      std::stringstream str;
      str.setf(ios::scientific);
      str.precision(9);
      str <<this->getTitle();
      str << this->DA;
      str << "\t";
      str << this->DB;
      str << this->ssr;
    }
}
