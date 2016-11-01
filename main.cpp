#include <iostream>
#include <cmath>
#include "coefficient_thomas.h"
#include "functions.h"
#include "iodata.h"
#include "simdata.h"

using namespace std;


int main(int argc, char *argv[])
{

  char filename[_MAX_FNAME];
  char file_path[_MAX_PATH];
  char dir[_MAX_DIR];
  char drive[_MAX_DRIVE];


  char* inputFilename;
  char* outputFilename;

  inputFilename="input";
  outputFilename="output";

  _splitpath(argv[0],drive,dir,filename,NULL);
  _makepath(file_path,drive,dir,inputFilename,".txt");

  cout << file_path << endl;
  //getchar();

  int rows = getRows(file_path)-1;
  //dataRW(io.id);
  //  io.id.DA=1.0E-8;
  //  io.id.DB=1.2E-8;
  //  io.id.t0=1;
  //  io.id.t1=2;
  //  io.id.time = io.id.seriesTime();
  //  io.id.theta0=0.3;
  //  io.id.theta1=-0.3;
  //  io.id.theta2=0.3;
  //io.id.current.resize(io.id.time.size(),);
  iodata io[rows];
  dataRW(file_path,"r",io);

  inputFilename="data";
  _makepath(file_path,drive,dir,inputFilename,".txt");


  simdata sd[rows] ;
  for(int i=0; i<rows; i++)
    {

      //io[i].id.readData(file_path);

      sd[i] = simdata(io[i].id);
      sd[i] = gridScale(io[i].id);

      simulation(sd[i]);
      //fitting(io[i],sd[i]);
      sd[i] = dimensionTrans(sd[i]);
      sd[i].setOdata(io[i]);


      _makepath(file_path,drive,dir,outputFilename,".txt");
      dataRW(file_path,"w",io[i].od);


    }
  //class iodata(subclass idata and odata) fitting();-> ssrFunc();
  //dataRW();




  return 0;
}
