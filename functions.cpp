
#include "functions.h"
double fitting(iodata& io,simdata &sd)
{

  int i=0;
  double step=0.5;
  double w=step;
  double ssr[]={0.0,0.0,0.0};
  iodata::idata id=io.id;
  sd=gridScale(id);
  std::vector<double> current=sd.current;
  std::vector<double> dtau=sd.dtau;

  sd.dtau=dtau;
  sd.current=current;
  simulation(sd);
  ssr[1]=io.ssrFunc(current,sd.current);

  sd.dtau.assign(sd.dtau.size(),sd.dtau.front()*(1+step/2));
  simulation(sd);
  ssr[2]=io.ssrFunc(current,sd.current);
  sd.dtau=dtau;

  sd.dtau.assign(sd.dtau.size(),sd.dtau.front() /(1+step/2));
  simulation(sd);
  ssr[0]=io.ssrFunc(current,sd.current);
  sd.dtau=dtau;


        while(i<70 )
          {
            if(ssr[0] <= ssr[1])
              {
                ssr[1]=ssr[0];
                sd.dtau.assign(sd.dtau.size(),sd.dtau.front() /(1+step));
                simulation(sd);
                ssr[0]=io.ssrFunc(current,sd.current);

                if(ssr[0] > ssr[1])
                  {
                    //sd.d*=(1+step);
                    sd.dtau.assign(sd.dtau.size(),sd.dtau.front() *(1+step));
                    step/=10;
                    simulation(sd);
                    ssr[0]=io.ssrFunc(current,sd.current);
                  }else if(ssr[0] == ssr[1]) break;
              }else if(ssr[2] <= ssr[1]){
                ssr[1]=ssr[2];
                sd.dtau.assign(sd.dtau.size(),sd.dtau.front() *(1+step));
                simulation(sd);
                ssr[2]=io.ssrFunc(current,sd.current);

                if(ssr[2] > ssr[1])
                  {
                    sd.dtau.assign(sd.dtau.size(),sd.dtau.front() /(1+step));
                    simulation(sd);
                    sd.setOdata(io);
                    ssr[2]=io.ssrFunc(current,sd.current);
                    step/=10;
                  }
              }else{break;}
            i++;
             cout<< i<< "'s fitting for DA, ssr is "<<ssr[1]<<endl;
          }
//  sd=gridScale(id);

//  simulation(sd);


//  //dtau = sd.dtau;


//  sd=dimensionTrans(sd);
//  sd.setOdata(io);
//  ssr[1]=io.ssrFunc();
//  double DA=id.DA;
//  double d=sd.d;


//    id.DA*=(1+step/2);
//    sd=gridScale(id);
//    simulation(sd);
//    sd=dimensionTrans(sd);
//    sd.setOdata(io);
//    ssr[2]=io.ssrFunc();
//    id.DA=DA;
//    step=w;


//    id.DA/=(1+step/2);
//    sd=gridScale(id);
//    simulation(sd);
//    sd=dimensionTrans(sd);
//    sd.setOdata(io);
//    ssr[0]=io.ssrFunc();
//    id.DA=DA;
//    step = w;

//      while(i<50 && id.DA>1E-8)
//        {
//          if(ssr[0] <= ssr[1])
//            {
//              ssr[1]=ssr[0];
//              id.DA/=(1+step);
//              //sd.d/=(1.0+step);
//              sd=gridScale(id);
//              simulation(sd);
//              sd=dimensionTrans(sd);
//              sd.setOdata(io);
//              ssr[0]=io.ssrFunc();


//              if(ssr[0] > ssr[1])
//                {
//                  //sd.d*=(1+step);
//                  id.DA*=(1+step);
//                  step/=10;
//                  sd=gridScale(id);
//                  simulation(sd);
//                  sd=dimensionTrans(sd);
//                  sd.setOdata(io);
//                  ssr[0]=io.ssrFunc();
//                }else if(ssr[0] == ssr[1]) break;
//            }else if(ssr[2] <= ssr[1]){
//              ssr[1]=ssr[2];
//              id.DA *= (1.0+step);
//              sd=gridScale(id);
//              simulation(sd);
//              sd=dimensionTrans(sd);
//              sd.setOdata(io);
//              ssr[2]=io.ssrFunc();
//              if(ssr[2] > ssr[1])
//                {
//                  id.DA/=(1+step);
//                  step/=10;
//                  sd=gridScale(id);
//                  simulation(sd);
//                  sd=dimensionTrans(sd);
//                  sd.setOdata(io);
//                  ssr[2]=io.ssrFunc();
//                }
//            }else{break;}
//          i++;
//           cout<< i<< "'s fitting for DA, ssr is "<<ssr[1]<<endl;
//        }
//      //fitting for C0
//      i=0;
//      step=w;
//      sd=gridScale(id);
//      simulation(sd);
//      sd=dimensionTrans(sd);
//      sd.setOdata(io);
//      ssr[1]=io.ssrFunc();
//      double A0=id.A0;


//        step =w;
//        id.A0*=(1+step/2);
//        sd=gridScale(id);
//        simulation(sd);
//        sd=dimensionTrans(sd);
//        sd.setOdata(io);
//        ssr[2]=io.ssrFunc();
//        id.A0=A0;

//        step =w;
//        id.A0/=(1+step/2);
//        sd=gridScale(id);
//        simulation(sd);
//        sd=dimensionTrans(sd);
//        sd.setOdata(io);
//        ssr[0]=io.ssrFunc();
//        id.A0=A0;

//          while(i<50 )
//            {
//              if(ssr[0] < ssr[1])
//                {
//                  ssr[1]=ssr[0];
//                  id.A0/=(1+step);
//                  //sd.d/=(1.0+step);
//                  sd=gridScale(id);
//                  simulation(sd);
//                  sd=dimensionTrans(sd);
//                  sd.setOdata(io);
//                  ssr[0]=io.ssrFunc();


//                  if(ssr[0] > ssr[1])
//                    {
//                      id.A0 *= (1+step);
//                      step/=10;
//                      sd=gridScale(id);
//                      simulation(sd);
//                      sd=dimensionTrans(sd);
//                      sd.setOdata(io);
//                      ssr[2]=io.ssrFunc();
//                    }else if(ssr[0] == ssr[1]) break;
//                }else if (ssr[2] < ssr[1]){
//                  ssr[1]=ssr[2];
//                  id.A0*=(1+step);
//                  sd=gridScale(id);
//                  simulation(sd);
//                  sd=dimensionTrans(sd);
//                  sd.setOdata(io);
//                  ssr[2]=io.ssrFunc();
//                  if(ssr[1] < ssr[2])
//                    {
//                      id.A0/=(1+step);
//                      step/=10;
//                      sd=gridScale(id);
//                      simulation(sd);
//                      sd=dimensionTrans(sd);
//                      sd.setOdata(io);
//                      ssr[2]=io.ssrFunc();
//                    }else{
//                      break;
//                    }
//                }else{
//                  break;
//                }
//              i++;
//               cout<< i<< "'s fitting for A0, ssr is "<<ssr[1]<<endl;
//            }


//  sd = simdata(io.id);
//  sd = dimensionTrans(io.id);
//  sd[i] = gridScale(io[i].id);
//  simulation(sd);
//  ssr[1] = io.ssrFunc(io.id.current,sd.current);
//  double d=sd.d;

//  sd.DA*=1.5;
//  simulation(sd);
//  ssr[2]=io.ssrFunc(io.id.current,sd.current);

//  sd.d = d;
//  sd.d /=1.5;
//  simulation(sd);
//  ssr[0]=io.ssrFunc(io.id.current,sd.current);

//  while(i<20)
//    {
//      if(ssr[0] < ssr[1])
//        {
//          ssr[0]=ssr[1];
//          io.od.ssr=ssr[0];
//          sd.d/=(1.0+step);
//          simulation(sd);
//          ssr[1]=io.ssrFunc(io.id.current,sd.current);
//          if(ssr[0] > ssr[1])
//            {
//              sd.d*=(1+step);
//              step/=10;
//            }else if(ssr[0] == ssr[1]) break;
//        }else{
//          ssr[0]=ssr[2];
//          io.od.ssr=ssr[0];
//          sd.d*=(1.0+step);
//          simulation(sd);
//          ssr[2]=io.ssrFunc(io.id.current,sd.current);
//          if(ssr[0] > ssr[2])
//            {
//              sd.d/=(1+step);
//              step/=10;
//            }
//        }
//      i++;
//       cout<< i<< "'s fitting, ssr is "<<ssr[0]<<endl;
//    }

  return io.od.ssr=ssr[1];
}
int getRows(char* filename)
{
  int Rows=0;
  ifstream file(filename);
  if(file.is_open())
    for(char str[256];!file.eof();Rows++)
      file.getline(str,256);
  else
    {
      cout << "Can not opening file !!!\n";
      system("pause");
      exit(0);
    }
  file.close();
  return Rows;

}

void writeDataPlot(char* filename, iodata::odata odata)
{
  //std::stringstream str;
  //str << odata.getPlotData();
  ofstream file(filename);
  file << odata.getPlotData() <<endl;
  file.close();
}

void dataRW(char* filename,const char* o,iodata* io)
{
  int i=0;
  if(*o=='r')
    {
      ifstream file(filename);

      if(!file.is_open())
        {cout << "Can not opening file !!!\n";system("pause");exit(1);}
      while(file.peek() != EOF)
        {
          char str[256];

          file.getline(str,256);
          if(i>0)
            io->id.setData(str);
          //return idata;
          i++;
        }
      file.close();
    }
}

void dataRW(char* filename, const char* o, iodata::odata& odata)
  {
      fstream file;

      if (*o=='w') {
          file.open(filename,ios::in);
          if(!file.is_open())
            {
              file.close();
              file.open(filename,ios::out);
              file << odata.getTitle() <<endl;
            }
          else
            {

//              file.close();
//              cout << "File Exist!! Will save data in append mode!!" <<endl;
//              file.open(filename,ios::app);
            }
        }else if (*o=='a'){
          file.open(filename,ios::app);
        }


      if(file.is_open())
        {
            file << odata.getPlotData() <<endl;
            //file << odata.getData() <<endl;
        }
      file.close();


}

void Calculate_ThomasCoefficient(std::vector<double> X,coefficient_Thomas& coe,double deltaT)
{
  /*To Use Thomas Algorithm properly
   * you MUST set the Thomas coefficients of index[0]
   * MANNUALY BEFORE
   * and index[N-1]
   * MANNUALY AFTER
   */


  //check size of input elements
  if(coe.size(coe.delta) != X.size())
    {
      std::cout << "Please check size of Delta AND X whether equal" <<std::endl;
      return;
    }else if(coe.isempty(coe.alpha) || coe.isempty(coe.beta) || coe.isempty(coe.gamma) || coe.isempty(coe.delta)){
      std::cout << "Please check if set the intital value on Alpha Beta Gamma or Delta" <<std::endl;
      return;
    }


  //Calculate coe :Alpha Beta Gamma Delata first;
  for(unsigned int i=1; i<X.size()-1;i++)
    {
      double deltaX_p = X[i+1] - X[i];
      double deltaX_m = X[i] - X[i-1];

      coe.setAlpha(-2  / (deltaX_m * deltaX_m + deltaX_p * deltaX_m));
      coe.setGamma(-2  / (deltaX_p * deltaX_p + deltaX_p * deltaX_m));
      coe.setBeta( - coe.getAlpha(i) - coe.getGamma(i) - 1/deltaT);
    }
  coe.setG_mod();
  coe.setDelta_mod();



}
void Calculate_ThomasCoefficient(std::vector<double> X,coefficient_Thomas& coe,double* factor,const simdata& sd)
{
  /*
   * To Use Thomas Algorithm properly
   * you MUST set the Thomas coefficients of index[0]
   * MANNUALY BEFORE
   * and index[N-1]
   * MANNUALY AFTER
   */


  //check size of input elements
  if(coe.size(coe.delta) != X.size())
    {
      std::cout << "Please check size of Delta AND X whether equal" <<std::endl;
      return;
    }else if(coe.isempty(coe.alpha) || coe.isempty(coe.beta) || coe.isempty(coe.gamma) || coe.isempty(coe.delta)){
      std::cout << "Please check if set the intital value on Alpha Beta Gamma or Delta" <<std::endl;
      return;
    }


  //Calculate coe :Alpha Beta Gamma Delata first;
  for(unsigned int i=1; i<X.size()-1;i++)
    {
      double deltaX_p = X[i+1] - X[i];
      double deltaX_m = X[i] - X[i-1];

      coe.setAlpha(factor[i]* sd.dtau[i] * -2/ (deltaX_m * deltaX_m + deltaX_p * deltaX_m));
      coe.setGamma(factor[i]* sd.dtau[i] * -2/ (deltaX_p * deltaX_p + deltaX_p * deltaX_m));
      coe.setBeta( - coe.getAlpha(i) - coe.getGamma(i) + 1);
    }

  coe.setAlpha(0.0);
  coe.setBeta(1.0);
  coe.setGamma(0.0);
  //coe.setLast(coe.delta,1.0);

  coe.setG_mod();
  coe.setDelta_mod();

}

void getConcentration(coefficient_Thomas& coe,std::vector<double> &C )
{
  if( coe.size(coe.g_mod) != C.size())
    {std::cout<<"length of C is not "
                "equal to coe.g_mod"<< std::endl;return;}
  int length = coe.size(coe.g_mod);
  C[length-1]=coe.getD_mod(length-1);
  for(int i=length-2; i>=0; i--)
    C[i]=coe.getD_mod(i)-coe.getG_mod(i)*C[i+1];
}

simdata gridScale(iodata::idata id)
{
  /*
   * grid way define:
   * grid r,z by uniform theta,Tau
   * t,A,B in a expand way
   */
  const double F=96500;
  const double R=8.31452;
  double FRT=F/R/id.T;

  simdata sd=dimensionTrans(id);
  int N=50;
  //vars following are all dimentionless var


  /*
       * tau stands time while dtau stands the time interval;
       * N is the uniform grid number at the beginning of potentional step
       * tau0 is the initial time, also means first potentional step start time
       * tau1 is the time when second potentional step start
       * tau1 = tau0 + internal(which is set to 2.0 as default)
       */
  sd.tau0=id.t0*id.DA/(id.rd*id.rd);
  sd.tau1=id.t1*id.DA/(id.rd*id.rd);
  double tau=0;
  sd.theta0= FRT *(id.theta0);
  sd.theta1= FRT *(id.theta1);
  sd.theta2= FRT *(id.theta2);
  /*
       * Time grid
       * normal grid before tau0
       * expand grid between tau0 and tau1
       * w is the expand coefficient
       */
  double w=0.03;
  double h=1;
  //double _time=0.0;
  //bool flag=sd.time.empty();
  bool flagt=sd.tau1<sd.tau0;
  sd.dtau.resize(0);
  int i;
  if(sd.potention.empty())
  {
    for( i=0; i<=N; i++)
      {
        sd.dtau.push_back(sd.tau0/N);
        sd.potention.push_back(sd.theta0);
        tau+=sd.dtau[i];
        //if(flag) sd.time.push_back(_time+=sd.dtau[i]);
      }

    if(flagt) sd.tau1=sd.tau0+4*sd.DA/(sd.rd*sd.rd);

    for( sd.dtau[i-1]=h ; tau<=sd.tau1; tau += sd.dtau[i],i++)
      {
        sd.dtau.push_back(sd.dtau[i-1]*(1+w));
        sd.potention.push_back(sd.theta1);
        //if(flag) sd.time.push_back(_time+=sd.dtau[i]);
      }
    if(!flagt)
      for (sd.dtau[i-1]=sd.tau0/N; tau<=(sd.tau1+sd.tau0); tau+= sd.dtau[i],i++)
        {
          sd.dtau.push_back(sd.dtau[i-1]*(1+w));
          sd.potention.push_back(sd.theta2);
          //if(flag) sd.time.push_back(_time+=sd.dtau[i]);
        }
  }else{
      for( i=1; i<id.time.size(); i++)
       sd.dtau.push_back(sd.time[i]-sd.time[i-1]);

//      if(sd.tau1<=sd.tau0) sd.tau1=sd.tau0+4*sd.DA/(sd.rd*sd.rd);

//      for( i=N; tau<sd.tau1; tau += sd.dtau[i],i++)
//          sd.dtau.push_back(sd.dtau[i-1]*(1+w));

//      if(sd.tau1>sd.tau0)
//        for (sd.dtau[i-1]=sd.tau0/N; tau<sd.tau1*2; tau+= sd.dtau[i],i++)
//          sd.dtau.push_back(sd.dtau[i-1]*(1+w));
    }
  /*
       * Spatial grid
       * move to simulation() in order to save time;
       */



  return sd;
}
simdata dimensionTrans(iodata::idata id)
{
  /*
   *  time,potention and current nead dimenless
   */
  const double F=96500;
  const double R=8.31452;
  double FRT=F/R/id.T;
  simdata sd(id);



  for(unsigned int i=0; i<id.time.size(); i++)
    {
      if(!id.time.empty())
        sd.time[i] = 4*id.time[i]*id.DA/(id.rd*id.rd);
      if(!id.current.empty())
        sd.current[i]=id.current[i]/(4*F*id.DA*id.A0);
      if(!id.potention.empty())
        sd.potention[i] = id.potention[i]*FRT;

    }
  if(id.time.size() != id.current.size())
    {std::cout<<"Please check input data! cols of time and current are not equal! " <<std::endl;}
  /*
       * theta0 is the initial potentional,
       * theta1 is the potentional after first potentional step,
       * theta2 is the seccond;
       */

  return sd;
}
simdata dimensionTrans(simdata sd)
{
  double tau=0;
  const double F=96500;
  const double R=8.31452;
  double FRT=F/R/sd.T;
  for(unsigned int i=0; i<sd.dtau.size(); i++)
    {
      if(!sd.R.empty())
      sd.r.push_back(sd.R[i]*sd.rd);
      if(!sd.Z.empty())
      sd.z.push_back(sd.Z[i]*sd.rd);
      if(!sd.time.empty())
      sd.time[i] *= (sd.rd*sd.rd/sd.DA);
      //sd.A.push_back(sd.A0*sd.a[i]);
      sd.current[i] *= (4*F*sd.DA*sd.A0);
      sd.potention[i] /= FRT;
    }
  return sd;
}
//void simulation(simdata& sd)
//{
//  const double Pi=M_PI;
//  double theta=0;
//  double Tau=0;
//  int N=500;//sd.dtau.size()+1;
//  double dTheta=1.0/N;
//  double dTau=1.0/N;
//  std::vector<std::vector<std::vector<double> > >
//      C(sd.dtau.size()+1,std::vector<std::vector<double> >(N,std::vector<double> (N)));
//  std::vector<double> a(N,0);
//  std::vector<double> sTau(N,0);

//  /* fix theta = 0  varries dTau first
//   *
//   * if theta=0, R varries while Z=0;
//   * if set Tau=0 as initial, concentration of A=0;
//   */
//  // the most easy way:

//  coefficient_Thomas &coe = *(new coefficient_Thomas());
//  //set initial value of delta

//  coe.asignValue("ALL",0,1);//NOT include delta
//  coe.asignValue(coe.delta,1,N);

//  coe.setFirst(coe.alpha,0.0);
//  coe.setFirst(coe.beta,1.0);
//  coe.setFirst(coe.gamma,0.0);
//  coe.setFirst(coe.delta,1);


//  /*
//   * set time count by dtau's(time's) size
//   * where it should be dtau.size()+1;
//   */
//  double time=0;
//  sd.current.resize(sd.dtau.size()+1,0);
//  /*
//   * Tau != tau,
//   * Tau stands an auxilatury var in the way of solve Diff Eqn;
//   */
//  for(unsigned int i=1; i<sTau.size(); i++)
//    sTau[i] = sTau[i-1] + dTau;


////std::vector<double> factor(coe.size(coe.delta),1.0);
//  double factor[N];
//  for(unsigned int k=0; k<sd.dtau.size(); time+=sd.dtau[k],k++)
//    {


//      theta=0;
//      for(unsigned int i=0; theta<1 ; theta+=dTheta,i++)
//        {
//          //        double delta[N];
//          //        delta[N]=1;

//          Tau=0;

//          if(time<sd.tau0)
//            {
//              coe.asignValue("ALL",0.0,1);

//              coe.setFirst(coe.alpha,0.0);
//              coe.setFirst(coe.beta,1.0);
//              coe.setFirst(coe.gamma,0.0);
//              coe.setFirst(coe.delta,1/(1+exp(-10)));
//             // coe.setFirst(coe.delta,1/(1+exp(-sd.potention.at(k))));
//              coe.setLast(coe.delta,1.0);

//              Calculate_ThomasCoefficient(sTau,coe,factor,sd);
//              getConcentration(coe,C[k][i]);


//              sd.current[k] += (C[k][i][1]-C[k][i][0])/(dTau) * dTheta;

//              std::cout<<"1k= " << k <<"\t i=" <<i <<std::endl;


//            }else if(time<sd.tau1 || time<sd.time.back() ){
//              /*
//               * First potention step
//               */
//              coe.asignValue("ALL",0.0,1);

//              coe.setFirst(coe.alpha,0.0);
//              coe.setFirst(coe.beta,1.0);
//              coe.setFirst(coe.gamma,0.0);

//              //coe.setFirst(coe.delta,1/(1+exp(0)));
//              coe.setFirst(coe.delta,1/(1+exp(sd.potention.at(k))));
//              coe.setLast(coe.delta,1);

//              Calculate_ThomasCoefficient(sTau,coe,factor,sd);
//              getConcentration(coe,C[k][i]);


//              sd.current[k] += (C[k][i][1]-C[k][i][0])/(dTau) * dTheta;

//              std::cout<<"2k= " << k <<"\t i=" <<i <<std::endl;

//            }else{
//              /*
//               * Seccond potention step
//               */
//              coe.asignValue("ALL",0.0,1);

//              coe.setFirst(coe.alpha,0.0);
//              coe.setFirst(coe.beta,1.0);
//              coe.setFirst(coe.gamma,0.0);

//              //coe.setFirst(coe.delta,1/(1+exp(-sd.potention.at(k))));
//              coe.setFirst(coe.delta,1);
//              coe.setLast(coe.delta,1);


//              Calculate_ThomasCoefficient(sTau,coe,factor,sd);
//              getConcentration(coe,C[k][i]);

//              sd.current[k] += (C[k][i][1]-C[k][i][0])/(dTau) * dTheta /(sd.d);

//              std::cout<<"3k= " << k <<"\t i=" <<i <<std::endl;

//            }

//          for(unsigned int j=0; Tau<1; Tau+=dTau,j++)
//            {
//              if(i>0 && j>0)
//                factor[j]=(4*cos(Pi/2*Tau)*cos(Pi/2*Tau)/Pi/Pi)/(theta*theta+tan(Pi/2*Tau)*tan(Pi/2*Tau));
//              //Get a series of factory to change EACH of alpah beta and delta
//              if(k>0)
//                coe.setValue(coe.delta,(coe.getAlpha(k)*C[k-1][i][j-1]+
//                    coe.getBeta(k)*C[k-1][i][j]+
//                    coe.getGamma(k)*C[k-1][i][j+1]),j);
//                //coe.setValue(coe.delta,C[k-1][i][j],j);
//              //set R & Z for each iteration
//              //sd.R.push_back(sqrt(1-theta*theta)*cos(Pi*Tau/2));
//              //sd.Z.push_back(theta*tan(Pi*Tau/2));
//            }





//          if(k==240 && i==337) system("pause");
//          //factor.~vector();
//          std::cout<<"what the hell!!! " <<std::endl;
//        }

//    }
//  //return sd;
//}
void simulation(simdata& sd)
{
  sd.current.resize(0);
  //Set number of treads to be used by OpenMP
  omp_set_num_threads(4);

  //Specify simulation parameters
  double theta_i = 20;
  double theta_v = -20;
  double sigma = 1000;

  double h0 = 1e-4;
  double omega = 1.08;
  double deltaTheta = 0.01;

  //Determine other parameters
  double deltaT = deltaTheta /sigma;
  double maxT = 2*fabs(theta_v - theta_i) / sigma;
  double maxR = 6* sqrt(maxT) +1;
  double maxZ = 6* sqrt(maxT);
  //int t =(int)(maxT/deltaT);
  int t= sd.dtau.size()+1;//number of timesteps;

  //Make Z grid;
  std::vector<double> Z;
  double h = h0;
  Z.push_back(0.0);
  while(Z.back() <= maxZ)
    {
      Z.push_back(Z.back()+h);
      h *= omega;
    }
  int m = Z.size(); // number of spacesteps(Z)

  //Make R grid
  std::vector<double> R;
  h=h0;

  R.push_back(0);
  while (R.back()<0.5)
    {
      R.push_back(R.back() + h);
      h *= omega;
    }
  R.back() = 0.5;

  for(int i=R.size() -2; i>=0; i--)
    {
      R.push_back(1- R[i]);
    }
  int n_e=R.size();//number spacestep over electrode

  h=h0;
  while(R.back() <= maxR)
    {
      R.push_back(R.back() +h);
      h *= omega;
    }
  int n= R.size(); //number of spacesteps(R);

  //Make concentration grids
  double** Ck = new double*[n];
  double** C_ = new double*[n];

  for(int i=0; i<n; ++i)
    {
      Ck[i] = new double[m];
      C_[i] = new double[m];
      //set initial concentration
      for(int j=0; j<m; j++)
        {
          Ck[i][j] = 1.0;
          C_[i][j] = 1.0;
        }
    }

  //Create and set Thomas coefficients
  std::vector<double> z_al(m,0.0),z_be(m,0.0),z_ga(m,0.0);
  std::vector<double> r_al(n,0.0),r_be(n,0.0),r_ga(n,0.0);
  std::vector<double> ga_modZ1(m,0.0), ga_modZ2(m,0.0);
  std::vector<double> ga_modR(n,0.0);

  for(int j=1; j<m-1; j++)
    {
      z_al[j] = 2.0 / ( (Z[j+1] - Z[j-1])*(Z[j] - Z[j-1]) );
      z_ga[j] = 2.0 / ( (Z[j+1] - Z[j-1])*(Z[j+1] - Z[j]) );
      z_be[j] = -z_al[j] - z_ga[j] - 2.0/sd.dtau[j];
    }

  for(int i=1; i<n-1; i++)
    {
      r_al[i] = ( 1.0 / (R[i+1]-R[i-1]) )
          * ( 2.0/(R[i] - R[i-1]) - 1.0/R[i] );
      r_ga[i] = ( 1.0 / (R[i+1]-R[i-1]) )
          * ( 2.0/(R[i+1] - R[i]) + 1.0/R[i] );
      r_be[i] = -r_al[i] - r_ga[i] - 2.0 / sd.dtau[i];
    }

  //Modify gamma coefficients for Z sweep
  ga_modZ1[0] =0; //electrode boundary condition
  for(int j=1; j<m-1; j++)
    {
      ga_modZ1[j] = z_ga[j] / (z_be[j] - ga_modZ1[j-1] * z_al[j]);
    }
  ga_modZ2[0] = -1; //no-flux boundary condition
  for(int j=1; j<m-1; j++)
    {
      ga_modZ2[j] = z_ga[j] / (z_be[j] - ga_modZ1[j-1] * z_al[j]);
    }

  //Modify gamma  coefficients for R sweep
  ga_modR[0]= -1; //boundary condition;
  for(int j=1; j<n-1; j++)
    {
      ga_modR[j] = r_ga[j] / (r_be[j] - ga_modR[j-1] * r_al[j]);
    }

  //Open file to output CV
  std::ofstream CV("CV_Output.txt");

  //BEGIN SIMULATION
  double Theta = theta_i;
  //	std::vector<double> Theta(t,10);
  //	for(int i=t/2; i<t; i++)
  //		Theta[i] = -10;
  for(int k=0; k<t; k++)
    {
     //k<t/2 ? Theta -= deltaTheta:Theta += deltaTheta;

//     Theta  = sd.potention[k];
      Theta = sd.potention[k];
      deltaT = sd.dtau[k];
      double flux[2] = {0.0, 0.0};
      //copy concentration grid
      for(int i=0; i<n; i++)
        for(int j =0; j<m; j++)
          C_[i][j] = Ck[i][j];

      //-----Z SWEEP--
#pragma omp parallel for
      for(int i=1; i<n-1; i++)
        {
          Ck[i][m-1] = 1.0;
          for(int j=1; j<m-1; j++)
            {
              Ck[i][j] = - C_[i-1][j] * r_al[i] - C_[i][j] * (-r_al[i] - r_ga[i]) - C_[i][j] * 2.0/deltaT- C_[i+1][j] * r_ga[i];
            }

          //Set surface deltas and pointer to gamma_mod
          std::vector<double>* ga_modZ;
          if(i < n_e){
              Ck[i][0]= 1.0 /(1.0+exp(-Theta));
              ga_modZ = &ga_modZ1;
            }else{
              Ck[i][0] = 0;
              ga_modZ = &ga_modZ2;
            }

          //Modify deltas
          for(int j=1; j<m-1; j++)
            Ck[i][j] = ( Ck[i][j] - Ck[i][j-1] * z_al[j] ) / ( z_be[j] - (*ga_modZ)[j-1] * z_al[j] );

          //solve by back substitutuio
          Ck[i][m-1] = 1.0;
          for(int j=m-2; j>=0; j--)
            Ck[i][j] = Ck[i][j]-(*ga_modZ)[j]*Ck[i][j+1];

        }
      //copy concentration grid
      for(int i=0; i<n; i++)
        for(int j =0; j<m; j++)
          C_[i][j] = Ck[i][j];
      //Output current

      for(int i=1; i<n_e; i++)
        {
          double J2 = (Ck[i][1] -Ck[i][0]) * R[i];
          double J1 = (Ck[i-1][1] -Ck[i-1][0]) * R[i-1];
          flux[0] -= (0.5/h0)*(J2+J1)*(R[i] - R[i-1]);
        }
//      if(k%2 == 1)
//        {
//          flux[k%2] -= (0.5/h0)*(J2+J1)*(R[i] - R[i-1]);
//          sd.current.push_back(0.5*(flux[0]+flux[1]));
//           CV << Theta <<"\t" << sd.current.back() << "\n";
//        }
 //
      std::cout << "interation iteratio" << std::endl;



		//--- R SWEEP---
#pragma omp parallel for
		for(int j=1; j<m-1; j++)
		{
			// set deltas
			Ck[0][j] = 0;
			Ck[n-1][j] = 1.0;
			for(int i=1; i<n-1; i++)
			{
				Ck[i][j] = - C_[i][j-1] * z_al[j]
				- C_[i][j] * (-z_al[j] - z_ga[j])
					- C_[i][j] * 2.0/deltaT
					- C_[i][j+1] * z_ga[j];
			}
			// modify deltas
			for(int i=1; i<n-1; i++)
			{
				Ck[i][j] = ( Ck[i][j] - Ck[i-1][j] * r_al[i] )
					/ ( r_be[i] - ga_modR[i-1] * r_al[i] );
			}
			//solve by back substitution
			for(int i=n-2; i>=0; i--) {
				Ck[i][j] = Ck[i][j] - ga_modR[i]*Ck[i+1][j];
			}
			std::cout << "j=" << j <<std::endl;
		}
		std::cout << "k=" << k <<std::endl;


		//-----Z SWEEP--

		for(int i=0; i<n; i++)
		  for(int j =0; j<m; j++)
		    C_[i][j] = Ck[i][j];

#pragma omp parallel for
                for(int i=1; i<n-1; i++)
                  {
                    Ck[i][m-1] = 1.0;
                    for(int j=1; j<m-1; j++)
                      {
                        Ck[i][j] = - C_[i-1][j] * r_al[i] - C_[i][j] * (-r_al[i] - r_ga[i]) - C_[i][j] * 2.0/deltaT- C_[i+1][j] * r_ga[i];
                      }

                    //Set surface deltas and pointer to gamma_mod
                    std::vector<double>* ga_modZ;
                    if(i < n_e){
                        Ck[i][0]= 1.0 /(1.0+exp(-Theta));
                        ga_modZ = &ga_modZ1;
                      }else{
                        Ck[i][0] = 0;
                        ga_modZ = &ga_modZ2;
                      }

                    //Modify deltas
                    for(int j=1; j<m-1; j++)
                      Ck[i][j] = ( Ck[i][j] - Ck[i][j-1] * z_al[j] ) / ( z_be[j] - (*ga_modZ)[j-1] * z_al[j] );

                    //solve by back substitutuio
                    Ck[i][m-1] = 1.0;
                    for(int j=m-2; j>=0; j--)
                      Ck[i][j] = Ck[i][j]-(*ga_modZ)[j]*Ck[i][j+1];

                  }
                //copy concentration grid
                for(int i=0; i<n; i++)
                  for(int j =0; j<m; j++)
                    C_[i][j] = Ck[i][j];
                //Output current
                for(int i=1; i<n_e; i++)
                  {
                    double J2 = (Ck[i][1] -Ck[i][0]) * R[i];
                    double J1 = (Ck[i-1][1] -Ck[i-1][0]) * R[i-1];
                    flux[1] -= (0.5/h0)*(J2+J1)*(R[i] - R[i-1]);
                  }
                //      if(k%2 == 1)
                //        {
                //          flux[k%2] -= (0.5/h0)*(J2+J1)*(R[i] - R[i-1]);
                //          sd.current.push_back(0.5*(flux[0]+flux[1]));
                //           CV << Theta <<"\t" << sd.current.back() << "\n";
                //        }



                //--- R SWEEP---
#pragma omp parallel for
                for(int j=1; j<m-1; j++)
                  {
                    // set deltas
                    Ck[0][j] = 0;
                    Ck[n-1][j] = 1.0;
                    for(int i=1; i<n-1; i++)
                      {
                        Ck[i][j] = - C_[i][j-1] * z_al[j]
                            - C_[i][j] * (-z_al[j] - z_ga[j])
                            - C_[i][j] * 2.0/deltaT
                            - C_[i][j+1] * z_ga[j];
                      }
                    // modify deltas
                    for(int i=1; i<n-1; i++)
                      {
                        Ck[i][j] = ( Ck[i][j] - Ck[i-1][j] * r_al[i] )
                            / ( r_be[i] - ga_modR[i-1] * r_al[i] );
                      }
                    //solve by back substitution
                    for(int i=n-2; i>=0; i--) {
                        Ck[i][j] = Ck[i][j] - ga_modR[i]*Ck[i+1][j];
                      }
                    std::cout << "j=" << j <<std::endl;
                  }
                std::cout << "k=" << k <<std::endl;


                sd.current.push_back(0.5*(flux[0]+flux[1]));
                CV << k <<"\t" << sd.current.back() << "\n";
                std::cout << "k=" << k <<std::endl;

    }

  CV.close();
  std::cout << "END" <<std::endl;

}
