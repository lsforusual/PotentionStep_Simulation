#include "simdata.h"

simdata::simdata(iodata::idata id)
{
  t0=id.t0;
  t1=id.t1;
  DA=id.DA;
  DB=id.DB;
  rd=id.rd;
  A0=id.A0;
  B0=id.B0;
  d=DA/DB;
  time=id.time;
  current=id.current;
  potention=id.potention;
//  if(!id.potention.size())
//    for(unsigned int i=0; i<id.potention.size(); i++)
//      potention[i]=id.potention[i];

}
simdata::simdata()
{}

void simdata::setOdata(iodata &io)
{
  int N=this->current.size();
  io.od.current.assign(N,-1);
  io.od.potention.assign(N,-1);
  io.od.time.assign(N,-1);
  io.od.DA=this->DA;
  io.od.DB=this->DB;
  for(int i=0; i<N; i++)
    {
      io.od.current[i]=this->current[i];
      io.od.time[i]=this->time[i];
      io.od.potention[i]=this->potention[i];
    }

}
