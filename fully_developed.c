#include "udf.h"
DEFINE_PROFILE(inlet_vel, t, i  )
{

real H=5e-5, x1,y;
real x[ND_ND];
real inlet_vel=17;
real Kn=0.02;
face_t f;

begin_f_loop(f,t)
{
C_CENTROID(x,f,t);
y=x[1];
F_PROFILE(f, t, i)= 2*inlet_vel *( 1- pow( ((y)/(H/2)),2 ) + (4*Kn) ) /  ( 1+8*Kn );

}
end_f_loop(f,t)

}