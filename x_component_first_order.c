/****************************************************************************
Slip_Thermal_BC.c
UDF – User Defined Function
For Slip velocity(Maxwell) + temperature jump (SmoluchowsKi)
*****************************************************************************/
#include "udf.h"  /*must be at the beginning of every UDF */
#include "sg.h"    /*must be at the beginning of this UDF */
#include "math.h"  /*must be at the beginning of this UDF */
#define not SCHEMEMFP /* shall MeanFreePath be scheme var.? */
#define SCHEMEUDRLXCOEFF   /* shall the under–relaxation Coefficient... */
#define SCHEMESpHR  /* shall SpHR ... */
#define SCHEMESIGMASQUARE /* shall sigma-square value... */
#define SCHEMEAMBPRESS  /* shall ambient pressure value... */
#define UNDERRLX    /* shall slip veloc be under-relaxated in Severe cases? */
#define WALLMOTION    /* shall wall-motion be regarded? */
#define sigma_square 1.35305239e-19   /* squared value of Sigma (molecule  diameter ) */
#define ambpress     101325  /* ambient  pressure */
#define TMAC         1.0   /* tangential momentum accomodation Coefficient */
#define ThAC         1.0   /* thermal accomodation Coefficient */
#define SpHR        1.4   /* specific heat ratio;  Air, Oxygen, Nitrogen */
#define UDRLXCOEFF  0.02   /* under-relaxation coefficient */
#define Boltzmann   1.3806505e-23   /* Boltzmann constant */
#define PI          3.14159265358979323846 /* number pi */
#define SQRT_2      1.41421356237309504880 /* sqrt(2) */
#define	knodson 	0.02
#define hydraulic_diameter 5.0e-5
#define Twall	320
/*
==========================================
    Velocity slip at wall boundaries, 
    a separate routine for every velocity coordinate 
    this is x-coordinate routine 
    with thermal creep term  
==========================================
*/
DEFINE_PROFILE(maxwell_slip_velocity_x,f_thread,index)
{
		face_t face;
		cell_t cell;
        Thread *c_thread;
        real slip, thcreep, dveloc;
        real normal_slip, tangential_slip, tangential_thcreep;
        real Coeff1[ND_ND], Coeff2[ND_ND],a;
        real u[ND_ND];
		real MeanFreePath;			//=6.8e-8;

	
		real y[ND_ND];
		real A[ND_ND];
		real dr0[ND_ND], es[ND_ND], ds, A_by_es;
		
		
		begin_f_loop(face,f_thread)
		{
		    F_CENTROID(y,face,f_thread);
			cell=F_C0(face,f_thread);
			c_thread=THREAD_T0(f_thread);	
			BOUNDARY_FACE_GEOMETRY( face,f_thread,A,ds,es,A_by_es,dr0) ;		
			ND_SET(u[0],u[1],u[2], C_U(cell,c_thread), C_V(cell,c_thread), C_W(cell,c_thread));
			a=NV_MAG(u);	
			ND_SET(F_UDMI(face,f_thread,0),   F_UDMI(face,f_thread,1),	F_UDMI(face,f_thread,2), NVD_DOT(u,1,0,0)/a ,	NVD_DOT(u,0,1,0)/a, 	NVD_DOT(u,0,0,1)/a  );
			ND_SET(F_UDMI(face,f_thread,3),   F_UDMI(face,f_thread,4),	F_UDMI(face,f_thread,5),  NVD_DOT(A,1,0,0)/A_by_es,	 NVD_DOT(A,0,1,0)/A_by_es,	NVD_DOT(A,0,0,1)/A_by_es) ;
		}
		end_f_loop(face,f_thread)
		
			// if (!Data_Valid_P())
			// {
			begin_f_loop(face,f_thread)
				{
		
				/* get cell and cell thread pointer */
                cell=F_C0(face,f_thread);
                c_thread=THREAD_T0(f_thread);
				
				/* compute mean free path at every position */

             MeanFreePath = Boltzmann * F_T(face,f_thread) / (sigma_square * (F_P(face,f_thread)+ambpress)* PI * SQRT_2 ) ;
	
			/* save the velocity coordinates into u[ND_ND] */
                ND_SET(u[0],u[1],u[2], F_U(face,f_thread), F_V(face,f_thread), F_W(face,f_thread)) ;
			/* save the transformation coefficients c_Mm:	c-11, c_21, c_31 into Coeff1[] */
                ND_SET(Coeff1[0],Coeff1[1],Coeff1[2], F_UDMI(face,f_thread,0), F_UDMI(face,f_thread,1), F_UDMI(face,f_thread,2)) ;
			/* save the transformation  coefficients c_Mm:	c_12, c_22, c_32 into Coeff2[]  */
                ND_SET(Coeff2[0],Coeff2[1],Coeff2[2], F_UDMI(face,f_thread,3), F_UDMI(face,f_thread,4), F_UDMI(face,f_thread,5)) ;
					
			/* evaluate the du/dl (tangential to surface) term in the local coord. System */
                tangential_slip= NVD_DOT(Coeff1, NV_DOT(Coeff1,C_U_G(cell,c_thread)), NV_DOT(Coeff1,C_V_G(cell,c_thread)), NV_DOT(Coeff1,C_W_G(cell,c_thread))) ;
					
			/* evaluate the du/dy (normal to surface ) term in the local coord. System */
                normal_slip=-1 * NVD_DOT(Coeff1, NV_DOT(Coeff2,C_U_G(cell,c_thread)), NV_DOT(Coeff2,C_V_G(cell,c_thread)), NV_DOT(Coeff2,C_W_G(cell,c_thread))) ;
					
			/* add theses values and multiply with MFP and TMAC */
                slip = ((2-TMAC)/TMAC) * MeanFreePath * (tangential_slip+normal_slip) ;
					
			/* evaluate the dT/dl (tangential to surface) term in the  local coord. System */
                tangential_thcreep=NV_DOT(Coeff1, C_T_G(cell,c_thread)) ;
                thcreep  = 0.75 * C_MU_L(cell,c_thread)/ (C_R(cell,c_thread)  *  C_T(cell,c_thread)) * tangential_thcreep ;
					

				dveloc = Coeff1[0]*  (slip + thcreep) ;	

                dveloc = (1-UDRLXCOEFF) *u[0] + UDRLXCOEFF * dveloc;
	
			/* boundary condition value is returned */
                F_PROFILE(face,f_thread,index) = dveloc;
				}
        end_f_loop(face,f_thread)

}

//A.3. Temperature Jump Boundary Condition Routine

/*
   =======================================
   Temperature change at wall boundaries
   =======================================
*/   
DEFINE_PROFILE(temperature_jump,f_thread,index)
{

        real MeanFreePath;		//=6.8e-8;
      
        face_t face;
        cell_t cell;
		Thread *c_thread;
		real Coeff2[ND_ND], a;
		real Prandtl, gamma, temp, normal;
		
        real u[ND_ND];
		real y[ND_ND];
		real A[ND_ND];
		real dr0[ND_ND], es[ND_ND], ds, A_by_es;
		
		
		begin_f_loop(face,f_thread)
		{
		    F_CENTROID(y,face,f_thread);
			cell=F_C0(face,f_thread);
			c_thread=THREAD_T0(f_thread);	
			BOUNDARY_FACE_GEOMETRY( face,f_thread,A,ds,es,A_by_es,dr0) ;		
			ND_SET(u[0],u[1],u[2], C_U(cell,c_thread), C_V(cell,c_thread), C_W(cell,c_thread));
			a=NV_MAG(u);	
			ND_SET(F_UDMI(face,f_thread,0),   F_UDMI(face,f_thread,1),	F_UDMI(face,f_thread,2), NVD_DOT(u,1,0,0)/a ,	NVD_DOT(u,0,1,0)/a, 	NVD_DOT(u,0,0,1)/a  );
			ND_SET(F_UDMI(face,f_thread,3),   F_UDMI(face,f_thread,4),	F_UDMI(face,f_thread,5),  NVD_DOT(A,1,0,0)/A_by_es,	 NVD_DOT(A,0,1,0)/A_by_es,	NVD_DOT(A,0,0,1)/A_by_es) ;
		}
		end_f_loop(face,f_thread)
		
/* call the routine that computes the transformation  coefficients, once per iteration should be enough */
       
/* loops over all faces in the thread passed in the DEFINE  macro argument */

        begin_f_loop(face,f_thread)
		{
		
		        cell=F_C0(face,f_thread);
				c_thread=THREAD_T0(f_thread);
				
                 MeanFreePath = Boltzmann * F_T(face,f_thread) / (sigma_square * (F_P(face,f_thread)+ambpress)*	PI * SQRT_2 ) ;
				
                Prandtl = (C_MU_L(cell,c_thread)*   C_CP(cell,c_thread))/C_K_L(cell,c_thread) ;
					
			    ND_SET(Coeff2[0],Coeff2[1],Coeff2[2], F_UDMI(face,f_thread,3),  F_UDMI(face,f_thread,4),  F_UDMI(face,f_thread,5)) ;
				/* evaluate the dT/dn (normal to surface) term in the  local coord. system */
                normal=NV_DOT(Coeff2,C_T_G(cell,c_thread));

                gamma=(2*SpHR)/(SpHR+1);
                temp=((2-ThAC)/ThAC) * gamma *  MeanFreePath/Prandtl * normal ;	
					
				 F_PROFILE(face,f_thread,index) =   F_T(face,f_thread)*(1-UDRLXCOEFF) + (Twall+temp)*UDRLXCOEFF ;

        }
        end_f_loop(f,thread)
}





