#define EXPORT_CODE extern "C" __declspec(dllimport)
#define CONVENTION __stdcall
#include "CoolPropLib.h"
#include <windows.h>
#undef EXPORT_CODE
#undef CONVENTION

#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <windows.h>

struct inputcontrol
{
char  coat_material[100];
char  channel_material[100];
char  coolname[100];
double mass_rate;
double ini_coolant_temp;
double ini_coolant_pressure;

};
struct gasinput 
{
double position;
double ma;
double totaltemperature;
double totalpressure;
};
struct section 
{
double coat_thickness;
double hotwall_thickness;
double channel_thickness;
double coldwall_thickness;
double channel_width;
double rib_width;
int nummesh_1;
int nummesh_2;
};
//GH3128
double CK_Coef_1[10][3] ={
		{100,515.1,11.3},
		{200,512.8,12.56},
		{300,533.4,14.24},
		{400,518.7,15.49},
		{500,515.2,16.75},
		{600,538.9,18.42},
		{700,537  ,19.68},
		{800,618.2,21.35},
		{900,628.1,23.02},
		{950,651  ,23.86}		                       
};
//GH3536
double CK_Coef_2[9][3] ={
		{100,372.6,13.38},
		{200,389.4,17.97},
		{300,456.4,20.27},
		{400,427.1,22.40},
		{500,452.2,24.62},
		{600,464.7,26.79},
		{700,515,29.05},
		{800,535.9,31.04},
		{900,561,33.44}
  };
//TA15
double CK_Coef_3[9][3] ={
		{100,545,8.8},
		{200,587,10.2},
		{300,628,10.9},
		{400,670,12.2},
		{500,712,13.8},
		{600,755,15.1},
		{700,838,16.8},
		{800,880,18.0},
		{900,922,19.7}
  };
//06Cr19Ni10
double CK_Coef_4[8][3] ={
		{100,502,16.3},
		{200,502,17.6},
		{300,502,18.8},
		{400,502,20.5},
		{500,502,21.8},
		{600,502,23.5},
		{700,502,24.7},
		{800,502,28.5}			
                       };
//T2-金属所
double CK_Coef_5[7][3] ={
		{25,250.7,0.5659},
		{400,309.1,0.5442},
		{600,337.3,0.5179},
		{800,378.6,0.5714},
		{1200,418.8,0.639},
		{1400,446.3,0.7981},
		{1600,535.1,1.4624}		
                          };
//T3-过程所
double CK_Coef_6[5][3] ={
		{20,350,0.75},
		{200,447,0.868},
		{600,568,0.866},
		{1000,548,0.705},
		{1400,645,0.949}		
                       };

//T1-陈剑
double CK_Coef_7[2][3] ={
		{20,4500,0.429},
		{1800,4500,0.429}	
                       };
//F1-C/C-SiC
double CK_Coef_8[8][3] ={
		{25,670,12.98},
		{200,1016,16.61},
		{400,1275,18.36},
		{600,1420,18.65},
		{800,1450,18.17},
		{1000,1450,21.523},
		{1400,1450,23.166},
		{1800,1450,25.822}			
                       };



//function to calculate the Taw
double findTaw(double Tstatic,double Ttotal,double Twall)
{
double delta = 9999.0;
double mu_hotgas=1.5E-5;
double cp_hotgas=2800;
double k_hotgas=0.04;
double pr_hotgas;
double Taw;
double Tfold,Tfnew;
Tfold = (Twall+Ttotal)/2.0;
while (delta>=1E-04)
{
	   
	mu_hotgas=(1.0E-6)*1.458*( pow(Tfold,1.5) )/(Tfold+110.4);
    cp_hotgas=10*(0.021144)*( pow((0.001*Tfold),9) )+9*(-0.45165)*( pow((0.001*Tfold),8) ) 
         +8*(3.83286)*( pow((0.001*Tfold),7) )+7*(-15.8513)*( pow((0.001*Tfold),6) ) 
         +6*(26.75992)*( pow((0.001*Tfold),5) )+5*(29.9455575)*( pow((0.001*Tfold),4) ) 
         +4*(-220.448)*( pow((0.001*Tfold),3) )+3*(388.667)*( pow((0.001*Tfold),2) ) 
         +2*(-206.916)*( pow((0.001*Tfold),1) )+1043.76;    
    k_hotgas=(1.0E-3)*2.648151*( pow(Tfold,1.5) )/( Tfold+( 245.4*( pow(10,(-12/Tfold)) ) ) );
	pr_hotgas=mu_hotgas*cp_hotgas/k_hotgas;
Taw = ( pow(pr_hotgas,(1.0/3.0)) )*(Ttotal-Tstatic)+Tstatic;
Tfnew = Twall+0.23*(Tstatic-Twall)+0.19*(Taw-Twall);
delta = pow((Tfnew-Tfold),2.0);
Tfold = Tfnew;
}
return Taw;
};

//function to calculate the hg0
double hotgas(double Taw, double Ts,double Twall,double Ps,double Ug)
{
double T_ref,Re_ref,St_ref;
double mu_hotgas_ref=1.5E-5;
double cp_hotgas_ref=2800;
double k_hotgas_ref=0.04;
double pr_hotgas_ref;
double Rg=287;
double L_ch=1.5;
double hg0;
T_ref=0.5* (Twall+Ts)+0.22* (Taw-Ts);
	mu_hotgas_ref=(1.0E-6)*1.458*( pow(T_ref,1.5) )/(T_ref+110.4);
    cp_hotgas_ref=10*(0.021144)*( pow((0.001*T_ref),9) )+9*(-0.45165)*( pow((0.001*T_ref),8) ) 
         +8*(3.83286)*( pow((0.001*T_ref),7) )+7*(-15.8513)*( pow((0.001*T_ref),6) ) 
         +6*(26.75992)*( pow((0.001*T_ref),5) )+5*(29.9455575)*( pow((0.001*T_ref),4) ) 
         +4*(-220.448)*( pow((0.001*T_ref),3) )+3*(388.667)*( pow((0.001*T_ref),2) ) 
         +2*(-206.916)*( pow((0.001*T_ref),1) )+1043.76;    
    k_hotgas_ref=(1.0E-3)*2.648151*( pow(T_ref,1.5) )/( T_ref+( 245.4*( pow(10,(-12/T_ref)) ) ) );
	pr_hotgas_ref=mu_hotgas_ref*cp_hotgas_ref/k_hotgas_ref;
	
//calulate Re number under the reference temperature 
Re_ref=Ps/Rg/T_ref*Ug*L_ch/mu_hotgas_ref;
St_ref=0.0287/(( pow(Re_ref,0.2) )*( pow(pr_hotgas_ref,0.4)) );
hg0=Ps/Rg/T_ref*Ug*St_ref*cp_hotgas_ref;
return hg0;
};


//function to calculate convection  Nusselt number  between wall and fluid 
double Nusselt(double p_liquid,double T_liquid,double Re_liquid,double Pr_liquid,double Pr_at_wall,double viscosity_liquid,double viscosity_at_wall,char cooling_material[100])
{
double Nusselt_wf;
double p_cr,T_cr;
double a[4];
//////////////////////////////////////////////////////
if ( (strcmp(cooling_material,"n-Decane")==0) ){
p_cr=2.11E6;T_cr=617.7;
a[0]=0.027;a[1]=0.8;a[2]=1.0/3;a[3]=0.14;
//a[0]=0.000045;a[1]=1.4;a[2]=0.4;a[3]=0.1;
//a[0]=0.0065;a[1]=0.89;a[2]=0.4;a[3]=0.1;
//Sieder-Tate换热关联式(适用于油等）：Nusselt_wf=0.027*pow(Re_liquid,0.8)*pow(Pr_liquid,1.0/3)*pow(viscosity_liquid/viscosity_at_wall,0.14);
Nusselt_wf=a[0]*pow(Re_liquid,a[1])*pow(Pr_liquid,a[2])*pow(viscosity_liquid/viscosity_at_wall,a[3]);
}
//////////////////////////////////////////////////////
if ( (strcmp(cooling_material,"IF97::water")==0) ){
p_cr=22.1E6,T_cr=647.15;
a[0]=0.021;a[1]=0.8;a[2]=0.43;a[3]=0.25;
//米海耶夫换热关联式（最符合水的换热关联式）：Nusselt_wf=0.021*pow(Re_liquid,0.8)*pow(Pr_liquid,0.43)*pow(Pr_liquid/Pr_at_wall,0.25);
Nusselt_wf=a[0]*pow(Re_liquid,a[1])*pow(Pr_liquid,a[2])*pow(Pr_liquid/Pr_at_wall,a[3]);
}

return Nusselt_wf;
}

//function to calculate convection coeffient hl's temperature gradient dhl_dT
double hl_T(double p_liquid,double T_liquid,double mass_rate,double equi_diameter,double viscosity_at_wall,double Pr_at_wall,double cp_liquid,double dcp_dT,double conduction_liquid,double dk_dT,double viscosity_liquid,double dmu_dT,char cooling_material[100])
{
double dhl_dT,p_cr,T_cr,Area;
double a[4],pi=3.1415926;

Area=pi*equi_diameter*equi_diameter/4.0;
//////////////////////////////////////////////////////
if ( (strcmp(cooling_material,"n-Decane")==0) ){
p_cr=2.11E6;T_cr=617.7;
a[0]=0.027;a[1]=0.8;a[2]=1.0/3;a[3]=0.14;
//a[0]=0.000045;a[1]=1.4;a[2]=0.4;a[3]=0.1;
//a[0]=0.0065;a[1]=0.89;a[2]=0.4;a[3]=0.1;
dhl_dT=(1.0-a[2])*pow(conduction_liquid,-a[2])* dk_dT*pow(cp_liquid,a[2])*pow(viscosity_liquid,-a[1]+a[2]+a[3]);
dhl_dT=dhl_dT+pow(conduction_liquid,1.0-a[2])*a[2]*pow(cp_liquid,a[2]-1.0)*dcp_dT*pow(viscosity_liquid,-a[1]+a[2]+a[3]);
dhl_dT=dhl_dT+pow(conduction_liquid,1.0-a[2])*pow(cp_liquid,a[2])*(-a[1]+a[2]+a[3])*pow(viscosity_liquid,-a[1]+a[2]+a[3]-1.0)*dmu_dT;
dhl_dT=dhl_dT*a[0]*pow(mass_rate,a[1])*pow(equi_diameter,a[1]-1)/pow(Area,a[1])/pow(viscosity_at_wall,a[3]);
}
//////////////////////////////////////////////////////
if ( (strcmp(cooling_material,"IF97::water")==0) ){
p_cr=22.1E6,T_cr=647.15;
a[0]=0.021;a[1]=0.8;a[2]=0.43;a[3]=0.25;
dhl_dT=(1.0-a[2]-a[3])*pow(conduction_liquid,-a[2]-a[3])* dk_dT*pow(cp_liquid,a[2]+a[3])*pow(viscosity_liquid,-a[1]+a[2]+a[3]);
dhl_dT=dhl_dT+pow(conduction_liquid,1.0-a[2]-a[3])*(a[2]+a[3])*pow(cp_liquid,a[2]+a[3]-1.0)*dcp_dT*pow(viscosity_liquid,-a[1]+a[2]+a[3]);
dhl_dT=dhl_dT+pow(conduction_liquid,1.0-a[2]-a[3])*pow(cp_liquid,a[2]+a[3])*(-a[1]+a[2]+a[3])*pow(viscosity_liquid,-a[1]+a[2]+a[3]-1.0)*dmu_dT;
dhl_dT=dhl_dT*a[0]*pow(mass_rate,a[1])*pow(equi_diameter,a[1]-1)/pow(Area,a[1])/pow(Pr_at_wall,a[3]);

}
return dhl_dT;
}

//function to calculate convection coeffient hl's wall temperature gradient dhl_dTw
double hl_Tw(double p_liquid,double T_liquid,double mass_rate,double equi_diameter,double cp_liquid,double conduction_liquid,double viscosity_liquid,double viscosity_at_wall,double Pr_at_wall,double dmu_wall_dTw,double dpr_wall_dTw,char cooling_material[100])
{
double dhl_dTw,p_cr,T_cr,Area;
double a[4],pi=3.1415926;

Area=pi*equi_diameter*equi_diameter/4.0;
//////////////////////////////////////////////////////
if ( (strcmp(cooling_material,"n-Decane")==0) ){
p_cr=2.11E6;T_cr=617.7;
a[0]=0.027;a[1]=0.8;a[2]=1.0/3;a[3]=0.14;
dhl_dTw=(-a[3])*pow(viscosity_at_wall,-a[3]-1.0)*dmu_wall_dTw;
dhl_dTw=dhl_dTw*a[0]*pow(mass_rate,a[1])*pow(equi_diameter,a[1]-1)/pow(Area,a[1])*pow(conduction_liquid,1.0-a[2])*pow(cp_liquid,a[2])*pow(viscosity_liquid,-a[1]+a[2]+a[3]);
}
//////////////////////////////////////////////////////
if ( (strcmp(cooling_material,"IF97::water")==0) ){
p_cr=22.1E6,T_cr=647.15;
a[0]=0.021;a[1]=0.8;a[2]=0.43;a[3]=0.25;
dhl_dTw=(-a[3])*pow(Pr_at_wall,-a[3]-1.0)*dpr_wall_dTw;
dhl_dTw=dhl_dTw*a[0]*pow(mass_rate,a[1])*pow(equi_diameter,a[1]-1)/pow(Area,a[1])*pow(conduction_liquid,1.0-a[2]-a[3])*pow(cp_liquid,a[2]+a[3])*pow(viscosity_liquid,-a[1]+a[2]+a[3]);
}


return dhl_dTw;
}

//function to calculate Re number
double Re(double mass_rate,double equi_diameter,double area,double viscosity_liquid)
{
double Re_number;
Re_number=mass_rate*equi_diameter/area/viscosity_liquid;
return Re_number;
}

//function to calculate rib effect
double rib_effect(double channelwidth,double ribwidth,double channelthickness,double hl, double kswall)
{
double ribeffect;
double ef1;
double ef2;
ef1=sqrt(2*hl/kswall*ribwidth)* channelthickness/ribwidth;
ef2=( exp(ef1)-exp(-ef1) )/( exp(ef1)+exp(-ef1) );
ribeffect=channelwidth/(channelwidth+ribwidth)+2*channelthickness/(channelwidth+ribwidth)*ef2/ef1;
return ribeffect;
}

//function to calculate materials' conductivity
double CondMat(double T,char solid_materal[100])
{
int k,l1,l2,m;
double CK_Coef[20][3]={0};
/////////////////////////////////////////
if ( (strcmp(solid_materal,"GH3128")==0) )
{
m=10;
for (l1=0;l1<m;l1++)
{
	for (l2=0;l2<3;l2++)
	{
	CK_Coef[l1][l2]=CK_Coef_1[l1][l2];
	}
}
}
/////////////////////////////////////////////
if ( (strcmp(solid_materal,"GH3536")==0) )
{
m=9;
for (l1=0;l1<m;l1++)
{
	for (l2=0;l2<3;l2++)
	{
	CK_Coef[l1][l2]=CK_Coef_2[l1][l2];
	}
}
}
/////////////////////////////////////////////
if ( (strcmp(solid_materal,"TA15")==0) )
{
m=9;
for (l1=0;l1<m;l1++)
{
	for (l2=0;l2<3;l2++)
	{
	CK_Coef[l1][l2]=CK_Coef_3[l1][l2];
	}
}
}
////////////////////////////////////////////
if ( (strcmp(solid_materal,"06Cr19Ni10")==0) )
{
m=8;
for (l1=0;l1<m;l1++)
{
	for (l2=0;l2<3;l2++)
	{
	CK_Coef[l1][l2]=CK_Coef_4[l1][l2];
	}
}
}
//////////////////////////////////////////
if ( (strcmp(solid_materal,"T2")==0) )
{
m=7;
for (l1=0;l1<m;l1++)
{
	for (l2=0;l2<3;l2++)
	{
	CK_Coef[l1][l2]=CK_Coef_5[l1][l2];
	}
}
}
/////////////////////////////////////////
if ( (strcmp(solid_materal,"T3")==0) )
{
m=5;
for (l1=0;l1<m;l1++)
{
	for (l2=0;l2<3;l2++)
	{
	CK_Coef[l1][l2]=CK_Coef_6[l1][l2];
	}
}
}
/////////////////////////////////////////
if ( (strcmp(solid_materal,"T1")==0) )
{
m=2;
for (l1=0;l1<m;l1++)
{
	for (l2=0;l2<3;l2++)
	{
	CK_Coef[l1][l2]=CK_Coef_7[l1][l2];
	}
}
}
/////////////////////////////////////////
if ( (strcmp(solid_materal,"F1")==0) )
{
m=8;
for (l1=0;l1<m;l1++)
{
	for (l2=0;l2<3;l2++)
	{
	CK_Coef[l1][l2]=CK_Coef_8[l1][l2];
	}
}
}
/////////////////////////////////////////

T=T-273.15;
   if (T<CK_Coef[0][0])
   {
  return CK_Coef[0][2];
   }
   if (T>CK_Coef[m-1][0])
   {
  return CK_Coef[m-1][2];
   }
   for (k=0;k<m-1;k++)
   {
    if (T>CK_Coef[k][0]&&T<CK_Coef[k+1][0])
   {
  return CK_Coef[k][2]+(T-CK_Coef[k][0])*(CK_Coef[k+1][2]-CK_Coef[k][2])/(CK_Coef[k+1][0]-CK_Coef[k][0]);
   }
   }

}

//function to calculate friction coeffient f
double f(double Ra,double equi_diameter,double Re_eff )
{
double f_old=1.0E-6;
double f_new=1.0E-6;
double delta=1;

while (delta>=1E-08)
{
f_new=	-2.0*log(Ra/3.7065/equi_diameter+2.5226/Re_eff/sqrt(f_old));
f_new=pow(1.0/f_new,2.0);
delta = pow((f_new-f_old),2.0);
f_old = f_new;
}
return f_new;
}

//function to calculate chemical reaction rate k according to Arrhenius equation (fuction works when 1050K >T>750K)
double k_Arrhenius(double p,double T)
{
double K_Arrhenius_Coef[3][3] ={
		{3.0E6,251.88E3,1.8*pow(10,15.2)},
		{4.5E6,251.88E3,1.8*pow(10,15.4)},
		{6.0E6,251.88E3,2.0*pow(10,15.5)}		
                          };
double k_Arrhenius_temp=0.0;
double temp;

if (750<=T) {
	if (p<=K_Arrhenius_Coef[0][0]){
	k_Arrhenius_temp=K_Arrhenius_Coef[0][2]*exp(-K_Arrhenius_Coef[0][1]/8.314/T);	
	}

	if (p>=K_Arrhenius_Coef[2][0]){
	k_Arrhenius_temp=K_Arrhenius_Coef[2][2]*exp(-K_Arrhenius_Coef[2][1]/8.314/T);	
	}

	if ( (p>K_Arrhenius_Coef[0][0])&&(p<=K_Arrhenius_Coef[1][0]) ){
    temp=K_Arrhenius_Coef[0][2]+(p-K_Arrhenius_Coef[0][0])*(K_Arrhenius_Coef[1][2]-K_Arrhenius_Coef[0][2])/(K_Arrhenius_Coef[1][0]-K_Arrhenius_Coef[0][0]);
	k_Arrhenius_temp=temp*exp(-K_Arrhenius_Coef[1][1]/8.314/T);	
	}

	if ( (p>K_Arrhenius_Coef[1][0])&&(p<K_Arrhenius_Coef[2][0]) ){
    temp=K_Arrhenius_Coef[1][2]+(p-K_Arrhenius_Coef[1][0])*(K_Arrhenius_Coef[2][2]-K_Arrhenius_Coef[1][2])/(K_Arrhenius_Coef[2][0]-K_Arrhenius_Coef[1][0]);
	k_Arrhenius_temp=temp*exp(-K_Arrhenius_Coef[2][1]/8.314/T);	
	}
           }
 return k_Arrhenius_temp; 
}

//fuction to calculate chemical heat sink according to temprature and pressure
double CHS(double p,double T)
{
double temp,temp1,temp2;
int i;
//p=3.0E6
double CHS_Coef_1[7][2] ={
		{750,0.0},
		{840,0.018E6},
		{900,0.1786E6},
		{960,0.786E6},
		{990,1.161E6},
		{1020,1.411E6},
		{1050,1.61E6}
};
//p=4.5E6
double CHS_Coef_2[7][2] ={
	{750,0},
	{840,0.018E6},
	{900,0.232E6},
	{960,0.875E6},
	{990,1.1964E6},
	{1020,1.411E6},
	{1050,1.518E6}
};

//p=6.0E6
double CHS_Coef_3[7][2] ={
	{750,0},
	{840,0.018E6},
	{900,0.2857E6},
	{960,0.875E6},
	{990,1.125E6},
	{1020,1.286E6},
	{1050,1.375E6}            
};
	if (p<=3.0E6){
		if(T<=CHS_Coef_1[0][0]){
		temp=0;
		}
		if(T>=CHS_Coef_1[6][0]){
		temp=CHS_Coef_1[6][1];
		}
		for (i=0;i<6;i++){
			if ( (T>CHS_Coef_1[i][0])&&(T<=CHS_Coef_1[i+1][0]) ){
			temp=CHS_Coef_1[i][1]+(T-CHS_Coef_1[i][0])*(CHS_Coef_1[i+1][1]-CHS_Coef_1[i][1])/(CHS_Coef_1[i+1][0]-CHS_Coef_1[i][0]);
			}
		}	
	}

	if (p>=6.0E6){
		if(T<=CHS_Coef_3[0][0]){
		temp=0;
		}
		if(T>=CHS_Coef_3[6][0]){
		temp=CHS_Coef_3[6][1];
		}
		for (i=0;i<6;i++){
			if ( (T>CHS_Coef_3[i][0])&&(T<=CHS_Coef_3[i+1][0]) ){
			temp=CHS_Coef_3[i][1]+(T-CHS_Coef_3[i][0])*(CHS_Coef_3[i+1][1]-CHS_Coef_3[i][1])/(CHS_Coef_3[i+1][0]-CHS_Coef_3[i][0]);
			}
		}	
	}

	if ( (p<=4.5E6)&&(p>3.0E6) ){
		if(T<=CHS_Coef_1[0][0]){
		temp=0;
		}
		if(T>=CHS_Coef_1[6][0]){
		temp=CHS_Coef_1[6][1]+(p-3.0E6)*(CHS_Coef_2[6][1]-CHS_Coef_1[6][1])/(4.5E6-3.0E6);
		}
		for (i=0;i<6;i++){
			if ( (T>CHS_Coef_1[i][0])&&(T<=CHS_Coef_1[i+1][0]) ){
			temp1=CHS_Coef_1[i][1]+(T-CHS_Coef_1[i][0])*(CHS_Coef_1[i+1][1]-CHS_Coef_1[i][1])/(CHS_Coef_1[i+1][0]-CHS_Coef_1[i][0]);
			temp2=CHS_Coef_2[i][1]+(T-CHS_Coef_2[i][0])*(CHS_Coef_2[i+1][1]-CHS_Coef_2[i][1])/(CHS_Coef_2[i+1][0]-CHS_Coef_2[i][0]);
			temp=temp1+(p-3.0E6)*(temp2-temp1)/(4.5E6-3.0E6);
			}
		}	
	}

	if ( (p<6.0E6)&&(p>4.5E6) ){
		if(T<=CHS_Coef_2[0][0]){
		temp=0;
		}
		if(T>=CHS_Coef_2[6][0]){
		temp=CHS_Coef_2[6][1]+(p-4.5E6)*(CHS_Coef_3[6][1]-CHS_Coef_2[6][1])/(6.0E6-4.5E6);
		}
		for (i=0;i<6;i++){
			if ( (T>CHS_Coef_2[i][0])&&(T<=CHS_Coef_2[i+1][0]) ){
			temp1=CHS_Coef_2[i][1]+(T-CHS_Coef_2[i][0])*(CHS_Coef_2[i+1][1]-CHS_Coef_2[i][1])/(CHS_Coef_2[i+1][0]-CHS_Coef_2[i][0]);
			temp2=CHS_Coef_3[i][1]+(T-CHS_Coef_3[i][0])*(CHS_Coef_3[i+1][1]-CHS_Coef_3[i][1])/(CHS_Coef_3[i+1][0]-CHS_Coef_3[i][0]);
			temp=temp1+(p-4.5E6)*(temp2-temp1)/(6.0E6-4.5E6);
			}
		}	
	}
return temp;
}

int main ( )
{	
	char coating_material[100];
	char channel_material[100];
	char cooling_material[100];
	struct inputcontrol  ictr;
	struct gasinput ginp[1000];
	struct section  sinp[1000];	
	FILE *in_file,*gfile,*sfile,*outfile;	
	double Tw[1000][100],T[1000],p[1000];	
	int i,j,n=0;
	double ave_ma,ave_totaltemperature,ave_totalpressure;
	double Taw,Ts,Ug,ps;
	double Rg=287;
	double gam=1.4;
	double pi=3.1415926;
	double deltax,deltay,delta_T,delta_Tw,delta_err,delta=-1E-3;

	//dimensions
	double area_old,area_new;
	double perimeter_old,perimeter_new;
	double equi_diameter_old,equi_diameter_new;
	//Properties
	double rho_old,rho_new,rho_wall_new;
	double H_old,H_new;
	double mu_new,mu_wall_new;
	double k_new,k_wall_new;
	double cp_new,cp_wall_new;
	double pr_new,pr_wall_new;
	double Re_new;
    double h_conv;
	//friction
	double Ra=3.2E-6,friction;
	double momentum_pressure_loss,friction_pressure_loss,ec_pressure_loss;
	//expansion or contraction coeffient
	double e_c_coeffient;

	//N-L
	double F1,F2;
	double dF1_dT,dF1_dTw,dF2_dT,dF2_dTw,detF;
	//gradient to T
	double dH_dT,dH_dp,dp_dT,drho_dT,dhl_dT,dcp_dT,dk_dT,dmu_dT;
	//gradient to Tw
	double dp_dTw,dhl_dTw,dmu_wall_dTw,dpr_wall_dTw;
	   //Properties vary with Temperature
	double H_new_1,rho_new_1,mu_new_1,k_new_1,cp_new_1;
	double mu_wall_new_1,k_wall_new_1,cp_wall_new_1,pr_wall_new_1;
	   //Properties vary with Pressure
	double H_new_2;	
	
	double hg0,q;
	double ks;

	double relax_number=0.7;
	double sum_H=0;

	//chemical
	double Z_conversion[1000];
	double Z_conversion_new_1,Z_conversion_new_2;
	double coff=0;
	double dqe_dT=0,qe_old=0,qe_new=0,qe_new_1=0,qe_new_2=0;
	double k_A_old=0,k_A_new=0,k_A__new_1=0,k_A__new_2=0;
	double he_old=0,he_new=0,he_new_1=0,he_new_2=0;
	//stay time
	double delta_t[1000];
	int inter_k=0;
	
//read from input_control.txt
in_file = fopen (".\\input_control.txt", "r") ; /* open the file */
if (in_file == NULL) /* if fopen() returns NULL to in_file */ 
{
printf("\nThe file cannot be opened.");
printf("\nPlease check that the file currently exists.");
}
else
printf("\n(1/3)The input control txt has been successfully opened for reading.");
while (!feof(in_file))
{
fscanf(in_file,"%s %s %s %lf %lf %lf",&ictr.coat_material,&ictr.channel_material,&ictr.coolname,&ictr.mass_rate,&ictr.ini_coolant_temp,&ictr.ini_coolant_pressure);
}
fclose(in_file);
// change into SI units
ictr.mass_rate=0.001*ictr.mass_rate;

//read from gas_input.txt
gfile = fopen (".\\gas_input.txt", "r") ; /* open the file */
if (gfile == NULL) /* if fopen() returns NULL to in_file */
{
printf("\nThe file cannot be opened.");
printf("\nPlease check that the file currently exists.");
}
else
printf("\n(2/3)The gas input txt has been successfully opened for reading.");
while (!feof(gfile))
{
//                                       position,         ma,         totaltemperature,         totalpressure
fscanf(gfile,"%lf %lf %lf %lf ",&ginp[n].position,&ginp[n].ma,&ginp[n].totaltemperature,&ginp[n].totalpressure);
n=n+1;
}
fclose(gfile);
n=0;
//read from section.txt
sfile = fopen (".\\section.txt", "r") ; /* open the file */
if (sfile == NULL) /* if fopen() returns NULL to in_file */
{
printf("\nThe file cannot be opened.");
printf("\nPlease check that the file currently exists.");
}
else
printf("\n(3/3)The section txt has been successfully opened for reading.");
while (!feof(sfile))
{
//                                               coat_thickness(mm),     hotwall_thickness(mm),     channel_thickness(mm),     coldwall_thickness(mm),    channel_width(mm),      rib_width(mm),       number of mesh(n1), 
fscanf(sfile,"%lf %lf %lf %lf %lf %lf %d %d",&sinp[n].coat_thickness,&sinp[n].hotwall_thickness,&sinp[n].channel_thickness,&sinp[n].coldwall_thickness,&sinp[n].channel_width,&sinp[n].rib_width,&sinp[n].nummesh_1,&sinp[n].nummesh_2);
// change into SI units
sinp[n].coat_thickness=0.001*sinp[n].coat_thickness;
sinp[n].hotwall_thickness=0.001*sinp[n].hotwall_thickness;
sinp[n].channel_thickness=0.001*sinp[n].channel_thickness;
sinp[n].coldwall_thickness=0.001*sinp[n].coldwall_thickness;
sinp[n].channel_width=0.001*sinp[n].channel_width;
sinp[n].rib_width=0.001*sinp[n].rib_width;
n=n+1;
}
fclose(sfile);
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//compare and diside to use which cooling material ,currently the  materials are those in  open source  "CoolProp"  library(122 kinds of pure substances) .
if ( (strcmp(ictr.coolname,"RP-3")==0)|| (strcmp(ictr.coolname,"rp-3")==0) )
{
//strcpy(cooling_material,"n-Decane[0.6]&n-Dodecane[0.4]");
strcpy(cooling_material,"n-Decane");
}

if ( (strcmp(ictr.coolname,"Water")==0)|| (strcmp(ictr.coolname,"water")==0)|| (strcmp(ictr.coolname,"WATER")==0) )
{
strcpy(cooling_material,"IF97::water");
}
//compare and  diside to use which channel material ,currently the  materials are GH3128/GH3536/TA15/06Cr19Ni10/T2/T3/T1
if ( (strcmp(ictr.channel_material,"J9")==0)|| (strcmp(ictr.channel_material,"j9")==0) )
{
strcpy(channel_material,"GH3128");
}

if ( (strcmp(ictr.channel_material,"J8")==0)|| (strcmp(ictr.channel_material,"j8")==0) )
{
strcpy(channel_material,"GH3536");
}

if ( (strcmp(ictr.channel_material,"J6")==0)|| (strcmp(ictr.channel_material,"j6")==0) )
{
strcpy(channel_material,"TA15");
}
if ( (strcmp(ictr.channel_material,"J2")==0)|| (strcmp(ictr.channel_material,"j2")==0) )
{
strcpy(channel_material,"06Cr19Ni10");
}

//F1-C/C-SiC
if ( (strcmp(ictr.channel_material,"F1")==0)|| (strcmp(ictr.channel_material,"f1")==0) )
{
strcpy(channel_material,"F1");
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//T3-过程所
if ( (strcmp(ictr.coat_material,"T3")==0)|| (strcmp(ictr.coat_material,"t2")==0) )
{
strcpy(coating_material,"T3");
}
//T2-金属所
if ( (strcmp(ictr.coat_material,"T2")==0)|| (strcmp(ictr.coat_material,"t2")==0) )
{
strcpy(coating_material,"T2");
}
//T1-陈剑
if ( (strcmp(ictr.coat_material,"T1")==0)|| (strcmp(ictr.coat_material,"t1")==0) )
{
strcpy(coating_material,"T1");
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//ini coolant temperature and wall temperature
for (i=0;i<n;i++) {
T[i]=ictr.ini_coolant_temp;p[i]=ictr.ini_coolant_pressure;Z_conversion[i]=0;delta_t[i]=0;
   for (j=0;j<sinp[i].nummesh_1+sinp[i].nummesh_2+1;j++) {
   Tw[i][j]=ictr.ini_coolant_temp+150;
   }
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
outfile = fopen (".\\out.txt", "w") ; /* open the file */
if ( (strcmp(cooling_material,"n-Decane")==0) ){
fprintf(outfile,"位置m 对流换热系数 稳态热流     油温K     油压Pa      热面温度K  中间温度K 盖板温度K  流速   停留时间s   裂解率   总热沉MJ/kg\n"); 
}

if ( (strcmp(cooling_material,"IF97::water")==0) ){
fprintf(outfile,"位置m 对流换热系数 稳态热流     水温K     水压Pa      热面温度K  中间温度K 盖板温度K  流速   停留时间s    总热沉MJ/kg\n"); 
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// calulate along the flow direction 
for (i=0;i<n-1;i++) {
//step size of position along the flow
deltax              =ginp[i+1].position-ginp[i].position;
//average Mach number between i and i+1
ave_ma              =0.5*(ginp[i].ma+ginp[i+1].ma);
//average Total Temperature between i and i+1
ave_totaltemperature=0.5*(ginp[i].totaltemperature+ginp[i+1].totaltemperature);
//average Total Pressure between i and i+1
ave_totalpressure   =0.5*(ginp[i].totalpressure+ginp[i+1].totalpressure);
//hot gas static temperature
Ts=ave_totaltemperature/(1.0+(gam-1)/2.0*ave_ma*ave_ma);
//hot gas static pressure
ps=ave_totalpressure*pow((1.0+(gam-1)/2.0*ave_ma*ave_ma),-gam/(gam-1));
//hot gas velosity
Ug=ave_ma*sqrt(gam*Rg*Ts);
//hot gas mass flux
//mass_flow_gas=ave_ma*ps*sqrt(gam/(Rg*Ts));
area_old=sinp[i].channel_thickness*sinp[i].channel_width;
perimeter_old=2.0*(sinp[i].channel_thickness+sinp[i].channel_width);
equi_diameter_old=sqrt(4*area_old/pi);

area_new=sinp[i+1].channel_thickness*sinp[i+1].channel_width;
perimeter_new=2.0*(sinp[i+1].channel_thickness+sinp[i+1].channel_width);
equi_diameter_new=sqrt(4*area_new/pi);

if (equi_diameter_new>=equi_diameter_old) {
e_c_coeffient=pow(1-pow(equi_diameter_old/equi_diameter_new,2),2);
}
if (equi_diameter_new<equi_diameter_old) {
e_c_coeffient=0.5-0.167*equi_diameter_new/equi_diameter_old-0.125*pow(equi_diameter_new/equi_diameter_old,2)-0.208*pow(equi_diameter_new/equi_diameter_old,3);
}
rho_old=PropsSI("D","P",p[i],"T",T[i],cooling_material);
H_old  =PropsSI("H","P",p[i],"T",T[i],cooling_material);

//qw*(channel_width+rib_width)=q*2*(channel_width+channel_thickness)
coff=perimeter_new/(sinp[i+1].channel_width+sinp[i+1].rib_width);

delta_err=100;p[i+1]=p[i];T[i+1]=T[i];Tw[i+1][sinp[i+1].nummesh_1+sinp[i+1].nummesh_2]=Tw[i][sinp[i].nummesh_1+sinp[i].nummesh_2];
while (delta_err>1.0E-4)
{
//keep T[i+1] into proper range [200,1650K]
	if ( (T[i+1]<200)||(T[i+1]>1650) ) {
	T[i+1]=T[i];
	}
//keep Tw[i+1][sinp[i+1].nummesh_1+sinp[i+1].nummesh_2 [200K ,2500K]
   if ( (Tw[i+1][sinp[i+1].nummesh_1+sinp[i+1].nummesh_2]>2500)||(Tw[i+1][sinp[i+1].nummesh_1+sinp[i+1].nummesh_2]<200) ){
   Tw[i+1][sinp[i+1].nummesh_1+sinp[i+1].nummesh_2]=Tw[i][sinp[i].nummesh_1+sinp[i].nummesh_2];
}

rho_new     =PropsSI("D","P",p[i+1],"T",T[i+1],cooling_material);
rho_wall_new=PropsSI("D","P",p[i+1],"T",Tw[i+1][sinp[i+1].nummesh_1+sinp[i+1].nummesh_2],cooling_material);
mu_new      =PropsSI("V","P",p[i+1],"T",T[i+1],cooling_material);
mu_wall_new =PropsSI("V","P",p[i+1],"T",Tw[i+1][sinp[i+1].nummesh_1+sinp[i+1].nummesh_2],cooling_material);
k_new       =PropsSI("L","P",p[i+1],"T",T[i+1],cooling_material);
k_wall_new  =PropsSI("L","P",p[i+1],"T",Tw[i+1][sinp[i+1].nummesh_1+sinp[i+1].nummesh_2],cooling_material);
cp_new      =PropsSI("C","P",p[i+1],"T",T[i+1],cooling_material);
cp_wall_new =PropsSI("C","P",p[i+1],"T",Tw[i+1][sinp[i+1].nummesh_1+sinp[i+1].nummesh_2],cooling_material);
pr_new=mu_new*cp_new/k_new;
pr_wall_new=mu_wall_new*cp_wall_new/k_wall_new;
//momentum  pressure loss 1/3
momentum_pressure_loss=2.0/(area_old+area_new)*pow(ictr.mass_rate,2)*(1.0/rho_new/area_new-1.0/rho_old/area_old);
p[i+1]=p[i]-momentum_pressure_loss;

//friction pressure loss  2/3
Re_new=Re(ictr.mass_rate,equi_diameter_new,area_new,mu_new);
friction=f(Ra,equi_diameter_new,Re_new*rho_wall_new/rho_new*mu_new/mu_wall_new );
friction_pressure_loss=friction*pow(ictr.mass_rate,2)/2/rho_new/pow(area_new,2)/equi_diameter_new*deltax;
p[i+1]=p[i+1]-friction_pressure_loss;

//expansion or contraction pressure loss 3/3
ec_pressure_loss=e_c_coeffient*pow(ictr.mass_rate,2)/2/rho_new/pow(area_new,2);
p[i+1]=p[i+1]-ec_pressure_loss;

H_new  =PropsSI("H","P",p[i+1],"T",T[i+1],cooling_material);
h_conv=k_new/equi_diameter_new*Nusselt(p[i+1],T[i+1],Re_new,pr_new,pr_wall_new,mu_new,mu_wall_new,cooling_material);
//calculate rib effect
ks=CondMat(Tw[i+1][sinp[i+1].nummesh_1+sinp[i+1].nummesh_2],channel_material);
h_conv=h_conv*rib_effect(sinp[i+1].channel_width,sinp[i+1].rib_width,sinp[i+1].channel_thickness,h_conv, ks);

q=h_conv*(Tw[i+1][sinp[i+1].nummesh_1+sinp[i+1].nummesh_2]-T[i+1]);
//qw*(channel_width+rib_width)=q*2*(channel_width+channel_thickness)
//calculate qw
q=q*coff;


//find Tw[i+1][0] according to q
for (j=0;j<sinp[i+1].nummesh_2;j++)
{
ks=CondMat(Tw[i+1][sinp[i+1].nummesh_1+sinp[i+1].nummesh_2-j],channel_material);
deltay=sinp[i+1].hotwall_thickness/sinp[i+1].nummesh_2;
Tw[i+1][sinp[i+1].nummesh_1+sinp[i+1].nummesh_2-j-1]=q*deltay/ks+Tw[i+1][sinp[i+1].nummesh_1+sinp[i+1].nummesh_2-j];
}
for (j=0;j<sinp[i+1].nummesh_1;j++)
{
ks=CondMat(Tw[i+1][sinp[i+1].nummesh_1-j],coating_material);
deltay=sinp[i+1].coat_thickness/sinp[i+1].nummesh_1;
Tw[i+1][sinp[i+1].nummesh_1-j-1]=q*deltay/ks+Tw[i+1][sinp[i+1].nummesh_1-j];
}
//
if ( (Tw[i+1][0]<0)||(Tw[i+1][0]>ave_totaltemperature) ){
Tw[i+1][0]=Tw[i][0];
}
Taw=findTaw(Ts,ave_totaltemperature,Tw[i+1][0]);
hg0=hotgas(Taw, Ts,Tw[i+1][0],ps,Ug);

delta_t[i+1]=rho_new*area_new*deltax/ictr.mass_rate;

H_new_1       =PropsSI("H","P",p[i+1],"T",T[i+1]+delta,cooling_material);
H_new_2       =PropsSI("H","P",p[i+1]+100*delta,"T",T[i+1],cooling_material);
rho_new_1     =PropsSI("D","P",p[i+1],"T",T[i+1]+delta,cooling_material);
mu_new_1      =PropsSI("V","P",p[i+1],"T",T[i+1]+delta,cooling_material);
mu_wall_new_1 =PropsSI("V","P",p[i+1],"T",Tw[i+1][sinp[i+1].nummesh_1+sinp[i+1].nummesh_2]+delta,cooling_material);
k_new_1       =PropsSI("L","P",p[i+1],"T",T[i+1]+delta,cooling_material);
k_wall_new_1  =PropsSI("L","P",p[i+1],"T",Tw[i+1][sinp[i+1].nummesh_1+sinp[i+1].nummesh_2]+delta,cooling_material);
cp_new_1      =PropsSI("C","P",p[i+1],"T",T[i+1]+delta,cooling_material);
cp_wall_new_1 =PropsSI("C","P",p[i+1],"T",Tw[i+1][sinp[i+1].nummesh_1+sinp[i+1].nummesh_2]+delta,cooling_material);

dH_dT=cp_new;
dH_dp=(H_new_2-H_new)/(100*delta);
drho_dT=(rho_new_1-rho_new)/(delta);
dcp_dT=(cp_new_1-cp_new)/(delta);
dk_dT=(k_new_1-k_new)/(delta);
dmu_dT=(mu_new_1-mu_new)/(delta);
pr_wall_new_1=cp_wall_new_1*mu_wall_new_1/k_wall_new_1;
dmu_wall_dTw=(mu_wall_new_1-mu_wall_new)/(delta);
dpr_wall_dTw=(pr_wall_new_1-pr_wall_new)/(delta);

dp_dT=2.0*pow(ictr.mass_rate,2)/area_new/(area_new+area_old)+friction*pow(ictr.mass_rate,2)*deltax/equi_diameter_new/2/pow(area_new,2)+e_c_coeffient*pow(ictr.mass_rate,2)/2/pow(area_new,2);
dp_dT=dp_dT/pow(rho_new,2)*drho_dT;

dhl_dT=hl_T(p[i+1],T[i+1],ictr.mass_rate,equi_diameter_new,mu_wall_new,pr_wall_new,cp_new,dcp_dT,k_new,dk_dT,mu_new,dmu_dT,cooling_material);

dp_dTw=0;
dhl_dTw=hl_Tw(p[i+1],T[i+1],ictr.mass_rate,equi_diameter_new,cp_new,k_new,mu_new,mu_wall_new,pr_wall_new,dmu_wall_dTw,dpr_wall_dTw,cooling_material);

// if oil temperature large than 750K，add chemical reaction heat absorption flux qe_new
if ( (T[i+1]>750)&&(strcmp(cooling_material,"n-Decane")==0) ) {
k_A_old=k_Arrhenius(p[i],T[i]);
k_A_new=k_Arrhenius(p[i+1],T[i+1]);
k_A__new_1=k_Arrhenius(p[i+1],T[i+1]+delta);
k_A__new_2=k_Arrhenius(p[i+1]+100*delta,T[i+1]);

he_old   =CHS(p[i],T[i]);
he_new   =CHS(p[i+1],T[i+1]);
he_new_1 =CHS(p[i+1],T[i+1]+delta);
he_new_2 =CHS(p[i+1]+100*delta,T[i+1]);
Z_conversion[i+1] =1.0-(1.0-Z_conversion[i])*exp(-k_A_new   *delta_t[i+1]);
Z_conversion_new_1=1.0-(1.0-Z_conversion[i])*exp(-k_A__new_1*delta_t[i+1]);
Z_conversion_new_2=1.0-(1.0-Z_conversion[i])*exp(-k_A__new_2*delta_t[i+1]);

qe_new  =ictr.mass_rate/(sinp[i+1].channel_width+sinp[i+1].rib_width)/deltax*(Z_conversion[i+1] *he_new  -Z_conversion[i]*he_old);
qe_new_1=ictr.mass_rate/(sinp[i+1].channel_width+sinp[i+1].rib_width)/deltax*(Z_conversion_new_1*he_new_1-Z_conversion[i]*he_old);
qe_new_2=ictr.mass_rate/(sinp[i+1].channel_width+sinp[i+1].rib_width)/deltax*(Z_conversion_new_2*he_new_2-Z_conversion[i]*he_old);
dqe_dT=(qe_new_1-qe_new)/delta+(qe_new_2-qe_new)/(100*delta)*dp_dT;
}

F1=H_new+pow(ictr.mass_rate,2)/2/pow(area_new,2)/pow(rho_new,2)-H_old-pow(ictr.mass_rate,2)/2/pow(area_old,2)/pow(rho_old,2)-perimeter_new*deltax/ictr.mass_rate*h_conv*(Tw[i+1][sinp[i+1].nummesh_1+sinp[i+1].nummesh_2]-T[i+1]);
F2=hg0*(Taw-Tw[i+1][0])-q-qe_new;

dF1_dT =dH_dT+dH_dp*dp_dT-pow(ictr.mass_rate,2)/pow(area_old,2)/pow(rho_new,3)*drho_dT-perimeter_new*deltax/ictr.mass_rate*dhl_dT*(Tw[i+1][sinp[i+1].nummesh_1+sinp[i+1].nummesh_2]-T[i+1])+perimeter_new*deltax/ictr.mass_rate*h_conv;
dF1_dTw=dH_dp*dp_dTw-perimeter_new*deltax/ictr.mass_rate*dhl_dTw*(Tw[i+1][sinp[i+1].nummesh_1+sinp[i+1].nummesh_2]-T[i+1])-perimeter_new*deltax/ictr.mass_rate*h_conv;
dF2_dT =-coff*dhl_dT*(Tw[i+1][sinp[i+1].nummesh_1+sinp[i+1].nummesh_2]-T[i+1])+coff*h_conv-dqe_dT;
dF2_dTw=-hg0-coff*h_conv;

detF=dF1_dT*dF2_dTw-dF1_dTw*dF2_dT;

delta_T =(F2*dF1_dTw-F1*dF2_dTw)/detF;
delta_Tw=(F1*dF2_dT-F2*dF1_dT)/detF;

T[i+1]=T[i+1]+relax_number*delta_T;
Tw[i+1][sinp[i+1].nummesh_1+sinp[i+1].nummesh_2]=Tw[i+1][sinp[i+1].nummesh_1+sinp[i+1].nummesh_2]+relax_number*delta_Tw;

inter_k=inter_k+1;
if (inter_k>100){
relax_number=0.7*exp(-inter_k/100);
}
delta_err=sqrt(delta_T*delta_T+delta_Tw*delta_Tw); 
}
inter_k=0;
relax_number=0.7;
//output total heat absorption values
sum_H=sum_H+(H_new+Z_conversion[i+1]*he_new)-(H_old+Z_conversion[i]*he_old);

if ( (strcmp(cooling_material,"n-Decane")==0) ){
fprintf(outfile, "%.2lf    %.2lf    %.2lf    %.2lf    %.2lf    %.2lf    %.2lf    %.2lf    %.2lf    %.6lf    %.3lf    %.4lf \n",ginp[i+1].position,h_conv,q,T[i+1],p[i+1],Tw[i+1][0],Tw[i+1][sinp[i+1].nummesh_1],Tw[i+1][sinp[i+1].nummesh_1+sinp[i+1].nummesh_2],ictr.mass_rate/rho_new/area_new,delta_t[i+1],Z_conversion[i+1],sum_H/1.0E6);
}

if ( (strcmp(cooling_material,"IF97::water")==0) ){
fprintf(outfile, "%.2lf    %.2lf    %.2lf    %.2lf    %.2lf    %.2lf    %.2lf    %.2lf    %.2lf    %.6lf    %.4lf \n",ginp[i+1].position,h_conv,q,T[i+1],p[i+1],Tw[i+1][0],Tw[i+1][sinp[i+1].nummesh_1],Tw[i+1][sinp[i+1].nummesh_1+sinp[i+1].nummesh_2],ictr.mass_rate/rho_new/area_new,delta_t[i+1],sum_H/1.0E6);
}
printf("\n%f",(i+1.0)/n*100); 
printf("percent!");
}
printf("\nDone!"); 
printf("\nPress any key to continue!"); 
getchar();
return 0;
}



