#define _CRT_SECURE_NO_WARNINGS
#include "stdafx.h"
#include "disEarth.h"
#include "MathConstant.h"
#include "EarthConstant.h"


#include <math.h>
#include "UserMath.h"
#include "UserMechanics.h"

#include <stdio.h>
#include <string.h>

#define min(a, b)  (((a) < (b)) ? (a) : (b))
int	ReadFileLine(FILE *inFile, char chrLineBuffer[82]);

inline unsigned short CEarthAgent::Index(const unsigned char n, const unsigned char m)
{
	return n*(n+1)/2+m;
}
void CEarthAgent::Attraction(const unsigned char iN,const unsigned char iM, const double Xi[3],double Acc[3])
{
	unsigned char N,M;
	N = iN;	M = iM;
	if (N>dimCS-2)  N = dimCS-2;  if (M>dimCS-2)  M = dimCS-2;
	if (M>N)  M = N;

	unsigned char i;
	double X[3];
    for(i=0;i<3;i++)	X[i]=Xi[i];

	/*double Rei[3][3];
    GetRei(GreenwichAngle,Rei);
    RotateFrame(Rei,X);*/
	
	RotateFrame(this->Rei,X);

	double Radius;
	Radius = GetLength(X);
	if (Radius <= EarthRadius) Radius = EarthRadius;
	unsigned char n,m;
	for(n=0;n<=N+1;n++)
	{
		for(m=0;m<=N+1;m++)
		{
			Pnm[n][m]	= 0.0;
		}
	}
	double ZR = X[2]/Radius;
	Pnm[0][0] = 1;
	Pnm[1][0] = ZR;
	for(n=2;n<=N;n++)
	{
		Pnm[n][0] = (2.0*n-1.0)/n*ZR*Pnm[n-1][0]-(n-1.0)/n*Pnm[n-2][0];
	}
	for(m=1;m<=N+1;m++)
	{
		Pnm[m][m] = Pnm[m-1][m-1]*(2.0*m-1);
		for(n=m+1;n<=N;n++)
		{
			Pnm[n][m] = (2.0*n-1.0)/(n-m)*ZR*Pnm[n-1][m]-(n+m-1.0)/(n-m)*Pnm[n-2][m];
		}
	}
	U[0] = 1.0;
	V[0] = 0.0;
	for(m=1;m<=N;m++)
	{
		U[m] = (U[m-1]*X[0]-V[m-1]*X[1])/Radius;
		V[m] = (U[m-1]*X[1]+V[m-1]*X[0])/Radius;
	}
	for(i=0;i<3;i++) Acc[i] = 0;
	double Temp ;
	for(n=1;n<=N;n++)
	{
		Temp = EarthAttractionCoefficient*pow((EarthRadius/Radius),n);
		unsigned char minNM = n;
		if ( minNM > M ) minNM = M;
		for(m=0;m<=minNM;m++)
		{
			double CSUV = Temp/pow(Radius,3)*(Cnm[n][m]*U[m]+Snm[n][m]*V[m]);
			for (i=0;i<3;i++)
			{
				Acc[i] -= CSUV*((n+m+1.0)*Pnm[n][m]*X[i]+X[2]*X[i]/Radius*Pnm[n][m+1]);
			}
			Acc[2] += CSUV*Pnm[n][m+1]*Radius;
			if ( m != 0)
			{
				Acc[0] += Temp/Radius/Radius*Pnm[n][m]*m*( Cnm[n][m]*U[m-1]+Snm[n][m]*V[m-1]);
				Acc[1] += Temp/Radius/Radius*Pnm[n][m]*m*(-Cnm[n][m]*V[m-1]+Snm[n][m]*U[m-1]);
			}
		}
	}
    for(i=0;i<3;i++)	Acc[i] -= EarthAttractionCoefficient/Radius/Radius/Radius*X[i];
    InverseRotateFrame(this->Rei,Acc);

}
void CEarthAgent::SpinVelocity(const double X[3], double V[3])
{
	double W[3];
	for(unsigned char i=0;i<3;i++) W[i] = 0.0;	W[2] = EarthSpinSpeed;
	cross(W,X,V);
}
void CEarthAgent::GetAirVelocity(const double R[3], const double V[3],double Vair[3])
{
	this->SpinVelocity (R,Vair);
	unsigned char i;
	for(i=0;i<3;i++) Vair[i] = V[i] - Vair[i];
}
double CEarthAgent::AirDensity(const double HeightM)
{
	double Height = HeightM/1.0E3;

    if (Height >= 1000.0) return 0.0;
    short i = 47-1;
    while ( (Height < C[i][0]) && (i>0) ) i--;
    double DetHeight;
    DetHeight = Height-C[i][0];
    if (DetHeight < 0.0 ) DetHeight = 0.0;
	double temp;
	temp =  C[i][2]+C[i][3]*DetHeight+C[i][4]*DetHeight*DetHeight
                   +C[i][5]*DetHeight*DetHeight*DetHeight ;
    return pow(10.0,temp);
}
double CEarthAgent::Temperature(const double HeightM)
{
	if (HeightM <= 0 ) return T[0][1];
	double Height = HeightM/1.0E3/(1.0+HeightM/EarthRadius);
    short i = 8-1;
    if (Height >= T[i][0]) return T[i][1];
    while ( (Height < T[i][0]) && (i>0) ) i--;
	double temp;
    temp =  T[i][1]+(Height-T[i][0])*(T[i+1][1]-T[i][1])/(T[i+1][0]-T[i][0]);
	return temp;
}

double	CEarthAgent::Pressure(const double Height)
{
	if ( Height > P[17][0] ) return 0.0;
	unsigned char i = 0;
	while ( ( i < 17 )&&( Height > P[i+1][0]) ) i++ ;
	double temp;
	temp = P[i][1] + ( Height - P[i][0] )/( P[i+1][0] - P[i][0] )*( P[i+1][1] - P[i][1] );
	return temp;
}

double CEarthAgent::MachVelocity(const double Height)
{
	return 20.04*sqrt(fabs(Temperature(Height)));
}

void CEarthAgent::MagneticDensity(const unsigned char iN,const unsigned char iM, const double Xi[3],double Acc[3])
{
	unsigned char N,M;
	N = iN;	M = iM;
	if (N>dimGH-1)  N = dimGH-1;  if (M>dimGH-1)  M = dimGH-1;
	if (M>N)  M = N;

	unsigned char n,m;

	double Xe[3];
	for(unsigned char j=0;j<3;j++)	Xe[j]=Xi[j];
    RotateFrame(this->Rei,Xe);

	double Radius = GetLength(Xe);
	double longitude = atan2(Xe[1],Xe[0]);
	double latitude = asin(fLimitInOne(Xe[2]/Radius));
    double srlon = sin(longitude);
    double srlat = sin(latitude);
    double crlon = cos(longitude);

	if ( PI/2-fabs(latitude) < 1.0E-6 )
		return;

	pp[0] = 1.0;
	sp[0] = 0.0;
	cp[0] = 1.0;
    sp[1] = srlon;
    cp[1] = crlon;

/* CONVERT FROM GEODETIC COORDS. TO SPHERICAL COORDS.  */
	double ct,st;
	ct = srlat;
	st = sqrt(fabs(1.0-(ct*ct)));

	if ( fabs(latitude) != PI/2) 
	{
		for (m=2; m<=M; m++) 
		{
			sp[m] = sp[1]*cp[m-1]+cp[1]*sp[m-1];
			cp[m] = cp[1]*cp[m-1]-sp[1]*sp[m-1];
		}
	}

/*   COMPUTE UNNORMALIZED ASSOCIATED LEGENDRE POLYNOMIALS AND DERIVATIVES VIA RECURSION RELATIONS */

	double	aor = EarthRadius/Radius;
	double	ar = aor*aor;
	double	br,bt,bp,bpp;
	br = bt = bp = bpp = 0.0;
	MPnm[0][0] = 1.0;
	MDPnm[0][0] = 0.0;
	for (n=1; n<=N; n++) 
	{
		ar = ar*aor;
		for (m=0;m<=min(M,n);m++) 
		{
			if ( fabs(latitude) != PI/2) 
			{
				if (n == m) 
				{
					MPnm[n][m] = st*MPnm[n-1][m-1];
					MDPnm[n][m] = st*MDPnm[n-1][m-1]+ct*MPnm[n-1][m-1];
				}
				if (n == 1 && m == 0) 
				{
					MPnm[n][m] = ct*MPnm[n-1][m];
					MDPnm[n][m] = ct*MDPnm[n-1][m]-st*MPnm[n-1][m];
				}
				if (n > 1 && n != m) 
				{
					if (m > n-2) MPnm[n-2][m] = 0.0;
					if (m > n-2) MDPnm[n-2][m] = 0.0;
					MPnm[n][m] = ct*MPnm[n-1][m]-k[n][m]*MPnm[n-2][m];
					MDPnm[n][m] = ct*MDPnm[n-1][m] - st*MPnm[n-1][m]-k[n][m]*MDPnm[n-2][m];
				}
			}

/*		ACCUMULATE TERMS OF THE SPHERICAL HARMONIC EXPANSIONS */

			double	par,temp1,temp2;
			par = ar*MPnm[n][m];
			if (m == 0) 
			{
				temp1 = Gnm[n][m]*cp[m];
				temp2 = Gnm[n][m]*sp[m];
			}
			else 
			{
				temp1 = Gnm[n][m]*cp[m]+Hnm[n][m]*sp[m];
				temp2 = Gnm[n][m]*sp[m]-Hnm[n][m]*cp[m];
			}
			bt -= ar*temp1*MDPnm[n][m];
			bp += (m*temp2*par);
			br += ((n+1)*temp1*par);
/*		SPECIAL CASE:  NORTH/SOUTH GEOGRAPHIC POLES */
			if (st == 0.0 && m == 1) 
			{
				if (n == 1) pp[n] = pp[n-1];
				else pp[n] = ct*pp[n-1]-k[n][m]*pp[n-2];
				double parp = ar*pp[n];
				bpp += (m*temp2*parp);
			}
		}
	}
    if (st == 0.0) bp = bpp;
    else bp /= st;

/*    ROTATE MAGNETIC VECTOR COMPONENTS FROM SPHERICAL TO GEODETIC COORDINATES */

    Acc[0] = -bt;
    Acc[1] = bp;
    Acc[2] = -br;

/*    COMPUTE DECLINATION (DEC), INCLINATION (DIP) AND TOTAL INTENSITY (TI) */

	double Rvi[3][3];
	double Alongitude = atan2(Xi[1],Xi[0]);;
	double Blatitude = asin(fLimitInOne(Xi[2]/Radius));
    GetRvi(Alongitude,Blatitude, Rvi);
    InverseRotateFrame(Rvi,Acc);

}
void CEarthAgent::InertiaToFixed(const double Xinertia[3],double xP[3])
{
	for(unsigned char i=0;i<3;i++)	xP[i] = Xinertia[i];
	RotateFrame(this->Rei,xP);
}
void CEarthAgent::FixedToInertia(const double xP[3],double Xinertia[3])
{
	for(unsigned char i=0;i<3;i++)	Xinertia[i] = xP[i];
	InverseRotateFrame(this->Rei,Xinertia);
}

void CEarthAgent::TrackToPosition(const struct SGeographyElements &HLB,double xP[3])
{
	double Xfixed[3];
	this->HLBToEarthFixedPosition(HLB,Xfixed);
	this->FixedToInertia(Xfixed,xP);
}
void CEarthAgent::HLBToEarthFixedPosition(const struct SGeographyElements &HLB,double xP[3])
{
	double Longitude = HLB.Longitude;
	double Latitude = HLB.Latitude ;
	double Altitude = HLB.Altitude ;

	double Eccentricity2 = 2*EarthCuvature-EarthCuvature*EarthCuvature;
	double RadiusN = EarthSemimajorAxis/sqrt(fabs(1-Eccentricity2*pow(sin(Latitude),2)));
	
	xP[0] = (RadiusN+Altitude)*cos(Latitude)*cos(Longitude);
	xP[1] = (RadiusN+Altitude)*cos(Latitude)*sin(Longitude);
	xP[2] = (RadiusN*(1-Eccentricity2)+Altitude)*sin(Latitude);
}
void CEarthAgent::EarthFixedPositionToHLB (const double X[3],struct SGeographyElements &HLB)
{

	double        Radius,Altitude,Elevation,longitude,latitude;

	Radius = GetLength(X);
	Elevation = asin(fLimitInOne(X[2]/Radius));
	longitude = atan2(X[1],X[0]);
	longitude = fLimitInPI(longitude);

	double RadiusN = EarthSemimajorAxis;
	Altitude = Radius - sqrt(EarthSemimajorAxis*EarthSemimajorAxis*(1.0-EarthCuvature));
	double Eccentricity2 = 2*EarthCuvature-EarthCuvature*EarthCuvature;
	if( (X[0]*X[0]+X[1]*X[1]) <= 1.0E-10 )
	{
		if (X[2]>0)
			latitude=PI/2;
		else
		    latitude=-PI/2;
        RadiusN = EarthSemimajorAxis*(1.0-EarthCuvature);
        Altitude = Radius-RadiusN;
	}
	else
	{	
		latitude = 0.0;
		for(int i=0;i<4;i++)
		{
			latitude = atan(X[2]/sqrt(X[0]*X[0]+X[1]*X[1])/(1-Eccentricity2*RadiusN/(RadiusN+Altitude)));
			RadiusN = EarthSemimajorAxis/sqrt(fabs(1-Eccentricity2*pow(sin(latitude),2)));
			Altitude = sqrt(X[0]*X[0]+X[1]*X[1])/cos(latitude)-RadiusN;
		}
	}
	HLB.Altitude = Altitude;
	HLB.Longitude = longitude;
	HLB.Latitude = latitude;
}
void CEarthAgent::Track(const double X[3],struct SGeographyElements &HLB)
{
	double Xfixed[3];
	this->InertiaToFixed(X,Xfixed);
	this->EarthFixedPositionToHLB(Xfixed,HLB);
	HLB.Longitude = fLimitInPI(HLB.Longitude);
}
void CEarthAgent::OrbitElementToRectangularElement (const struct SOrbitalElements &xa,double R[3],double V[3])
{

	double a,e,i,Q,w,f;
	a	 = xa.SemimajorAxis ;
	e	 = xa.Eccentricity ;
	i	 = xa.Inclination ;
	Q	 = xa.AscendingNode ;
	w	 = xa.Argument ;
	f	 = xa.TrueAnormaly ;

	double sini,cosi,sinq,cosq,sinw,cosw,sinf,cosf;
	sini = sin(i);   	cosi = cos(i);
	sinq = sin(Q);   	cosq = cos(Q);
	sinw = sin(w);   	cosw = cos(w);
	sinf = sin(f);   	cosf = cos(f);

	double u1,u2,v1,v2,p1,p2;
	u1   =  cosq*cosw-sinq*cosi*sinw;
	u2   = -cosq*sinw-sinq*cosi*cosw;
	v1   =  sinq*cosw+cosq*cosi*sinw;
	v2   = -sinq*sinw+cosq*cosi*cosw;
	p1   =  sini*sinw;
	p2   =  sini*cosw;

	double tmp,x0,y0;
	tmp  =  a*(1.0-e*e)/(1.0+e*cosf);
	x0   =  tmp*cosf;
	y0   =  tmp*sinf;

	R[0] = u1*x0+u2*y0;
	R[1] = v1*x0+v2*y0;
	R[2] = p1*x0+p2*y0;

	double un = EarthAttractionCoefficient;
	tmp  = sqrt(fabs(un/a/(1.0-e*e)));
	x0	 = tmp*sinf;
	y0	 = tmp*(e+cosf);

	V[0] = -u1*x0+u2*y0;
	V[1] = -v1*x0+v2*y0;
	V[2] = -p1*x0+p2*y0;

}
void CEarthAgent::RectangularElementToOrbitElement(const double R[3],const double V[3], struct SOrbitalElements &xa)
{
	double   un = EarthAttractionCoefficient;
	double r[3],v[3],h[3],n[3],e[3];
	unsigned char i;
	for(i=0;i<3;i++)
	{
		r[i] = R[i];
		v[i] = V[i];
		e[i] = 0;
	}
	cross(r,v,h);
	e[2] =1;
	cross(e,h,n);

	double rr = GetLength(r);
	double v2 = dot(v,v);

	double tmp1 = v2 - un/rr;
	double tmp2 = dot(r,v);
	for(i=0;i<3;i++) e[i] = (tmp1*r[i]-tmp2*v[i])/un;

	xa.SemimajorAxis = un*rr/(2.0*un-rr*v2);
	xa.Eccentricity = GetLength(e);
	xa.Inclination = acos(fLimitInOne(h[2]/GetLength(h)));

	if( xa.Inclination <= 1.0E-10 )
	{
		xa.AscendingNode = 0;
		n[0]  = 1.0;	n[1] = 0.0;	n[2] = 0.0;
	}
	else
	{
		xa.AscendingNode = atan2(n[1],n[0]);
	}
	if ( xa.Eccentricity > 1.0E-10 )
	{
		tmp1 = dot(n,e)/GetLength(n)/xa.Eccentricity;
		xa.Argument = acos(fLimitInOne(tmp1));
		if(e[2]<0) xa.Argument = -xa.Argument;
	}
	else
	{
		xa.Argument = 0.0;
	}
	tmp1 = dot(n,r)/rr/GetLength(n);
	xa.TrueAnormaly = acos(fLimitInOne(tmp1));
	if(r[2]<0) xa.TrueAnormaly = -xa.TrueAnormaly;
	xa.TrueAnormaly -= xa.Argument;

}
void CEarthAgent::RectangularElementToOrbitVector(const double R[3],const double V[3],struct SOrbitalVectors &xa)
{
	double   un = EarthAttractionCoefficient;
	double r[3],v[3],h[3],e[3];
	unsigned char i;
	for(i=0;i<3;i++)
	{
		r[i] = R[i];
		v[i] = V[i];
		e[i] = 0;
	}
	cross(r,v,h);

	double rr = GetLength(r);
	double v2 = dot(v,v);

	double tmp1 = v2 - un/rr;
	double tmp2 = dot(r,v);
	for(i=0;i<3;i++) e[i] = (tmp1*r[i]-tmp2*v[i])/un;

	xa.SemimajorAxis = un*rr/(2.0*un-rr*v2);
	double ee = GetLength(e);
	double w = atan2(e[1],e[0]);
	xa.EccentricityX = ee*cos(w);
	xa.EccentricityY = ee*sin(w);

	double ii = acos(fLimitInOne(h[2]/GetLength(h)));
	double Q = atan2(h[1],h[0])+0.5*PI;
	xa.InclinationX = ii*cos(Q);
	xa.InclinationY = ii*sin(Q);

	xa.TrueAnormaly = atan2(r[1],r[0]);

}
void CEarthAgent::GetFixedToLocalMatrix(const struct SGeographyElements &HLB,double REe[3][3])
{
    double sinlong,coslong,sinlat,coslat;
    sinlong = sin(HLB.Longitude);
    coslong = cos(HLB.Longitude);
    sinlat  = sin(HLB.Latitude);
    coslat  = cos(HLB.Latitude);
    REe[0][0] = -coslong*sinlat;
    REe[0][1] = -sinlong*sinlat;
    REe[0][2] =  coslat;
    REe[1][0] = -sinlong;
    REe[1][1] =  coslong;
    REe[1][2] =  0;
    REe[2][0] = -coslat*coslong;
    REe[2][1] = -coslat*sinlong;
    REe[2][2] = -sinlat;
}

void CEarthAgent::GetInertiaToFixedMatrix(double Rei[3][3])
{
	for(unsigned char i=0;i<3;i++)
		for(unsigned char j=0;j<3;j++)
			Rei[i][j] =  this->Rei[i][j];
}

void CEarthAgent::GetInertiaToLocalMatrix(const struct SGeographyElements &HLB,double REi[3][3])
{
	double Rei[3][3],REe[3][3];
	this->GetInertiaToFixedMatrix(Rei);
	this->GetFixedToLocalMatrix(HLB,REe);
	MatrixProduct(3,3,3,&REe[0][0],&Rei[0][0],&REi[0][0]);
}
void CEarthAgent::QuaternionToLocalEularAngle(const struct SGeographyElements &HLB, const double Q[4], struct SEularAngle &Eular)
{
	double RbE[3][3],REi[3][3],Rbi[3][3];
	GetRbi(Q,Rbi); 
	this->GetInertiaToLocalMatrix(HLB,REi);
	MatrixProductTranspose(3,3,3,&Rbi[0][0],&REi[0][0],&RbE[0][0]);
	EularXYZExtract(RbE,&Eular);
}
void CEarthAgent::LocalEularAngleToQuaternion(const struct SGeographyElements &HLB, const struct SEularAngle &Eular, double Q[4])
{
	double RbE[3][3],REi[3][3],Rbi[3][3];
	EularXYZToMatrix(&Eular,RbE);
	this->GetInertiaToLocalMatrix(HLB,REi);
	MatrixProduct(3,3,3,&RbE[0][0],&REi[0][0],&Rbi[0][0]);
	QuaternionExtract(Rbi,Q );
}
void CEarthAgent::InertiaToGPS(const double R[3], const double V[3],struct SGeographyElements &HLB, struct SSphereElements &VAE)
{
	this->Track(R,HLB);
	double Vair[3];
	this->GetAirVelocity(R,V,Vair);
	double REi[3][3];
	this->GetInertiaToLocalMatrix(HLB,REi);
	RotateFrame(REi,Vair);
	RectAngularToSphere(Vair,&VAE);
}
void CEarthAgent::GPSToInertia(const struct SGeographyElements &HLB, const struct SSphereElements &VAE,double R[3],double V[3])
{
	this->TrackToPosition(HLB,R);
	double Vair[3];
	SphereToRectAngular(&VAE,Vair);
	double REi[3][3];
	this->GetInertiaToLocalMatrix(HLB,REi);
	InverseRotateFrame(REi,Vair);
	this->SpinVelocity(R,V);
	for(unsigned char i=0;i<3;i++) V[i] += Vair[i];
}

CEarthAgent::CEarthAgent()
{

	g10	= -30.339E-6;
	g11	=  -2.123E-6;
	h11	=   5.738E-6;

	C[0][0]=0.				;C[0][1]=.122500E+01	;C[0][2]=.881361E-01;
	C[0][3]=-.419340E-01	;C[0][4]=-.235728E-03	;C[0][5]=-.585315E-04;
	C[1][0]=2.				;C[1][1]=.100660E+01	;C[1][2]=.285693E-02;
	C[1][3]=-.435793E-01	;C[1][4]=-.586917E-03	;C[1][5]=.148469E-04;
	C[2][0]=4.				;C[2][1]=.819350E+00	;C[2][2]=-.865305E-01;
	C[2][3]=-.457488E-01	;C[2][4]=-.497835E-03	;C[2][5]=-.455277E-04;
	C[3][0]=6.				;C[3][1]=.660110E+00	;C[3][2]=-.180384E+00;
	C[3][3]=-.482865E-01	;C[3][4]=-.771001E-03	;C[3][5]=.106618E-03;
	C[4][0]=8.				;C[4][1]=.525790E+00	;C[4][2]=-.279188E+00;
	C[4][3]=-.500911E-01	;C[4][4]=-.131292E-03	;C[4][5]=-.452379E-03;
	C[5][0]=10.				;C[5][1]=.413510E+00	;C[5][2]=-.383514E+00;
	C[5][3]=-.560448E-01	;C[5][4]=-.284557E-02	;C[5][5]=.132107E-03;
	C[6][0]=12.				;C[6][1]=.311940E+00	;C[6][2]=-.505929E+00;
	C[6][3]=-.658417E-01	;C[6][4]=-.205292E-02	;C[6][5]=.436527E-03;
	C[7][0]=14.				;C[7][1]=.227860E+00	;C[7][2]=-.642332E+00;
	C[7][3]=-.688151E-01	;C[7][4]=.566240E-03	;C[7][5]=-.120855E-03;
	C[8][0]=16.				;C[8][1]=.166470E+00	;C[8][2]=-.778664E+00;
	C[8][3]=-.680004E-01	;C[8][4]=-.158892E-03   ;C[8][5]=.515678E-04;
	C[9][0]=18.				;C[9][1]=.121650E+00	;C[9][2]=-.914888E+00;
	C[9][3]=-.680172E-01	;C[9][4]=.150514E-03	;C[9][5]=-.911512E-04;
	C[10][0]=20.			;C[10][1]=.889100E-01	;C[10][2]=-.105105E+01;
	C[10][3]=-.685089E-01	;C[10][4]=-.396393E-03  ;C[10][5]=.517994E-04;
	C[11][0]=25.			;C[11][1]=.400840E-01	;C[11][2]=-.139703E+01;
	C[11][3]=-.685879E-01	;C[11][4]=.380598E-03	;C[11][5]=-.359416E-04;
	C[12][0]=30.			;C[12][1]=.184100E-01	;C[12][2]=-.173495E+01;
	C[12][3]=-.674775E-01	;C[12][4]=-.158526E-03	;C[12][5]=.307354E-04;
	C[13][0]=35.			;C[13][1]=.846340E-02	;C[13][2]=-.207246E+01;
	C[13][3]=-.667576E-01	;C[13][4]=.302506E-03	;C[13][5]=.218860E-05;
	C[14][0]=40.			;C[14][1]=.399570E-02	;C[14][2]=-.239841E+01;
	C[14][3]=-.635684E-01	;C[14][4]=.335335E-03	;C[14][5]=.121259E-04;
	C[15][0]=45.			;C[15][1]=.196630E-02	;C[15][2]=-.270635E+01;
	C[15][3]=-.593056E-01	;C[15][4]=.517224E-03	;C[15][5]=.118082E-04;
	C[16][0]=50.			;C[16][1]=.102690E-02	;C[16][2]=-.298847E+01;
	C[16][3]=-.532478E-01	;C[16][4]=.694347E-03	;C[16][5]=-.657844E-04;
	C[17][0]=55.			;C[17][1]=.568100E-03	;C[17][2]=-.324558E+01;
	C[17][3]=-.512382E-01	;C[17][4]=-.292419E-03	;C[17][5]=-.829393E-07;
	C[18][0]=60.			;C[18][1]=.309680E-03	;C[18][2]=-.350909E+01;  
	C[18][3]=-.541686E-01	;C[18][4]=-.293663E-03  ;C[18][5]=.144155E-06;
	C[19][0]=65.			;C[19][1]=.163210E-03	;C[19][2]=-.378725E+01;
	C[19][3]=-.570944E-01	;C[19][4]=-.291501E-03	;C[19][5]=-.144392E-04;
	C[20][0]=70.			;C[20][1]=.828290E-04	;C[20][2]=-.408182E+01;
	C[20][3]=-.610923E-01	;C[20][4]=-.508089E-03  ;C[20][5]=.946301E-05;
	C[21][0]=75.			;C[21][1]=.399210E-04	;C[21][2]=-.439880E+01;
	C[21][3]=-.654635E-01	;C[21][4]=-.366144E-03  ;C[21][5]=.116343E-04;
	C[22][0]=80.			;C[22][1]=.184580E-04	;C[22][2]=-.473382E+01;
	C[22][3]=-.682524E-01	;C[22][4]=-.191629E-03	;C[22][5]=-.422518E-04;
	C[23][0]=85.			;C[23][1]=.821960E-05	;C[23][2]=-.508515E+01;  
	C[23][3]=-.733375E-01	;C[23][4]=-.825405E-03  ;C[23][5]=.479193E-04;
	C[24][0]=90.			;C[24][1]=.341600E-05	;C[24][2]=-.546648E+01;
	C[24][3]=-.779976E-01	;C[24][4]=-.106615E-03  ;C[24][5]=.561803E-05;
	C[25][0]=100.			;C[25][1]=.560400E-06	;C[25][2]=-.625150E+01;
	C[25][3]=-.784445E-01	;C[25][4]=.619256E-04	;C[25][5]=.168843E-04;
	C[26][0]=110.			;C[26][1]=.970800E-07	;C[26][2]=-.701287E+01;
	C[26][3]=-.721407E-01	;C[26][4]=.568454E-03	;C[26][5]=.241761E-04;
	C[27][0]=120.			;C[27][1]=.222200E-07	;C[27][2]=-.765326E+01;  
	C[27][3]=-.535188E-01	;C[27][4]=.129374E-02	;C[27][5]=-.296655E-04;
	C[28][0]=130.			;C[28][1]=.815200E-08	;C[28][2]=-.808874E+01;  
	C[28][3]=-.365437E-01	;C[28][4]=.403772E-03	;C[28][5]=-.289208E-05;
	C[29][0]=140.			;C[29][1]=.383100E-08	;C[29][2]=-.841669E+01;
	C[29][3]=-.293359E-01	;C[29][4]=.317010E-03	;C[29][5]=-.442685E-05;
	C[30][0]=150.			;C[30][1]=.207600E-08	;C[30][2]=-.868277E+01;
	C[30][3]=-.243238E-01	;C[30][4]=.184205E-03	;C[30][5]=-.144722E-05;
	C[31][0]=160.			;C[31][1]=.123300E-08	;C[31][2]=-.890904E+01;
	C[31][3]=-.210738E-01	;C[31][4]=.140788E-03	;C[31][5]=-.137461E-05;
	C[32][0]=170.			;C[32][1]=.781500E-09	;C[32][2]=-.910707E+01;
	C[32][3]=-.186705E-01	;C[32][4]=.995495E-04	;C[32][5]=-.677454E-06;
	C[33][0]=180.			;C[33][1]=.519400E-09	;C[33][2]=-.928450E+01;
	C[33][3]=-.168827E-01	;C[33][4]=.792259E-04	;C[33][5]=-.593218E-06;
	C[34][0]=190.			;C[34][1]=.358100E-09	;C[34][2]=-.944600E+01;
	C[34][3]=-.154761E-01	;C[34][4]=.614293E-04	;C[34][5]=-.381119E-06;
	C[35][0]=200.			;C[35][1]=.254100E-09	;C[35][2]=-.959500E+01;  
	C[35][3]=-.143619E-01	;C[35][4]=.499958E-04	;C[35][5]=-.249568E-06;
	C[36][0]=220.			;C[36][1]=.136700E-09	;C[36][2]=-.986423E+01;
	C[36][3]=-.126615E-01	;C[36][4]=.350217E-04	;C[36][5]=-.154281E-06;
	C[37][0]=240.			;C[37][1]=.785800E-10	;C[37][2]=-.101047E+02;  
	C[37][3]=-.114458E-01	;C[37][4]=.257648E-04	;C[37][5]=-.925137E-07;
	C[38][0]=260.			;C[38][1]=.474200E-10	;C[38][2]=-.103240E+02;  
	C[38][3]=-.105262E-01	;C[38][4]=.202140E-04	;C[38][5]=-.774691E-07;
	C[39][0]=280.			;C[39][1]=.297100E-10	;C[39][2]=-.105271E+02;  
	C[39][3]=-.981064E-02	;C[39][4]=.155659E-04	;C[39][5]=-.650883E-07;
	C[40][0]=300.			;C[40][1]=.191600E-10	;C[40][2]=-.107176E+02;  
	C[40][3]=-.926611E-02	;C[40][4]=.116606E-04	;C[40][5]=-.249222E-07;
	C[41][0]=400.			;C[41][1]=.280200E-11	;C[41][2]=-.115525E+02;
	C[41][3]=-.768166E-02	;C[41][4]=.418393E-05	;C[41][5]=-.388723E-08;
	C[42][0]=500.			;C[42][1]=.521500E-12	;C[42][2]=-.122827E+02;  
	C[42][3]=-.696149E-02	;C[42][4]=.301776E-05	;C[42][5]=.447753E-08;
	C[43][0]=600.			;C[43][1]=.113700E-12	;C[43][2]=-.129442E+02;  
	C[43][3]=-.622361E-02	;C[43][4]=.436102E-05	;C[43][5]=.998741E-08;
	C[44][0]=700.			;C[44][1]=.306900E-13	;C[44][2]=-.135130E+02;  
	C[44][3]=-.505179E-02	;C[44][4]=.735724E-05	;C[44][5]=-.124098E-10;
	C[45][0]=800.			;C[45][1]=.113600E-13	;C[45][2]=-.139446E+02;
	C[45][3]=-.358071E-02	;C[45][4]=.735352E-05	;C[45][5]=-.104955E-07;
	C[46][0]=900.			;C[46][1]=.575900E-14	;C[46][2]=-.142397E+02;
	C[46][3]=-.242487E-02	;C[46][4]=.420487E-05	;C[46][5]=-.833682E-08;


	T[0][0] = 0				;T[0][1] = 288.15;
	T[1][0] = 11.0			;T[1][1] = 216.65;
	T[2][0] = 20.0			;T[2][1] = 216.65;
	T[3][0] = 32.0			;T[3][1] = 228.65;
	T[4][0] = 47.0			;T[4][1] = 270.65;
	T[5][0] = 51.0			;T[5][1] = 270.65;
	T[6][0] = 71.0			;T[6][1] = 214.65;
	T[7][0] = 86.0			;T[7][1] = 186.95;

	P[0][0] = 0E+03             ;P[0][1] = .101325E+06;
	P[1][0] = 1E+03				;P[1][1] = .898763E+05;
	P[2][0] = 2E+03				;P[2][1] = .795014E+05;
	P[3][0] = 3E+03				;P[3][1] = .701211E+05;
	P[4][0] = 4E+03				;P[4][1] = .616604E+05;
	P[5][0] = 5E+03				;P[5][1] = .540483E+05;
	P[6][0] = 6E+03				;P[6][1] = .472176E+05;
	P[7][0] = 7E+03				;P[7][1] = .411053E+05;
	P[8][0] = 8E+03				;P[8][1] = .356516E+05;
	P[9][0] = 9E+03				;P[9][1] = .308007E+05;
	P[10][0] = 10E+03			;P[10][1] = .264999E+05;
	P[11][0] = 20E+03			;P[11][1] = .552932E+04;
	P[12][0] = 30E+03			;P[12][1] = .119703E+04;
	P[13][0] = 40E+03			;P[13][1] = .287144E+03;
	P[14][0] = 50E+03			;P[14][1] = .797788E+02;
	P[15][0] = 60E+03			;P[15][1] = .219586E+02;
	P[16][0] = 70E+03			;P[16][1] = .522088E+01;
	P[17][0] = 80E+03			;P[17][1] = .105247E+01;


	Cnm[0][0] =	0.0;
	Cnm[1][0] = 0.0;
	Cnm[1][1] = 0.0;
	Snm[0][0] =	0.0;
	Snm[1][0] = 0.0;
	Snm[1][1] = 0.0;

	unsigned char n,m;

	for(n=0;n<dimCS;n++)
	{
		for(m=0;m<dimCS;m++)
		{
			Cnm[n][m] = 0.0;
			Snm[n][m] = 0.0;
		}
	};

	FILE	*fp;
	fp = fopen("JGM3.grv","r");
	char	chrLineBuffer[82];
	unsigned short		usLine = 0;
	unsigned short		nn,mm;
	double	tmp1,tmp2;
	do
	{
		ReadFileLine(fp,chrLineBuffer);	usLine++;
		if(usLine > 21)
		{	
			sscanf_s(chrLineBuffer,"%hd %hd %lf %lf",&nn,&mm,&tmp1,&tmp2);
			Cnm[nn][mm] = tmp1;
			Snm[nn][mm] = tmp2;

		}

	}while(usLine < 2575);
	fclose(fp);
	double	det,NormalizeCoefficient;
	for(n=0;n<dimCS;n++)
	{
		for(m=0;m<n+1;m++)
		{
			if (m == 0)
			{
				det = 0;
			} 
			else
			{
				det = 1;
			}
			NormalizeCoefficient = sqrt(FactorialDivision(n+m,n-m)/(1+det)/(2*n+1));
			Cnm[n][m] = Cnm[n][m]/NormalizeCoefficient;
			Snm[n][m] = Snm[n][m]/NormalizeCoefficient;
		}		
	};

	for(n=0;n<dimGH;n++)
	{
		for(m=0;m<dimGH;m++)
		{
			Gnm[n][m] = 0.0;
			Hnm[n][m] = 0.0;
		}
	};

	Gnm[1][0] = -29682E-9;
	Gnm[1][0] = -30.339E-6;
	Gnm[1][1] = -1789E-9;
	Gnm[1][1] =  -2.123E-6;
	Gnm[2][0] = -2197E-9;
	Gnm[2][1] = 3074E-9;
	Gnm[2][2] = 1685E-9;

	Gnm[3][0] = 1329E-9;
	Gnm[3][1] = -2268E-9;
	Gnm[3][2] = 1249E-9;
	Gnm[3][3] = 769E-9;
	Gnm[4][0] = 941E-9;

	Gnm[4][1] = 782E-9;
	Gnm[4][2] = 291E-9;
	Gnm[4][3] = -421E-9;
	Gnm[4][4] = 116E-9;
	Gnm[5][0] = -210E-9;

	Gnm[5][1] = 352E-9;
	Gnm[5][2] = 237E-9;
	Gnm[5][3] = -122E-9;
	Gnm[5][4] = -167E-9;
	Gnm[5][5] = -26E-9;

	Gnm[6][0] = 66E-9;
	Gnm[6][1] = 64E-9;
	Gnm[6][2] = 65E-9;
	Gnm[6][3] = -172E-9;
	Gnm[6][4] = 2E-9;

	Gnm[6][5] = 17E-9;
	Gnm[6][6] = -94E-9;
	Gnm[7][0] = 78E-9;
	Gnm[7][1] = -67E-9;
	Gnm[7][2] = 1E-9;

	Gnm[7][3] = 29E-9;
	Gnm[7][4] = 4E-9;
	Gnm[7][5] = 8E-9;
	Gnm[7][6] = 9E-9;
	Gnm[7][7] = -2E-9;

	Gnm[8][0] = 24E-9;
	Gnm[8][1] = 4E-9;
	Gnm[8][2] = -1E-9;
	Gnm[8][3] = -9E-9;
	Gnm[8][4] = -14E-9;

	Gnm[8][5] = 4E-9;
	Gnm[8][6] = 5E-9;
	Gnm[8][7] = 0E-9;
	Gnm[8][8] = -7E-9;
	Gnm[9][0] = 4E-9;

	Gnm[9][1] = 9E-9;
	Gnm[9][2] = 1E-9;
	Gnm[9][3] = -12E-9;
	Gnm[9][4] = 9E-9;
	Gnm[9][5] = -4E-9;
	

	Gnm[9][6] = -2E-9;
	Gnm[9][7] = 7E-9;
	Gnm[9][8] = 0E-9;
	Gnm[9][9] = -6E-9;
	Gnm[10][0] = -3E-9;

	Gnm[10][1] = -4E-9;
	Gnm[10][2] = 2E-9;
	Gnm[10][3] = -5E-9;
	Gnm[10][4] = -2E-9;
	Gnm[10][5] = 4E-9;

	Gnm[10][6] = 3E-9;
	Gnm[10][7] = 1E-9;
	Gnm[10][8] = 3E-9;
	Gnm[10][9] = 3E-9;
	Gnm[10][10] = 0E-9;

	Hnm[1][0] = 0E-9;
	Hnm[1][1] = 5318E-9;
	Hnm[1][1]=   5.738E-6;
	Hnm[2][0] = 0E-9;
	Hnm[2][1] = -2356E-9;
	Hnm[2][2] = -425E-9;

	Hnm[3][0] = 0E-9;
	Hnm[3][1] = -263E-9;
	Hnm[3][2] = 302E-9;
	Hnm[3][3] = -406E-9;
	Hnm[4][0] = 0E-9;

	Hnm[4][1] = 262E-9;
	Hnm[4][2] = -232E-9;
	Hnm[4][3] = 98E-9;
	Hnm[4][4] = -301E-9;
	Hnm[5][0] = 0E-9;

	Hnm[5][1] = 44E-9;
	Hnm[5][2] = 157E-9;
	Hnm[5][3] = -152E-9;
	Hnm[5][4] = -64E-9;
	Hnm[5][5] = 99E-9;

	Hnm[6][0] = 0E-9;
	Hnm[6][1] = -16E-9;
	Hnm[6][2] = 77E-9;
	Hnm[6][3] = 67E-9;
	Hnm[6][4] = -57E-9;

	Hnm[6][5] = 4E-9;
	Hnm[6][6] = 28E-9;
	Hnm[7][0] = 0E-9;
	Hnm[7][1] = -77E-9;
	Hnm[7][2] = -25E-9;

	Hnm[7][3] = 3E-9;
	Hnm[7][4] = 22E-9;
	Hnm[7][5] = 16E-9;
	Hnm[7][6] = -23E-9;
	Hnm[7][7] = -3E-9;

	Hnm[8][0] = 0E-9;
	Hnm[8][1] = 12E-9;
	Hnm[8][2] = -20E-9;
	Hnm[8][3] = 7E-9;
	Hnm[8][4] = -21E-9;

	Hnm[8][5] = 12E-9;
	Hnm[8][6] = 10E-9;
	Hnm[8][7] = -11E-9;
	Hnm[8][8] = -10E-9;
	Hnm[9][0] = 0E-9;

	Hnm[9][1] = -19E-9;
	Hnm[9][2] = 15E-9;
	Hnm[9][3] = 11E-9;
	Hnm[9][4] = -7E-9;
	Hnm[9][5] = -7E-9;
	
	Hnm[9][6] = 9E-9;
	Hnm[9][7] = 7E-9;
	Hnm[9][8] = -8E-9;
	Hnm[9][9] = 1E-9;
	Hnm[10][0] = 0E-9;

	Hnm[10][1] = 2E-9;
	Hnm[10][2] = 1E-9;
	Hnm[10][3] = 3E-9;
	Hnm[10][4] = 6E-9;
	Hnm[10][5] = -4E-9;

	Hnm[10][6] = 0E-9;
	Hnm[10][7] = -2E-9;
	Hnm[10][8] = 3E-9;
	Hnm[10][9] = -1E-9;
	Hnm[10][10] = -6E-9;

	double Snorm[dimGH][dimGH];
    Snorm[0][0] = 1.0;

    for (n=1; n<=dimGH-1; n++) 
    {
		Snorm[n][0] = Snorm[n-1][0]*(double)(2*n-1)/(double)n;
		for (m=0;m<=n;m++) 
		{
			k[n][m] = (double)(((n-1)*(n-1))-(m*m))/(double)((2*n-1)*(2*n-3));
			if (m > 0) 
			{
				double flnmj;
				if ( m== 1 )	flnmj = (double)((n-m+1)*2)/(double)(n+m);
				else			flnmj = (double)((n-m+1)*1)/(double)(n+m);
				Snorm[n][m] = Snorm[n][m-1]*sqrt(fabs(flnmj));
			}
			Gnm[n][m] = Snorm[n][m]*Gnm[n][m];
			Hnm[n][m] = Snorm[n][m]*Hnm[n][m];
		}
    }
    k[1][1] = 0.0;

	return ;

}

inline double CEarthAgent::FactorialDivision(unsigned short a,unsigned short b)
{
	if(a < b)		return 0;
	if(a == b)		return 1;
	unsigned short		i;
	double				result = 1;
	for (i=a;i>b;i--)
	{
		result = result*i;
	}
	return result;

}