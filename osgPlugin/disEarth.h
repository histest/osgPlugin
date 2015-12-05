#ifndef __DISTRIBUTE_EARTH_H
#define __DISTRIBUTE_EARTH_H

#include "UserMechanics.h"

#define dimCS 71
#define dimGH 11

struct SOrbitalElements
{
	double SemimajorAxis;
	double Eccentricity;
	double Inclination;
	double AscendingNode;
	double Argument;
	double TrueAnormaly;
};
struct SOrbitalVectors
{
	double SemimajorAxis;
	double EccentricityX;
	double EccentricityY;
	double InclinationX;
	double InclinationY;
	double TrueAnormaly;
};

struct SGeographyElements
{
	double Altitude;
	double Longitude;
	double Latitude;
};

class CEarthAgent
{
public:
	double			Rei[3][3];
private:
	double			Cnm[dimCS][dimCS],Snm[dimCS][dimCS];
	double			Gnm[dimGH][dimGH],Hnm[dimGH][dimGH];
	double			k[dimGH][dimGH];
	double			C[47][6];
	double			T[8][2];
	double			P[18][2];
	double			g10,g11,h11;
	inline unsigned short	Index(const unsigned char n, const unsigned char m);
private:
	double			Pnm[dimCS][dimCS+1],U[dimCS],V[dimCS];
	double			sp[dimGH],cp[dimGH],pp[dimGH];
	double			MPnm[dimGH][dimGH],MDPnm[dimGH][dimGH];

public:
				CEarthAgent();
	void 		Attraction(const unsigned char N,const unsigned char M, const double X[3],double Acceleration[3]);
	void		SpinVelocity(const double R[3], double V[3]);
	void		GetAirVelocity(const double R[3], const double V[3],double Vair[3]);
	double		AirDensity(const double height);
	double		MachVelocity(const double height);
	double		Temperature(const double height);
	double		Pressure(const double height);
	void		MagneticDensity(const unsigned char N, const unsigned char M, const double X[3], double B[3]);
	void		OrbitElementToRectangularElement(const struct SOrbitalElements &xOrbit,double R[3],double V[3]);
	void		RectangularElementToOrbitElement(const double R[3],const double V[3],struct SOrbitalElements &xOrbit);
	void		RectangularElementToOrbitVector(const double R[3],const double V[3],struct SOrbitalVectors &xx);
	void 		Track(const double R[3],struct SGeographyElements &Tracker);
	void 		TrackToPosition(const struct SGeographyElements &Tracker,double xP[3]);
	void 		EarthFixedPositionToHLB(const double R[3],struct SGeographyElements &HLB);
	void 		HLBToEarthFixedPosition(const struct SGeographyElements &HLB,double xP[3]);
	void 		InertiaToFixed(const double Xinertia[3],double xPosition[3]);
	void 		FixedToInertia(const double XPosition[3],double Xinertia[3]);
	void		GetFixedToLocalMatrix(const struct SGeographyElements &HLB,double REe[3][3]);
	void		GetInertiaToLocalMatrix(const struct SGeographyElements &HLB,double REi[3][3]);
	void		GetInertiaToFixedMatrix(double Rei[3][3]);
	void		QuaternionToLocalEularAngle(const struct SGeographyElements &HLB, const double Q[4], struct SEularAngle &Eular);
	void		LocalEularAngleToQuaternion(const struct SGeographyElements &HLB, const struct SEularAngle &Eular, double Q[4]);
	void		InertiaToGPS(const double R[3], const double V[3],struct SGeographyElements &HLB, struct SSphereElements &VAE);
	void		GPSToInertia(const struct SGeographyElements &HLB, const struct SSphereElements &VAE,double R[3],double V[3]);

	inline double FactorialDivision(unsigned short a,unsigned short b);

};

#endif
