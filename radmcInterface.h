//radmcInterface.h
//handles outputting the ensemble average density into the input files for RADMC3d

#include <cassert>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <limits>
#include <ostream>
#include <random>
#include <string>
#include <utility>
#include <vector>

//the density ensemble is a subset of the actual ensemble.
//The full ensemble is defined in driver.cpp, and contains more data as well 
//as the entire mcmc interface, which don't need here.
//The data here is a subset of the actual ensemble.coords.
//The coords will be what we need to generate files like in andrea's param2radmc
//script, so sigma0gas, rc P, S, and h0gas, for both gas and dust.
//We'll need to calculate mgas from sigma0gas and the others
//The other parameters, rmin, rmax, etc, we'll hardcode to match what the grid
//in modelEval will be using.
struct densityEnsemble{
public:
	densityEnsemble():coords(), Sigma0Avg(0), rcAvg(0), PAvg(0), SAvg(0), h0Avg(0), mStarAvg(0),
			Sigma0AvgDust(0), rcAvgDust(0), PAvgDust(0), SAvgDust(0), h0AvgDust(0){	}
	densityEnsemble(std::vector<std::vector<double> > coords, double Sigma0Avg, double rcAvg, double PAvg, double SAvg, 
		double h0Avg, double mStarAvg, double Sigma0AvgDust, double rcAvgDust, double PAvgDust, double SAvgDust, double h0AvgDust):
		coords(coords), Sigma0Avg(Sigma0Avg), rcAvg(rcAvg), PAvg(PAvg), SAvg(SAvg), h0Avg(h0Avg), mStarAvg(mStarAvg),
		Sigma0AvgDust(Sigma0AvgDust), rcAvgDust(rcAvgDust), PAvgDust(PAvgDust), SAvgDust(SAvgDust), h0AvgDust(h0AvgDust){}
	densityEnsemble(const densityEnsemble& other):
		coords(other.coords), Sigma0Avg(other.Sigma0Avg), rcAvg(other.rcAvg), PAvg(other.PAvg), 
		SAvg(other.SAvg), h0Avg(other.h0Avg), mStarAvg(other.mStarAvg), Sigma0AvgDust(other.Sigma0AvgDust), 
		rcAvgDust(other.rcAvgDust), PAvgDust(other.PAvgDust), SAvgDust(other.SAvgDust), h0AvgDust(other.h0AvgDust)	{}
		
	//The special constructor calculates the averages based on an input ensemble vector
	densityEnsemble(std::vector<std::vector<double> >& coords):coords(coords){
		double sig0=0;	double sig0Dust=0;
		double rc=0;	double rcDust=0;
		double P=0; 	double PDust=0;
		double S=0;		double SDust=0;
		double h0=0;	double h0Dust=0;
		double mStar=0;	
		for(int i=0;i<coords.size();i++){
			sig0+=coords[i][0];
			rc+=coords[i][1];
			P+=coords[i][2];
			S+=coords[i][3];
			h0+=coords[i][4];
			mStar+=coords[i][5];
			sig0Dust+=coords[i][6];
			rcDust+=coords[i][7];
			PDust+=coords[i][8];
			SDust+=coords[i][9];
			h0Dust+=coords[i][10];
		}
		sig0/=(double)coords.size();	sig0Dust/=(double)coords.size();
		rc/=(double)coords.size();		rcDust/=(double)coords.size();
		P/=(double)coords.size();		PDust/=(double)coords.size();
		S/=(double)coords.size();		SDust/=(double)coords.size();
		h0/=(double)coords.size();		h0Dust/=(double)coords.size();
		mStar/=(double)coords.size();	
		Sigma0Avg=pow(10,sig0);			Sigma0AvgDust=pow(10,sig0Dust);
		rcAvg=rc;						rcAvgDust=rcDust;
		PAvg=P;							PAvgDust=PDust;
		SAvg=S;							SAvgDust=SDust;
		h0Avg=h0;						h0AvgDust=h0Dust;
		mStarAvg=mStar;
	}
	
	void outputToRADMC(const unsigned int nr, const unsigned int ntheta){
		//TODO: Move these out of here
		const double pi=4*atan(1);
		const double AU=1.496e13;
		const double amu=1.670539040e-24;//g
		const double gravConst=6.67259e-8;
		//first I'm going to hardcode the radmc disk dimensions and bin numbers and things
		double rmin=0.1*AU;
		double rmax=600*AU;
		//unsigned int nr=600;
		
		//radmc can be set to assume midplane symmetry, which I don't generally assume,
		//but have literally always used, so we'll use it here and go from 60-90°
		double thetamin=50.0*pi/180.0;
		double thetamax=pi/2;
		//unsigned int ntheta=300;
		
		double numberRatio=1e-5;//number of 12CO molecules per H2
		double dustRatio=1e-2;//***MASS*** ratio of dust to H2
	
		//first we need to define the grid, in amr_grid.inp
		//first the radmc header stuff
		std::ofstream gridfile("radmc/amr_grid.inp");
		gridfile.precision(10);
		gridfile << 1   << std::endl;   // iformat
		gridfile << 0   << std::endl;   // Grid style
		gridfile << 100 << std::endl;   // coordinate system; if <100 the system is cartesian, 100 <= coord < 200 the system is spherical, 
	    	                   // 200 <= coord < 300 the system is cylindrical
		gridfile << 0   << std::endl;   // grid info: 1 = verbose, 0 = standard
		gridfile << 1   << " " << 1 << " " << 0 << std::endl; // include dimension x, y, z: 1 = active, 0 = not active
		gridfile << nr-1  << " " << ntheta-1 << " " << 1 << std::endl; // number of grid cells along r, theta, and phi
		
		//now the actual bin numbers
		//first r, then theta
		for(int i=0;i<nr;i++){
			gridfile << rmin*pow(rmax/rmin, (double)i/(double)(nr-1)) << std::endl;
		}
		for(int i=0;i<ntheta;i++){
			gridfile << thetamin + (double)i*((thetamax-thetamin)/(double)(ntheta-1)) << std::endl;
		}
		gridfile << 0 << std::endl;
		gridfile << 2*pi << std::endl; //and phi is just 0 to 2pi
		gridfile.close();
		
		//now we need to specify the co density file, and gas velocity file
		std::ofstream molFile("radmc/numberdens_co.inp");
		std::ofstream dustFile("radmc/dust_density.inp");
		std::ofstream velFile("radmc/gas_velocity.inp");
		molFile.precision(10);
		velFile.precision(10);		
		
		//the molFile we'll also give the parameters to, so that we can tell what parameter set generated a given temperature
		molFile << "#" << Sigma0Avg << "\t" << rcAvg << "\t" << PAvg << "\t" << SAvg << "\t" << h0Avg << "\t" << mStarAvg << std::endl;
	    	molFile << 1 << std::endl;
	    	molFile <<  (nr-1)*(ntheta-1) << std::endl; //number of bins total

		dustFile << 1 << std::endl;
		dustFile <<  (nr-1)*(ntheta-1) << std::endl; //number of bins total
		dustFile << 1 << std::endl;

	    	velFile << 1 << std::endl;
	    	velFile << (nr-1)*(ntheta-1) << std::endl;
	
		//numberdens_co.inp is the number density in each bin
		//gas_velocity.inp is the velocity in each bin, which for this case
		//is just keplerian
		double maxHeight=5;	//max number of scale heights
		for(int i=0;i<ntheta-1;i++){
			for(int j=0;j<nr-1;j++){
				double rLowerBinEdge=rmin*pow(rmax/rmin, (double)j/(double)(nr-1));
                		double rUpperBinEdge=rmin*pow(rmax/rmin, (double)(j+1)/(double)(nr-1));
				double tLowerBinEdge=thetamin + (double)i*((thetamax-thetamin)/(double)(ntheta-1));
				double tUpperBinEdge=thetamin + (double)(i+1)*((thetamax-thetamin)/(double)(ntheta-1));
				double rSphere=(rUpperBinEdge+rLowerBinEdge)/2.0;
				double theta=(tUpperBinEdge+tLowerBinEdge)/2.0;
				double rCyl=rSphere*sin(theta);
				double z=rSphere*cos(theta);
				
				double surf=Sigma0Avg*pow(rCyl/rcAvg,-PAvg);//* exp(-(pow(rCyl/rcAvg, 2-PAvg)));; //surface mass density
				double h=h0Avg*pow(rCyl/rcAvg,SAvg);	//scale height
				
				double massDens=surf/sqrt(2*pi)/h * exp(-z*z/2.0/h/h);
				
				//std::cout << tLowerBinEdge << "\t" << tUpperBinEdge << "\t" << (double)j*((thetamax-thetamin)/(double)(ntheta-1)) << "\t" << theta << std::endl;
				//std::cout << rCyl/AU << "\t" << z/AU << "\t" << h/AU << "\t" << surf << "\t" << massDens << std::endl;
				dustFile /*<< rSphere/AU << "\t" << theta << "\t"*/ << massDens*0.01 << std::endl;
				if(z<maxHeight*h){
					molFile << massDens/28.0/amu*numberRatio << std::endl;
				} else{
					molFile << 0.0 << std::endl;
				}
				velFile << 0.0 << "\t" << 0.0 << "\t" << sqrt(gravConst*mStarAvg/rCyl) << std::endl;
			}
		}
		molFile.close();
		velFile.close();
		dustFile.close();
		std::cout << "printed amr_grid.inp, numberdens_co.inp, and gas_velocity.inp" << std::endl;
		
	}
	
private:
	std::vector<std::vector<double> > coords;
	double Sigma0Avg, rcAvg, PAvg, SAvg, h0Avg, mStarAvg;
	double Sigma0AvgDust, rcAvgDust, PAvgDust, SAvgDust, h0AvgDust;
};
