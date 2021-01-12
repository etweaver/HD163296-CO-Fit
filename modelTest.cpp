//modelTest.cpp
#include <cassert>
#include <chrono>
#include <cmath>
#include <complex>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <ostream>
#include <random>
#include <utility>
#include <vector>
#include "vcl/vectorclass.h"
#include "vcl/vectormath_trig.h"

#include "diskPhysics.h"
#include "geometry.h"
#include "grid.h"
#include "image.h"
#include "ParameterSet.h"
#include "radmcInterface.h"

struct newDensity {
	double Sigma0;
	double rc;
	double h0;
	double P;//density index
	double S;//scale height index
	double p1, p2, p3, p4; //ring positions (au)
	double d1, d2, d3, d4; //ring depths (0 to 1)
	double w1, w2, w3, w4; //ring widths (au)

	newDensity(): Sigma0(0), rc(0), h0(0), P(0), S(0), p1(0), p2(0), p3(0), p4(0), d1(0), d2(0), d3(0), d4(0), w1(0), w2(0), w3(0), w4(0) { }
	newDensity(const newDensity& other): Sigma0(other.Sigma0), rc(other.rc), h0(other.h0),
	P(other.P), S(other.S), p1(other.p1), p2(other.p2), p3(other.p3), p4(other.p4), d1(other.d1), d2(other.d2), 
	d3(other.d3), d4(other.d4), w1(other.w1), w2(other.w2), w3(other.w3), w4(other.w4) {  }
	newDensity(double Sigma0, double rc, double h0, double P, double S, double p1, double p2, double p3, double p4,
	double d1, double d2, double d3, double d4, double w1, double w2, double w3, double w4):
	Sigma0(Sigma0), rc(rc), h0(h0), P(P), S(S), p1(p1), p2(p2), p3(p3), p4(p4), d1(d1), d2(d2), d3(d3), 
	d4(d4), w1(w1), w2(w2), w3(w3), w4(w4) {	}

	newDensity& operator= (const newDensity& other){
		Sigma0=other.Sigma0; rc=other.rc; h0=other.h0; P=other.P; S=other.S;
		p1=other.p1; p2=other.p2; p3=other.p3; p4=other.p4; d1=other.d1; d2=other.d2; 
		d3=other.d3; d4=other.d4; w1=other.w1; w2=other.w2; w3=other.w3; w4=other.w4;
		return *this;
	}

	double surfaceMassDensity(const double r) const{
		return (Sigma0*pow(r/rc, -P)) * exp(-(pow(r/rc, 2-P)));
	}

	double scaleHeight(const double r) const{
		return (h0*pow(r/rc,S));
	}

	double operator()(double r, double theta, double phi) const{
		double r_cyl=r*sin(theta);
		double z=r*cos(theta);
		double h=scaleHeight(r_cyl);
		double ring1=(1-d1*(gaussianNotNorm(r_cyl,p1,w1)));
		double ring2=(1-d2*(gaussianNotNorm(r_cyl,p2,w2)));
		double ring3=(1-d3*(gaussianNotNorm(r_cyl,p3,w3)));
		double ring4=(1-d4*(gaussianNotNorm(r_cyl,p4,w4)));

		return(((surfaceMassDensity(r_cyl))/(sqrt(2*pi)*h)) * exp(-z*z/(2*h*h))*ring1*ring2*ring3*ring4);
	}
	
	//for a given point in the disk, return the column density ABOVE that position
	//needed to calculate the CO photodisociation
	//Here, the z that matters is the absolute value, since we want this symmetric
	//about the midplane;
	double colDensAbove(double r, double theta, double phi) const{
		double r_cyl=r*sin(theta);
		double z=abs(r*cos(theta));
		double result=surfaceMassDensity(r_cyl)/2;
		double arg=z/(sqrt(2)*scaleHeight(r_cyl));
		result *= erfc(arg);
		return result;
	}
};

//given a final set of parameters, print the file to a fits file and make the uv tables
void printModel(const grid<newDensity>& g1, image& im, const double PA, const std::vector<uvPoint>& data, ThreadPool& pool){
	im.propagate(g1,grid<newDensity>::continuum, pool);
	
	double totalFlux1=0;
	for(int i=0; i<im.hpix; i++){
		for(int j=0; j<im.vpix; j++){
			totalFlux1+=im.data[0][i][j];
		}
	}
	//im.data[0][1400][850]=0.05; //indeces are data[f][y][x]
	std::cout << "total flux: " << totalFlux1 << " jy" << std::endl;
	im.printToFits("model.fits");
	fourierImage FFTs=FFTDifferent(im);
	//FFTs.printToFits("fft.fits");
	double result=0;
	
	//std::ofstream UVoutfile("uvtableSim.txt");
	std::ofstream UVoutfile1("uvtableData.txt");
	std::ofstream UVoutfile2("uvtableResid.txt");
	//UVoutfile1.precision(12);
	//UVoutfile2.precision(12);
	//double inc=(90-22.004)*(pi/180); double PA=(110.3-90)*pi/180;
	for(int i=0;i<data.size();i++){
		std::pair<double,double> interpPair = FFTs.interp(0,data[i].u, data[i].v);
		UVoutfile1 << /*data[i].u << "\t" << data[i].v << "\t" <<*/ interpPair.first << "\t" << interpPair.second << std::endl;
		double diffReal=data[i].real-interpPair.first;
		double diffImag=data[i].imaginary-interpPair.second;
		UVoutfile2 << diffReal << "\t" << diffImag << std::endl;
		//temporary: need to deproject the uv data.
		//double rotU=data[i].u*cos(PA)-data[i].v*sin(PA);
		//double rotV=data[i].u*sin(PA)+data[i].v*cos(PA);
		//rotV*=cos(inc);
		//UVoutfile << sqrt(rotU*rotU + rotV*rotV) << "\t" << sqrt(data[i].u*data[i].u + data[i].v*data[i].v) << "\t" << interpPair.first << "\t" << interpPair.second << std::endl;
		//UVoutfile << data[i].u << "\t" << data[i].v << "\t" << sqrt(data[i].u*data[i].u+data[i].v*data[i].v) << "\t" << interpPair.first << "\t" << interpPair.second << std::endl;
		result+=(diffReal*diffReal + diffImag*diffImag)*data[i].weight;
	}
	UVoutfile1.close();
	UVoutfile2.close();
}

double chiSquaredStupid(image& model, const image& obs, beam& bm, unsigned int index, const double dx, const double dy){
	//first we need to convolve and offset the model
	fourierImage ffts=FFTDifferent(model);
	ffts.offset(dx,dy);
	for(int f=0;f<ffts.frequencies.size();f++){
		for(int i=0; i<ffts.vpix;i++){
			for(int j=0; j<ffts.hpix;j++){
				ffts.realPart[f][i][j]*=bm.realPart[i][j];
				ffts.imaginaryPart[f][i][j]*=bm.realPart[i][j];
			}
		}	
	}
	//ffts.printToFits("FFTs.fits");
	fourierImage back=backTransform(ffts);
	
	image final=model;
	final.data=back.realPart;
	std::string name="modelConvolved";
	name += std::to_string(index);
	name+=".fits";
	final.printToFits(name);
	
	//now to actually count the chi^2
	double sum=0;
	for(int i=0;i<final.hpix;i++){
		for(int j=0;j<final.vpix;j++){
			double realPoint=obs.data[index][i][j];
			double modelPoint=final.data[0][i][j];
			sum+=(realPoint-modelPoint)*(realPoint-modelPoint);
		}
	}
	//and just for good measure, the residuals
	/*image resid=obsFixed;
	for(int i=0;i<resid.hpix;i++){
		for(int j=0;j<resid.vpix;j++){
			resid.data[0][i][j]-=final.data[0][i][j];
		}
	}
	resid.printToFits("resid.fits");*/
	return sum;
}

int main(int argc, char* argv[]){
	//TODO: Figure out why fitsExtract reads the image in upside down
	//I know this is terrible
	//flip it here before we do anything with it
	//gotta stitch the fits files together, and take care of flipping the inputs vertically
	image obsPart1=fitsExtract("frame10.fits",{2.30538e11},101*3.0857e+18);//~-3.2km/s
	image obsPart2=fitsExtract("frame15.fits",{2.30538e11},101*3.0857e+18);//~-1.6km/s
	image obsPart3=fitsExtract("frame20.fits",{2.30538e11},101*3.0857e+18);//~0km/s
	image obsPart4=fitsExtract("frame25.fits",{2.30538e11},101*3.0857e+18);//~1.6km/s
	image obsPart5=fitsExtract("frame30.fits",{2.30538e11},101*3.0857e+18);//~3.2km/s
	image obs(obsPart1.vpix, obsPart1.hpix, obsPart1.width, obsPart1.height, {0,0,0,0,0}, 0, 0, obsPart1.distance);
	for(int i=0;i<obs.vpix;i++){
		for(int j=0;j<obs.hpix;j++){
			obs.data[0][i][j]=obsPart1.data[0][obs.vpix-i-1][j];
			obs.data[1][i][j]=obsPart2.data[0][obs.vpix-i-1][j];
			obs.data[2][i][j]=obsPart3.data[0][obs.vpix-i-1][j];
			obs.data[3][i][j]=obsPart4.data[0][obs.vpix-i-1][j];
			obs.data[4][i][j]=obsPart5.data[0][obs.vpix-i-1][j];
		}
	}
	std::cout << "data read in" << std::endl;
	
	//read in the fit parameters from the command line
	if(argc != 11){
		std::cout << "error: arguments" << std::endl;
		exit(1);
	}
	double logSigGas, rcgas, Pgas, Sgas, h0gas, logSigDust, rcdust, Pdust, h0dust, Sdust; //don't need mStar and dF here
	double mStar=2.05; double deltaF=-37000; //hardcode these for now
	logSigGas=strtod(argv[1],NULL);
	rcgas=strtod(argv[2],NULL);
	Pgas=strtod(argv[3],NULL);
	h0gas=strtod(argv[4],NULL);
	Sgas=strtod(argv[5],NULL);
	logSigDust=strtod(argv[6],NULL);
	rcdust=strtod(argv[7],NULL);
	Pdust=strtod(argv[8],NULL);
	h0dust=strtod(argv[9],NULL);
	Sdust=strtod(argv[10],NULL);
	
	//the radmc interface assumes it's getting a whole ensemble to average.
	//here the "ensemble" is one point, so it looks a little dumb but works
	std::vector<std::vector<double> > ensemble;
	std::vector<double> ensemblePoint;
	ensemblePoint.push_back(logSigGas);
	ensemblePoint.push_back(rcgas*AU);
	ensemblePoint.push_back(Pgas);
	ensemblePoint.push_back(Sgas);
	ensemblePoint.push_back(h0gas*AU);
	ensemblePoint.push_back(mStar*mSun);
	ensemblePoint.push_back(logSigDust);
	ensemblePoint.push_back(rcdust*AU);
	ensemblePoint.push_back(Pdust);
	ensemblePoint.push_back(Sdust);
	ensemblePoint.push_back(h0dust*AU);
	ensemble.push_back(ensemblePoint);
	
	densityEnsemble radmcEnsemble(ensemble);
	radmcEnsemble.outputToRADMC(110,100);
	system("./calcTemp.sh");
	std::cout << "temperature model generated" << std::endl;
	
	//now to actually make the model image
	ThreadPool pool(50);
	double deltaV[5]={-3.2,-1.6,0,1.6,3.2};
	
	double inc, PA, p1, p2, p3, d1, d2, d3, w1, w2, w3;
	inc=0.768; PA=2.307;
	p1=48.23*AU; p2=85.37*AU; p3=98.89*AU;
	d1=0.99; d2=0.96; d3=-0.904916;
	w1=8.5*AU; w2=5.67*AU; w3=3.246*AU;
	
	vect pos(10000*AU*cos(inc),0,10000*AU*sin(inc));
	newDensity gasdens(pow(10, logSigGas), rcgas*AU, h0gas*AU, Pgas, Sgas, 0, 0, 0,0, 0, 0, 0, 0, 0, 0, 0, 0);
	newDensity dustDens(pow(10,logSigDust),rcdust*AU,h0dust*AU,Pdust,Sdust,p1,p2,p3,0,d1,d2,d3,0,w1,w2,w3,0);
	std::ifstream gridStream("radmc/amr_grid.inp");
	std::ifstream temperatureStream("radmc/dust_temperature_phi0.ascii");
	grid<newDensity> g(0.01*AU,1000*AU,60*(pi/180),120*(pi/180),0, 6.28318530717,mStar*mSun, false, true, gridStream,
		temperatureStream, "radmc/dust_density.inp","radmc/dustopac.txt", gasdens, true, 0);
	g.dustDens=dustDens; //dust structure
	g.dens=gasdens;
	g.freezeout=true;
	image im(500, 500, 1312.990851*AU, 1312.990851*AU, {0,0,0,0}, inc, PA, 101*3.0857e+18);
	im.position=pos;
	
	beam bm(obs,(0.104/(2*log(2))), (0.095/(2*log(2))),-80.216*pi/180);
	
	unsigned int fIndex=3;
	
	unsigned int nFreqs=1;
	std::vector<double> globalFrequencies;
	double centfreq=2.30538e11;
	double frequency=centfreq+(deltaV[fIndex]*1e5/c*centfreq);
    
	std::vector<double> frequencies;
	double startfreq=frequency;
	double freqStep=1e5;
	//double freqStep=0;
	frequencies.push_back(startfreq-(freqStep*0.875));
	//frequencies.push_back(startfreq-(freqStep*0.625));
	frequencies.push_back(startfreq-(freqStep*0.375));
	//frequencies.push_back(startfreq-(freqStep*0.125));
	//frequencies.push_back(startfreq+(freqStep*0.125));
	frequencies.push_back(startfreq+(freqStep*0.375));
	//frequencies.push_back(startfreq+(freqStep*0.625));
	frequencies.push_back(startfreq+(freqStep*0.875));
	std::cout.precision(10);
	for(int i=0;i<frequencies.size();i++){ 
		frequencies[i]+=deltaF;
		//std::cout << frequencies[i] << std::endl;
	}
	im.frequencies=frequencies;
	
	im.propagate(g,grid<newDensity>::normal, pool);

	std::vector<double> finalfreq; finalfreq.push_back(centfreq);
	image finalIm(im.vpix, im.hpix, im.width, im.height, finalfreq, inc, PA, im.distance);
	for(int i=0; i<im.hpix; i++){
		for(int j=0; j<im.vpix; j++){
			finalIm.data[0][j][i]=0;
			for(int f=0; f<frequencies.size(); f++){
				finalIm.data[0][j][i]+=im.data[f][j][i]/frequencies.size();
			}
		}
	}
	double chi2=chiSquaredStupid(finalIm,obs,bm,fIndex,0,0);
	std::cout << "chi squared: " << chi2 << std::endl;
	
	return 0;
}
