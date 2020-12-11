//modelEval.cpp
//given a set of parameters, make and evaluate a model disk

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
#include "HTTPRequests.h"
#include "image.h"
#include "ParameterSet.h"
#include "worker.h"

struct diskParams{
	double inc, PA, logSig0, rc, P, h0dust, Sdust, h0gas, Sgas, p1,p2,p3, d1, d2, d3, w1, w2, w3, mStar, deltaF, dx, dy;
	diskParams()=default;
	diskParams(double inc, double PA, double logSig0,double rc,double P,double h0dust,double Sdust,
	double h0gas,double Sgas, double p1, double p2, double p3, double d1,
	double d2,double d3,double w1,double w2,double w3,double mStar,double deltaF,double dx,double dy):
	inc(inc), PA(PA), logSig0(logSig0), rc(rc), P(P), h0dust(h0dust), Sdust(Sdust), 
	h0gas(h0gas), Sgas(Sgas), p1(p1), p2(p2), p3(p3), d1(d1), d2(d2), 
	d3(d3), w1(w1), w2(w2), w3(w3), mStar(mStar), deltaF(deltaF), dx(dx), dy(dy)	{	}
};

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
		
	/*double scaleHeight(const double r) const{
	double T100=14.2243;
	double prefactor=sqrt(kboltzmann*T100/gravConst/1.12/mSun/3.819239518e-24)/pow(100*AU,-0.25);
	return prefactor*pow(r,1.25);
	}*/

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

struct uvTable{
	unsigned int nfreqs;
	double** u;
	double** v;
	double** real;
	double** imag;
	double** weight;
	std::vector<std::unique_ptr<double[]> > u_store;
	std::vector<std::unique_ptr<double[]> > v_store;
	std::vector<std::unique_ptr<double[]> > real_store;
	std::vector<std::unique_ptr<double[]> > imag_store;
	std::vector<std::unique_ptr<double[]> > weight_store;
	unsigned int size; //size per channel

	uvTable()=default;

	uvTable(uvTable&& other)=default;

	double* makeAligned(size_t vectorSize, double* unalignedAddress){
		const int alignBy = vectorSize*sizeof(double);
		double* alignedAddress = (double*)(((size_t)unalignedAddress + alignBy - 1) & (-alignBy));
		return alignedAddress;
	}

	void readFile(std::string inFileName, int nchans, unsigned int nLines){
		const size_t vectorSize=4;
		const size_t padding=vectorSize-1;
		nfreqs=nchans;
		u=new double*[nchans];
		v=new double*[nchans];
		real=new double*[nchans];
		imag=new double*[nchans];
		weight=new double*[nchans];
		u_store.resize(nchans);
		v_store.resize(nchans);
		real_store.resize(nchans);
		imag_store.resize(nchans);
		weight_store.resize(nchans);
		size=nLines;
		
		for(int f=0;f<nchans; f++){
			u_store[f].reset(new double[nLines+2*padding]);
			v_store[f].reset(new double[nLines+2*padding]);
			real_store[f].reset(new double[nLines+2*padding]);
			imag_store[f].reset(new double[nLines+2*padding]);
			weight_store[f].reset(new double[nLines+2*padding]);

			u[f]=makeAligned(vectorSize,u_store[f].get());
			v[f]=makeAligned(vectorSize,v_store[f].get());
			real[f]=makeAligned(vectorSize,real_store[f].get());
			imag[f]=makeAligned(vectorSize,imag_store[f].get());
			weight[f]=makeAligned(vectorSize,weight_store[f].get());

			std::ifstream infile(inFileName);
			std::string line;
			double tempu, tempv, tempReal, tempImag, tempWeight;
			for(int i=0; i<nLines; i++){
				std::getline(infile, line);
				std::stringstream stream;
				stream << line;
				stream >> tempu >> tempv >> tempReal >> tempImag >> tempWeight;
				u[f][i]=tempu;
				v[f][i]=tempv;
				real[f][i]=tempReal;
				imag[f][i]=tempImag;
				weight[f][i]=tempWeight;
				//std::cout << tempu << "\t" << tempv << "\t" << tempReal << "\t" << tempImag << "\t" << tempWeight << std::endl;
			}
			infile.close();
			if(nLines%vectorSize){
				for(int i=0;i<vectorSize-nLines%vectorSize;i++){
					u[f][i+nLines]=v[f][i+nLines]=real[f][i+nLines]=imag[f][i+nLines]=weight[f][i+nLines]=0;
				}
			}
		}
		//don't forget to adjust the size per channel based on the additional padding
		//I'm assuming that all channels have the same number of points
		if(nLines%vectorSize)
			size+=vectorSize-nLines%vectorSize;
		
	}
};

//given a final set of parameters, print the file to a fits file and make the uv tables
void printModel(const grid<newDensity>& g1, image& im, const double PA, const std::vector<uvPoint>& data, ThreadPool& pool){
	im.propagate(g1,-PA,grid<newDensity>::continuum, pool);
	
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

double chiSquaredAVX(image& model, const uvTable& data, unsigned int index, const double dx, const double dy){
	fourierImage FFT=FFTDifferent(model);
	FFT.offset(dx,dy);
	double result=0;
	int totalNumber=data.size;

	//std::cout << "starting chi^2 " << index << "\t" << dx << "\t" << dy << std::endl;

	for(int i=0;i<data.size;i+=4){
		Vec4d uPointVec(0);
		Vec4d vPointVec(0);
		Vec4d realVec(0);
		Vec4d imagVec(0);
		Vec4d weightVec(0);
		uPointVec.load_a(data.u[index]+i);
		vPointVec.load_a(data.v[index]+i);
		realVec.load_a(data.real[index]+i);
		imagVec.load_a(data.imag[index]+i);
		weightVec.load_a(data.weight[index]+i);
			
		std::pair<Vec4d,Vec4d> interpPair = FFT.interpAVX(0,uPointVec, vPointVec);
		Vec4d diffReal=realVec-interpPair.first;
		Vec4d diffImag=imagVec-interpPair.second;

		Vec4d resultVec=(diffReal*diffReal + diffImag*diffImag)*weightVec;

		for(int i=0;i<resultVec.size();i++){
			//std::cout << uPointVec[i] << "\t" << vPointVec[i] << "\t" << realVec[i] << "\t" << imagVec[i] << "\t" << weightVec[i] << "\t\t" << diffReal[i] <<std::endl;
			result+=resultVec[i];
		}
	}

	return result;
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

class modelWorker : public Worker{
public:
	modelWorker(std::string serverAddress, std::string serverTempFileURL, image& obs):
		Worker(serverAddress),
		dens(-1, -1, -1, -1, -1, -1, -1, -1, 0, -1, -1, -1, 0, -1, -1, -1, 0),	//-1 is just a placeholder for things we'll put in once we get the parameters
		g(0.01*AU,1000*AU,60*(pi/180),120*(pi/180),0, 6.28318530717, 0, false, true, "diskmodels/HD163296/amr_grid.inp",
			"diskmodels/HD163296/dust_temperature_phi0.ascii", "diskmodels/HD163296/dust_temperature_phi0.ascii",
			"diskmodels/HD163296/dustopac.txt", dens, false, 0),
		im(500, 500, 1312.990851*AU, 1312.990851*AU , origin, origin, 2.30538e+11, 8,{0,0,0,0,0,0,0,0}, 101*3.0857e+18),
		obs(obs),
		bm(obs,(0.104/(2*log(2))), (0.095/(2*log(2))),-80.216*pi/180)
	{
		//the grid has separate dust and gas density structures, but it didn't used to, and the constructors don't yet reflect this
		//now the two structures are nearly the same, so we only need to change a few things
		g.dustDens=g.dens;
		g.dustDens.h0=-1; g.dustDens.S=-1;
		im.RA = 269.0886636250;
		im.DEC = -21.9562676325;
		vect pos(-1,0,-1);
		im.position=pos;
		tempFileURL=serverTempFileURL;
	}
			
protected:
	virtual double compute(const std::vector<double>& params) override {
		ThreadPool pool(1); //this will be taken out once things are working
		double deltaV[5]={-3.2,-1.6,0,1.6,3.2};

		//first we need to extract the full set of disk parameters from the input
		unsigned int chainIndex, fIndex;
		//some are hardcoded for now:
		chainIndex=(unsigned int)params[0]; fIndex=(unsigned int)params[1];
		
		std::cout << "work item " << chainIndex << "\t" << fIndex << std::endl;

		double logSig0gas=params[2]; double rcgas=params[3]; double  Pgas=params[4]; 
		double  h0gas=params[5]; double Sgas=params[6]; double mStar=params[7]; 
		double deltaF=params[8]; double logSig0dust=params[9]; double rcdust=params[10]; 
		double  Pdust=params[11]; double  h0dust=params[12]; double Sdust=params[13];
		
		double inc, PA, p1, p2, p3, d1, d2, d3, w1, w2, w3;
		inc=-0.768; PA=2.307;
		p1=48.23*AU; p2=85.37*AU; p3=98.89*AU;
		d1=0.99; d2=0.96; d3=-0.904916;
		w1=8.5*AU; w2=5.67*AU; w3=3.246*AU;
		
		vect pos(10000*AU*cos(inc),0,10000*AU*sin(inc));
		
		//we're almost ready to do the imaging, but we need to get the temperature file
                //from the server
                httpRequests::Response resp=httpRequests::httpGet(tempFileURL);
                if(resp.status!=200){
                        std::cout << "Error: could not retrieve temperature file" << std::endl;
			std::cout << resp.status << std::endl;
                        exit(1);
                }
		std::istringstream temperatureStream(resp.body);
		std::ifstream gridStream("radmc/amr_grid.inp");
	
		newDensity gasdens(pow(10, logSig0gas), rcgas*AU, h0gas*AU, Pgas, Sgas, 0, 0, 0,0, 0, 0, 0, 0, 0, 0, 0, 0);
		newDensity dustDens(pow(10,logSig0dust),rcdust*AU,h0dust*AU,Pdust,Sdust,p1,p2,p3,0,d1,d2,d3,0,w1,w2,w3,0);
		grid<newDensity> g(0.01*AU,1000*AU,60*(pi/180),120*(pi/180),0, 6.28318530717,mStar*mSun, false, true, gridStream,
			temperatureStream, "radmc/dust_density.inp","radmc/dustopac.txt", gasdens, true, 0);
		g.dustDens=dustDens; //dust structure
		g.dens=gasdens;
		g.freezeout=true;
		image im(500, 500, 1312.990851*AU, 1312.990851*AU , origin, origin, 2.30538e+11, 4,{0,0,0,0}, 101*3.0857e+18);
		im.position=pos;
				
		unsigned int nFreqs=1;
		std::vector<double> globalFrequencies;
		double frequency=im.centfreq+(deltaV[fIndex]*1e5/c*im.centfreq);
        
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
		im.freqbins=frequencies.size();
		
		im.propagate(g,PA,grid<newDensity>::normal, pool);
		
		std::vector<double> finalfreq; finalfreq.push_back(im.centfreq);
		image finalIm(im.vpix, im.hpix, im.width, im.height , origin, origin, im.centfreq, 1, finalfreq, im.distance);
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
		return chi2*25;
	}

private:
	newDensity dens;
	grid<newDensity> g;
	image im, finalImg;
	image obs;
	beam bm;
	double chi2;
	std::string tempFileURL;
};

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
	image obs(obsPart1.vpix, obsPart1.hpix, obsPart1.width, obsPart1.height, obsPart1.position, obsPart1.target, 
		2.30538e11, 5, {0,0,0,0,0}, obsPart1.distance);
	for(int i=0;i<obs.vpix;i++){
		for(int j=0;j<obs.hpix;j++){
			obs.data[0][i][j]=obsPart1.data[0][obs.vpix-i-1][j];
			obs.data[1][i][j]=obsPart2.data[0][obs.vpix-i-1][j];
			obs.data[2][i][j]=obsPart3.data[0][obs.vpix-i-1][j];
			obs.data[3][i][j]=obsPart4.data[0][obs.vpix-i-1][j];
			obs.data[4][i][j]=obsPart5.data[0][obs.vpix-i-1][j];
		}
	}
	
	std::string serverAddress="obelix.rice.edu:15002";
	std::string serverTempFileURL="obelix.rice.edu:15001/radmc/dust_temperature_phi0.ascii";
	modelWorker testWorker(serverAddress, serverTempFileURL, obs);
	testWorker.run();
	
	return 0;
}
