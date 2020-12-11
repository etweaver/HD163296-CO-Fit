#include "driver.h"
#include "ParameterSet.h"
#include <cstdlib>
#include <random>
#include <sstream>
#include <string>
#include <condition_variable>

#include "radmcInterface.h"

struct moveProposal{
	std::vector<double> coordinates;
	double logAcceptanceFactor;

	moveProposal()=default;
	moveProposal(std::size_t n):coordinates(n){}
	moveProposal(moveProposal&& mp):
	coordinates(std::move(mp.coordinates)),logAcceptanceFactor(mp.logAcceptanceFactor){}
};

struct ensembleMember{
	std::vector<double> coords;
	double currentValue;
	std::vector<double> proposedCoords;
	std::vector<double> partialValues;
	double logAcceptanceFactor;
	
	ensembleMember()=default;
	ensembleMember(unsigned int nFrequencies):partialValues(nFrequencies,0) {	}

	void proposeJump(moveProposal mp){
		proposedCoords=std::move(mp.coordinates);
		logAcceptanceFactor=mp.logAcceptanceFactor;
	}

	template<typename RNG>
	bool decideJump(std::uniform_real_distribution<double>& acceptDist, RNG& rng, ParameterSet& params){
		double proposedValue=0; 
		for(auto v : partialValues)
			proposedValue+=v;
		double logRatio=(currentValue-proposedValue)+logAcceptanceFactor; 
		bool accept=(logRatio>=0) || (logRatio>log(acceptDist(rng)));
		if(!params.inBounds(proposedCoords)){
			accept=false;
			std::cout << "step rejected: out of bounds" << std::endl;
		}
		if(accept){
			std::cout << "move accepted: old value: " << currentValue << ", new value: " << proposedValue << ", logAcceptanceFactor: " << logAcceptanceFactor << std::endl;
			std::copy(proposedCoords.begin(),proposedCoords.end(),coords.begin());
			currentValue=proposedValue;
		}
		return accept;
	}
};

struct stretchMove{

	//square root distribution
	struct sqrtDist{
		double range; //extends from 1/r to r
		double norm;
		mutable std::uniform_real_distribution<double> uniform; //mutable is needed to override a const later because the
		//uniform distribution has no const call operator. Ugly but necessary.

		sqrtDist(double range): range(range),norm(1/((sqrt(range)-sqrt(1/range)))),uniform(norm*sqrt(1/range),norm*sqrt(range)){
			if(range<=1)
				throw std::domain_error("square_root_distribution requires range>1");
		}

		template<typename RNG>
		double operator()(RNG& rng) const{
			double v=uniform(rng)/(norm);
			return(v*v);
		}
	};

	sqrtDist jumpDist;
	stretchMove():jumpDist(2){}

	template<typename RNG>
	moveProposal operator()(const std::vector<double>& coords, const std::vector<ensembleMember>& ensemble, RNG& rng) const {
		assert(!coords.empty());
		assert(ensemble.size() > 1);
		moveProposal mp(coords.size());

		//choose a random member of the ensemble, but not the one we are already using
		std::uniform_int_distribution<int> idxDist(0,ensemble.size()-1);
		int idx=idx=idxDist(rng);;
		unsigned int maxTrials=100; //if we need to keep trying
		//std::cout << coords.size() << "\t" << ensemble.size() << "\t" << idx << std::endl;
		while(std::equal(coords.begin(),coords.end(),ensemble[idx].coords.begin())){
			idx=idxDist(rng);
			if(!--maxTrials)
				throw std::runtime_error("StretchMove failed too many times to find a distinct ensemble member. "
					"Ensmeble may have become degenerate.");
		}

		//jump distance
		double z=jumpDist(rng);

		//direction
		for(std::size_t i=0; i<coords.size(); i++)
			mp.coordinates[i]=ensemble[idx].coords[i]+z*(coords[i]-ensemble[idx].coords[i]);
		//and the penalty associated with going there
		mp.logAcceptanceFactor=(coords.size()-1)*log(z);

		return mp;
	}
};

template<typename RNG>
double randInRange(double min, double max, RNG& rng){
	std::uniform_real_distribution<double> dist(min,max);
	return dist(rng);
}

template<typename RNG, class Jumper>
class modelingServer : public WorkServerImpl{
public:
	modelingServer()=default;
	modelingServer(unsigned int nSamples, unsigned int nFrequencies, unsigned int ensembleSize,
	const ParameterSet& params, std::size_t rngSeed, bool readFromFile):
	nSamples(nSamples),nFrequencies(nFrequencies),ensembleSize(ensembleSize),counter(0),
	ensembleCounter(0),accepted(0),rejected(0),params(params),rng(rngSeed),firstGeneration(true),
	processor([this](){this->processResultsImpl();}){
		if(readFromFile){
			std::string lineIn;
			firstGeneration=false;
			std::ifstream statefile("finalState.txt");
			int numRead=0;
			while(!statefile.eof() && numRead < ensembleSize){
				ensembleMember member(nFrequencies);
				std::getline(statefile, lineIn);
				if(lineIn[0]=='#')
					continue;
				if(lineIn[0]=='\n')
					continue;
				std::vector<double> ensembleMember;
				double logSigGas, rcgas, Pgas, Sgas, h0gas, mStar, dF, logSigDust, rcdust, Pdust, h0dust, Sdust, like;
				std::stringstream stream;
				stream << lineIn;
				stream >> logSigGas >> rcgas >> Pgas >> Sgas >> h0gas >> mStar >> dF >> logSigDust >> rcdust >> Pdust >> h0dust >> Sdust >> like;

				member.coords = {logSigGas, rcgas, Pgas, Sgas, h0gas, mStar, dF, logSigDust, rcdust, Pdust, h0dust, Sdust};
				member.currentValue=like;
				ensemble.push_back(member);
				numRead++;
				//for(auto i : ensembleMember)		
				//std::cout << i << "\t";
				//std::cout << std::endl;
			}
			statefile.close();
			std::cout << "Ensemble loaded from file" << std::endl;
		}else{
			for(int i=0;i<ensembleSize;i++){
				ensembleMember member(nFrequencies);
				member.coords.push_back(randInRange(-1.25,-1.15,rng));//logSig0
				member.coords.push_back(randInRange(375,385,rng));//rc
				member.coords.push_back(randInRange(0.5,0.55,rng));//P
				member.coords.push_back(randInRange(130,140,rng));//h0gas
				member.coords.push_back(randInRange(1.0,1.1,rng));//Sgas
				member.coords.push_back(randInRange(2.17,2.19,rng));//mStar
				member.coords.push_back(randInRange(3.8e4,4e4,rng));//deltaF*/
				member.coords.push_back(randInRange(-0.2,0,rng));//logSig0dust
				member.coords.push_back(randInRange(80,82,rng));//rcdust
				member.coords.push_back(randInRange(0.3,0.4,rng));//Pdust
				member.coords.push_back(randInRange(1.5,2,rng));//h0dust
				member.coords.push_back(randInRange(0.1,0.5,rng));//Sdust
				ensemble.push_back(member);
			}
		}
		/*std::cout << "starting ensemble: " << std::endl;
		for(int i=0;i<ensembleSize;i++){
		for(int j=0;j<4;j++){
		std::cout << ensemble[i].coords[j] << "\t";	
		}
		std::cout << std::endl;
		}*/
		generateWorkBlock();
		std::cout << "first work block ready" << std::endl;
	}
		
protected:
	virtual void processResult(WorkItem item, double value) override{
		std::lock_guard<std::mutex> lock(resultProcessingLock);
		std::cout << "placing result " << item << " in processing queue" << std::endl;
		resultsToProcess.push(std::make_pair(item,value)); //TODO: maybe value should just be a member of WorkItem?
		needToProcessResults.notify_one();
	}
    
	void processResultsImpl(){
		while(!done()){
			//wait for something to process
			std::unique_lock<std::mutex> lock(resultProcessingLock);
			if(resultsToProcess.empty())
				needToProcessResults.wait(lock);
			
			//grab it
			if(resultsToProcess.empty())
				continue;
			auto [item,value]=resultsToProcess.front();
			resultsToProcess.pop();
			lock.unlock(); //done accessing queue, release lock
			
			//process it
			std::cout << "processing result " << item << std::endl;
			unsigned int chainIndex=(unsigned int)item.parameters[0];
			unsigned int freqIndex=(unsigned int)item.parameters[1];
			//std::cout << "got indeces: " << chainIndex << "\t" << freqIndex << std::endl;
			if(chainIndex > ensembleSize)
				throw std::runtime_error("Error: Bad chain index");
			if(freqIndex > nFrequencies)
				throw std::runtime_error("Error: Bad frequency index");
			//to construct the ensemble from each result, we need to go to
			//the correct spot for the member piece, and add its value to the 
			//value list. The coordinates are already there from when the proposed
			//ensemble was generated.
			ensemble[chainIndex].partialValues[freqIndex]=value;
			ensembleCounter++;
			std::cout << "ensemble counter: " << ensembleCounter.load() << std::endl;

			//if that was the end of a block, do the block processing, then generate a new block
			if(ensembleCounter.load()==ensembleSize*nFrequencies){
				std::cout << "processing finished ensemble" << std::endl;
				if(firstGeneration){
					for(int i=0;i<ensembleSize;i++){
						ensemble[i].currentValue=0;
						for(auto v : ensemble[i].partialValues)
							ensemble[i].currentValue+=v;
					}
					firstGeneration=false;
				}
				else{
					unsigned int acceptedThisRun=0;
					unsigned int rejectedThisRun=0;
					std::uniform_real_distribution<double> acceptDist(0,1);
					for(int i=0;i<ensembleSize;i++){
						if(ensemble[i].decideJump(acceptDist,rng,params)){
							acceptedThisRun++;
						}else{
							rejectedThisRun++;
						}
					}
					std::cout << "acceptedThisRun: " << acceptedThisRun << std::endl;
					std::cout << "rejectedThisRun: " << rejectedThisRun << std::endl;
					accepted+=acceptedThisRun;
					rejected+=rejectedThisRun;
					counter++;
					std::ofstream outfile;
					outfile.open("samples.txt", std::ofstream::app);
					for(int i=0;i<ensembleSize;i++){
						for(int j=0;j<ensemble[i].coords.size();j++){
							outfile << ensemble[i].coords[j] << "\t";
						}
						outfile << ensemble[i].currentValue << std::endl;
					}
					outfile.close();
				}
				ensembleCounter=0;
				std::cout << "done processing finished ensemble" << std::endl;
				
				generateWorkBlock();
			}
		}
	}

	virtual bool done() const override{
		return counter.load()>=nSamples;
	}
			
	//I think that this should generate ensembleSize*nFrequencies work units
	//starting from an ensemble generated by generateNextStep(). Each needs 
	//its ensemble number and frequency number and then to be sent out.
	void generateWorkBlock(){
		std::cout << "generating work block" << std::endl;
		//this is where we will use the ensemble to generate the radmc files
		//later this may move as the code is restructured to fit this step better,
		//for now, it's easiest to stop other processes and use the server computer 
		//to generate the new temperature. Ideally, this could be done in parallel
		//with the fitting on the workers, but for now we'll pause to generate the 
		//temperatures
		
		//first we ned to pack the parameters we need into a vector so that they
		//can be processed into the averages
		std::vector<std::vector<double> > subEnsemble;
		for(int i=0;i<ensembleSize;i++){
			std::vector<double> coords;
			coords.push_back(ensemble[i].coords[7]);
			coords.push_back(ensemble[i].coords[8]);
			coords.push_back(ensemble[i].coords[9]);
			coords.push_back(ensemble[i].coords[11]);
			coords.push_back(ensemble[i].coords[10]);
			coords.push_back(ensemble[i].coords[5]);
			//coords.push_back(-2.108641);
			//coords.push_back(100*AU);
			//coords.push_back(-0.1);
			//coords.push_back(1.4);
			//coords.push_back(11*AU);
			//coords.push_back(2.05*mSun);
			subEnsemble.push_back(coords);
		}
		densityEnsemble radmcEnsemble(subEnsemble);
		radmcEnsemble.outputToRADMC(400,200);
		
		//radmc is all fortran, so there isn't much of a way to run it from here besides system()
		system("./calcTemp.sh");
		std::cout << "temperature model generated" << std::endl;

		//now to generate the next actual work block
		std::vector<WorkItem> workBlock;
		for(int i=0;i<ensembleSize;i++){
			auto proposal=jumper(ensemble[i].coords,ensemble,rng);
			ensemble[i].proposeJump(std::move(proposal));
			for(int j=0; j<nFrequencies; j++){
				WorkItem w({(double)i,(double)j});
				for(int k=0; k<ensemble[i].coords.size(); k++)
					w.parameters.push_back(firstGeneration ? ensemble[i].coords[k] : ensemble[i].proposedCoords[k]);
				workBlock.push_back(w);
			}
		}
		addWork(std::move(workBlock));
	}

private:
	unsigned int nSamples;
	unsigned int nFrequencies;
	unsigned int ensembleSize;
	std::atomic<unsigned int> counter;		//from 0 to nSamples
	std::atomic<unsigned int> ensembleCounter;	//from 0 to ensembleSize*nFrequencies. Just for tracking when each ensemble finishes
	std::vector<ensembleMember> ensemble;
	unsigned int accepted;				//steps that were accepted
	unsigned int rejected;
	ParameterSet params;
	RNG rng;
	Jumper jumper;
	bool firstGeneration;
	//list of all states of the ensemble over time
	//	std::vector<std::vector<std::vector<double>>> ensembleHistory;
    std::mutex resultProcessingLock;
    std::condition_variable needToProcessResults;
	std::queue<std::pair<WorkItem,double>> resultsToProcess;
	std::thread processor;
};

int main(int argc, char* argv[]){
	std::mt19937 rng(137);

	auto randInRange=[&rng](const double min, const double max){
		double range = max-min;
		double randNum=rng()/(double)rng.max();
		return min+randNum*range;
	};

	ParameterSet params;
	/*params.addParameter("inclination");
	params.setParameterLowerLimit("inclination",0); params.setParameterUpperLimit("inclination",pi/2);
	params.addParameter("PA");
	params.setParameterLowerLimit("PA",0); params.setParameterUpperLimit("PA",2*pi);*/
	//gas structure parameters
	params.addParameter("logSigma0");
	params.setParameterLowerLimit("logSigma0",-3); params.setParameterUpperLimit("logSigma0",2);
	params.addParameter("rc"); //units of au
	params.setParameterLowerLimit("rc",50); params.setParameterUpperLimit("rc",500);
	params.addParameter("P"); //unitless
	params.setParameterLowerLimit("P",-1); params.setParameterUpperLimit("P",1.999);
	params.addParameter("h0gas"); //au
	params.setParameterLowerLimit("h0gas",25); params.setParameterUpperLimit("h0gas",500);
	params.addParameter("Sgas"); //unitless
	params.setParameterLowerLimit("Sgas",-1); params.setParameterUpperLimit("Sgas",2);
	params.addParameter("mStar"); //mSun
	params.setParameterLowerLimit("mStar",2.05); params.setParameterUpperLimit("mStar",2.35);
	params.addParameter("deltaF"); //Hz
	params.setParameterLowerLimit("deltaF",-1e7); params.setParameterUpperLimit("deltaF",1e7);
	//dust structure parameters
	params.addParameter("logSigma0dust");
	params.setParameterLowerLimit("logSigma0dust",-3); params.setParameterUpperLimit("logSigma0dust",2);
	params.addParameter("rcdust"); //units of au
	params.setParameterLowerLimit("rcdust",50); params.setParameterUpperLimit("rcdust",500);
	params.addParameter("Pdust"); //unitless
	params.setParameterLowerLimit("Pdust",-1); params.setParameterUpperLimit("Pdust",1.999);
	params.addParameter("h0dust"); //au
	params.setParameterLowerLimit("h0dust",0.01); params.setParameterUpperLimit("h0dust",25);
	params.addParameter("Sdust"); //unitless
	params.setParameterLowerLimit("Sdust",-1); params.setParameterUpperLimit("Sdust",1.25);
	/*params.addParameter("r1"); //au
	params.setParameterLowerLimit("r1",30); params.setParameterUpperLimit("r1",65);
	params.addParameter("r2"); //au
	params.setParameterLowerLimit("r2",70); params.setParameterUpperLimit("r2",100);
	params.addParameter("r3"); //au
	params.setParameterLowerLimit("r3",120); params.setParameterUpperLimit("r3",180);
	params.addParameter("d1"); //unitless
	params.setParameterLowerLimit("d1",0); params.setParameterUpperLimit("d1",1);
	params.addParameter("d2"); //unitless
	params.setParameterLowerLimit("d2",0); params.setParameterUpperLimit("d2",1);
	params.addParameter("d3"); //unitless
	params.setParameterLowerLimit("d3",0); params.setParameterUpperLimit("d3",1);
	params.addParameter("w1"); //au
	params.setParameterLowerLimit("w1",0); params.setParameterUpperLimit("w1",25);
	params.addParameter("w2"); //au
	params.setParameterLowerLimit("w2",0); params.setParameterUpperLimit("w2",25);
	params.addParameter("w3"); //au
	params.setParameterLowerLimit("w3",0); params.setParameterUpperLimit("w3",50);*/
	/*params.addParameter("dx"); //pixels
	params.setParameterLowerLimit("dx",-25); params.setParameterUpperLimit("dx",25);
	params.addParameter("dy"); //pixels
	params.setParameterLowerLimit("dy",-25); params.setParameterUpperLimit("dy",25);*/
	
	modelingServer<std::mt19937,stretchMove> testServer(1200,5,40,params,137,true);
	testServer.setWatchdogInterval(std::chrono::seconds(2));
	testServer.setMaxWorkerSilenceTime(std::chrono::seconds(4));
	testServer.run();

	return 0;
}
