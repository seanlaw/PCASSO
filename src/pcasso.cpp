//Sean M. Law

#include "Molecule.hpp"
#include "Analyze.hpp"
#include "Trajectory.hpp"

#include <iostream>
#include <ctime>
#include <fstream>
#include <cstdlib>
#include <limits>
#include <sstream>

void usage(){
  std::cerr << "Usage:   pcasso [-options] <PDBfile>" << std::endl;
  std::cerr << "Options: [-predict | -features]" << std::endl;
	std::cerr << "         [-trj TRAJfile]" << std::endl;
	std::cerr << "         [-skip frames] [-start frame] [-stop frame]" << std::endl;
//	std::cerr << "         [-dssp dsspFile]" << std::endl;
//	std::cerr << "         [-trial]" << std::endl;
  std::cerr << std::endl;
  exit(0);
}

int main (int argc, char **argv){


  int i;
	unsigned int j;
  std::vector<std::string> pdbs;
  std::string currArg;
  std::string dssp;
	PcassoOutEnum out;
	std::vector<std::string> trajs;
	int start=0;
	int stop=std::numeric_limits<int>::max();
	int skip=0;
	bool startFlag=false;
	unsigned int itrj;
	std::ifstream trjin;
	Trajectory *ftrjin;
	unsigned int nframe;
	AnalyzePcasso *anin;

  dssp.clear();
	out=PREDICT;
	trajs.clear();
	ftrjin=NULL;
	nframe=0;
	anin=NULL;

  for (i=1; i<argc; i++){
    currArg=argv[i];
    if (currArg.compare("-h") == 0 || currArg.compare("-help") == 0){
      usage();
    }
    else if (currArg.compare("-dssp") == 0){
      currArg=argv[++i];
      dssp=currArg;
    }
		else if (currArg.compare("-predict") == 0 || currArg.compare("-prediction") == 0){
			out=PREDICT;
		}
		else if (currArg.compare("-features") == 0 || currArg.compare("-feature") == 0){
			out=FEATURES;
		}
		else if (currArg.compare("-trj") == 0 || currArg.compare("-traj") == 0){
			currArg=argv[++i];
			trajs.push_back(currArg);
		}
		else if (currArg.compare("-skip") == 0){
      currArg=argv[++i];
      std::stringstream(currArg) >> skip;
    }
    else if (currArg.compare("-start") == 0){
      currArg=argv[++i];
      std::stringstream(currArg) >> start;
      start--;
      startFlag=true;
    }
    else if (currArg.compare("-stop") == 0){
      currArg=argv[++i];
      std::stringstream(currArg) >> stop;
    }
		else if (currArg.compare(0,1,"-") == 0){
      std::cerr << "Warning: Skipping unknown option \"" << currArg << "\"" << std::endl;
    }
    else{
      pdbs.push_back(currArg);
    }
  }

  if (pdbs.size() == 0){
    std::cerr << std::endl << "Error: Please provide an input file" << std::endl << std::endl;
    usage();
  }

	Molecule *mol=NULL;

	if (trajs.size() > 0){
		if (pdbs.size() > 1){
			std::cerr << std::endl << "Warning: Only the first PDB structure is used for trajectory analysis" << std::endl << std::endl;
		}
		//Trajectory analysis
		anin=new AnalyzePcasso;
		anin->addSel(":.CA");
		anin->setOutType(out);
		mol=Molecule::readPDB(pdbs.at(0));
		anin->preAnalysis(mol, "");
		mol->selAll();
		//Process trajectories
		for (itrj=0; itrj< trajs.size(); itrj++){
			trjin.open(trajs.at(itrj).c_str(), std::ios::binary);
			if (trjin.is_open()){
				ftrjin=new Trajectory;
      	ftrjin->setMolecule(mol);
				if (ftrjin->findFormat(trjin) == true){
					ftrjin->readHeader(trjin);
        	if (skip > 0 && startFlag == false){
          	start=skip;
        	}
					//Loop through desired frames
        	for (i=start; i< ftrjin->getNFrame() && i< stop; i=i+1+skip){
          	ftrjin->readFrame(trjin, i);
          	nframe++;
						//Analyze PCASSO
						anin->runAnalysis();					
					}
				}
				else{
					std::cerr << "Warning: Skipping unknown trajectory format \"";
        	std::cerr << trajs.at(itrj) << "\"" << std::endl;
				}
				if (ftrjin != NULL){
					delete ftrjin;
				}
			}
			trjin.close();
		}
	}
  else if (dssp.length() > 0){
    if (pdbs.size() > 1){
      std::cerr << "Warning: \"-dssp\" option can only be used with a single PDB input" << std::endl;
      std::cerr << "Warning: Only the first PDB structure was processed" << std::endl;
    }
    //Only process first structure!
    mol=Molecule::readPDB(pdbs.at(0));
    std::cerr << "Processing file \"" << pdbs.at(0) << "..." << std::endl;
		mol->pcasso(dssp, out); //Makes temporary clone with C-alpha only, and analyzes it
    delete mol;
  }
  else{
		//Placed here for efficiency; construct trees once instead of calling mol->pcasso()
		anin=new AnalyzePcasso;
		anin->addSel(":.CA");
		anin->setOutType(out);
		
	  for (j=0; j< pdbs.size(); j++){
  	  mol=Molecule::readPDB(pdbs.at(j));
		  std::cerr << "Processing file \"" << pdbs.at(j) << "..." << std::endl;
		
			/*
			std::clock_t start;
			double duration;
			start=std::clock();
			*/

			//mol->pcasso("", out); //Removed for efficiency; avoid re-constructing trees

			anin->clearMol();
			anin->preAnalysis(mol, "");

			anin->runAnalysis();

			/*
			duration = (std::clock() - start) / (double) CLOCKS_PER_SEC;
			std::cerr << "* " << mol->getNAtom() << " " << duration << std::endl;
			*/

		  delete mol;
    }
	}
	if (anin != NULL){
		delete anin;
	}

  return 0;
}
