//Sean M. Law

#include "Analyze.hpp"

#include "Molecule.hpp"
#include "Chain.hpp"
#include "Residue.hpp"
#include "Atom.hpp"
#include "Vector.hpp"
#include "Constants.hpp"
#include "DTree.hpp"
#include "PCASSO.hpp"
#include "Misc.hpp"

#include <fstream>
#include <iomanip>

Analyze::Analyze (){
	sel.clear();
	mol.clear();
	nframe=0;
	verbose=false;
}

Analyze::~Analyze (){
	//Do nothing
}

AnalyzePcasso::AnalyzePcasso (std::string delim){
	std::vector<std::string> tokens;

	pout=PREDICT;
	t.clear();
	t.resize(PCASSO::getNTree());
	for (unsigned int i=0; i< PCASSO::getNTree(); i++){
		t.at(i)=new DTree;
		Misc::splitStr(Misc::trim(PCASSO::getTree(i)), " \t", tokens, false);
		t.at(i)->genDTree(tokens, delim);
	}
}

void Analyze::addSel(const std::string& selin){
	this->sel.push_back(selin);
}

std::string Analyze::getSel(const int& element){
	return this->sel.at(element);
}

unsigned int Analyze::getNSel(){
	return this->sel.size();
}

void Analyze::addMol(Molecule* molin){
	this->mol.push_back(molin);
}

void Analyze::setMol(const int& element, Molecule* molin){
	this->mol.at(element)=molin;
}

void Analyze::clearMol(){
	this->mol.clear();
}

void Analyze::resizeNMol(const int sizein){
  this->mol.resize(sizein);
}

void Analyze::readTopology(Molecule* molin, std::string topin){
  if (topin.length() > 0){
    molin->readTopology(topin);
    molin->setMass();
    molin->setCharge();
  }
}

void Analyze::setupMolSel(Molecule* molin){
	Molecule* tmpmol;

  molin->select(this->getSel(0));
  tmpmol=molin->copy();
  this->addMol(tmpmol);
}

Molecule* Analyze::getMol(const int& element){
	return this->mol.at(element);
}

unsigned int Analyze::getNMol(){
	return this->mol.size();
}

std::vector<std::vector<double> >& Analyze::getFDataVec(){
	return fdata;
}

void Analyze::setInput(const std::string& fin){
  ifile=fin;
}

std::string Analyze::getInput(){
  return ifile;
}

void Analyze::setNFrame(const unsigned int nframein){
	nframe=nframein;
}

unsigned int Analyze::getNFrame(){
	return nframe;
}

bool Analyze::getVerbose(){
  return verbose;
}

void Analyze::setVerbose(bool verbosein){
  verbose=verbosein;
}


//All preAnalysis Functions
void Analyze::preAnalysis(){
  //Do nothing
}

void Analyze::preAnalysis(Molecule* molin, std::string topin){
  this->readTopology(molin, topin);
	this->setupMolSel(molin);
}

void AnalyzePcasso::preAnalysis(Molecule* molin, std::string fin){
	this->setupMolSel(molin);
	this->setInput(fin);
	
	//Resize FData in Analyze::pcasso
}


//All runAnalysis Functions

void AnalyzePcasso::runAnalysis(){
	std::ifstream dsspFile;
  std::istream* dsspinp;
  std::string line;
  std::vector<std::string> dssp;
  std::vector<std::string> s;
  unsigned int natom;

	Analyze::pcasso(this->getMol(0), this->getFDataVec()); //PCASSO features get stored in the second argument (a 2-D vector double)

	natom=0;

  //Read DSSP file first
  if (this->getInput().length() > 0){
    dsspFile.open(this->getInput().c_str(), std::ios::in);
    dsspinp=&dsspFile;
    while (dsspinp->good() && !(dsspinp->eof())){
      getline(*dsspinp, line);
      Misc::splitStr(line, " \t", s, false);
      if (s.size() > 0){
        if (s.at(s.size()-2).compare(0,1,"E") == 0 || s.at(s.size()-2).compare(0,1,"B") == 0){
          dssp.push_back("E");
        }
        else if (s.at(s.size()-2).compare(0,1,"H") == 0 || s.at(s.size()-2).compare(0,1,"I") == 0 || s.at(s.size()-2).compare(0,1,"G") == 0){
          dssp.push_back("H");
        }
        else if (s.at(s.size()-2).compare(0,1,"-") == 0 || s.at(s.size()-2).compare(0,1,"S") == 0 || s.at(s.size()-2).compare(0,1,"T") == 0){
          dssp.push_back("C");
        }
        else{
          std::cerr << "Warning unrecognized DSSP classification \"" << s.at(s.size()-2) << "\" has been set to \"X\"" << std::endl;
          dssp.push_back("X");
        }
      }
    }
    if (dsspFile.is_open()){
      dsspFile.close();
    }
  }

  std::vector<std::vector<double> > &feat=this->getFDataVec();
	std::map<std::string, unsigned int> vote; 
	std::map<std::string, unsigned int>::iterator iter;
	std::string tmpClass;
	std::string maxClass;
	unsigned int maxVote;
	bool majority;
	unsigned int ntree;

	ntree=PCASSO::getNTree();

	if (this->getOutType() == PREDICT || this->getOutType() == PREDICTION){
		for (unsigned int i=0; i< feat.size(); i++){
			vote.clear();
			maxVote=0;
			majority=false;
			for (unsigned int j=0; j< ntree && majority == false; j++){
				tmpClass=t.at(j)->getDTreeClass(feat.at(i));
				if (vote.find(tmpClass) != vote.end()){
					vote.at(tmpClass)++;
					if (vote.at(tmpClass) > maxVote){
						//Find majority vote, method adapted from openCV
						//Is pseudo-random since it depends on the order of the trees
						maxVote=vote.at(tmpClass);
						maxClass=tmpClass;
						if (static_cast<float>(maxVote)/ntree > 0.5){ //Unsigned integer division!
							majority=true;
						}
					}
				}
				else{
					vote.insert(std::pair<std::string, unsigned int>(tmpClass,1));
				}
			}
			if (this->getVerbose() == true && this->getMol(0)->getNAtom() == feat.size()){
				std::cout << this->getMol(0)->getAtom(i)->getSummary() << " ";
				std::cout << maxClass; //Print majority vote
				if (this->getNFrame() > 0){
					std::cout << " " << this->getNFrame();
				}
				std::cout << std::endl;
			}
			else{
				std::cout << maxClass << std::endl; //Print majority vote
			}
		}
	}
	else if (this->getOutType() == FEATURES || this->getOutType() == FEATURE){
    //Print features
    for (unsigned int i=0; i < feat.size(); i++){
      if (this->getInput().length() > 0 && natom < dssp.size()){
        std::cout << dssp.at(natom) << ",";
      }
      else{
        std::cout << "X" << ",";
      }

      for (unsigned int j=0; j< feat.at(i).size(); j++){
        std::cout << feat.at(i).at(j);
        if (j < feat.at(i).size()-1){
          std::cout << ",";
        }
      }

      std::cout << std::endl;
      natom++;
    } //Loop through chains

    if (this->getInput().length() > 0 && dssp.size() != natom){
      std::cerr << "Warning: DSSP (" << dssp.size() << ") and NATOM (" << natom << ") mismatch" << std::endl;
    }
	}
	else{
		std::cerr << "Warning: Unrecognized PCASSO output type" << std::endl;
	}

}

//All postAnalysis functions

void Analyze::postAnalysis(){
	//Do nothing
}

void AnalyzePcasso::postAnalysis(){
	while (!t.empty()){
		t.back()->delDTree();
		t.pop_back();
	}
}

//Basic analysis functions

Vector Analyze::centerOfGeometry(Molecule *mol, bool selFlag){
  Vector cog=Vector(0.0, 0.0, 0.0);

  for (unsigned int i=0; i< mol->getAtmVecSize(); i++){
    if (selFlag == true && mol->getAtom(i)->getSel() == false){
      continue;
    }
    cog+=mol->getAtom(i)->getCoor();
  }

  if (selFlag == true){
    cog/=mol->getNAtomSelected();
  }
  else{
    cog/=mol->getNAtom();
  }

  return cog;
}

double Analyze::distance (const Vector& u, const Vector& v){
	Vector d=u-v;
	return d.norm();
}

double Analyze::angle (const Vector& u, const Vector& v, const Vector& w){
  double angle;
  Vector dx, dy;
  double dp, nx, ny;

  dx=u-v;
  dy=w-v;

  nx=dx.norm();
  ny=dy.norm();

  dp=dx.dot(dy);

  angle=acos(dp/(nx*ny));

  return angle/PI*180.0;
}

double Analyze::dihedral (const Vector& t, const Vector& u, const Vector& v, const Vector& w) {
  double dihedral;
  Vector dx, dy, dz, p1, p2, p3;
  double np1, np2, dp1, dp2, ts;

  dx=t-u;
  dy=u-v;
  dz=w-v; //This is correct!

  p1=dx.cross(dy);

  np1=p1.norm();
  p1=p1/np1;

  p2=dz.cross(dy);
  np2=p2.norm();
  p2=p2/np2;

  dp1=p1.dot(p2); //Dot product

  ts=1.0-dp1*dp1;
  ts=(ts<0.0)?0.0:sqrt(ts);
  dihedral=PI/2.0-atan2(dp1,ts);

  p3=p1.cross(p2);

  dp2=p3.dot(dy); //Dot product

  if (dp2 > 0.0){
    dihedral=-dihedral;
  }

  return dihedral/PI*180.0;
}

double Analyze::distance (Molecule* sel1, Molecule* sel2, bool selFlag){
	return Analyze::distance(Analyze::centerOfGeometry(sel1,selFlag), Analyze::centerOfGeometry(sel2,selFlag));
}

double Analyze::angle (Molecule* sel1, Molecule* sel2, Molecule* sel3, bool selFlag){
	return Analyze::angle(Analyze::centerOfGeometry(sel1,selFlag), Analyze::centerOfGeometry(sel2,selFlag), Analyze::centerOfGeometry(sel3,selFlag));
}

double Analyze::dihedral (Molecule* sel1, Molecule* sel2, Molecule* sel3, Molecule* sel4, bool selFlag){
	return Analyze::dihedral(Analyze::centerOfGeometry(sel1,selFlag), Analyze::centerOfGeometry(sel2,selFlag), Analyze::centerOfGeometry(sel3,selFlag), Analyze::centerOfGeometry(sel4,selFlag));
}

void Analyze::pairwiseDistance(Molecule *mol, std::vector<std::vector<double> >& pdin){
	std::vector<Atom*>::iterator ai;
	std::vector<Atom*>::iterator aj;
	unsigned int natom;
	unsigned int aiInx;
	bool flag;

	natom=mol->getAtmVecSize();

	pdin.clear();
	pdin.resize(natom);

	for (ai=mol->getAtmVec().begin(); ai != mol->getAtmVec().end(); ++ai){
		aiInx=(*ai)->getAtmInx();
		pdin.at(aiInx).resize(natom);
		pdin.at(aiInx).at(aiInx)=0.0; //Zero diagonal
	}

	for (ai=mol->getAtmVec().begin(); ai != mol->getAtmVec().end(); ++ai){
		aiInx=(*ai)->getAtmInx();
		if ((*ai)->getX() < 9999.9){
			flag=true;
		}
		else{
			flag=false;
		}
		//Lower Triangle
		for (aj=mol->getAtmVec().begin(); aj != ai; ++aj){
      if (flag && (*aj)->getX() < 9999.9){
        pdin.at(aiInx).at((*aj)->getAtmInx())=Analyze::distance((*ai)->getCoor(), (*aj)->getCoor());
      }
      else{
        pdin.at(aiInx).at((*aj)->getAtmInx())=9999.9;
      }
    }
	
		//Upper Triangle
		for (aj=ai+1; aj != mol->getAtmVec().end(); ++aj){
			if (flag && (*aj)->getX() < 9999.9){
        pdin.at(aiInx).at((*aj)->getAtmInx())=Analyze::distance((*ai)->getCoor(), (*aj)->getCoor());
      }
      else{
        pdin.at(aiInx).at((*aj)->getAtmInx())=9999.9;
      }
		}
	}

}

void Analyze::allAnglesDihedrals(Molecule *mol, std::vector<std::vector<double> >& anglesin){
	unsigned int i, j;
	Chain *c;
	Atom *atmI, *iMinusTwo, *iMinusOne, *iPlusOne, *iPlusTwo, *iPlusThree;
	unsigned int size;
	std::vector<double> angles;

	anglesin.resize(mol->getAtmVecSize());

	for (i=0; i< mol->getChnVecSize(); i++){
    c=mol->getChain(i);
    size=c->getAtmVecSize();
    for (j=0; j< c->getAtmVecSize(); j++){
      atmI=c->getAtom(j);
			iMinusTwo=NULL;
			iMinusOne=NULL;
      iPlusOne=NULL;
      iPlusTwo=NULL;
      iPlusThree=NULL;

			if (j > 1 && atmI->getResId()-2 == c->getAtom(j-2)->getResId()){
        iMinusTwo=c->getAtom(j-2);
      }
      else{
        if (j > 1 && atmI->getResId() == c->getAtom(j-2)->getResId() && (atmI->getICode().compare(0,1,c->getAtom(j-1)->getICode(),0,1) != 0)){
          iMinusTwo=c->getAtom(j-2);
        }
      }

			if (j > 0 && atmI->getResId()-1 == c->getAtom(j-1)->getResId()){
        iMinusOne=c->getAtom(j-1);
      }
      else{
        if (j > 0 && atmI->getResId() == c->getAtom(j-1)->getResId() && (atmI->getICode().compare(0,1,c->getAtom(j-1)->getICode(),0,1) != 0)){
          iMinusOne=c->getAtom(j-1);
        }
      }

			if (j+1 < size && atmI->getResId()+1 == c->getAtom(j+1)->getResId()){
				iPlusOne=c->getAtom(j+1);
			}
			else{
				if (j+1 < size && atmI->getResId() == c->getAtom(j+1)->getResId() && (atmI->getICode().compare(0,1,c->getAtom(j+1)->getICode(),0,1) != 0)){
					iPlusOne=c->getAtom(j+1);
				}
			}

			if (j+2 < size && atmI->getResId()+2 == c->getAtom(j+2)->getResId()){
				iPlusTwo=c->getAtom(j+2);
			}
		 	else{
        if (j+2 < size && atmI->getResId() == c->getAtom(j+2)->getResId() && (atmI->getICode().compare(0,1,c->getAtom(j+2)->getICode(),0,1) != 0)){
          iPlusTwo=c->getAtom(j+2);
        }
      }

			if (j+3 < size && atmI->getResId()+3 == c->getAtom(j+3)->getResId()){
        iPlusThree=c->getAtom(j+3);
      }
			else{
        if (j+3 < size && atmI->getResId() == c->getAtom(j+3)->getResId() && (atmI->getICode().compare(0,1,c->getAtom(j+3)->getICode(),0,1) !=0)){
          iPlusThree=c->getAtom(j+3);
        }
      }
			
			angles.clear();
			angles.resize(3, 9999.9);
			if (iMinusOne != NULL && iPlusOne != NULL){
				//Get Angle
				angles.at(0)=Analyze::angle(iMinusOne->getCoor(), atmI->getCoor(), iPlusOne->getCoor());
			}
			if (iPlusOne != NULL && iPlusTwo != NULL && iPlusThree != NULL){
				//Get Dihedral
				angles.at(1)=Analyze::dihedral(atmI->getCoor(), iPlusOne->getCoor(), iPlusTwo->getCoor(), iPlusThree->getCoor());
			}
			if (iMinusTwo != NULL && iPlusTwo != NULL){
				//Get Wide Angle (i-2, i, i+2)
				angles.at(2)=Analyze::angle(iMinusTwo->getCoor(), atmI->getCoor(), iPlusTwo->getCoor());
			}
			anglesin.at(atmI->getAtmInx()).resize(3);
			anglesin.at(atmI->getAtmInx())=angles;
    }
  }
}

void Analyze::pcasso(Molecule* mol, std::vector<std::vector<double> > &fdataIO){
	Chain *c;
	Atom *ai, *aj, *ak;
	unsigned int i, j, start;
	double defVal;
	double iMinus6;//For non-local contacts
	double iPlus6; //For non-local contacts
	unsigned minx;
	unsigned pinx;
	int diffResId;
	std::vector<std::vector<double> > caPairDist; //Ca-Ca Distances
	std::vector<std::vector<double> > pcPairDist; //Pc-Pc Distances
	std::vector<std::vector<double> > caAngles; //Ca-Ca Angle/Diehdral/Wide
	std::vector<std::vector<double> > pcAngles; //Pc-Pc Angle/Dihedral/Wide
	Molecule *camol, *pcmol;
	unsigned int natom;
	double dist;

	defVal=9999.9;
	camol=NULL;
	pcmol=NULL;

	mol->storeSel();
	mol->select(":.CA");
	camol=mol->clone(true,true); //Copy selection, keep original
	pcmol=mol->clone(true,true);
	mol->recallSel(); //Restore original selection
	mol->eraseSel();

	camol->assignAtmInx();
	pcmol->assignAtmInx();

	//Analyze all C-alpha first
	Analyze::pairwiseDistance(camol, caPairDist);
	Analyze::allAnglesDihedrals(camol, caAngles);

	natom=camol->getAtmVecSize();

	for (unsigned int ichain=0; ichain < camol->getChnVecSize(); ichain++){
		c=camol->getChain(ichain);
		for (unsigned int iatom=0; iatom < c->getAtmVecSize(); iatom++){
			ai=c->getAtom(iatom);
			ai->clearData();
			i=ai->getAtmInx();
			//i-5, i-4, i-3, i-2, i-1, i+1, i+2, i+3, i+4, i+5 Distances
			//Deal with unsigned int subtraction from zero
			if (iatom == 0){
				ai->addData(defVal);
				ai->addData(defVal);
				ai->addData(defVal);
				ai->addData(defVal);
				ai->addData(defVal);
				start=iatom-0;
			}
			else if (iatom == 1){
				//std::cout << defVal << " ";
				ai->addData(defVal);
        ai->addData(defVal);
        ai->addData(defVal);
        ai->addData(defVal);
				start=iatom-1;
			}
			else if (iatom == 2){
				ai->addData(defVal);
        ai->addData(defVal);
        ai->addData(defVal);
				start=iatom-2;
			}
			else if (iatom == 3){
				ai->addData(defVal);
        ai->addData(defVal);
				start=iatom-3;
			}
			else if (iatom == 4){
				ai->addData(defVal);	
				start=iatom-4;
			}
			else{
				start=iatom-5;
			}
			for (j=start; j<= iatom+5; j++){
				aj=c->getAtom(j);
				if (aj == NULL){
					ai->addData(defVal);
				}
				else if (j == iatom){
					//Distance == 0
					continue;
				}
				else{
					ai->addData(caPairDist.at(i).at(aj->getAtmInx()));
				}
			}

			//Angles and Dihedrals
			for (j=0; j< caAngles.at(i).size(); j++){
				ai->addData(caAngles.at(i).at(j));
			}
		
			//Shortest non-local contact distance, >= i+6 and <= i-6
			iPlus6=1E10;
			iMinus6=1E10;
			pinx=natom;
			minx=natom;

			for (j=0; j< natom; j++){
				aj=camol->getAtom(j);

				if (ai == aj){
					continue;
				}
				
				dist=caPairDist.at(i).at(j);

				//i+6
				//Assess distance first to avoid unnecessary string comparison
				if (dist < iPlus6){
					if (ai->getChainId().compare(aj->getChainId()) != 0){
						//atom i and atom j are on different chains
						pinx=j;
						iPlus6=dist;
					}
					else{
						//atom i and atom j are on the same chain
						diffResId=aj->getResId() - ai->getResId();
						if (diffResId >= 6){
						  pinx=j;
							iPlus6=dist;
						}
						else{
							if (diffResId == 0 && (Misc::atoi(aj->getICode()) - Misc::atoi(ai->getICode()) >= 6)){
           	 		pinx=j;
								iPlus6=dist;
          		}
						}
					}
				}

				//i-6
				//Assess distance first to avoid unnecessary string comparison
				if (dist < iMinus6){
					if (ai->getChainId().compare(aj->getChainId()) != 0){
					//atom i and atom j are on different chains
            minx=j;
						iMinus6=dist;
          }
					else{
						//atom i and atom j are on the same chain
						diffResId=aj->getResId() - ai->getResId();
						if (diffResId <= -6){
            	minx=j;
							iMinus6=dist;
          	}
						else{
							if (diffResId == 0 && (Misc::atoi(aj->getICode()) - Misc::atoi(ai->getICode()) <= -6)){
            		minx=j;
								iMinus6=dist;
          		}
						}
					}
				}
			}

			int k;
			unsigned int q;
			int max;
			max=10;
			if (iatom == 0){
				for (k=0; k< max; k++){
					ai->addData(defVal);
				}
				start=0;
			}
			else{
				start=iatom-1;
			}
			for (q=start; q<= iatom+1; q++){
				if (q < c->getAtmVecSize()){
					ak=c->getAtom(q);
					for (k=static_cast<int>(pinx)-2; k<=static_cast<int>(pinx)+2; k++){
						if (k >= 0 && k < static_cast<int>(natom)){
							aj=camol->getAtom(k);
							ai->addData(caPairDist.at(ak->getAtmInx()).at(aj->getAtmInx()));
						}
						else{
							ai->addData(defVal);
						}
					}
					for (k=static_cast<int>(minx)-2; k<=static_cast<int>(minx)+2; k++){
        		if (k >= 0 && k < static_cast<int>(natom)){
							aj=camol->getAtom(k);
							ai->addData(caPairDist.at(ak->getAtmInx()).at(aj->getAtmInx()));
        		}
        		else{
          		ai->addData(defVal);
        		}	
					}
				}
				else{
					for (k=0; k< max; k++){
						ai->addData(defVal);
					}
				}
			}

		} //Loop through atoms
	}//Loop through chains

	//Analyze all pseudocenter
	pcmol->modPseudoCenter();
	Analyze::pairwiseDistance(pcmol, pcPairDist);
	Analyze::allAnglesDihedrals(pcmol, pcAngles);

	natom=pcmol->getAtmVecSize();

	for (unsigned int ichain=0; ichain < pcmol->getChnVecSize(); ichain++){
		c=pcmol->getChain(ichain);
		for (unsigned int iatom=0; iatom < c->getAtmVecSize(); iatom++){
			ai=camol->getChain(ichain)->getAtom(iatom); //From C-alpha
			ak=c->getAtom(iatom);
			i=ak->getAtmInx();
			//i-5, i-4, i-3, i-2, i-1, i+1, i+2, i+3, i+4, i+5 Distances
			//Deal with unsigned int subtraction from zero
			if (iatom == 0){
				ai->addData(defVal);
				ai->addData(defVal);
				ai->addData(defVal);
				ai->addData(defVal);
				ai->addData(defVal);
				start=iatom-0;
			}
			else if (iatom == 1){
				//std::cout << defVal << " ";
				ai->addData(defVal);
        ai->addData(defVal);
        ai->addData(defVal);
        ai->addData(defVal);
				start=iatom-1;
			}
			else if (iatom == 2){
				ai->addData(defVal);
        ai->addData(defVal);
        ai->addData(defVal);
				start=iatom-2;
			}
			else if (iatom == 3){
				ai->addData(defVal);
        ai->addData(defVal);
				start=iatom-3;
			}
			else if (iatom == 4){
				ai->addData(defVal);	
				start=iatom-4;
			}
			else{
				start=iatom-5;
			}
			for (j=start; j<= iatom+5; j++){
				aj=c->getAtom(j);
				if (aj == NULL){
					ai->addData(defVal);
				}
				else if (j == iatom){
					//Distance == 0
					continue;
				}
				else{
					ai->addData(pcPairDist.at(i).at(aj->getAtmInx()));
				}
			}

			//Angles and Dihedrals
			for (j=0; j< pcAngles.at(ak->getAtmInx()).size(); j++){
				ai->addData(pcAngles.at(ak->getAtmInx()).at(j));
			}
		
			//Shortest non-local contact distance, >= i+6 and <= i-6
			iPlus6=1E10;
			iMinus6=1E10;
			pinx=natom;
			minx=natom;
			for (j=0; j< natom; j++){
				aj=pcmol->getAtom(j);
				if (ak == aj){
					continue;
				}

				dist=pcPairDist.at(i).at(j);

				//i+6
				//Assess distance first to avoid unnecessary string comparison
				if (dist < iPlus6){
					if (ak->getChainId().compare(aj->getChainId()) != 0){
						//atom i and atom j are on different chains
						pinx=j;
						iPlus6=dist;
					}
					else{
						//atom i and atom j are on the same chain
						diffResId=aj->getResId() - ak->getResId();
						if (diffResId >= 6){
						  pinx=j;
							iPlus6=dist;
						}
						else{
							if (diffResId == 0 && (Misc::atoi(aj->getICode()) - Misc::atoi(ak->getICode()) >= 6)){
           	 		pinx=j;
								iPlus6=dist;
          		}
						}
					}
				}

				//i-6
				//Assess distance first to avoid unnecessary string comparison
				if (dist < iMinus6){
					if (ak->getChainId().compare(aj->getChainId()) != 0){
						//atom i and atom j are on different chains
            minx=j;
						iMinus6=dist;
          }
					else{
						//atom i and atom j are on the same chain
						diffResId=aj->getResId() - ak->getResId();
						if (diffResId <= -6){
            	minx=j;
							iMinus6=dist;
          	}
						else{
							if (diffResId == 0 && (Misc::atoi(aj->getICode()) - Misc::atoi(ak->getICode()) <= -6)){
            		minx=j;
								iMinus6=dist;
          		}
						}
					}
				}
			}

			int k;
			unsigned int q;
			int max;
			max=10;
			if (iatom == 0){
				for (k=0; k< max; k++){
					ai->addData(defVal);
				}
				start=0;
			}
			else{
				start=iatom-1;
			}
			for (q=start; q<= iatom+1; q++){
				if (q < c->getAtmVecSize()){
					ak=c->getAtom(q);
					for (k=static_cast<int>(pinx)-2; k<=static_cast<int>(pinx)+2; k++){
						if (k >= 0 && k < static_cast<int>(natom)){
							aj=pcmol->getAtom(k);
							ai->addData(pcPairDist.at(ak->getAtmInx()).at(aj->getAtmInx()));
						}
						else{
							ai->addData(defVal);
						}
					}
					for (k=static_cast<int>(minx)-2; k<=static_cast<int>(minx)+2; k++){
        		if (k >= 0 && k < static_cast<int>(natom)){
							aj=pcmol->getAtom(k);
							ai->addData(pcPairDist.at(ak->getAtmInx()).at(aj->getAtmInx()));
        		}
        		else{
          		ai->addData(defVal);
        		}	
					}
				}
				else{
					for (k=0; k< max; k++){
						ai->addData(defVal);
					}
				}
			}

		} //Loop through atoms
	}//Loop through chains

	//Always clear and resize!
	fdataIO.clear();
	fdataIO.resize(camol->getNAtom());

	//Store features
	natom=0; //This is needed since natom is used for other things above
	for (unsigned int ichain=0; ichain < camol->getChnVecSize(); ichain++){
    c=camol->getChain(ichain);
    for (unsigned int iatom=0; iatom < c->getAtmVecSize(); iatom++){
      ai=c->getAtom(iatom);
			fdataIO.at(natom).reserve(3*ai->getDataSize());
	
			//Store S(i)
      for (j=0; j< ai->getDataSize(); j++){
				fdataIO.at(natom).push_back(ai->getDataPoint(j));
      }

			//Store S(i-1)
			if (iatom > 0){
				//Not first atom of chain
				aj=c->getAtom(iatom-1);
				for (j=0; j< ai->getDataSize(); j++){
					if (aj != NULL){
						fdataIO.at(natom).push_back(aj->getDataPoint(j));
					}
					else{
						fdataIO.at(natom).push_back(defVal);
					}
				}
			}
			else{
				//Store S(i-1) which is all default values
				for (j=0; j< ai->getDataSize(); j++){
					fdataIO.at(natom).push_back(defVal);
				}
			}

			//Store S(i+1)
			aj=c->getAtom(iatom+1);
			for (j=0; j< ai->getDataSize(); j++){
        if (aj != NULL){
					fdataIO.at(natom).push_back(aj->getDataPoint(j));
        }
        else{
					fdataIO.at(natom).push_back(defVal);
        }
      }
			natom++;
		} //Loop through atoms
	} //Loop through chains


	if (camol != NULL){
		delete camol;
	}
	if (pcmol != NULL){
		delete pcmol;
	}
}

void AnalyzePcasso::setOutType(PcassoOutEnum pin){
	pout=pin;
}

PcassoOutEnum AnalyzePcasso::getOutType(){
 return pout;
}

