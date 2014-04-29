//Sean M. Law
//Aaron T. Frank

/*
This file is part of PCASSO.

PCASSO is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

PCASSO is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with PCASSO.  If not, see <http://www.gnu.org/licenses/>.

This file is part of MoleTools.

MoleTools is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

MoleTools is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with MoleTools.  If not, see <http://www.gnu.org/licenses/>.
*/

#include "Molecule.hpp"

#include "Chain.hpp"
#include "Residue.hpp"
#include "Atom.hpp"
#include "PDB.hpp"
#include "Select.hpp"
#include "Analyze.hpp"
#include "Misc.hpp"

#include <limits>

Molecule::Molecule (){
  chnVec.clear();
  resVec.clear();
  atmVec.clear();
  copyFlag=false;
  storedSel.clear();
  remarks.clear();
  iCodeFlag=false;
  year=0;
  exp.clear();
}

Molecule::~Molecule (){
  Chain *c;
  Residue *r;
  Atom *a;

  for (unsigned int i=0; i< this->getChnVecSize(); i++){
    c=this->getChain(i);
    for (unsigned int j=0; j< c->getResVecSize(); j++){
      r=c->getResidue(j);
      if (this->getCopyFlag() == false){
        //Only destruct atoms if molecule is NOT a copy
        for (unsigned int k=0; k< r->getAtmVecSize(); k++){
          a=r->getAtom(k);
          delete a;
        }
      }
      delete r;
    }
    delete c;
  }
}

Molecule* Molecule::readPDB (const std::string ifile, const int model, const std::string format, const bool hetFlag){
  if (ifile.length() == 0){
    std::cerr << "Error: PDB file \"" << ifile << "\" cannot be found" << std::endl;
    return new Molecule;
  }
  else{
    Molecule* mol=PDB::readPDB(ifile, model, format, hetFlag);
    mol->format();
    return mol;
  }
}

Molecule* Molecule::readPDB (const std::string ifile, const std::string format, const bool hetFlag){
  //Overloaded function
  return Molecule::readPDB(ifile, 0, format, hetFlag);
}

std::string Molecule::writePDB(bool selFlag, bool print, bool chnFlag){

  std::ostringstream out;

  PDB::writePDBFormat(this, out, selFlag, chnFlag);

  if (print == true){
    std::cout << out.str();
  }
  
  return out.str();
}

std::string Molecule::writePDB(bool selFlag, bool print){
  return this->writePDB(selFlag, print, false);
}

std::string Molecule::writePDB(bool chnFlag){
  return this->writePDB(true, true, chnFlag);
}

std::string Molecule::writePDB(){
  return this->writePDB(true, true, false);
}

void Molecule::addAtom(Atom* atmEntry) {
  if (atmEntry->getAtmNum()){
    atmVec.push_back(atmEntry);
  }
}

Molecule* Molecule::clone (bool selFlag, bool keep){
  //Deep copy
  Molecule *mol=new Molecule;
  Chain *c, *chnEntry;
  Residue *r, *resEntry;
  Atom *a, *atmEntry;

  c=NULL;
  r=NULL;
  a=NULL;
  chnEntry=NULL;
  resEntry=NULL;
  atmEntry=NULL;
  mol->setCopyFlag(false); //Not a copy
  
  for (unsigned int i=0; i< this->getChnVecSize(); i++){
    c=this->getChain(i);
    chnEntry=new Chain; //Deleted later if no residues/atoms

    for (unsigned int j=0; j< c->getResVecSize(); j++){
      r=c->getResidue(j);
      resEntry=new Residue; //Deleted later if no atoms

      for (unsigned int k=0; k< r->getAtmVecSize(); k++){
        a=r->getAtom(k);
        if(selFlag == true && a->getSel() == false){
          continue;
        }
        //Add each selected atom
        atmEntry=new Atom; //This is necessary!
        atmEntry->clone(a); //Clone Atom
        mol->addAtom(atmEntry);
        resEntry->addAtom(atmEntry);
        chnEntry->addAtom(atmEntry);
      }
      if (resEntry->getAtmVecSize() > 0){
        //Add each residue with selected atoms
        mol->addResidue(resEntry);
        chnEntry->addResidue(resEntry);
      }
      else{
        delete resEntry;
      }
    }
    if (chnEntry->getResVecSize() > 0){
      //Add each chain with selected atoms
      mol->addChain(chnEntry);
    }
    else{
      delete chnEntry;
    }
  }

  if (keep == false){
    delete this;
  }

  mol->setICodeFlag(mol->checkICode());

  return mol;
}

Molecule* Molecule::copy (bool selFlag){
  //Shallow copy of selected atom pointers.
  //Chains and Residues are required to be created new.
  //Useful when you need to make multiple selections
  Molecule *mol=new Molecule;
  Chain *c, *chnEntry;
  Residue *r, *resEntry;
  Atom *a;

  c=NULL;
  r=NULL;
  a=NULL;
  chnEntry=NULL;
  resEntry=NULL;
  mol->setCopyFlag(true); //Is a copy, do not destruct atoms!

  for (unsigned int i=0; i< this->getChnVecSize(); i++){
    c=this->getChain(i);
    chnEntry=new Chain;

    for (unsigned int j=0; j< c->getResVecSize(); j++){
      r=c->getResidue(j);
      resEntry=new Residue;

      for (unsigned int k=0; k< r->getAtmVecSize(); k++){
        a=r->getAtom(k);
        if(selFlag == true && a->getSel() == false){
          continue;
        }
        //Add each selected atom
        mol->addAtom(a); //Copy atom pointer, not clone
        resEntry->addAtom(a);
        chnEntry->addAtom(a);
      }
      if (resEntry->getAtmVecSize() > 0){
        //Add each residue with selected atoms
        mol->addResidue(resEntry);
        chnEntry->addResidue(resEntry);
      }
      else{
        delete resEntry;
      }
    }
    if (chnEntry->getResVecSize() > 0){
      //Add each chain with selected atoms
      mol->addChain(chnEntry);
    }
    else{
      delete chnEntry;
    }     
  } 

  mol->setICodeFlag(mol->checkICode());

  return mol;
}

void Molecule::cat (Molecule* catmol, bool selFlag, bool keep){
  //Deep copy
  Chain *c, *chnEntry;
  Residue *r, *resEntry;
  Atom *a, *atmEntry;

  c=NULL;
  r=NULL;
  a=NULL;
  chnEntry=NULL;
  resEntry=NULL;
  atmEntry=NULL;
  this->setCopyFlag(false); //Not a copy

  for (unsigned int i=0; i< catmol->getChnVecSize(); i++){
    c=catmol->getChain(i);
    chnEntry=new Chain; //Deleted later if no residues/atoms

    for (unsigned int j=0; j< c->getResVecSize(); j++){
      r=c->getResidue(j);
      resEntry=new Residue; //Deleted later if no atoms

      for (unsigned int k=0; k< r->getAtmVecSize(); k++){
        a=r->getAtom(k);
        if(selFlag == true && a->getSel() == false){
          continue;
        }
        //Add each selected atom
        atmEntry=new Atom; //This is necessary!
        atmEntry->clone(a); //Clone Atom
        this->addAtom(atmEntry);
        resEntry->addAtom(atmEntry);
        chnEntry->addAtom(atmEntry);
      }
      if (resEntry->getAtmVecSize() > 0){
        //Add each residue with selected atoms
        this->addResidue(resEntry);
        chnEntry->addResidue(resEntry);
      }
      else{
        delete resEntry;
      }
    }
    if (chnEntry->getResVecSize() > 0){
      //Add each chain with selected atoms
      this->addChain(chnEntry);
    }
    else{
      delete chnEntry;
    }
  }

  if (keep == false){
    delete catmol;
  }

}

Atom* Molecule::getAtom(const unsigned int& element){
  if (element >= atmVec.size()){
    return NULL;
  }
  else{
    return atmVec.at(element);
  }
}

unsigned int Molecule::getAtmVecSize(){
  return atmVec.size();
}

std::vector<Atom*>& Molecule::getAtmVec(){
  return atmVec;
}

std::vector<Atom*> Molecule::getAtmVecClone(){
  return atmVec;
}

void Molecule::addChain(Chain* chnEntry){
  if(chnEntry->getChainId().size()){ //Check if empty string
    chnVec.push_back(chnEntry);
  }
}

void Molecule::addResidue(Residue* resEntry){
  if (resEntry->getResId()){
    resVec.push_back(resEntry);
  }
}

Chain* Molecule::getChain(const unsigned int& element){
  if (element >= chnVec.size()){
    return 0;
  }
  else{
    return chnVec.at(element);
  }
}

unsigned int Molecule::getChnVecSize(){
  return chnVec.size();
}

unsigned int Molecule::getResVecSize(){
  return resVec.size();
}

Residue* Molecule::getResidue(const unsigned int& element){
  if (element >= resVec.size()){
    return NULL;
  }
  else{
    return resVec.at(element);
  }
}

void Molecule::readTopology(const std::string& topin){
  this->toppar.readTopology(topin);
  this->setMass();
  this->setCharge();
}

void Molecule::readParameter(const std::string& prmin){
  std::cerr << "Warning: Molecule::readParameter() has not been implemented yet!" << std::endl;
}

void Molecule::setMass(){
  Atom* a;
  for (unsigned int i=0; i< this->getAtmVecSize(); i++){
    a=this->getAtom(i);
    a->setMass(this->toppar.getMass(Misc::trim(a->getResName()), Misc::trim(a->getAtmName())));
  }
}

void Molecule::setCharge(){
  Atom* a;
  for (unsigned int i=0; i< this->getAtmVecSize(); i++){
    a=this->getAtom(i);
    a->setCharge(this->toppar.getCharge(Misc::trim(a->getResName()), Misc::trim(a->getAtmName())));
  }
}

void Molecule::selAll(){
  for (unsigned int i=0; i< this->getNAtom(); i++){
    this->getAtom(i)->setSel(true);
  }
}

void Molecule::deselAll(){
  for (unsigned int i=0; i< this->getNAtom(); i++){
    this->getAtom(i)->setSel(false);
  }
}

void Molecule::select(std::string sel, bool dieFlag, bool verbose){
  Select::makeSel(this, sel, dieFlag, verbose);
}

unsigned int Molecule::getNAtom(){
  return atmVec.size();
}

unsigned int Molecule::getNAtomSelected(){
  unsigned int natom=0;
  for (unsigned int i=0; i< this->getAtmVecSize(); i++){
    if (this->getAtom(i)->getSel() == true){
      natom++;
    }
  }
  return natom;
}

void Molecule::setCopyFlag(bool copyFlagIn){
  copyFlag=copyFlagIn;
}

bool Molecule::getCopyFlag(){
  return copyFlag;
}

void Molecule::storeSel(std::string key){
  std::vector<bool> selFlags;
  Chain *chn;
  Residue *res;
  Atom *atm;

  for (unsigned int i=0; i< this->getChnVecSize(); i++){
    chn=this->getChain(i);
    for (unsigned int j=0; j< chn->getResVecSize(); j++){
      res=chn->getResidue(j);
      for (unsigned int k=0; k< res->getAtmVecSize(); k++){
        atm=res->getAtom(k);
        selFlags.push_back(atm->getSel());
      }
    }
  }

  storedSel[key]=selFlags;
}

void Molecule::recallSel(std::string key){
  Chain *chn;
  Residue *res;
  Atom *atm;
  unsigned int n;
  std::vector<bool>* storedSelVal;

  n=0;
  storedSelVal=&storedSel[key];
  
  for (unsigned int i=0; i< this->getChnVecSize(); i++){
    chn=this->getChain(i);
    for (unsigned int j=0; j< chn->getResVecSize(); j++){
      res=chn->getResidue(j);
      for (unsigned int k=0; k< res->getAtmVecSize(); k++){
        atm=res->getAtom(k);
        atm->setSel(storedSelVal->at(n));
        n++;
      }
    }
  }
}


void Molecule::eraseSel(std::string key){
  storedSel.erase(key);
}

void Molecule::zeroCoor(){
  Chain *chn;
  Residue *res;
  Atom *atm;
  Coor coor;

  coor=Coor(0.0, 0.0, 0.0);

  for (unsigned int i=0; i< this->getChnVecSize(); i++){
    chn=this->getChain(i);
    for (unsigned int j=0; j< chn->getResVecSize(); j++){
      res=chn->getResidue(j);
      for (unsigned int k=0; k< res->getAtmVecSize(); k++){
        atm=res->getAtom(k);
        atm->setCoor(coor);
      }
    }
  }
  
}

void Molecule::modPseudoCenter(){
  Atom *atmEntry;
  Atom *lastAtom;

  for (unsigned int i=0; i< this->getChnVecSize(); i++){
    atmEntry=NULL;
    lastAtom=NULL;
    for (unsigned int j=0; j< this->getChain(i)->getAtmVecSize(); j++){
      atmEntry=this->getChain(i)->getAtom(j);

      if(lastAtom !=NULL){
        if (lastAtom->getResId()+1 == atmEntry->getResId()){
          //Compute pseudocenter for last atom
          lastAtom->setCoor((lastAtom->getCoor()+atmEntry->getCoor())/2.0);
        }
        else if (lastAtom->getResId() == atmEntry->getResId() && (lastAtom->getICode().compare(0,1,atmEntry->getICode(),0,1) != 0)){
          //Compute pseudocenter for last atom
          lastAtom->setCoor((lastAtom->getCoor()+atmEntry->getCoor())/2.0);
        }
        else{
          //No i+1 neighbor, modify coordinates
          lastAtom->setCoor(Coor(9999.9, 9999.9, 9999.9));
        }
      }
      lastAtom=atmEntry;
    }
    //No i+1 neighbor for last atom in chain, modify coordinates
    lastAtom->setCoor(Coor(9999.9, 9999.9, 9999.9));
  }
}

void Molecule::pcasso (std::string dsspin, PcassoOutEnum out){
  AnalyzePcasso* anin=new AnalyzePcasso;
  
  anin->addSel(":.CA");
  anin->setOutType(out);

  anin->preAnalysis(this, dsspin);

  anin->runAnalysis();

  delete anin;
}

//Virtual Functions

void Molecule::format(){
  //Do nothing
}

void MoleculeCHARMM::format(){
  Chain *chn;
  Residue *res;
  Atom *atm;
  Residue *lastRes;
  Residue *nextRes;
  int nUNK; //Number of unknown residues
  
  lastRes=NULL;
  nextRes=NULL;
  nUNK=0;

  for (unsigned int i=0; i< this->getChnVecSize(); i++){
    chn=this->getChain(i);
    for (unsigned int j=0; j< chn->getResVecSize(); j++){
      if (j>1){
        lastRes=chn->getResidue(j-1);
      }
      else{
        lastRes=NULL;
      }
      if (j< chn->getResVecSize()-1){
        nextRes=chn->getResidue(j+1);
      }
      else{
        nextRes=NULL;
      }
      res=chn->getResidue(j);
      if (res->getResName().compare(0,3,"UNK") == 0){
        nUNK++;
      }
      for (unsigned int k=0; k< res->getAtmVecSize(); k++){
        atm=res->getAtom(k);
        //if (selFlag == true && atm->getSel() == false){
        //  continue;
        //}
        //Perform Formatting
        if (atm->getResName().compare("HIE") == 0){
          atm->setResName("HSE");
        }
        else if (atm->getResName().compare("HID") == 0){
          atm->setResName("HSD");
        }
        else if (atm->getResName().compare("HIP") == 0){
          atm->setResName("HSP");
        }
        else if (atm->getResName().compare("CYX") == 0){
          std::cerr << "Warning: " << atm->getSummary() << " has a disulfide bond" << std::endl;
          this->addRemark("Warning: ");
          atm->setResName("CYS");
        }
        else if (atm->getResName().compare("AHE") == 0 && lastRes != NULL){
          atm->setResName(lastRes->getResName());
          atm->setResId(lastRes->getResId());
        }
        else if (atm->getResName().compare("NME") == 0 && lastRes != NULL){
          atm->setResName(lastRes->getResName());
          atm->setResId(lastRes->getResId());
        }
        else if (atm->getResName().compare("ACE") == 0 && nextRes != NULL){
          atm->setResName(nextRes->getResName());
          atm->setResId(nextRes->getResId());
        }
        else if (atm->getResName().compare("FOR") == 0 || atm->getResName().compare("CSO") == 0 || atm->getResName().compare("CME") == 0){
          std::cerr << "Warning: " << atm->getSummary();
          std::cerr << " has no matching residue name in CHARMM" << std::endl;
          this->addRemark("Warning: ");
          //Do nothing
        }
        else{
          //Do nothing
        }
      } 
    }
  }
  if (nUNK > 1){
    std::cerr << "Warning: More than one (" << nUNK << ") UNK residues found" << std::endl;
  }
}

void Molecule::addRemark(const std::string& remin){
  this->remarks+=remin;
}

void Molecule::clearRemark(){
  this->remarks.clear();
}

std::string Molecule::getRemark(){
  return remarks;
}

bool Molecule::checkICode(){
  for (unsigned int i=0; i< this->getAtmVecSize(); i++){
    if (this->getAtom(i)->getICode().length() != 0 && this->getAtom(i)->getICode().compare(0,1," ") != 0){
      return true;
    }
  }
  return false;
}

void Molecule::setICodeFlag(bool iCodeFlagIn){
  iCodeFlag=iCodeFlagIn;
}

bool Molecule::getICodeFlag(){
  return iCodeFlag;
}

void Molecule::assignAtmInx(){
  unsigned int natom;

  natom=0;

  for (unsigned int i=0; i< this->getAtmVecSize(); i++){
    this->getAtom(i)->setAtmInx(natom);
    natom++;
  }
}

void Molecule::resetAtmInx(){
  for (unsigned int i=0; i< this->getAtmVecSize(); i++){
    this->getAtom(i)->setAtmInx(std::numeric_limits<unsigned int>::max());
  }
}

void Molecule::setYear(const unsigned int& yearin){
  year=yearin;
}

unsigned int Molecule::getYear(){
  return year;
}

void Molecule::setExp(const std::string& expin){
  exp=expin;
}

std::string Molecule::getExp(){
  return exp;
}

