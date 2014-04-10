//Sean M. Law

#include "Chain.hpp"

#include "Residue.hpp"
#include "Atom.hpp"

Chain::Chain (){
  resVec.clear();
  atmVec.clear();
//  sel=true;
}

void Chain::reset(){
  resVec.clear();
  atmVec.clear();
//  sel=true;
}

void Chain::addResidue(Residue* resEntry){
  if (resEntry->getResId()){
    resVec.push_back(resEntry);
  }
}

void Chain::addAtom(Atom* atmEntry){
  if (atmEntry->getAtmNum()){
    atmVec.push_back(atmEntry);
  }
}

Atom* Chain::getAtom (const unsigned int& element){
	if (element >= atmVec.size()){
		return NULL;
	}
	else{
  	return atmVec.at(element);
	}
}

Residue* Chain::getResidue (const unsigned int& element){
	if (element >= resVec.size()){
		return NULL;
	}
	else{
  	return resVec.at(element);
	}
}

std::string Chain::getChainId(){
  return this->getAtom(0)->getChainId();
}

unsigned int Chain::getAtmVecSize(){
  return atmVec.size();
}

unsigned int Chain::getResVecSize(){
  return resVec.size();
}

//void Chain::setSel(bool selin){
//  sel=selin;
//}

//bool& Chain::getSel(){
//  return sel;
//}

void Chain::selAll(){
//  sel=true;
  unsigned int i;
//  for (i=0; i< this->getResVecSize(); i++){
//    this->getResidue(i)->setSel(true);
//  }
  for (i=0; i< this->getAtmVecSize(); i++){
    this->getAtom(i)->setSel(true);
  }
}

void Chain::deselAll(){
//  sel=false;
  unsigned int i;
//  for (i=0; i< this->getResVecSize(); i++){
//    this->getResidue(i)->setSel(false);
//  }
  for (i=0; i< this->getAtmVecSize(); i++){
    this->getAtom(i)->setSel(false);
  }  
}
