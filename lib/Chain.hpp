//Sean M. Law

#ifndef CHAIN_H
#define CHAIN_H

#include <vector>
#include <string>

//Forward Declaration
class Residue;
class Atom;

class Chain {
  private:
    std::vector<Residue *> resVec; //Vector of residue pointers
    std::vector<Atom *> atmVec; //Vector of atom pointers
//    bool sel;

  public:
    Chain();

    void reset();

    void addResidue(Residue* resEntry);
    void addAtom(Atom* atmEntry);

    //Get Chain info
    Atom* getAtom(const unsigned int& element);
    Residue* getResidue(const unsigned int& element);
    std::string getChainId();
    unsigned int getAtmVecSize();
    unsigned int getResVecSize();
//    void setSel(bool selin);
//    bool& getSel();
    void selAll();
    void deselAll();
};

#endif
