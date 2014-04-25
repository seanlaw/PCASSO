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

#ifndef MOLECULE_H
#define MOLECULE_H

#include "Prmtop.hpp"
#include "Coor.hpp"
#include "Enum.hpp"

#include <vector>
#include <map>

//Forward Declaration
class Chain;
class Residue;
class Atom;

//Base class
class Molecule {
  private:
    std::vector<Chain*> chnVec;
    std::vector<Residue*> resVec;
    std::vector<Atom*> atmVec;
    bool copyFlag; //This molecule is a copy if true
    std::map< std::string, std::vector<bool> > storedSel;
    std::string remarks;
    bool iCodeFlag;
    Prmtop toppar;
    unsigned int year;
    std::string exp;

  public:
    Molecule(); //Constructor
    virtual ~Molecule();
    static Molecule* readPDB (const std::string ifile, const int model=0, const std::string format="", const bool hetFlag=true);
    static Molecule* readPDB (const std::string ifile, const std::string format, const bool hetFlag=true);
//    std::string writePDB (bool selFlag=true, bool print=true, bool chnFlag=false);
    std::string writePDB (bool selFlag, bool print, bool chnFlag);
    std::string writePDB (bool selFlag, bool print);
    std::string writePDB (bool chnFlag);
    std::string writePDB ();
    Molecule* clone(bool selFlag=true, bool keep=true);
    Molecule* copy(bool selFlag=true);
    void cat (Molecule* catmol, bool selFlag=true, bool keep=true);
    void addAtom(Atom* atmEntry);
    Atom* getAtom(const unsigned int& element); 
    void addChain(Chain* chnEntry);
    void addResidue(Residue* resEntry);
    Chain* getChain(const unsigned int& element);
    unsigned int getChnVecSize();
    unsigned int getResVecSize();
    unsigned int getAtmVecSize();
    std::vector<Atom*>& getAtmVec();
    std::vector<Atom*> getAtmVecClone();
    Residue* getResidue(const unsigned int& element);
    void readTopology(const std::string& topin);
    void readParameter(const std::string& prmin);
    void setMass();
    void setCharge();
    void selAll();
    void deselAll();
    void select(std::string sel);
    unsigned int getNAtom();
    unsigned int getNAtomSelected(); //Determining this on the fly is a good safety measure
    void setCopyFlag(bool copyFlagIn=false);
    bool getCopyFlag();
    void storeSel(std::string key="tmp");
    void recallSel(std::string key="tmp");
    void eraseSel(std::string key="tmp");
    void zeroCoor();
    void addRemark(const std::string& remin);
    void clearRemark();
    std::string getRemark();
    bool checkICode();
    void setICodeFlag(bool iCodeFlagIn=false);
    bool getICodeFlag();
    void assignAtmInx();
    void resetAtmInx();
    void setYear(const unsigned int& yearin);
    unsigned int getYear();
    void setExp(const std::string& expin);
    std::string getExp();

    void modPseudoCenter();
    void pcasso (std::string dsspin="", PcassoOutEnum out=PREDICT);

    //Virtual functions
    virtual void format();
};

class MoleculeCHARMM: public Molecule{
  public:
    void format();
};

#endif
