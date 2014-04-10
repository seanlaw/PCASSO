//Sean M. Law
#ifndef PDB_H
#define PDB_H

#include <map>
#include <string>

class Molecule;
class Residue;
class Atom;

class PDB {
  private:
    std::map<std::string, int> chnMap;
		std::string format; //Output format

  public:
    PDB();
    static void writePDBFormat (Molecule* mol, std::ostringstream &out, bool selFlag=true, bool chnFlag=false);
    static Molecule* readPDB (const std::string ifile, const int model=0, const std::string format="", const bool hetFlag=true, const bool remFlag=false);
    Atom* processAtomLine (std::string line, Atom* lastAtom);
		static std::string formatCHARMMResName (Atom* atmEntry);
		static int formatCHARMMResId(Atom* atmEntry, Residue* lastRes, Residue* nextRes);
};

#endif
