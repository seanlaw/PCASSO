//Sean M. Law

#ifndef ANALYZE_H
#define ANALYZE_H

#include "Enum.hpp"

#include <vector>
#include <string>

//Forward Declaration
class Molecule;
class Vector;
class DTree;

//Abstract base class (cannot create instance of it!)
class Analyze {
	private:
    //Since this is an abstract base class
    //all members need to be accessed via a function or passed
    //directly to the analysis function or set the members as protected
    std::vector<std::string> sel;
		std::vector<Molecule*> mol;
		std::vector<double> tdata; //Time dependent data, maybe for averaging
		std::vector<std::vector<double> > fdata; //Frame data, cleared after each frame
    std::vector<unsigned int> modes;
		int ndata; //Total number of datapoints
		bool resel; //Re-do selection for each analysis, not implemented yet
    std::string ifile;
    std::string ofile;

	public:
		Analyze ();
		virtual ~Analyze ();
		void addSel(const std::string& selin);
		std::string getSel(const int& element);
		unsigned int getNSel();
    void addMol(Molecule* molin);
		void setMol(const int& element, Molecule* molin);
		void clearMol();
		void resizeNMol(const int sizein);
		Molecule* getMol(const int& element);
		unsigned int getNMol();
		void setNData(const int& ndatain);
		int& getNData();
		std::vector<double>& getTDataVec();
		std::vector<std::vector<double> >& getFDataVec();
    void setInput(const std::string& fin);
    std::string getInput();
    void setOutput (const std::string& fin);
    std::string getOutput();
    void addModes(const std::vector<unsigned int>& modesin);
    std::vector<unsigned int>& getModes();
		
		//Virtual functions
    virtual void readTopology(Molecule* molin, std::string topin="");
		virtual void setupMolSel(Molecule* molin);
		virtual void preAnalysis(Molecule* molin, std::string topin=""); 
    virtual void preAnalysis();
		virtual void runAnalysis() =0; //Pure virtual function
		virtual void postAnalysis();

		//Analysis functions
		static Vector centerOfGeometry(Molecule* mol, bool selFlag=true);
		static double distance (const Vector& u, const Vector& v);
		static double angle (const Vector& u, const Vector& v, const Vector& w);
    static double dihedral (const Vector& t, const Vector& u, const Vector& v, const Vector& w);
		static double distance (Molecule* sel1, Molecule* sel2, bool selFlag=true);
		static double angle (Molecule* sel1, Molecule* sel2, Molecule* sel3, bool selFlag=true);
	  static double dihedral (Molecule* sel1, Molecule* sel2, Molecule* sel3, Molecule* sel4,bool selFlag=true);
		static void pairwiseDistance(Molecule *mol, std::vector<std::vector<double> >& pdin);
		static void allAnglesDihedrals(Molecule *mol, std::vector<std::vector<double> >& anglesin);
		static void pcasso(Molecule* mol, std::vector<std::vector<double> > &fdataIO);
    //static std::vector<double> gyration(Molecule* mol);
};

//Derived classes

class AnalyzePcasso: public Analyze {
	private:
		PcassoOutEnum pout;
		std::vector<DTree *> t;

	public:
		AnalyzePcasso(std::string delim=":");
		void setOutType(PcassoOutEnum pin);
		PcassoOutEnum getOutType();
		void preAnalysis(Molecule* molin, std::string fin="");
		void runAnalysis();
		void postAnalysis();
};

#endif
