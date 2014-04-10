//Sean M. Law

#ifndef PCASSO_H
#define PCASSO_H

#include <string>

class PCASSO {
	private:
		//Don't put std::vector<std::string> t here!
		//Can't instantiate object when using static call!

	public:
		static std::string getTree (unsigned int elem);
		static unsigned int getNTree();

};

#endif
