//Sean M. Law

#ifndef DTREE_H
#define DTREE_H

#include <string>
#include <vector>
#include <map>

struct DTreeNode{
  double key_value;
	unsigned int inx; //1-D vector<double> index, no response column!
	DTreeNode *left;
	DTreeNode *right;
	std::string cls;

	DTreeNode* getDTreeNodeLeft();
	DTreeNode* getDTreeNodeRight();

};


class DTree {
	private:

		void delDTree(DTreeNode *leaf);
		DTreeNode* root;
		std::map<std::string, std::string> classMap;
		
		void addDTree(double key, DTreeNode *leaf, unsigned int index, std::string classin="");
		std::string getDTreeClass(DTreeNode *leaf, const std::vector<double> &fin);
		void genDTree(DTreeNode *&leaf, std::vector<std::string> &t, unsigned int &inx, std::string delim=":");

	public:
		DTree();
		~DTree();

		void delDTree();

		void addDTree(double key, unsigned int index, std::string classin="");

		void genDTree(std::vector<std::string> &t, std::string delim=":");

		DTreeNode* getDTreeRoot();
		std::string getDTreeClass(const std::vector<double> &fin);

		unsigned int getClassSize();
		void addClass(std::string valin, std::string classin);
		void delClass(std::string classin);
		std::string getClass(std::string valin);
};

#endif


