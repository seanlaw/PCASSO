//Sean M. Law

#ifndef BTREE_H
#define BTREE_H

struct BTreeNode{
  int key_value;
	BTreeNode *left;
	BTreeNode *right;
};


class BTree {
	private:
		void delBTree(BTreeNode *leaf);
		void addBTree(int key, BTreeNode *leaf);
		BTreeNode *searchBTree(int key, BTreeNode *leaf);
		BTreeNode *root;

	public:
		BTree();
		~BTree();

		void delBTree();
		void addBTree(int key);
		BTreeNode *searchBTree(int key);
};

#endif
