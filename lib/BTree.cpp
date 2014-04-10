//Sean M. Law

#include "BTree.hpp"

#include <cstdlib>

BTree::BTree(){
	root=NULL;
}

BTree::~BTree(){
	delBTree();
}

void BTree::delBTree(BTreeNode *leaf){
	if (leaf != NULL){
		delBTree(leaf->left);
		delBTree(leaf->right);
		delete leaf;
	}
}

void BTree::addBTree(int key, BTreeNode *leaf){
	if (key < leaf ->key_value){
		if (leaf->left != NULL){
			addBTree(key, leaf->left);
		}
		else{
			leaf->left=new BTreeNode;
			leaf->left->key_value=key;
			leaf->left->left=NULL; //Sets left child of child node to NULL
			leaf->left->right=NULL; //Sets right child of child node to NULL
		}
	}
	else{
		if (leaf->right != NULL){
			addBTree(key, leaf->right);
		}
		else{
			leaf->right=new BTreeNode;
			leaf->right->key_value=key;
			leaf->right->left=NULL; //Sets left child of child node to NULL
			leaf->right->right=NULL; //Sets right child of child node to NULL
		}
	}
}

BTreeNode* BTree::searchBTree(int key, BTreeNode *leaf){
	if (leaf != NULL){
		if (key == leaf->key_value){
			return leaf;
		}
		else if (key < leaf->key_value){
			return searchBTree(key, leaf->left);
		}
		else{
			return searchBTree(key, leaf->right);
		}
	}
	else{
		return NULL;
	}
}

void BTree::addBTree(int key){
	//Public version of the addBTree function
	//Takes care of when the root node has not been initialized
	if (root != NULL){
		addBTree(key, root);
	}
	else{
		root=new BTreeNode;
		root->key_value=key;
		root->left=NULL;
		root->right=NULL;
	}
}

BTreeNode* BTree::searchBTree(int key){
	//Public version of the searchBTree function
	//Starts search recursion starting at root node
	return searchBTree(key, root);
}

void BTree::delBTree(){
	//Public version of the delBTree function
	//Starts delete recursion starting at root node
  delBTree(root);
}

