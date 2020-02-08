/**@file solver.h
 * @brief Header file for TSP Solver.
 *       
 * Detailed description follows here.
 * @author     : Abhinav Modi
 * @created on : Sep 24, 2019
 * @copyright  : This code is developed for Homework 
 *				 assignments for the course CMSC818B. 
 *				 Do not use without citation.
 */

#ifndef SOLVER_H_
#define SOLVER_H_

#include<stdio.h>
#include<iostream>
#include<vector>
#include<algorithm>
#include<string>
#include<typeinfo>
#include<stack>
#include<stdlib.h>
#include<fstream>
#include<cstdlib>
#include<iomanip>

using std :: cout;
using std :: cin;
using std :: endl;
using std :: vector;
using std :: stack;

/**
 * Structure for reading the data from the provided .tsp file
 */
struct Node {
	int key;
	int xCord;
	int yCord;
};


/**
 * Structure for creating the Minimum Spanning Tree
 */
struct Tree {
	int *parent; int *rank;
	int n;

	//Constructor
	Tree(int n){
		//Memory allocation
		this->n = n;
		parent = new int[n+1];
		rank = new int[n+1];

		//Initially all vertices are parts of
		//different trees and have rank 0
		for (int i=0; i <= n; i++) {
			rank[i] = 0;

			//every element is parent of itself
			parent[i] = i;
		}
	}

	//Utility functions for the tree
	//Find the parent tree of vertex 'v1'
	int findTree(int v1){
		/*Form a linked list of parents of the tree to
		 * find the subtree connected to the tree*/
		if (v1 != parent[v1])
			parent[v1] = findTree(parent[v1]);
		return parent[v1];
	}

	//Merge the subtrees collected to form a larger tree
	void mergeTrees(int tree1, int tree2){
		tree1 = findTree(tree1), tree2 = findTree(tree2);

		if (rank[tree1] > rank[tree2])
			parent[tree2] = tree1;
		else //If rank[tree1] <= rank[tree2]
			parent[tree1] = tree2;

		if (rank[tree1] == rank[tree2])
			rank[tree2]++;
	}
};

/**
 * Creating a shortcut for an integer pair
 */
typedef std::pair<int, int> nodePair;


class TspSolver{
	//Generate Graph
	// <edge_cost,<source, destination>>
	vector<std::pair<double, nodePair>> edges;

	/**
     * @brief: Utilityy function to add edges to the graph
	 * @param: node1 - index of the first node in the edge
	 * @param: node2 - index of the second node in the edge
	 * @param: distance - length of the edge
     * @return: None
     */
	void addEdge(int node1, int node2, double distance) {
		edges.push_back(std::pair<double, nodePair> \
						(distance, nodePair(node1, node2)));
	}

	/**
     * @brief: Checks if the input edge is the last edge in the edges list
     * @param: node1 - index of the first node in the edge
	 * @param: node2 - index of the second node in the edge
	 * @param: distance - length of the edge
     * @return: bool - false if the edge is the last in the list true otherwise
     */
	bool checkEdge(int node1, int node2, double distance) {
		if (std::find(edges.begin(), edges.end(), std::pair<double, nodePair> \
					(distance, nodePair(node2, node1))) != edges.end()) {
			return false;
		}
		else {return true;}
	}

public:
	vector<Node> nodes;			
	vector<std::pair<double, nodePair>> tree;
	vector<int> tour;
	stack<int> visited;

	//Define Constructor for the class
	TspSolver(vector<Node> &nodeList){ nodes = nodeList;};

	double calcDistance(Node &node1, Node &node2) {
	  double distance = 0.0;
	  distance = std::pow(std::pow((node1.xCord - node2.xCord),2) \
	                      + std::pow((node1.yCord - node2.yCord),2),0.5);
	  return distance;
	}

	/**
     * @brief: generate graph using all nodes
     * @param: None
     * @return: None
     */
	void genGraph() {
		double distance_ij;
		int it = nodes.size();
		for (int i = 0; i < it-1; i++){
			for (int j = 0; j < it; j++){
				if (i !=j){
					distance_ij = calcDistance(nodes[i], nodes[j]);
					if (!checkEdge(nodes[i].key,nodes[j].key,distance_ij)){
						continue;}

					//If node does not exist in the edges add it to the vector
					addEdge(nodes[i].key, nodes[j].key, distance_ij);
				};
			}
		}
		cout << "Number of edges in the graph are " << edges.size() << endl;
	};

	/**
     * @brief: generate MST based on Kruskal's algorithm
     * @param: None
     * @return: None
     */
	void genMST() {

		//Sort edges in increasing order on basis of cost
		std::sort(edges.begin(), edges.end());

		//Initialize MST as an empty Tree
		Tree mst(nodes.size());

		//Iterate through all sorted edges
		vector< std::pair<double,nodePair>> :: iterator currEdge;
		for (currEdge=edges.begin(); currEdge!=edges.end(); currEdge++) {
			int node1 = currEdge->second.first;
			int node2 = currEdge->second.second;
			//Find the tree corresponding to the nodes
			int tree1 = mst.findTree(node1);
			int tree2 = mst.findTree(node2);

			if (tree1 != tree2) {
				//Add the current edge to the mst
				mst.mergeTrees(tree1, tree2);
				tree.push_back(std::pair<double, nodePair> \
							  (currEdge->first,nodePair(node1, node2)));
			}
		}
	}


	/**
     * @brief: Finds the next node in the tree in a depth-wise manner
     * @param: currNode - integer index of the current node
     * @return: int - integer index of the next node
     */
	int nextNode(int currNode) {
	  int nextNode = currNode;
	  std::sort(tree.begin(), tree.end());
	  vector< std::pair<double,nodePair>> :: iterator search=tree.begin();
	  int len = tree.size();
	  for (search; search!=tree.end();search++) {
	    int x = search->second.first;
	    int y = search->second.second;

	    if (x == currNode) {
	      nextNode = y;
	      tree.erase(search);
	      return nextNode;
	    }
	    if (y == currNode) {
	      nextNode = x;
	      tree.erase(search);
        return nextNode;
	    }
	  }
	  return nextNode;
	}

	/**
     * @brief: checks if the node was previously visited
     * @param: node - integer index of the current node
     * @return: bool - false if the node was not visited true otherwise
     */
	bool ifVisited(int node) {
	  for (int i=0; i<tour.size();i++) {
	    if (tour[i] == node) {
	      return true;
	    }
	  }
	  return false;
	}

	/**
     * @brief: nodeTraverse is the main function which traverses all nodes in the graph
	 * 		   to find the minimum length tour
     * @param: start - integer key of the start node
     * @return: None
     */
	void nodeTraverse(int start) {

	  int currNode=0, childNode = 0;
	  std::sort(tree.begin(), tree.end());
	  vector< std::pair<double,nodePair>>  fullTree;
	  std::copy(tree.begin(),tree.end(), back_inserter(fullTree));
	  currNode = start;
	  tour.push_back(currNode);
	  visited.push(currNode);
	  while (tree.size() != 0) {

	    childNode = nextNode(currNode);

	    if (ifVisited(childNode)) {
	      visited.pop();
        currNode = visited.top();
	    }
	    else {
        tour.push_back(childNode);
        visited.push(childNode);
        currNode = childNode;
	    }
	  }
	  tour.push_back(start);
	}

	/**
     * @brief: Calculates the length of the final tour generated by the solver
     * @param: None
     * @return: euclidean length of the tour
     */
	double tourLength() {
	  int numVertices = tour.size();
	  double tourLen = 0.0;
	  Node c1, c2;
	  for (int i=0; i < numVertices - 1; i++) {
	    //Find city1 (x,y) coord
	    vector<Node>::iterator city1 = nodes.begin();
	    std::advance(city1, tour[i]-1);

	    //Find city1 (x,y) coord
      vector<Node>::iterator city2 = nodes.begin();
      std::advance(city2, tour[i+1]-1);

      //Calculate distance between the two cities
      tourLen += calcDistance(*city1, *city2);
	  }
	  cout << "Length of tour is: " << tourLen << endl;
	  return tourLen;
	}

    /**
     * @brief: prints the edges of the minimium spanning tree generated 
     * @param: None
     * @return: None
   	 */
	void printTree() {
	  double sum=0;
    vector< std::pair<double,nodePair>> :: iterator iter;
    for (iter=tree.begin(); iter!=tree.end(); iter++) {
      cout << iter->second.first << " - " << iter->second.second << endl;
      sum += iter->first;
    }
    cout << "Length of MST "<< sum << endl;
	}

    /**
	 * @brief: prints the final tour predicted by the algorithm
	 * @param: None
	 * @return: None
	 */
	void printTour(){
	  vector<int>::iterator n;
	  for (n=tour.begin(); n != tour.end(); n++) {
	    cout << *n << endl;
	  }
	}
};

/**
 * Data handling class for reading and writing .tsp files 
 */
class DataHandling {
 private:
  std::string name; //Name of the file

 public:

  /**
   * @brief: readGraph funtion is to read from .tsp file and 
   * 		 save the data in a vector of Nodes
   * @param: ff - ifstream object to read input stream from the file
   * @return: vector of nodes
   */
  vector<Node> readGraph(std::ifstream& ff){
    vector<Node> nodes;
    std::string line;
    int x=0, y=0, key=0;

    std::cout << "data stream" << endl;

    std::getline(ff,line);
    name = line.substr(line.find(":")+2);
    cout << name << endl;

    for (int i=0; i<4;i++){
      std::getline(ff,line);
    }

    cout << line << endl;
    int count=0;
    while (std::getline(ff, line)) {
      if (line.compare("EOF")!=0){

        ff >> key >> x >> y;

        //Create a new node for each line
        Node currNode;
        currNode.key = key;
        currNode.xCord = x;
        currNode.yCord = y;

        //Push it in the vector
        nodes.push_back(currNode);
        count++;
      };
    }
    nodes.pop_back();
    return nodes;
  };

  /**
   * @brief: genOutputFile funtion is to write the generated tour to a .tsp file  
   * @param: solver - TspSolver
   * @return: vector of nodes
   */
  void genOutputFile(TspSolver* solver) {
    //generate output file with name "name"
    std::ofstream outFile;
    std::string filename = name+".out.tour";
    outFile.open(filename);
    outFile << "NAME : " << name+".out.tour" << endl;
    outFile << "COMMENT : " << "Tour for " + name + ".tsp"
            << " (Length " << solver->tourLength() << ")" << endl;
    outFile << "TYPE : TOUR" << endl;
    outFile << "DIMENSION : " << solver->nodes.size() << endl;
    outFile << "TOUR_SECTION" << endl;
    vector<int>::iterator n;
    for (n=solver->tour.begin(); n != solver->tour.end(); n++) {
      outFile << *n << endl;
    }
    outFile.close();
  }


  vector<Node> inputreadGraph(std::ifstream& ff){
    vector<Node> nodes;
    std::string line;
    int x=0, y=0, key=0;

    std::cout << "data stream" << endl;

    std::getline(ff,line);
    name = line.substr(line.find(":")+2);
    cout << line << endl;
    for (int i=0; i<4;i++){
      std::getline(ff,line);
    }

    cout << line << endl;
    int count=0;
    while (std::getline(ff, line)) {
      if (line.compare("EOF")!=0){

        ff >> key >> x >> y;
        //Create a new node for each line
        Node currNode;
        currNode.key = key;
        currNode.xCord = x;
        currNode.yCord = y;

        //Push it in the vector
        nodes.push_back(currNode);
        count++;
      };
    }
    nodes.pop_back();
    return nodes;
  };
};

#endif /* SOLVER_H_ */
