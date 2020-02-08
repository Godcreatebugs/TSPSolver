/**@file Wrapper.cpp
 * @brief Main file for TSP Solver.
 *       
 * Detailed description follows here.
 * @author     : Abhinav Modi
 * @created on : Sep 24, 2019
 * @copyright  : This code is developed for Homework 
 *				 assignments for the course CMSC818B. 
 *				 Do not use without citation.
 */

#include<chrono>
#include "solver.h"
using namespace std :: chrono;

int main(int argc, char* argv[]){
	//Initialize vector to save Nodes
	vector<Node> nodes;
	DataHandling files;

	std::ifstream datafile {argv[1]}; //first arg is filename
	nodes = files.readGraph(datafile);
	datafile.close();

	TspSolver tsp(nodes);

	tsp.genGraph();
	tsp.genMST();

	//initialize start node
	int start = 1;
	tsp.printTree();
	auto start_time = high_resolution_clock::now();

	tsp.nodeTraverse(start);

	auto end_time = high_resolution_clock::now();

	auto duration = duration_cast<microseconds>(end_time - start_time);

	cout << "Time taken for finding the tour is: "
	     << duration.count() << " micro secs" << endl;

	tsp.printTour();
	tsp.tourLength();
	files.genOutputFile(&tsp);

	return 0;
}
