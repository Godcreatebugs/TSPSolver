# CMSC818B Homework 1 - Travelling Salesman Problem
The Travelling Salesman Problem is to find the shortest possible route that visits every city in a given list of cities and returning to the origin for the tour to be complete. It is an NP-hard problem in combinatorial optimization. This implementation of the Travelling Salesman Problem is based on Kruskal's agorithm to find minimum spanning trees.

## Algorithm Overview
For this project I have followed MST algorithm which requires a weighted graph as an input and then tries to find a tour which visits each vertex at least once.

The basic approach is as follows:
1. Find the Minimum Spanning Tree(MST).
2. Traverse the edges in a depth-first fashion.
3. When going up the tree, skip an already visited vertex and add a shortcut to the next unvisited one.

## Compiling and running the code
- The code is present in the Code folder
- To build the code, open the `Code` folder in terminal and type:
```
g++ -std=c++11 Wrapper.py -o tsp
```
- Once the code is built use the following command to run the code:
```
./tsp <path to .tsp file>
```
eg.,
```
./tsp ../eil51.tsp
```
This will generate the output file containing the tour with name eil51.out.tour
