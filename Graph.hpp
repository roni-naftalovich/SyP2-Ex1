// Roni Naftalovich 319049060, roni100400100400@gmail.com

#ifndef GRAPH_HPP
#define GRAPH_HPP

#include <iostream>
#include <vector>
#include <stdexcept>

using namespace std;

namespace ariel
{
    class Graph
    {
    private:
        vector<vector<int>> adjMatrix;
        size_t numVertices;
        bool isDirected;
        bool isWeighted;

    public:
        // Constructor that getting  if the graph is directed or not
        Graph(bool isDirected = false); // default constructor is not directed

        size_t getNumVertices() const;
        
        //getter to adjMatrix
        vector<vector<int>> getAdjMatrix() const;
        
        //getter to isDirected
        bool getisDirected() const;


        // Function to load graph from an adjacency matrix, like a set function
        void loadGraph(vector<vector<int>> &matrix); // passing it not by reference because we don't want to change it the orginal matrix

        // Function to print the adjacency matrix
        void printGraph();

        //getter to isWeighted
        bool getisWeighted();

        void isDirectedAndWeighted(vector<vector<int>> &matrix); 
    };
} 

#endif // GRAPH_HPP