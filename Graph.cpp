// Roni Naftalovich 319049060, roni100400100400@gmail.com

#include <iostream>
#include <vector>
#include "Graph.hpp"

using namespace std;
using namespace ariel;

Graph::Graph(){
    this->isDirected = false; //temporarly
    this->isWeighted = false; //temporarly
    this->numVertices = 0;
}

Graph::Graph(bool isDirected)
{
    this->isDirected = isDirected;
    this->isWeighted = false; //temporarly
    this->numVertices = 0;
}

size_t Graph::getNumVertices() const
{
    return this->numVertices;
}

vector<vector<int>> Graph::getAdjMatrix() const
{
    return this->adjMatrix;
}

bool Graph::getisDirected() const
{
    return this->isDirected;
}


void Graph::loadGraph(vector<vector<int>> &matrix)
{
    if (matrix.size() !=  matrix[0].size() || matrix.size()==0) {
            throw std::invalid_argument("Graphs must be of the same size.");
            return;
        }

    this->numVertices = matrix.size();
    this->adjMatrix = matrix;
    for (size_t i = 0; i < adjMatrix.size(); i++)
    {
        for (size_t j = 0; j < adjMatrix[i].size(); j++)
        {
            if (i == j && adjMatrix[i][j] != 0)
            {
               cerr << "Invalid adjacency matrix - diagonal must be zero "<< endl;
               return;
            }
        }
    }
    this->isDirectedAndWeighted(matrix);
}

void Graph::printGraph()
{
    if (this->isDirected==true)
    {
        cout << "The graph is directed, containing " << this->numVertices << " vertices" << endl;
    }
    else
    {
        cout << "The graph is undirected, containing " << this->numVertices << " vertices" << endl;
    }
}

bool Graph::getisWeighted(){
    return this->isWeighted;
}

void Graph::isDirectedAndWeighted(vector<vector<int>> &matrix) // for graph being undirected it first need to be symmetric
{
    for (size_t i = 0; i < matrix.size(); i++)
    {
        for (size_t j = 0; j < matrix[i].size(); j++)
        {
            if (this->isDirected== false && matrix[i][j] != matrix[j][i])
            {
                this->isDirected= true;
            }
            if(this->isWeighted == false && matrix[i][j] !=0 && matrix[i][j]!=1)
            {
                this->isWeighted=true;
            }
        }
    }
    return;
}