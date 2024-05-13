// Roni Naftalovich 319049060, roni100400100400@gmail.com
#ifndef ALGORITHMS_HPP
#define ALGORITHMS_HPP
#include <iostream>
#include <vector>
#include <queue>
#include <stack>
#include <unordered_map>
#include <algorithm> 
#include <climits>
#include "Graph.hpp"





using namespace std;

namespace ariel
{

    class Algorithms
    {
    public:
        static bool isConnected(Graph &g);  

        static string shortestPath(Graph &g, int start, int end);

        static bool isContainsCycle(Graph &g);

        static void isBipartite(Graph &g);

        static void negativeCycle(Graph &g);
    };
}
#endif