// Roni Naftalovich 319049060, roni100400100400@gmail.com

#include <iostream>
#include <vector>
#include "Graph.hpp"
#include "Algorithms.hpp"
#include <string>
#include <climits>
#include <limits>
#include <queue>
#include <stack>
#include <algorithm>
#include <sstream>


using namespace std;
using namespace ariel;



std::vector<size_t> DFS(Graph &g, size_t vertex, std::vector<int> &verticesColors) {
    std::vector<size_t> visitedVertices; // Vector to store visited vertices
    std::stack<size_t> stack; // Stack for DFS traversal

    // Push the starting vertex onto the stack
    stack.push(vertex);

    // Mark the starting vertex as visited
    verticesColors[vertex] = 1; // Mark as gray (being visited)

    // Perform DFS traversal
    while (!stack.empty()) {
        size_t currentVertex = stack.top();
        stack.pop();
        visitedVertices.push_back(currentVertex); // Mark the current vertex as visited

        // Get the adjacency list of the current vertex
        std::vector<int> adjacencyList = g.getAdjMatrix()[currentVertex];

        // Visit adjacent vertices
        for (size_t i = 0; i < adjacencyList.size(); ++i) {
            if (adjacencyList[i] != 0 && verticesColors[i] == 0) { // If the vertex is adjacent and not visited
                stack.push(i); // Push the adjacent vertex onto the stack
                verticesColors[i] = 1; // Mark as gray (being visited)
            }
        }
    }

    return visitedVertices;
}

std::vector<std::vector<size_t>> DFSforest(Graph &g) {
    size_t V = g.getNumVertices();
    std::vector<int> Colors(V, 0); // Initialize all vertices as 0  means not visited)
    std::vector<std::vector<size_t>> forest; // DFS forest output

    for (size_t i = 0; i < V; i++) {
        if (Colors[i] == 0) { // i not visited yet
            std::vector<size_t> tree = DFS(g, i, Colors); // Start DFS visit from the vertex
            forest.push_back(tree); // Push the visit tree into the forest
        }
    }

    return forest;
}

bool directedIsConnected(Graph &g)
{
   size_t numVertices = g.getNumVertices();
    if (numVertices == 0)
        return true; // An empty graph is considered connected

    for (size_t i = 0; i < numVertices; ++i) {
        std::vector<bool> visited(numVertices, false);
        std::queue<int> q;
        q.push(i); // Start BFS from vertex i
        visited[i] = true;

        while (!q.empty()) {
            int u = q.front();
            q.pop();

            for (size_t j=0; j< numVertices; j++) {
            if(g.getAdjMatrix()[static_cast<size_t>(u)][j]!=0){
                if (!visited[j]) {
                    q.push(j);
                    visited[j] = true;
                }
            }}
        }

        if (std::find(visited.begin(), visited.end(), false) != visited.end())
            return false; // If any vertex is not reachable, graph is not connected
    }

    return true;
}



bool undirectedIsConnected(Graph &g)
{
    size_t numVertices = g.getNumVertices();
    if (numVertices == 0)
        return true; // An empty graph is considered connected

    std::vector<bool> visited(numVertices, false);
    std::queue<int> q;

    q.push(0); // Start BFS from vertex 0
    visited[0] = true;

    while (!q.empty()) {
        int u = q.front();
        q.pop();

        for (size_t j=0; j< numVertices; j++) {
            if(g.getAdjMatrix()[static_cast<size_t>(u)][j]!=0){
                if (!visited[j]) {
                    q.push(j);
                    visited[j] = true;
                }
            }
        }}

    return std::all_of(visited.begin(), visited.end(), [](bool v) { return v; });
}


bool dfsForCycleDetection(size_t vertex, int parent, Graph& graph, vector<bool>& visited, vector<int>& cycle)
{
    visited[vertex] = true;
    for (size_t i = 0; i < graph.getNumVertices(); ++i) {
        if (graph.getAdjMatrix()[vertex][i] > 0) { // Check if there is an edge between the vertices
            if (!visited[i]) {
                cycle.push_back(vertex); // Add current vertex to cycle
                if (dfsForCycleDetection(i, vertex, graph, visited, cycle)) {
                    return true;
                }
                cycle.pop_back(); // Remove current vertex if no cycle found in subtree
            } 
            else if (i != parent) {
                cycle.push_back(vertex); // Add current vertex to cycle
                cycle.push_back(i); // Add neighbor to cycle
                return true; // Cycle detected
            }
        }
    }
    return false; // No cycle found
}

bool isContainCylcleDirected(Graph &graph)
{
    size_t n = graph.getNumVertices();
    vector<bool> visited(n, false);

    for (size_t v = 0; v < n; ++v) {
        if (!visited[v]) {
            vector<int> cycle;
            if (dfsForCycleDetection(v, -1, graph, visited, cycle)) {
                cout << "Cycle found: ";
                for (size_t i = 0; i < cycle.size() ; ++i) {
                    cout << cycle[i];
                    if (i < cycle.size() - 1) {
                         cout << "->";
        }
    }
    cout << endl;
}
                
                return true; // Cycle detected
            }
        }
    


    return false; // No cycle found
}

bool dfsCycleDetectionUnDirected(Graph &g, size_t v, int parent, std::vector<bool> &visited) {
    visited[v] = true;

    // Go through all adjacent vertices of v
    for (size_t i = 0; i < g.getNumVertices(); ++i) {
        // If an adjacent vertex is not visited yet, then recur for it
        if (g.getAdjMatrix()[v][i]) { 
            if (!visited[i]) {
                if (dfsCycleDetectionUnDirected(g, i, v, visited))
                    return true;
            }
            // If an adjacent vertex is visited and is not the parent of the current vertex,
            // then there is a cycle in the graph
            else if (i != parent)
                return true;
        }
    }
    return false;
}

bool isContaionCylcleUnDirected(Graph &g) {
    size_t V = g.getNumVertices();
    std::vector<bool> visited(V, false);

    // Call the recursive helper function to detect cycle in different DFS trees
    for (size_t i = 0; i < V; ++i) {
        if (!visited[i]) {
            if (dfsCycleDetectionUnDirected(g, i, -1, visited))
                return true;
        }
    }
    return false;
}




std::string BFS(Graph &g, int start, int end) {
    size_t numVertices = g.getNumVertices();
    std::queue<int> q;
    std::unordered_map<int, int> parent;
    std::vector<bool> visited(numVertices, false);

    q.push(start);
    visited[static_cast<size_t>(start)] = true;

    while (!q.empty()) {
        int u = q.front();
        q.pop();

        if (u == end) {
            // Reconstruct path
            std::stringstream path;
            path << end;
            int current = end;
            while (current != start) {
                current = parent[static_cast<size_t>(current)];
                path << " -> " << current;
            }
            return path.str();
        }

                    for (size_t j=0; j< numVertices; j++) {
                        if(g.getAdjMatrix()[static_cast<size_t>(u)][j]!=0){
                            if (!visited[j]) {
                                q.push(j);
                                parent[j] = u;
                                visited[j] = true;
                            }
                        }
                    }}

    return "No path exist";
}

string dijkstra(Graph &g, int start, int end) {
    size_t numVertices = g.getNumVertices();
    size_t startIndex = static_cast<size_t>(start);
    size_t endIndex = static_cast<size_t>(end);

    vector<int> dist(numVertices, numeric_limits<int>::max());
    vector<size_t> prev(numVertices, numeric_limits<size_t>::max());
    priority_queue<pair<int, size_t>, vector<pair<int, size_t>>, greater<pair<int, size_t>>> pq;

    dist[startIndex] = 0;
    pq.emplace(0, startIndex);

    while (!pq.empty()) {
        size_t u = pq.top().second;
        pq.pop();

        if (u == endIndex) {
            break;
        }

        for (size_t v = 0; v < numVertices; v++) {
            if (g.getAdjMatrix()[u][v] != 0) {
                int weight = g.getAdjMatrix()[u][v];

                if (dist[u] + weight < dist[v]) {
                    dist[v] = dist[u] + weight;
                    prev[v] = u;
                    pq.emplace(dist[v], v);
                }
            }
        }
    }

    if (dist[endIndex] == numeric_limits<int>::max()) {
        return "No path exists";
    } else {
        stringstream pathStream;
        vector<size_t> path;
        for (size_t at = endIndex; at != numeric_limits<size_t>::max(); at = prev[at]) {
            path.push_back(at);
        }
        reverse(path.begin(), path.end());

        for (size_t i = 0; i < path.size(); ++i) {
            pathStream << path[i];
            if (i < path.size() - 1) {
                pathStream << " -> ";
            }
        }
        return pathStream.str();
    }


    if (dist[static_cast<size_t>(end)] == numeric_limits<int>::max()) {
        return "No path exists";
    } else {
        stringstream pathStream;
        vector<size_t> path;
        for (size_t at = static_cast<size_t>(end); at != -1; at = prev[at]) {
            path.push_back(at);
        }
        reverse(path.begin(), path.end());

        for (size_t i = 0; i < path.size(); ++i) {
            pathStream << path[i];
            if (i < path.size() - 1) {
                pathStream << " -> ";
            }
        }
        return pathStream.str();
    }}


string bellmanFord(Graph &g, int start, int end)
{
    size_t numVertices = g.getNumVertices();
    vector<int> distance(numVertices, INT_MAX);
    vector<int> predecessor(numVertices, -1);
    distance[(size_t)start] = 0;

    vector<vector<int>> adjacencyMatrix = g.getAdjMatrix();
    // Running |n-1| iterations
    for (size_t i = 0; i < numVertices - 1; i++)
    {
        for (size_t j = 0; j < numVertices; j++)
        {
            for (size_t k = 0; k < numVertices; k++)
            {
                // Relax
                if (adjacencyMatrix[j][k] != 0 && distance[j] != INT_MAX && distance[j] + adjacencyMatrix[j][k] < distance[k])
                {
                    distance[k] = distance[j] + adjacencyMatrix[j][k];
                    predecessor[k] = j;
                }
            }
        }
    }

    // Check for negative-weight cycles, by doing anthoer iteration and finding if we can still relax the path
    for (size_t j = 0; j < numVertices; j++)
    {
        for (size_t k = 0; k < numVertices; k++)
        {
            if (adjacencyMatrix[j][k] != 0 && distance[j] != INT_MAX && distance[j] + adjacencyMatrix[j][k] < distance[k])
            {
                return "Negative cycle detected, cannot define the shortest path";
            }
        }
    }

    // If there is no path from start to end
    if (distance[(size_t)end] == INT_MAX)
    {
        return "No path exist";
    }

    // Build the shortest path from end to start
    string path = to_string(end);
    for (size_t v = (size_t)end; v != start; v = (size_t)predecessor[v])
    {
        path.insert(0, to_string((int)predecessor[v]) + " -> ");
    }

    return path;
}




bool Algorithms::isConnected(Graph &g)
{
    if (g.getNumVertices() == 0)
    {
        return true;
    }

    if (g.getisDirected()) // if the graph is directed
    {
        return directedIsConnected(g);
    }
    return undirectedIsConnected(g);
}

bool Algorithms::isContainsCycle(Graph &g)
{
    if (g.getisDirected()) 
    {
        return isContainCylcleDirected(g);
    }
    return isContaionCylcleUnDirected(g);
}

// According to what is returning  form graph kind we'll use the correct algorithm
string Algorithms::shortestPath(Graph &g, int start, int end)
{
    if (start < 0 || start >= g.getNumVertices() || end < 0 || end >= g.getNumVertices())
    {
        return "Illegal start & end values, cannot define shortest path";
    }

    if(g.getisWeighted()== false){
        return BFS(g, start, end);
    }else if(g.getisDirected()== false && g.getisWeighted()==true && !negativeCycle(g)){
        return dijkstra(g, start, end);
    }else{
        return bellmanFord(g, start, end);
    }}

// To check negative cycle in a graph we just need to run bellman-ford algorithm and che what have returned
bool Algorithms::negativeCycle(Graph &g){
            size_t numVertices = g.getNumVertices();
            vector<int> dist(numVertices, INT_MAX);
            vector<size_t> parent(numVertices, 0);
            dist[0] = 0; // Starting vertex
            parent[0]= 0;
            size_t cycleVertex= INT_MAX;
            bool flag= false;

            // Relax all edges V-1 times
            for (int i = 0; i < numVertices - 1; i++)
            {
                for (size_t u = 0; u < numVertices; u++)
                {
                    for (size_t v = 0; v < numVertices; v++)
                    {
                        int weight = g.getAdjMatrix()[u][v];
                        if (dist[u] != INT_MAX && (weight != 0 || u == v) && dist[u] + weight < dist[v])
                        {
                            dist[v] = dist[u] + weight;
                            parent[v]= u;
                        }
                    }
                }
            }

            // Check for negative weight cycle
            for (size_t u = 0; u < numVertices; u++)
            {
                for (size_t v = 0; v < numVertices; v++)
                {
                    int weight = g.getAdjMatrix()[u][v];
                    if (dist[u] != INT_MAX && (weight != 0 || u == v) && dist[u] + weight < dist[v])
                    {
                        flag= true;
                        cout<< "Negative cycle detected :" << endl;
                         // Negative cycle detected, backtrack to find a vertex in the cycle
                        cycleVertex = v;
                        for (int j = 0; j < numVertices; j++) {
                            cycleVertex = parent[cycleVertex]; // Backtrack
                        }
                        std::string cycleVertices = std::to_string(cycleVertex);
                        size_t nextVertex = parent[cycleVertex];
                        while (nextVertex != cycleVertex) {
                            cycleVertices += " -> " + std::to_string(nextVertex);
                            nextVertex = parent[nextVertex];
                        }
                          cycleVertices += " -> " + std::to_string(cycleVertex);
                          cout<< cycleVertices <<endl;
                          return flag;
                    }
                }
            }

            cout<< "No negative cycle detected."<<endl;
            return flag;
        }
    


bool Algorithms::isBipartite(Graph &g) {
    // Get the number of vertices in the graph
    size_t numVertices = g.getNumVertices();
    
    // Initialize a color vector with all vertices marked as unvisited (value 4)
    std::vector<size_t> color(numVertices, 2);
    
    // Vector to store vertices belonging to each of the two groups
    std::vector<std::vector<size_t>> groups(2);
    
    // Boolean flag to check if the graph is bipartite
    bool isBipartite = true;

    // Iterate through all vertices
    for (size_t i = 0; i < numVertices; ++i) {
        // If the vertex is unvisited
        if (color[i] == 2) {
            // Initialize a queue for BFS and push the current vertex
            std::queue<int> q;
            q.push(i);
            
            // Assign the starting vertex the color 0
            color[i] = 0;

            // Perform BFS
            while (!q.empty() && isBipartite) {
                size_t u = static_cast<size_t>(q.front());
                q.pop();
                
                // Add the current vertex to its corresponding group
                groups[color[u]].push_back(u);

                // Check all vertices for adjacency
                for (size_t v = 0; v < numVertices; v++) {
                    // Check if there's an edge between vertex u and vertex v
                    if (g.getAdjMatrix()[u][v] != 0) {
                        // If the neighbor vertex v is unvisited
                        if (color[v] == 2) {
                            // Assign the opposite color to the neighbor
                            color[v] = 1 - color[u];
                            // Push the neighbor to the queue for further processing
                            q.push(v);
                        } else if (color[v] == color[u]) {
                            // If the neighbor has the same color as the current vertex, the graph is not bipartite
                            isBipartite = false;
                            break;
                        }
                    }
                }
            }
        }
    }

    // Output the result
    if (!isBipartite) {
        std::cout << "Graph is not bipartite." << std::endl;
    } else {
        std::cout << "Graph is bipartite." << std::endl;
        std::cout << "Group 0: ";
        for (size_t v : groups[0]) {
            std::cout << v << " ";
        }
        std::cout << std::endl;
        std::cout << "Group 1: ";
        for (size_t v : groups[1]) {
            std::cout << v << " ";
        }
        std::cout << std::endl;
    }

    return isBipartite;
}