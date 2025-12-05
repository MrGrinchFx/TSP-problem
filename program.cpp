#include <chrono>
#include <algorithm>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <string>
#include <tuple>
#include <vector>

// max double value
const double INF = std::numeric_limits<double>::max();

struct GraphData {
  std::vector<std::vector<double>> distMatrix; // The lookup table for distances
  std::map<int, int> idToIndex; // Map FileID -> Array Index (0..N)
  std::vector<int> indexToId;   // Map Array Index -> FileID (for output)
  int numNodes;
};

GraphData loadGraphFromFile(const std::string &filename) {
  std::ifstream file(filename);
  GraphData graph;

  if (!file.is_open()) {
    std::cerr << "Error: Could not open file " << filename << std::endl;
    exit(1);
  }

  std::string line;

  std::getline(file, line);
  std::getline(file, line);

  std::vector<std::tuple<int, int, double>> edges;
  std::map<int, int> uniqueNodes;
  int nodeCounter = 0;

  int u_id, v_id;
  double weight;

  while (file >> u_id >> v_id >> weight) {
    if (uniqueNodes.find(u_id) == uniqueNodes.end()) {
      uniqueNodes[u_id] = nodeCounter++;
    }
    if (uniqueNodes.find(v_id) == uniqueNodes.end()) {
      uniqueNodes[v_id] = nodeCounter++;
    }
    edges.emplace_back(u_id, v_id, weight);
  }

  graph.numNodes = uniqueNodes.size();
  graph.idToIndex = uniqueNodes;
  graph.indexToId.resize(graph.numNodes);

  // create reverse map for printing later
  for (auto const &[id, index] : uniqueNodes) {
    graph.indexToId[index] = id;
  }

  // initialize with infinity
  graph.distMatrix.resize(graph.numNodes,
                          std::vector<double>(graph.numNodes, INF));

  for (const auto &edge : edges) {
    int u = uniqueNodes[std::get<0>(edge)];
    int v = uniqueNodes[std::get<1>(edge)];
    double w = std::get<2>(edge);

    // TSP Graphs are undirected (symmetric)
    graph.distMatrix[u][v] = w;
    graph.distMatrix[v][u] = w;
  }

  // Distance to self is 0
  for (int i = 0; i < graph.numNodes; i++) {
    graph.distMatrix[i][i] = 0;
  }

  return graph;
}

// calculate path length using the matrix
double getPathCost(const std::vector<int> &path,
                   const std::vector<std::vector<double>> &matrix) {
  double total = 0;
  for (size_t i = 0; i < path.size(); ++i) {
    int u = path[i];
    int v = path[(i + 1) % path.size()]; // wrap around
    if (matrix[u][v] == INF) {
      return INF; // invalid path check
    }
    total += matrix[u][v];
  }
  return total;
}

// intial path construction
void nearestNeighbor(const GraphData &graph, int startNode,
                     std::vector<int> &path) {
  int n = graph.numNodes;
  path.clear();
  path.reserve(n);
  std::vector<bool> visited(n, false);

  int current = startNode;
  path.push_back(current);
  visited[current] = true;

  for (int i = 0; i < n - 1; ++i) {
    int bestNext = -1;
    double minDist = INF;

    const std::vector<double> &row = graph.distMatrix[current];

    for (int next = 0; next < n; ++next) {
      if (!visited[next]) {
        double d = row[next];
        if (d < minDist) {
          minDist = d;
          bestNext = next;
        }
      }
    }

    if (bestNext != -1) {
      current = bestNext;
      path.push_back(current);
      visited[current] = true;
    } else {
      break;
    }
  }
}

// post optimization
void twoOpt(std::vector<int> &path, const GraphData &graph,
            long long &evalCount) {
  int n = path.size();
  bool improved = true;

  // true implies node u's edges weren't improved last time,
  // so we skip checking u until one of its neighbors changes.
  std::vector<bool> dontLook(graph.numNodes, false);

  while (improved) {
    improved = false;
    for (int i = 0; i < n - 1; ++i) {
      int u = path[i];

      // if this node didn't improve last time, skip it
      if (dontLook[u])
        continue;

      bool improvement_found_for_u = false;

      for (int j = i + 1; j < n; ++j) {
        if (j == n - 1 && i == 0)
          continue; // skip wrap-around edge

        int v = path[(i + 1) % n];
        int x = path[j];
        int y = path[(j + 1) % n];

        evalCount++;

        double dist_uv = graph.distMatrix[u][v];
        double dist_xy = graph.distMatrix[x][y];
        double dist_ux = graph.distMatrix[u][x];
        double dist_vy = graph.distMatrix[v][y];

        // optimization check
        if (dist_ux + dist_vy < dist_uv + dist_xy) {
          std::reverse(path.begin() + i + 1, path.begin() + j + 1);
          improved = true;
          improvement_found_for_u = true;

          // reset bits for the nodes involved in the swap
          // their edges changed so we must check them again in future
          dontLook[u] = false;
          dontLook[v] = false;
          dontLook[x] = false;
          dontLook[y] = false;
        }
      }

      // if we checked all possible swaps for u and found nothing,
      // mark it so we don't check it again until a neighbor changes.
      if (!improvement_found_for_u) {
        dontLook[u] = true;
      }
    }
  }
}

std::tuple<GraphData, double, std::vector<int>> solve(std::string filename) {
  GraphData graph = loadGraphFromFile(filename);
  int N = graph.numNodes;

  double bestCost = INF;
  std::vector<int> bestTour;
  long long evalCount = 0;

  std::vector<int> currentPath;
  currentPath.reserve(N);

  auto start = std::chrono::high_resolution_clock::now();

  for (int startNodeIndex = 0; startNodeIndex < N; ++startNodeIndex) {

    nearestNeighbor(graph, startNodeIndex, currentPath);

    twoOpt(currentPath, graph, evalCount);

    double currentCost = getPathCost(currentPath, graph.distMatrix);

    if (currentCost < bestCost) {
      bestCost = currentCost;
      bestTour = currentPath; // copy only if better
    }
  }

  auto end = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> elapsed = end - start;

  std::cerr << "Evaluated Cycles: " << std::scientific << (double)evalCount
            << "\n";
  std::cerr << "Final Cost: " << std::fixed << std::setprecision(2) << bestCost
            << "\n";
  std::cerr << "Total Time: " << elapsed.count() << "s\n";

  return {graph, bestCost, bestTour};
}

int main(int argc, char **argv) {

  std::string euclideanFile = "TSP_1000_euclidianDistance.txt";
  std::string randomFile = "TSP_1000_randomDistance.txt";

  auto [euclidGraph, euclidCost, euclidTour] = solve(euclideanFile);

  auto [randomGraph, randomCost, randomTour] = solve(randomFile);

  std::string outputFilename = "solution_919858479.txt";
  std::ofstream outFile(outputFilename);

  if (!outFile.is_open()) {
    std::cerr << "Could not open " << outputFilename << " for writing.\n";
    return 1;
  }

  for (size_t i = 0; i < euclidTour.size(); ++i) {
    outFile << euclidGraph.indexToId[euclidTour[i]] << ", ";
  }

  outFile << euclidGraph.indexToId[euclidTour[0]] << "\n";

  for (size_t i = 0; i < randomTour.size(); ++i) {
    outFile << randomGraph.indexToId[randomTour[i]] << ", ";
  }
  outFile << randomGraph.indexToId[randomTour[0]];

  outFile.close();

  return 0;
}
