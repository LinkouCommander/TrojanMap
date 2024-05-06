# EE538 Final Project - TrojanMap Final Report
**Chih-Cheng Hsieh**\
**Hung-Chih Lin**

video presentation link: https://www.youtube.com/watch?v=3IhQXKd003I

## Item 1: Autocomplete

#### Algorithm description: 
This code first converts the input string to lowercase, then iterates over each nodes in the map and converts the node's name to lowercase, using case-insensitive comparison to find matching location names, and adds them to the result list to return.
 
#### Time complexity: O(M * N)
m: length of string `name`\
n: length of node list `data`

#### Runtime: 1 ms

#### Discussion: 
The use of case-insensitive comparison makes the code more user-friendly. The `transform()` function easily implements the conversion between uppercase and lowercase, enabling seamless case conversion.

<!----------------------------------------------------------------------------------------------------------->
## Item 2: Find the location
Including Item 2-1: Find the place's coordinates in the Map (Phase 1) & Item 2-2: Check Edit Distance Between Two Location Names (Phase 2)

#### Algorithm description: 
* In Item 2-1, we first call `GetID()` with string `name` as input to find its id, then we call `GetLat()` & `GetLon()` with id as input to find latitude and longitude of target.
* In Item 2-2, we calculates the edit distance between two strings `a` and `b` using dynamic programming. It initializes a 2D vector `dp` to store the distances between substrings of `a` and `b`. Then, it iterates over through `dp`, filling it based on the minimum number of operations required to transform one substring into another.

#### Time complexity: O(M * N)
Item 2-1: O(N)\
Item 2-2: O(M * N)

#### Runtime: 3 ms

#### Discussion: 
* Item 2-1: It's important to note that we cannot directly find the latitude and longitude using the target's name because the key in `data` is id. Therefore, we must first find the id corresponding to this name.
* Item 2-2: For me, this simple dynamic programming algorithm is a great opportunity for practice.

<!----------------------------------------------------------------------------------------------------------->
## Item 3: Find all location categories

#### Algorithm description: 
Function `GetAllCategories()` iterates through each node in `data`, retrieves the `attribute_list` of each node, and adds unique attributes to a vector `vec`.

#### Time complexity: O(M * N)

#### Runtime: --

#### Discussion: --

<!----------------------------------------------------------------------------------------------------------->
## Item 4: Get all locations of a category

#### Algorithm description: 
Function `GetAllLocationsFromCategory()` iterates through each node in the `data` and checks if the specified category exists in the attributes of each location. If the category exists, the location ID is added to the result vector. 

#### Time complexity: O(M * N)

#### Runtime: --

#### Discussion: --

<!----------------------------------------------------------------------------------------------------------->
## Item 5: Get location matching regular expression

#### Algorithm description: 
Function `GetLocationRegex()` call function `regex_match` to ensure that the input regular expression is correct.

#### Time complexity: O(N)

#### Runtime: --

#### Discussion: 
I've never encountered regular expressions before and not quite sure what they're used for.

<!----------------------------------------------------------------------------------------------------------->
## Item 6: CalculateShortestPath

#### Algorithm description:
* **Dijkstra:** The function `CalculateShortestPath_Dijkstra()` implements the Dijkstra algorithm to find the shortest path between two locations. We initialize a priority queue `pq` to store the nodes to be visited, and unlike BFS and DFS, we prioritize visiting nodes with the shortest distance in `pq`. We also create an unordered_map `visited` to keep track of visited nodes and another unordered_map `node_info` to store every node's previous visited node. We start from the `location1_name` node, select the element on the top of `pq` for visiting, and if we reach the `location2_name` node during the process, we break out of the while loop. We then use `node_info` to find all nodes between `location1_name` and `location2_name` and store them for return.
* **Bellman Ford:** The function `CalculateShortestPath_Bellman_Ford()` implements the Bellman-Ford algorithm to find the shortest path between two locations. It initializes a queue `q` to store the nodes to be visited and creates two unordered_maps: `revisit` to keep track of visited nodes and `node_info` to store each node's shortest path from the starting node and the previous visited node. Starting from the location1_name node, it selects the top element from `q` for visiting. Whenever a node is processed, it updates the shortest path to its neighbor nodes if either the neighbor node is not in `node_info` or the newly calculated shortest path (sum of current visit node's shortest path in `node_info` and the distance between current node and neighbor node) is shorter than the neighbor node's shortest path stored in `node_info`. If the algorithm reaches the `location2_name` node during the process, it breaks out of the while loop. Finally, it uses the information stored in `node_info` to find all nodes between `location1_name` and `location2_name` and stores them for return.

#### Time complexity: 
* **Dijkstra:** O(E+VlogV)\
  We use priority queue (heap) to save node.\
  V: number of vertices\
  E: number of edges
* **Bellman Ford:** O(M * N)

#### Runtime:
| Point A to Point B | Dijkstra | Bellman Ford |
|--------------------|----------|------------- |
| Target to Ralphs | 20 ms | 937 ms |
| Ralphs to Honda | 37 ms | 462 ms |
| Honda to Cava | 27 ms | 858 ms |
| Cava to 66 | 28 ms | 670 ms |

#### Discussion:
If we only look at time complexity, Dijkstra's algorithm with heap will always outperform Bellman-Ford. The experimental results, where Dijkstra significantly outperforms Bellman-Ford, also corroborate this point.

<!----------------------------------------------------------------------------------------------------------->
## Item 7: Cycle Detection

#### Algorithm description: 
The function `Cycle Detection()` implements a BFS algorithm to detect cycles in a graph. Each time we visit a node, we add it to an unordered_map `checked` with its previous node as the value. While examining the neighbor nodes, if a neighbor is already present in `checked` and its value in `checked` is not the current visited node, we conclude that the graph contains a cycle and return true. If the code is not interrupted by this condition, we conclude that there is no cycle and return false.

#### Time complexity: O(M + N)
M: number of vertices\
N: number of edges

#### Runtime: 0 ms

#### Discussion:
We noticed that directly determining cycle based on whether a node is visited or not would include the parent node, resulting in the program falsely detecting a cycle. This could lead to incorrect judgments.

<!----------------------------------------------------------------------------------------------------------->
## Item 8: Topological Sort

#### Algorithm description: 
Function `DeliveringTrojan()` implements a algorithm to solve dependency problems. This algorithm used is trying to calculate the total number of prior nodes for each node, including prior nodes' prior nodes. It then prioritizes nodes without dependency restrictions and fills them in sequentially. First, we create an unordered map `dependencies_map` with value as vector type, where nodes in value prior to key. Then, we create a 2D vector `dependencies_vec` to store these vectors and sort them based on their sizes from shortest to longest. Next, we iterate over `dependencies_vec` and store elements into a vector `result` and return it.

#### Time complexity: O(M * N)
M: number of vectors in dependencies_vec\
N: number of elements in vectors

#### Runtime:
| Number of nodes | Runtime |
|--------------------|----------|
| 3 nodes | 0 ms |
| 5 nodes | 0 ms |
| 10 nodes | 0 ms |

#### Discussion:
This algorithm is not a typical solution using DFS, but it still achieves a runtime of 0 ms within 10 nodes.

<!----------------------------------------------------------------------------------------------------------->
## Item 9: Traveling salesman problem

#### Algorithm description: 
* **Brute-force:** Function `TravelingTrojan_Brute_force()` implements a brute-force algorithm to solve the Traveling Salesman Problem (TSP). It initializes necessary data structures such as the vector `temp` to store the current permutation, integer `dist` to store the shortest path for the current permutation, and vector `records` to store the globally optimal permutation along with its shortest path. Then, it calls `generatePermutations_Backtracking` to do recursion task. Each time `generatePermutations_Backtracking` is called, it first checks if `temp` contains all locations. If so, it updates the shortest path, then backtracks out of the recursion. If not, it iterates through the elements in the input vector `location_ids`, updating `dist`, pushing the element into `temp`, and calling the next level of `generatePermutations_Backtracking`.
* **Brute-force with backtracking:** Function TravelingTrojan_Backtracking() is basically identical to `TravelingTrojan_Brute_force()`, but it includes a new backtracking condition. If we discover that `temp` hasn't included all locations and its total distance `dist` has already exceeded the globally optimal solution, then we terminate the recursion prematurely.
* **2-opt:** Function `TravelingTrojan_2opt()` implements the 2-opt algorithm to solve the traveling salesman problem (TSP). Firstly, we copy the given location IDs into `int_ids` and add the starting point to the end to form a cyclic path. Then, we calculate the initial length of this cyclic path and set it as the first element of `records`. Next, we use the 2-opt algorithm to optimize this cyclic path, continuously searching for ways to reduce the total path length by swapping two adjacent subpaths. After each swap, we update the first element of records to reflect the new total path length and add the new path to the second element of `records`. The algorithm terminates when no further improvement is found, meaning no shorter path is discovered, and returns the optimal solution.

#### Time complexity:
* **Brute-force:** O(N!)
* **Brute-force with backtracking:** O(N!)
* **2-opt:** O(N^2)

#### Runtime:
| Number of nodes | Brute-force | Backtracking | 2-opt |
|-----------------|-------------|--------------|-------|
| 5 nodes | 0 ms | 0 ms | 0 ms |
| 7 nodes | 35 ms | 18 ms | 0 ms |
| 9 nodes | 3334 ms | 447 ms | 0 ms |

#### Discussion:
Although both brute-force algorithms with and without backtracking have a time complexity of O(N!), the performance of the brute-force algorithm with backtracking is far superior to that without it, as observed from the results. Moreover, the 2-opt algorithm consistently outperforms the previous two, despite the very small chance of not finding the optimal solution or getting stuck in an infinite loop.

<!----------------------------------------------------------------------------------------------------------->
## Item 10: Find Nearby

#### Algorithm description: 
Function `FindNearby` iterates over data and checks if each node meets our requirements. For each node examined, we first check if the node has the attribute `attributesName`, then we call `CalculateDistance()` to calculate the distance between the node being examined and the origin node `name`, and determine if this distance is less than `r`. If these conditions are all satisfied, we add it into vector `vec`.

#### Time complexity: O(nlogn)
The time complexity of traversing each node in the `data` is O(n), and the time complexity of sorting is O(nlogn). Therefore, the total time complexity is O(N + NlogN), which is O(NlogN), where N is the total number of nodes.

#### Runtime: 4 ms

#### Discussion:
Due to the requirement that the order of locations should be from nearest to farthest, to achieve efficient sorting, our `vec` also stores the distance, and we call `sort()` to implement quicksort.

<!----------------------------------------------------------------------------------------------------------->
## Item 11: Find Path to Visit All Places

#### Algorithm description: 
Function `TrojanPath` implements Dijkstra and Brute-force algorithms. Essentially, we use Brute-force as the framework for this algorithm, replacing the function that calculates the distance between two points in Brute-force with `CalculateShortestPath_Dijkstra()`. Other parts of the algorithm remain almost identical to Brute-force.

#### Time complexity: O(N! * (E+VlogV))
N: number of location\
V: number of vertices\
E: number of edges

#### Runtime:
| Number of nodes | Runtime |
|--------------------|----------|
| 3 nodes | 308 ms |
| 5 nodes | 8401 ms |
| 6 nodes | 54209 ms |

#### Discussion:
The primary reason for the inefficiency of this algorithm comes from brute-force algorithm, with a time complexity of O(N!). Consequently, if the number of locations to be calculated exceeds 6, the code's runtime becomes excessively long. We guess that replacing the brute-force algorithm with dynamic programming approach might lead to better performance.

<!----------------------------------------------------------------------------------------------------------->
## Item 12: Check Exist of Path with Constrain

#### Algorithm description: 
Function `Queries()` implements the DFS algorithm to solve the gas tank problem. The only change we make in `Queries()` is that we add an additional condition while pushing neighbors into the stack, which is whether the distance is shorter than the size of the gas tank.

#### Time complexity: O(m+n)
m: number of vertices\
n: number of edges

#### Runtime: 11 ms (result of example in README.md)

#### Discussion:
This algorithm is a simple variant of DFS.
