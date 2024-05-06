#include "trojanmap.h"
// using namespace std;

//-----------------------------------------------------
// TODO: Students should implement the following:
//-----------------------------------------------------
/**
 * GetLat: Get the latitude of a Node given its id. If id does not exist, return
 * -1.
 *
 * @param  {std::string} id : location id
 * @return {double}         : latitude
 */
double TrojanMap::GetLat(const std::string &id) {
  auto it = data.find(id);
  if(it != data.end()) return it->second.lat;

  return -1;
}

/**
 * GetLon: Get the longitude of a Node given its id. If id does not exist,
 * return -1.
 *
 * @param  {std::string} id : location id
 * @return {double}         : longitude
 */
double TrojanMap::GetLon(const std::string &id) {
  auto it = data.find(id);
  if(it != data.end()) return it->second.lon;

  return -1;
}

/**
 * GetName: Get the name of a Node given its id. If id does not exist, return
 * "NULL".
 *
 * @param  {std::string} id : location id
 * @return {std::string}    : name
 */
std::string TrojanMap::GetName(const std::string &id) {
  auto it = data.find(id);
  if(it != data.end()) return it->second.name;

  return "NULL";
}

/**
 * GetNeighborIDs: Get the neighbor ids of a Node. If id does not exist, return
 * an empty vector.
 *
 * @param  {std::string} id            : location id
 * @return {std::vector<std::string>}  : neighbor ids
 */
std::vector<std::string> TrojanMap::GetNeighborIDs(const std::string &id) {
  auto it = data.find(id);
  if(it != data.end()) return it->second.neighbors;

  return {};
}

/**
 * GetID: Given a location name, return the id.
 * If the node does not exist, return an empty string.
 * The location name must be unique, which means there is only one node with the name.
 *
 * @param  {std::string} name          : location name
 * @return {std::string}               : id
 */
std::string TrojanMap::GetID(const std::string &name) {
  std::string res = "";

  string lname = name;
  transform(name.begin(), name.end(), lname.begin(), [](unsigned char c){ return std::tolower(c); });

  for(const auto& node : data) {
    string n = node.second.name;
    if(n.size() == 0) continue;
    string ln = n;

    transform(n.begin(), n.end(), ln.begin(), [](unsigned char c){ return std::tolower(c); });

    if(lname == ln) return node.first;
  }

  return res;
}

/**
 * GetPosition: Given a location name, return the position. If id does not
 * exist, return (-1, -1).
 *
 * @param  {std::string} name          : location name
 * @return {std::pair<double,double>}  : (lat, lon)
 */
std::pair<double, double> TrojanMap::GetPosition(std::string name) {
  std::pair<double, double> results(-1, -1);

  string id = GetID(name);
  if(id.size() == 0) return results;

  results.first = GetLat(id);
  results.second = GetLon(id);

  return results;
}

/**
 * CalculateEditDistance: Calculate edit distance between two location names
 * @param  {std::string} a          : first string
 * @param  {std::string} b          : second string
 * @return {int}                    : edit distance between two strings
 */
// int CED_dp(string a, string b, size_t i, vector<vector<int>>& dp) {
//   if(i >= a.size()) return b.size() - i;
//   if(i >= b.size()) return a.size() - i;

//   if(dp[i][b.size()] != -1) dp[i][b.size()];

//   if(a[i] == b[i]) return CED_dp(a, b, i + 1, dp);

//   string inss = b.insert(i, 1, a[i]);
//   string dels = b.erase(i, 1);
//   string reps = b;
//   b[i] = a[i];

//   int ins = CED_dp(a, inss, i + 1, dp) + 1;
//   int del = CED_dp(a, dels, i, dp) + 1;
//   int rep = CED_dp(a, reps, i + 1, dp) + 1;

//   return dp[i][b.size()] = min(ins, min(del, rep));

// }
int TrojanMap::CalculateEditDistance(std::string a, std::string b) {
  int m = a.size();
  int n = b.size();

  vector<vector<int>> dp(m + 1, vector<int>(n + 1, 0));

  for(int i = 0; i <= m; i++) {
    for(int j = 0; j <= n; j++) {
      if(i == 0) dp[i][j] = j;
      else if(j == 0) dp[i][j] = i;
      else if(a[i - 1] == b[j - 1]) dp[i][j] = dp[i - 1][j - 1];
      else dp[i][j] = min(dp[i - 1][j - 1] + 1, min(dp[i - 1][j] + 1, dp[i][j - 1] + 1));
    }
  }

  return dp[m][n];
}


/**
 * FindClosestName: Given a location name, return the name with the smallest edit
 * distance.
 *
 * @param  {std::string} name          : location name
 * @return {std::string} tmp           : the closest name
 */
std::string TrojanMap::FindClosestName(std::string name) {
  std::string tmp = ""; // Start with a dummy word
  int mini = INT_MAX;

  for(const auto& node : data) {
    if(node.second.name.size() == 0) continue;
    int minDis = CalculateEditDistance(name, node.second.name);
    // cout << node.second.name << ": " << minDis << "\n";
    if(minDis < mini) {
      mini = minDis;
      tmp = node.second.name;
    }
  }
  return tmp;
}

/**
 * Autocomplete: Given a parital name return all the possible locations with
 * partial name as the prefix. The function should be case-insensitive.
 *GetName(node.first)
 * @param  {std::string} name          : partial name
 * @return {std::vector<std::string>}  : a vector of full names
 */
std::vector<std::string> TrojanMap::Autocomplete(std::string name) {
  std::vector<std::string> results;
  int l = name.size();

  if(l == 0) return {};

  string lname = name;
  transform(name.begin(), name.end(), lname.begin(), [](unsigned char c){ return std::tolower(c); });

  for(const auto& node : data) {
    string n = node.second.name;
    if(n.size() == 0) continue;
    string ln = n;

    transform(n.begin(), n.end(), ln.begin(), [](unsigned char c){ return std::tolower(c); });
    if(lname == ln.substr(0, l)) results.push_back(n);
  }

  return results;
}

/**
 * GetAllCategories: Return all the possible unique location categories, i.e.
 * there should be no duplicates in the output.
 *
 * @return {std::vector<std::string>}  : all unique location categories
 */
std::vector<std::string> TrojanMap::GetAllCategories() {
  vector<string> vec;
  for(const auto& node : data) {
    auto attribute_list = node.second.attributes;
    if(attribute_list.size() == 0) continue;

    for(auto& att : attribute_list) {
      auto it = find(vec.begin(), vec.end(), att);
      if(it == vec.end()) vec.push_back(att);
    }
  }
  return vec;
}

/**
 * GetAllLocationsFromCategory: Return all the locations of the input category (i.e.
 * 'attributes' in data.csv). If there is no location of that category, return
 * (-1, -1). The function should be case-insensitive.
 *
 * @param  {std::string} category         : category name (attribute)
 * @return {std::vector<std::string>}     : ids
 */
std::vector<std::string> TrojanMap::GetAllLocationsFromCategory(
    std::string category) {
  std::vector<std::string> res;
  for(const auto& node : data) {
    auto attribute_list = node.second.attributes;
    if(attribute_list.size() == 0) continue;

    auto it = attribute_list.find(category);
    if(it != attribute_list.end()) res.push_back(node.first);
  }
  return res;
}

/**
 * GetLocationRegex: Given the regular expression of a location's name, your
 * program should first check whether the regular expression is valid, and if so
 * it returns all locations that match that regular expression.
 *
 * @param  {std::regex} location name      : the regular expression of location
 * names
 * @return {std::vector<std::string>}     : ids
 */
std::vector<std::string> TrojanMap::GetLocationRegex(std::regex location) {
  vector<string> vec = {};

  for(const auto& node : data) {
    string n = node.second.name;
    if(n.size() == 0) continue;
    if(regex_match(n, location)) vec.push_back(node.first);
  }
  return vec;
}

/**
 * CalculateDistance: Get the distance between 2 nodes.
 * We have provided the code for you. Please do not need to change this function.
 * You can use this function to calculate the distance between 2 nodes.
 * The distance is in mile.
 * The distance is calculated using the Haversine formula.
 * https://en.wikipedia.org/wiki/Haversine_formula
 * 
 * @param  {std::string} a  : a_id
 * @param  {std::string} b  : b_id
 * @return {double}  : distance in mile
 */
double TrojanMap::CalculateDistance(const std::string &a_id,
                                    const std::string &b_id) {
  // Do not change this function
  Node a = data[a_id];
  Node b = data[b_id];
  double dlon = (b.lon - a.lon) * M_PI / 180.0;
  double dlat = (b.lat - a.lat) * M_PI / 180.0;
  double p = pow(sin(dlat / 2), 2.0) + cos(a.lat * M_PI / 180.0) *
                                           cos(b.lat * M_PI / 180.0) *
                                           pow(sin(dlon / 2), 2.0);
  double c = 2 * asin(std::min(1.0, sqrt(p)));
  return c * 3961;
}

/**
 * CalculatePathLength: Calculates the total path length for the locations
 * inside the vector.
 * We have provided the code for you. Please do not need to change this function.
 * 
 * @param  {std::vector<std::string>} path : path
 * @return {double}                        : path length
 */
double TrojanMap::CalculatePathLength(const std::vector<std::string> &path) {
  // Do not change this function
  double sum = 0;
  for (int i = 0; i < int(path.size()) - 1; i++) {
    sum += CalculateDistance(path[i], path[i + 1]);
  }
  return sum;
}

/**
 * CalculateShortestPath_Dijkstra: Given 2 locations, return the shortest path
 * which is a list of id. Hint: Use priority queue.
 *
 * @param  {std::string} location1_name     : start
 * @param  {std::string} location2_name     : goal
 * @return {std::vector<std::string>}       : path
 */
std::vector<std::string> TrojanMap::CalculateShortestPath_Dijkstra(
    std::string location1_name, std::string location2_name) {
  std::vector<std::string> path;
  if(location1_name == "" || location2_name == "") return path;

  typedef pair<double, string> idp;
  priority_queue<idp, vector<idp>, greater<idp>> pq;
  unordered_map<string, bool> visited;
  unordered_map<string, pair<double, string>> node_info;

  string id1 = GetID(location1_name);
  string id2 = GetID(location2_name);

  pq.push(make_pair(0, id1));
  node_info[id1] = make_pair(0, "-1");

  while(!pq.empty()) {
    string tp = pq.top().second;
    double cur_dist = node_info[tp].first;
    pq.pop();

    if(tp == id2) {
      visited[id2] = true;
      break;
    }

    if(visited[tp]) continue;

    visited[tp] = true;
    for(const auto& ele : GetNeighborIDs(tp)) {
      double dist_to_nb = CalculateDistance(tp, ele);
      double total_dist = dist_to_nb + cur_dist;

      if(node_info.count(ele) == 0 || total_dist < node_info[ele].first) {
        node_info[ele] = make_pair(total_dist, tp);
        pq.push(make_pair(total_dist, ele));
      }
    }
  }

  if(!visited[id2]) return path;

  string nd = id2;
  while(nd != "-1") {
    path.push_back(nd);
    nd = node_info[nd].second;
  }
  reverse(path.begin(), path.end());

  return path;
}

/**
 * CalculateShortestPath_Bellman_Ford: Given 2 locations, return the shortest
 * path which is a list of id. Hint: Do the early termination when there is no
 * change on distance.
 *
 * @param  {std::string} location1_name     : start
 * @param  {std::string} location2_name     : goal
 * @return {std::vector<std::string>}       : path
 */
std::vector<std::string> TrojanMap::CalculateShortestPath_Bellman_Ford(
    std::string location1_name, std::string location2_name) {
  std::vector<std::string> path;
  if(location1_name == "" || location2_name == "") return path;

  queue<string> q;

  unordered_map<string, bool> revisit;
  unordered_map<string, pair<double, string>> node_info;

  string id1 = GetID(location1_name);
  string id2 = GetID(location2_name);

  q.push(id1);
  node_info[id1] = {0, "-1"};

  while(!q.empty()) {
    string cur = q.front();
    double cur_dist = node_info[cur].first;
    q.pop();

    revisit[cur] = false;
    for(const auto& ele : GetNeighborIDs(cur)) {
      double dist_to_nb = CalculateDistance(cur, ele);
      double total_dist = dist_to_nb + cur_dist;

      if(node_info.count(ele) == 0 || total_dist < node_info[ele].first) {
        node_info[ele] = make_pair(total_dist, cur);
        if(!revisit[ele]) {
          revisit[ele] = true;
          q.push(ele);
        }
        
      }
    }
  }

  string nd = id2;
  while(nd != "-1") {
    path.push_back(nd);
    nd = node_info[nd].second;
  }
  reverse(path.begin(), path.end());

  return path;
}

/**
 * Traveling salesman problem: Given a list of locations, return the shortest
 * path which visit all the places and back to the start point.
 *
 * @param  {std::vector<std::string>} input : a list of locations needs to visit
 * @return {std::pair<double, std::vector<std::vector<std::string>>} : a pair of total distance and the all the progress to get final path, 
 *                                                                      for example: {10.3, {{0, 1, 2, 3, 4, 0}, {0, 1, 2, 3, 4, 0}, {0, 4, 3, 2, 1, 0}}},
 *                                                                      where 10.3 is the total distance, 
 *                                                                      and the first vector is the path from 0 and travse all the nodes and back to 0,
 *                                                                      and the second vector is the path shorter than the first one,
 *                                                                      and the last vector is the shortest path.
 */
// Please use brute force to implement this function, ie. find all the permutations.
void  TrojanMap::generatePermutations(vector<string>& locations, vector<string> temp, pair<double, vector<vector<string>>>& records) {
  if(temp.size() == locations.size()) {
    temp.push_back(temp[0]);
    double temp_dist = CalculatePathLength(temp);
    // cout << "temp_dist: " << temp_dist << "\n";
    if(temp_dist < records.first) {
      records.first = temp_dist;
      records.second.push_back(temp);
      // cout << "dist: " << records.first << "\n";
    }
    return;
  }

  for(auto e : locations) {
    
    if(find(temp.begin(), temp.end(), e) != temp.end()) continue;
    // cout << e << " ";
    vector<string> new_temp = temp;
    new_temp.push_back(e);
    generatePermutations(locations, new_temp, records);
    // cout << "\n";
  }
}
std::pair<double, std::vector<std::vector<std::string>>> TrojanMap::TravelingTrojan_Brute_force(
                                    std::vector<std::string> location_ids) {
  std::pair<double, std::vector<std::vector<std::string>>> records;

  if(location_ids.size() == 0) {
    std::vector<std::vector<std::string>> path;
    return make_pair(0, path);
  }

  vector<string> temp;

  vector<string> int_ids =  location_ids;
  int_ids.push_back(int_ids[0]);
  double int_dist = CalculatePathLength(int_ids);

  records.first = int_dist;
  records.second.push_back(int_ids);

  generatePermutations(location_ids, temp, records);

  // cout << vec.size() << " "<< vec[0].size();
  // for(int i = 0; i < records.second.size(); i++) {
  //   for(auto e : records.second[i]) {
  //     cout << e << " ";
  //   }
  //   cout << "\n";
  // }

  return records;
}

// Please use backtracking to implement this function
void  TrojanMap::generatePermutations_Backtracking(vector<string>& locations, vector<string> temp, pair<double, vector<vector<string>>>& records, double dist) {
  if(dist >= records.first) return;

  if(temp.size() == locations.size()) {
    temp.push_back(temp[0]);
    double temp_dist = CalculatePathLength(temp);
    // cout << "temp_dist: " << temp_dist << "\n";
    if(temp_dist < records.first) {
      records.first = temp_dist;
      records.second.push_back(temp);
      // cout << "dist: " << records.first << "\n";
    }
    return;
  }

  for(auto e : locations) {
    if(find(temp.begin(), temp.end(), e) != temp.end()) continue;

    double cur_dist = dist;
    if(!temp.empty()) cur_dist += CalculateDistance(temp.back(), e);

    vector<string> new_temp = temp;
    new_temp.push_back(e);

    generatePermutations_Backtracking(locations, new_temp, records, cur_dist);
  }
}
std::pair<double, std::vector<std::vector<std::string>>> TrojanMap::TravelingTrojan_Backtracking(
                                    std::vector<std::string> location_ids) {
  std::pair<double, std::vector<std::vector<std::string>>> records;

  if(location_ids.size() == 0) {
    std::vector<std::vector<std::string>> path;
    return make_pair(0, path);
  }

  vector<string> temp;

  vector<string> int_ids =  location_ids;
  int_ids.push_back(int_ids[0]);
  double int_dist = CalculatePathLength(int_ids);

  records.first = int_dist;
  records.second.push_back(int_ids);

  generatePermutations_Backtracking(location_ids, temp, records, 0);

  return records;
}

// Hint: https://en.wikipedia.org/wiki/2-opt
void swap_2opt(vector<string>& path, int i, int j) {
  reverse(path.begin()+i, path.begin()+j+1);
  return;
}
std::pair<double, std::vector<std::vector<std::string>>> TrojanMap::TravelingTrojan_2opt(
      std::vector<std::string> location_ids){
  std::pair<double, std::vector<std::vector<std::string>>> records;

  if(location_ids.size() == 0) {
    std::vector<std::vector<std::string>> path;
    return make_pair(0, path);
  }

  vector<string> int_ids =  location_ids;
  int_ids.push_back(int_ids[0]);
  double int_dist = CalculatePathLength(int_ids);
  cout << int_dist;
  
  records.first = int_dist;
  records.second.push_back(int_ids);

  int l = int_ids.size();
  bool foundImprovement = true;
  while(foundImprovement) {
    foundImprovement = false;

    for(int i = 1; i < l - 2; i++) {
      for(int j = i + 1; j < l - 1; j++) {
        double lengthDelta = -CalculateDistance(int_ids[i], int_ids[i - 1]) - CalculateDistance(int_ids[j], int_ids[j + 1])
                         + CalculateDistance(int_ids[i], int_ids[j + 1]) + CalculateDistance(int_ids[i - 1], int_ids[j]);
        // cout << "lengthDelta: " << lengthDelta << "\n";

        if(lengthDelta < 0) {
          swap_2opt(int_ids, i, j);
          records.second.push_back(int_ids);
          records.first += lengthDelta;
          // cout << "records.first: " << records.first << "\n";
          foundImprovement = true;
        }
      }
    }
  }
  return records;
}

// This is optional
std::pair<double, std::vector<std::vector<std::string>>> TrojanMap::TravelingTrojan_3opt(
      std::vector<std::string> location_ids){
  std::pair<double, std::vector<std::vector<std::string>>> records;
  return records;
}

/**
 * Given CSV filename, it read and parse locations data from CSV file,
 * and return locations vector for topological sort problem.
 * We have provided the code for you. Please do not need to change this function.
 * Example: 
 *   Input: "topologicalsort_locations.csv"
 *   File content:
 *    Name
 *    Ralphs
 *    KFC
 *    Chick-fil-A
 *   Output: ['Ralphs', 'KFC', 'Chick-fil-A']
 * @param  {std::string} locations_filename     : locations_filename
 * @return {std::vector<std::string>}           : locations
 */
std::vector<std::string> TrojanMap::ReadLocationsFromCSVFile(
    std::string locations_filename) {
  std::vector<std::string> location_names_from_csv;
  std::fstream fin;
  fin.open(locations_filename, std::ios::in);
  std::string line, word;
  getline(fin, line);
  while (getline(fin, word)) {
    location_names_from_csv.push_back(word);
  }
  fin.close();
  return location_names_from_csv;
}

/**
 * Given CSV filenames, it read and parse dependencise data from CSV file,
 * and return dependencies vector for topological sort problem.
 * We have provided the code for you. Please do not need to change this function.
 * Example: 
 *   Input: "topologicalsort_dependencies.csv"
 *   File content:
 *     Source,Destination
 *     Ralphs,Chick-fil-A
 *     Ralphs,KFC
 *     Chick-fil-A,KFC
 *   Output: [['Ralphs', 'Chick-fil-A'], ['Ralphs', 'KFC'], ['Chick-fil-A', 'KFC']]
 * @param  {std::string} dependencies_filename     : dependencies_filename
 * @return {std::vector<std::vector<std::string>>} : dependencies
 */
std::vector<std::vector<std::string>> TrojanMap::ReadDependenciesFromCSVFile(
    std::string dependencies_filename) {
  std::vector<std::vector<std::string>> dependencies_from_csv;
  std::fstream fin;
  fin.open(dependencies_filename, std::ios::in);
  std::string line, word;
  getline(fin, line);
  while (getline(fin, line)) {
    std::stringstream s(line);
    std::vector<std::string> dependency;
    while (getline(s, word, ',')) {
      dependency.push_back(word);
    }
    dependencies_from_csv.push_back(dependency);
  }
  fin.close();
  return dependencies_from_csv;
}

/**
 * DeliveringTrojan: Given a vector of location names, it should return a
 * sorting of nodes that satisfies the given dependencies. If there is no way to
 * do it, return a empty vector.
 *
 * @param  {std::vector<std::string>} locations                     : locations
 * @param  {std::vector<std::vector<std::string>>} dependencies     : prerequisites
 * @return {std::vector<std::string>} results                       : results
 */
bool compareLength(const vector<string>& vec1, const vector<string>& vec2) {
  return vec1.size() < vec2.size();
}
std::vector<std::string> TrojanMap::DeliveringTrojan(
    std::vector<std::string> &locations,
    std::vector<std::vector<std::string>> &dependencies) {
  std::vector<std::string> result;
  stack<string> stk;
  unordered_map<string, bool> visited;

  unordered_map<string, vector<string>> dependencies_map;
  for(const auto& dep : dependencies) {
    dependencies_map[dep[1]].push_back(dep[0]);
  }

  vector<vector<string>> dependencies_vec;
  for(auto& node : locations) {
    vector<string> vec;
    if(dependencies_map.count(node) == 0) {
      vec.push_back(node);
      dependencies_vec.push_back(vec);
    }
    else {
      vec.push_back(node);
      if(dependencies_map[node].size() > 0) vec.insert(vec.begin() + 1, dependencies_map[node].begin(), dependencies_map[node].end());
      dependencies_vec.push_back(vec);
    }
  }

  // for(auto e : dependencies_vec) {
  //   for(auto v : e) {
  //     cout << v << ", ";
  //   }
  //   cout << "\n";
  // }

  sort(dependencies_vec.begin(), dependencies_vec.end(), compareLength);

  for(auto& ele : dependencies_vec) {
    for(auto& e : ele) {
      if(visited[e] == true) continue;
      visited[e] = true;
      stk.push(e);
    }
  }

  while(!stk.empty()) {
    string s = stk.top();
    stk.pop();
    result.push_back(s);
  }

  reverse(result.begin(), result.end());

  return result;     
}

/**
 * inSquare: Give a id retunr whether it is in square or not.
 *
 * @param  {std::string} id            : location id
 * @param  {std::vector<double>} square: four vertexes of the square area
 * @return {bool}                      : in square or not
 */
bool TrojanMap::inSquare(std::string id, std::vector<double> &square) {
  double lat = GetLat(id);
  double loc = GetLon(id);

  if(loc >= square[0] && loc <= square[1] && lat <= square[2] && lat >= square[3]) return true;

  return false;
}


/**
 * GetSubgraph: Give four vertexes of the square area, return a list of location
 * ids in the squares
 *
 * @param  {std::vector<double>} square         : four vertexes of the square
 * area
 * @return {std::vector<std::string>} subgraph  : list of location ids in the
 * square
 */
std::vector<std::string> TrojanMap::GetSubgraph(std::vector<double> &square) {
  // include all the nodes in subgraph
  std::vector<std::string> subgraph;

  for(const auto& node : data) {
    string id = node.first;
    if(inSquare(id, square)) subgraph.push_back(id);
  }

  return subgraph;
}

/**
 * Cycle Detection: Given four points of the square-shape subgraph, return true
 * if there is a cycle path inside the square, false otherwise.
 *
 * @param {std::vector<std::string>} subgraph: list of location ids in the
 * square
 * @param {std::vector<double>} square: four vertexes of the square area
 * @return {bool}: whether there is a cycle or not
 */
bool TrojanMap::CycleDetection(std::vector<std::string> &subgraph, std::vector<double> &square) {
  // cout << subgraph.size();
  unordered_map<string, string> checked;
  queue<string> q;

  for(auto& location : subgraph) {
    if(checked.count(location) != 0) continue;

    q.push(location);
    checked[location] = "";

    while(!q.empty()) {
      string cur = q.front();
      q.pop();

      // cout << cur << ": ";
      auto neighbor_list = GetNeighborIDs(cur);
      if(neighbor_list.empty()) continue;

      for(const auto& nb : neighbor_list) {
        // cout << nb << ", ";
        if(!inSquare(nb, square)) continue;

        if(checked.count(nb) != 0) {
          if(nb != checked[cur]) return true;
        }
        else {
          checked[nb] = cur;
          q.push(nb);
        }
      }
      // cout << "\n";
    }
  }
  return false;
}

/**
 * FindNearby: Given a class name C, a location name L and a number r,
 * find all locations in class C on the map near L with the range of r and
 * return a vector of string ids
 *
 * @param {std::string} className: the name of the class
 * @param {std::string} locationName: the name of the location
 * @param {double} r: search radius
 * @param {int} k: search numbers
 * @return {std::vector<std::string>}: location name that meets the requirements
 */
std::vector<std::string> TrojanMap::FindNearby(std::string attributesName, std::string name, double r, int k) {
  std::vector<std::string> res;
  vector<pair<double, string>> vec;

  int count = 0;
  for(const auto& node : data) {
    string node_name = node.second.name;
    if(node_name.empty() || node_name == name) continue;

    auto attribute_list = node.second.attributes;
    if(attribute_list.empty() || find(attribute_list.begin(), attribute_list.end(), attributesName) == attribute_list.end()) continue;

    double dist = CalculateDistance(node.first, GetID(name));

    if(dist <= r) {
      vec.push_back({dist, node.first});
      count++;
    }
    if(count == k) break;
  }

  auto compare = [](const std::pair<double, std::string>& a, const std::pair<double, std::string>& b) {
    return a.first < b.first;
  };
  std::sort(vec.begin(), vec.end(), compare);

  for(auto& e : vec) {
    res.push_back(e.second);
  }
  return res;
}

/**
 * Shortest Path to Visit All Nodes: Given a list of locations, return the shortest
 * path which visit all the places and no need to go back to the start point.
 *
 * @param  {std::vector<std::string>} input : a list of locations needs to visit
 * @return {std::vector<std::string> }      : the shortest path
 */
// Please use backtracking to implement this function
void  TrojanMap::generatePermutations_TrojanPath(vector<string>& locations, vector<string> temp, pair<double, vector<string>>& records, double dist) {
  if(dist >= records.first) return;

  if(temp.size() == locations.size()) {
    // cout << "dist: " << dist << "\n";
    if(dist < records.first) {
      records.first = dist;
      records.second = temp;
      // cout << "records.first: " << records.first << "\n";
    }
    return;
  }

  for(auto e : locations) {
    if(find(temp.begin(), temp.end(), e) != temp.end()) continue;

    double cur_dist = dist;
    if(!temp.empty()) {
      auto temp_path = CalculateShortestPath_Dijkstra(temp.back(), e);
      double d = CalculatePathLength(temp_path);
      cur_dist += d;
    }
    vector<string> new_temp = temp;
    new_temp.push_back(e);
    // cout << e << " ";
    generatePermutations_TrojanPath(locations, new_temp, records, cur_dist);
    // cout << "\n";
  }
}
std::vector<std::string> TrojanMap::TrojanPath(
      std::vector<std::string> &location_names) {
  std::vector<std::string> res;

  if(location_names.size() == 0) return res;

  vector<string> temp;
  pair<double, vector<string>> records;

  vector<string> intitial_names =  location_names;
  double initial_dist = 0;
  for(int i = 0; i < intitial_names.size() - 1; i++) {
    auto initial_path = CalculateShortestPath_Dijkstra(intitial_names[i], intitial_names[i + 1]);
    double temp_dist = CalculatePathLength(initial_path);
    initial_dist += temp_dist;
  }

  records.first = initial_dist;
  records.second = intitial_names;

  generatePermutations_TrojanPath(location_names, temp, records, 0);

  for(int i = 0; i < records.second.size() - 1; i++) {
    auto path = CalculateShortestPath_Dijkstra(records.second[i], records.second[i + 1]);
    if(i != records.second.size() - 2)
      path.pop_back();
    res.insert(res.end(), path.begin(), path.end());
  }
  return res;
}

/**
 * Given a vector of queries, find whether there is a path between the two locations with the constraint of the gas tank.
 *
 * @param  {std::vector<std::pair<double, std::vector<std::string>>>} Q : a list of queries
 * @return {std::vector<bool> }      : existence of the path
 */
std::vector<bool> TrojanMap::Queries(const std::vector<std::pair<double, std::vector<std::string>>>& q) {
    std::vector<bool> ans(q.size(), false);

    for(int i = 0; i < q.size(); i++) {
      double fuel = q[i].first;
      string start = q[i].second[0], destination = q[i].second[1];
      if(fuel == 0 || start.empty() || destination.empty()) continue;

      string start_ID = GetID(start), end_ID = GetID(destination);
      if(start_ID.empty() || end_ID.empty()) continue;

      unordered_map<string, bool> visited;
      stack<string> s;

      s.push(start_ID);

      while(!s.empty()) {
        string temp = s.top();
        s.pop();

        if(temp == end_ID) {
          ans[i] = true;
          break;
        }

        if(visited[temp]) continue;

        visited[temp] = true;
        auto neighbor_list = GetNeighborIDs(temp);
        for(auto& nb : neighbor_list) {
          double dist = CalculateDistance(temp, nb);
          if(dist <= fuel && !visited[nb]) s.push(nb);
        }
      }
    }
    return ans;
}

/**
 * CreateGraphFromCSVFile: Read the map data from the csv file
 * We have provided the code for you. Please do not need to change this function.
 */
void TrojanMap::CreateGraphFromCSVFile() {
  // Do not change this function
  std::fstream fin;
  fin.open("src/lib/data.csv", std::ios::in);
  std::string line, word;

  getline(fin, line);
  while (getline(fin, line)) {
    std::stringstream s(line);

    Node n;
    int count = 0;
    while (getline(s, word, ',')) {
      word.erase(std::remove(word.begin(), word.end(), '\''), word.end());
      word.erase(std::remove(word.begin(), word.end(), '"'), word.end());
      word.erase(std::remove(word.begin(), word.end(), '{'), word.end());
      word.erase(std::remove(word.begin(), word.end(), '}'), word.end());
      if (count == 0)
        n.id = word;
      else if (count == 1)
        n.lat = stod(word);
      else if (count == 2)
        n.lon = stod(word);
      else if (count == 3)
        n.name = word;
      else {
        word.erase(std::remove(word.begin(), word.end(), ' '), word.end());
        if (isalpha(word[0])) n.attributes.insert(word);
        if (isdigit(word[0])) n.neighbors.push_back(word);
      }
      count++;
    }
    data[n.id] = n;
  }
  fin.close();
}
