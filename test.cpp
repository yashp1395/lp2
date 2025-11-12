//o/1 knapsack
#include <iostream>
#include <vector>
#include <algorithm>
using namespace std;

// Function to solve 0/1 Knapsack using DP
void Knapsack(vector<int> Profit, vector<int> Weight, int n, int Capacity) {
    // DP table initialization: (n+1) x (Capacity+1), all values set to 0
    vector<vector<int>> dp(n + 1, vector<int>(Capacity + 1, 0));
    vector<int>itemsInKnapsack; // Stores indices of selected items

    // Fill DP table using bottom-up approach
    for (int i = 1; i <= n; i++) {
        for (int j = 0; j <= Capacity; j++) {
            if (Weight[i - 1] <= j) {
                // Option to include or exclude the item
                dp[i][j] = max(dp[i - 1][j], dp[i - 1][j - Weight[i - 1]] + Profit[i - 1]);
            } else {
                // Can't include the item
                dp[i][j] = dp[i - 1][j];
            }
        }
    }

    // Print the DP table
    cout << "\nDP Table:\n\n";
    cout << "\t"; // Top-left empty cell for column header
    for (int j = 0; j <= Capacity; j++) {
        cout << j << "\t"; // Column headers = capacity values
    }
    cout << endl;

    for (int i = 0; i <= n; i++) {
        cout << "i=" << i << "\t"; // Row header = item index
        for (int j = 0; j <= Capacity; j++) {
            cout << dp[i][j] << "\t";
        }
        cout << endl;
    }

    // Maximum profit from DP table
    cout << "\nFor given capacity, maximum profit of items stored in knapsack is = " << dp[n][Capacity] << endl;

    // Backtrack to find items included in knapsack
    int remainingCapacity = Capacity;
    for (int i = n; i >= 1; i--) {
        if (dp[i][remainingCapacity] != dp[i - 1][remainingCapacity]) {
            itemsInKnapsack.push_back(i - 1); // Store 0-based index
            remainingCapacity -= Weight[i - 1]; // Decrease remaining capacity
        }
    }

    // Reverse to maintain original order
    reverse(itemsInKnapsack.begin(), itemsInKnapsack.end());

    // Print selected items (converted to 1-based index)
    cout << "Items included in the knapsack are (1-based index): ";
    for (auto item : itemsInKnapsack) {
        cout << item + 1 << " ";
    }
    cout << endl;
}

int main() {
    int n, Capacity;
    
    // Input: number of items
    cout << "Enter number of items: ";
    cin >> n;

    vector<int> Profit(n), Weight(n);

    // Input: profit and weight of each item
    cout << "Enter Profit and Weight of each item:\n";
    for (int i = 0; i < n; i++) {
        cout << "Enter Profit of item " << i + 1 << ": ";
        cin >> Profit[i];
        cout << "Enter Weight of item " << i + 1 << ": ";
        cin >> Weight[i];
    }

    // Input: capacity of knapsack
    cout << "Enter Capacity of Knapsack: ";
    cin >> Capacity;

    // Print input summary
    cout << "\nItems: ";
    for (int i = 1; i <= n; i++) {
        cout << i << " ";
    }
    cout << "\nProfits: ";
    for (auto p : Profit) {
        cout << p << " ";
    }
    cout << "\nWeights: ";
    for (auto w : Weight) {
        cout << w << " ";
    }
    cout << endl;

    // Call Knapsack function
    Knapsack(Profit, Weight, n, Capacity);
    return 0;
}


// -------------------------------------------------------------------------------------------------------------------------------------
fractional knap
#include <bits/stdc++.h>
using namespace std;

// Structure to store item details
struct Item {
    int id;         // Item ID
    int profit;     // Profit of item
    int weight;     // Weight of item
    float ratio;    // Profit/Weight ratio

    // Constructor to initialize item and calculate ratio
    Item(int id, int profit, int weight) {
        this->id = id;
        this->profit = profit;
        this->weight = weight;
        this->ratio = float(profit) / weight; // Calculate profit/weight ratio
    }
};

// Partition function for quicksort based on ratio (descending)
int partition(vector<Item>& items, int low, int high) {
    float pivot = items[high].ratio; // Choose pivot as last element's ratio
    int i = low - 1;

    for (int j = low; j < high; j++) {
        if (items[j].ratio > pivot) { // If current ratio > pivot, swap
            i++;
            swap(items[i], items[j]);
        }
    }

    swap(items[i + 1], items[high]); // Place pivot in correct position
    return i + 1;
}

// Recursive quicksort to sort items by ratio in descending order
void quickSort(vector<Item>& items, int low, int high) {
    if (low < high) {
        int pi = partition(items, low, high); // Partition index
        quickSort(items, low, pi - 1);        // Sort left subarray
        quickSort(items, pi + 1, high);       // Sort right subarray
    }
}

// Function to perform fractional knapsack using P/W ratio
float fractionalKnapsack(vector<Item> items, int W, float &totalWeightUsed, vector<int> &itemsTaken) {
    quickSort(items, 0, items.size() - 1); // Sort items by P/W ratio

    float profit = 0;
    totalWeightUsed = 0;

    for (auto& item : items) {
        if (W == 0) break; // Stop if knapsack is full

        if (item.weight <= W){
            profit += item.profit;
            W -= item.weight;
            totalWeightUsed += item.weight;
            itemsTaken.push_back(item.id); // Record taken item
        } else {
            // Take fractional part of item
            profit += item.ratio * W;
            totalWeightUsed += W;
            itemsTaken.push_back(item.id); // Record partial item taken
            W = 0;
        }
    }

    return profit; // Return total profit
}

int main() {
    // Sample input data: IDs, Profits, and Weights
    vector<int> obj = {1, 2, 3, 4, 5, 6, 7};
    vector<int> profit = {5, 10, 15, 7, 8, 9, 4};
    vector<int> wt = {1, 3, 5, 4, 1, 3, 2};

    int W;
    cout << "Enter capacity of knapsack: ";
    cin >> W; // User input for knapsack capacity
    cout << endl;

    // Build vector of Item objects
    vector<Item> items;
    for (int i = 0; i < obj.size(); i++) {
        items.push_back(Item(obj[i], profit[i], wt[i]));
    }

    float total = 0;
    vector<int> itemsTaken;

    // Call fractional knapsack algorithm
    float maxProfit = fractionalKnapsack(items, W, total, itemsTaken);

    // Display results
    cout << "Maximum profit: " << maxProfit << endl;
    cout << "Total weight used: " << total << endl;
    cout << "Items taken (object IDs): ";
    for (int it : itemsTaken) {
        cout << it << " ";
    }
    cout << endl;

    return 0;
}




// -------------------------------------------------------------------------------------------------------------------------------------

bellmanford
#include <bits/stdc++.h>
using namespace std;

struct Edge {
    int u, v, wt; // from, to, weight
};

vector<int> bellman_ford(int V, vector<Edge> &edges, int src, vector<int> &parent) {  // 🔹 added parent vector
    vector<int> dist(V, 1e8);  // initially all distances = infinity
    dist[src] = 0;
    parent.assign(V, -1); // 🔹 initialize all parents as -1

    // Relax all edges (V-1) times
    for (int i = 0; i < V - 1; i++) {
        for (auto e : edges) {
            if (dist[e.u] != 1e8 && dist[e.u] + e.wt < dist[e.v]) {
                dist[e.v] = dist[e.u] + e.wt;
                parent[e.v] = e.u; // 🔹 store parent
            }
        }
    }

    // Check for negative weight cycle
    for (auto e : edges) {
        if (dist[e.u] != 1e8 && dist[e.u] + e.wt < dist[e.v]) {
            cout << "Graph contains a negative weight cycle!\n";
            return {-1};
        }
    }

    return dist;
}

int main() {
    int V, E;
    cout << "Enter number of vertices and edges: ";
    cin >> V >> E;

    vector<Edge> edges(E);
    cout << "Enter edges (u v wt):\n";
    for (int i = 0; i < E; i++) {
        cin >> edges[i].u >> edges[i].v >> edges[i].wt;
    }

    int src, dest;
    cout << "Enter source vertex: ";
    cin >> src;
    cout << "Enter destination vertex: "; // 🔹 new input
    cin >> dest;

    vector<int> parent; // 🔹 for path tracking
    vector<int> result = bellman_ford(V, edges, src, parent);

    if (result[0] != -1) {
        cout << "Shortest distances from source " << src << ":\n";
        for (int i = 0; i < V; i++)
            cout << "Vertex " << i << " : " << result[i] << "\n";

        // 🔹 Print shortest path to destination
        cout << "\nShortest path from " << src << " to " << dest << ": ";
        vector<int> path;
        for (int v = dest; v != -1; v = parent[v]) path.push_back(v);
        reverse(path.begin(), path.end());
        for (int v : path) cout << v << " ";
        cout << "\nTotal cost = " << result[dest] << "\n";
    }

    return 0;
}

-------------------------------------------------------------------------------------------------------------------------------------
//nqueens
#include<bits/stdc++.h>  // Import all standard libraries
using namespace std;
int total=0;  // To count total number of solutions

// Hash maps to keep track of attacks by queens
unordered_map<int,bool>rowCheck;         // Check if a row already has a queen
unordered_map<int,bool>leftUpperCheck;   // Check diagonal (upper-left to lower-right)
unordered_map<int,bool>leftLowerCheck;   // Check diagonal (lower-left to upper-right)


// Function to print the current board configuration
void printBoard(vector<vector<char>>& board,int n){
    for(int i=0;i<n;i++){
        for(int j=0;j<n;j++){
            cout<<board[i][j]<<" ";  // Print each cell
        }
        cout<<endl;
    }
    cout<<endl<<endl;
}

// Function to check if it's safe to place a queen at (row, col)
bool isSafe(int row,int col,vector<vector<char>>& board,int n){
    if(rowCheck[row]) return false;                   // Row already occupied
    if(leftUpperCheck[n-1+row-col]) return false;     // Left-upper diagonal occupied
    if(leftLowerCheck[row+col]) return false;         // Left-lower diagonal occupied

    return true;  // Safe to place queen
}

// Recursive function to solve N-Queens problem
void solve(vector<vector<char>>& board,int col,int n){
    // Special case for n=1
    if(n==1){
        cout<<'Q'<<endl;
        return;
    }
    
    // If all queens are placed, print solution
    if(col>=n){
        printBoard(board,n);
        total++;   // Count this arrangement
        return ;
    }

    // Try placing queen in each row of the current column
    for(int row=0;row<n;row++){
        if(isSafe(row,col,board,n)){   // Check safety before placing
            // Place the queen
            board[row][col]='Q';
            rowCheck[row]=true;
            leftUpperCheck[n-1+row-col]=true;
            leftLowerCheck[row+col]=true;

            // Recurse for next column
            solve(board,col+1,n);

            // Backtrack: remove queen and reset checks
            board[row][col]='-';
            rowCheck[row]=false;
            leftUpperCheck[n-1+row-col]=false;
            leftLowerCheck[row+col]=false;
        }
    }
}


int main(){
    int n;
    cout<<"Enter the number of Queen : ";
    cin>>n;

    // For n=2 and n=3, no solution exists (except n=1 special case)
    if(n<=3 && n!=1) cout<<"Cant place Queen"<<endl;
    else{
        // Initialize the chessboard with '-'
        vector<vector<char>>board(n,vector<char>(n,'-'));

        // Start solving from column 0
        cout<<"The following are the N-Queen arrangements : "<<endl<<endl;
        solve(board,0,n);
        
        if(n!=1) cout<<"Total number of possibilities to place the queens is : "<< total<<endl;
        else cout<<"Total number of possibilities to place the queens is : 1" <<endl;
        // Print total solutions
    }
    return 0;
}

// -------------------------------------------------------------------------------------------------------------------------------------
// subset
#include <iostream>
using namespace std;

bool foundSolution = false;

void findSubsetSums(int index, int n, int set[], int targetSum, int subset[], int subsetSize) {
    // Base case: if targetSum is zero and subset is non-empty, print the subset in reverse order
    if (targetSum == 0 && subsetSize > 0) {
        foundSolution = true;
        cout << "[ ";
        for (int i = subsetSize - 1; i >= 0; i--) {  // Print in reverse order
            cout << subset[i] << " ";
        }
        cout << "] ";
        return;
    }

    // If all elements are considered and no solution found, return
    if (index == n) return;

    // Recur without including the current element
    findSubsetSums(index + 1, n, set, targetSum, subset, subsetSize);

    // Include the current element and recur with reduced target
    subset[subsetSize] = set[index];
    findSubsetSums(index + 1, n, set, targetSum - set[index], subset, subsetSize + 1);
}

int main() {
    int n, targetSum;

    // Input number of elements in the set
    cout << "Enter the number of elements in the set: ";
    cin >> n;

    int set[n];

    // Input the elements of the set
    cout << "Enter the elements of the set: ";
    for (int i = 0; i < n; i++) {
        cin >> set[i];
    }

    // Input the target sum
    cout << "Enter the target sum: ";
    cin >> targetSum;

    int subset[n];  // Array to store current subset

    cout << "Subsets that sum to " << targetSum << ": ";
    findSubsetSums(0, n, set, targetSum, subset, 0);

    // If no solution is found
    if (!foundSolution) {
        cout << "No subset found that sums to " << targetSum << endl;
    }

    return 0;
}
// --------------------------------------------------------------------------------------------------------------------------------------
