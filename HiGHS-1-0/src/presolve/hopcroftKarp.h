#include <vector>
#include <queue>
#include <climits>
#define hZERO 0
#define hINF INT_MAX
struct matching{
    int numMatched;
    std::vector<int> pairL;
    std::vector<int> pairR;
};
// Hopcroft Karp Algorithm implemented as a class
class hopcroftKarp{
public:
    // Let side vertices
    int m;
    // Right side vertices
    int n;
    struct matching swap;
    // Adjacency list for right side vertices
    std::vector<int> adj;
    std::vector<int> adjStart;
    // Distance from u to u'
    std::vector<int> dist;
    // Set up
    void setUp(int mNodes, int nNodes, std::vector<int> conn,
                 std::vector<int> connStart);
    // BFS for augmenting paths
    bool bfs();
    // DFS for augmenting path beginning with vertex u
    bool dfs(int u);
    // hopcroftKarp algo
    matching hopKarp();
};