#include <hopcroftKarp.h>

struct matching hopcroftKarp::hopKarp(){
    swap.pairL.assign(m + 1, hZERO);
    swap.pairR.assign(n + 1, hZERO);
    dist.assign(m + 1, hINF);
    swap.numMatched = 0;
    while (bfs())
        for (int u = 1; u <= m; ++u)
            if (!swap.pairL[u] && dfs(u)) ++swap.numMatched;
    return swap;
}

bool hopcroftKarp::bfs(){
    std::queue<int> nodeQ;
    for (int u = 1; u <= m; ++u){
        if (!swap.pairL[u]){
            dist[u] = 0;
            nodeQ.push(u);
        }
        else 
            dist[u] = hINF;
    }
    dist[hZERO] = hINF;
    while (!nodeQ.empty()){
        int u = nodeQ.front();
        nodeQ.pop();
        if (dist[u] < dist[hZERO]){
            for (int i = adjStart[u - 1]; i < adjStart[u]; ++i){
                if (dist[swap.pairR[adj[i]]] == hINF){
                    dist[swap.pairR[adj[i]]] = dist[u] + 1;
                    nodeQ.push(swap.pairR[adj[i]]);
                }
            }
        }
    }
    return (dist[hZERO] != hINF);
}

bool hopcroftKarp::dfs(int u){
    if (u){
        for (int i = adjStart[u - 1]; i < adjStart[u]; ++i){
            if (dist[swap.pairR[adj[i]]] == dist[u] + 1){
                if (dfs(swap.pairR[adj[i]])){
                    swap.pairR[adj[i]] = u;
                    swap.pairL[u] = adj[i];
                    return true; 
                }
            }
        }
        dist[u] = hINF;
        return false;
    }
    return true;
}

void hopcroftKarp::setUp(int mNodes, int nNodes, std::vector<int> conn,
                            std::vector<int> connStart){
    m = mNodes;
    n = nNodes;
    adj = conn;
    adjStart = connStart;
}