#include "Solver.hpp"
#include <iomanip>
#include <algorithm>
#include <limits>
#include <cstdlib>
#include <ctime>
#include <cstdint>

using namespace std;

#define QTD_VIZINHANCAS 3
// #define VALIDATE_ROUTES

typedef enum
{
    SWAP = 1,
    TWO_OPT,
    REINSERTION
} Vizinhança;

struct Saving
{
    double value;
    int fromStation;
    int toStation;
};
struct SavingKey { int a, b; };

// PRNG rápido (xorshift32)
static uint32_t __rng_state = 0x9e3779b9u;
static inline uint32_t fastRand()
{
    uint32_t x = __rng_state;
    x ^= x << 13;
    x ^= x >> 17;
    x ^= x << 5;
    __rng_state = x;
    return x;
}
static inline int fastRandBound(int n) { return (int)(fastRand() % (uint32_t)n); }
static constexpr int RCL_MAX = 32; // limite top-k da RCL para acelerar a construção

vector<Route> VND(vector<Route> routes, const Instance *instance);
vector<Route> performSwap(vector<Route> routes, const Instance *instance);
vector<Route> performTwoOpt(vector<Route> routes, const Instance *instance);
vector<Route> performReinsertion(vector<Route> routes, const Instance *instance);
vector<Route> GRASPConstruction(const Instance *instance);

static inline int computeRouteCost(const vector<int> &st, const vector<vector<int>> &dist)
{
    int cost = 0;

    for (int i = 0; i < (int)st.size() - 1; ++i)
        cost += dist[st[i]][st[i + 1]];

    return cost;
}

// KNN: listas de k-vizinhos por estação para podar vizinhanças
static std::vector<std::vector<int>> gKNN;
static std::vector<std::vector<uint8_t>> gIsKNN; // matriz booleana i->j
static int gKNN_k = 0;
static int gKNN_n = 0;
static double gAvgEdge = 0.0; // média global das arestas

static double computeAvgEdge(const Instance *ins)
{
    const auto &d = ins->get_distancias();
    int n = ins->get_qtdEstacoes();
    long long cnt = 0, sum = 0;
    for (int i = 0; i <= n; ++i)
        for (int j = 0; j <= n; ++j)
            if (i != j) { sum += d[i][j]; ++cnt; }
    return cnt ? (double)sum / (double)cnt : 0.0;
}

static void buildKNN(const Instance *ins, int k)
{
    int n = ins->get_qtdEstacoes();
    if (gKNN_n == n && gKNN_k == k && (int)gKNN.size() == n + 1)
        return;
    gKNN_n = n; gKNN_k = k;
    gKNN.assign(n + 1, {});
    gIsKNN.assign(n + 1, std::vector<uint8_t>(n + 1, 0));
    const auto &d = ins->get_distancias();
    std::vector<int> cand; cand.reserve(n - 1);
    for (int i = 1; i <= n; ++i)
    {
        cand.clear();
        for (int j = 1; j <= n; ++j) if (j != i) cand.push_back(j);
        std::nth_element(cand.begin(), cand.begin() + std::min(k, (int)cand.size()) - 1, cand.end(), [&](int a, int b){ return d[i][a] < d[i][b]; });
        int take = std::min(k, (int)cand.size());
        gKNN[i].reserve(take);
        std::partial_sort(cand.begin(), cand.begin() + take, cand.end(), [&](int a, int b){ return d[i][a] < d[i][b]; });
        for (int t = 0; t < take; ++t) { int nb = cand[t]; gKNN[i].push_back(nb); gIsKNN[i][nb] = 1; }
    }
}
static inline bool isKnn(int i, int j)
{
    if (i <= 0 || j <= 0) return true; // ignore depósito nas podas
    if (i >= (int)gIsKNN.size() || j >= (int)gIsKNN[i].size()) return true; // segurança
    return gIsKNN[i][j] != 0;
}

static inline double totalCost(const vector<Route> &routes, const Instance *ins)
{
    double cost = 0;
    const auto &d = ins->get_distancias();

    for (auto &route : routes)
    {
        for (int i = 0; i + 1 < (int)route.stations.size(); ++i)
        {
            cost += d[route.stations[i]][route.stations[i + 1]];
        }
    }

    return cost;
}

static inline int routeCapacity(const vector<int> &stations, const Instance *ins)
{
    const auto &demandas = ins->get_demandas();
    int load = 0, maxLoad = 0;

    for (int i = 1; i < (int)stations.size() - 1; ++i)
    {
        if (demandas[stations[i] - 1] < 0)
            load += abs(demandas[stations[i] - 1]);
    }

    maxLoad = load;
    for (int i = 1; i < (int)stations.size() - 1; ++i)
    {
        load += demandas[stations[i] - 1];
        maxLoad = max(maxLoad, load);
    }
    return maxLoad;
}

static inline bool validateQuiet(const vector<Route> &routes, const Instance *ins)
{
    vector<bool> vis(ins->get_qtdEstacoes() + 1, false);

    vis[0] = true;
    const auto &demandas = ins->get_demandas();

    for (auto &route : routes)
    {
        if (route.stations.size() < 3 || route.stations.front() != 0 || route.stations.back() != 0)
            return false;

        int init = 0;
        for (int i = 1; i < (int)route.stations.size() - 1; ++i)
        {
            int station = route.stations[i];

            if (station == 0 || vis[station])
                return false;
            if (demandas[station - 1] < 0)
                init += abs(demandas[station - 1]);
        }

        int load = init;
        for (int i = 1; i < (int)route.stations.size() - 1; ++i)
        {
            int station = route.stations[i];
            vis[station] = true;
            load += demandas[station - 1];

            if (load < 0 || load > ins->get_capVeiculos())
                return false;
        }
    }

    for (int s = 1; s <= ins->get_qtdEstacoes(); ++s)
    {
        if (!vis[s])
            return false;
    }

    return true;
}

bool validateRoutes(const vector<Route> &routes, const Instance *ins)
{
    vector<bool> vis(ins->get_qtdEstacoes() + 1, false);

    vis[0] = true;
    const auto &dm = ins->get_demandas();

    cout << "\n=== VALIDAÇÃO DAS ROTAS ===\n";
    for (int route = 0; route < (int)routes.size(); ++route)
    {
        auto &r = routes[route];

        if (r.stations.size() < 3 || r.stations.front() != 0 || r.stations.back() != 0)
            return false;

        int init = 0;
        for (int i = 1; i < (int)r.stations.size() - 1; ++i)
        {
            if (dm[r.stations[i] - 1] < 0)
                init += abs(dm[r.stations[i] - 1]);
        }

        int load = init;
        for (int i = 1; i < (int)r.stations.size() - 1; ++i)
        {
            int s = r.stations[i];

            if (s == 0)
                return false;
            if (vis[s])
                return false;

            vis[s] = true;
            load += dm[s - 1];
            if (load < 0)
                return false;
            if (load > ins->get_capVeiculos())
                return false;
        }
    }

    for (int s = 1; s <= ins->get_qtdEstacoes(); ++s)
        if (!vis[s])
            return false;

    return true;
}

vector<Route> VND(vector<Route> routes, const Instance *ins)
{
    vector<Vizinhança> ord = {SWAP, TWO_OPT, REINSERTION};
    int k = 0;
    while (k < QTD_VIZINHANCAS)
    {
        vector<Route> cand;
        switch (ord[k])
        {
        case SWAP:        cand = performSwap(routes, ins); break;
        case TWO_OPT:     cand = performTwoOpt(routes, ins); break;
        case REINSERTION: cand = performReinsertion(routes, ins); break;
        }
        if (totalCost(cand, ins) < totalCost(routes, ins))
        {
            routes = std::move(cand);
            k = 0;
        }
        else
        {
            ++k;
        }
    }
    return routes;
}

vector<Route> performSwap(vector<Route> routes, const Instance *ins)
{
    vector<Route> bestRoutes = routes;
    const auto &d = ins->get_distancias();
    double base = totalCost(routes, ins), bestC = base;
    const double STRONG = std::max(1.0, 0.25 * gAvgEdge);

    // intra
    for (int route = 0; route < (int)routes.size(); ++route)
    {
        auto &stations = routes[route].stations;
        if (stations.size() <= 3) continue;
        for (int i = 1; i < (int)stations.size() - 1; ++i)
        {
            for (int j = i + 1; j < (int)stations.size() - 1; ++j)
            {
                int a = stations[i-1], b = stations[i];
                int c = stations[i+1];
                int d1 = stations[j-1], e = stations[j], f = stations[j+1];
                int delta;
                if (i + 1 < j)
                {
                    delta = -d[a][b] - d[b][c] - d[d1][e] - d[e][f]
                            + d[a][e] + d[e][c] + d[d1][b] + d[b][f];
                }
                else // i+1 == j (adjacentes)
                {
                    delta = -d[a][b] - d[b][e] - d[e][f]
                            + d[a][e] + d[e][b] + d[b][f];
                }
                if (delta >= 0) continue;
                if (!(isKnn(b, e) || isKnn(e, b)) && (double)(-delta) < STRONG) continue;
                auto tmp = stations; std::swap(tmp[i], tmp[j]);
                int cap = routeCapacity(tmp, ins); if (cap > ins->get_capVeiculos()) continue;
                double cand = base + delta;
                if (cand < bestC)
                {
                    bestC = cand;
                    auto nb = routes;
                    nb[route].stations = std::move(tmp);
                    nb[route].capacity = cap;
                    bestRoutes = std::move(nb);
                    return bestRoutes; // first-improvement
                }
            }
        }
    }

    // inter
    for (int a = 0; a < (int)routes.size(); ++a)
    {
        for (int b = a + 1; b < (int)routes.size(); ++b)
        {
            if (routes[a].stations.size() <= 2 || routes[b].stations.size() <= 2) continue;
            for (int ia = 1; ia < (int)routes[a].stations.size() - 1; ++ia)
            for (int ib = 1; ib < (int)routes[b].stations.size() - 1; ++ib)
            {
                int a0 = routes[a].stations[ia-1];
                int x  = routes[a].stations[ia];
                int a1 = routes[a].stations[ia+1];
                int b0 = routes[b].stations[ib-1];
                int y  = routes[b].stations[ib];
                int b1 = routes[b].stations[ib+1];
                int deltaA = -d[a0][x] - d[x][a1] + d[a0][y] + d[y][a1];
                int deltaB = -d[b0][y] - d[y][b1] + d[b0][x] + d[x][b1];
                int delta = deltaA + deltaB;
                if (delta >= 0) continue;
                if (!(isKnn(x,y) || isKnn(y,x)) && (double)(-delta) < STRONG) continue;
                auto ta = routes[a].stations, tb = routes[b].stations; std::swap(ta[ia], tb[ib]);
                int capA = routeCapacity(ta, ins); if (capA > ins->get_capVeiculos()) continue;
                int capB = routeCapacity(tb, ins); if (capB > ins->get_capVeiculos()) continue;
                double cand = base + delta;
                if (cand < bestC)
                {
                    bestC = cand;
                    auto nb = routes;
                    nb[a].stations = std::move(ta);
                    nb[b].stations = std::move(tb);
                    nb[a].capacity = capA;
                    nb[b].capacity = capB;
                    bestRoutes = std::move(nb);
                    return bestRoutes; // first-improvement
                }
            }
        }
    }
    return bestRoutes;
}

vector<Route> performTwoOpt(vector<Route> routes, const Instance *ins)
{
    vector<Route> best = routes;
    const auto &d = ins->get_distancias();
    double base = totalCost(routes, ins), bestC = base;
    const double STRONG = std::max(1.0, 0.25 * gAvgEdge);
    for (int r = 0; r < (int)routes.size(); ++r)
    {
        auto &st = routes[r].stations;
        if (st.size() <= 4)
            continue;
        // int old = computeRouteCost(st, d);
        for (int i = 1; i < (int)st.size() - 2; ++i)
            for (int j = i + 1; j < (int)st.size() - 1; ++j)
            {
                int a = st[i - 1], b = st[i], c = st[j], e = st[j + 1];
                int delta = -d[a][b] - d[c][e] + d[a][c] + d[b][e];
                if (delta >= 0)
                    continue;
                if (!(isKnn(a, c) || isKnn(c, a) || isKnn(b, e) || isKnn(e, b)) && (double)(-delta) < STRONG) continue;
                auto tmp = st;
                reverse(tmp.begin() + i, tmp.begin() + j + 1);
                int cap = routeCapacity(tmp, ins);
                if (cap > ins->get_capVeiculos())
                    continue;
                double cand = base + delta; // delta-only evaluation
                if (cand < bestC)
                {
                    bestC = cand;
                    auto nb = routes;
                    nb[r].stations = move(tmp);
                    nb[r].capacity = cap;
                    best = move(nb);
                    return best; // first-improvement
                }
            }
    }
    return best;
}

vector<Route> performReinsertion(vector<Route> routes, const Instance *ins)
{
    vector<Route> best = routes;
    const auto &d = ins->get_distancias();
    double base = totalCost(routes, ins), bestC = base; // intra
    const double STRONG = std::max(1.0, 0.25 * gAvgEdge);
    for (int r = 0; r < (int)routes.size(); ++r)
    {
        auto &st = routes[r].stations;
        if (st.size() <= 3)
            continue;
        for (int i = 1; i < (int)st.size() - 1; ++i)
        {
            int node = st[i];
            for (int j = 1; j < (int)st.size() - 1; ++j)
            {
                if (i == j)
                    continue;
                int a = st[i-1];
                int c = st[i+1];
                // destino
                int dest = (j < i ? j : j - 1);
                int p = st[dest - 1];
                int s = st[dest];
                int delta = -d[a][node] - d[node][c] - d[p][s] + d[a][c] + d[p][node] + d[node][s];
                if (delta >= 0) continue;
                if (!(isKnn(node, p) || isKnn(p, node) || isKnn(node, s) || isKnn(s, node)) && (double)(-delta) < STRONG) continue;
                auto tmp = st; tmp.erase(tmp.begin() + i); tmp.insert(tmp.begin() + dest, node);
                int cap = routeCapacity(tmp, ins); if (cap > ins->get_capVeiculos()) continue;
                double cand = base + delta;
                if (cand < bestC)
                {
                    bestC = cand;
                    auto nb = routes;
                    nb[r].stations = move(tmp);
                    nb[r].capacity = cap;
                    best = move(nb);
                    return best; // first-improvement
                }
            }
        }
    }
    // inter (move)
    for (int a = 0; a < (int)routes.size(); ++a)
        for (int b = 0; b < (int)routes.size(); ++b)
        {
            if (a == b)
                continue;
            for (int ia = 1; ia < (int)routes[a].stations.size() - 1; ++ia)
            {
                int node = routes[a].stations[ia];
                for (int jb = 1; jb < (int)routes[b].stations.size(); ++jb)
                {
                    int a0 = routes[a].stations[ia-1];
                    int x  = routes[a].stations[ia];
                    int a1 = routes[a].stations[ia+1];
                    int pb = routes[b].stations[jb-1];
                    int sb = routes[b].stations[jb];
                    int deltaA = -d[a0][x] - d[x][a1] + d[a0][a1];
                    int deltaB = -d[pb][sb] + d[pb][x] + d[x][sb];
                    int delta = deltaA + deltaB;
                    if (delta >= 0) continue;
                    if (!(isKnn(x, pb) || isKnn(pb, x) || isKnn(x, sb) || isKnn(sb, x)) && (double)(-delta) < STRONG) continue;
                    auto ta = routes[a].stations, tb = routes[b].stations; ta.erase(ta.begin() + ia); tb.insert(tb.begin() + jb, x);
                    int capB = routeCapacity(tb, ins); if (capB > ins->get_capVeiculos()) continue;
                    bool removeA = (int)ta.size() == 2;
                    int capA = 0; if (!removeA) { capA = routeCapacity(ta, ins); if (capA > ins->get_capVeiculos()) continue; }
                    double cand = base + delta;
                    if (cand < bestC)
                    {
                        bestC = cand;
                        if (!removeA)
                        {
                            auto nb = routes;
                            nb[a].stations = move(ta);
                            nb[a].capacity = capA;
                            nb[b].stations = move(tb);
                            nb[b].capacity = capB;
                            best = move(nb);
                            return best; // first-improvement
                        }
                        else
                        {
                            vector<Route> nb;
                            nb.reserve(routes.size() - 1);
                            for (int r = 0; r < (int)routes.size(); ++r)
                                if (r != a)
                                    nb.push_back(routes[r]);
                            int idxB = b - (a < b ? 1 : 0);
                            nb[idxB].stations = move(tb);
                            nb[idxB].capacity = capB;
                            best = move(nb);
                            return best; // first-improvement
                        }
                    }
                }
            }
        }
    return best;
}

vector<Route> GRASPConstruction(const Instance *ins)
{
    vector<Route> routes;
    const auto &dm = ins->get_demandas();
    const auto &d = ins->get_distancias();
    float alpha = ins->alpha;

    // 1) Inicializa rotas singletons com ordem embaralhada
    vector<int> order; order.reserve(ins->get_qtdEstacoes());
    for (int i = 1; i <= ins->get_qtdEstacoes(); ++i) order.push_back(i);
    for (int i = (int)order.size() - 1; i > 0; --i) { int j = fastRandBound(i + 1); std::swap(order[i], order[j]); }
    for (int s : order) routes.push_back({{0, s, 0}, std::abs(dm[s - 1])});

    // Mapeia endpoints atuais de cada rota (first/last) e índice da rota que contém um endpoint
    auto getEndpoints = [&](const Route &r){ vector<int> ep; if (r.stations.size()>2){ ep.push_back(r.stations[1]); ep.push_back(r.stations[(int)r.stations.size()-2]); } return ep; };
    auto routeIndexOfEndpoint = [&](int station){
        for (int k = 0; k < (int)routes.size(); ++k){ if (routes[k].stations.size()>2){ int f=routes[k].stations[1]; int l=routes[k].stations[(int)routes[k].stations.size()-2]; if (f==station || l==station) return k; } }
        return -1; };

    // 2) Heap global de savings entre endpoints (max-heap)
    struct HeapItem{ double val; int a; int b; };
    auto cmp = [](const HeapItem &x, const HeapItem &y){ return x.val < y.val; }; // max-heap
    std::vector<HeapItem> heap;
    heap.reserve(routes.size()*routes.size());

    auto pushSaving = [&](int a, int b){ if (a==b) return; double s = (double)d[0][a] + (double)d[0][b] - (double)d[a][b]; heap.push_back({s,a,b}); };
    // inicial
    for (int i = 0; i < (int)routes.size(); ++i)
        for (int j = i+1; j < (int)routes.size(); ++j)
        {
            if (routes[i].capacity + routes[j].capacity > ins->get_capVeiculos()) continue;
            auto Ci = getEndpoints(routes[i]);
            auto Cj = getEndpoints(routes[j]);
            for (int si : Ci) for (int sj : Cj) pushSaving(si, sj);
        }
    std::make_heap(heap.begin(), heap.end(), cmp);

    // 3) Itera mesclagens usando heap + invalidação preguiçosa com threshold RCL
    while (!heap.empty() && routes.size()>1)
    {
        // Reconstrói RCL por threshold a partir do topo atual
        int RCLlim = std::min(128, std::max(32, (int)routes.size()));
        std::vector<HeapItem> rcl; rcl.reserve(RCLlim);
        // peek top válido e obtém threshold dinâmico
        // Tenta algumas extrações até achar um topo que leve a rotas distintas e ainda válidas por capacidade
        HeapItem top{}; bool foundTop=false; std::vector<HeapItem> popped;
        while (!heap.empty()){
            std::pop_heap(heap.begin(), heap.end(), cmp); auto cur = heap.back(); heap.pop_back();
            popped.push_back(cur);
            int ri = routeIndexOfEndpoint(cur.a); int rj = routeIndexOfEndpoint(cur.b);
            if (ri==-1 || rj==-1 || ri==rj) { continue; }
            if (routes[ri].capacity + routes[rj].capacity > ins->get_capVeiculos()) { continue; }
            top = cur; foundTop = true; break;
        }
        if (!foundTop) break;
        double bestS = top.val;
        double worstS = bestS - 1e6; // aproximação: sem revarrer todo heap; suficiente para threshold
        double thr = bestS - alpha * (bestS - worstS);
        rcl.push_back(top);
        // Recolhe mais candidatos do heap dentro do threshold e válidos
        int take = RCLlim - 1;
        while (take>0 && !heap.empty()){
            std::pop_heap(heap.begin(), heap.end(), cmp); auto cur = heap.back(); heap.pop_back();
            popped.push_back(cur);
            if (cur.val < thr) { continue; }
            int ri = routeIndexOfEndpoint(cur.a); int rj = routeIndexOfEndpoint(cur.b);
            if (ri==-1 || rj==-1 || ri==rj) continue;
            if (routes[ri].capacity + routes[rj].capacity > ins->get_capVeiculos()) continue;
            rcl.push_back(cur); --take;
        }

        // escolhe um da RCL
        auto chosen = rcl[fastRandBound((int)rcl.size())];
        // Recoloca itens popados que não foram escolhidos e ainda são potencialmente válidos
        auto push_back_heap = [&](const HeapItem &h){ heap.push_back(h); std::push_heap(heap.begin(), heap.end(), cmp); };
        for (auto &h : popped)
        {
            if (h.a==chosen.a && h.b==chosen.b && h.val==chosen.val) continue;
            int ri_test = routeIndexOfEndpoint(h.a); int rj_test = routeIndexOfEndpoint(h.b);
            if (ri_test==-1 || rj_test==-1 || ri_test==rj_test) continue;
            if (routes[ri_test].capacity + routes[rj_test].capacity > ins->get_capVeiculos()) continue;
            push_back_heap(h);
        }
        int ri = routeIndexOfEndpoint(chosen.a);
        int rj = routeIndexOfEndpoint(chosen.b);
        if (ri==-1 || rj==-1 || ri==rj) continue; // segurança
        if (routes[ri].capacity + routes[rj].capacity > ins->get_capVeiculos()) continue;

        // realiza mescla ri <- rj, orientando pelos endpoints escolhidos (a deve ser o último interior de ri, b o primeiro de rj)
        auto ra = routes[ri].stations;
        auto rb = routes[rj].stations;
        if ((int)ra.size() > 2)
        {
            int fi = ra[1];
            int li = ra[(int)ra.size()-2];
            if (chosen.a == fi)
                std::reverse(ra.begin()+1, ra.end()-1);
        }
        if ((int)rb.size() > 2)
        {
            int fj = rb[1];
            int lj = rb[(int)rb.size()-2];
            if (chosen.b == lj)
                std::reverse(rb.begin()+1, rb.end()-1);
        }
        // agora ra termina em 'a' e rb começa em 'b'
        auto ns = ra; ns.pop_back();
        for (int t = 1; t < (int)rb.size(); ++t) ns.push_back(rb[t]);
        int newCap = routes[ri].capacity + routes[rj].capacity;

        // remove rj e ri mantendo índices
        int a = std::max(ri, rj), b = std::min(ri, rj);
        routes.erase(routes.begin() + a);
        routes.erase(routes.begin() + b);
        routes.push_back({ns, newCap});

        // Atualiza savings: gera savings do novo endpoint com todas as rotas atuais
        int newIdx = (int)routes.size()-1;
        auto Cn = getEndpoints(routes[newIdx]);
        for (int k = 0; k < (int)routes.size()-1; ++k)
        {
            if (routes[k].capacity + routes[newIdx].capacity > ins->get_capVeiculos()) continue;
            auto Ck = getEndpoints(routes[k]);
            for (int si : Cn) for (int sj : Ck)
            {
                double s = (double)d[0][si] + (double)d[0][sj] - (double)d[si][sj];
                heap.push_back({s, si, sj});
                std::push_heap(heap.begin(), heap.end(), cmp);
            }
        }
    }
    return routes;
}

void printRoute(const Route &route)
{
    cout << "[";
    for (int i = 0; i < (int)route.stations.size(); ++i)
    {
        cout << route.stations[i];
        if (i + 1 < (int)route.stations.size())
            cout << " -> ";
    }
    cout << "] capacidade: " << route.capacity << ")\n";
}
void printAllRoutes(const vector<Route> &routes, const Instance *ins)
{
    cout << "=== ROTAS ENCONTRADAS ===\n";
    cout << "número total de rotas: " << routes.size() << "\n";
    cout << "custo total: " << totalCost(routes, ins) << "\n\n";
    for (int i = 0; i < (int)routes.size(); ++i)
    {
        cout << "rota " << (i + 1) << ": ";
        printRoute(routes[i]);
    }
    cout << endl;
}

void Solver::Solve(Instance *ins)
{
    cout << "\n=== Iniciando GRASP ===\n";
    srand((unsigned)time(nullptr));
    gAvgEdge = computeAvgEdge(ins);
    // k dinâmico: limita por n e mantém mínimo razoável para qualidade
    int knn_k = std::min(ins->get_qtdEstacoes(), std::max(32, ins->get_qtdEstacoes()/3));
    buildKNN(ins, knn_k); // k-vizinhos para podas
    vector<Route> bestSol;
    double bestCost = numeric_limits<double>::max();
    double opt = ins->valorOtimo;
    for (int it = 0; it < ins->repetitions; ++it)
    {
        cout << "\n--- Iteração " << (it + 1) << " ---\n";
        auto routes = GRASPConstruction(ins);
        double c1 = totalCost(routes, ins);
        cout << "Custo após construção: " << c1 << "\n";
        auto improved = VND(routes, ins);
        double c2 = totalCost(improved, ins);
        cout << "Custo após VND: " << c2 << "\n";
        if (validateQuiet(improved, ins))
        {
            if (c2 < bestCost)
            {
                bestCost = c2;
                bestSol = improved;
                cout << "*** NOVA MELHOR SOLUÇÃO VÁLIDA! (Gap: " << fixed << setprecision(2) << (((bestCost - opt) / opt) * 100.0) << "%) ***\n";
            }
        }
        if (abs(bestCost - opt) < 0.01)
        {
            cout << "*** ÓTIMO ENCONTRADO! Parando execução. ***\n";
            break;
        }
    }
    cout << "\n=== RESULTADO FINAL GRASP ===\n"
         << fixed << setprecision(2);
    cout << "--------------------------------------------------------------------------\n";
    cout << "| Instancia | Custo GRASP + VND     | Valor Otimo |   Gap (%)   |\n";
    cout << "--------------------------------------------------------------------------\n";
    cout << "|" << setw(10) << ins->instanceName << " |" << setw(22) << bestCost << " |" << setw(12) << opt << " |" << setw(11) << (((bestCost - opt) / opt) * 100.0) << " |\n";
    cout << "--------------------------------------------------------------------------\n";
    cout << "\n=== MELHOR SOLUÇÃO ENCONTRADA ===\n";
    printAllRoutes(bestSol, ins);
}