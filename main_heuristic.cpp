#include "bounds.hpp"
#include "local_search.hpp"
#include "local_search_edgew.hpp"
#include "reductions.hpp"
#include <algorithm>
#include <cmath>
#include <cstring>
#include <fstream>
#include <signal.h>
#include <sstream>
#include <unistd.h>

using namespace std;

volatile sig_atomic_t tle = 0;

void term(int signum) {
    tle = 1;
}

/*
double Ts = 0.15 + ((1.0 - (std::min((double)g.size(), 100000.0) / 100000.0)) * 0.2)
Ts = 0.15 + ((1.0 - (std::min((double)ls.get_best().popcount(), 100000.0) / 100000.0)) * 0.15);
*/

bitvector search_until_signal(sparse_graph g, size_t cost = 0) {
    double Ts = 0.35, Tf = 0.08, T = Ts, alfa = 0.999;

#ifdef VERBOSE
    std::cout << "Ts " << Ts << ", Tf " << Tf << std::endl;
#endif

    local_search ls(g, T, 0);
    uint32_t imp = g.size(), no_imp = 0;
    ls.greedy_one_zero_swaps(g);

    while (!tle) {
        ls.set_temperature(T);
        ls.search(g, 10000);

        if (imp > ls.get_best().popcount()) {
            imp = ls.get_best().popcount();
            no_imp = 0;
        } else if (no_imp++ > 3) {
            T *= alfa;
            if (T < Tf) {
                ls.shuffle_solution(g);
                ls.greedy_one_zero_swaps(g);
                T = Ts;
            }
            no_imp = 0;
        }
#ifdef VERBOSE
        std::cout << "\x1b[2K" << ls.get_best().popcount() + cost << " " << ls.get_current().popcount() + cost << " " << T << '\r' << std::flush;
#endif
    }
    return ls.get_best();
}

void solve_heuristic(sparse_graph &g) {
    graph_search gs(g.size());
    bitvector res(g.size());
    reduction_engine re;

#ifdef VERBOSE
    cout << "Starting reductions" << endl;
#endif

    reduce_graph(g, res, gs, re, true);

    if (g.active_vertices().popcount() == 0) {
        re.unfold_graph(g, 0, res);
#ifdef VERBOSE
        cout << res.popcount() << endl;
#else
        for (auto u : res) {
            cout << u + 1 << endl;
        }
#endif
        return;
    }

    size_t total_non_pi_degree = 0, total_pi_degree = 0;
    for (auto u : g.active_vertices()) {
        total_non_pi_degree += g.out_degree_non_pi(u);
        total_pi_degree += g.pi_degree(u);
    }

#ifdef VERBOSE
    cout << "Start size: " << g.size() << ", Reduced size: " << g.active_vertices().popcount() << ", Components: " << gs.SCC.size() << ", Solution cost offset: " << res.popcount() << endl;

    cout << "Average non-pi degree: " << (double)total_non_pi_degree / (double)g.active_vertices().popcount() << ", Average pi degree: " << ((double)total_pi_degree / 2.0) / (double)g.active_vertices().popcount() << endl;
#endif

    sparse_graph _g(g, g.active_vertices());

    auto _res = search_until_signal(_g, res.popcount());
    for (auto u : _res) {
        res.set(_g.original_label(u));
    }

    re.unfold_graph(g, 0, res);

#ifdef VERBOSE
    cout << "\x1b[2K" << res.popcount() << endl;
#else
    for (auto u : res) {
        cout << u + 1 << endl;
    }
#endif
}

int main(int, char **) {

    struct sigaction action;
    memset(&action, 0, sizeof(struct sigaction));
    action.sa_handler = term;
    sigaction(SIGTERM, &action, NULL);

    sparse_graph g;
    cin >> g;
    solve_heuristic(g);

    return 0;
}