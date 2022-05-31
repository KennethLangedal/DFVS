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

    // ofstream fs("scripts/plot_data");
    // for (auto u : g.active_vertices()) {
    //     for (auto v : g.out(u)) {
    //         fs << u << " " << v << std::endl;
    //     }
    // }

    local_search ls(g, 0.2, 3);
    ls.greedy_one_zero_swaps();
    ls.greedy_one_zero_swaps();

    double avg_move_cost = ls.get_average_move_cost(g);
    double T_offset = (std::max(std::min(avg_move_cost, 7.0), 0.0) / 7.0) * 0.1;
    double Ts_offset = (1.0 - ((double)ls.get_best_cost() / (double)g.size())) * 0.3;

    double Ts = 0.15 + Ts_offset + T_offset, Tf = 0.05 + T_offset, T = Ts, alfa = 0.999;

#ifdef VERBOSE
    std::cout << "Ts " << Ts << ", Tf " << Tf << std::endl;
#endif

    uint32_t imp = g.size();

#ifdef VERBOSE
    std::cout << "Average move cost: " << ls.get_average_move_cost(g) << std::endl;
#endif

    uint64_t it = 0;
    int32_t two_one_imps_left = 0;

    while (!tle) {
        ls.set_temperature(T);
        ls.search(g, g.size() * (1.0 - ((T - Tf) / (Ts - Tf))));

        if (imp > ls.get_best_cost()) {
            imp = ls.get_best_cost();
        } else {
            T *= alfa;
            if (T < Tf) {
                ls.greedy_one_zero_swaps();
                ls.greedy_one_zero_swaps_dfs(g);
                if (g.size() < 10000) {
                    while (ls.check_every_two_one_swap(g)) {
                        two_one_imps_left++;
                        ls.greedy_one_zero_swaps();
                        ls.greedy_one_zero_swaps_dfs(g);
                    }
                }
                // if (ls.get_current_cost() > ls.get_best_cost())
                //     ls.return_to_best(g);

                // ls.random_walk(10);
                T = Ts;
            }
        }

#ifdef VERBOSE
        if ((it++ & ((1 << 4) - 1)) == 0)
            std::cout << "\x1b[2K" << ls.get_best_cost() + cost << " " << ls.get_current_cost() + cost << " " << T << " " << two_one_imps_left << '\r' << std::flush;
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

    reduce_graph(g, res, gs, re, false);

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

    cout << "\x1b[2K";

#ifdef _VERBOSE
    cout << "\x1b[2K" << res.popcount() << endl;
#else
    for (auto u : res) {
        cout << u + 1 << "\n";
    }
    cout << flush;
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