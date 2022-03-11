#include "solver.hpp"
#include <optional>
#include <fstream>

size_t find_dense_branch_node(const graph &g, const bitvector &nodes) {
    size_t best = g.size(), best_cnt;
    for (size_t u : nodes) {
        size_t u_cnt = g.out_degree(u) + g.in_degree(u);
        if (u_cnt > 0 && (best == g.size() || u_cnt > best_cnt)) {
            best = u;
            best_cnt = u_cnt;
        }
    }

    return best;
}

bitvector solve(graph &g, const bitvector &nodes, graph_search &gs, size_t d) {
    size_t t0 = g.get_timestamp();
    bitvector fvs(g.size());
    std::vector<bitvector> SCC;
    reduce_graph(g, fvs, nodes, gs, SCC);
    if (nodes.intersection_size(g.active_vertices()) == 0) {
        g.unfold_graph(t0, fvs);
        return fvs;
    }
    if (d > 10) {
        std::cout << g.active_vertices().popcount() << std::endl;
        std::ofstream fs("scripts/plot_data");
        g.print_edgelist(fs, g.active_vertices());
        exit(0);
    }
    // if (popcount(remaining) < 30) {
    //     std::cout << popcount(remaining) << std::endl;
    //     std::ofstream fs("scripts/plot_data");
    //     g.print_edgelist(fs, remaining);
    //     exit(0);
    // }
    // for each CC
    for (auto &&CC : SCC) {
        if (CC.popcount() == 0) {
            continue;
        }
        size_t u = find_dense_branch_node(g, CC);
        size_t t1 = g.get_timestamp();

        // exclude u
        bitvector fvs_ex_u(g.size());
        exclude_from_fvs(g, fvs_ex_u, gs, u);
        fvs_ex_u |= solve(g, CC, gs, d + 1);
        g.unfold_graph(t1, fvs_ex_u);

        // include u
        bitvector fvs_inc_u(g.size());
        add_to_fvs(g, fvs_inc_u, gs, u);
        fvs_inc_u |= solve(g, CC, gs, d + 1);
        g.unfold_graph(t1, fvs_inc_u);

        if (fvs_inc_u.popcount() < fvs_ex_u.popcount()) {
            fvs |= fvs_inc_u;
        } else {
            fvs |= fvs_ex_u;
        }
    }

    g.unfold_graph(t0, fvs);

    return fvs;
}

branch_info analysis(graph &g, graph_search &gs) {
    bitvector fvs(g.size());
    std::vector<bitvector> SCC;
    reduce_graph(g, fvs, g.active_vertices(), gs, SCC);
    branch_info bi;
    if (g.active_vertices().popcount() == 0) {
        return bi;
    }
    // size_t t0 = g.get_timestamp();
    // visit(g.active_vertices(), [&](size_t u) {
    //     exclude_from_fvs(g, fvs, gs, u);
    //     reduce_graph(g, fvs, g.active_vertices(), gs, SCC);
    //     size_t exc_num = popcount(g.active_vertices());
    //     bi.avg_exc += exc_num;
    //     bi.best_exc = std::min(bi.best_exc, exc_num);
    //     bi.worst_exc = std::max(bi.worst_exc, exc_num);
    //     g.unfold_graph(t0, fvs);

    //     add_to_fvs(g, fvs, gs, u);
    //     reduce_graph(g, fvs, g.active_vertices(), gs, SCC);
    //     size_t inc_num = popcount(g.active_vertices());
    //     bi.avg_inc += inc_num;
    //     bi.best_inc = std::min(bi.best_inc, inc_num);
    //     bi.worst_inc = std::max(bi.worst_inc, inc_num);
    //     g.unfold_graph(t0, fvs);
    // });
    // bi.avg_exc /= popcount(g.active_vertices());
    // bi.avg_inc /= popcount(g.active_vertices());
    return bi;
}