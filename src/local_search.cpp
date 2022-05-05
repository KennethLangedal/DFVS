#include "local_search.hpp"
#include <algorithm>
#include <iostream>
#include <math.h>

void local_search::_search_step(const graph &g) {
    size_t best_move = g.size(), best_i;
    double best_score;
    bool best_pos;

    bool found = false;
    size_t t = 0;
    while (!found && t++ < 50000) {
        tmp0.fill();
        tmp0.set_and_not(tmp0, _config_flag);
        size_t u = *tmp0.begin().advance(_dist_int(_reng) % tmp0.popcount());
        bool pos = _dist_real(_reng) > 0.5;
        auto [i, score] = _move_score(g, u, pos);
        if (score <= 0 || expl(-score / T) > _dist_real(_reng)) {
            found = true;
            _apply_move(g, u, i, pos);
            if (_N - _config_flag.popcount() < _best) {
                _best = _N - _config_flag.popcount();
                _best_dfvs.fill();
                _best_dfvs.set_and_not(_best_dfvs, _config_flag);
                imp = true;
            }
        } else if (best_move == g.size() || score < best_score) {
            best_move = u;
            best_i = i;
            best_score = score;
            best_pos = pos;
        }
    }
    if (!found) {
        _apply_move(g, best_move, best_i, best_pos);
        if (_N - _config_flag.popcount() < _best) {
            _best = _N - _config_flag.popcount();
            _best_dfvs.fill();
            _best_dfvs.set_and_not(_best_dfvs, _config_flag);
            imp = true;
        }
    }
}

// void local_search::_search_step(const graph &g) {
//     size_t best_move = g.size(), best_i;
//     double best_score;
//     bool best_pos;

//     double n = _N - _config_flag.popcount();
//     double n_draw = (uint64_t)((1.0 - T) * _N);
//     tmp0.fill();
//     tmp0.set_and_not(tmp0, _config_flag);
//     for (size_t u : tmp0) {
//         if (n_draw / n >= _dist_real(_reng)) {
//             bool pos = _dist_real(_reng) > 0.5;
//             auto [i, score] = _move_score(g, u, pos);
//             if (best_move == g.size() || best_score > score) {
//                 best_move = u;
//                 best_i = i;
//                 best_score = score;
//                 best_pos = pos;
//             }
//             n_draw -= 1;
//         }
//         n -= 1;
//     }
//     _apply_move(g, best_move, best_i, best_pos);
//     if (_N - _config_flag.popcount() < _best) {
//         _best = _N - _config_flag.popcount();
//         _best_dfvs.fill();
//         _best_dfvs.set_and_not(_best_dfvs, _config_flag);
//         imp = true;
//         // std::cout << _best << " " << _best << " " << T << " " << std::flush << '\r';
//     }
// }

std::pair<size_t, double> local_search::_move_score(const graph &g, size_t u, bool pos) const {
    double res = -1;
    size_t i;
    tmp0.set_and(g.out(u), _config_flag);
    tmp1.set_and(g.in(u), _config_flag);
    if (pos) {
        i = _config_c;
        for (size_t v : tmp0) {
            i = std::min(_pos_in_config[v], i);
        }
        for (size_t v : tmp1) {
            if (_pos_in_config[v] >= i) {
                ++res;
            }
        }
    } else {
        i = 0;
        for (size_t v : tmp1) {
            i = std::max(_pos_in_config[v] + 1, i);
        }
        for (size_t v : tmp0) {
            if (_pos_in_config[v] < i) {
                ++res;
            }
        }
    }
    return {i, res};
}

void local_search::_apply_move(const graph &g, size_t u, size_t i, bool pos) {
    for (size_t j = 0; j < _N; ++j) {
        if (_pos_in_config[j] >= i) {
            ++_pos_in_config[j];
        }
    }
    _config_flag.set(u);
    _pos_in_config[u] = i;
    ++_config_c;

    if (pos) {
        tmp0.set_and(g.in(u), _config_flag);
        for (size_t v : tmp0) {
            if (_pos_in_config[v] > i) {
                _config_flag.reset(v);
            }
        }
    } else {
        tmp0.set_and(g.out(u), _config_flag);
        for (size_t v : tmp0) {
            if (_pos_in_config[v] < i) {
                _config_flag.reset(v);
            }
        }
    }
}

local_search::local_search(size_t N, size_t seed)
    : _pos_in_config(N, N), _config_flag(N), _best_dfvs(N), tmp0(N), tmp1(N), _reng(seed), _dist_int(), _dist_real(0.0, 1.0), _best(N), _N(N), _config_c(0),
      l1(200), l2(N * 5), T(0.6), a(0.99) {
    _best_dfvs.fill();
}

void local_search::search(const graph &g) {
    size_t l1_it = 0;
    T = 0.4;
    while (l1_it < l1) {
        size_t l2_it = 0;
        imp = false;
        while (l2_it < l2) {
            _search_step(g);
            ++l2_it;
        }
        if (!imp) {
            T = T * a;
#ifdef VERBOSE
            std::cout << _best << " " << _N - _config_flag.popcount() << " " << T << " " << std::flush << '\r';
#endif
            ++l1_it;
        }
    }
}

const bitvector &local_search::get_best() const {
    return _best_dfvs;
}