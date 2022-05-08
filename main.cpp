#include <fstream>
#include <iostream>

// #define HEURISTIC
#define EXACT

#include "solver.hpp"

using namespace std;

#ifdef EXACT

int main(int, char **)
{
    graph g;
    cin >> g;

    // cout << g.size() << endl;

    auto fvs = solve(g);

    cout << endl
         << "Solution size " << fvs.popcount() << endl;

    // for (size_t v : fvs) {
    //     cout << v + 1 << endl;
    // }

    return 0;
}

#endif

#ifdef HEURISTIC

#include "local_search.hpp"

int main(int, char **)
{
    graph g;
    cin >> g;

    local_search ls(g.size(), 0);

    for (size_t i = 0; i < 10; ++i)
        ls.search(g);

    cout << endl
         << ls.get_best().popcount() << endl;

    // for (size_t u : ls.get_best()) {
    //     cout << u + 1 << endl;
    // }

    return 0;
}

#endif