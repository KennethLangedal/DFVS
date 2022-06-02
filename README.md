# Directed Feedback Vertex Set: PACE 2022

This repository contains both an exact and heuristic solver. However, only the heuristic solver is relevant for the PACE challenge 2022.

The solver uses reduction rules to preprocess the graph, followed by a simulated annealing procedure that repeats until the SIGTERM is received. At the end of each cooling schedule, the heuristic also greedily looks for one-one and two-one swaps before starting a new cooling cycle.

There are no external dependencies. Compile using the included CMake file or manually using the corresponding main file, **main_exact.cpp** or **main_heuristic.cpp**. The only requirement is **c++ 17**.

For information on the input-output format, see [PACE 2022 - Track Details](https://pacechallenge.org/2022/tracks/).