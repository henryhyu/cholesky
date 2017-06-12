# Cholesky Solver for a Tree or General Graph

### libraries:
- http://people.sc.fsu.edu/~jburkardt/c_src/csparse/csparse.html
- http://math.nist.gov/MatrixMarket/mmio-c.html
- matplotlib, numpy
- https://github.com/json-c/json-c

### File description
* Solver.c
  * reads given matrix, performs csparse's cs_cholsol()
  * to run: make; ./solve {file.mtx};
* grapher/bargraph.py
  * creates various bargraphs
* Solver.c
  * outputs json
  
### TODO:
* more performance analysis?
* improve performance?

### Notes:
- perf stat -x, -e cache-misses,cache-references ./exec
- Tree Plus Edges Problem (not dense enough, will stay sparse most of the time, won’t have dense sub matrices)
- Existing solvers are going to fail because they are dealing with denser problems
- Parallelizing Cholesky
-

- LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/local/lib/
- export LD_LIBRARY_PATH