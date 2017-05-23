# Cholesky Solver for a Tree or General Graph

- Tree Plus Edges Problem (not dense enough, will stay sparse most of the time, won’t have dense sub matrices)
- Existing solvers are going to fail because they are dealing with denser problems
- Parallelizing Cholesky

use ```pip install numpy``` and ```pip install scipy```

### libraries:
http://people.sc.fsu.edu/~jburkardt/c_src/csparse/csparse.html
http://math.nist.gov/MatrixMarket/mmio-c.html

### I need to learn how to write README
* Solver.c
	* reads given matrix, performs csparse's cs_cholsol()
	* to run: make; ./solve {file.mtx};
	
### TODO:
* more performance analysis?
* improve performance?


#### Todo 5/23/17:
- Write scripts to run Kevin's .mtx files
- Spit out results into a raw data file
- From those raw data files, parse and use relevant data for output or potential data plots
- Presentation next Tuesday on a projector to run our code and show results
