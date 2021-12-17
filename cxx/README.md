## C++ version
Here we describe the procdure for compiling and running the C++ implementation which includes error and permutation analysis.
### Dependencies

In addition to Nupack [1,2], the C++ implementation has several dependencies [3]. Required components are bundled into the repository (at the expense of some bloat) in the `./vendor/` directory except for Python-dev which includes the [Python/C API](https://docs.python.org/3/c-api/index.html) used to call Nupack 4.0 functions, which needs to be installed on the user's machine.

#### Linux

On Linux Python can be installed with development packages with the following command
```
>$ sudo apt install python-dev
```
(replacing `apt` with your package manager) or through a dedicated distribution such as Anaconda with command
```
>$ conda install -c conda-forge python-devtools
```

#### MacOS

On MacOS the development packages can be installed with Anaconda in the above manner or (for example) through [homebrew](https://brew.sh/). Once brew is installed, python (and developer headers etc.) can be installed by entering
```
brew install python
```
into a terminal. If python is already installed, it may not have the development headers and so you can type 
```
brew reinstall python
```
to install them.

### Compiling

Once the dependencies are installed, a `CMakeLists.txt` script in the `./cxx/` directory allows for compilation using the standard [CMake](https://cmake.org/) procedure, namely
```
>$ mkdir build
>$ cd build
>$ cmake .. 
>$ make
```
from the `./cxx/` directory. The executable `DNA_bind` will appear in the created `./cxx/build/` directory.

The above procedures have been tested on Ubuntu 20.04 with gcc version 9.3.0, and on MacOS BigSur version 11.6 with clang version 12.0.5.

### Running

The executable runs assuming it is is being calling from the `cxx` directory (i.e. with command `>$ ./build/DNA_bind` from the `cxx` directory). If not you need to provide a command line argument with value equal to the path to the directory `<repo root>/cxx/nupack4/` (containing the required Python script) found in the repository root, relative to the current working directory.

[1] [www.nupack.org](www.nupack.org)

[2] M.E. Fornace, N.J. Porubsky, and N.A. Pierce (2020). A unified dynamic programming framework for the analysis of interacting nucleic acid strands: enhanced models, scalability, and speed. [ACS Synth Biol, 9:2665-2678](https://pubs.acs.org/doi/abs/10.1021/acssynbio.9b00523), 2020. 

[3] [Python/C API](https://docs.python.org/3/c-api/index.html). [OptimLib](https://www.kthohr.com/optimlib.html). [Armadillo](http://arma.sourceforge.net/). [Boost](https://www.boost.org/).