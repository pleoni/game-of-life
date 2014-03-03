# ~Game of Life~

Purpose of this project is to make a useful performance measurement tool based on the "Game of Life" (http://en.wikipedia.org/wiki/Conway%27s_Game_of_Life), suitable on HPC systems.

The code is currently targeted at CINECA's Eurora cluster, which is equipped with Nvidia Tesla K20 and Intel Xeon Phi accelerators, and can be compiled to exploit the presence of any of the two.
K20 GPGPUs are programmed via OpenACC directives, which are currently best supported by the Portland C compiler; on the other hand, the Intel C compiler is required to target MIC accelerators. OpenMP is used to provide multi-threading in a host-only or MIC environment, while MPI is used for inter-process communication on all platforms. Because of the different possible combinations, the Makefile has a few different targets which allow to select the correct compiler and to enable OpenMP or OpenACC support as needed.

On the Eurora cluster, the `module` command is used to load the appropriate environment and libraries for the chosen compiler, and this is what the Makefile script does prior to compilation and execution. If your system doesn't support such command, or if you want to use other compilers or libraries than the default ones, you can adapt the Makefile as needed. If you have manually setup the correct environment, you can ignore any "module: command not found" error. Please note that an MPI implementation is required.

## Compilation

To compile, type `make <target>`, where &lt;target&gt; is one of the following:

* `accpgi`: compile with the Portland C compiler, enabling OpenACC directives
* `omppgi`: compile with the Portland C compiler, enabling OpenMP directives
* `serpgi`: compile with the Portland C compiler, in "serial" mode (neither OpenACC nor OpenMP are enabled)
* `ompgnu`: compile with the GNU C compiler, enabling OpenMP directives
* `sergnu`: compile with the GNU C compiler, in "serial" mode
* `ompmic`: compile with the Intel C compiler, enabling OpenMP directives

Other compilation targets are available which generate multiple binaries in one run:

* `pgi`: compile with the Portland C compiler, in all supported acceleration modes
* `gnu`: compile with the GNU C compiler, in all supported acceleration modes
* `mic`: compile with the Intel C compiler, in all supported acceleration modes
* `all`: compile with all available compilers, in all supported acceleration modes

## Running

If you followed the above step, you should have (at least) an executable named `life_hpc2_<target>`. Simply run it in one of these two ways:

    ./life_hpc2_<target> <PARAMETERS>  # or
    mpirun -np N life_hpc2_<target> <PARAMETERS>

where N is the number of MPI processes you want to launch, and the `<PARAMETERS>` are:

> `-r#` - number of rows
> `-c#` - number of columns
> `-s#` - number of steps of the simulation
> `-n#` - number of artificial extra computations per cell, per step
> `-t#` - number of OpenMP threads for each MPI process (when compiled with OpenMP support)
> `-G#` - ID number of the GPU to which to offload the calculations (when compiled with openACC support)
> `-d#` - debug mode (0 = concise output, 1 = verbose output (default), 2 = display interactively)
> `-f<NAME>`: output file (mainly for debugging)

## Testing

You can test for code correctness by typing:

    make tests

which will perform a short, standard test run on all versions of the code, and compare their outputs with a reference output file named 'test_reference.dat'. (You will get an error for each version of the code which you have not compiled, but you can ignore them). You can also run single-target tests by typing, for example:

    make test_accpgi

Please refer to the Makefile for all other available targets.
