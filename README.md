# DNA-bind

A simple model for capturing DNA hybridisation rates, fit to experimental data, as appearing in [1]. There are two implementations, one in C++ and one in Python, detailed below.

## Data

Program output containing all model data for the supplied sequences and experimental rates, using Nupack 4.0.0.21 [2,3], for bind lengths of 1-4 bases, on debug setting 0 can be found in the [data directory](./data).

## Nupack dependency

Both implementations rely on [Nupack version 4.0](http://www.nupack.org/) [2,3] which can be downloaded [here](http://www.nupack.org/downloads). Installation instructions for MacOS/Linux can be found [here](https://docs.nupack.org/start/#maclinux-installation).

## Python version

A simple implementation which performs only point-estimate optimisations is contained in the script `dna_bind.py` in the `./python/` directory of the repository root. It should be run using `python3`. I.e.

```
>$ python3 ./python/dna_bind.py
```
This file acts as the public interface where you can adjust parameters, e.g. `bindLength`. The file imports functions from `dna_bind_functions.py` which in contrast should generally be left alone.

## C++ version

The full model is provided as a C++ implementation, including error and permutation analysis, as well as the ability for users to design similar models by deriving sub-classes. Details are provided in its dedicated [README](./cxx/README.md) file in its self contained `./cxx/` directory.

## References

[1] S. Hertel and R.E. Spinney et al. (2021). Mechanisms underlying sequence-dependent DNA hybridisation rates in the absence of secondary structure. (to appear)

[2] [www.nupack.org](www.nupack.org)

[3] M.E. Fornace, N.J. Porubsky, and N.A. Pierce (2020). A unified dynamic programming framework for the analysis of interacting nucleic acid strands: enhanced models, scalability, and speed. [ACS Synth Biol, 9:2665-2678](https://pubs.acs.org/doi/abs/10.1021/acssynbio.9b00523), 2020. 
