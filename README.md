# DFT-D4 damping parameter optimization

This project provides a driver for optimization of damping parameters in the DFT-D4 method.


## Usage

This project is built with the Fortran package manager ([fpm](https://github.com/fortran-lang/fpm)).
To run a parameter optimization out of this source directory use

```
git clone https://github.com/dftd4/dftd4-fitset
fpm run --profile release --flag -fopenmp -- -C dftd4-fitset ./data.csv
```

The data set must contain the dispersion energies for each reaction in Hartree:

```csv
S22x5/01-0.9, S22x5/01-A, S22x5/01-B, 1.0007611865e-03
S22x5/01-1.0, S22x5/01-A, S22x5/01-B, 1.5228237266e-03
S22x5/01-1.2, S22x5/01-A, S22x5/01-B, 1.6586059147e-03
S22x5/01-1.5, S22x5/01-A, S22x5/01-B, 1.2297590834e-03
S22x5/01-2.0, S22x5/01-A, S22x5/01-B, 6.2420992500e-04
...
```

For more information checkout the project help page

```
fpm run -- --help
```
