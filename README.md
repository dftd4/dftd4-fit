# DFT-D4 damping parameter optimization

This project provides a driver for optimization of damping parameters in the DFT-D4 method.

## Installation

To use this project the following dependencies are required

- [`dftd4`](https://github.com/dftd4/dftd4) version 3.0.0 or newer,
  for evaluation of the DFT-D4 dispersion correction
- [`nlopt`](https://nlopt.readthedocs.io) version 2.0.0 or newer,
  for the optimization of the damping parameters
- [`nlopt-f`](https://github.com/grimme-lab/nlopt-f),
  to provide Fortran bindings for `nlopt`
- [`minpack`](https://github.com/fortran-lang/minpack),
  for the optimization of the damping parameters
- [`mctc-lib`](https://github.com/grimme-lab/mctc-lib),
  for reading geometry inputs when creating the data set

### Meson

Create a new build with

```
meson setup _build --prefix ~/.local
```

The Fortran compiler can be adjusted by setting the `FC` environment variable.
Meson will find installed dependencies automatically or fetch them automatically if they are not available.

To build the project use

```
meson compiler -C _build
```

Finally, you install the binary with

```
meson install -C _build
```

### Fortran package manager

This project is built with the Fortran package manager ([fpm](https://github.com/fortran-lang/fpm)).
Fpm will handle the dependencies for `dftd4`, `minpack`, `nlopt-f` and `mctc-lib`, only `nlopt` has to be provided on the system.

To create the binaray use

```
fpm install --profile release --flag -fopenmp
```

This will install `dftd4-fit` to `~/.local/bin` on Unix and `%APPDATA%\local\bin` on Windows.
Make sure the respective directories are in your `PATH` or adjust the installation prefix with the `--prefix` option.

## Usage

```
git clone https://github.com/dftd4/dftd4-fitset
dftd4-fit -C dftd4-fitset ./data.csv
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

For an overview over all command line arguments use the `--help` argument or checkout the [`dftd4-fit(1)`](man/dftd4-fit.1.adoc) manpage.

## License

This project is free software: you can redistribute it and/or modify it under
the terms of the Lesser GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This project is distributed in the hope that it will be useful,
but without any warranty; without even the implied warranty of
merchantability or fitness for a particular purpose. See the
Lesser GNU General Public License for more details.

Unless you explicitly state otherwise, any contribution intentionally
submitted for inclusion in this project by you, as defined in the
Lesser GNU General Public license, shall be licensed as above, without any
additional terms or conditions.
