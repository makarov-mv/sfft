Implementation of dimension independent sparse fourier transform
# Requirements
To compile and run tests you will need:
* CMAKE >=3.8
* FFTW3

To make plots you will additionally need:
* python 3.5
* matplotlib

# Compilation
Open repository directory in terminal and run:
```
mkdir build
cd build
cmake -DCMAKE_BUILD_TYPE=Release ..
make
```

If you want to add `sfft v.1-v.2` to graphs, compile them by running:

```
cd sfft-v1v2
mkdir build
cd build
cmake -DCMAKE_BUILD_TYPE=Release ..
make
```

Now you can run one of the `test_*` executables in `build` directory.

# Making graphs

After compilation, you can run `plot_experiments_vs_sparsity.sh`,
`plot_experiments_vs_fftw.sh` and `plot_experiments_vs_sfft.sh` to make graphs. It may take several hours.
Don't forget to run `run_sfft.sh` in `sfft-v1v2` directory before running `plot_experiments_vs_sfft.sh`.