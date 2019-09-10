Implementation of dimension independent sparse fourier transform
# Requirements
To compile and run tests you will need:
- CMAKE >=3.8
- FFTW3

To make plots you will additionally need:
- python 3.5
- matplotlib

# Compilation
Open repository directory in terminal and run:
```
mkdir build
cd build
cmake -DCMAKE_BUILD_TYPE=Release ..
make
```

Now you can run one of the `test_*` executables in `build` directory.

# Making graphs

After compilation, you can use `run_experiment.py` to make graphs. You can use calls in `plot_experiments.sh`
as examples.

Run `plot_experiments.sh` to plot several basic graphs.