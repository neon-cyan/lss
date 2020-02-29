# Langevin system simulator

A small diffusion / Brownian motion simulation and data analysis suite. All the maths & graphs are based on [this paper](https://aapt.scitation.org/doi/suppl/10.1119/1.4772632)

## Set up

First, build the `lang` executable by running

```bash
make all
```
[Has been tested on gcc version 9.2.1]

For the data analysis it is recommended to set up a new python3/conda virtual environment, then install the required python packages with

```bash
pip install -r requirements.txt
```

Lastly start the jupyter server with


```bash
jupyter notebook
```
Lastly, navigate your way to `Data analysis.ipynb` - the notebook is self-explanatory

## Running a simulation

Some system parameters like the temperature[=300K], viscosity [=0.001Ns/m²], particle mass[=11pg] and radius[=1μm] are hard-coded. Any changes made in the `lang2.c` file will require re-compilation with `make all`.

[TODO] : Make more parameters changeable at runtime

Below are all the parameters that can be modified at runtime and the 3 simulations one can run

Free diffusion

```bash
./lang freediff [deltaT (double) = 0.1]  [NSteps (int) = 500] 
```

Brownian motion
```bash
./lang browndiff [deltaT (double) = 0.1]  [NSteps (int) = 500] 
```

Optical trap
```bash
./lang opttrap [deltaT (double) = 0.1]  [NSteps (int) = 500] [kx (double)]_[ky (double)]_[kz (double)]
```

## Output format

The program generates tab-separated lists of values for each experiment according to the following schema:

Free diffusion

| Position (x_i) | Scaled weight (W_i) |  Time [s] |
|----------------|:-------------------:|----------:|

Brownian motion

| Position with inertia (x_i) | Position without inertia (x_i) |  Time [tau] |
|-----------------------------|:------------------------------:|------------:|

Brownian motion in an optical trap

| Position (x_i) | Position (y_i) | Position (z_i) | Time [s]|
|----------------|:--------------:|:--------------:|--------:|

For a examples, take a look at the SAMPLES folder

## Data analysis

The jupyter notebook is organised such that every experiment has a corresponding class. The constructor can be called with the _generated_ file name automatically populating the class parameters. 
