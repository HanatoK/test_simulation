 
# WTM-eABF example with anisotropic diffusion

This is a prototypical molecular dynamic simulation program of anisotropic diffusion. Specifically, After compilation, the program will run simulations of an atom of 12 g/mol on the following Berezhkovskii-Szabo potential:

$$
\beta U(x,y) = \beta U(x) + \Omega^2 (x-y)^2/2,
$$

where $U(x)$ is defined as

$$
\beta U(x) =
\begin{cases}
-\omega^2 x_0^2 / 4 + \omega^2 (x+x_0)^2/2, &x < -x_0 / 2\\
-\omega^2 x^2 / 2, &-x_0 / 2 < x < x_0/2\\
-\omega^2 x_0^2 / 4 + \omega^2 (x-x_0)^2/2, & x_0 / 2 < x,
\end{cases}
$$

where $x_0 = 2.2 \mathrm{\mathring{A}}$, $\omega^2 = 4 \mathrm{\mathring{A}}^{-2}$ and $\Omega=1.01\omega^2$. The temperature is 300 K.

The program uses Langevin thermostats and BAOAB integrator. Assuming $\gamma_x$ and $\gamma_y$ are the friction coefficients along $X$ and $Y$ of the thermostats, there are six simulations:

1. Unbiased simulation with $\gamma_x = 10.0$ and $\gamma_y = 10.0$. The trajectory is dumped to `XYZ_10_10.traj`.
2. Unbiased simulation with $\gamma_x = 100.0$ and $\gamma_y = 10.0$. The trajectory is dumped to `XYZ_100_10.traj`.
3. Unbiased simulation with $\gamma_x = 10.0$ and $\gamma_y = 100.0$. The trajectory is dumped to `XYZ_10_100.traj`.
4. WTM-eABF simulation with $\gamma_x = 10.0$ and $\gamma_y = 10.0$. The trajectory is dumped to `XYZ_10_10_b.traj` and WTM-eABF related data are saved to files starting with `bias_10_10*`.
5. WTM-eABF simulation with $\gamma_x = 100.0$ and $\gamma_y = 10.0$. The trajectory is dumped to `XYZ_100_10_b.traj` and WTM-eABF related data are saved to files starting with `bias_100_10*`.
6. WTM-eABF simulation with $\gamma_x = 10.0$ and $\gamma_y = 100.0$. The trajectory is dumped to `XYZ_10_100_b.traj` and WTM-eABF related data are saved to files starting with `bias_10_100*`.

## Compilation

You need a C++20 compatible compiler and the `fmt` package. To compile, simply run the following commands in the root directory of the source code:
```
mkdir build
cd build/
cmake ../ -DCMAKE_BUILD_TYPE=release
make -j4
```

## Running

Just run `./TestSimulation`. The log is dumped to `stdout`. After simulation, you can compile the `abf_integrate` from [Colvars](https://github.com/Colvars/colvars/tree/master/colvartools) and run
```
abf_integrate bias_100_10.czar.grad
abf_integrate bias_10_10.czar.grad
abf_integrate bias_10_100.czar.grad
```
to integrate the PMFs. Copy the resultant PMFs to `scripts/` and run `./draw_fes_2D_2.py` to visualize them.