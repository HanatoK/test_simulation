 
# WTM-eABF example with anisotropic diffusion

This is a prototypical molecular dynamic simulation program of anisotropic diffusion. Specifically, After compilation, the program will run simulations of an atom of 12 g/mol on the following Berezhkovskii-Szabo potential [1]:

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

The program uses Langevin thermostats and BAOAB integrator. Assuming $\gamma_x$ and $\gamma_y$ are the friction coefficients along $X$ and $Y$ of the thermostats, there are nine simulations:

1. Unbiased simulation with $\gamma_x = 10.0$ and $\gamma_y = 10.0$. The trajectory is dumped to `XYZ_10_10.traj`.
2. Unbiased simulation with $\gamma_x = 100.0$ and $\gamma_y = 10.0$. The trajectory is dumped to `XYZ_100_10.traj`.
3. Unbiased simulation with $\gamma_x = 10.0$ and $\gamma_y = 100.0$. The trajectory is dumped to `XYZ_10_100.traj`.
4. Two-dimensional WTM-eABF simulation along $X$ and $Y$ with $\gamma_x = 10.0$ and $\gamma_y = 10.0$. The trajectory is dumped to `XYZ_10_10_b.traj` and WTM-eABF related data are saved to files starting with `bias_10_10*`.
5. Two-dimensional WTM-eABF simulation along $X$ and $Y$ with $\gamma_x = 100.0$ and $\gamma_y = 10.0$. The trajectory is dumped to `XYZ_100_10_b.traj` and WTM-eABF related data are saved to files starting with `bias_100_10*`.
6. Two-dimensional WTM-eABF simulation along $X$ and $Y$ with $\gamma_x = 10.0$ and $\gamma_y = 100.0$. The trajectory is dumped to `XYZ_10_100_b.traj` and WTM-eABF related data are saved to files starting with `bias_10_100*`.
7. One-dimensional WTM-eABF simulation along the pathCV (PCV) [2] $s$ with $\gamma_x = 10.0$ and $\gamma_y = 10.0$. Results are saved to files starting with `PCV_bias_10_10*`.
8. One-dimensional WTM-eABF simulation along the pathCV (PCV) [2] $s$ with $\gamma_x = 100.0$ and $\gamma_y = 10.0$. Results are saved to files starting with `PCV_bias_100_10*`.
9. One-dimensional WTM-eABF simulation along the pathCV (PCV) [2] $s$ with $\gamma_x = 10.0$ and $\gamma_y = 100.0$. Results are saved to files starting with `PCV_bias_10_100*`.

The two components of PCV, $s(x,y)$ and $z(x,y)$, are defined as

$$
s(x,y)=\frac{1}{N-1}\frac{\sum_{i=0}^{N-1} i \exp\left(-\lambda \left((x-x_i)^2+(y-y_i)^2\right) \right)}{\sum_{i=0}^{N-1} \exp\left(-\lambda \left((x-x_i)^2+(y-y_i)^2\right) \right)},
$$

$$
z(x,y)=-\frac{1}{\lambda} \ln \left(\sum_{i=0}^{N-1} \exp\left( -\lambda (x-x_i)^2+(y-y_i)^2 \right)\right),
$$

where $N$ is the number of images of the string, and $(x_i, y_i)$ is the image position of $i$-th image. $s(x,y)$ and $z(x,y)$ represent the progress along the string and the distance perpendicular to it, respectively. All PCV simulations above use the string in `data/path_new3.txt` by default.

## Compilation

You need a C++20 compatible compiler and the `fmt` package. To compile, simply run the following commands in the root directory of the source code:
```
mkdir build
cd build/
cmake ../ -DCMAKE_BUILD_TYPE=release
make -j4
```

## Running

Just run the executables `TestUnbiased1`, `TestUnbiased2`, `TestUnbiased3`, `TestBiased1`, `TestBiased2`, `TestBiased3`, `TestPCVBiased1`, `TestPCVBiased2` and `TestPCVBiased3`. The logs are dumped to `stdout`.

After simulation, you can compile the `abf_integrate` from [Colvars](https://github.com/Colvars/colvars/tree/master/colvartools) and run
```
abf_integrate bias_100_10_300000000.czar.grad
abf_integrate bias_10_10_300000000.czar.grad
abf_integrate bias_10_100_300000000.czar.grad
```
to integrate the 2D PMFs. Copy the resultant PMFs to `scripts/` and run `./draw_fes_2D_2.py` to visualize them. For 1D PMFs along PCV $s$, it is easy to integrate them by the trapezoidal rule.

## Reference

[1] [Berezhkovskii, A.; Szabo, A. One-Dimensional Reaction Coordinates for Diffusive Activated Rate Processes in Many Dimensions. The Journal of Chemical Physics 2005, 122 (1), 014503.](https://doi.org/10.1063/1.1818091)

[2] [Branduardi, D.; Gervasio, F. L.; Parrinello, M. From A to B in Free Energy Space. The Journal of Chemical Physics 2007, 126 (5), 054103.]( https://doi.org/10.1063/1.2432340)