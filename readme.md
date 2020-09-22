[![Double bubbles sans toil and trouble: discrete circulation-preserving vortex sheets for soap films and foams](http://www.cs.columbia.edu/cg/doublebubbles/title.jpg)](http://www.cs.columbia.edu/cg/doublebubbles/)

SoapFilm3D is an open source project for the physical simulation of soap films and bubbles. It is
cross-platform (Mac OS X, Linux, Windows, and more), and licensed under the Mozilla Public License
v. 2.0.

We would like to hear from you if you appreciate this work.

The repository contains the test program accompanying the paper by
[Da et al.](http://www.cs.columbia.edu/cg/doublebubbles/), including the multimaterial mesh-based
surface tracking library LosTopo. It also includes the curvature computation proposed in the article
by [Fei et al](https://arxiv.org/abs/1910.06402). 

The 

This program is built by standard procedures using CMAKE (http://www.cmake.org).

Dependencies
--------------------
The following external libraries are required:

* Eigen (http://eigen.tuxfamily.org)
* OpenGL and GLUT (http://www.opengl.org/resources/libraries/glut/)
* LAPACK (http://www.netlib.org/lapack/)
* BLAS (http://www.netlib.org/blas/)
* libPNG (https://libpng.sourceforge.io/)
* zLib (https://www.zlib.net/)
* GLEW (http://glew.sourceforge.net/)

On Mac OS X or Linux-based systems, most of the dependencies are either included, or can be easily
installed with Homebrew (https://brew.sh) or the APT package handling utility. 

On Windows you may need manually download and install some of them.

Compilation
-----------------
SoapFilm3D has been tested with Clang (under Mac OS X), and GCC 4.8+ (under Linux). Currently the
MSVC compiler is partially supported without the FMMTL module.

To compile SoapFilm3D, you'll need CMake on Mac OS X or Linux.

1. make a directory, say, `build`, with `mkdir build`, enter the *build* directory, type `cmake ..`
2. Optionally you can adjust the options with `ccmake .`
3. type `make` to compile the code. For speeding up the compilation process through multithreading
   you may use `make -j`.

On Windows:

1. open CMake-GUI, enter the correct directory for source code and build. Then click *Configure*,
   choose your installed version of the Microsoft Visual Studio.
2. after configuration you may find several libraries not found (with notifications of errors),
   check the *Advanced* box and *specify those missing header path and libraries manually*. For
   example, if Eigen is missing, then please specify the `EIGEN3_INCLUDE_DIR` to the path of directory
   we provided. For the ones we have not provided, you need to download and compile them, and then
   specify the missing directories to the path containing your headers or compiled libraries. Please
   make sure you have picked the libraries corresponding to the architecture you have selected (say,
   32-bit libraries for x86, and 64-bit libraries for x64).
3. click generate after fixing all missing variables to generate your Visual Studio solution.
4. open the Visual Studio solution and compile the code.
5. before running the demo, all the compiled dynamic linking libraries (DLLs) for your dependencies
   should be accessible from your PATH environment variable that can be changed in system settings,
   or you may simply copy them into your System32 (x64) or SysWOW64 (x86) directories.

Since collecting and compiling the dependencies could be tricky for Windows, we provided some
compiled libraries and DLLs for your convenient. You may find them in the zip file under the bin
folder, which also contains a compiled executable (under Visual Studio 2015 and Windows 10). 

Run the Demo
--------------------
Running the executables without command line arguments will display usage
information. All the data files for testing are located in the assets folder.
Since FMMTL is not supported on Windows, running can be slower than expected.

Please run the code using the project folder that contains everything (assets, code, etc.) as the
working directory. 

## Simulation Options

The software accepts different options that must be stored in a csv file with spaces as delimiter.

We present here a quick description for each possible options and their defaults. For
a more detailed documentation see [simulation\_options.md](simulation_options.md).

| Key                                         | Use                                                                              | Default        |
| ---                                         | ---                                                                              | -------        |
| `scene`                                     | Chose the scene to simulate (See [scenes.md](scenes.md))                         | "T1"           |
| `load-dir`                                  | Used for the scene `load` (See [scenes.md](scenes.md))                           | ""             |
| `load-increment`                            | Used for the scene `load` (See [scenes.md](scenes.md))                           | 1              |
| `time-step`                                 | The time step used in the simulation                                             | 0.01           |
| `simulation-time`                           | The total duration of the simulation                                             | 1.0            |
| `implicit-integration`                      | Enable implicit integration                                                      | 0              |
| `pbd-implicit`                              | Choses between Newton's Method (0) or Position Based Dynamics implicit solve (1) | 0              |
| `RK4-velocity-integration`                  | Do implicit integration using a Runge-Kutta method of order 4                    | 0              |
| `smoothing-coef`                            | Weight of the circulation smoothing                                              | 0.0            |
| `damping-coef`                              | Damping coefficient (`e^{-a} where a is the coefficient of the damping force)    | 1.0            |
| `sigma`                                     | The ratio of surface tension over density                                        | 1.0            |
| `regularization-coefficient`                | Controls the regularization parameter of the Biot-Savart law                     | .5             |
| `mean-curvature-computation`                | Choses the method used for computing the mean curvature                          | "average-area" |
| `fast-summation`                            | Chooses method for the computation of the Biot-Savart law                        | "fmmtl"        |
| `winding-beta`                              | The beta parameter of the Barnes-Hut method number method                        | 4.0            |
| `winding-expansion-order`                   | The order of the expansion used in the Barnes-Hut method                         | 2              |
| `fmmtl-theta`                               | The theta parameter of the Fast Multipole Method                                 | 0.5            |
| `fmmtl-minimum-cell-size`                   | Minimum number of elements in the octree's leaves of the Fast Multipole Method   | 128            |
| --------                                    | -------                                                                          | ----           |
| `output-png`                                | Enables output of rendering of the scene                                         | 1              |
| `output-png-every-n-frames`                 | Time step period of rendering output (0 is equivalent to 1)                      | 0              |
| `output-mesh`                               | Enables output of mesh binary files of the scene (in .rec)                       | 0              |
| `output-mesh-every-n-frames`                | Time step period of mesh binary file output (0 is equivalent to 1)               | 0              |
| `output-obj`                                | Enables output of the mesh text file output (in .obj)                            | 0              |
| `output-obj-every-n-frames`                 | Time step period of mesh text file output (0 is equivalent to 1)                 | 0              |
| ---------                                   | ---------                                                                        | -----          |
| `remeshing-resolution`                      | Target edge length                                                               | 0.1            |
| `remeshing-iterations`                      | Number of remeshing iterations to do at each step                                | 1              |
| `initial-remeshing-iterations`              | Number of remeshing iterations to do before starting the simulation              | 1              |
| ----------                                  | ------------                                                                     | -------        |
| `lostopos-collision-epsilon-fraction`       | Lostopos parameter (See LosTopos)                                                | 1e-4           |
| `lostopos-merge-proximity-epsilon-fraction` | Lostopos parameter (See LosTopos)                                                | 0.02           |
| `lostopos-perform-smoothing`                | Lostopos parameter (See LosTopos)                                                | 0              |
| `lostopos-max-volume-change-fraction`       | Lostopos parameter (See LosTopos)                                                | 1e-4           |
| `lostopos-min-triangle-angle`               | Lostopos parameter (See LosTopos)                                                | 3.0            |
| `lostomos-max-triangle-angle`               | Lostopos parameter (See LosTopos)                                                | 177.0          |
| `lostopos-large-triangle-angle-to-split`    | Lostopos parameter (See LosTopos)                                                | 160            |
| `lostopos-min-triangle-area-fraction`       | Lostopos parameter (See LosTopos)                                                | 0.02           |
| `lostopos-t1-transisition-enabled`          | Lostopos parameter (See LosTopos)                                                | 1              |
| `lostopos-smooth-subdivision`               | Lostopos parameter (See LosTopos)                                                | 0              |
| `lostopos-allow-non-manifold`               | Lostopos parameter (See LosTopos)                                                | 1              |
| `lostopos-allow-topology-changes`           | Lostopos parameter (See LosTopos)                                                | 1              |
| --------                                    | ----------                                                                       | -----          |
| `mesh-size-n`                               | Depends on the scene (See [scenes.md](scenes.md))                                | 2              |
| `mesh-size-m`                               | Depends on the scene (See [scenes.md](scenes.md))                                | 2              |
| `foam-burst-interval`                       | Option of the scene `foam` (See [scenes.md](scenes.md)                           | 20.0           |
| `foam-burst-start`                          | Option of the scene `foam` (See [scenes.md](scenes.md)                           | 10.0           |

## Keybindings

| Key            | Effect                                                                                        |
| -------        | ------                                                                                        |
| `q`/`Q`/Escape | quits the simulation                                                                          |
| ` `            | Starts/Stops the simulation                                                                   |
| `s`/`S`        | Make a step in the simulation                                                                 |
| `m`            | Change the rendering mode to the next one                                                     |
| `M`            | Change the rendering mode the the previous one                                                |
| `v`            | Start vertices selection                                                                      |
| `V`            | Stop vertices selection                                                                       |
| `e`            | Start edge selection                                                                          |
| `E`            | Stop edge selection                                                                           |
| `f`            | Start face selection                                                                          |
| `F`            | Stop face selection                                                                           |
| `n`/`N`        | Show selected primitive information                                                           |
| `]`            | Load the 1st mesh forward from the current one (Only works when the scene type is `load`)     |
| `[`            | Load the 1st mesh backward from the current one (Only works when the scene type is `load`)    |
| `}`            | Load the 10th mesh forward from the current one (Only works when the scene type is `load`)    |
| `{`            | Load the 10th mesh backward from the current one (Only works when the scene type is `load`)   |
| `.`            | Load the 100th mesh forward from the current one (Only works when the scene type is `load`)   |
| `,`            | Load the 100th mesh backward from the current one (Only works when the scene type is `load`)  |
| `>`            | Load the 1000th mesh forward from the current one (Only works when the scene type is `load`)  |
| `<`            | Load the 1000th mesh backward from the current one (Only works when the scene type is `load`) |
| `o`            | Save current mesh to `mesh.obj`                                                               |
| `a`/`A`        | Automatically load the next mesh                                                              |

Citations
--------------------
Da, Fang, et al. "Double bubbles sans toil and trouble: Discrete circulation-preserving vortex sheets for soap films and foams." ACM Transactions on Graphics (TOG) 34.4 (2015): 149. (http://www.cs.columbia.edu/cg/doublebubbles/)

Da, Fang, Christopher Batty, and Eitan Grinspun. "Multimaterial mesh-based surface tracking." ACM Trans. Graph. 33.4 (2014): 112-1. (http://www.cs.columbia.edu/cg/multitracker/)

Fei, Yun (Raymond), et al. "Addressing Troubles with Double Bubbles: Convergence and Stability at Multi-Bubble Junctions." arXiv preprint arXiv:1910.06402 (2019). (https://arxiv.org/abs/1910.06402)
