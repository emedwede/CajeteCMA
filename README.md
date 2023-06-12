# A Graph Grammar Simulator for the Cortical Microtubule Array (CMA)

A serial implementation for a graph grammar simulator featured in the paper titled, "Approximate Simulation of Cortical Microtubule Models using Dynamical Graph Grammars". This repo is intended to be a proof of concept for the approximate simulation algorithm and not for production use. The backend is powered by an early version of YAGL: Yet another Graph Library.


# Requirements
The minimum CXX version required to compile is C++17. The code has been tested on Pop!_OS 20.04 LTS, but most versions of linux should be able to compile and run the code. We also use the state of the art ODE solver [SUNDIALS](https://github.com/LLNL/sundials), and a build of the latest version is required to be provided in the build script for CMAKE to find. Every other requirement is packaged in. The build instructions are fairly simple, but require a little bit of work. First build [SUNDIALS](https://github.com/LLNL/sundials), then read and modify the builder script in the scripts directory and simply run ./scripts/builder.sh The build will be placed out of source wherever the user defines.

# Usage and Output

To use the main simulation code, navigate to your build directory and then to the simulators folder. Use the command `mt_dgg_simulator configuration/[experiment_name].json` to run the experiment with the defined settings. Feel free to modify parameters in the config files, but there is little error checking so changes may lead to undefined behavior. 

The results of the simulations are saved as graphs in the [VTK](https://kitware.github.io/vtk-examples/site/VTKFileFormats/) file format. The system graph file, simulation_step_xx.vtu, is saved at the checkpoint frequency defined in the config file. The cell complex for the simulation domain is saved as factory_geoplex.vtu. Both types of data files are saved in [experiment_name]_results folder locally created in the executable directory. 

**Please be aware that simulation step files can get rather large, and frequent checkpointing can add up! So, saving several thousand steps for several simulations is a bad idea unless you have the storage to spare!**

To visualize any results, we highly recommend using [Paraview](https://www.paraview.org/download/) and loading the simulation step files and geoplex files twice. On the first load you can apply a surface with edges filter and on the second load you can apply the points filter. Then, color the points by node_type.

# Experiments and Results

In the paper we feature three experiments. The first, experiment 1, tests network formation. Experiment 2 tests alignment. Experiment 3 tests runtime performance on five different domain subdivisions: 1x1, 2x1, 2x2, 4x4, and 8x8. The configuration files for experiments can be found in the simulators/configurations folder. 

We have also included the jupyter notebooks to produce most of the plots in the paper. To run these, python 3 is required along with numpy, matplotlib, and scipy. These can be found in the folder simulators/notebooks.
