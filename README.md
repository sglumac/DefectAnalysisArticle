### DefectAnalysisArticle

The code for the article Defect Analysis of a Non-iterative Co-simulation.

The algorithm code and code used to produce figure data is available in *src/Experiments*.
The models used in the experiments are implemented according to [FMI 2.0 standard](https://fmi-standard.org/)
- *src/OscillatorOmega2Tau*
- *src/Tau2Oscillator*

#### Quickstart

Run *produce_figures.bat*. This script will
- compile Experiments.exe
- and run it to produce all figures
The figures will be stored in *pdf* files found in the *out* folder.

Prerequisites:
- [Python with Numpy and Matplotlib](https://anaconda.org/anaconda/python)
- [CMake](https://cmake.org/)

The code has been tested on Windows with Visual Studio 2019.
