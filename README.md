# Data Assimilation (DA) for Subsurface Biogeochemical Research (SBR)

[![N|Solid](https://upload.wikimedia.org/wikipedia/en/thumb/1/17/Pacific_Northwest_National_Laboratory_logo.svg/200px-Pacific_Northwest_National_Laboratory_logo.svg.png)](https://www.pnnl.gov/)

DA-SBR ueses ensemble Kalman filter (EnKF) based DA methods to assimilate observed field data at Hanford site (e.g. temperature, hydaulic heads...). This repository provides the entire workflow for the implementation of DA. Objectives include but not limited to:

  - Estimate dynamic hydraulic properties
  - Infer subdaily hydraulic exchange flux

# Contents
# 
| Workflow | Contents |
| ------ | ------ |
| 1. [Pre-Processing](https://github.com/lovingckw/DA-SBR/blob/de59203c7399c596e7b24f05a147e28e63c7e23f/pre-processing.ipynb) | 1.1 Read HDF5 file generated by PFLOTRAN |
|                   | 1.2 Extract geometry and properties |
|                   | 1.3 Generate inputs for 1D hydro-thermal simulation |
| 2. [Data Assimilation](https://github.com/lovingckw/DA-SBR/blob/master/Data%20Assimilation.ipynb) | 2.1 Data input |
|                      | 2.2 PFLOTRAN input configuration |
|                      | 2.3 EnKF implementation |
| 3. [Post-Processing](https://github.com/lovingckw/DA-SBR/blob/de59203c7399c596e7b24f05a147e28e63c7e23f/post-processing.ipynb) | 3.1 Plot hydraulic conductivity  |
|                    | 3.2 Plot the vertical hydraulic exchange flux (HEF) |
|                    | 3.3 Compare the calculated HEF with that from temperature method|


# Installation and Configuration

1. The workflow is written in [Jupyter notebook](http://jupyter.org/) which supports both Python and R. A recommended distribution of Jupyter notebook is [Anaconda](https://www.anaconda.com/download/).
  (1) To start Jupyter notebook on Mac/Linux after installation of Anaconda, typing the following command in terminal:
    ```sh
    jupyter notebook
    ```
    (2) To start Jupyter notebook on Windows, just click the desktop icon. 
2. The numerical model for subsurface flow at Hanford site is built using [PFLOTRAN](http://www.pflotran.org/), a massively parallel reactive flow and transport model for describing surface and subsurface processes.
3. The workflow is adapted to the supercomputers at [National Energy Research Scientific Computing Center (NERSC)](http://www.nersc.gov/).

![Estimated flux](https://github.com/lovingckw/DA-SBR/blob/master/Doc/temp/Picture1.jpg)
![Estimated perm](https://github.com/lovingckw/DA-SBR/blob/master/Doc/temp/Picture4.jpg)
