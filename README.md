# Data Assimilation (DA) for Subsurface Biogeochemical Research (SBR)

[![N|Solid](https://upload.wikimedia.org/wikipedia/en/thumb/1/17/Pacific_Northwest_National_Laboratory_logo.svg/200px-Pacific_Northwest_National_Laboratory_logo.svg.png)](https://www.pnnl.gov/)

DA-SBR ueses ensemble Kalman filter (EnKF) based DA methods to assimilate observed field data at Hanford site (e.g. temperature, hydaulic heads...). This repository provides the entire workflow for the implementation of DA. Objectives include but not limited to:

  - Estimate dynamic hydraulic properties
  - Infer subdaily hydraulic exchange flux

# Contents
# 
| Workflow | Contents |
| ------ | ------ |
| 1. Preprocessing | procesing the raw temperatue data, hydraulic head data |
| 2. Input | model parameters, batch job configuration |
| 3. Processing | conducting data assimilation using different DA methods |
| 4. Postprocessing | display the results |


# Installation and Configuration

1. The workflow is written in [Jupyter notebook](http://jupyter.org/) which supports both Python and R. A recommended distribution of Jupyter notebook is [Anaconda](https://www.anaconda.com/download/).
    (1) To start Jupyter notebook on Mac/Linux after installation of Anaconda, typing the following command in terminal:
    ```sh
    jupyter notebook
    ```
    (2) To start Jupyter notebook on Windows, just click the desktop icon.
2. The numerical model for subsurface flow at Hanford site is built using [PFLOTRAN](http://www.pflotran.org/), a massively parallel reactive flow and transport model for describing surface and subsurface processes
3. The workflow is adapted to the supercomputers at [National Energy Research Scientific Computing Center (NERSC)](http://www.nersc.gov/)

Estimated flux:
[!alt text](https://www.flickr.com/photos/164278148@N08/shares/5Fn8B6)
