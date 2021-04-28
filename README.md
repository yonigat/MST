# MST
MST - multiplexed scattering tomography performs 3D reconstruction of mediums using a scattering tomography algorithm while multiplexing sources.
Using this code the user can determine the scanner setup (i.e. the formation of detectors) and the sources location, orientation and FOV.
The code is devided into 2 parts:
1. Python code (located in scattering_tomography_src) where the phantom, setup and optimization are defined and run.
2. Geant4 simulation to perform forward model and calculate the gradient of the inverse problem. It can be run alone or through the Python part.

##  Setup
This section contains explantions on how to setup the enviroment to run the code.

### Conda enviroment
Set up the conda enviroment by:
1. Install conda [from here](https://docs.conda.io/projects/conda/en/latest/user-guide/install/linux.html).
2. In terminal: `conda env create -f environment.yml`

### Installing Geant4
1. Download files from [Geant4 site](http://geant4.web.cern.ch/support/download).
2. Place downloaded file in new directory (~/Geant4_workspace for example) 
3. Extract file using `tar -xvf geant4.10.06.p01.tar.gz`. 
4. Create 2 directories: geant4.10.06.p01-build and geant4.10.06.p01-install 
5. `cd` into build dir 
6. Run command: `cmake -DCMAKE_BUILD_TYPE=RelWithDebInfo -DGEANT4_BUILD_MULTITHREADED=ON -DGEANT4_INSTALL_DATA=ON -DGEANT4_USE_QT=ON -DGEANT4_USE_RAYTRACER_X11=ON -DGEANT4_BUILD_TLS_MODEL=global-dynamic -DCMAKE_INSTALL_PREFIX=/home/yonatangat/Geant4_workspace/geant4.10.06.p01-install /home/yonatangat/Geant4_workspace/geant4.10.06.p01`
7. Run command: `make -jN` (N number of threads).
8. Run connand: `make install`

### Building G4 simulator
The dir `ring_sim` contains only the source code for the radiative transfer simulator. To build it we use cmake. The code is built such that it can
be imported as project into Eclipse where you can debug it.
1. In the repository dir: `mkdir ring_sim_build`.
2. `cd` into `mkdir ring_sim_build`.
3. Run `cmake  -G"Eclipse CDT4 - Unix Makefiles" -DGeant4_DIR=/path/to/G4/install/dir ../ring_sim -DCMAKE_BUILD_TYPE=Debug`.
4. Run command: `make -jN` (N number of threads).

## Code Structure

