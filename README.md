# SCAN2CAD
AAU Bachelor Project

# Introduction
This project uses a RANSAC based appoach to construct a CAD model of an object from a point cloud.

# Build Instructions
The project uses code from the Point Cloud Library (https://pointclouds.org/). Made with the 1.10 version for Visual Studio 19 (available at https://github.com/PointCloudLibrary/pcl/releases). Due to time constraints on the project, the build process is slightly inconvenient.

Recommended build instructions:

1. Download and install the appropriate versions of Visual Studio and the PCL.
2. Clone the repository, build from Visual Studio
3. A copy of OpenNI2.dll must be provided for IO functionality in the PCL. The OpenNI2 SDK contains all needed dlls. (https://structure.io/openni)
4. The program should now build and run

# Variables
The test object is an I-beam, defined by three variables:

testbeam_position: A 3D vector (x, y, z) pointing to the start location for the beam
testbeam_size: A 3D vector (x, y, z) defining the size in each direction. 
testbeam_noiselevel: The max error introduced to the simulated scan points (evenly distributed with mean 0, range \[-noiselevel/2; +noiselevel/2\].

The simulated scan points of the object are even distributed on a 3D cubic grid.

Other than that, the relevant variables for the different parts are explained in the source code.
