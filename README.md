# Latent RANSAC
See paper in https://arxiv.org/abs/1802.07045

## Depencencies

* Config++ (note - this DLL need to be in PATH)
* Lapack
* OpenGV
* Eigen (needed by OpenGV)


## Basis of our code

Since we based our on (and compared to) USAC, we forked this repository from cr333/usac-cmake who ported it to camke format. 

The original USAC source code (in branch 'usac') is available from:
http://www.cs.unc.edu/~rraguram/usac/

CMake/Mac porting by Christian Richardt (http://richardt.name).

FindOpenGV.cmake taken from  https://github.com/alicevision/AliceVision/blob/develop/src/cmake/FindOpenGV.cmake