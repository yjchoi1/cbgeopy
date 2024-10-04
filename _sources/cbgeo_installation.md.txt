# CB-geo MPM installation
For the original installation guide, refer to the CB-geo website [here](https://mpm.cb-geo.com/#/). 
This document deal with the latest successful compilation on Stampede3 system on TACC.  

### Stampede3
```shell
cd $WORK
module load boost
module load hdf5
module load intel
# export LD_LIBRARY_PATH=$SWR_LD_LIBRARY_PATH:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=/work2/01197/semeraro/stampede3/VTK/lib64:$LD_LIBRARY_PATH

git clone https://gitlab.com/libeigen/eigen.git
cd $WORK && git clone https://github.com/KaHIP/KaHIP && \
   cd KaHIP && sh ./compile_withcmake.sh
cd $WORK && git clone https://github.com/cb-geo/mpm.git
git clone https://github.com/cb-geo/mpm-benchmarks.git benchmarks
export CC=icx
export CXX=icpx

# Set the HDF5 environment variables
export HDF5_ROOT=$TACC_HDF5_DIR
export HDF5_INCLUDE_DIRS=$TACC_HDF5_INC

cd mpm && mkdir build && cd build

cmake -DMPM_BUILD_TESTING=Off \
      -DBOOST_ROOT=$TACC_BOOST_DIR \
      -DBOOST_INCLUDE_DIRS=$TACC_BOOST_INC \
      -DCMAKE_BUILD_TYPE=Release \
      -DEIGEN3_INCLUDE_DIR=$WORK/eigen \
      -DKAHIP_ROOT=$WORK/KaHIP \
      -DHDF5_ROOT=$HDF5_ROOT \
      -DHDF5_INCLUDE_DIRS=$HDF5_INCLUDE_DIRS \
      -DVTK_ROOT=/work2/01197/semeraro/stampede3/VTK/lib64 \
      ..

make -j8
```
Note:
> * The compilation will step at around 45%, but MPM runs.
> * VTKs are not generated as outputs.
> * When set `-DMPM_BUILD_TESTING=Off`, the build goes to 100%, but still no VTK as outputs.  

### Frontera
Please refer to https://mpm.cb-geo.com/#/user/run/hpc/tacc-frontera.
```shell
# (on update)
```