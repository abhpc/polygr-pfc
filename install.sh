#! /bin/bash



wget http://www.fftw.org/fftw-3.3.10.tar.gz
tar -vxf fftw-3.3.10.tar.gz
cd fftw-3.3.10
module load mpi/2021.9.0 mkl/2023.1.0
./configure --enable-mpi --prefix="/home-beegfs/apps/devt/fftw-3.3.10mpi"
make -j 10
make install