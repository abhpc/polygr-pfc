# Modeling polycrystalline graphene for molecular dynamics simulations via phase field method <!-- omit in toc -->

Some codes in this program are forked from https://github.com/petenez/pfc. However, the author stated that he has leaved academia and had limited time to prepare these tools for publication. Since more codes have been added, I created this repository for easily generate initial configuration of polycrystalline graphene in LAMMPS data format.

# Table of Contents <!-- omit in toc -->
- [1 System requirements](#1-system-requirements)
  - [1.1 Install FFTW-MPI](#11-install-fftw-mpi)
  - [1.2 Install OpenJDK 13.0.1](#12-install-openjdk-1301)
  - [1.3 OVITO 2.9.0](#13-ovito-290)
- [2 Build phase field program](#2-build-phase-field-program)
- [3 Example](#3-example)


## 1 System requirements

These packages below are required for running this program:
```bash
mkl/2023.1.0  mpi/2021.9.0  gcc/7.5.0  openjdk/13.0.1  ovito/2.9.0
```

### 1.1 Install FFTW-MPI

```bash
wget http://www.fftw.org/fftw-3.3.10.tar.gz
tar -vxf fftw-3.3.10.tar.gz
cd fftw-3.3.10
module load mpi/2021.9.0
./configure --enable-mpi --prefix="/home-beegfs/apps/devt/fftw-3.3.10mpi"
make -j 10
make install
```
### 1.2 Install OpenJDK 13.0.1

### 1.3 OVITO 2.9.0




## 2 Build phase field program

```bash
git clone https://github.com/petenez/pfc
cd pfc/src
mpicc pfc.c -lfftw3_mpi -lfftw3 -lm -Ofast -Wall -I/opt/devt/fftw3-mpi/include -L/opt/devt/fftw3-mpi/lib -o pfc
```

## 3 Example