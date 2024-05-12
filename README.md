# Modeling polycrystalline graphene for molecular dynamics simulations via phase field method <!-- omit in toc -->

Some codes in this program are forked from https://github.com/petenez/pfc. However, the author stated that he has leaved academia and had limited time to prepare these tools for publication. Since more codes have been added, I created this repository for easily generate initial configuration of polycrystalline graphene in LAMMPS data format.

# Table of Contents <!-- omit in toc -->
- [1 System requirement](#1-system-requirement)
  - [Install FFTW-MPI](#install-fftw-mpi)
- [2 Build phase field program](#2-build-phase-field-program)
- [3](#3)
    - [2.下载安装PFC程序](#2下载安装pfc程序)
      - [编译安装](#编译安装)
    - [3.输入文件解析](#3输入文件解析)


## 1 System requirement

### Install FFTW-MPI



```bash
# wget http://www.fftw.org/fftw-3.3.8.tar.gz
# tar -vxf fftw-3.3.8.tar.gz
# cd fftw-3.3.8
# ./configure --enable-mpi --prefix="/opt/devt/fftw3-mpi"
# make -j 10
# make install
```



如果已经安装Intel Parallel Studio的话，可以直接使用mpicc（对应C编译器为系统gcc）；如果没有的话，可以安装openmpi或者mpich，这里推荐使用mpich：

```bash
# wget http://www.mpich.org/static/downloads/3.3.1/mpich-3.3.1.tar.gz
# tar -vxf mpich-3.3.1.tar.gz
# cd mpich-3.3.1/
# ./configure --prefix="/opt/devt/mpich3/"
# make -j 10
# make install
```
安装好以后，记得将/opt/devt/mpich3/bin添加到系统的PATH中，方便后面使用。

下载FFTW3源代码并安装：




## 2 Build phase field program



## 3


#### 2.下载安装PFC程序

本项目使用Hirvonen等人开发的相场程序进行建模，源项目地址：https://github.com/petenez/pfc

使用该方法务必引用参考文献：
```
Hirvonen et al., Physical Review B 94, 035414 (2016)
```

##### 编译安装
```
$ git clone https://github.com/petenez/pfc
$ cd pfc/src
$ mpicc pfc.c -lfftw3_mpi -lfftw3 -lm -Ofast -Wall -I/opt/devt/fftw3-mpi/include -L/opt/devt/fftw3-mpi/lib -o pfc
```
编译成功后，在src目录下会生成可执行文件pfc。

#### 3.输入文件解析
