### 使用相场方法建模多晶石墨烯

#### 1.建模PFC程序

源项目地址：https://github.com/petenez/pfc

参考文献引用：
```
Hirvonen et al., Physical Review B 94, 035414 (2016)
```

#### 2.安装MPI和FFTW-MPI

如果已经安装Intel Parallel Studio的话，可以直接使用mpicc（对应C编译器为系统gcc）；如果没有的话，可以安装openmpi或者mpich，这里推荐使用mpich：

```
# wget http://www.mpich.org/static/downloads/3.3.1/mpich-3.3.1.tar.gz
# tar -vxf mpich-3.3.1.tar.gz
# cd mpich-3.3.1/
# ./configure --prefix="/opt/devt/mpich3/"
# make -j 10
# make install
```

下载FFTW3源代码并安装：

```
# wget http://www.fftw.org/fftw-3.3.8.tar.gz
# tar -vxf fftw-3.3.8.tar.gz
# cd fftw-3.3.8
# ./configure --enable-mpi --prefix="/opt/devt/fftw3-mpi"
# make -j 10
# make install
```
