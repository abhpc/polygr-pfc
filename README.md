# DCU 加速的多晶石墨烯 PFC 程序

本项目用相场晶体方法（PFC）生成多晶石墨烯，并转换为 LAMMPS 可读取的原子初始构型。程序基于 [Petri Hirvonen 的 PFC 工具](https://github.com/petenez/pfc)，保留原物理方程、离散算法和双精度计算。

计算分为两个阶段：`step1` 生成并生长晶粒，`step2` 读取第一阶段末态进行结构弛豫。DCU 负责 FFT 和场更新；原子转换沿用 Java Coordinator 和 OVITO，灰度图片由 Python 导出。

![约100×100 nm、400个初始晶核的多晶石墨烯生长动画](docs/images/step1-growth-preview-white.gif)

展示案例：约 100×100 nm、400 个初始晶核的第一阶段生长过程，从第 0 步演化至第 10000 步，每 1000 步取一帧，共 11 帧循环播放。

编译与运行示例使用 Linux、Intel MPI 和单节点 4 张 DCU。使用资源调度系统时，按所在站点的规则申请计算资源。Slurm 部署示例见 [Host.md](Host.md)。

## 1. 获取源码并设置环境

需要 GCC、GNU Make、Intel MPI、DTK/HIP 和 hipFFT。后处理需要 Java 13、OVITO 2.9.0，以及 Python 3.11、NumPy 和 Pillow。已验证的工具链为 GCC 13.1.0、Intel MPI 2021.9.0、DTK 26.04，设备为 `gfx906`。

```bash
git clone https://github.com/abhpc/polygr-pfc.git
cd polygr-pfc
export PFC_ROOT="$PWD"
export FFTW_PREFIX="$PFC_ROOT/.local/fftw-3.3.10"
```

先加载或初始化已安装的编译器、MPI 和 DTK。没有使用环境模块时，按实际安装位置设置下面的路径；`/path/to/...` 需要替换：

```bash
export I_MPI_ROOT=/path/to/intel/mpi
export DTKROOT=/path/to/dtk
export HIP_PATH="$DTKROOT/hip"
export PATH="$I_MPI_ROOT/bin:$HIP_PATH/bin:$DTKROOT/bin:$PATH"
export LD_LIBRARY_PATH="$FFTW_PREFIX/lib:$I_MPI_ROOT/lib/release:$DTKROOT/lib:$HIP_PATH/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
```

Makefile 使用 Intel MPI 的 `include` 和 `lib/release` 布局。其他安装布局需对应修改 `MPI_FLAGS`；更换 MPI 实现时，FFTW 和最终程序必须使用同一 MPI ABI。`PFC_GPU_ARCH` 应与实际 DCU 型号一致。

## 2. 编译并安装 FFTW

已有兼容的 MPI 版 FFTW 时，将 `FFTW_PREFIX` 指向该安装目录，并更新 `LD_LIBRARY_PATH`，即可跳过本步。

```bash
mkdir -p "$PFC_ROOT/.local/sources" "$PFC_ROOT/.build"
curl --fail --location https://www.fftw.org/fftw-3.3.10.tar.gz \
  -o "$PFC_ROOT/.local/sources/fftw-3.3.10.tar.gz"
tar -xzf "$PFC_ROOT/.local/sources/fftw-3.3.10.tar.gz" -C "$PFC_ROOT/.build"
cd "$PFC_ROOT/.build/fftw-3.3.10"
CFLAGS='-O3 -march=native -ffp-contract=off -Wall' MPICC=mpicc CC=gcc \
  ./configure --prefix="$FFTW_PREFIX" --enable-mpi --enable-shared \
  --enable-sse2 --enable-avx --enable-avx2
make -j 8
make install
cd "$PFC_ROOT"
```

上述 SIMD 选项对应已验证的 x86_64 CPU。其他 CPU 架构需按硬件调整 FFTW 配置。

## 3. 编译 PFC

```bash
cd "$PFC_ROOT"
export PFC_GPU_ARCH=gfx906
make -j 8 bin/polygrpfc
```

生成的程序为 `bin/polygrpfc`。主机代码使用 `-O3 -march=native -ffp-contract=off -Wall`，设备使用 FP64，不启用 fast-math 或浮点收缩。程序链接本地 `hipfft` 和 MPI 版 FFTW，不混用其他 MPI ABI 的 `hipfftMp`。

需要编译原 CPU 对照程序时，执行：

```bash
make bin/polygrpfc_cpu_ref
```

该目标从固定提交 `c86efb51cf03086fb98897ab8abfe2704a62096e` 提取原始 `src/pfc.c`，生成 `bin/polygrpfc_cpu_ref`。因此需要保留包含该提交的 Git 历史。

## 4. 准备 100×100 nm、400 个初始晶核的输入

使用新的运行目录，避免与已有结果重名。下面使用 4240² 网格、固定随机种子、每阶段 10000 步，每 1000 步保存一次：

```bash
mkdir -p "$PFC_ROOT/runs"
mkdir "$PFC_ROOT/runs/graphene-100nm-400"
cd "$PFC_ROOT/runs/graphene-100nm-400"

sed -e 's@xabc@34789847@g;s@WW@4240@g;s@HH@4240@g;s@NN@400@g;s@RR@2@g' \
  "$PFC_ROOT/example/step1.in" > step1.in
sed -e 's@WW@4240@g;s@HH@4240@g;s@step1-t-10000.dat@step1-10000.dat@g' \
  "$PFC_ROOT/example/step2.in" > step2.in
```

`400` 是初始晶核数，晶核半径为 2 个晶格常数；生长后保留的晶粒数可能变化。沿用原转换尺度 `0.7 0.7 7.3 2.46` 时，4240² 网格的名义尺寸约为 100.017534×100.017534 nm。

原输入标签及含义保持不变：

| 标签 | 含义 |
|---|---|
| `S` | 随机种子 |
| `O` | 科学日志间隔、密度保存间隔 |
| `A` | 网格宽、高 |
| `I` | 初始化方式；第二阶段用 `I 2` 读取第一阶段结果 |
| `M` | 原物理模型参数 |
| `R` | 步数、网格间距、时间步长、盒优化间隔 |

修改步数时，同步修改两份输入的 `R` 行、第二阶段重启文件名，以及后处理命令中的末态文件和绘图区间。运行目录中的输入可直接编辑，源码目录的模板保持原样。

## 5. 运行两阶段 PFC

单节点使用 32 个 MPI rank 和 4 张可见 DCU：

```bash
export OMP_NUM_THREADS=1
export I_MPI_PIN=1
export I_MPI_PIN_DOMAIN=core

mpirun -bootstrap fork -np 32 -ppn 32 "$PFC_ROOT/bin/polygrpfc" step1 \
  > step1.stdout 2> step1.stderr
```

确认第一阶段成功完成后，运行第二阶段：

```bash
mpirun -bootstrap fork -np 32 -ppn 32 "$PFC_ROOT/bin/polygrpfc" step2 \
  > step2.stdout 2> step2.stderr
```

程序读取当前目录的 `step1.in` 和 `step2.in`；第二阶段读取 `step1-10000.dat`。运行过程中保留原六列科学输出，并记录每个实际推进步的计算时间。

这里的 `-bootstrap fork`、`-ppn` 是 Intel MPI 参数。当前验证保留 32 个真实科学 rank 的初始化和局部统计，每卡由一个执行进程聚合成员的 FFT 与场更新。不要仅为匹配卡数把进程数改成 4，也不要给每个 rank 独占绑定一张卡。

跨节点时，已验证的配置为两节点、每节点 16 rank 和 4 张卡，总计 32 rank。将 `node-a`、`node-b` 替换为可用节点名；两个节点需要相同的运行环境，并能以相同路径访问程序和运行目录：

```bash
mpirun -hosts node-a,node-b -np 32 -ppn 16 "$PFC_ROOT/bin/polygrpfc" step1 \
  > step1.stdout 2> step1.stderr
```

第一阶段成功完成后，再运行：

```bash
mpirun -hosts node-a,node-b -np 32 -ppn 16 "$PFC_ROOT/bin/polygrpfc" step2 \
  > step2.stdout 2> step2.stderr
```

跨节点命令用于新的运行目录，替代上面的单节点命令。节点启动和环境传播按本站 Intel MPI 配置处理；主机与进程布局参数见 [Intel MPI 文档](https://www.intel.com/content/www/us/en/developer/articles/technical/controlling-process-placement-with-the-intel-mpi-library.html)。

## 6. 转换 LAMMPS data

在 DCU 运行目录中，使用 Java Coordinator 提取第二阶段末态的原子坐标和邻居信息：

```bash
java -jar "$PFC_ROOT/java/coordinator.jar" step2-10000.dat 4240 4240 \
  0.7 0.7 7.3 2.46 step2-10000.xy step2-10000.nh \
  > conversion.stdout 2> conversion.stderr
```

整理为 XYZ，再使用 OVITO 2.9.0 自带的 `ovitos` 导出 LAMMPS data：

```bash
awk 'END {print NR}' step2-10000.xy > atoms.count
{
  cat atoms.count
  printf '\n'
  awk '{print $1+1,1,$2,$3,0}' step2-10000.xy
} > graphene.xyz

ovitos "$PFC_ROOT/source/xyz2lmp.py" >> conversion.stdout 2>> conversion.stderr
sed -i 's@0.0 0.0 zlo zhi@-10.0 10.0 zlo zhi@g' graphene.data
```

`xyz2lmp.py` 从当前目录读取 `graphene.xyz`，写出 `graphene.data`；最后一条命令沿用原流程，将 z 方向盒边界设为 −10～10 Å。Java、XYZ 和 OVITO 转换会写固定名称的文件，应在独立运行目录执行。

保留现有 `java/coordinator.jar`，不要直接用接口不同的新版 `Coordinator.java` 重编替换。其他 OVITO 版本需要适配其 Python API 和输出格式。

## 7. 批量导出灰度 PNG

创建独立 Python 环境，安装已验证版本的 NumPy 和 Pillow：

```bash
python3.11 -m venv "$PFC_ROOT/.venv"
"$PFC_ROOT/.venv/bin/python" -m pip install numpy==1.24.3 Pillow==9.4.0
```

回到 DCU 运行目录导出两个阶段的全部保存帧：

```bash
cd "$PFC_ROOT/runs/graphene-100nm-400"
"$PFC_ROOT/.venv/bin/python" "$PFC_ROOT/source/plotter.py" \
  'step1-#.dat' 'step1-#.png' 4240 4240 0 1000 10000
"$PFC_ROOT/.venv/bin/python" "$PFC_ROOT/source/plotter.py" \
  'step2-#.dat' 'step2-#.png' 4240 4240 0 1000 10000
```

`#` 按 `起始步、间隔、终止步` 展开，路径模板需要引号。每阶段导出 11 张图片；图片逐帧归一化为 8 位灰度并上下翻转，常量场输出黑图。已有图片不会被覆盖。

单帧导出的命令为：

```bash
"$PFC_ROOT/.venv/bin/python" "$PFC_ROOT/source/plotter.py" \
  step2-10000.dat step2-final-10000.png 4240 4240
```

求解器、普通 Python 和 OVITO 自带 Python 使用各自的环境；先完成求解，再切换后处理环境，避免不同环境的库路径混入 MPI/HIP 进程。

## 8. 输出文件

| 文件 | 内容 |
|---|---|
| `step1.in`、`step2.in` | 本次输入 |
| `step1-STEP.dat`、`step2-STEP.dat` | 各时刻密度场 |
| `step1-STEP.png`、`step2-STEP.png` | 灰度图片 |
| `step1.out`、`step2.out` | 原六列科学记录 |
| `step1.stdout`、`step2.stdout` | 求解器输出及逐步计时 |
| `step1.stderr`、`step2.stderr` | 求解器错误输出 |
| `step2-10000.xy`、`step2-10000.nh` | 末态原子坐标和邻居信息 |
| `graphene.xyz`、`graphene.data` | 原子坐标中间文件和 LAMMPS 初始构型 |

密度文件保持无表头、每行一个数、x 方向优先的原格式。新名称为 `PREFIX-STEP.dat` 和 `PREFIX-STEP.png`。求解器拒绝覆盖已有输出；旧重启路径先按原名查找，缺失时将 `step1-t-10000.dat` 映射到 `step1-10000.dat`。

## 验证结果

固定 100×100 nm 性能算例使用 4240² 网格、100 个初始晶核、每阶段 2000 步。原 CPU 两阶段共 458.16 秒，单节点四卡 DCU 共 83.12 秒，加速比 5.51 倍；两侧末态均通过原容差对照，转换后均为 388638 个原子。

计时覆盖两阶段完整 PFC 求解，包括 MPI 启动、初始化、FFT 规划、推进、通信、日志和数据落盘，排除图片导出及 LAMMPS 转换。上述实测对应固定硬件和工作量，与本文 400 晶核、每阶段 10000 步的运行示例不同。算例、容差与适用范围见 [validation.md](specs/001-dcu-pfc-acceleration/validation.md)。

## 环境模块与 Slurm 部署

需要通过 `module load pfc/2026.1` 管理环境时，先按本站模块名调整 `slurm/modulefiles/pfc/2026.1`，再从仓库根目录注册：

```bash
cd "$PFC_ROOT"
module use "$PWD/slurm/modulefiles"
module load pfc/2026.1
```

`slurm/` 提供批作业示例，脚本中的 `Debug` 分区和资源声明来自特定集群的部署约定，部署到其他环境时按站点配置调整。模块与模板用法见 [Host.md](Host.md)。个人 SSH、账户、私有配置和运行记录集中保存在被 Git 忽略的 `yhmsr/` 中。
