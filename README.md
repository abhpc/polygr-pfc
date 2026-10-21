# DCU 加速的多晶石墨烯 PFC 程序

本项目用相场晶体方法（PFC）生成多晶石墨烯，再转换为 LAMMPS 可读取的原子初始构型。程序基于 [Petri Hirvonen 的 PFC 工具](https://github.com/petenez/pfc)，保留原物理方程、离散算法和双精度计算。

计算分为两个阶段：`step1` 生成并生长晶粒，`step2` 读取第一阶段末态进行结构弛豫。DCU 负责 FFT 和场更新；原子转换沿用 Java Coordinator 和 OVITO，灰度图片由 Python 导出。

编译、测试、计算及后处理均通过 Slurm 作业执行，分区固定为 `Debug`。登录节点用于传文件、编辑、加载环境和提交作业。公开文件不记录个人账户、SSH 地址或个人安装路径；本地集群配置与运行记录放在被 Git 忽略的 `yhmsr/` 中。

## 1. 准备集群环境

将仓库放到计算节点可访问的共享存储中，进入仓库根目录：

```bash
cd /path/to/polygr-pfc
module av
```

需要 MPI、GCC、DTK/HIP、hipFFT 和 MPI 版 FFTW。后处理还需要 Java、OVITO 2.9.0，以及安装了 NumPy 和 Pillow 的 Python。

仓库中的 modulefile 使用以下模块名作为部署示例：

```text
gcc/13.1.0
mpi/2021.9.0
dtk/26.04
```

先用 `module av` 查询本站环境，再修改 `slurm/modulefiles/pfc/2026.1` 中的模块名。该构建使用 Intel MPI；`I_MPI_ROOT`、`DTKROOT` 和 `HIP_PATH` 由工具链模块提供。`slurm/polygr.slm` 和 `slurm/convert.sh` 中的 Python、Java、OVITO 模块名也需要与本站一致。

随后注册并加载 PFC 模块：

```bash
module use "$PWD/slurm/modulefiles"
module load pfc/2026.1
mkdir -p "$PFC_ROOT/runs"
cd "$PFC_ROOT/runs"
```

模块从自身在仓库中的位置确定 `PFC_ROOT`，设置程序和 FFTW 库的搜索路径。无需写入个人绝对路径；模块加载本身不进行计算。新登录会话需重新执行 `module use` 和 `module load`。更多部署说明见 [Host.md](Host.md)。

## 2. 首次安装 FFTW

已有兼容的 MPI 版 FFTW 时，可以修改 modulefile 中的 `FFTW_PREFIX` 及库路径并跳过本步。默认安装位置为 `$PFC_ROOT/.local/fftw-3.3.10/`。

先准备官方源码包；下载和传文件可以在登录节点完成：

```bash
mkdir -p "$PFC_ROOT/.local/sources"
curl --fail --location https://www.fftw.org/fftw-3.3.10.tar.gz \
  -o "$PFC_ROOT/.local/sources/fftw-3.3.10.tar.gz"
sbatch "$PFC_ROOT/slurm/fftw.slm"
```

安装作业会解压、配置、编译并安装 FFTW。目标安装目录已存在时不覆盖；已有安装直接用于后续编译。

提交成功仅表示作业进入队列。用提交时返回的编号查询，确认状态为 `COMPLETED` 且 `ExitCode` 为 `0:0` 后再继续：

```bash
sacct -j JOB_ID --format=JobID,State,Elapsed,ExitCode
```

## 3. 编译 PFC

在集群提交：

```bash
sbatch "$PFC_ROOT/slurm/build.slm"
```

脚本执行 `make bin/polygrpfc`，生成 `$PFC_ROOT/bin/polygrpfc`。默认目标设备为 `gfx906`，其他 DCU 型号需修改 `build.slm` 中的 `PFC_GPU_ARCH`。

主机代码使用 `-O3 -march=native -ffp-contract=off -Wall`，设备使用 FP64，不启用 fast-math 或浮点收缩。FFT 链接本地 `hipfft` 和 MPI 版 FFTW；不要混用不同 MPI ABI 的库。确认编译作业完成后再提交计算。

## 4. 提交两阶段计算

```bash
sbatch "$PFC_ROOT/slurm/polygr.slm"
```

脚本在 `runs/JobID/dcu/` 中复制并填写输入文件，然后依次运行两个阶段、转换 LAMMPS data、导出 PNG。资源在 `#SBATCH` 中指定：`Debug` 分区、单节点、32 个 MPI rank、4 张 DCU。账户和 QoS 使用集群默认值；需要显式配置时，在自己的作业副本中填写。

两阶段的实际求解命令为：

```bash
mpirun -bootstrap fork -np 32 -ppn 32 polygrpfc step1 > step1.stdout 2> step1.stderr
mpirun -bootstrap fork -np 32 -ppn 32 polygrpfc step2 > step2.stdout 2> step2.stderr
```

这些命令只能放在 Slurm 作业内执行。程序读取当前目录的 `step1.in` 和 `step2.in`，第二阶段按输入中的路径读取第一阶段末态。脚本按顺序执行，命令出错即停止，不包含额外检测、测试或外层计时命令；程序自身保留逐步计算时间日志。

当前验证使用总计 32 rank 的划分，每卡由一个执行进程聚合成员的 FFT 和场更新。初始化和局部统计依赖 rank 划分，不要仅为匹配卡数把进程数改成 4，也不要添加每 rank 独占一张卡的绑定。

跨节点时，已验证的配置为两节点、每节点 16 rank 和 4 张卡，总计仍为 32 rank。相应修改 `#SBATCH --nodes=2`、`#SBATCH --ntasks-per-node=16`，以及两条求解命令的 `-bootstrap slurm -np 32 -ppn 16`。`Debug` 分区保持不变。

## 5. 设置尺寸、晶核数和步数

修改 `slurm/polygr.slm` 中的参数：

| 参数 | 默认值 | 含义 |
|---|---:|---|
| `W`、`H` | 1000、1000 | 薄片尺寸，单位 Å |
| `N` | 100 | 初始晶核数 |
| `R` | 2 | 晶核半径，单位为晶格常数 |

工作流沿用旧换算 `int(1024/243*W)`，默认网格为 4213²。要使用已验证的约 100×100 nm 网格，直接将 `PFC_WIDTH` 和 `PFC_HEIGHT` 设为 `4240`；原转换尺度下名义尺寸为 100.017534 nm。模拟 400 个初始晶核时，将 `N` 设为 `400`；生长后保留下来的晶粒数可能变化。

物理参数和推进步数在 `example/step1.in`、`example/step2.in` 中：

| 标签 | 含义 |
|---|---|
| `S` | 随机种子；默认工作流随机生成第一阶段种子 |
| `O` | 科学日志间隔、密度保存间隔 |
| `A` | 网格宽、高 |
| `I` | 初始化方式；第二阶段用 `I 2` 读取第一阶段结果 |
| `M` | 原物理模型参数 |
| `R` | 步数、网格间距、时间步长、盒优化间隔 |

脚本中的晶核半径 `R` 与输入文件的 `R` 标签含义不同。输入文件注释给出了各参数顺序。

默认每阶段 10000 步，每 1000 步保存一次。修改步数时，同步修改两份输入的 `R` 行、第二阶段重启文件名，以及 `polygr.slm` 中的末态转换文件和绘图区间。固定随机种子后才能复现同一初始构型。

## 6. 查看输出

在登录节点查询作业：

```bash
squeue --me
sacct -j JOB_ID --format=JobID,State,Elapsed,ExitCode
```

Slurm 总日志为提交目录中的 `JobID.out` 和 `JobID.err`，计算结果在 `$PFC_ROOT/runs/JobID/dcu/`：

| 文件 | 内容 |
|---|---|
| `step1.in`、`step2.in` | 本次输入 |
| `step1-STEP.dat`、`step2-STEP.dat` | 各时刻密度场 |
| `step1-STEP.png`、`step2-STEP.png` | 灰度图片 |
| `step1.out`、`step2.out` | 原六列科学记录 |
| `step1.stdout`、`step2.stdout` | 求解器输出及逐步计时 |
| `step1.stderr`、`step2.stderr` | 求解器错误输出 |
| `step2-10000.xy`、`step2-10000.nh` | 默认末态的原子坐标和邻居信息 |
| `graphene.data` | LAMMPS 原子初始构型 |

新文件名使用 `PREFIX-STEP.dat` 和 `PREFIX-STEP.png`。密度文件仍无表头、每行一个数、x 方向优先，已有输出不会被覆盖。旧重启路径先按原名查找，缺失时将 `step1-t-10000.dat` 映射到 `step1-10000.dat`。

灰度图片逐帧归一化为 8 位并上下翻转，常量场输出黑图。独立批量导出的接口如下，也需放入 Slurm 脚本：

```text
python source/plotter.py 'step1-#.dat' 'step1-#.png' W H START INCREMENT END
```

原子转换保持 `0.7 0.7 7.3 2.46` 尺度，沿用 `java/coordinator.jar`。不要用接口不同的新版 `Coordinator.java` 直接重编替换。

## 验证结果

固定 100×100 nm 算例（4240²，每阶段 2000 步）已通过原容差 CPU/DCU 数值对照：原 CPU 两阶段共 458.16 秒，单节点四卡 DCU 共 83.12 秒，加速比 5.51 倍。计时包含两阶段完整 PFC 求解，排除 LAMMPS 转换和图片导出；两侧转换后均为 388638 个原子。

该结果对应已验证的硬件和固定算例，不代表任意尺寸与环境都能达到相同加速比。算例、容差与验证范围见 [validation.md](specs/001-dcu-pfc-acceleration/validation.md)。日常运行脚本不重复执行验收。
