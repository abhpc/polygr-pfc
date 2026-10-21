# Research: DCU 双精度 PFC 加速

本文记录实现采用的技术选择与科学不变量；冻结输入见第10节，已验证结果见 [validation.md](validation.md)。

## 1. 依赖选择

采用HIP/本地hipFFT、Intel MPI以及FFTW-MPI双精度库。
已验证组合为 `dtk/26.04`、`gcc/13.1.0`、`mpi/2021.9.0` 和FFTW 3.3.10。
HIP头文件版本为6.3.26113，hipFFT为1.0.14，rocFFT为1.0.28。
安装前缀分别由 `DTKROOT`、`HIP_PATH`、`I_MPI_ROOT` 和 `FFTW_PREFIX` 指定，不固定站点绝对路径。

FFTW从[官方3.3.10归档](https://www.fftw.org/fftw-3.3.10.tar.gz)在Slurm Debug批作业内构建，
使用 `--enable-mpi --enable-shared --enable-sse2 --enable-avx --enable-avx2` 和第6节共同CFLAGS。
CPU参考与DCU程序使用同一FFTW安装；后者仅用其准确的科学slab划分查询。
目标CPU须支持所选指令集，不能以标量库或不明预编译依赖削弱CPU基准。

所有编译、测试与计算仅在目标Slurm集群的Debug批作业内执行，不使用其他分区或回退。
日常 `.slm` 仅包含顺序构建/运行命令，分区与资源写在 `#SBATCH` 中；
不自动执行测试、探针、证据生成或外层计时。账户、QoS、安装位置与资源上限按目标站点配置。

## 2. 不使用集群 hipFFT-Mp

**Decision**：直接链接 `-lhipfft`，不链接 `libhipfftMp.so`，不使用把它自动加入的
`hipfft_LIBRARIES` CMake 汇总变量。

**Evidence**：集群 `libhipfftMp.so` 虽存在分布式接口，却引用 `ompi_mpi_comm_world`、
`ompi_mpi_double` 等 OpenMPI 对象；本地 `libhipfft.so` 无此 MPI 符号依赖，且导出
`hipfftPlanMany`、`hipfftExecD2Z/Z2D/Z2Z`。厂商版本可能有回移，不能仅凭版本号判断能力。
进一步静态核对 `hipfft → rocfft → galaxyhip` 的未定义符号，均未发现 MPI/OpenMPI
对象依赖；动态加载和设备兼容性仍留作批作业准入，不能由 ELF 阅读断言运行成功。

**Rationale**：避免将 OpenMPI 句柄传给 Intel MPI。采用本地 FFT 加显式 MPI 转置，
既支持跨节点，也不假定当前 MPI 能收发设备指针。

**Alternatives considered**：切换整个项目到 OpenMPI，或另装支持当前 MPI 的分布式 FFT，
均比使用已安装的本地变换增加依赖与验证范围。

参考：[hipFFT API](https://rocm.docs.amd.com/projects/hipFFT/en/latest/reference/hipfft-api-usage.html)、
[rocFFT 分布式变换约束](https://rocm.docs.amd.com/projects/rocFFT/en/latest/how-to/distributed-transforms.html)。
上游文档只解释接口；本计划的 ABI 判断来自集群实际 ELF/头文件读取。

## 3. 原 MPI 科学划分保持不变

**Decision**：CPU/DCU 验收均为 32 个真实 MPI rank。单节点 DCU 每卡 8 rank；
跨节点保持总 rank 数 32，每节点 16 rank、4 DCU，每卡 4 rank。
继续调用原 `fftw_mpi_local_size_2d_transposed(H, W/2+1, ...)` 确定 slab。

**Rationale**：`seed_rngs()` 的随机序列与 `read_state()` 的局部均值/极差重标度依赖 rank。
原 CPU 是纯 MPI，设置 `OMP_NUM_THREADS=8` 不能让 4 个 rank 自动使用 32 个核心。
相同真实 rank 数及 slab 保持科学含义，CPU 同时占满 32 个物理核心。

**Alternatives considered**：4 rank CPU 是被削弱的基准；把重启改为全局归一化改变算法；
不引入虚拟科学rank；同卡真实rank的FFT批次交由每卡最低rank执行，
其余真实rank仍负责原初始化、局部统计与fp累加。该通信/执行聚合不改变科学分解。
设备共享和主机转置是明确性能风险，不能据设计宣称达到 5 倍。

参考：[FFTW 划分查询](https://www.fftw.org/fftw3_doc/MPI-Data-Distribution-Functions.html)。
该接口返回实际起点、长度及分配容量，不能用简单整除替代，也不能漏掉空 slab。

## 4. 分布式 FFT 与常驻设备数据

**Decision**：科学实场仍按原y slab，科学频场仍按原kx slab；同卡成员的设备数据由最低rank
执行进程聚合，明确保存各原slab的global x/y映射。每卡一个执行进程、本地FP64批量变换与
GPU打包/解包之间，同节点使用注册共享主机缓冲，跨节点使用页锁定缓冲及 `MPI_Alltoallv`。

正变换：本地W方向 `D2Z` → 按目标kx slab打包 → 共享主机区或Alltoallv → H2D →
排为 `[local_kx][global_ky]` → H 方向 `Z2Z` 正变换；逆变换反向执行。
保留原 `2*(W/2+1)` 实数行跨度，物理 W 点之外的 padding 不参与计算或输出。
变换不自动归一化，只保留原 `scale_P/scale_Q/B` 位置的 `1/(W*H)`。

**Rationale**：不要求 GPU-aware MPI；不整除、奇数尺寸及空 slab 也使用实际元数据。
主要密度/频场、A/B 与 FFT 工作区常驻设备，逐步不下载完整场；必要转置暂存和定期输出除外。
周期刷新在GPU执行，次数与原调度不变；同卡区块直接D2D复制。
MPI 以双精度分量计数，检查乘法和 `int` 计数溢出，不静默截断。

**Alternatives considered**：每卡复制全局场不是真正跨节点并行；只加速逐点循环而保留 CPU
FFT 不能合理承载本目标；GPU-aware MPI 快路径没有可用性证据，首版不增加可选路径。

参考：[FFTW 实数据 MPI 布局](https://www.fftw.org/fftw3_doc/Multi_002ddimensional-MPI-DFTs-of-Real-Data.html)。

## 5. 最小 C/HIP 衔接与物理不变量

**Decision**：保留 `src/pfc.c` 的 C 主流程，增加窄 C 接口和一个 HIP C++ 翻译单元。
`PFC_DCU` 构建仅替换数组配置、步推进和诊断FFT的执行路径；A/B保留原C计算并上传，
随机初始化、重启局部统计、主机能量累加、盒优化公式和调度尽量保留原代码。
原 CPU 基准从固定提交提取原始 `src/pfc.c`，不以修改后的源码冒充原程序。

| 原代码 | 迁移必须保持 |
|---|---|
| `seed_rngs/polycrystalline_state/read_state` | libc 随机序列、严格最近晶核比较、局部重启统计 |
| `update_AB` | 原 kx/ky 索引规则；`L==0` 分支；B 中已有的归一化 |
| `step` | 非线性项、谱更新及 `t%100==0` 的推进前重建 |
| `fp` | 恢复实场，同时重新建立归一化持久谱场 P；不是无副作用查询 |
| `optimize` | 五次有序采样、二次插值、盒位移截断和 d 调整 |
| `relax` | 优化 → 诊断 → 保存 → 推进；终止时刻仍可执行优化 |

GPU 完成诊断 FFT/谱算子后，下载原密度与算子结果，由原主机逐行公式累加自由能与密度，
保持 32 rank 的原 MPI 归约。暂不改为 GPU 树形归约，以免给盒优化五个相近能量的差分
额外引入舍入来源。代价只发生在原诊断/优化时刻，并计入完整时间。

**Alternatives considered**：整文件改成 C++ 会牵涉 C99 复数；复制一套完整求解器容易使
物理实现分叉；通用后端框架和可重现求和库均非必要。

## 6. 浮点策略与数值门槛

**Decision**：原 CPU 与 DCU 主机共同采用 GCC 13.1、
`-O3 -march=native -ffp-contract=off -Wall`；设备采用 `-O3 -ffp-contract=off`、FP64，
不启用 fast-math、降精度或额外自动缩放。构建时针对已分配设备记录并选择真实 gfx ISA。

**Rationale**：O3 与原生指令集提供优化，关闭收缩/快速重结合有利于保留科学累加及非有限值检查。
这与历史 README 的 `-Ofast` 配方不同；本次结论只针对冻结的共同浮点策略，不宣称是所有
编译选项中的最快 CPU，也不进行选项扫描。基准仍为未改写源码、32 核完整优化的纯 CPU。

冻结规格上限：密度 `atol=1e-6, rtol=1e-5`；能量/均密度/盒尺寸
`atol=1e-8, rtol=1e-6`。原 `%e` 序列化约 7 位有效数字，正常量级样本存在约末位的量化；
该事实不是上述容差的推导或 FP64 证明。盒优化差分可能放大舍入，故保留原主机累加，
把 `dx/dy` 和非有限值纳入硬门槛；若条件敏感而失败，应修正执行路径或报告未通过，
不得改物理公式、关闭优化或事后放宽容差。采用上限，实施前不得根据运行结果反向挑选阈值。

## 7. Python 灰度导出与已有转换链

**Decision**：`source/plotter.py` 使用已安装的 `anaconda3/2023.07`：包元数据表明
Python 3.11.3、NumPy 1.24.3、Pillow 9.4.0。仅在后处理子 Shell 中加载，
避免将 Anaconda 库路径混入 HIP/MPI 求解环境；不新增 Matplotlib。

读取无表头单列 FP64 文本并验证长度/有限性；`reshape(H,W)` 后上下翻转。
逐帧 min/max 线性映射并四舍五入到 0..255，常量场全黑，输出 PNG 的 8 位 L 模式。
保留旧批处理的输入/输出模板、显式 W/H 和起点/间隔/终点模式，但明确写出文件扩展名。

**Rationale**：库已安装，适合真实数据格式；不需要更重的绘图库。
静态检查表明 bundled JAR 与 `source/` 新版 Java 的表头/参数接口不一致，
继续调用现有 `java/coordinator.jar`，只更新路径，不重编替换它或迁移其他 Java 工具。

**Alternatives considered**：按新版 Java 源码要求加尺寸表头会破坏实际求解器契约；
2026.07 的 Python 3.14 环境不是本任务所需升级。

参考：[NumPy 1.24 loadtxt](https://numpy.org/doc/1.24/reference/generated/numpy.loadtxt.html)、
[Pillow 图像模式与坐标](https://pillow.readthedocs.io/en/stable/handbook/concepts.html)。

## 8. 文件命名与历史输入

**Decision**：新密度/PNG 为安全前缀加整数 step；C/Python 按相同文件名字节规则将
ASCII 白名单之外的字节替换为 `_`。输入目录与原路径不改写。
先读取精确旧路径；仅对缺失的已知 `-t-<step>.dat` 做到新名称的同目录回退。
不写旧名副本，不静默覆盖已有独立运行数据；输出与基准分别使用独立目录。

**Rationale**：兼容未修改的第二阶段引用，且避免 C 按字节、Python 按 Unicode 字符
产生不同输出前缀。具体文件及冲突规则见 contracts/。

**Alternatives considered**：批量重命名旧数据或强迫用户改输入不满足兼容；保留两套新旧文件
会引入歧义和重复 I/O。

## 9. 逐步时间与完整性能

**Decision**：每个实际 `t→t+1` 推进记为 `step=t+1`，设备完成后关闭该 rank 的墙钟区间。
缓存每步时间，在阶段结束统一做 MPI_MAX 并由 rank 0 输出，避免新增每步集合通信。
盒优化单独记，包括终止时刻的优化。原六列科学 `.out` 保留，新增记录使用独立
`PFC_TIMING` 前缀写入作业标准输出。

外部 `/usr/bin/time` 分别包住两个阶段的完整 MPI 启动，求和得到 T；
原 CPU 源码不加计时代码。新增日志成本、FFT 规划和 I/O 都计入加速门槛。

**Alternatives considered**：只计异步 launch、只计 rank 0、用内部晚启动计时替代完整耗时，
或每步调用 `fp()`，均不符合已确认的验收范围/物理调度。

## 10. 冻结的最小验证与作业边界

**Decision**：只使用下列组，不扫描网格/参数/编译选项，也不重复运行生产时长。

| 组 | 网格、晶核 | 每阶段步数 | O: 诊断/写场 | 资源 |
|---|---|---|---|---|
| V1 | `256×256, N=16,R=2` | `128/128` | `64/128` | 原 CPU 与 DCU 各一次，单节点 32 rank，DCU 4 卡 |
| V2 | 复用 V1 | `128/128` | `64/128` | 2 节点，总 32 rank，每节点 4 卡，复用 V1 CPU |
| V3 | `4240×4240,N=100,R=2` | `2000/2000` | `1000/2000` | 同一独占单节点依次运行原 CPU 与 4 DCU，32 rank |

共同种子 `34789847`；初始化 `dx=dy=.7,l0=7.2552,p=.25,A=-.1`。
第一阶段模型 `(-.18,1,0,1,1)`、dt=.2、T_optimize=0；第二阶段读各自第一阶段末态，
p=-.072、A=.1、模型 `(-.15,1,.87,1,0)`、dt=1、T_optimize=100。
其他参数和局部重标度保持原模板；CPU/DCU 输入字节相同，第二阶段均保留旧 `-t-` 引用，
DCU 路径由兼容规则解析。V1 经过内部 t=100 的刷新与一次盒优化；V3 包含 t=2000 的终态优化。

**Rationale**：保留转换尺度dx=dy=.7、lPFC=7.3、lA=2.46 Å；每格为.7*2.46/7.3 Å，
N=ceil(1000*7.3/(.7*2.46))=4240，对应名义100.017534×100.017534 nm。
这是不小于100nm目标尺寸的最小整数网格，不为FFT友好性额外补点。
该名义周期长度不等于OVITO根据有限原子坐标自动推断的包围盒，原转换链不改。
128步覆盖短程分支，2000步保持冻结终点；种子、模型、N/R、步数和输出间隔均固定。
V1 较多晶核降低局部 slab 量化成常量的风险；若仍退化必须失败，不改归一化掩盖。

所有构建、比较、灰度边界检查、原子转换和性能解析均在Debug分区的Slurm作业，
脚本头部显式声明 `#SBATCH --partition=Debug`，不得回退。
资源申请按目标站点配置，不能以调度器的内存声明推断实际物理内存或显存。
CPU/DCU 性能配对在同一个节点顺序执行，直接消除 CPU 型号不一致；不计构建与后处理时间。

## 11. 结论与边界

实现保留32个真实科学rank，每卡一个执行进程聚合同卡FFT批次；A/B由原C循环和主机exp计算后上传，
保留L==0分支和原相消式，不替换为expm1。
已验证的冻结工作量达到完整5倍目标，详见[匿名结果](validation.md)。
硬件、FFT后端及MPI环境变化后仍需保持相同科学输入与容差，不得据既有结果推断新平台通过。
不增加全量回归、参数扫描、额外物理演化或日常脚本中的调试验证步骤。
