# Implementation Plan: DCU 双精度 PFC 加速与批量灰度导出

**Branch**: `master`
**Date**: 2026-10-20 | **Spec**: [spec.md](spec.md)

**Input**: `specs/001-dcu-pfc-acceleration/spec.md`。

## Summary

以现有 C/MPI 求解器为主流程，通过窄 C 接口调用 HIP FP64 内核及本地 hipFFT；
用页锁定主机暂存和MPI转置连接设备与节点，不依赖有ABI冲突的hipFFT-Mp。
保持原物理公式、真实 MPI slab、主机初始化/能量累加和盒优化；CPU/DCU 对照都使用 32 rank。
单节点 4 卡每卡 8 rank，跨节点验证为 2 节点、每节点 16 rank 和 4 卡。
新增 Python 灰度导出，统一安全文件名并兼容旧重启引用。

验收按完整两阶段求解耗时 `T_CPU/T_4DCU >= 5`，包含初始化、FFT 规划、日志及 I/O。
只运行固定 V1/V2/V3，验证冻结终点，不追加生产规模作业。规划完成不代表数值或加速已实测。

## Technical Context

**Language/Version**: 现有 C（C11）主流程；HIP C++14 翻译单元；Python 3.11.3 灰度/验收工具。

**Primary Dependencies**: `gcc/13.1.0`、`mpi/2021.9.0`、`dtk/26.04` 的 HIP 6.3.26113 /
本地 hipFFT 1.0.14（头文件值）、FFTW 3.3.10 MPI；后处理采用 `anaconda3/2023.07`
的 NumPy 1.24.3 / Pillow 9.4.0；保留 bundled `java/coordinator.jar`、OpenJDK 13.0.1 和 OVITO 2.9.0。

**Storage**: 集群共享文件系统上的输入、文本密度、PNG、坐标、LAMMPS 数据及作业日志。
无数据库。源码、依赖安装、二进制和每次运行目录分开。

**Testing**: 最小 assert/数值对照工具，不引入测试框架；所有构建、语法检查、比较与图片验证
仅在经 目标 Slurm 集群 提交的 Slurm 批作业。冻结输入和资源见 [research.md](research.md) 第 10 节。

**Target Platform**: 目标Slurm集群的 `Debug` 分区，验证资源为每节点32物理CPU核、4 DCU。
账户、模块安装位置和内存申请按目标站点配置，不使用登录节点硬件代替计算节点信息。
日常脚本以 `#SBATCH --partition=Debug` 明确分区，只包含顺序构建/运行命令。

**Project Type**: 科学计算命令行工具及 Slurm 批处理流程。

**Performance Goals**: 单节点 32 核原 CPU 与相同节点 4 DCU 比较，固定 V3 的完整两阶段耗时
至少 5 倍；跨节点要求正确完成，不附加未经要求的线性扩展门槛。

**Constraints**: FP64、不改物理/离散算法、不降低 CPU 利用、保持输入格式和原诊断副作用；
新名称为安全前缀加整数 step；禁止本机/登录节点测试、全量回归和参数扫描。

**Scale/Scope**: V1 为 256²、每阶段 128 步；V2 复用 V1，32 rank 跨两节点；
V3 为100×100 nm：4240²、名义100.017534 nm、每阶段2000步。灰度迁移仅 Plotter，不重写其他 Java 分析/原子转换工具。

## Constitution Check

*Gate 在 Phase 0 研究前和 Phase 1 设计后各审阅一次；下表是设计符合性，不是运行验收结果。*

| 关卡 | 研究前 | 设计后 | 依据 |
|---|---|---|---|
| DCU 范围 | 通过 | 通过 | C/HIP、构建/调度和验证直接服务于加速；Python绘图和新命名属于功能范围 |
| 物理不变 | 通过 | 通过 | 保留原算子/归一化、32 rank 初始化、主机 fp 累加和盒优化，详见研究第 3–6 节 |
| 兼容 | 通过 | 通过 | 原标签/字段、旧重启文件精确优先/同目录回退、旧 JAR 与转换参数保持 |
| 集群执行 | 通过 | 通过 | 所有编译与计算仅在Slurm Debug批作业执行，无其他分区或回退 |
| 证据与最小验证 | 通过 | 通过 | 基准提交、固定算例、容差上限、外部计时与作业证据有明确契约 |

不存在需要豁免的关卡。
真实库加载、设备计算、数值与 5 倍指标均属于后续运行关卡，未通过不得声称完成。

## Project Structure

### Documentation (this feature)

```text
specs/001-dcu-pfc-acceleration/
├── spec.md
├── plan.md
├── research.md
├── data-model.md
├── quickstart.md
├── contracts/
│   ├── cli.md
│   └── data-and-logs.md
└── checklists/requirements.md
```

实施任务见 `tasks.md`，匿名验收结果见 `validation.md`。

### Source Code (repository root)

实现文件位置如下，是否通过验收以任务清单及实际集群证据为准：

```text
src/
├── pfc.c                 # 共享 C 主流程，添加 PFC_DCU 分支与必要兼容/计时
├── pfc_dcu.h             # 必要的窄 C/HIP 接口
└── pfc_dcu.hip.cpp       # HIP 内核、本地 FFT、打包/转置支持
source/plotter.py         # Python 灰度工具
Makefile                  # C 对象、HIP 对象及明确的本地 FFT 链接
slurm/
├── fftw.slm              # 依赖构建
├── build.slm             # 程序构建
└── polygr.slm            # 适配原工作流、新命名及 Python 绘图
tests/
├── minimal.py            # 最小比较/断言，只有批作业调用
└── fixtures/
    ├── v1/step1.in, step2.in
    └── v3/step1.in, step2.in
```

**Structure Decision**: 沿用 C、source/、slurm/ 组织。只增加一个 HIP 翻译单元和必要头文件，
不增加通用后端框架。`bin/`、`.local/`、运行数据和日志作为生成物排除出源码提交。
CPU 参考从固定 git blob 导出到构建工作区，不从已加 DCU 分支的工作树推定“原 CPU”。

## Implementation Design

### C 主流程与 HIP 边界

- C 保留输入解析、libc 播种/初始化、原局部重启统计、科学诊断输出与盒优化。
  用 `PFC_DCU` 构建选择执行路径；设备失败明确退出，不静默使用 CPU 继续。
- HIP 接口限定为配置/销毁、主机初态上传、原C计算的A/B系数上传、一步推进、诊断变换与场下载。
  GPU 状态由 C 持有的单个句柄管理，参数/尺寸作为值传入，不移动整个 C99 复数程序到 C++。
- 初态上传一次；推进时主要实场、谱场和 A/B 常驻设备。输出/原诊断需要时才同步到主机。
  多个原 R 指令块沿用已有状态，不得在第二个块用过期主机数组覆盖 GPU 场。
- 原 `update_AB()` 的C表达式、主机exp及FP64顺序保持，按原频slab计算后上传，
  避免主机/设备指数函数差异在B的相消式中引入额外舍入。原 `fp()` 的主机累加与 MPI 归约保持；GPU 提供算子结果，并保留恢复 q、重建 P 的副作用。
  原 `optimize()` 五次采样与标量算术保留，包括终止时刻优化及 d 更新。

### FFT、归一化与通信

按原 FFTW 查询得到所有 rank 的 `lH/lh0/lW/lw0`，不平均重分局部初始化区域。
在保留32个真实科学rank的基础上合并同卡执行，降低同卡多进程FFT开销。
每卡由最低本地rank作为执行进程，聚合同卡成员的原实/频slab；其余成员继续执行原初始化、
局部重启统计、主机fp累加及全世界归约。局部集合只在原上传/诊断/输出/A-B更新位置发生。
推进时只由每卡执行进程计算聚合批次，其他成员在下一次原诊断/输出位置等待，不增加逐步计时集合。
本地D2Z/Z2D沿x，Z2Z沿y；转置用GPU pack/unpack。
单节点执行进程用MPI共享主机接收区，各GPU直接D2H推送其他卡的目标片段，目标GPU从本地接收区H2D；
本卡块直接在同一设备内复制，省去主机往返及Alltoallv的CPU复制。跨节点继续使用页锁定主机缓冲及执行进程间Alltoallv。
消息块按64字节对齐，原消息计数不变，pack/unpack映射跳过填充；MPI窗口、注册、容量及同步必须核查。
不依赖HIP IPC；PFC_TRANSFER明确记录host-shared或host-mpi实际模式。
源执行组i到目标组j的量为 `2*sum(lH_i)*sum(lW_j)` 个MPI_DOUBLE；用明确global x/y映射处理
组内不连续的原科学slab，逆向用反向布局，不将局部重标度改成按执行组统计。
检查尺寸/容量与 MPI 计数溢出；零批次跳过本地 FFT 计划但仍参加集合通信。

```mermaid
flowchart LR
  R[原 y slab 实场] --> X[本地 x D2Z]
  X --> S[GPU 打包]
  S --> M[单节点共享主机区或跨节点Alltoallv]
  M --> U[按 kx 解包]
  U --> Y[本地 y Z2Z]
  Y --> K[原转置 kx slab 谱场]
```

逆变换反向执行。所有 FFT 均不自动缩放，B、scale_P/Q 中的原归一化保持位置与次数。
保留原奇/偶 ky 索引规则、L==0 分支、周期边界和 t%100 重建，不能趁迁移改公式。
周期刷新保持在GPU执行。
仅链接本地 `-lhipfft`，明确排除 OpenMPI ABI 的 hipFFT-Mp。

### 设备与主机资源

用 MPI 共享内存通信域确定 node-local rank，在 Slurm 可见的 4 个设备中按 local_rank%4 选择。
不设置每 task 独占一张 GPU，也不取消调度器可见性限制；启动记录 global/local rank、主机、
设备身份与绑定，并记录leader/member角色及实际执行进程。每节点4个执行进程驱动4卡，
仍有32个真实科学rank；不得将host科学划分缩成4rank或全局重启统计。
CPU 基准 32 rank、每核一个；单节点 DCU 32 rank 同样绑定 32 核；V2 每节点 16 rank。
DCU进程读取GPU sysfs local_cpulist，在各节点原MPI分配CPU集合的并集中选唯一核心，
优先对应GPU的NUMA核心，不扩大分配范围；设备日志记录实际cpu，验证器核查无同节点重叠。
每节点先由local rank0读取设备PCI信息并广播，其他GPU执行进程在首次HIP调用前绑定，
避免已创建运行线程留在其他科学rank的核心；实际选卡后再次核对PCI身份。
主机诊断保持全局 rank 的原顺序/分块，V2 仅改变进程所在节点。
单节点在已有Debug分配内使用Intel MPI fork启动，同作业CPU/DCU相同；跨节点继续Slurm启动。
参考[Intel Hydra文档](https://www.intel.com/content/www/us/en/docs/mpi-library/developer-reference-linux/2021-15/hydra-environment-variables.html)。
构建和运行脚本必须在 `#SBATCH` 中显式指定Debug，不切换分区。
日常脚本不加入测试、元数据生成、外层计时或探针。

### 输入、文件与绘图

原 `polygrpfc INPUT_PREFIX` 仍读取 `INPUT_PREFIX.in`；资源不加入科学输入。
新名称、路径映射及碰撞处理见 [CLI 契约](contracts/cli.md) 和
[数据契约](contracts/data-and-logs.md)。保留原无表头 `%e` 文本，不增加尺寸表头。

Python 使用现有 NumPy/Pillow，按 `reshape(H,W)` 加上下翻转、逐帧线性灰度、8 位 PNG。
继续调用 bundled coordinator JAR 和 OVITO；保留既有 `0.7,0.7,7.3,2.46` 转换尺度参数，
只更新文件路径，不把发现的原流程尺度问题作为本次物理修复。
不同工具的模块环境在后处理子 Shell 中隔离，避免 Anaconda 改变求解器动态库搜索路径。

### 计时与数值验收

每次真实推进结束同步设备并记录执行leader墙钟耗时；成员计时缓冲保持0，每个R块结束一次MPI_MAX数组归约，
rank 0 输出逐步记录。无需实时流式日志；原 `.out` 六列不变。
盒优化、初始化和内部总时长分别记录，终态无推进时不得添加虚构 step。
外部计时分别覆盖两次完整 MPI 启动并求和，避免内部时间遗漏进程启动/FFT 规划/退出成本。

按规格冻结容差上限逐元素检查两阶段末态，同时对比原均密度、自由能和盒尺寸；
保持原 host 累加及共同浮点策略，减少盒拟合差分中的额外舍入来源。阈值不是 FP64 证明；
若敏感性导致超限，必须失败，不能改优化策略或事后放宽。
CPU/DCU 主机均用 `-O3 -march=native -ffp-contract=off -Wall`，设备 `-O3 -ffp-contract=off`；
CPU 基准不加计时代码。具体数据精度和计时字段见契约。
FFTW 依赖自身也固定为同一原生 O3 构建及 MPI/SSE2/AVX/AVX2 配置，先在构建节点检查支持，
保留构建证据；不接受来源不明或标量化的 FFTW 库作为 CPU 性能基准。

## Frozen Validation and Delivery

研究第 10 节冻结全部科学参数；测试输入作为只读 fixtures，CPU/DCU 两份内容相同。

| 验证 | 必须证明 | 不额外执行 |
|---|---|---|
| 构建准入 | 实际 CPU/设备/ISA、Intel MPI 与本地 FFT 链接、FP64 构建、Python/旧 JAR 可用 | 不运行上游全量套件，不在登录节点做任何可执行检查 |
| V1 | 256² 两阶段128步、旧重启引用、数值/最终转换、命名、每步日志、至多3帧及最小坏数据检查 | 不新增模拟来测试绘图 |
| V2 | 同一256²算例，32rank跨2节点，节点设备均实际工作，复用V1 CPU末态 | 不做节点数/扩展率扫描 |
| V3 | 4240²两阶段2000步，同节点32核CPU/4DCU各一次，末态合格且总时间比≥5 | 不先额外跑一次代表性算例，不追加生产时长 |

构建准入与实际检查只能在相应 Slurm 作业中；预先核对 fixture 哈希、输出目录和资源，
不得在观察性能后改变网格、步数、物理参数、输出间隔或 CPU 配额。
运行目录采用 `PFC_RUN_ROOT/SLURM_JOB_ID`，参考与候选输出独立。
失败只允许有记录的定向排错/重跑，不自动扩大测试矩阵；未达标保留未完成状态。

实施顺序：先完成 C/HIP 最小接口及构建，再完成准确 FFT/原物理推进、命名兼容和计时；
随后接入Python/原转换链，并完成冻结V1/V2/V3的独立验收；日常脚本不自动执行验收。
实现检查项见 `tasks.md`，实际结果见 `validation.md`。

## Complexity Tracking

无宪章违规。窄 C 接口、本地 FFT 的 MPI 转置和最小验收工具均为已要求功能所必需。
不增加多后端框架、虚拟科学 rank、GPU-aware MPI 可选路径、GPU 树形能量归约或绘图库选型开关。
32 rank 共享设备、暂存转置、诊断传输和原串行文本 I/O 可能限制完整5倍目标；
只有冻结的集群配对能判定达标，不能用设计审阅替代运行证据。
