# Data and Logs Contract

日期：2026-10-20。[规格](../spec.md)规定验收；[研究](../research.md)记录算法与资源选择。
本契约适用于新 DCU 程序及 Python 绘图；原 CPU 参考保持原样。

## 输入、场布局与命名

- 原输入字段、注释、解析顺序和科学含义保持兼容；原路径不做重命名。
- 密度文件保持无表头、每行一个 `%e` 实数，共 `W*H` 个样本；行内 x 优先、随后 y，
  对应 `q[y*W+x]`。MPI padding 不写出，不添加 rank 分片或尺寸行。
- 新场与图片分别为 `<safe-prefix>-<step>.dat`、`<safe-prefix>-<step>.png`；
  `step` 是不补零的非负十进制整数。科学诊断仍写 `<safe-prefix>.out`。
- 只规范化 basename 的前缀，目录字节保持不变。原文件名字节中 ASCII
  `A-Z a-z 0-9 _ -` 原样保留，其他每一个字节替换为一个 `_`；前缀不得为空。
  C 与 Python 使用相同规则，Python 按 `os.fsencode()` 得到文件名字节后处理，
  不能按 Unicode 字符逐个替换而与 C 产生不同名称。
- 每个调用先拒绝独立运行已有的目标文件及已知的规范化冲突，不生成旧名副本。
  创建时仍须排除并发覆盖；写入/关闭错误使作业失败，部分文件不能算完整结果。
- 同一求解器调用内多个 `R` 块可按原行为再次写入本调用自己创建的同一步文件；
  此例外不允许覆盖启动前已存在的结果，不改变原科学更新或保存调度。
  两个阶段或独立调用必须使用不同安全前缀或独立目录。

重启与绘图读取的解析顺序固定为：

1. 尝试输入指定的精确路径；存在时使用它，即使对应新名称也存在。
2. 仅当路径不存在，且 basename 确切符合 `<prefix>-t-<step>.dat` 时，
   将该前缀按上述字节规则规范化，并在同一目录尝试 `<safe-prefix>-<step>.dat`。
3. 回退必须记录原始路径和实际路径；其他形式不猜测、不扫描目录、不移动或复制旧文件。
   求解器仍执行原重启局部统计与重标度，不能把读取兼容改成新的归一化算法。

## 灰度像素

Python 必须验证单列实数、样本数与有限性；空数据、表头、复数、额外列、截断、
NaN/Inf 或非法值均失败。每帧使用该帧自身的有限最小值 `m`、最大值 `M`：

```text
M == m: gray = 0
M > m:  gray = floor(255 * (q - m)/(M - m) + 0.5), 范围 0..255
PNG(x, H-1-y) = gray(q[y*W+x])
```

输出尺寸严格为 `(W,H)`、8 位 `L` 灰度，不转置 x/y、不增加坐标轴或边框。
计算使用 FP64，最终量化为 8 位；相对原灰度约定最多允许 1 灰阶差。
具体单帧、批量和错误行为见 [CLI 契约](cli.md)。

## 科学诊断与逐步计时

原科学 `.out` 保持六列数值格式及原诊断频率：

```text
step elapsed_integer_seconds dx dy mean_free_energy mean_density
```

新记录仅由 rank 0 写作业标准输出，以 `PFC_TIMING` 区分，固定使用秒：

```text
PFC_TIMING kind=init phase=step1 seconds=0.123456789
PFC_TIMING kind=step phase=step1 block=1 step=1 seconds=0.012345678 cumulative_seconds=0.023456789
PFC_TIMING kind=optimize phase=step2 block=1 step=100 seconds=0.034567890
PFC_TIMING kind=total phase=step1 seconds=1.234567890
```

- `phase` 为安全 basename 前缀；`block` 为本输入中 `R` 块的 1 起始顺序。
  每块 T 次推进恰有 T 条 `kind=step`，step 为 `1..T`；初态 0 无推进耗时记录。
- `seconds` 为完成的该次推进墙钟时间，包含设备完成、FFT、转置和 MPI 等待。
  每个执行leader缓存本块各步时间，成员缓冲保持0，块结束后统一`MPI_MAX`，由rank0按步数输出。
  优化及累计时间同样只由执行leader记录，避免将成员等待此前GPU推进的时间归到下一次操作。
  不新增每步集合通信，不承诺运行中实时输出每一步。
- `cumulative_seconds` 为自本 `R` 块开始至该步完成的墙钟累计值，同样取最慢 rank，
  包含期间原有优化、诊断和写场；它不等于仅把各步 `seconds` 相加。
- `kind=optimize` 覆盖每次原优化调用，包括原调度发生的终止时刻优化；不计入推进条数。
  `kind=init` 报告初始化成本；`kind=total` 覆盖该求解器调用的内部总时间。
  时间必须有限、非负；累计值按 step 不递减。记录不能触发新增 `fp()` 或物理更新。
- `kind=total` 不替代外部启动计时。完整性能验收由 `/usr/bin/time` 分别包住两次 MPI
  启动，含启动、初始化、FFT 规划、全部计算/优化/通信、I/O、计时输出和退出，
  两阶段时间相加为 `T_CPU` 或 `T_4DCU`；排除编译、排队、绘图和原子后处理。

## 验收证据

执行进程的`PFC_PROFILE`补充累计转置次数、等待GPU完成及D2H、MPI、H2D和已有同步耗时，
单位为秒。计时仅包围既有操作，不添加同步或物理调用；各rank等待可重叠，不能将这些值
跨rank相加充当总耗时，也不能把ready_download解释为纯DMA时间。验收仍只使用外部完整PFC计时。

证据记录须关联原 CPU 提交、输入与初始场、种子、真实 MPI/slab 分解、CPU/DCU 资源、
模块/库/编译选项、可选调度作业身份、受管集群的实际 `SLURM_JOB_PARTITION`、退出状态、日志、误差和外部两阶段时间。
本会话代理的受管集群作业必须在脚本头部声明 `#SBATCH --partition=Debug`，不得使用分区回退。公共用户可按README直接编译运行，Slurm/Debug不是公共使用前提。
日常 `.slm` 不包含测试、探针、证据/元数据生成或外层计时；这些约束不改变求解器逐步计时。
只有规格冻结的 V1/V2/V3；两阶段末态按冻结容差通过后，才能报告 `T_CPU/T_4DCU >= 5`。
失败、未运行或资料不完整的项目保持未验证，不以日志格式正确替代数值或性能通过。
