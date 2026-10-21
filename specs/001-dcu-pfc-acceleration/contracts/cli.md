# CLI Contract

日期：2026-10-20。约束来源：[规格](../spec.md)；技术选择见[研究](../research.md)。
以下程序命令仅在目标集群的Slurm Debug批作业内执行，禁止其他分区或回退。
安装路径和模块名按站点配置；从仓库根目录注册项目模块并提交作业：

```bash
module use "$PWD/slurm/modulefiles"
module load pfc/2026.1
sbatch slurm/build.slm
# 等待构建作业成功完成后再提交计算。
sbatch slurm/polygr.slm
```

脚本头部使用 `#SBATCH --partition=Debug`；日常脚本不包含测试、探针、元数据生成或外层计时。

## 求解器

```text
bin/polygrpfc INPUT_PREFIX
bin/polygrpfc_cpu_ref INPUT_PREFIX
```

- `INPUT_PREFIX` 保持原义：读取原样的 `INPUT_PREFIX.in`，含目录时保持目录不变。
  不增加物理输入字段，不要求改写既有 `S/O/A/I/M/R` 标签、顺序或参数。
- `bin/polygrpfc` 是以 `PFC_DCU` 构建的 C 主程序及 HIP 后端，必须使用已分配 DCU。
  MPI/设备配置来自 Slurm 与启动环境，不能静默回退为 CPU。
- `bin/polygrpfc_cpu_ref` 从规格冻结提交的原 `src/pfc.c` 构建，保留原 CLI、科学行为、
  旧输出名称和诊断；不能为了本契约的计时或命名修改它。两种程序使用隔离运行目录。
- 新程序输出、安全前缀、重启回退与计时遵循[数据与日志契约](data-and-logs.md)。
  同一输入内多个 `R` 块仍按原顺序执行，不把输入限制为仅有一个 `R`。
- 缺少参数、无法读取输入、非法设备绑定、分配/通信/设备/I/O 失败均须以非零状态退出；
  MPI 作业不得在某 rank 已失败后继续报告成功。错误说明包含相关路径或资源及原因。

## Python 灰度导出

```text
python source/plotter.py INPUT OUTPUT W H [START INCREMENT END]
```

| 参数 | 契约 |
|---|---|
| `INPUT` | 单帧的显式 `.dat` 路径，或批量模式中恰含一个 `#` 的 `.dat` 路径模板 |
| `OUTPUT` | 显式 `.png` 路径；批量模式恰含一个 `#`，展开后采用规范的前缀与步数 |
| `W H` | 正整数，必须与原场宽高一致；不从数据中猜测表头或尺寸 |
| `START INCREMENT END` | 必须三个同时提供；非负整数起止、正间隔且 `START <= END` |

单帧示例：

```text
python source/plotter.py step2-128.dat step2-128.png 256 256
```

批量示例：

```text
python source/plotter.py 'step1-#.dat' 'step1-#.png' 256 256 0 128 128
python source/plotter.py 'step1-t-#.dat' 'step1-#.png' 256 256 0 128 128
```

- 批量按 `START + k*INCREMENT <= END` 的整数顺序展开；不依赖目录字典排序。
  单帧不接受 `#`。输出 basename 必须含非空前缀、连字符和非负整数 step，
  再按数据契约规范化；步数不补零。单帧输出 step 与已知输入帧的 step 一致。
- 输入先尝试精确路径，仅在已知旧格式路径缺失时尝试约定的新名映射。
  图片不能覆盖已有文件；同一批次的规范化目标若重名，必须失败。
- 批量缺帧告警后继续；单帧缺失、全部缺帧、非法区间、损坏数据或输出失败为非零退出。
  数据损坏不能作为“缺帧”跳过，也不能生成部分合法数据的图片。
- 成功汇总输出 `PFC_PLOT generated=<整数> missing=<整数>`；警告与错误包含路径和原因。
  成功至少生成一张图；不将已存在图片计为本次成功生成。
- 使用研究冻结的 Python/NumPy/Pillow 环境，每次处理一帧；灰度导出不依赖 Java。

## 既有原子转换

继续调用仓库的 `java/coordinator.jar`，保留原有参数次序：

```text
java -jar java/coordinator.jar step2-128.dat 256 256 0.7 0.7 7.3 2.46 \
  step2-128.xy step2-128.nh
```

仅将原脚本引用的密度、XY、邻居文件路径改为新名；后续 XYZ/OVITO/LAMMPS 链路保持原义。
不能用接口不同的 `source/Coordinator.java` 重编产物替换现有 JAR。
图片导出和原子转换均在完整求解耗时的计时区间之外，但仍在批作业内运行。
