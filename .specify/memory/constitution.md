<!--
Sync Impact Report
- Version change: 1.2.1 → 1.2.2 (PATCH: clarify public use versus managed-agent execution scope)
- Principle IV: Debug-only Slurm is a managed-agent/private-deployment rule, not a public prerequisite.
- Public users may compile and run directly using README shell commands.
- Unchanged: this session's agent execution boundary, FP64, physics, tolerances and simple .slm examples.
- Updated: AGENTS.md; .agents/skills/speckit-{implement,tasks}/SKILL.md;
  .specify/templates/{plan,spec,tasks}-template.md; feature scope, contracts and quickstart.
- Added/removed principles: none; principle IV scope clarified; no pending propagation.
- Historical anonymous Debug validation facts retained; no solver changes or execution.
-->
# polygr-pfc Constitution

## Core Principles

### I. 仅限 DCU 加速

本项目通过相场晶体（PFC）方法生成多晶石墨烯结构，并转换为 LAMMPS 可用的原子初始构型。
所有代码改写 MUST 仅服务于既有 PFC 程序的 DCU 加速及其必要的构建、调度、兼容与验证。
允许迁移计算内核、FFT、归约、内存管理及 MPI/DCU 数据交换；每项修改 MUST 说明与加速的直接关系。
无关重构、新物理模型、新业务功能及可视化或格式转换工具的独立改造 MUST NOT 混入迁移。

### II. 物理算法与公式不变

原有 CPU 实现及输入模板 MUST 作为可追溯的科学基准。
自由能泛函、演化方程、守恒与非守恒动力学、线性谱算子、非线性项、周期边界、
FFT 归一化、时间推进方法、计算盒优化准则、初始化语义、参数含义及单位 MUST 保持不变。
不得通过改变时间步长、网格、迭代次数、物理参数或收敛标准来制造加速结果。

执行位置、内存布局、并行划分与 FFT 后端可以改变，但 MUST 实现同一离散数学运算，
保持原有双精度数值表示；MUST NOT 引入降精度、近似物理项或替代积分算法。
并行归约和 FFT 实现带来的浮点舍入差异 MUST 通过预先规定的容差验证，不能借此接受公式变化。
发现原有物理算法问题时 MUST 单独记录，不能在 DCU 移植中顺带修正。

### III. 两阶段工作流与数据兼容

MUST 保持第一阶段晶粒生长、第二阶段读取前一阶段状态并弛豫的顺序与科学含义。
既有命令行调用、输入格式、状态文件及重启语义、坐标与邻居数据、LAMMPS 输出约定 MUST 兼容。
DCU 资源选择只能增加必要的执行配置，MUST NOT 改变输入中物理参数的解释。
密度场到原子坐标及 LAMMPS 初始构型的转换链路 MUST 保持可用，
以避免求解器加速后产生下游工具无法使用的结果。

### IV. 区分公共使用与受管代理执行

公共用户 MUST 能按 README 中的直接 Shell 命令编译、运行和后处理。
Slurm、Debug 分区及私有站点调度规则 MUST NOT 成为公共安装或使用前提。

本会话代理在受管集群执行依赖/程序编译、语法或可执行检查、测试、数值回归、
性能测试、试算及转换验证时，MUST 仍经该集群提交 Slurm 批作业，在分配的计算节点执行。
代理使用的 `.slm` MUST 显式声明 `#SBATCH --partition=Debug`；
MUST NOT 改用其他分区、分区列表或自动回退，也 MUST NOT 在本机或登录节点直接执行上述任务。
这项执行限制属于受管代理及私有部署，不限制公共用户自己的执行环境。

提供的 `.slm` 示例 MUST 使用 `#!/bin/bash -e` 及顺序命令，资源由 `#SBATCH` 声明。
MUST NOT 重新加入 `set -euo pipefail`、`${VAR:?}`、分区/节点/进程数的重复Shell检查，
或用Shell分支、包装函数隐藏作业流程。
日常 `.slm` MUST 仅保留构建/运行命令，MUST NOT 调用tests、
evidence/metadata生成、`/usr/bin/time`、设备探针或调试验证；验证材料与日常运行路径分开。
既有科学与性能验收记录 MUST 保留；公开摘要不包含账户、主机、作业编号或部署路径。
移除Shell外层计时封装 MUST NOT 删除求解器原科学逐步计时日志。

代理本机仅可用于编辑、静态文本审阅、文件同步及受管集群作业管理。
代理DCU作业 MUST 申请实际可用设备并遵循该部署的Debug资源配置；
公共用户自行提供兼容的MPI/DCU环境，公共命令不得强制依赖私有调度设施。

### V. 数值等价与加速证据

求解器迁移 MUST 保留可追溯的原 CPU 基准版本，并使用相同输入、固定随机种子和相同初始场
进行对照；MPI 进程布局及其对初始化和归约的影响 MUST 记录并控制。
规格与计划 MUST 在运行前规定误差指标、绝对/相对容差及其依据，禁止根据失败结果放宽容差。
验证 MUST 覆盖单步或短程密度场、平均密度与自由能、两阶段衔接及最终原子构型转换；
守恒性检查 MUST 匹配各阶段实际采用的动力学，不能要求非守恒阶段保持平均密度不变。
非有限值、超限误差、重启或转换不兼容 MUST 阻止验收。

性能结论 MUST 在数值验证通过后，以相同科学工作量的 CPU/DCU 对照为依据。
MUST 区分内核时间与包含初始化、FFT、通信、数据传输和 I/O 的求解器总时间，
记录硬件、软件版本、编译选项、输入、资源配置、退出状态与日志路径；使用调度器时记录作业身份。
没有实际运行证据时 MUST 明确标记为未验证，不能宣称已正确运行或获得加速。

## Scope and Execution Constraints

- `src/pfc.c` 是当前两阶段 MPI/FFTW 求解器入口；`source/` 包含其他求解器及 Java 工具。
  计划 MUST 确认实际修改范围，不得将整个目录中的工具默认纳入重写。
- DCU 编程环境与 FFT 库的具体选择 MUST 在目标环境核实后进入计划；本宪章不预设库版本、
  设备型号、设备数量或多卡能力，也不代表 DCU 后端已经实现。
- 构建和运行 MUST 使用独立工作目录，保留基准、输入及日志；MUST NOT 将生成数据和无关二进制
  混入源码提交。运行目录 MUST 与源码分离，保留既有输入和结果。
- 通用部署与模块配置见 `README.md`；代理操作遵循 `AGENTS.md`。

## Development Workflow and Quality Gates

1. **规格关卡**：列出加速目标、范围外事项、物理不变量、兼容接口、数值容差与可测量性能目标。
2. **计划关卡**：明确 CPU 基准、DCU 工具链、运算到原公式的映射、资源和编译/验证方案；
   代理执行遵守受管Debug批作业边界，公共用户可使用直接Shell命令；研究前及设计后逐项通过 Constitution Check。
3. **任务关卡**：将构建、CPU/DCU 数值对照、两阶段转换回归和性能测量列为迁移必需任务。
   文档变更仅需静态审阅；不因模板示例将迁移验证标为可选。
4. **实现关卡**：每项代码改写说明加速目的及保持的物理运算；代理的构建和所有测试只通过受管Debug批作业。
5. **验收关卡**：审阅最小必要差异及实际执行证据，数值与兼容性先通过，再报告性能。
   无法运行或证据不足时保持任务未完成，明确阻塞原因。

## Governance

本宪章优先于仓库中冲突的模板、历史命令示例与开发惯例；执行时仍须遵守用户的明确指令。
每份规格、计划、任务清单及 PR MUST 核查五项原则，记录数值影响、执行位置与验证证据。
不得使用 Complexity Tracking、性能收益或测试便利性豁免 DCU 范围、物理不变和受管代理执行要求。

修订 MUST 记录动机、条款变化、影响及同步文件；改变既有约束 MUST 先取得项目维护者明确授权，
再正式修订宪章及全部受影响模板和指引。版本采用语义化版本：原则删除或不兼容重定义升 MAJOR，
新增原则或实质扩展升 MINOR，文字澄清与非语义修正升 PATCH。
每次修订 MUST 更新顶部 Sync Impact Report、版本与最后修订日期，保留首次批准日期。

**Version**: 1.2.2 | **Ratified**: 2026-10-20 | **Last Amended**: 2026-10-21
