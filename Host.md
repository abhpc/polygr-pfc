# Slurm 集群部署示例

本文说明可选的环境模块和 Slurm 部署方式，直接编译、运行命令见 [README.md](README.md)。个人 SSH 地址、端口、账户、站点路径、模块清单和作业记录保存在本地 `yhmsr/` 目录；该目录已被 `.gitignore` 忽略，不随仓库发布。

## 查询模块

在目标集群的登录会话中执行：

```bash
module av
module av dtk
module av mpi
module list
module show dtk/26.04
```

`module av` 是 `module avail` 的缩写。模块名和版本由集群管理员提供；本文中的版本是部署示例，不能代替本站模块清单。

## 加载 PFC

先在 `slurm/modulefiles/pfc/2026.1` 中设置本站的 GCC、Intel MPI 和 DTK 模块名，再从仓库根目录执行：

```bash
module use "$PWD/slurm/modulefiles"
module load pfc/2026.1
```

modulefile 使用 Environment Modules 提供的 `ModulesCurrentModulefile` 定位自身，再计算项目根目录，设置 `PFC_ROOT`、`FFTW_PREFIX`、程序及运行库路径，并设置单线程 MPI 环境。工具链需提供 `I_MPI_ROOT`、`DTKROOT` 和 `HIP_PATH`；Makefile 不保存站点安装路径。变量含义见 [Modules 官方文档](https://modules.readthedocs.io/en/stable/modulefile.html#modulescurrentmodulefile)。

此 modulefile 应留在仓库的 `slurm/modulefiles/pfc/2026.1`。若复制到集中管理的模块目录，管理员需将其中的 `root` 改为实际共享安装目录。每个新登录会话先注册模块目录，再加载 PFC；也可把 `module use` 写入自己的 Shell 登录配置，个人配置不应提交到仓库。

后处理分别使用 Python、Java 和 OVITO 环境。按本站模块名修改 `slurm/polygr.slm` 和 `slurm/convert.sh`；保留原 Coordinator JAR 及 OVITO 2.9.0 的转换接口。

## 获取 Slurm 参考模板

提供 `slmscript` 的集群可以先查看并复制参考脚本：

```bash
slmscript ls
slmscript cat slm.vasp-dcu
mkdir -p ~/pfc-slurm-template
cd ~/pfc-slurm-template
slmscript cp slm.vasp-dcu
```

`slmscript cp` 只把模板复制到当前目录，不提交作业。VASP 模板仅用于参考资源声明和模块加载方式；PFC 应使用本仓库的 `.slm` 脚本。模板中的账户、分区、路径和应用命令必须按项目修改。

## 提交示例

仓库附带脚本的 `#SBATCH --partition=Debug` 来自特定集群的部署约定。公共程序可直接编译和运行；使用 Slurm 的站点需要按实际调度规则配置分区、账户和资源。

公开脚本不写个人账户或 QoS，使用集群默认值；需要指定时，保存在自己的作业副本中。单节点示例申请 32 个 MPI rank、4 张 DCU。不同 DCU 架构需修改 `build.slm` 的 `PFC_GPU_ARCH`；跨节点配置见 README。

资源写在 `#SBATCH` 中，命令按顺序写在 `.slm` 中，提交时只需：

```bash
sbatch /path/to/polygr-pfc/slurm/build.slm
sbatch /path/to/polygr-pfc/slurm/polygr.slm
```

先确认构建完成，再提交计算。登录节点的使用遵循所在站点的规定。完整安装与直接运行步骤见 [README.md](README.md)。
