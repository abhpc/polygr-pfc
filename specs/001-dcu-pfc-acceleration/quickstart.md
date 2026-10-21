# Quickstart

公共用户的依赖准备、直接Shell编译、MPI运行和后处理命令见[项目 README](../../README.md)。
不要求Slurm、Debug分区或私有集群模块。可选调度部署见[Host.md](../../Host.md)。

科学输入、输出兼容与验收标准分别见 [spec.md](spec.md)、[contracts/](contracts/) 和
[匿名验证摘要](validation.md)。历史实测在Slurm Debug节点完成，不构成公共运行前提。

本会话代理的编译、语法检查、测试、计算及后处理仍只能在受管集群Debug批作业内执行，
禁止本机或登录节点直接执行。提供的日常 `.slm` 保持顺序构建/运行命令，
不自动调用测试、探针、证据生成或外层计时。
