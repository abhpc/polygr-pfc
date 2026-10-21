# Quickstart

构建、依赖模块、批处理资源和日常运行命令见[项目 README](../../README.md)。
科学输入、输出兼容与验收标准分别见 [spec.md](spec.md)、[contracts/](contracts/) 和
[匿名验证摘要](validation.md)。

所有编译、测试和计算仅在目标集群的 Slurm Debug 批作业内执行。
日常 `.slm` 只保留顺序构建/运行命令，不自动调用测试、探针、证据生成或外层计时。
