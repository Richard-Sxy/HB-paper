# Overleaf 双机写作

- JHD 论文以已上传并编译成功的 Overleaf 项目为主稿。两台电脑用浏览器登录同一个个人账号，打开同一项目地址。免费版即可完成双机写作，无需购买会员、下载模板或配置 API。
- 让 Codex 修改正文时，在当前电脑已登录的浏览器中打开该项目，核对文件树中的 `paper.tex` 和待改段落，再按用户要求编辑。项目地址可从已编译成功的编辑页地址栏复制。
- 换电脑后，先打开刚编辑的源文件，确认上一处改动已经出现，再继续写。Overleaf 自动保存源码；需要查看排版时，点击 **Recompile**。

**首次双向验证（待实测）**：在电脑 A 的 `paper.tex` 顶部单独添加一行 `% sync-test`。电脑 B 看到后删掉这一行，再回到 A 确认它已消失。若改动没有出现，先核对两端的项目地址、打开的文件和网络连接。

# SCNet 连接与操作范围

- 在本机 PowerShell 中运行 `ssh scnet`，通过已配置的 SSH 别名连接 SCNet 华东三区（乌镇）。
- **远程文件操作仅限 `~/sxy/Paper/` 及其子目录。** 登录后先运行 `cd "$HOME/sxy/Paper"`，失败就停止操作。
- 在同一个 SSH 会话中，进入该目录后先运行 `of2112` 加载 OpenFOAM v2112。用 `echo "$WM_PROJECT_VERSION"` 确认输出为 `v2112` 后，再执行 OpenFOAM 相关命令。
- 读取、列出、创建、修改、删除及传输文件时，涉及的远程路径都必须位于该目录内。禁止通过 `..` 或符号链接访问目录外的文件，包括 `~/sxy/` 下的其他目录。
- **即使在 `~/sxy/Paper/` 内，也禁止执行高风险命令或代码。** 包括递归强制删除、通配符批量删除或覆盖、提权、修改系统配置、执行来源不明的脚本，以及可能影响其他用户或集群服务的操作。无法确认操作安全时，停止执行并说明原因。

# 启用脚本的执行方式

- mpirun -np 16(并行核数) interPhaseChangeDyMFoam_cor(这个是求解器的名称) -parallel | tee log
- 作业上传之后可能会等待排队，这个不需要管，之后会自己执行计算
- 在计算完成后，执行 reconstructPar 构建时间场，然后使用命令 touch case.foam 创建一个 foam 文件后续给 Paraview 能阅读到。

# 算例计算和下载细节

- 算例计算的结果和报错在文件夹下的 log 文件能看到，运行项目需要没有任何报错，如果有报错就需要重点排查报错的问题。

# 脚本提交指令模板

命令行提交指令会在命令行界面运行，偶尔会出现问题，用脚本提交指令就稳定运行：
#!/bin/bash
#SBATCH -J [这边放入作业名，一般用修改的参数命名，下划线隔开，不要加多余的字符]
#SBATCH -p wzhcnormal
#SBATCH -N [这边放入计算节点数量]
#SBATCH --ntasks-per-node=[这边放入申请CPU核数]
#SBATCH --cpus-per-task=[每个任务申请一个CPU]
module purge
module load compiler/devtoolset/7.3.1
module load mpi/openmpi/4.1.1/gcc-7.3.1
source /work/home/acrkaz0qxz/OpenFOAM/OpenFOAM-v2112/etc/bashrc

srun --mpi=pmix_v3 [这边放入求解器名称] -parallel 2>&1
这边可以放入其他指令，比方说如下
reconstructPar [这边放入这个指令的构建参数]
