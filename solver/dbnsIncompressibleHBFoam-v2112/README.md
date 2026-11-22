# dbnsIncompressibleHBFoam

OpenFOAM-v2112（OpenCFD / ESI）的不可压缩谐波平衡求解器。水翼和螺旋桨使用同一个程序，相位数、边界条件及求解控制由各自的算例字典指定。

## 适用范围与算例设置

本程序用谐波平衡法计算具有已知周期的单相不可压缩流动。网格、运动、边界条件和初场由 OpenFOAM 算例提供。支持层流，以及 `kEpsilon`、`kOmegaSST` 两种 RAS 模型。

| 设置位置 | 作用 |
| --- | --- |
| `system/fvSolution` 中的 `harmonicBalance` | `instantsNumber` 指定相位数，默认 3。`selectedPeriod` 指定采样周期。相位数须能表示所选谐波。 |
| `constant/HBProperties` | 指定 HB 区域及谐波。此实现要求活动 HB 区域覆盖全网格（`allMesh true`），使用笛卡尔速度分量（`cylCoords false`）。 |
| `system/controlDict` | `startTime`、`endTime` 和 `deltaT` 管理外层伪时间推进。物理振荡周期由 HB 设置给定。 |

选择 `pressureSchur` 预条件器时，压力增量边界支持 `fixedValue`、`zeroGradient`、`empty` 和耦合边界。可选调试开关 `stageFoilAudit` 默认关闭。它专用于历史水翼算例的两侧 AMI 和指定局部单元取样，常规算例保持默认值即可。

同一程序曾运行三相位层流水翼（周期 1.25 s）和九相位 `kEpsilon` 螺旋桨（周期约 0.039767 s）。两份已有运行记录结束时均报告 `converged=0`。其他工况的收敛和精度应按该工况的物理残差、载荷和参考结果核实。

## 目录

| 路径 | 内容 |
| --- | --- |
| `applications/dbnsIncompressibleHBFoam/` | 最终求解器源码及其 `Make/` 配置，入口为 `dbnsIncompressibleHBFoam.C` |
| `src/ICSFoam/` | 配套库源文件及其 `Make/` 配置，已包含本版本使用的兼容性修改 |
| `platforms/linux64GccDPInt32Opt/` | 已编译的 `bin/dbnsIncompressibleHBFoam` 和 `lib/libICSFoam.so` |
| `activate.sh` | 在当前 Bash 终端中启用本包的程序和库 |
| `Allwmake` | 在本包内先编译库，再编译求解器 |
| `LICENSE` | GNU GPL v3 文本 |

## 运行

在 Linux 或 WSL 的 Bash 中解压，先加载本机 OpenFOAM-v2112 环境。下面的环境脚本路径适用于此次验证的 Ubuntu 安装，其他安装请换成实际路径。

~~~bash
tar -xzf dbnsIncompressibleHBFoam-v2112.tar.gz
cd dbnsIncompressibleHBFoam-v2112
source /usr/lib/openfoam/openfoam2112/etc/bashrc
source ./activate.sh
command -v dbnsIncompressibleHBFoam
dbnsIncompressibleHBFoam -help
~~~

随包程序适用于兼容的 Linux x86-64、GCC、双精度、32 位标签 OpenFOAM-v2112 环境。程序路径应指向解压目录里的 `platforms/`。更换系统、编译器或 MPI 环境时，在目标机加载 v2112 环境后，进入包根目录运行：

~~~bash
./Allwmake
source ./activate.sh
~~~

成功后在同一终端进入已配置好的 HB 算例。以下示例使用已经分为四个子域、具有完整各相位初场或重启场的算例：

~~~bash
cd /path/to/decomposed-hb-case
mpirun -np 4 dbnsIncompressibleHBFoam -parallel
~~~

算例的网格、各相位场和控制字典由使用者提供。运行日志会显示各伪时间迭代及外层接受情况。

## 验算与退出状态

本包已在 WSL Ubuntu 20.04 的 OpenFOAM-v2112 环境从源码编译，并使用随包程序和库在四核三相位水翼算例上从伪时间 0.3 续算至 0.4。外层接受全步（α = 1），最大物理动量残差从 7631.31006329 降至 6677.81329644。三个相位的力与力矩在日志输出精度下均复现了原运行记录，四个子域均写出了各相位的速度和压力场。

程序以 `IHBFinished converged=1` 和退出码 `0` 表示达到设定的完整物理收敛判据。本次单步验算到达指定终点时记录 `IHBFinished converged=0`，对应退出码 `2`。提前停滞也可能以 `2` 退出，可用日志中的 `IHBStagnated` 区分。

## 常见报错与排查

先看失败发生在哪一步。在算例目录中，分区失败看 `log.decomposePar`，求解开始后看求解日志。以下命令假定求解日志名为 `log.hb`。并行日志行首的 `[2]` 等数字是 MPI 进程编号。

~~~bash
# 找到第一处致命错误，并显示它之后的文件名、字典条目或相位号
grep -n -m 1 -A 15 -E 'FOAM FATAL|error while loading shared libraries|command not found' log.hb
# 已经进入迭代时，查看最近的试探步、物理残差和退出状态
grep -E 'Iteration time =|IHBFlowTrial|IHBFlowAccepted|IHBFlowRejected|IHBPhysicalResidual|IHBTurbulenceResidual|IHBStagnated|IHBFinished' log.hb | tail -n 35
~~~

### 启动、分区与初始场

| 看到的提示或现象 | 原因与处理 |
| --- | --- |
| `Load the OpenFOAM-v2112 environment before sourcing activate.sh.` 或 `wmake: command not found` | 在当前 Bash 里加载目标机 OpenCFD/ESI v2112 的 `etc/bashrc`，用 `echo "$WM_PROJECT_VERSION $WM_OPTIONS"` 查看版本和编译选项，然后在包根目录运行 `source ./activate.sh`。 |
| `dbnsIncompressibleHBFoam: command not found`、`libICSFoam.so` 找不到或出现 `undefined symbol` | 用 `command -v dbnsIncompressibleHBFoam` 确认程序路径。找到程序后运行 `ldd "$(command -v dbnsIncompressibleHBFoam)"`，检查 `libICSFoam.so` 是否来自此包、其余 OpenFOAM 库是否来自当前 v2112 安装。如果编译环境、MPI 或 ABI 不同，按“运行”一节在目标机重新构建。 |
| 启动后提示找不到 `processorN`、网格或 `subTimeLevelK/U`、`p` 等文件 | 先看 `log.decomposePar` 第一处错误。四核运行应有 `processor0` 至 `processor3`，与 `numberOfSubdomains 4` 和 `mpirun -np 4` 一致。再按实际起始时刻检查每个核的 `processorN/<起始时刻>/subTimeLevelK/`：三相位的 K 为 0–2，九相位为 0–8。层流需有 `U,p`；本文的螺旋桨还需有 `k,epsilon,nut`。从零起算应选择 `startFrom startTime`、`startTime 0`；从检查点起算则选择所有相位都写全的时刻。 |
| `phiAbsolute must be present for every HB phase and every rank` 或 `Invalid phiAbsolute field class` | 某些核或相位混用了不同的初始通量，或文件类型有误。冷启动用干净的初场目录，让所有核和相位都由求解器建立绝对面通量；要复现续算，则从同一检查点恢复整套 `phiAbsolute`，类型为 `surfaceScalarField`。日志中的 `IHBInitialAbsoluteFlux read=0` 是正常的冷启动提示。 |
| `Expected kinematic pressure p` 或 `Expected volume flux` | 核对场文件的 `dimensions`。本程序的 `p` 是运动学压力 `[0 2 -2 0 0 0 0]`，`U` 为 `[0 1 -1 0 0 0 0]`，`phi` 和 `phiAbsolute` 为体积通量 `[0 3 -1 0 0 0 0]`。从其他求解器迁移时应转换物理量及其边界值，单改量纲标签会留下错误的数值。 |

### 相位、周期与模型配置

| 报错关键字 | 应核对的文件与设置 |
| --- | --- |
| `HB requires instantsNumber >= 2` | 检查 `system/fvSolution` 中的 `harmonicBalance.instantsNumber` 和实际生成的相位数。本文水翼用 3 相位、1 阶谐波，螺旋桨用 9 相位、4 阶谐波；相位目录须与配置对应。 |
| `cycles over the sampled period` 或 `use selectedPeriod = 2*pi*n/omega` | 核对 `constant/HBProperties` 的角频率（rad/s）和 `system/fvSolution` 的 `selectedPeriod`（s）。**实际螺旋桨示例**：角频率为 158 rad/s，把一转周期截成 `0.039767` 后，日志报告 `1.00000011 cycles` 并退出。该工况采用 `0.039766995615060675`。编辑后重新读取字典，确认有效数字保留下来。 |
| `Mesh motion is not periodic over the HB sample period` | HB 周期对应的首末网格位置未通过检查。对照 `constant/dynamicMeshDict` 的运动函数、轴、角频率和 HB 周期。水翼的 `omega=5.026548245743669 rad/s` 对应 `T=1.25 s`；先确认运动定义和物理周期，再判断网格偏差。 |
| `set allMesh true` 或 `cylCoords true is not supported` | 检查 `constant/HBProperties` 中每个活动 HB 区域。本实现要求 `allMesh true;` 和 `cylCoords false;`。 |
| `Use the same RAS model in every HB phase` 或 `Supported HB RAS models` | 核对根目录及各相位的 `turbulenceProperties`。水翼用 `laminar`；螺旋桨用 `RAS`、`kEpsilon`、`turbulence on`。也支持 `kOmegaSST`，其第二个输运场为 `omega`。同一算例各相位的模型和场文件应一致。 |
| `pressureSchur only supports` | 日志会指出相位、边界名称和实际压力边界类型。检查该相位的 `p`。此预条件器支持 `fixedValue`、`zeroGradient`、`empty` 和耦合边界；修改边界类型前应确认它符合该边界的物理条件。 |
| `stageNewtonSteps requires fullStageCorrectors>0 and flowLineSearch` | 在 `system/fvSolution` 检查阶段修正设置。启用内层 Newton（`stageNewtonSteps > 0`）时，须同时启用外层 `flowLineSearch` 并设置正的 `fullStageCorrectors`。 |

### 已经进入迭代，但更新被拒绝或计算停止

| 日志标志 | 如何判断和处理 |
| --- | --- |
| `IHBFlowTrial`、`IHBFlowAccepted` | `alpha 1` 试全步，小于 1 表示缩步。比较同一试探步的 `merit` 和 `target`，并看随后的 `IHBFlowAccepted`。某个试探步超出目标后，后续更小的步长仍可能被接受。 |
| `IHBFlowRejected` | 外层线搜索没有找到可接受步，程序恢复上一流场。对照这一轮每次试探的 `alpha`、`merit`、`target`、`momentum`、`continuity`、`flux` 和 `pressureStep`。若内层出现 `IHBStageNewtonAccepted`，它只确认阶段修正；外层接受情况要看 `IHBFlowAccepted`。 |
| `IHBStagnated repeated minimum-Co rejected steps` | 伪 Courant 数已到下限，连续拒步次数达到 `flowRejectedStepLimit`。旧水翼调试算例曾在连续 3 次拒步后停止。此时应回看**首次拒步**前后的动量、连续性和压力增量；单独延长 `endTime` 无法改变当前的更新方向。 |
| `IHBFinished converged=0`，随后 MPI 报退出码 `2` | 求解器已执行最终写出，完整物理收敛判据尚未满足。先看末尾是否有 `IHBStagnated`，再核对最后的伪时间、`endTime` 和 `stopAt`。短跑到达设定终点也会出现此结果。若每步仍被接受、物理残差持续下降，可以从完整检查点继续；若反复拒步，按上一行查原因。 |
| 只出现 `MPI_ABORT`、非零退出或作业被终止，找不到 `IHBFinished` | 从日志开头定位第一处 `FOAM FATAL` 或库加载错误及其后的具体信息。若求解日志中没有原因，再查看作业标准错误和调度器的内存、时限记录。单凭 MPI 的最后一行无法确定故障类别。 |

## 来源

本版本沿用 ICSFoam 的配套库，保留各源文件的署名。上游 ICSFoam 按 GPL v3 或更新版本分发，源自以 HiSA 为起点的工作。上游推荐引用 [Computer Physics Communications 论文](https://doi.org/10.1016/j.cpc.2023.108673) 和 [Machines 论文](https://doi.org/10.3390/machines10040279)。
