# HB 求解器内的螺旋桨固壁声学计算

`dbnsIncompressibleHBFoam` 与 `dbnsFullyImplicitHBFoam` 均有声学入口。将示例字典作为算例的 `system/hbAcousticsDict` 后，求解器达到自身的收敛判据时，在同一个 OpenFOAM 进程内调用 `libAcoustics.so`。按伪迭代数停止的未收敛算例不会启动声学重建。字典缺省或 `enabled false` 时按原流程计算流场。声学是从流场到观察点的单向计算。

## 传给 libAcoustics 的数据

1. 从当前相位网格的 `p.boundaryField()` 取同一材料面的压力，保留平均载荷。运动学压力乘以 `rhoInf` 转成 Pa，绝对压力直接以 Pa 读取，两者最后都减去以 Pa 表示的环境压力 `pInf`。
2. 从当前 ICSFoam 全网格 HB 区域读取频率和矩阵 `E`。先把面压力投影到保留的模态，再按指定物理时间合成桨面压力。程序核对每个相位的面中心、面积向量和给定刚体运动，包括求解器计算网格通量时的微小时间偏移。
3. 每个声源时刻用旋转后的面中心、从固体指向流体的面积向量和刚体速度计算声源到观察点的距离及到达时刻。对不可渗透壁面取流体与壁面相同的法向速度，压力承担应力载荷。适配器把这些数据传给上游 `Farassat1AFormulation::observerAcousticPressure`，由库求厚度和载荷源项、插值到接收时刻并做 MPI 求和。本地补丁仅为适配类开放两个上游类的数据访问权限。
4. 程序根据观察距离、声速和差分启动所需步数预热，随后输出完整的接收记录。观察点坐标为惯性系坐标，时间列以秒表示。启用叶频分析时，对该记录做离散傅里叶投影。

本实现覆盖单个刚体恒角速度、全湿固壁、一个全网格 HB 区域及均匀静止的无界声传播介质。选定的多个 patch 合并为一个积分面；实际桨的叶片、桨毂和截面覆盖需结合同一套物理模型核对。可压缩流场的静压同样可作为壁面输入。可渗透面、空化体积源及有均匀对流的声传播需要相应的声源和传播模型。[libAcoustics v2112 源码](https://github.com/unicfdlab/libAcoustics/tree/digitef-dev-2112) · [Abergo 等的 HB 桨面重建与 FW–H](https://doi.org/10.1016/j.ast.2024.109343)

## 在 OpenFOAM v2112 中构建

下面的 Bash 命令从本项目根目录执行，使用与两个求解器一致的 OpenFOAM v2112、编译器及 MPI 环境。本地可使用已准备好的上游工作副本。其他机器可检出 `unicfdlab/libAcoustics` 的 `digitef-dev-2112` 分支，并调整第一条构建命令的源码目录。构建脚本从固定提交 `1ac0140e4dcfe8c8add29293c793fc57d73854ec` 导出源码，在临时目录应用补丁。

```bash
source /usr/lib/openfoam/openfoam2112/etc/bashrc
apps="$PWD/solver/ICSFoam/applications/solvers"
ac="$apps/HBacoustics"
ihb="$apps/dbnsIncompressibleHBFoam"
bash "$ac/build-lib.sh" "$PWD/05_开发/.libAcoustics-v2112-reference" "$ac/platforms/$WM_OPTIONS/lib"
bash "$ac/build-solvers.sh"
export LD_LIBRARY_PATH="$ac/platforms/$WM_OPTIONS/lib:$ihb/platforms/$WM_OPTIONS/lib:$LD_LIBRARY_PATH"
export PATH="$ac/platforms/$WM_OPTIONS/bin:$PATH"
```

`build-solvers.sh` 复用本项目不可压缩求解器兼容版的 ICSFoam 头文件和共享库。初次构建 ICSFoam 时可先运行 `bash "$ihb/build.sh"`。若这套兼容库和头文件装在其他位置，分别设置 `ICSFOAM_LIB_DIR` 和 `ICSFOAM_BUILD_SRC`，并将实际库目录加入 `LD_LIBRARY_PATH`。

## 配置与结果

从 `hbAcousticsDict.example` 取字典模板，将 `wallPatches`、`origin`、带方向的角速度向量 `omega`、介质参数和观察点按实际工况填写。默认输出名含最终伪迭代时刻。可以用 `outputName` 区分多套观察点配置。

```bash
cp "$ac/hbAcousticsDict.example" "$CASE/system/hbAcousticsDict"
cd "$CASE"
dbnsIncompressibleHBFoam
```

可压缩算例在算例目录启动 `dbnsFullyImplicitHBFoam`。并行运行沿用相应 HB 算例原有的分区、相位场和 `-parallel` 启动方式。字典由求解器结束时主动读取，勿将声学对象再加入 `controlDict.functions`。

不可压求解器也可从已有检查点单独计算声学。将 `system/controlDict` 中的起始时刻选为目标检查点，保留求解时的基础网格、相位配置和每个相位的 `p`。动网格还需保留该时刻各相位的 `polyMesh/points`。在 `system/hbAcousticsDict` 中启用声学，给不同的观察点配置设置不同的 `outputName`。串行运行 `dbnsIncompressibleHBFoam -replayHBAcoustics`。已分区的算例使用原分区数，例如四分区运行：

```bash
mpirun -np 4 dbnsIncompressibleHBFoam -parallel -replayHBAcoustics
```

查看日志中的 `IHBAcousticReplay checkpoint` 核对读入时刻。重放直接使用指定检查点的压力。判断流场是否收敛应查原流动日志的 `IHBFinished converged=`。

- `acousticData/<outputName>-time.dat` 给出接收时间和各观察点的总脉动声压，单位 Pa。参考零点 `pInf` 必须与输入压力的定义一致，绝对空气静压与水动力表压通常取不同值。
- 可选 `tones { bladeCount 5; orders 3; referencePressure 1e-6; }` 中的叶数须按实际桨填写。它在 `acousticData/<outputName>-bpf.dat` 给出阶次、频率、总声压 RMS、相对于指定参考声压的声压级和相对于接收记录首点的相位。水声参考通常取 1 µPa，空气声通常取 20 µPa。
- `samplesPerPeriod` 和 `periods` 使用流场 HB 基本周期。一转对应几个 HB 周期由流场设置决定。评价指向性时可配置多个观察点；求叶频时记录应含整数次叶通过，且采样应解析所需的最高频率。增加采样数只细化现有流动谐波的时间重建。

## 运行受控声源检验

```bash
python3 -B "$ac/smoke_test.py"
```

脚本在临时目录调用真实求解器及动态库，检查静止面元的解析偶极子声压、旋转面元与独立 Python Farassat 1A 实现的总声压、串并行结果、绝对压力的单位转换及原生叶频输出。可压缩测试在壁面规定常量绝对压力，并对这个测试放宽归一化残差阈值，以单独检验终止回调和静压基准。生产算例沿用其流场收敛判据。
