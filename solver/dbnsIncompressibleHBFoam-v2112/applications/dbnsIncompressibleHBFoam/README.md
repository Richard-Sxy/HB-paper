# dbnsIncompressibleHBFoam

## 一次螺旋桨计算，核心 .H 分别做了什么？

九个相位怎么建立？| createMeshes.H updateMeshes.H | 前者读取 instantsNumber=9，创建九套网格并确定采样时刻；后者把各个相位网格移到对应位置，计算运动所需的网格通量
程序从哪里拿到压力、速度和湍流场？| initialise.H createFields.H | initialise.H 读取控制参数和容限，并引入 createFields.H。后者为每个相位建立 p、'U'、通量和湍流模型，还会检查续算所需的 ohiAbsolute 是否齐全。
单个相位的方程怎么样写出来？| preparePhaseEquations.H、'residualsUpdate.H' | 当前建立动量对流项、粘性项和压力矩阵块；后者计算这个相位当前还差多少才满足方程。之前查到的连续性缺陷，就在 residualsUpdate.H:21 的 fvc::div(phiAbsolute[K])
九个相位在哪里连成一个问题？| outerLoop.H | 他先组装各相位内部的空间方程，然后再第 314 行调用 HB.addBlock(hbVelocity, J, K)，把相位 K 的速度对相位 J 的动量贡献加入系统。
内层 Newton 在修正什么？| fullStageDefectCorrection.H、stageNewtonKrylov.H、stagePressureSchur.H | 第一个文件定义候选压力、速度、通量增量对应的阶段残差；第二个求能降低该残差的 Newton 方向；第三个提供当前选用的压力 Schur 预条件器，帮助线性迭代求这个方向。
算出了方向，这一步能走多远？| applyFlowLineSearch.H | 它用同一个 alpha 更新九个相位的候选场，再重新计算物理指标。先打印试探结果，之后判断是否接受。日志中的 IHBFlowTrail 就来自这边。
螺旋桨的 k 和 epsilon 在哪里更新？| HBKEpsilon.H、correctHBTurbulence.H、correctPhaseImplicitTurbulence.H | 第一个提供 k-ε 模型的 HB 方程；第二个按字典选择更新方法
最后如何宣布收敛？| evaluateFlowResiduals.H、evaluatePhysicalResiduals.H | 前者检查动量、连续性、压力参考和通量一致性；后者再合并湍流残差。

## src/ICSFoam/ 下面的文件有什么用？

HBZone.H 定义一个 HB 区域保存哪些数据，包括频率、相位变换矩阵和时间微分算子 D。
coupleMatrix.H 定义整个耦合系统能存哪些矩阵块。相位 J 的动量方程和相位 K 的速度之间的联系。
blockFvMatix.H 定义单个有限体积矩阵块怎么样保存单元对角、相邻单元及边界贡献。

## 求解器的所有文件
### 可执行算法步骤

HBKEpsilon.H | k-ε 模型，湍动能k、耗散率ε、空间输运、壁面条件
HBKOmegaSST.H | k 和 omega 的相位输运与残差计算接口
HBRotatingFrame.H | dynamicMeshDict 指定转速
HBTurbulenceStages.H | 湍流阶段计算工具
applyFlowLineSearch.H | 外层流场更新能不能接受。同一个 alpha 更新所有相位的候选压力、速度、绝对面通量
boundLocalTimeStep.H | 用密度、动量、能量等变量限制局部伪时间更新，不可压缩没有直接包含它
calculateScalarTransport.H | 实现额外标量的全相位输运，不可压缩变量没有直接包含它
compressibleCourantNo.H | 将声速与相对网格速度计入可压缩库朗数。不可压缩应用计算伪时间尺度使用 setCoAndDelta.H
correctHBTurbulence.H | 湍流更新的分支入口。读取 hbTurbulenceMethod，进入全隐式/显式更新
correctPhaseImplicitTurbulence.H | phaseImmplicit 设置，进入这里
createFields.H | 创建并读取各相位的场。同时检查 phiAbsolute 是否齐全
createMeshes.H | 建立 HB 相位网格。fvSolution 读取 instantsNumber，随后创建相位网格、确定采样时刻
evaluateFlowResiduals.H |  当前流场还差多少能满足方程。日志的 IHBPhysicalResidual
evaluatePhysicalResiduals.H | 流场残差基础加入湍流残差