压力—速度—面通量联合 Newton 原型（2026-11-14）

状态
  已实现可选路径；默认关闭。尚未在 OpenFOAM 中完成编译或运行。
  本机没有配置好的 v2112 工具链；本地源码 instantList.H 引用的
  instant.H 缺失（同目录仅找到 Instant.H），源码头文件语法检查止于此。
  不将脚本生成成功、静态检查或备份校验视为数值验证。

修改主目录
  ICSFoam/applications/solvers/dbnsIncompressibleHBFoam
  01_源码/dbnsIncompressibleHBFoam 和 OpenFOAM-v2112/ICSFoam 未同步修改。

算子
  每相位未知量：运动学压力 p、绝对速度 U、绝对面体积通量 q。
  连续性：V div(q)；动量：空间离散(p,U,q-meshPhi)+原有 HB/参考系项；
  通量闭合：q-Q(p,U,q)。每次试探均从当前 q 组装动量系数，重新计算
  rAU=1/A，而不是使用 flowFrozenResponse。湍流系数在本次流动求解中冻结。
  仅动量加入 V/deltaTau*(U-U0)。最终验收仍检查原物理残差。

  jointFlowStage.H 保存和恢复全部相位状态，计算三组阶段残差。
  jointFlowNewtonKrylov.H 对完整阶段残差作有限差分 Jv，包含 HB 项的作用；
  尚未将 HB 的解析导数单独拆出。使用现有 Newton/右预条件 FGMRES 算法。
  p/U 使用现有块 Jacobi 或压力 Schur 预条件；q 使用无量纲单位块。
  当前预条件器是分块对角原型，尚未实现面通量消元或闭合块 Schur 补。
  因此不能假设其对大型螺旋桨算例一定更快或更稳。

  q 增量以 |Sf|*Uref 归一化；p 以 Uref^2、U 以 Uref 归一化。
  单元残差采用体积加权，面残差采用面积加权，耦合面两侧各记一半面积。
  固定速度/empty 边界的通量由边界条件确定，不增加独立 q 自由度。
  AMI、cyclic、processor 上保留独立的本地通量迹；闭合目标沿用原生相位
  接口操作。这不是一个新的保守 AMI 插值方案。必须进一步验证接口两侧的
  通量守恒、并行一致性和闭合收敛；本原型未通过这些数值检查。
  保留所有连续性行，参考压力以原有 gauge 权重加入对应行；最终分别检查
  完整连续性和参考值，不允许通过丢弃参考单元连续性掩盖全局不守恒。

开启方法（添加/修改 system/fvSolution 的 pseudoTime 条目，勿重复键）
  jointFlowOperator         true;
  flowLineSearch            true;
  stageZeroInitialGuess     true;
  fullStageCorrectors       1;
  stageNewtonSteps          3;
  stageKrylovDirections     30;
  stageKrylovCycles         3;
  stageKrylovRelTol         0.1;
  stageNewtonRelTol         0.1;
  stageFiniteDifferenceEpsilon 1e-7;
  stagePreconditioner       blockJacobi; // 或 pressureSchur
  flowFacePseudoTime        false;
  stageAdvectiveMomentum    false;
  gaugeReplacesContinuity   false;
  pressureFluxCorrectors    0;
  massConsistentMomentum    false;
  jointFlowJacobianAudit    false; // 小算例诊断时设 true

  这些数值仅为起始配置，并非螺旋桨验证参数。保留原算例的 HB 频率、运动、
  边界、湍流、物理容差、尺度和其他必要设置。关闭 jointFlowOperator 恢复
  原求解路径；若此前修改了其他控制参数，也应恢复对应的基线配置。

编译与小算例检查（在具备 OpenFOAM v2112 的 Linux 环境执行）
  bash build.sh
  source env.sh
  python3 Make/checkJointFlow.py --binary "$FOAM_USER_APPBIN/dbnsIncompressibleHBFoam"

  脚本在临时目录生成 32 单元、3 相位、静止网格的层流测试。均匀速度为
  (1,0.2,0.1)，施加小的非均匀压力扰动。它要求新路径确实运行、残差打包
  范数一致、方向导数审计通过、至少一个阶段步和物理步被接受。
  这是三次外迭代的步进检查，允许求解器返回未完全收敛的退出码 2；
  不是螺旋桨、AMI、MPI、动态网格或湍流验证。
  --prepare-only 仅生成输入，不调用 OpenFOAM。本地只运行了此模式。

日志
  IHBJointStage：连续性/动量/闭合的阶段 RMS 和最大值。
  IHBJointNewtonNormCheck：打包向量范数与阶段综合指标一致性。
  IHBJointLinearBlocks：独立重算 Jv 后的三组线性残差。
  IHBJointTaylor：扰动步长、Taylor 余项和相对线性项误差。
    在光滑分支、合适步长范围内应呈二阶余项；有限差分噪声会产生底限。
    边界/限幅器分支切换时不能强行要求二阶，需先定位。
  IHBFlowAccepted/Rejected：最终物理验收结果。

后续验证
  先在 Linux 编译并运行小测试，再检查非正交网格、周期面、AMI 和 MPI。
  保持同一个初场、频率、网格和物理容差，比较旧/新路径的物理残差、
  线搜索接受率、总耗时、推力、扭矩和质量守恒。不得只比较 GMRES 残差。
