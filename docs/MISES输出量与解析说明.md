# MISES 输出量与解析说明（草稿）

本文区分三层：原始 MISES polar 输出、改进 polar 可执行程序增加的输出、Python 数据库派生量。改进 polar 程序的完整计算源码尚未作为本文依据提供。因此，只将可由输出与后处理规则核实的关系写成确定公式；无法确认的内部算法明确标出，不以通用公式代替具体实现。

## 1. 原始 polar 量

`Sinl = tan(β1)`、`Sout = tan(β2)`，其中 β 为入口、出口流角，Python 用 `atan` 转为 `beta1_deg`、`beta2_deg`，转角为 `β1−β2`。`Minl`、`Mout` 是入口/出口 Mach。`omega` 为 polar 的总损失系数，`omega_v`（原头 `omega_V`）为粘性损失分量；二者不可互换。`DF` 是 polar 给出的扩散因子。解析程序直接读取这些字段，并用反正切转换角度；后处理**不重算** `omega`、`omega_v` 或 `DF`；其精确压力基准、归一化公式须由相应 MISES/改进 polar 源码核实，不能据通用教材公式冒充本库的精确定义。

## 2. 改进 polar 新增量

基础版 polar 只有前 18 列；改进可执行程序在后面增加 `Hte_*`、`Wmax`、`W1`、`W2`、`Deq`、尾缘边界层和尾迹等列。数据库读取程序根据实际表头识别新增字段。这些量不是旧版原始 MISES polar 的直接输出；数据库不能声称所有版本都具备。

| 字段 | 当前可确认定义 | 测量位置 / 注意 |
|---|---|---|
| `W1`, `W2` | 改进 polar 的入口、出口相对速度；`Wmax_over_W1 = Wmax/W1` | 速度采用同一程序归一化，绝对 SI 速度尺度需源码核查 |
| `Wmax_top/bottom` | 改进 polar 分别给出的上/下表面最大相对速度 | 取极值的索引范围和是否覆盖尾迹尚未核实 |
| `W_te_top/bottom` | 上/下表面尾缘位置速度 | 尾缘表面站位 `ITEB` |
| `Wmax` | 改进 polar 的单个最大速度字段 | 与物理吸力面选择可能不一致，见下文 |
| `Deq` | 改进 polar 输出的等效扩散比，数值关系约为 `Wmax/W2` | 原始取值使用的 `Wmax` 选择需核查；不可与派生的 `Deq_physical` 混用 |
| `theta_te_top/bottom` | 上/下表面尾缘动量厚度 `θ = ∫(ρu/ρ_eu_e)(1−u/u_e)dn` | 原生值按 MISES 参考长度 `Lref` 无量纲化，站位 `ITEB` |
| `delta_star_te_top/bottom` | 上/下表面尾缘位移厚度 `δ* = ∫[1−ρu/(ρ_eu_e)]dn` | 同上 |
| `H_te_top/bottom` | 各侧形状因子 `δ*/θ` | 同上 |
| `theta_wake_exit`, `delta_star_wake_exit`, `H_wake_exit` | 尾迹出口动量厚度、位移厚度及其比值 `H=δ*/θ` | 当前程序标记站位为 `ITEB+NWAK−2`，不是两侧尾缘值的简单和 |

上表积分式是边界层量的物理定义；改进可执行程序的离散积分、边界位置及速度归一化需要其 Fortran 源码确认。后处理对厚度的单位换算为 `theta_m = theta_raw × Lref`、`theta/c = theta_raw × Lref/chord_m`，位移厚度同理；合并尾缘形状因子为 `(δ*_top+δ*_bottom)/(θ_top+θ_bottom)`。

### Wmax 的未解决问题

后处理 用 `D_velocity_term = 1−W2/W1`、`D_loading_signed = DF−D_velocity_term` 判断反载；正载选 `Wmax_top`，反载选 `Wmax_bottom`，并算 `Deq_physical = Wmax_physical_suction/W2`。另存两侧最大值的 envelope 供审计。代码明确标记 `surface_selected_location_not_audited`，表示**虽然选了表面，尚未核查最大速度发生的位置**，不能视为 Wmax 问题已解决。应保留 `Wmax_raw_polar`、两侧最大值和质量状态，待核查后再用于物理关联式。

## 3. 几何弦长、TSRAT 和哨兵值

示例中的 `chord` 是 MISES 的 `m′–θ` 计算平面几何弦长，单位为该平面的坐标单位，并非米；所述设计的米制弦长由 `baseline_b2b_chord_m × (design_chord / baseline_mises_plane_chord)` 求得。`b2b_chord_m_per_mises_plane_unit` 是基准物理弦长除以基准平面弦长，`input_length_unit_to_m` 是原始输入长度单位换算到米的系数。

本文所述转子热力学流程计算 `tsrat = T_ref / Tt,rel`，默认 `T_ref=110 K`，`Tt,rel` 是转子相对坐标系入口总温；然后写入 MISES IDAT 第四条记录的第 7 个浮点数并回读验证。这一实现关系可以确认；IDAT 参数在 MISES 求解器中的完整理论定义仍须对照源手册/源码。

`separation_fraction_top/bottom` 对应改进 polar 表头 `Ssep_top/bot`。当前数据可见 `-1`；解析代码仅按数值读取，没有把它解释为真实分离位置，也没有定义该哨兵的完整含义。规范应把 `-1` 当作**非物理/未提供值**处理，不纳入统计；正值的确切弧长归一化和其他潜在哨兵需改进 polar 源码确认。`Cf_min_top/bottom` 的负值可能是物理摩擦系数，不应一概当哨兵。后处理的缺失扩展字段使用 `NaN`，不是 `-1`。`wake_exit_data_status=not_available_without_rerun` 表示已有结果缺少尾迹出口数据。

## 4. 其余 polar 输出量

| 字段 | 含义与约定 |
|---|---|
| `Pinl_Po1`, `Pout_Po1` | 入口、出口静压相对入口总压 `Po1` 的比值；其参考状态要与损失系数定义区分。 |
| `Re_million`, `Tu_percent` | polar 报告的 Reynolds 数（单位为百万）和入口湍流强度百分数。 |
| `Xtr_top/bottom` | 上、下表面转捩位置；位置坐标和归一化尺度需由求解器版本说明。 |
| `rVt1`, `delta_rVt` | 旋转流道的角动量相关输出；旋转半径、速度归一化和符号需由所用流道模型说明。 |
| `cl`, `Phi`, `Psi` | 叶型升力系数、流量系数、负荷系数；速度及长度参考尺度需由所用 polar 版本说明。 |
| `separation_fraction_top/bottom`, `Cf_min_top/bottom` | 表面分离位置相关值和最小摩擦系数。`-1` 的分离值不可作为真实分离位置；负摩擦系数可能是真实结果。 |
| `omega_passage` | 扩展程序报告的叶栅通道损失量；不能仅按名称与总损失相加，其差异和归一化方式需由程序定义。 |
| `W_wake_exit` | 尾迹出口速度，使用与其他 `W` 字段相同的程序速度单位。 |
| `density_inlet/outlet`, `density_te_suction/passage` | 入口、出口、尾缘吸力面与通道侧密度；若求解器采用无量纲密度，须给出参考密度。 |
| `passage_width`, `trailing_edge_thickness` | 扩展程序内部使用的通道宽度和尾缘厚度；几何坐标及单位基准需要记录。 |
| `two_theta_over_passage_width`, `te_thickness_over_passage_width` | `2θ/passage_width` 与尾缘厚度/通道宽度；`θ` 的选取站位需记录。 |
| `blockage_over_remaining_open_width`, `blockage_fraction` | 相对剩余开放宽度和相对总宽度的阻塞指标；两者分母不同。 |
| `velocity_te_suction_over_outlet`, `velocity_te_passage_over_outlet` | 尾缘吸力侧或通道侧速度与出口速度之比；`*_minus_one` 为相应比值减 1。 |
| `te_mean_pressure_coefficient_proxy` | 尾缘平均压力系数代理量；它是代理指标，不能不经定义直接当作实测压力系数。 |

这些扩展字段的精确离散公式、参考状态或站位尚未逐项从计算源码核对。在发布用于建模的数值表时，应同时提供程序版本、字段可用性以及速度、密度和长度的归一化约定。
