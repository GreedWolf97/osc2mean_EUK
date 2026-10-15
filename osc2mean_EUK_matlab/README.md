# osc2mean_EUK — MATLAB

[返回仓库首页](../README.md) · [Python 版本](../osc2mean_EUK-python/README.md)

基于 MATLAB 的瞬时轨道根数转平均轨道根数程序。先使用 Eckstein–Ustinov 方法迭代求解 J2 平均根数，再减去 Kaula 一阶非 J2 短周期改正。

## 快速运行

使用 MATLAB 打开本目录 `osc2mean_EUK_matlab`，运行：

```matlab
verifyOsc2Mean
```

已在 MATLAB R2024b 上验证。示例读取 `orbit_input.csv` 和 `GGM03C.gfc`，绘制瞬时、J2 平均及 EUK 平均半长轴，并打印样本数、半长轴峰峰值、J2 迭代残差和迭代次数。结果保留在工作区，不自动写入结果文件。

## 输入数据

`orbit_input.csv` 包含 18 个测试点。最后一项为平近点角 M，读取后无需再次转换。历元按文本读取并按 9 位小数秒格式解析，避免自动日期识别造成精度损失。

| 列 | 含义 | 单位 |
| --- | --- | --- |
| epoch | 历元，格式 `yyyy-MM-dd HH:mm:ss.SSSSSSSSS` | 日期时间 |
| a_km | 半长轴 | km |
| e | 偏心率 | 无量纲 |
| i_deg | 轨道倾角 | deg |
| RAAN_deg | 升交点赤经 | deg |
| omega_deg | 近地点幅角 | deg |
| M_deg | 平近点角 | deg |

替换数据时请沿用列名、单位、角度定义及一致的参考坐标系和时间约定。程序不会转换参考坐标系或时间尺度。

## 函数调用

```matlab
opts = struct('Lmax', 30, 'qMax', 4, 'gN', 4);
[opts.C, opts.S] = read_GGM03C_gfc('GGM03C.gfc', opts.Lmax);
coe_osc = [6907.58807407627; 0.00134444185009642; ...
           97.4982057684827; 0.0135806901159646; ...
           70.6898082174351; 289.455537785476];
epoch = datetime('2027-07-08 08:00:01.054932961', ...
                 'InputFormat', 'yyyy-MM-dd HH:mm:ss.SSSSSSSSS');
[coe_mean, info] = osc2mean_EUK(coe_osc, epoch, opts);
```

输入和输出均为列向量 `[a; e; i; RAAN; omega; M]`，半长轴单位 km，角度单位 deg。`info` 包含 J2 迭代收敛标志、位置和速度残差、迭代历史、J2 平均根数与 Kaula 改正量。`info.converged` 只表示 J2 反演收敛。

## 文件组成

| 文件 | 用途 |
| --- | --- |
| verifyOsc2Mean.m | 示例入口、绘图及收敛检查 |
| orbit_input.csv | 18 点瞬时轨道根数，末项为 M |
| osc2mean_EUK.m | J2 平均根数与 Kaula 改正的组合入口 |
| osc2mean_EU.m | J2 平均根数迭代反演 |
| eu_J2_short_period.m | Eckstein–Ustinov J2 短周期摄动 |
| coe2EUK_deg.m、EUK2coe_deg.m | 经典根数与 EU 根数互转 |
| coe2rv_EU.m、rv2coe_EU.m | 轨道根数与位置速度互转 |
| KaulaShortPeriod.m | 非 J2 的 Kaula 短周期改正 |
| KaulaFTable_EU.m、KaulaF_EU.m | 倾角函数及其导数 |
| KaulaG_EU.m | 偏心率函数及其导数 |
| gmst_EU.m | 格林尼治恒星时角 |
| read_GGM03C_gfc.m | 重力场系数读取 |
| GGM03C.gfc | GGM03C 重力场，保留原文件及模型说明 |

示例参数为 `Lmax = 30`、`qMax = 4`、`gN = 4`。重力场原文件保留，可在其支持范围内调整阶数；提高截断参数会增加计算量。

## 适用范围

本程序使用 J2 与一阶重力场短周期模型，不包含大气阻力、日月引力或太阳光压等摄动。Kaula 经典根数表达要求 `0 < e < 1` 且轨道非赤道；圆轨道及接近奇异条件的轨道需另行处理。示例为低偏心率近地轨道，其他轨道应独立验证截断参数与精度。

GGM03C 的模型说明和参考文献信息见 `GGM03C.gfc` 文件头。程序相关理论参考：*Precise mean orbital elements determination for LEO monitoring and maintenance*。
## 验证记录

在 MATLAB R2024b 中完成全部 18 点运行：全部 J2 反演均以 5 次迭代收敛。示例 EUK 平均半长轴峰峰值约 1.367967 m。该数值描述示例样本的剩余变化，不代表相对外部真值的误差。

