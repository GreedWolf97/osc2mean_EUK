# osc2mean_EUK

**基于 Eckstein–Ustinov（EU）与 Kaula 理论的瞬时轨道根数转平均轨道根数工具，提供 MATLAB 和 Python 两套实现。**

本项目面向低偏心率近地轨道（LEO），从给定历元的六个瞬时经典轨道根数出发，先迭代消除 J2 短周期摄动，再结合 GGM03C 地球重力场计算非 J2 谐项改正，输出六个 EUK 平均根数。可用于分析轨道参数中的周期性变化，为轨道监测和维持研究提供平均根数输入。

## 实现功能

- 支持单历元瞬时根数转换，以及示例脚本中的多历元逐点处理。
- 输出 EU 阶段的 J2 平均根数、Kaula 改正量与最终 EUK 平均根数。
- 通过重构瞬时位置、速度的残差判断 J2 反演收敛，返回迭代次数及残差；MATLAB 版另保留完整迭代历史。
- 包含经典根数与 EU 元素互转、根数与位置速度互转、恒星时计算、重力场读取，以及 Kaula 倾角函数和偏心率函数。
- 附带 18 点输入数据，可绘制半长轴对比并输出峰峰值、迭代残差与迭代次数。

## 示例结果

示例包含 2027-07-08 08:00:01 至 08:44:52 的 18 个历元，轨道倾角约为 97.5°。采用 `Lmax=30`、`qMax=4`、`gN=4`，Python 版本实测结果如下：

| 处理阶段 | 半长轴峰峰值 |
| --- | ---: |
| 瞬时根数（Osculating） | 18762.640122 m |
| EU 阶段 J2 平均根数（J2-Mean） | 169.809141 m |
| EUK 平均根数（EUK-Mean） | 1.367967 m |

![18 点示例的半长轴对比及 EUK 平均半长轴剩余变化](docs/images/semi-major-axis.png)

第一幅图对比三个处理阶段，第二幅图展示 EUK 平均半长轴相对这 18 个样本均值的偏差。峰峰值为样本最大值减最小值，说明该示例中的周期性波动得到明显抑制。**1.37 m 是样本内的剩余波动，并非相对外部真值的误差，也不代表所有轨道均能达到相同效果。**

本次复现环境为 Python 3.13.11、NumPy 2.5.5、Matplotlib 3.11.2。全部 18 点的 J2 反演均在第 5 次迭代收敛，最大位置残差为 `7.020372e-11 km`，最大速度残差为 `1.845240e-13 km/s`。

MATLAB 版已有说明记录了 R2024b 的验证结果；本次未重新完成 MATLAB 运行。两版输入分别使用 Excel 日期和 CSV 微秒时间文本，历元解析精度及浮点运算可能造成结果差异，应以实际运行输出为准。

## 方法

### 1. EU：迭代反演 J2 平均根数

将输入经典根数转换为 EU 元素：

$$
\mathbf{x}=[a,h,l,i,\Omega,\lambda],\qquad
h=e\sin\omega,\quad l=e\cos\omega,\quad\lambda=\omega+M.
$$

以瞬时根数作为平均根数初值，求解：

$$
\mathbf{x}_{\mathrm{osc}}=
\mathbf{x}_{\mathrm{J2mean}}+
\Delta\mathbf{x}_{\mathrm{J2}}(\mathbf{x}_{\mathrm{J2mean}}).
$$

每次迭代由当前平均根数计算 EU 短周期改正，重构瞬时根数并转换为位置速度。将输入瞬时状态与重构状态之差加到当前平均位置速度上，再转回 EU 元素继续迭代。代码在半长轴改正中额外保留纯 J2 二阶项。

默认最多迭代 20 次，位置残差容限为 `1e-9 km`，速度残差容限为 `1e-12 km/s`，两项同时满足才判定收敛。这些容限约束的是模型内部的重构残差。

### 2. Kaula：非 J2 重力谐项改正

读取 GGM03C 的完全归一化球谐系数，转换为计算所需系数，并在 J2 平均根数处计算一阶 Kaula 改正：

- 通过倾角函数及其导数、偏心率函数及其数值导数构造各谐项。
- 使用历元对应的格林尼治恒星时计算地球自转相位，并用 J2 长期进动率构造谐项频率。
- 跳过已由 EU 处理的 `l=2, m=0` 项，避免重复改正。
- 累加截断范围内的短周期及 m-daily 项；保留 `m=0` 且 `l-2p+q=0` 的带谐长期、长周期变化；跳过绝对频率小于 `1e-16 rad/s` 的项。

最终按经典根数逐项相减：

$$
\mathbf{c}_{\mathrm{EUKmean}}=
\mathbf{c}_{\mathrm{J2mean}}-
\Delta\mathbf{c}_{\mathrm{Kaula}}.
$$

输出的升交点赤经、近地点幅角和平近点角归算至 `[0, 360)` 度。这里的平均根数由上述摄动模型和截断方式定义，并非对输入时间序列做滑动平均。

## 部署与运行

本仓库以本地脚本和函数库方式运行。MATLAB 与 Python 版本可独立使用，无需同时安装。

### 获取代码

```sh
git clone https://github.com/GreedWolf97/osc2mean_EUK_matlab.git
cd osc2mean_EUK_matlab
```

以下终端命令均从仓库根目录执行。仓库名与内部 MATLAB 子目录名相同，运行时注意当前所在位置。

### Python

环境：Python 3.10 或更新版本；依赖为 `numpy>=1.24`、`matplotlib>=3.7`。

建议创建独立虚拟环境：

```sh
python -m venv .venv
```

Windows PowerShell：

```powershell
.\.venv\Scripts\python.exe -m pip install -r osc2mean_EUK-python/requirements.txt
.\.venv\Scripts\python.exe osc2mean_EUK-python/verifyOsc2Mean.py
```

Linux / macOS：

```sh
.venv/bin/python -m pip install -r osc2mean_EUK-python/requirements.txt
.venv/bin/python osc2mean_EUK-python/verifyOsc2Mean.py
```

运行后终端打印统计信息并显示两幅子图。也可在 PyCharm 中选择安装好依赖的解释器，直接运行 `verifyOsc2Mean.py`。脚本按自身路径寻找数据。

无图形界面的环境可先进入 `osc2mean_EUK-python`，使用对应虚拟环境的 Python 执行：

```sh
python -c "import matplotlib; matplotlib.use('Agg'); from verifyOsc2Mean import main; means, j2means = main(show=False)"
```

`main` 返回两个 `N×6` 数组：EUK 平均根数与 J2 平均根数。示例默认不写出结果文件。

### MATLAB

已有项目说明记录的验证版本为 MATLAB R2024b。在 MATLAB 中将当前文件夹设为仓库根目录，运行：

```matlab
cd osc2mean_EUK_matlab
verifyOsc2Mean
```

脚本读取本目录的 `orbit_input.xlsx` 和 `GGM03C.gfc`，绘制半长轴对比并打印统计信息。输出保留在工作区：`coe_means` 为 EUK 平均根数，`coe_J2means` 为 J2 平均根数，另有 `r_errors`、`v_errors` 和 `iterations`。示例不自动保存结果文件。

## 输入与函数接口

输入、输出的经典根数顺序统一为：

```text
[a, e, i, RAAN, omega, M]
```

| 参数 | 含义 | 单位 |
| --- | --- | --- |
| `a` | 半长轴 | km |
| `e` | 偏心率 | 无量纲 |
| `i` | 轨道倾角 | deg |
| `RAAN` / `Omega` | 升交点赤经 | deg |
| `omega` | 近地点幅角 | deg |
| `M` | 平近点角 | deg |

**最后一项必须是平近点角 M，不能直接填入真近点角。**

- Python：替换 `data/checkpoints_18.csv`，保留表头 `epoch,a_km,e,i_deg,raan_deg,omega_deg,M_deg`；历元采用如 `2027-07-08 08:00:01.054932` 的 ISO 时间文本。
- MATLAB：替换 `orbit_input.xlsx`，第一行为表头，A 列为历元，B:G 列按上述顺序排列；代码按列位置读取，时间格式须能由 MATLAB 自动识别。

程序不进行参考坐标系转换或完整的时间尺度转换。示例历元无时区；Python 恒星时函数对带时区的 `datetime` 会先转换到 UTC。自定义数据须保持参考坐标系、时间约定及单位一致。

### Python 单点调用

在 `osc2mean_EUK-python` 目录下执行：

```python
from datetime import datetime
from osc2mean import read_GGM03C_gfc, osc2mean_EUK

C, S = read_GGM03C_gfc('data/GGM03C.gfc', 30)
opts = dict(Lmax=30, qMax=4, gN=4, C=C, S=S)
coe = [6907.5880740762695, 0.0013444418500964154,
       97.49820576848269, 0.013580690115964586,
       70.68980821743506, 289.45553778547554]
epoch = datetime.fromisoformat('2027-07-08 08:00:01.054932')

mean, info = osc2mean_EUK(coe, epoch, opts)
print(mean)
print(info['converged'], info['iter'])
```

MATLAB 对应接口为（在 MATLAB 子目录中调用）：

```matlab
opts = struct('Lmax', 30, 'qMax', 4, 'gN', 4);
[opts.C, opts.S] = read_GGM03C_gfc('GGM03C.gfc', opts.Lmax);
coe_osc = [6907.58807407627; 0.00134444185009642; ...
           97.4982057684827; 0.0135806901159646; ...
           70.6898082174351; 289.455537785476];
epoch = datetime(2027, 7, 8, 8, 0, 1.055);
[coe_mean, info] = osc2mean_EUK(coe_osc, epoch, opts);
```

MATLAB 的输入、输出为 6 元素列向量；Python 的 `mean` 为一维 NumPy 数组。`info` 包含收敛标志、迭代次数、位置速度残差、`coe_J2mean` 与 `dKaula`。**`info.converged` 仅表示 EU 阶段 J2 反演收敛，不表示最终 Kaula 改正的精度已经验证。**

## 参数设置

| 参数 | 示例值 | 含义 |
| --- | ---: | --- |
| `Lmax` | 30 | 重力场最大次数，同时包含各次内的所有阶 |
| `qMax` | 4 | Kaula 求和范围 `q ∈ [-qMax, qMax]` |
| `gN` | 4 | 偏心率函数级数的保留项数 |
| `maxIter` | 20 | EU 反演最大迭代次数 |
| `tol_r` | `1e-9` | 位置残差容限，km |
| `tol_v` | `1e-12` | 速度残差容限，km/s |

前 3 项由示例显式设置，其余为 EU 函数默认值。修改 MATLAB 示例中的 `opts`，或 Python 示例顶部的 `LMAX, QMAX, GN` 即可调整截断参数。两版均允许通过 `opts` 自定义 `mu`、`Re`、`J2`。提供的 GGM03C 原文件最高为 360 次，但较高截断值会增加计算量，须结合轨道条件验证数值稳定性和截断误差。

## 目录结构

```text
.
├── README.md
├── docs/images/semi-major-axis.png   # Python 示例生成的结果图
├── osc2mean_EUK_matlab/
│   ├── verifyOsc2Mean.m             # 18 点验证与绘图入口
│   ├── osc2mean_EUK.m               # EU + Kaula 组合入口
│   ├── osc2mean_EU.m                # J2 平均根数迭代反演
│   ├── eu_J2_short_period.m         # EU 短周期摄动
│   ├── KaulaShortPeriod.m           # 非 J2 谐项改正
│   ├── KaulaF*.m / KaulaG_EU.m      # 倾角与偏心率函数
│   ├── orbit_input.xlsx             # MATLAB 输入数据
│   └── GGM03C.gfc                   # 重力场系数
└── osc2mean_EUK-python/
    ├── osc2mean.py                  # 算法与辅助函数
    ├── verifyOsc2Mean.py            # 18 点验证与绘图入口
    ├── requirements.txt
    └── data/
        ├── checkpoints_18.csv
        └── GGM03C.gfc
```

更多版本细节见 [MATLAB 说明](osc2mean_EUK_matlab/README.md) 和 [Python 说明](osc2mean_EUK-python/README.md)。

## 适用范围与参考资料

当前实现主要针对低偏心率近地轨道，使用 EU J2 模型及一阶非 J2 重力谐项改正，未包含大气阻力、日月引力及太阳光压等摄动，也未提供轨道预报或机动控制功能。Kaula 经典根数表达要求 `0 < e < 1` 且轨道非赤道；接近圆轨道、赤道轨道或谐项共振时，需额外验证数值稳定性及模型适用性。

仓库附有理论资料：

- [Precise mean orbital elements determination for LEO monitoring and maintenance](Precise%20mean%20orbital%20elements%20determination%20for%20LEO%20monitoring%20and%20maintenance.pdf)
- [Theory of Satellite Geodesy (1966)](Theory%20of%20Satellite%20Geodesy%281966%29.pdf)
- GGM03C 的模型信息、归一化方式及参考文献见所附 `.gfc` 文件头。
