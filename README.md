# osc2mean_EUK

瞬时轨道根数转平均轨道根数的 MATLAB 与 Python 实现。先使用 Eckstein–Ustinov（EU）方法迭代求解 J2 平均根数，再减去 Kaula 一阶非 J2 短周期改正，得到 EUK 平均根数。

## 目录结构

```text
.
├── osc2mean_EUK_matlab/       # MATLAB 实现
│   ├── verifyOsc2Mean.m      # 示例入口
│   ├── orbit_input.xlsx      # 18 点输入
│   └── GGM03C.gfc           # 重力场系数
└── osc2mean_EUK-python/      # Python 实现
    ├── osc2mean.py          # 算法及辅助函数
    ├── verifyOsc2Mean.py    # 示例入口
    ├── requirements.txt
    └── data/
        ├── checkpoints_18.csv
        └── GGM03C.gfc
```

两个版本可独立运行，各自保留一份内容相同的重力场文件。详细说明见 [MATLAB README](osc2mean_EUK_matlab/README.md) 和 [Python README](osc2mean_EUK-python/README.md)。

## MATLAB 运行

在 MATLAB 中将当前文件夹切换到仓库根目录，执行：

```matlab
cd osc2mean_EUK_matlab
verifyOsc2Mean
```

MATLAB 版本已在 R2024b 中验证。示例读取本目录数据，绘制半长轴变化并打印收敛信息，结果保留在工作区，不自动保存结果文件。

## Python 运行

使用 Python 3.10 或更新版本，在仓库根目录的终端执行：

```sh
python -m pip install -r osc2mean_EUK-python/requirements.txt
python osc2mean_EUK-python/verifyOsc2Mean.py
```

依赖为 NumPy 和 Matplotlib，也可在 PyCharm 中直接运行 `verifyOsc2Mean.py`。脚本按自身位置读取数据，无需修改绝对路径。运行后打印计算摘要并显示图窗。

## 输入约定

两版均提供 18 个测试点。根数顺序为 `[a, e, i, RAAN, omega, M]`：半长轴单位 km，偏心率无量纲，角度单位 deg。**最后一个根数已经是平近点角 M，无需再次转换。**

| 内容 | MATLAB Excel | Python CSV |
| --- | --- | --- |
| 文件 | `orbit_input.xlsx` | `data/checkpoints_18.csv` |
| 读取规则 | 按 A:G 列位置读取，不依赖列名 | 按 CSV 列名读取，升交点赤经为 `raan_deg` |
| 历元显示示例 | `2027/7/8 8:00:01.055` | `2027-07-08 08:00:01.054932` |
| 时间处理 | Excel 日期单元格，显示到毫秒；代码自动识别日期 | Python datetime 使用微秒精度 |

MATLAB 表格第一行为表头，A 列为历元，B:G 列依次为六个轨道根数；替换数据时保持列顺序和单位，日期须能被 MATLAB 自动识别。Python CSV 则需沿用原有列名和日期格式。Excel 的显示位数不等于底层日期数值的精度。示例历元不带时区；程序不自动完成坐标系或时间尺度转换。比较两版结果时，应统一历元精度、输入数值和参数，并允许浮点运算带来的微小差异。

## 参数与结果

默认截断参数为 `Lmax = 30`、`qMax = 4`、`gN = 4`。MATLAB 参数在示例的 `opts` 中设置，Python 参数在示例顶部的 `LMAX, QMAX, GN` 中设置。

示例输出区分瞬时根数（Osculating）、EU 阶段结果（J2-Mean）与 Kaula 改正后的结果（EUK-Mean）。18 点示例的半长轴峰峰值参考如下：

| 结果 | 半长轴峰峰值 |
| --- | --- |
| Osculating | 约 18762.640122 m |
| J2-Mean | 约 169.809141 m |
| EUK-Mean | 约 1.367967 m |

峰峰值为样本最大值减最小值，用于复现示例，不代表相对外部真值的误差。`info.converged` 仅表示 J2 反演收敛。

## 适用范围与资料

本程序面向低偏心率近地轨道，包含 J2 与一阶重力场短周期改正，不包含大气阻力、日月引力及太阳光压等摄动。Kaula 经典根数表达要求 `0 < e < 1` 且轨道非赤道；其他轨道需独立验证模型适用性和截断参数。

理论参考：*Precise mean orbital elements determination for LEO monitoring and maintenance*。GGM03C 模型说明和参考文献信息保留在 `.gfc` 文件头中。