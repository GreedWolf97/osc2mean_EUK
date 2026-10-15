# osc2mean_EUK — Python

[返回仓库首页](../README.md) · [MATLAB 版本](../osc2mean_EUK_matlab/README.md)

基于 Eckstein–Ustinov 与 Kaula 方法，将瞬时轨道根数转换为平均轨道根数。

## 运行

使用 Python 3.10 或更新版本，在本目录执行：

```sh
python -m pip install -r requirements.txt
python verifyOsc2Mean.py
```

也可在 PyCharm 中直接运行 `verifyOsc2Mean.py`。脚本按自身位置查找输入数据。依赖为 NumPy 和 Matplotlib，版本下限见 `requirements.txt`。

运行后计算 18 个点，打印半长轴峰峰值、J2 迭代残差及迭代次数，并绘制半长轴对比与 EUK 平均半长轴相对样本均值的偏差；不自动保存结果文件。

无图形界面的环境中可在本目录执行：

```sh
python -c "import matplotlib; matplotlib.use('Agg'); from verifyOsc2Mean import main; means, j2means = main(show=False)"
```

## 文件与输入

| 文件 | 用途 |
| --- | --- |
| osc2mean.py | EU 与 Kaula 算法及辅助函数 |
| verifyOsc2Mean.py | 18 点示例入口及绘图 |
| requirements.txt | Python 依赖 |
| data/checkpoints_18.csv | 18 个点的历元和瞬时根数 |
| data/GGM03C.gfc | GGM03C 重力场系数 |

CSV 列名为 `epoch,a_km,e,i_deg,raan_deg,omega_deg,M_deg`。根数顺序为 `[a, e, i, RAAN, omega, M]`，半长轴单位 km，角度单位 deg，最后一项已经是平近点角 M。

历元使用 ISO 格式，例如 `2027-07-08 08:00:01.054932`，精度为微秒。示例不带时区，按给定时刻计算，不自动平移时区。若直接调用函数传入带时区的 datetime，恒星时函数会先转换到 UTC；程序不完成其他时间尺度或参考坐标系转换。

## 参数与函数接口

主程序顶部可修改截断参数，默认 `LMAX=30`、`QMAX=4`、`GN=4`。提高截断参数会增加计算量。

在本目录可直接调用：

```python
from datetime import datetime
from osc2mean import read_GGM03C_gfc, osc2mean_EUK

C, S = read_GGM03C_gfc('data/GGM03C.gfc', 30)
opts = dict(Lmax=30, qMax=4, gN=4, C=C, S=S)
coe = [6907.5880740762695, 0.0013444418500964154,
       97.498205768482691, 0.013580690115964586,
       70.689808217435058, 289.45553778547554]
epoch = datetime.fromisoformat('2027-07-08 08:00:01.054932')
mean, info = osc2mean_EUK(coe, epoch, opts)
```

`mean` 是包含 6 个根数的一维 NumPy 数组。`info` 包含 `converged`、`iter`、`r_error`（km）、`v_error`（km/s）、`coe_J2mean` 和 `dKaula`。收敛标志仅针对 J2 反演。Python 的诊断字典不包含 MATLAB 版本提供的完整迭代历史。

J2-Mean 是 EU 阶段结果，EUK-Mean 是继续经过 Kaula 改正的结果。默认 18 点示例半长轴峰峰值参考分别为 169.809141 m 和 1.367967 m。峰峰值表示样本最大值减最小值，不是相对外部真值的误差。

模型适用范围、两版输入差异与理论参考见[仓库首页](../README.md)。