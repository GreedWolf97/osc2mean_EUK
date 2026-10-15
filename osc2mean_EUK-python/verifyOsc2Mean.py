"""18 点瞬时根数转平根测试；在 PyCharm 中直接运行本文件。"""
from pathlib import Path
from datetime import datetime
import csv
import numpy as np
import matplotlib.pyplot as plt
from osc2mean import read_GGM03C_gfc, osc2mean_EUK

BASE = Path(__file__).resolve().parent
LMAX, QMAX, GN = 30, 4, 4  # 重力场阶次、q 的截断范围、偏心率级数项数


def main(show=True):
    # CSV 最后一列已经是平近点角 M，单位为度，无需再次转换。
    with (BASE / 'data/checkpoints_18.csv').open(encoding='utf-8-sig', newline='') as f:
        rows = list(csv.DictReader(f))
    epochs = [datetime.fromisoformat(row['epoch']) for row in rows]
    fields = ['a_km', 'e', 'i_deg', 'raan_deg', 'omega_deg', 'M_deg']
    inputs = np.array([[float(row[key]) for key in fields] for row in rows])
    C, S = read_GGM03C_gfc(BASE / 'data/GGM03C.gfc', LMAX)
    opts = dict(Lmax=LMAX, qMax=QMAX, gN=GN, C=C, S=S)

    # 先用 EU 去除 J2 摄动，再用 Kaula 去除其余重力谐项。
    means, j2means, residuals, iterations = [], [], [], []
    for index, (coe, epoch) in enumerate(zip(inputs, epochs), 1):
        mean, info = osc2mean_EUK(coe, epoch, opts)
        if not info['converged']:
            raise RuntimeError(f'第 {index} 个点的 J2 迭代未收敛')
        means.append(mean)
        j2means.append(info['coe_J2mean'])
        residuals.append([info['r_error'], info['v_error']])
        iterations.append(info['iter'])
    means, j2means = np.array(means), np.array(j2means)
    residuals = np.array(residuals)

    # 峰峰值为各点半长轴的最大值减最小值；km 乘 1000 转为 m。
    print(f'Samples: {len(rows)}')
    print(f'Osculating a peak-to-peak: {1000*np.ptp(inputs[:, 0]):.6f} m')
    print(f'J2-mean a peak-to-peak: {1000*np.ptp(j2means[:, 0]):.6f} m')
    print(f'EUK-mean a peak-to-peak: {1000*np.ptp(means[:, 0]):.6f} m')
    print(f'Max J2 residual: {max(residuals[:, 0]):.6e} km, '
          f'{max(residuals[:, 1]):.6e} km/s')
    print(f'J2 iterations: {min(iterations)} to {max(iterations)}')

    # 对比瞬时根数、EU 平根和 EUK 平根，再单独显示米级剩余变化。
    fig, axes = plt.subplots(2, 1, figsize=(10, 8), sharex=True)
    axes[0].plot(epochs, inputs[:, 0], 'b-', label='Osculating')
    axes[0].plot(epochs, j2means[:, 0], 'g--', label='J2-Mean')
    axes[0].plot(epochs, means[:, 0], 'r-', label='EUK-Mean')
    axes[0].set_ylabel('Semi-major axis (km)')
    axes[0].legend()
    axes[1].plot(epochs, 1000*(means[:, 0]-means[:, 0].mean()), 'o-')
    axes[1].set(xlabel='Epoch', ylabel='EUK mean a deviation (m)')
    for ax in axes:
        ax.grid()
    fig.autofmt_xdate()
    fig.tight_layout()
    if show:
        plt.show()
    else:
        plt.close(fig)
    return means, j2means


if __name__ == '__main__':
    main()
