clc; clear;
projectDir = fileparts(mfilename('fullpath'));
addpath(projectDir);

opts = struct();

% 重力场参数
opts.Lmax = 30;
opts.qMax = 4;
opts.gN = 4; % 级数保留项，e约为0.001，4项就够

% 读取 GGM03C 重力场系数
gfcFile = fullfile(projectDir, 'GGM03C.gfc');
[C, S] = read_GGM03C_gfc(gfcFile, opts.Lmax);
opts.C = C;
opts.S = S;

% ===== 读取输入文件：按列位置读，不依赖列名 =====
inputFile = fullfile(projectDir, 'orbit_input.xlsx');

% 1. 读取数值部分（第2~7列，包含表头），然后去掉第1行表头
raw_nums = readmatrix(inputFile, 'Range', 'B:G');
coe_osc = raw_nums(2:end, :); % 剔除表头，得到纯数值矩阵 n x 6
n = size(coe_osc, 1);
assert(n > 0, 'Input file contains no data rows.');

% 2. 读取时间部分（第1列），去掉表头
% readtable 默认会把第一行作为列名，所以 inputTable{:,1} 直接就是数据
inputTable = readtable(inputFile, 'Range', 'A:A');
epochStrs = string(inputTable{:, 1});

% 3. 自动识别时间格式
times = datetime(epochStrs);
% ===== 读取输入文件结束 =====

a_osc = coe_osc(:, 1);
coe_means = zeros(n, 6);
coe_J2means = zeros(n, 6);
r_errors = zeros(n, 1);
v_errors = zeros(n, 1);
iterations = zeros(n, 1);

for j = 1:n
    [coe_mean, info] = osc2mean_EUK(coe_osc(j, :).', times(j), opts);
    assert(info.converged, 'J2 inversion did not converge at point %d.', j);
    coe_J2means(j, :) = info.coe_J2mean.';
    r_errors(j) = info.r_error;
    v_errors(j) = info.v_error;
    iterations(j) = info.iter;
    coe_means(j, :) = coe_mean.';
end

% 画图
plot(times, a_osc, 'b-', 'LineWidth', 1.5); hold on;
plot(times, coe_means(:, 1), 'r-', 'LineWidth', 1.5);
plot(times, coe_J2means(:, 1), 'g--', 'LineWidth', 1);
xlabel('Epoch');
ylabel('Semi-major Axis (km)');
legend('Osculating', 'EUK-Mean', 'J2-Mean', 'Location', 'best');
title('Semi-major Axis: Osculating vs EUK-Mean');
grid on;
fprintf('Samples: %d\n', numel(a_osc));
fprintf('Osculating a peak-to-peak: %.6f m\n', 1000 * (max(a_osc) - min(a_osc)));
fprintf('J2-mean a peak-to-peak: %.6f m\n', 1000 * (max(coe_J2means(:, 1)) - min(coe_J2means(:, 1))));
fprintf('EUK-mean a peak-to-peak: %.6f m\n', 1000 * (max(coe_means(:, 1)) - min(coe_means(:, 1))));
fprintf('Max J2 residual: %.6e km, %.6e km/s\n', max(r_errors), max(v_errors));
fprintf('J2 iterations: %d to %d\n', min(iterations), max(iterations));