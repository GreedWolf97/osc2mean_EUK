
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

% 输入第六项已经是平近点角 M，直接参与计算。
inputFile = fullfile(projectDir, 'orbit_input.csv');
importOpts = detectImportOptions(inputFile);
importOpts = setvartype(importOpts, 'epoch', 'string');
input = readtable(inputFile, importOpts);
times = datetime(string(input.epoch), 'InputFormat', 'yyyy-MM-dd HH:mm:ss.SSSSSSSSS');
coe_osc = input{:, {'a_km', 'e', 'i_deg', 'RAAN_deg', 'omega_deg', 'M_deg'}};
n = height(input);
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