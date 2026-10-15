function [coe_mean, info] = osc2mean_EUK(coe_osc, epoch, opts)
% =========================================================================
% OSC2MEAN_EUK
% 瞬时轨道根数 -> 完全平均轨道根数 (EUK + Kaula)
% 输入:
%   coe_osc = [a; e; i; RAAN; omega; M]  瞬时根数 (km, deg)
%   epoch   : datetime
% 输出:
%   coe_mean = [a; e; i; RAAN; omega; M]  完全平均根数
% =========================================================================
%% EUK 瞬转平
% Step 1: J2 平均
[coe_J2mean, info] = osc2mean_EU(coe_osc, epoch, opts);
if ~info.converged
    warning('osc2mean_EUK:NoConvergence', 'J2 inversion did not converge after %d iterations.', info.iter);
end

% Step 2: Kaula 短周期改正
dKaula = KaulaShortPeriod(coe_J2mean, epoch, opts);

% Step 3: 平均 = J2平均 - Kaula改正
coe_mean = coe_J2mean - dKaula;

% 角度归算
coe_mean(4:6) = mod(coe_mean(4:6), 360);
info.coe_J2mean = coe_J2mean;
info.dKaula = dKaula;

end
