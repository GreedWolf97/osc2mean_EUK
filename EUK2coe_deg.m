function coe = EUK2coe_deg(x)
% =========================================================================
% EUK2COE - EUK 非奇异元素 -> 经典轨道根数 (角度版本)
%
% 输入: x = [a, h, l, i, RAAN, lambda]
%      角度单位为 [deg]
% 输出: coe = [a, e, i, RAAN, omega, M]
%      角度为度
% =========================================================================

a      = x(1);
h      = x(2);
l      = x(3);
inc    = x(4);
RAAN   = x(5);
lambda = x(6);

% 偏心率
e = hypot(h, l);

% 近地点幅角 [deg] (使用 atan2d)
if e > 1e-15
    omega = atan2d(h, l);
else
    omega = 0;
end

% 平近点角 [deg]
% M = mod(lambda - omega, 360);
M = lambda - omega; 

% 输出
coe = zeros(6, 1);
coe(1) = a;
coe(2) = e;
coe(3) = inc;
coe(4) = mod(RAAN, 360);
coe(5) = mod(omega, 360);
coe(6) = M;

end