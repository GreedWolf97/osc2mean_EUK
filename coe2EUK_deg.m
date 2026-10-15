function x = coe2EUK_deg(coe)
% =========================================================================
% COE2EUK_DEG - 经典轨道根数 -> EUK 非奇异元素 (角度版本)
%
% 输入: coe = [a, e, i, RAAN, omega, M]
%      角度单位为 [deg]
% 输出: x = [a, h, l, i, RAAN, lambda]
%      h = e * sin(omega)
%      l = e * cos(omega)
%      lambda = omega + M (角度为度)
% =========================================================================

a      = coe(1);
e      = coe(2);
i      = coe(3);
RAAN   = coe(4);
omega  = coe(5);
M      = coe(6);

x = zeros(6, 1);
x(1) = a;
x(2) = e * sind(omega);
x(3) = e * cosd(omega);
x(4) = i;
x(5) = RAAN;
% x(6) = mod(omega + M, 360);
x(6) = omega + M; 

end