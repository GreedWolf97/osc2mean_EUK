function theta = gmst_EU(epoch)
% Compute GMST angle in radians

JD = juliandate(epoch);
T = (JD - 2451545.0) / 36525;

GMST_sec = mod(67310.54841 + (876600 * 3600 + 8640184.812866) * T + ...
               0.093104 * T^2 - 6.2e-6 * T^3, 86400);

theta = 2 * pi * GMST_sec / 86400;

end