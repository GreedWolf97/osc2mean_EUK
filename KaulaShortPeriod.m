function dcoe = KaulaShortPeriod(coe, epoch, opts)
% KaulaShortPeriod: Kaula一阶短周期摄动（所有非J2谐项）
% 输入:  coe = [a; e; i; Omega; omega; M] (km, deg), epoch (datetime)
% 输出: dcoe = [da; de; di; dOmega; domega; dM] (km, deg)
% 注: J2 (l=2, m=0) 已排除，由 Eckstein-Ustinov 理论处理

%% 1. 解析输入参数
coe = coe(:);
if numel(coe) ~= 6
    error('coe 必须包含6个元素。');
end

a      = coe(1);          % 半长轴 (km)
e      = coe(2);          % 偏心率
inc    = coe(3);          % 轨道倾角 (deg)
Omega  = coe(4);          % 升交点赤经 (deg)
omega  = coe(5);          % 近地点幅角 (deg)
M      = coe(6);          % 平近点角 (deg)

%% 2. 选项参数设置
% 物理常数
if ~isfield(opts, 'mu'), opts.mu = 398600.4418; end
if ~isfield(opts, 'Re'), opts.Re = 6378.1363; end
if ~isfield(opts, 'J2'), opts.J2 = 1.08262668e-3; end
mu = opts.mu;
Re = opts.Re;
J2 = opts.J2;
if e <= 0 || e >= 1 || abs(sind(inc)) < 1e-12
    error('KaulaShortPeriod:SingularElements', ...
        'Classical Kaula corrections require 0 < e < 1 and non-equatorial inclination.');
end

% 计算参数
if ~isfield(opts, 'Lmax')
    opts.Lmax = size(opts.C, 1) - 1;
end

if ~isfield(opts, 'qMax')
    opts.qMax = 5;
end

if ~isfield(opts, 'gN')
    opts.gN = 8;
end

if ~isfield(opts, 'gDerivativeStep')
    opts.gDerivativeStep = max(1e-8, 1e-4 * max(e, 1e-3));
end

if ~isfield(opts, 'gravityNormalized')
    opts.gravityNormalized = true;
end

Lmax = opts.Lmax;
qMax = opts.qMax;
gN   = opts.gN;

%% 3. 检查系数矩阵大小
if size(opts.C, 1) < Lmax + 1 || size(opts.C, 2) < Lmax + 1
    error('opts.C 维度小于 Lmax+1。');
end
if size(opts.S, 1) < Lmax + 1 || size(opts.S, 2) < Lmax + 1
    error('opts.S 维度小于 Lmax+1。');
end

%% 4. 系数归一化（如需要）
C = zeros(Lmax + 1, Lmax + 1);
S = zeros(Lmax + 1, Lmax + 1);

for l = 0:Lmax
    for m = 0:l
        Cbar = opts.C(l + 1, m + 1);
        Sbar = opts.S(l + 1, m + 1);
        
        if opts.gravityNormalized
            % 归一化因子
            Nlm = sqrt((2 * l + 1) * (2 - (m == 0)) * factorial(l - m) / factorial(l + m));
            C(l + 1, m + 1) = Cbar * Nlm;
            S(l + 1, m + 1) = Sbar * Nlm;
        else
            C(l + 1, m + 1) = Cbar;
            S(l + 1, m + 1) = Sbar;
        end
    end
end

%% 5. 地球自转（格林尼治恒星时）
theta = gmst_EU(epoch);

%% 6. 平均进动椭圆参数
i_rad   = deg2rad(inc);
n       = sqrt(mu / a^3);           % 平均角速度 (rad/s)
p_orbit = a * (1 - e^2);            % 半通径 (km)

%% 7. J2长期项变化率
Omega_dot = -1.5 * J2 * n * (Re / p_orbit)^2 * cos(i_rad);
omega_dot =  0.75 * J2 * n * (Re / p_orbit)^2 * (4 - 5 * sin(i_rad)^2);
M_dot     = n - 0.75 * J2 * n * (Re / p_orbit)^2 * sqrt(1 - e^2) * (3 * sin(i_rad)^2 - 2);
theta_dot = 7.29211514670698e-5;    % 地球自转角速度 (rad/s)

%% 8. Kaula倾角函数
[Ftable, dFtable] = KaulaFTable_EU(Lmax, inc);

%% 9. Kaula偏心率函数
Gtable  = zeros(Lmax + 1, Lmax + 1, Lmax + 1, 2*qMax + 1);
dGtable = zeros(Lmax + 1, Lmax + 1, Lmax + 1, 2*qMax + 1);

for l = 2:Lmax
    for p = 0:l
        [Gvec, dGvec] = KaulaG_EU(l, p, e, qMax, gN, opts.gDerivativeStep);
        % 将同一 (l,p) 的 G 向量复制到所有 m (0:l)
        for m = 0:l
            Gtable(l+1, m+1, p+1, :) = Gvec;
            dGtable(l+1, m+1, p+1, :) = dGvec;
        end
    end
end

%% 10. 初始化摄动累加器
d_a       = 0;
d_e       = 0;
d_i       = 0;
d_Omega   = 0;
d_omega   = 0;
d_M       = 0;

%% 11. Kaula级数主循环
for l = 2:Lmax
    for m = 0:l
        % 跳过 J2（由 Eckstein-Ustinov 处理）
        if l == 2 && m == 0
            continue;
        end
        
        Clm = C(l + 1, m + 1);
        Slm = S(l + 1, m + 1);
        
        if Clm == 0 && Slm == 0
            continue;
        end
        
        for p = 0:l
            F     = Ftable(l + 1, m + 1, p + 1);
            dF    = dFtable(l + 1, m + 1, p + 1);
            l_2p  = l - 2 * p;
            
            for q = -qMax:qMax
                G  = Gtable(l + 1, m + 1, p + 1, q + qMax + 1);
                dG = dGtable(l + 1, m + 1, p + 1, q + qMax + 1);
                
                % 谐波频率
                l_2pplusq = l - 2 * p + q;
                % Retain short-period and m-daily terms, but preserve
                % zonal long-period and secular variations in mean elements.
                if l_2pplusq == 0 && m == 0
                    continue;
                end
                psi_dot   = l_2p * omega_dot + l_2pplusq * M_dot + m * (Omega_dot - theta_dot);
                
                % 跳过长周期项（非周期项）
                if abs(psi_dot) < 1e-16
                    continue;
                end
                
                % 谐波相位
                psi = l_2p * deg2rad(omega) + l_2pplusq * deg2rad(M) + m * (deg2rad(Omega) - theta);
                
                % S_lmpq 和 S_tilde_lmpq
                if mod(l - m, 2) == 0
                    S_lmpq  =  Clm * cos(psi) + Slm * sin(psi);
                    S_tilde = Clm * sin(psi) - Slm * cos(psi);
                else
                    S_lmpq  = -Slm * cos(psi) + Clm * sin(psi);
                    S_tilde = -Slm * sin(psi) - Clm * cos(psi);
                end
                
                % 摄动计算（角度单位为弧度）
                common_a = mu * Re^l / (n * a^(l + 2));
                common_e = mu * Re^l / (n * a^(l + 3) * e);
                common_i = mu * Re^l / (n * a^(l + 3));
                denom    = sqrt(1 - e^2) * sin(i_rad);
                
                % 半长轴摄动
                Delta_a = common_a * 2 * F * G * l_2pplusq * S_lmpq / psi_dot;
                
                % 偏心率摄动
                % Integrate the Lagrange equation (paper section 2.1).
                Delta_e = common_e * sqrt(1 - e^2) * F * G * (sqrt(1 - e^2) * l_2pplusq - l_2p) * S_lmpq / psi_dot;
                
                % 倾角摄动
                Delta_i = common_i / denom * F * G * (l_2p * cos(i_rad) - m) * S_lmpq / psi_dot;
                
                % 升交点赤经摄动
                Delta_Omega = common_i / denom * dF * G * S_tilde / psi_dot;
                
                % 近地点幅角摄动
                term1 = sqrt(1 - e^2) / e * F * dG;
                term2 = cot(i_rad) / sqrt(1 - e^2) * dF * G;
                Delta_omega = common_i * (term1 - term2) * S_tilde / psi_dot;
                
                % 平近点角摄动
                term3 = -((1 - e^2) / e) * dG + 2 * (l + 1) * G;
                Delta_M = common_i * term3 * F * S_tilde / psi_dot;
                
                % 累加摄动
                d_a     = d_a + Delta_a;
                d_e     = d_e + Delta_e;
                d_i     = d_i + Delta_i;
                d_Omega = d_Omega + Delta_Omega;
                d_omega = d_omega + Delta_omega;
                d_M     = d_M + Delta_M;
            end
        end
    end
end

%% 12. 将角度摄动转换为度
dcoe = [d_a; d_e; rad2deg(d_i); rad2deg(d_Omega); rad2deg(d_omega); rad2deg(d_M)];

end
