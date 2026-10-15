function [Gvec, dGvec] = KaulaG_EU(l, p, e, qMax, N, varargin)
% KaulaG_EU 计算Kaula偏心率函数 G_lpq(e) 及其导数
%
% 基于 Kaula/Tisserand 公式:
%   G_lpq(e) = (-1)^|q| (1+beta^2)^l beta^|q| 
%              * sum_{k=0}^\infty P_lpqk Q_lpqk beta^(2k)
%   其中 beta = e / (1 + sqrt(1-e^2))
%
% 输入参数:
%   l     : Kaula degree (非负整数)
%   p     : Kaula index, 0 <= p <= l (整数)
%   e     : 偏心率, 0 <= e < 1
%   qMax  : q 范围 [-qMax, ..., qMax] (非负整数)
%   N     : 级数截断项数 (正整数, 默认50)
%
% 输出:
%   Gvec  : [G_l,p,-qMax ... G_l,p,qMax]
%   dGvec : 对e的数值导数
%
% Compute all Hansen coefficients; period selection belongs to the caller.

%% =========================================================================
% 1. 参数检查
% =========================================================================

if nargin < 5 || isempty(N)
    N = 50;
end

if nargin < 4 || isempty(qMax)
    qMax = 5;
end

if ~isscalar(l) || l < 0 || l ~= floor(l)
    error('l 必须是非负整数。');
end

if ~isscalar(p) || p < 0 || p > l || p ~= floor(p)
    error('p 必须满足 0 <= p <= l，且为整数。');
end

if ~isscalar(e) || e < 0 || e >= 1
    error('偏心率 e 必须满足 0 <= e < 1。');
end

if ~isscalar(qMax) || qMax < 0 || qMax ~= floor(qMax)
    error('qMax 必须是非负整数。');
end

if ~isscalar(N) || N < 1 || N ~= floor(N)
    error('N 必须是正整数。');
end

%% =========================================================================
% 2. 初始化
% =========================================================================

nq = 2*qMax + 1;
qList = -qMax:qMax;

Gvec  = zeros(1, nq);
dGvec = zeros(1, nq);

%% =========================================================================
% 3. e = 0 特殊处理
% =========================================================================

if e == 0
    for iq = 1:nq
        q = qList(iq);
        
        % Circular limit: G_lpq(0) = delta(q,0).
        if q == 0
            Gvec(iq) = 1;
        else
            Gvec(iq) = 0;
        end
        dGvec(iq) = NaN;
    end
    
    % 使用单边差分计算 e=0 处的导数
    h = 1e-7;
    for iq = 1:nq
        q = qList(iq);
        if q == 0
            dGvec(iq) = 0;
        else
            Gp = compute_G_single(l, p, h, q, N);
            dGvec(iq) = (Gp - 0) / h;
        end
    end
    return;
end

%% =========================================================================
% 4. 正常偏心率 e > 0
% =========================================================================

for iq = 1:nq
    q = qList(iq);
    
    % 计算 G_lpq(e)
    Gvec(iq) = compute_G_single(l, p, e, q, N);
    
    % 数值计算 dG/de (中心差分)
    h = max(1e-8, 1e-5 * max(e, 1e-3));
    if ~isempty(varargin)
        h = varargin{1};
        validateattributes(h, {'numeric'}, {'scalar','positive','finite'});
    end
    ePlus  = min(e + h, 1 - 1e-10);
    eMinus = max(e - h, 0);
    
    Gplus  = compute_G_single(l, p, ePlus,  q, N);
    Gminus = compute_G_single(l, p, eMinus, q, N);
    
    if ePlus ~= eMinus
        dGvec(iq) = (Gplus - Gminus) / (ePlus - eMinus);
    else
        dGvec(iq) = 0;
    end
end

end


%% =========================================================================
% 子函数: 计算单个 G_lpq(e)
% =========================================================================

function G = compute_G_single(l, p, e, q, N)
% 计算单个 Hansen 系数 G_lpq(e)

%% -------------------------------------------------------------------------
% 1. 零频项检查
% -------------------------------------------------------------------------

%% -------------------------------------------------------------------------
% 2. e = 0
% -------------------------------------------------------------------------

if e == 0
    G = (q == 0);
    return;
end

%% -------------------------------------------------------------------------
% 3. 计算 beta
% -------------------------------------------------------------------------

sqrtTerm = sqrt(1 - e^2);
beta = e / (1 + sqrtTerm);
beta2 = beta^2;

%% -------------------------------------------------------------------------
% 4. 确定 p' 和 q' (Kaula 定义)
% -------------------------------------------------------------------------

if p <= l/2
    pPrime = p;
    qPrime = q;
else
    pPrime = l - p;
    qPrime = -q;
end

%% -------------------------------------------------------------------------
% 5. 计算 alpha 和 A
% -------------------------------------------------------------------------

alpha = l - 2*pPrime + qPrime;

% 使用恒等式 e/beta = 1 + sqrt(1-e^2) 提高数值稳定性
A = 0.5 * alpha * (1 + sqrtTerm);

%% -------------------------------------------------------------------------
% 6. 级数求和
% -------------------------------------------------------------------------

sumG = 0;

for k = 0:(N-1)
    
    % ---- 确定 h 值 ----
    if qPrime > 0
        hP = k + qPrime;
        hQ = k;
    else
        hP = k;
        hQ = k - qPrime;
    end
    
    % ---- 计算 P_lpqk ----
    Psum = 0;
    for r = 0:hP
        n = hP - r;
        C1 = generalized_binomial(2*pPrime - 2*l, n);
        if C1 == 0
            continue;
        end
        Psum = Psum + C1 * (-1)^r / factorial(r) * A^r;
    end
    
    % ---- 计算 Q_lpqk ----
    Qsum = 0;
    for r = 0:hQ
        n = hQ - r;
        C2 = generalized_binomial(-2*pPrime, n);
        if C2 == 0
            continue;
        end
        Qsum = Qsum + C2 / factorial(r) * A^r;
    end
    
    % ---- 累加 ----
    termK = Psum * Qsum * beta2^k;
    sumG = sumG + termK;
    
    % ---- 提前收敛判断 ----
    if k > 5 && abs(termK) < 1e-15 * max(1, abs(sumG))
        break;
    end
end

%% -------------------------------------------------------------------------
% 7. 最终组合
% -------------------------------------------------------------------------

G = (-1)^abs(q) * (1 + beta2)^l * beta^abs(q) * sumG;

% 极小数清零
if abs(G) < 1e-15
    G = 0;
end

end


%% =========================================================================
% 子函数: 广义二项式系数
% =========================================================================

function C = generalized_binomial(a, n)
% 计算广义二项式系数 C(a,n) = a(a-1)...(a-n+1)/n!
%
% 支持 a 为负整数的情况 (原代码错误地将负数情况置零)
%
% 示例:
%   C(-4,0) = 1
%   C(-4,1) = -4
%   C(-4,2) = 10
%   C(-4,3) = -20
%   C(-4,4) = 35

if n < 0 || n ~= floor(n)
    C = 0;
    return;
end

if n == 0
    C = 1;
    return;
end

if a == 0
    C = 0;
    return;
end

% 正整数 a: C(a,n) = 0 当 n > a
if a > 0 && a == floor(a) && n > a
    C = 0;
    return;
end

% 递推计算: C(a,n) = C(a,n-1) * (a-n+1)/n
C = 1;
for j = 1:n
    C = C * (a - j + 1) / j;
end

end
