function [F, dF] = KaulaF_EU(l, m, p, i)
% Kaula inclination function and its derivative
% Based on direct formula (Kaula 1966)
%
% Inputs:
%   l, m, p - integers (0 <= m <= l, 0 <= p <= l)
%   i       - inclination angle in radians
% Outputs:
%   F       - inclination function F_lmp(i)
%   dF      - derivative dF_lmp / di

% --- Validate inputs ---
if l < 2 || m > l || p < 0 || p > l
    F = 0; dF = 0;
    return;
end

s = sin(i/2);
c = cos(i/2);

% --- Common factor ---
% factor = (l+m)! / (2^l * p! * (l-p)!)
factor = factorial(l + m) / (2^l * factorial(p) * factorial(l - p));

% --- Summation limits ---
J = min(l - m, 2*l - 2*p);
j0 = max(0, l - 2*p - m);

F_sum = 0;
dF_sum = 0;

for j = j0:J
    % Binomial coefficients
    bin1 = factorial(2*l - 2*p) / (factorial(j) * factorial(2*l - 2*p - j));
    bin2 = factorial(2*p) / (factorial(l - m - j) * factorial(2*p - l + m + j));
    
    % Check if binomial arguments are valid
    if l - m - j < 0 || 2*p - l + m + j < 0
        continue;
    end
    
    % Coefficient
    coeff = (-1)^j * bin1 * bin2;
    
    % Exponents
    exp_s = m - l + 2*p + 2*j;
    exp_c = 3*l - m - 2*p - 2*j;
    
    % Term and derivative term
    term = coeff * s^exp_s * c^exp_c;
    F_sum = F_sum + term;
    
    % Derivative of term w.r.t. i
    % d/di [s^exp_s * c^exp_c] = term * (exp_s/2 * cot(i/2) - exp_c/2 * tan(i/2))
    dTerm_di = term * (exp_s/2 * cot(i/2) - exp_c/2 * tan(i/2));
    dF_sum = dF_sum + dTerm_di;
end

% --- Apply sign factor ---
sign_factor = (-1)^ceil((l - m)/2);

F  = sign_factor * factor * F_sum;
dF = sign_factor * factor * dF_sum;

end