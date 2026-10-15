function [coe_mean, info] = osc2mean_EU(coe_osc, epoch, opts)
% TODosc2mean_EU: Osculating COE -> J2-mean COE (Eckstein-Ustinov)
% 
% Input:  coe_osc = [a; e; i; RAAN; omega; M]  (km, deg)
% Output: coe_mean = [a; e; i; RAAN; omega; M] (km, deg)
%         info: structure with convergence info
%
% Reference: "Precise mean orbital elements determination for LEO monitoring"
%
% Note: Perturbations are evaluated in EU space; updates use PV residuals.
%       Delta = Osculating - Mean (in EU space)

%% Default options
if nargin < 3, opts = struct(); end
if ~isfield(opts, 'mu'), opts.mu = 398600.4418; end          % Earth gravitational parameter [km^3/s^2]
if ~isfield(opts, 'Re'), opts.Re = 6378.1363; end            % Earth equatorial radius [km]
if ~isfield(opts, 'J2'), opts.J2 = 1.08262668e-3; end        % J2 coefficient
if ~isfield(opts, 'maxIter'), opts.maxIter = 20; end         % Maximum iterations
if ~isfield(opts, 'tol_r'), opts.tol_r = 1e-9; end           % Position tolerance [km]
if ~isfield(opts, 'tol_v'), opts.tol_v = 1e-12; end          % Velocity tolerance [km/s]
if ~isfield(opts, 'verbose'), opts.verbose = false; end      % Print iteration info

%% Initial osculating state (Keplerian -> PV)
coe_osc = coe_osc(:);
if numel(coe_osc) ~= 6
    error('coe_osc must contain 6 elements: [a; e; i; RAAN; omega; M]');
end

[r_osc, v_osc] = coe2rv_EU(coe_osc, opts.mu);

%% Convert osculating Keplerian elements to Equinoctial (EU) elements
% EU elements: x = [a; h; l; i; Omega; lambda]
% where: h = e*sin(omega), l = e*cos(omega), lambda = omega + M
xEU_osc = coe2EUK_deg(coe_osc);

%% Initial guess: mean = osculating (in EU space)
xEU_mean = xEU_osc;

%% History for convergence monitoring
history = struct();
history.iter = zeros(opts.maxIter, 1);
history.a = zeros(opts.maxIter, 1);
history.e = zeros(opts.maxIter, 1);
history.i = zeros(opts.maxIter, 1);
history.Omega = zeros(opts.maxIter, 1);
history.omega = zeros(opts.maxIter, 1);
history.M = zeros(opts.maxIter, 1);
history.r_err = zeros(opts.maxIter, 1);
history.v_err = zeros(opts.maxIter, 1);

%% ============================================================
% Iterative solution: find mean elements such that:
%   x_osc = x_mean + dx(x_mean)
% where dx(x_mean) is the J2 short-period perturbation
% =============================================================
converged = false;
for iter = 1:opts.maxIter
    % Step 1: Compute short-period perturbations in EU space
    % dx = Delta = Osculating - Mean (in EU space)
    dxEU = eu_J2_short_period(xEU_mean, opts);
    
    % Step 2: Reconstruct osculating elements in EU space
    x_oscEU_calc = xEU_mean + dxEU;
    
    % Step 3: Convert reconstructed EU -> Keplerian -> PV
    coe_osc_calc = EUK2coe_deg(x_oscEU_calc);
    [r_calc, v_calc] = coe2rv_EU(coe_osc_calc, opts.mu);
    
    % Step 4: Compute residuals in PV space (true osculating - reconstructed)
    dr = r_osc - r_calc;
    dv = v_osc - v_calc;
    r_err = norm(dr);
    v_err = norm(dv);
    
    % Step 5: Store history (convert EU mean back to Keplerian for plotting)
    coe_mean_cur = EUK2coe_deg(xEU_mean);
    history.iter(iter) = iter;
    history.a(iter) = coe_mean_cur(1);
    history.e(iter) = coe_mean_cur(2);
    history.i(iter) = coe_mean_cur(3);
    history.Omega(iter) = coe_mean_cur(4);
    history.omega(iter) = coe_mean_cur(5);
    history.M(iter) = coe_mean_cur(6);
    history.r_err(iter) = r_err;
    history.v_err(iter) = v_err;
    
    % Step 6: Print progress
    if opts.verbose
        fprintf('Iter %3d: dr = %.12e km, dv = %.12e km/s, a_mean = %.12f km\n', ...
                iter, r_err, v_err, coe_mean_cur(1));
    end
    
    % Step 7: Check convergence
    if r_err <= opts.tol_r && v_err <= opts.tol_v
        converged = true;
        break;
    end
    
    % Do not update an untested state after the final residual evaluation.
    if iter == opts.maxIter
        break;
    end
    % Step 8: Update mean state using PV correction
    % Convert current mean EU -> Keplerian -> PV
    coe_mean_cur = EUK2coe_deg(xEU_mean);
    [r_mean_cur, v_mean_cur] = coe2rv_EU(coe_mean_cur, opts.mu);
    
    % Apply correction in PV space
    r_mean_new = r_mean_cur + dr;
    v_mean_new = v_mean_cur + dv;
    
    % Convert corrected PV back to Keplerian, then to EU for next iteration
    coe_mean_new = rv2coe_EU(r_mean_new, v_mean_new, opts.mu);
    xEU_mean = coe2EUK_deg(coe_mean_new);
end

%% Trim history arrays to actual iteration count
history.iter = history.iter(1:iter);
history.a = history.a(1:iter);
history.e = history.e(1:iter);
history.i = history.i(1:iter);
history.Omega = history.Omega(1:iter);
history.omega = history.omega(1:iter);
history.M = history.M(1:iter);
history.r_err = history.r_err(1:iter);
history.v_err = history.v_err(1:iter);

%% Final output: convert converged mean EU -> Keplerian
coe_mean = EUK2coe_deg(xEU_mean);

%% Build info structure
info = struct();
info.iter = iter;
info.r_error = r_err;
info.v_error = v_err;
info.history = history;
info.converged = converged;
info.x_mean_final = xEU_mean;
info.x_osc = xEU_osc;
info.r_osc = r_osc;
info.v_osc = v_osc;
[r_mean, v_mean] = coe2rv_EU(coe_mean, opts.mu);
info.r_mean = r_mean;
info.v_mean = v_mean;

if ~converged && opts.verbose
    fprintf('\nWARNING: Maximum iterations (%d) reached without convergence.\n', opts.maxIter);
end

end
