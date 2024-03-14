%% a* line search and supporting functions from [12]
% Implementation of:
%   - T. Hitchcox and J. R. Forbes, "Mind the gap: Norm-aware adaptive
%       robust loss for multivariate least-squares problems," IEEE Robotics
%       and Automation Letters, vol. 7, pp. 7116–7123, 3 Jul. 2022.

function a_star = getAStar(epsilon, n_e, tau)
    % make normalized histogram q_k(\epsilon) from the residuals.  use
    % \epsilon <= \tau when finding optimal a_star to avoid influence of
    % very large outliers.
    n_bins             = 50;
    [ q_k, bin_edges ] = histcounts(epsilon(epsilon <= tau), n_bins, 'normalization', 'pdf');
    bin_centers        = bin_edges(1:end-1) + (bin_edges(2) - bin_edges(1))/2;
    
    % setup
    options.max_iter    = 15;   % iterations of newton's method
    options.min_delta   = 1e-3; % min change for convergence
    options.tau_0       = 0.5;  % initial line search shrink factor
    options.tau         = 0.5;  % shrink by this amount each ls iteration
    options.max_ls_iter = 15;   % number of line search iterations
    
    % objective function handle
    objfunc = @(a) objFuncA(bin_centers, q_k, a, n_e);
    
    % confine a* to this range, expect \tilde{\epsilon} to be <= tau/2
    a_min = 0.01;
    a_max = (tau/2)/sqrt(n_e-1);
    
    % initialize to expected value
    a_0 = 1;   
    
    % solve using Newton's method with backtracking basic line search
    a_star = newtonsMethodBounded(objfunc, a_0, a_min, a_max, options);
end

% evaluate objective function and objective function gradient
function [ L, dL_da ] = objFuncA(epsilon, q_k, a, n_e)
    p     = p_mb(epsilon, a, n_e);
    L     = sum((q_k.*(p - q_k)).^2);
    dL_da = 2*sum(q_k.^2.*(p - q_k).*dp_da(epsilon, a, n_e));
end

% Maxwell-Boltzmann speed dist., with shape parameter 'a' and n_e DoF.
% from Laurendeau2005 sec. 15.2, (15) in Hitchcox2022.
function p = p_mb(epsilon, a, n_e)
    num   = epsilon.^(n_e - 1).*exp(-(epsilon.^2)./(2*(a^2)));
    denom = (a^n_e)*2^(n_e/2 - 1)*gamma(n_e/2);
    p = num./denom;
end

% evaluate gradient of the Maxwell-Boltzmann speed function, from
function grad_p = dp_da(epsilon, a, n_e)
    num    = -2^(1 - n_e/2).*(n_e*(a^2) - epsilon.^2).*epsilon.^(n_e-1).* ...
             exp(-(epsilon.^2)./(2*(a^2)));
    denom  = (a^(n_e + 3))*gamma(n_e/2);
    grad_p = num./denom;
end