%% functions for computing alphaStar from [11]
% Implementation of:
%   - N. Chebrolu, T. Labe, O. Vysotska, J. Behley, and C. Stachniss,
%       "Adaptive robust kernels for non-linear least squares problems,"
%       IEEE Robotics and Automation Letters, vol. 6, pp. 2240-2247, 2 Apr.
%       2021

% alpha* line search and supporting functions
function alpha_star = getAlphaStar(epsilon, lower_bound, upper_bound)
    % setup
    options.max_iter    = 15;   % iterations of newton's method
    options.min_delta   = 1e-3; % min change for convergence
    options.tau_0       = 0.1;  % initial line search shrink factor
    options.tau         = 0.5;  % shrink by this amount each ls iteration
    options.max_ls_iter = 15;   % number of line search iterations

    % objective function handle
    objfunc = @(alpha) objFuncAlpha(epsilon, alpha, lower_bound, upper_bound);

    % confine alpha* to this range
    alpha_min = -10;
    alpha_max = 2;

    % initialize to a typical value
    alpha_0 = -2;

    % solve using Newton's method with backtracking basic line search
    alpha_star = newtonsMethodBounded(objfunc, alpha_0, alpha_min, alpha_max, options);
end

% evaluate objective function and objective function gradient
function [ Lambda, dLambda_dalpha ] = objFuncAlpha(epsilon, alpha, tau_min, tau_max)
    N       = numel(epsilon);
    z_tilde = integral(@(epsilon) exp(-adaptiveLoss(epsilon, alpha)), tau_min, tau_max); % (11)
    Lambda  = N.*log(z_tilde) + sum(adaptiveLoss(epsilon, alpha));                       % (12b)
    int_n   = integral(@(epsilon) -exp(-adaptiveLoss(epsilon, alpha)).* ...
              drho_dalpha(epsilon, alpha), tau_min, tau_max);
    dLambda_dalpha = (N./z_tilde).*int_n + sum(drho_dalpha(epsilon, alpha));
end

% evaluate the gradient of the generalized loss function, from (12)
function grad_rho = drho_dalpha(epsilon, alpha)
    % deal with singularities.  don't expect alpha == -Inf in practice.
    if abs(alpha) < 1e-10
        % alpha \approx 0
        grad_rho = (drho_dalpha(epsilon, 1e-9) + drho_dalpha(epsilon, -1e-9))/2;
    elseif abs(alpha - 2) < 1e-10
        % alpha \approx 2
        grad_rho = (drho_dalpha(epsilon, 2-1e-9) + drho_dalpha(epsilon, 2+1e-9))/2;
    else
        b        = (epsilon.^2)./abs(alpha - 2) + 1;
        term_1   = (abs(alpha - 2)/alpha).*(b.^(alpha/2)).*((1/2).*log(b) ...
                   - (alpha.*(epsilon.^2).*(alpha - 2))./(2.*b.*abs(alpha - 2)^3));
        term_2   = -(abs(alpha - 2)/(alpha^2)).*(b.^(alpha/2) - 1);
        term_3   = ((alpha - 2).*(b.^(alpha/2) - 1))./(alpha*abs(alpha - 2));
        grad_rho = term_1 + term_2 + term_3;
    end
end

% the adaptive robust loss function (7)
function rho = adaptiveLoss(epsilon, alpha)
    switch alpha
        case 2
            rho = (epsilon.^2)./2;
        case 0
            rho = log((epsilon.^2)./2 + 1);
        case -Inf
            rho = 1 - exp(-(epsilon.^2)./2);
        otherwise
            rho = (abs(alpha - 2)/alpha).*...
                (((epsilon.^2)./abs(alpha - 2) + 1).^(alpha/2) - 1);
    end
end

