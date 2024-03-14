function w = amb(epsilon, tau)
%AMB Adaptive Maxwell-Boltzmann by Hitchcox2022
% Implementation of:
%   - T. Hitchcox and J. R. Forbes, "Mind the gap: Norm-aware adaptive
%       robust loss for multivariate least-squares problems," IEEE Robotics
%       and Automation Letters, vol. 7, pp. 7116–7123, 3 Jul. 2022.
%
% Inputs:
% - epsilon: residuals e' * W * e [N x 1]
% - tau: parameter to compute alphaStar
%
% Editor: Kyungmin John Jung
% Date: 2024-03-12
% Lab: DECAR Group
% Institution: McGill University

% find a* using Newton's method with line search
n_e = 3;
aStar = getAStar(epsilon, n_e, tau);

% mode of distribution
epsilon_tilde = aStar * sqrt(n_e-1);

% shift residuals
e_shift = epsilon - epsilon_tilde;  % shifted residuals
idx_xi  = e_shift >= 0;
xi      = e_shift(idx_xi);          % xi is shifted errors >= 0
nu      = tau - epsilon_tilde;      % shifted critical value

% now find alpha* using Newton's method with line search
alphaStar = getAlphaStar(xi, 0, nu);

% assign the weights.  for e_shift < 0, w = 1.
w         = ones(numel(epsilon), 1);
w(idx_xi) = adaptiveWeights(xi, alphaStar);
end

%% functions for adaptive weight from Barron2019
% the adaptive robust weight function (13)
function w = adaptiveWeights(epsilon, alpha)
    switch alpha
        case 2
            w = ones(1, numel(epsilon));
        case 0
            w = 1./((epsilon.^2)./2 + 1);
        case -Inf
            w = exp(-(epsilon.^2)./2);
        otherwise
            w = ((epsilon.^2)./abs(alpha - 2) + 1).^(alpha/2 - 1);
    end
end

