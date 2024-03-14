function w = gnc_amb(residuals, params)
%GNC_AMB GNC weights update for the AMB cost function
%   The proposed GNC with Hitchcox's Adaptive Maxwell-Boltzmann (GNC-AMB)
%
% Implementaion of
%   - K. Jung, T. Hitchcox, and J. R. Forbes, "An Adaptive Graduated
%       Nonconvexity Loss Function for Robust Nonlinear Least-Squares
%       Solutions," https://arxiv.org/abs/2305.06869
%
% Editor: Kyungmin John Jung
% Date: 2024-03-12
% Lab: DECAR Group
% Institution: McGill University

residual_shift = residuals - params.epsilon_tilde;
idx_xi = residual_shift >= 0;
xi = residual_shift(idx_xi);

w = ones(length(residuals), 1);
w_tmp = ones(length(xi), 1);
mu = params.mu;
alpha = computeAlpha(mu ,params.alphaStar);
for k = 1:length(xi)
    if alpha <= 2
        switch alpha
            case 2
               w_tmp(k) = 1;
            case 0
               w_tmp(k) = 1 / (0.5 * xi(k)^2 + 1);
            case -inf
               w_tmp(k) = exp(-0.5 * xi(k)^2);
            otherwise
               w_tmp(k) = (xi(k)^2 / abs(alpha - 2) + 1)^(alpha/2 - 1);
        end
        assert(w_tmp(k) >= 0 && w_tmp(k) <= 1, ...
            'weight %g is wrong!', w_tmp(k));
    else
        error("Invalid mu value!");
    end
end
w(idx_xi) = w_tmp;

end

function alpha = computeAlpha(mu, alpha_star)
    % alpha = (alpha_star + 2 * mu - 2) / mu; % Example 1
    % alpha = (alpha_star * exp(1 / -mu) + 2 * exp(-mu)); % Example 2
    alpha = (alpha_star * mu + 2) / (mu + 1); % Example 3
end


