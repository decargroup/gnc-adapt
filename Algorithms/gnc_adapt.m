function w = gnc_adapt(residuals, params)
%GNC_ADAPT GNC weights update for the ADAPT cost function
%   The proposed GNC with Barron's Adaptive loss (GNC-ADAPT)
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

mu = params.mu;
alpha = computeAlpha(mu ,params.alphaStar);
w = ones(length(residuals), 1);
for k = 1:length(residuals)
    if alpha <= 2
        switch alpha
            case 2
               w(k) = 1;
            case 0
               w(k) = 1 / (0.5 * residuals(k)^2 + 1);
            case -inf
               w(k) = exp(-0.5 * residuals(k)^2);
            otherwise
               w(k) = (residuals(k)^2 / abs(alpha - 2) + 1)^(alpha/2 - 1);
        end
        assert(w(k) >= 0 && w(k) <= 1, 'weight %g is wrong!', w(k));
    else
        error("Invalid mu value!");
    end
end

end

function alpha = computeAlpha(mu, alpha_star)
    % alpha = (alpha_star + 2 * mu - 2) / mu; % Example 1
    % alpha = (alpha_star * exp(1 / -mu) + 2 * exp(-mu)); % Example 2
    alpha = (alpha_star * mu + 2) / (mu + 1); % Example 3
end

