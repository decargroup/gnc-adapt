function params = gncParams(rlf, residuals, tau, alpha)
%GNC_PARAMS Initialize GNC parameters
%   Parameters include
%   - mu: convexity parameter
%   - c: continuation constant
%   - @update: function to increase convexity
%
% Inputs:
% - rlf: loss function to be incorporated with GNC
% - residuals: current residuals
% - alpha: the optimal shape parameter
%
% Editor: Kyungmin John Jung
% Date: 2024-03-12
% Lab: DECAR Group
% Institution: McGill University

switch upper(functions(rlf).function)
    case 'GNC_TLS'
        params.mu = max(1 / (2 * max(residuals) - 1), 1e-3);
        params.update = @(x) x * 1.4;
        params.converged = @(x) false;
    case 'GNC_ADAPT' % Example 2 or 3 (Increasing mu)
        params.mu = max(1 / (2 * max(residuals) - 1), 1e-3);
        params.update = @(x, c) x * 1.4;
        params.converged = @(x) x > 1e3;
        if isempty(alpha)
            params.alphaStar = getAlphaStar(residuals, -tau, tau);
        else
            params.alphaStar = alpha;
        end
        params.tau = tau;
    case 'GNC_AMB'
        params.mu = max(1 / (2 * max(residuals) - 1), 1e-3);
        params.update = @(x, c) x * 1.4;
        params.converged = @(x) x > 1e3;
        aStar = getAStar(residuals, 3, tau);
        params.epsilon_tilde = aStar * sqrt(3 - 1);
        epsilon_shift = residuals - params.epsilon_tilde;
        idx_xi = epsilon_shift >= 0;
        nu     = tau - params.epsilon_tilde;
        params.alphaStar = getAlphaStar(epsilon_shift(idx_xi), 0, nu);
        params.tau = tau;
    otherwise
        error('GNC parameters not Implemented')
end

end

