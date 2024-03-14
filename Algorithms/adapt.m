function w = adapt(epsilon, tau)
%ADAPT Chebrolu's Adaptive kernel [11]
% Implementation of:
%   - N. Chebrolu, T. Labe, O. Vysotska, J. Behley, and C. Stachniss,
%       "Adaptive robust kernels for non-linear least squares problems,"
%       IEEE Robotics and Automation Letters, vol. 6, pp. 2240-2247, 2 Apr.
%       2021
%
% Inputs:
% - epsilon: residuals e' * W * e [N x 1]
% - tau: parameter to compute alphaStar
%
% Editor: Kyungmin John Jung
% Date: 2024-03-12
% Lab: DECAR Group
% Institution: McGill University

% get alpha star
alphaStar = getAlphaStar(epsilon, -tau, tau);

% assign the weights
w = adaptiveWeights(epsilon, alphaStar);
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

