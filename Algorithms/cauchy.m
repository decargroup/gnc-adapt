function w = cauchy(epsilon, ~)
%CAUCHY weights from Cauchy robust loss function
%
% Inputs:
% - epsilon: residuals e' * W * e [N x 1]
%
% Editor: Kyungmin John Jung
% Date: 2024-03-12
% Lab: DECAR Group
% Institution: McGill University

w = 1./((1/2).*epsilon.^2 + 1);
end

