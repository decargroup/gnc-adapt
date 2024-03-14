function w = gm(epsilon, ~)
%GM weights from Geman-McClure robust loss function
%
% Inputs:
% - epsilon: residuals e' * W * e [N x 1]
%
% Editor: Kyungmin John Jung
% Date: 2024-03-12
% Lab: DECAR Group
% Institution: McGill University

w = ( 1 ./ (epsilon.^2 ./ 4 + 1) ).^2;
end

