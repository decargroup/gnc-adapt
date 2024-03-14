function w = gnc_tls(residuals, params)
%GNC_TLS GNC weights update for the TLS cost function
% Implementation of
%   - H. Yang, P. Antonante, V. Tzoumas, and L. Carlone, "Graduated non-
%       convexity for robust spatial perception: From non-minimal solvers
%       to global outlier rejection," IEEE Robotics and Automation Letters,
%       vol. 5, pp. 1127–1134, 2 Apr. 2020.
%
% Editor: Kyungmin John Jung
% Date: 2024-03-12
% Lab: DECAR Group
% Institution: McGill University

mu = params.mu;
barc2 = sqrt(chi2inv(0.9973, 3));
th1 = (mu+1) / mu * barc2;
th2 = (mu) / (mu+1) * barc2; % th1 > th2
w = ones(length(residuals), 1);
for k = 1:length(residuals)
    if residuals(k) - th1 >= 0
        w(k) = 0;
    elseif residuals(k) - th2 <= 0
        w(k) = 1;
    else
        w(k) = sqrt( barc2*mu*(mu+1) / residuals(k) ) - mu;
        assert(w(k)>= 0 && w(k) <=1, 'weight %g is wrong!', w(k));
    end
end
end

