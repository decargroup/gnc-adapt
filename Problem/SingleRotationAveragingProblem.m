classdef SingleRotationAveragingProblem
    %SINGLEROTATIONAVERAGINGPROBLEM
    %   Generation and solution of the rotation averaging problem
    %
    % Editor: Kyungmin John Jung
    % Date: 2024-03-12
    % Lab: DECAR Group
    % Institution: McGilll University

    properties
        N
        outlierRatio
        noiseSigma
        outlierIndices
        inlierIndices
        R_gt
        R_measurements
    end

    methods
        function obj = SingleRotationAveragingProblem(problem_params)
            %SINGLEROTATIONAVERAGINGPROBLEM
            %   Constructor generates the problem.
            %
            %   Inputs:
            %   - problem_params: parameters that contains
            %       - N: number of measurements
            %       - outlierRatio: outlier ratio between 0 and 1
            %       - noiseSigma: noise in degrees

            % check mandatory fields in problem parameters
            mandatory_fields = {'N', 'outlierRatio', 'noiseSigma'};
            if ~all(isfield(problem_params, mandatory_fields))
                error('Wrong problem parameters!');
            end

            % parse problem parameters
            N = problem_params.N;
            outlierRatio = problem_params.outlierRatio;
            noiseSigma = deg2rad(problem_params.noiseSigma);
            covInv = 1 / noiseSigma.^2 * eye(3);

            % generate a random groundtruth rotation
            R_gt = randomRotation();

            % generate N noisy measurements
            R_measurements = zeros(3, 3, N);
            for i=1:N
                while true
                    axis = normalize(randn(3, 1), 'norm');
                    rad = noiseSigma * randn();
                    rotation = axang2rotm([axis', rad]);
                    R_measurements(:, :, i) = R_gt * rotation;
                    dphi = rotm2rotvec(R_gt' * R_measurements(:, :, i));
                    if (dphi' * covInv * dphi) <= chi2inv(0.9973, 3)
                        break;
                    end
                end
            end

            % generate outliers
            numOutliers = round(N * outlierRatio);
            for i=1:numOutliers
                while true
                    R_tmp = randomRotation();
                    dphi = rotm2rotvec(R_gt' * R_tmp);
                    if (dphi' * covInv * dphi) > chi2inv(0.9973, 3)
                        R_measurements(:, :, i) = R_tmp;
                        break;
                    end
                end
            end

            obj.N = N;
            obj.outlierRatio = outlierRatio;
            obj.noiseSigma = noiseSigma;
            obj.outlierIndices = 1:numOutliers;
            obj.inlierIndices = numOutliers+1:N;
            obj.R_gt = R_gt;
            obj.R_measurements = R_measurements;
        end

        function resultInfo = solve(obj, varargin)
            %SOLVE 
            %   Solve the problem in a nonlinear least-squares formulation.
            %
            %   Options:
            %   - R_init: initial rotation estimate (default: eye(3))
            %   - lossFunction: robust loss function (default: @l2)
            %   - inliers: known inlier indices (default: [])
            %   - maxIter: maximum GN iteration (default: 50)
            %   - tau: parameter for computing alphaStar (default: 40)
            %   - costThreshold: threshold to check convergence
            %                     (default: 1e-6)
            %   - use_gnc: flag to use GNC or not (default: false)
            %   - alphaStar: optimal shape parameter (default: [])
            %   - verbosity: print optimization process (default: false)

            loss_functions = { ...
                'l2', 'cauchy', 'gm', 'adapt', 'amb', ...
                'gnc_adapt', 'gnc_amb', 'gnc_tls'
            };

            % Parse optional arguments
            params = inputParser;
            params.CaseSensitive = false;
            params.addParameter('R_init', eye(3), ...
                @(x) all(ismembertol(x' * x, eye(3), 1e-6), 'all'));
            params.addParameter('lossFunction', @l2, ...
                @(x) any(strcmpi(functions(x).function, loss_functions)))
            params.addParameter('inliers', [], ...
                @(x) isnumeric(x) && isvector(x) && all(floor(x)==x));
            params.addParameter('maxIter', 50, ...
                @(x) isPositiveIntegerValuedNumeric(x));
            params.addParameter('tau', 40, ...
                @(x) isPositiveIntegerValuedNumeric(x) && x > 1);
            params.addParameter('costThreshold', 1e-6, ...
                @(x) isnumeric(x) && isscalar(x))
            params.addParameter('alphaStar', [], @(x) x < 2);
            params.addParameter('verbose', false, @(x) islogical(x));
            params.parse(varargin{:});

            R_init = params.Results.R_init;
            rlf = params.Results.lossFunction;
            inliers = params.Results.inliers;
            maxIter = params.Results.maxIter;
            tau = params.Results.tau;
            costThreshold = params.Results.costThreshold;
            alphaStar = params.Results.alphaStar;
            verbosity = params.Results.verbose;

            % Check the algorithm if it runs GNC
            use_gnc = false;
            if contains(upper(functions(rlf).function), 'GNC')
                use_gnc = true;
            end

            % Get measurements to use
            if isempty(inliers)
                inliers = 1:obj.N;
            end
            R_meas = obj.R_measurements(:, :, inliers);

            % measurement information matrix
            W = 1 / obj.noiseSigma.^2 * eye(3);

            % print
            if verbosity
                fprintf('\n===================================== %s =====================================\n', upper(functions(rlf).function))
                fprintf('------------------------------------------------------------------------------------------\n')
                if use_gnc
                    fprintf(' itr \t | \t obj \t | \t delobj  |   delobj/obj  | \t dx \t | \t sumw \t | \t mu \t | \n');
                else
                    fprintf(' itr \t | \t obj \t | \t delobj  |   delobj/obj  | \t dx \t | \t sumw \t | \n');
                end
            end

            % Run GN
            resultInfo = struct();
            iter = 0;
            dx = 0;
            relCostDiff = 0;
            absCostDiff = 0;
            prevCost = 0; 
            R_updated = R_init;
            R_history = R_init;
            while iter < maxIter
                % evaluate residuals and jacobians
                e = rotm2rotvec(mtimesx(R_updated', R_meas));
                J = -leftJacobianInverse(e);

                % formulate A and b matrices
                A = num2cell(J, [1, 2]);
                A = vertcat(A{:});
                b = reshape(e, [], 1);

                % compute cost
                e = reshape(e, 3, 1, []);
                residuals = squeeze(mtimesx(e, 'T', mtimesx(W, e)));
                currentCost = sum(residuals);
                costHistory(iter+1) = currentCost;

                % get weights
                if use_gnc
                    if iter == 0
                        gnc_params = gncParams(rlf, residuals, tau, alphaStar);
                    end
                    weights = rlf(residuals, gnc_params);
                else
                    weights = rlf(residuals, tau);
                end

                % check convergence
                if iter == 0
                    if currentCost == 0
                        resultInfo.status = 'AbsoluteCost';
                        break;
                    end
                else
                    absCostDiff = prevCost - currentCost;
                    relCostDiff = absCostDiff / prevCost;
                    if norm(dx) < 1e-6
                        resultInfo.status = 'SmallIncrement';
                        break;
                    end
                    if abs(absCostDiff) < costThreshold
                        resultInfo.status = 'AbsoluteCostDiff';
                        break;
                    end
                    if abs(relCostDiff) < costThreshold
                        resultInfo.status = 'RelativeCostDiff';
                        break;
                    end
                end

                if verbosity
                    if use_gnc
                        fprintf('%2i \t | %3.4e \t | %3.4e \t | %3.4e \t | %3.4e \t | %3.4e \t | %1.3e \t |\n',...
                            iter, currentCost, absCostDiff, max(abs(relCostDiff), 0), vecnorm(dx), sum(weights), gnc_params.mu);
                    else
                        fprintf('%2i \t | %3.4e \t | %3.4e \t | %3.4e \t | %3.4e \t | %3.4e \t |\n',...
                            iter, currentCost, absCostDiff, max(abs(relCostDiff), 0), vecnorm(dx), sum(weights));
                    end
                end

                % solve
                p = colamd(repelem(weights, 3, 1) .* A);
                Ap = A(:,p);
                dx = -(Ap' * Ap) \ (Ap' * (repelem(weights, 3, 1) .* b));

                inverseIdx = zeros(length(p), 1);
                for i=1:length(p)
                    inverseIdx(i) = find(p == i);
                end
                dx = dx(inverseIdx);

                % update states
                dx = full(reshape(dx, 3, []));
                R_updated = R_updated * rotvec2mat3d(dx);

                % update mu
                if use_gnc
                    if ~gnc_params.converged(gnc_params.mu)
                        gnc_params.mu = gnc_params.update(gnc_params.mu);
                    end
                end

                % increment
                iter = iter + 1;
                prevCost = currentCost;
                R_history(:, :, iter+1) = R_updated;
            end

            if iter == maxIter
                resultInfo.status = 'MaxIterations';
            end

            if verbosity
                absCostDiff = prevCost - currentCost;
                relCostDiff = absCostDiff / prevCost;
                if use_gnc
                    fprintf('%2i \t | %3.4e \t | %3.4e \t | %3.4e \t | %3.4e \t | %3.4e \t | %1.3e \t |\n',...
                            iter, currentCost, absCostDiff, max(abs(relCostDiff), 0), vecnorm(dx), sum(weights), gnc_params.mu);
                else
                    fprintf('%2i \t | %3.4e \t | %3.4e \t | %3.4e \t | %3.4e \t | %3.4e \t |\n',...
                            iter, currentCost, absCostDiff, max(abs(relCostDiff), 0), vecnorm(dx), sum(weights));
                end
            end

            resultInfo.iter = iter+1;
            resultInfo.R_optimized = R_updated;
            resultInfo.R_history = R_history;
            resultInfo.costHistory = costHistory;
            resultInfo.w_optimized = weights;
        end
    end
end

function R = randomRotation(varargin)
    params = inputParser;
    params.CaseSensitive = false;
    params.addParameter('RotationBound', 2*pi, @(x) 0 <= x && x <= 2*pi);
    params.parse(varargin{:});
    rotationBound = params.Results.RotationBound;

    angle = rotationBound*rand - rotationBound/2;
    axis  = randn(3, 1);
    axis  = axis / norm(axis);
    R     = axang2rotm([axis' angle]);
end

function J = leftJacobianInverse(phis)
phis  = reshape(phis, 3, 1, []);
angle = sqrt(mtimesx(phis, 'T', phis));

% when angle is small, compute Jinv using Taylor series expansion
tf_small = angle <= 1e-6; % small angle approximation
A        = zeros(1, 1, size(phis, 3));

% Taylor series expansion for small angle, see (163)
t2           = angle(tf_small).^2;
A(tf_small)  = (1/12) .* (1 + t2./60 .* (1 + t2./42.*(1 + t2./40)));

% in the regular case, compute using (125)
big_angle    = angle(~tf_small);
A(~tf_small) = 1./big_angle.^2 .*...
    (1 - ( big_angle.*sin(big_angle) ./ (2.*(1 - cos(big_angle))) ));

% create return using (125), `V' in Eade
J = eye(3) - 1/2 .* wedge(squeeze(phis)) ...
               + A .* mtimesx(wedge(squeeze(phis)), wedge(squeeze(phis)));
end

function elements_so3 = wedge(phis)
p1 = reshape(phis(1,:), 1, 1, []);
p2 = reshape(phis(2,:), 1, 1, []);
p3 = reshape(phis(3,:), 1, 1, []);
z  = zeros(1, 1, size(phis, 2));

% construct the element of the Lie algebra.
elements_so3 = [ z, -p3, p2  ;
                 p3, z, -p1  ;
                -p2, p1, z  ];
end