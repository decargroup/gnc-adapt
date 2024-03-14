%% Example: Single Rotation Averaging Problem
% Application of GNC-ADAPT to nonlinear least-squares problem compared to
% other state-of-the-art methods
%
% Editor: Kyungmin John Jung
% Date: 2024-03-12
% Lab: DECAR Group
% Institution: McGilll University

clc; clear; close all; restoredefaultpath;
addpath(genpath('./Algorithms'));
addpath(genpath('./mtimesx'));
addpath(genpath('./Problem'));
addpath(genpath('./utils'));

% seed random
rng(100);

% Problem generation
problem_params.N = 1000; % number of measurements available
problem_params.outlierRatio = 0.5; % outlier ratio
problem_params.noiseSigma = 5; % in degrees
problem = SingleRotationAveragingProblem(problem_params);

%% Run experiment

% To disable the execution of any algorithm feel free to comment out the
% relative section. `displayResults` will show only available results.

% initial condition
initialNoise = 90; % degrees
rotAxis = normalize(randn(3, 1), 'norm');
rotAngle = deg2rad(initialNoise) * randn();
R_init = problem.R_gt * axang2rotm([rotAxis', rotAngle]);

% solve
results.cauchy = problem.solve( ...
    'lossFunction', @cauchy, ...
    'R_init', R_init ...
);

results.gm = problem.solve( ...
    'lossFunction', @gm, ...
    'R_init', R_init ...
);

results.adapt = problem.solve( ...
    'lossFunction', @adapt, ...
    'R_init', R_init ...
);

results.amb = problem.solve( ...
    'lossFunction', @amb, ...
    'R_init', R_init ...
);

results.gnc_cauchy = problem.solve( ...
    'lossFunction', @gnc_adapt, ...
    'alphaStar', 0, ... % Cauchy
    'R_init', R_init ...
);

results.gnc_gm = problem.solve( ...
    'lossFunction', @gnc_adapt, ...
    'alphaStar', -2, ... % GM
    'R_init', R_init ...
);

results.gnc_adapt = problem.solve( ...
    'lossFunction', @gnc_adapt, ...
    'R_init', R_init ...
);

results.gnc_amb = problem.solve( ...
    'lossFunction', @gnc_amb, ...
    'R_init', R_init ...
);

results.gnc_tls = problem.solve( ...
    'lossFunction', @gnc_tls, ...
    'R_init', R_init ...
);

%% Display Results
displayResults(results, problem);