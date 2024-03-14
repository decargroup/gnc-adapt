<h1 style="padding-left: 1em; height: 70px;"> GNC-ADAPT </h1>


This repository contains the MATLAB implementation of **GNC (Graduated Non-Convexity)** on **ADAPT (Adaptive Loss)** by Barron described in the following paper:

- K. Jung, T. Hitchcox, and J. R. Forbes. "An Adaptive Graudated Nonconvexity Loss Function for Robust Nonlinear Least-Squares Solutions," 2024. arXiv:2305.06869

```bibtex
@misc{Jung2023Adaptive,
      title={An Adaptive Graduated Nonconvexity Loss Function for Robust Nonlinear Least Squares Solutions}, 
      author={Kyungmin Jung and Thomas Hitchcox and James Richard Forbes},
      year={2024},
      eprint={2305.06869},
      archivePrefix={arXiv},
      primaryClass={cs.RO}
}
```

## Quick-start

Open MATLAB and run
```matlab
example.m
```

## GNC-ADAPT Example

The code below demonstrates how to use GNC-ADAPT to solve a nonlinear least-squares problem.
The example problem used here is a **single rotation averaging problem**.
A user can modify the fields in `problem_params` to generate a problem with different parameters.
- N: number of measurements
- outlierRatio: the amount of outliers in the measurements
- noiseSigma: the standard deviation of the Gaussian noise in the measurements

A problem is then generated by executing
```matlab
problem = SingleRotationAveragingProblem(problem_params);
```

The initial guess for the rotation matrix is generated by adding a random rotation to the ground truth rotation matrix.
```matlab
initialNoise = 90; % degrees
rotAxis = normalize(randn(3, 1), 'norm');
rotAngle = deg2rad(initialNoise) * randn();
R_init = problem.R_gt * axang2rotm([rotAxis', rotAngle]);
```

Then, the problem is solved using the method of the user's choice (e.g., GNC-ADAPT).
```matlab
results.gnc_adapt = problem.solve( ...
    'lossFunction', @gnc_adapt, ...
    'alphaStar', [], ...
    'R_init', R_init ...
);
```

The problem can be solved by other loss functions such as
| Method     | Function Handle | alphaStar |
|------------|:---------------:|----------:|
| Cauchy     | @cauchy         | -         |
| GM         | @gm             | -         |
| ADAPT      | @adapt          | -         |
| AMB        | @amb            | -         |
| GNC-Cauchy | @gnc_adapt      | 0         | 
| GNC-GM     | @gnc_adapt      | -2        |
| GNC-ADAPT  | @gnc_adapt      | -         |
| GNC-AMB    | @gnc_amb        | -         |
| GNC-TLS    | @gnc_tls        | -         |


The `solve` method has the following optional parameters:
- `R_init`: the initial guess for the rotation matrix
- `lossFunction`: a function handle to the loss function to be used
- `inliers`: a logical vector indicating the inliers
- `maxIter`: the maximum number of iterations
- `tau`: the noise bound parameter to estimate the optimal shape parameter
- `costThreshold`: the cost threshold to stop the iterations
- `alphaStar`: the value of optimal shape parameter when it is known
- `verbose`: a flag to print the progress of the algorithm

## Acknowledgments

This work was partially funded by:

- Voyis Imaging Inc.
- Natural Sciences and Engineering Research Council of Canada (NSERC)
- McGill Engineering Doctoral Award (MEDA)

## License

[BSD License](LICENSE.BSD)

