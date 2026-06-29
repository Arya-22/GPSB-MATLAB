# GPSB-MATLAB
MATLAB implementation of Generalized Proximal Subgradient-Based (GPSB) algorithms (GPSB1 &amp; GPSB2) for non-smooth composite convex and weakly convex optimization problems.

    This repository provides the official MATLAB implementation for the Generalized Proximal Subgradient-Based (GPSB) algorithms (GPSB1 and GPSB2). The GPSB framework is designed to solve non-smooth composite convex and weakly convex optimization problems.
    The framework successfully generalizes several classical subgradient-based methods, including Dual Averaging (DA), Subgradient Method (SM), Regularized Proximal Subgradient Method (RPSB), and Incremental Proximal Subgradient Method (IPSB).
    📊 Numerical Experiments
        This repository includes code to reproduce the two main numerical experiments presented in the paper, strictly adhering to the reproducibility protocols (e.g., fixed random seeds, multiple independent trials, and identical parameter setups across benchmarks).
        Experiment 1: Elastic Net Regression (Convex Setting)
            <b>Problem:</b> High-dimensional regularized regression using the Elastic Net penalty (combining $\ell_1$ and $\ell_2$ regularizations).
            Characteristics: Non-smooth composite convex optimization.
            <b>Inner Solver:</b> The strongly convex subproblems involving non-smooth $\ell_1$ terms are efficiently solved using the Fast Iterative Shrinkage-Thresholding Algorithm (FISTA).
            Baselines: RPSB, Dual Averaging (DA), Proximal Subgradient Method (PSM), and standard Subgradient Method (SM).

        Experiment 2: Robust Sparse Recovery with Heavy-tailed Noise (Weakly Convex Setting)
            <b>Problem:</b> Signal recovery under heavy-tailed noise using a robust loss function blending $\ell_1$ and $\ell_2$ characteristics, combined with an $\sigma$ penalty.
            Characteristics: Non-smooth composite weakly convex optimization.
            Data Generation: Ground-truth signals are strictly sparse, and the measurement noise is generated using a heavy-tailed Student's t-distribution to properly evaluate the robust model.
            <b>Metrics:</b> Evaluates optimization quality and recovery accuracy, including objective values and reconstruction errors ($\|x_k - x_{true}\|_2$) [2].


    ⚙️ Requirements
        MATLAB (Tested on MATLAB R2020b or later)
        Statistics and Machine Learning Toolbox (required for data generation, e.g., mvnrnd, trnd)

