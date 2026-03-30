# Efficient Methods for Solving Least Squares Problems

This repository accompanies a paper on **efficient numerical methods for solving least squares problems**, a fundamental task that appears across **data science, AI, and engineering**.


##  Why Least Squares Matters

Many real-world problems do not have exact solutions due to:
- Noise in data
- Overdetermined systems (more equations than unknowns)
- Measurement errors

In such cases, we seek an **approximate solution** that minimizes the error:

\[
\min_x \|Ax - b\|_2^2
\]

This is known as a **least squares problem**.


##  Applications Leading to Least Squares Problems

Least squares is not just a mathematical abstraction—it is the backbone of countless applications.

###  Data Science & Machine Learning
- **Linear Regression**  
  The most basic predictive model is directly formulated as a least squares problem.

- **Feature Fitting & Trend Analysis**  
  Fitting curves or models to data (polynomial regression, splines)

- **Dimensionality Reduction**  
  Methods like PCA rely on solving least squares subproblems.

- **Model Calibration**  
  Adjusting parameters to best fit observed data


###  Artificial Intelligence
- **Training Neural Networks (locally)**  
  Some layers or approximations reduce to least squares problems.

- **Computer Vision**
  - Image reconstruction
  - Camera pose estimation
  - 3D reconstruction (structure from motion)

- **Natural Language Processing**
  - Word embedding approximations
  - Matrix factorization methods


###  Engineering
- **Signal Processing**
  - Filtering and denoising signals
  - Fourier and wavelet approximations

- **Control Systems**
  - System identification (estimating system parameters)

- **Robotics**
  - Localization and mapping (SLAM)
  - Sensor fusion


###  Scientific Computing
- **Inverse Problems**
  Recovering unknown causes from observed effects (often ill-posed)

- **Physics & Simulation**
  Parameter estimation in models governed by differential equations

- **Geophysics & Medical Imaging**
  Reconstruction problems (e.g., MRI, CT scans)


##  Why Efficiency is Crucial

While least squares problems are easy to state, solving them efficiently is **non-trivial in practice**.

### 1. Scale of Modern Data
- Datasets can contain **millions or billions of samples**
- Matrices are often:
  - Large
  - Sparse
  - Ill-conditioned

Naive methods quickly become computationally infeasible.


### 2. Numerical Stability
- Directly solving normal equations:
  \[
  A^T A x = A^T b
  \]
  can lead to **severe numerical instability**

- Efficient methods must:
  - Avoid loss of precision
  - Handle ill-conditioned systems robustly


### 3. Real-Time Constraints
- Applications like robotics, finance, and online systems require:
  - Fast updates
  - Incremental solutions


### 4. Energy and Resource Constraints
- Large-scale computation impacts:
  - Energy consumption
  - Memory usage

Efficient algorithms reduce both computational cost and environmental impact.


##  Role of Numerical Linear Algebra

Efficient least squares solvers rely on advanced techniques from **numerical linear algebra**:

- **QR Decomposition**  
  Stable and widely used for solving least squares problems

- **Singular Value Decomposition (SVD)**  
  Provides robustness and insight into rank-deficient systems

- **Iterative Methods**
  - Conjugate Gradient (for normal equations)
  - LSQR and related algorithms

- **Randomized Algorithms**
  - Sketching methods for large-scale problems
  - Subsampling techniques

These methods trade off:
- Accuracy
- Speed
- Memory usage



##  Resources

For more details, refer to the full paper and implementation:

-  **Paper**: [Read the full paper](https://doi.org/10.1002/nla.70017)
-  **MATLAB Code**: [View the implementation and the codes of the example projects inside the paper](https://github.com/AminAPNA/Efficient-Methods-for-Solving-Least-Squares-Problems/tree/main/matlab_codes)
