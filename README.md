# The-Corridor-Method [![DOI](https://zenodo.org/badge/902868851.svg)](https://doi.org/10.5281/zenodo.14446619)

## Authors
[Ahmet Agaoglu](https://github.com/Ahmet-Agaoglu), Namik Ciblak

## Description
This repository contains the MATLAB implementation of a novel method for line fitting under bounded-error constraints. It features the Corridor Method for accurately delineating the feasible parameter set (FPS) and provides solutions for estimating regression parameters for cases with unknown error distributions and normally distributed errors.

## Files Included
The repository includes a MATLAB function:
- **CorridorMethod.m**: This function determines the FPS and computes:
  - The centroid of the FPS, which serves as the optimal estimate in cases where no distributional information about the errors is available.
  - The optimal solution for normally distributed errors.

Additionally, the function visualizes the FPS along with the centroid (blue circle), the Ordinary Least Squares solution (black star), and the optimal solution (red circle) for normally distributed errors.

### Function Inputs:
- x: Independent variable values.
- y: Dependent variable values.
- Lb: Lower bounds of the dependent variable values.
- Ub: Upper bounds of the dependent variable values.

### Function Outputs:
- V: Vertices of the FPS.
- cntrd: Centroid of the FPS (used for unknown error distributions).
- optSol: Optimal solution for normally distributed errors.

## Installation
1. Clone or download the repository to your system.
2. Set the downloaded folder as the MATLAB working directory.

### Dependencies
- MATLAB R2014 or later.
- Required Toolboxes:
  - Statistics and Machine Learning Toolbox

## Usage
To run the code, call the function `CorridorMethod` with the required inputs.

### Example
x = [1, 2, 3, 4, 5];
y = [2.1, 4.2, 5.9, 7.8, 9.9];
Lb = y - 0.5;
Ub = y + 0.5;

[V, cntrd, optSol] = CorridorMethod(x, y, Lb, Ub);

## Reproducibility
To reproduce the results from the accompanying study:

1. Prepare synthetic or real-world datasets with known bounds on measurement errors.
2. Run the `CorridorMethod` function on these datasets.
3. Compare the outputs (FPS, centroid, and optimal point) with baseline methods such as Ordinary Least Squares, Constrained Linear Least Squares, and Regularized Chebyshev Center.

## Citation
If you use this code in your research, please cite:
@article{Agaoglu2024, title={A Fast and Accurate Geometric Method for Line Fitting with Bounded Errors}, author={Ahmet Agaoglu and Namik Ciblak}, journal={Under Review}, year={2024}, doi={}}


## License
This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.
