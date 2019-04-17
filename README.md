# Experiments for geometry of SaS problems

*Short-and-sparse* (SaS) inverse problems are frequently encountered in signal processing tasks across science and engineering, and are frequently posed as nonconvex optimization problems. Regardless of the specific formulation, the geometry of all such optimization problems are governed similarly by two main aspects of the SaS model: *problem symmetries* and *sparsity-coherence tradeoff*.

This repository contains experiments on the *short-and-sparse blind deconvolution* (SaS-BD) problem, designed to demonstrate the geometry of SaS inverse problems based on the two crucial aspects.


## Experiments

For a number of formulations, show that:

- `Todo:` Negative curvature exists at balanced points of shift space, pointing toward the shifts of **a_0**.
- `Todo:` Directions orthogonal to shifts of **a_0** in shift-space exhibit positive curvature.
- `Todo:` Negative curvature decreases when shift-coherence increases.


## Utils
We will use `sbd-tools` from `https://github.com/sbdsphere/sas-geometry.git`

- `Todo:` Tool to find saddle-point between shifts of **a_0**.
- `Todo:` Tool to determine negative / positive curvature.
- `Todo:` Projector onto shift space / tube.
- `Todo:` Solvers for different formulations.


# SaS-BD formulations

Here the tools for each formulation is written into a MATLAB class. Each formulation should allow for one to

- Initialize based on either **y**, (**a_0**, **x_0**), or something similar
- Probe function values *f*(**a**) for different **a**.
- Find the nearest critical point for given **a**.
- Get information about positive and negative curvature.
- Take descent steps using methods similar to
  - (accelerated) gradient descent
  - curvilinear search (to have more control over direction in saddle points)

The goal in `sas-geometry` is not recovery performance, but to get information about  the landscape geometry near shift-space.
