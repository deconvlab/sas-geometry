# SaS-BD formulations

Here the tools for each formulation is written into a MATLAB class. Each formulation should allow for one to

- Initialize based on either **y**, (**a_0**, **x_0**), or something similar
- Probe function values *f*(**a**) for different **a**.
- Find the nearest critical point for given **a**.
- Get information about positive and negative curvature.
- Take (accelerated) gradient descent steps or something similar.



The goal in `sas-geometry` is not recovery performance, but to get information about  the landscape geometry near shift-space.