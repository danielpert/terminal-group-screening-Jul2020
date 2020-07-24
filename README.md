# new-terminal-group-screening

### extension of https://github.com/daico007/ChBE4830-modified-project2.git,
### https://github.com/summeraz/monolayer_screening,
### https://github.com/summeraz/terminal_group_screening, and
### https://github.com/summeraz/terminal_groups_mixed

This signac project involves the molecular simulation of two functionalized monolayers shearing against one another at different normal loads (5, 15, and 25 nN), an amorphous silica surface, and different random seeds like the original study (__J.Chem.TheoryComput.2020, 16, 1779−1793__). A chain length of 17 carbons is used.

However, this project is an extension of the previous perojects in that is uses the following new terminal groups:
* cyclohexyl
* nitroso
* acetylene

This workspace is meant to provide simulation data for these groups not screened previously, and not fed to the machine learning model. (https://github.com/daico007/tribology-machine-learning) Using this data, along with RDKit to determine the molecular descriptors, this will allow for the generalization ability of the random forest and neural networks to new data.

In this project, the 81 monolayer systems being studied include all 27 combinations of the three new terminal groups, with each combination being tested with 3 different random seeds to provide variation in the positioning of the alkylsilane chains on the silica surface.

Each monolayer system will be simulated at 3 different normal loads for a total of 273 simulations. The coefficient of friction and adhesive force for each normal load will be calculated by averaging the values from each random seed, and linear regression is used to determine the COF and F0 for each system.
