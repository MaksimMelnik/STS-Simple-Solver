# STSe-01D-Solver
Universal State-to-State modeling code for strong non-equilibrium. Allows to simulate 1D and 0D chemical-electronic-vibrational relaxation. It's important feature is accounting of electronic excitation. Only for atoms and diatomic molecules.
## Objectives
The code aims to **simulate simple 0D and 1D problems using the detailed state-to-state approach** [1]. This approach considers each vibrational state of molecules as separate gas species. That provides an accurate instrument for the description of strong non-equilibrium in gas mixtures with detailed vibrational and chemical kinetics. The code is flexible in a choice of accounted processes, models for its description, problems formulations and mixture compositions. The described enhances and simplifies the **validation of models and model combinations comparing with experimental data and other models**. The last one is the second important feature of the code.
## Restrictions
The kinetic scheme includes following restrictions. Particles under consideration have to be diatomic molecules, atoms or electrons. Polyatomic molecules are not implemented. The problem has to be 0D, 1D or reducible to a system of ODEs, since only the ODE solver is available. STS for rotational degrees of freedom is not available in the present product.
## Requirements
The code was tested on MATLAB 2020a [2] and newest versions. 
Symbolic Math Toolbox [3] is required for some test cases such as 'O2O_SW_Shatalov.m'.
To simulate test cases with free electrons LoKI-B tool [4, 5] is required.
## Structure
The code is divided into four sections: /data, /examples, /plots and /src. /data contains information about particles, chemical reactions and experimental studies. All resource files designated to evaluate kinetic values, coefficients and functions are stored in /src. /examples provides prepared test cases for validation based on experimental studies as well as theoretical test cases for wider investigations. The rest folder /plots includes scripts for plotting the results against experimental studies, their analysis and validation.
## How to use?
Install MATLAB and run one of the /examples (better to start with 'O2O_SW_Shatalov.m'). You also may use /src functions, /plots and everything else for your purposes separately.
## Where is used?
- [M Yu Melnik and E V Kustova 2022 J. Phys.: Conf. Ser. 2308 012014] (DOI.org/10.1088/1742-6596/2308/1/012014)

- [M Melnik and E Kustova 2021 J. Phys.: Conf. Ser. 1959 012034] (DOI.org/10.1088/1742-6596/1959/1/012034)
## References
[1] [E. Nagnibeda, E. Kustova, Nonequilibrium Reacting Gas Flows. Kinetic Theory of Transport and Relaxation Processes, Springer-Verlag, Berlin, Heidelberg, 2009.](https://doi.org/10.1007/978-3-642-01390-4).

[2] [MATLAB website.] (https://www.mathworks.com/products/matlab.html)

[3] [Symbolic Math Toolbox] (https://www.mathworks.com/products/symbolic.html)

[4] [Tejero A et al "The LisbOn KInetics Boltzmann solver" 2019 Plasma Sources Sci. Technol. 28 043001] (https://doi.org/10.1088/1361-6595/ab0537)

[5] [Tejero A et al "On the quasi-stationary approach to solve the electron Boltzmann equation in pulsed plasmas" 2021 Plasma Sources Sci. Technol. 30 065008] (https://doi.org/10.1088/1361-6595/abf858)
[available as open-access papers]