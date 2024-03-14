# Quantum Monte Carlo

This program demonstrates the effectiveness of Quantum Monte Carlo methods in determining the ground states of small quantum systems, in this case by determining the ground state energy of a helium atom.

## Description



In principle the ground state energy of a system can be determined by solving Schrodinger's equation given the appropriate Hamiltonian, where the state and energy level come from the eigenvalue equation as the eigenstate and eigenvalues respectively. In practice anything more complicated than a hydrogen atom becomes an n-body problem, typically lacking a closed form solution. Even after constructing a good approximation for the wave function, the multidimensional integral necessary for determining the energy level remain problematic and computationally heavy. Instead, it is possible to find the energy level of a system with a given trial wave function stochastically, and apply a variational method to correct the wave function. This method sometimes referred to as as a Variational Monte Carlo (VMC) allows us to find the ground state as the solution to an optimization problem, namely finding the minimum energy for the system given a wave function with a variational parameter.

The ground state energy will be the energy expectation value associated with the trial wave function for the system, and as such can be calculated by sampling the configuration space of the system according to its probability density function. This is achieved by taking an initially random homogeneous spherical distribution of electron pair positions, and updating the points according to a Metropolis-Hastings algorithm until they approach an energy equilibrium. Continuing to update the points this way, after initially reaching equilibrium, results in a Markov Chain of position data values which obey the expected probability distribution. This process, known as a Markov Chain Monte Carlo (MCMC), is particularly useful for the evaluation of multidimensional integrals. The average energy value for the system as determined by the electron positions in the Markov Chain gives the energy expectation value, or ground level energy.

The following Hamiltonian defines the system

```math
H = {-\hbar \over 2m} (\nabla_1^2 + \nabla_2^2) - {e^2 \over 4 \pi \epsilon_0}({2 \over r_1} + {2 \over r_2} - {1 \over |\vec{r_1} - \vec{r_2}|})
```

### The Trial Wave function

In place of an analytic solution for the wave function, we can construct a trial wave function as an approximation. A typical construction involves neglecting the electron-electron Coulomb interaction, and gives us a solution in the form of a product between two hydrogen-like wave functions with a nuclear charge of 2e. 

```math
\psi_{trial} = {8 \over \pi a^3} e^{-2(r_1+r_2) \over a}
```
A common trick to improve the trial wave function involves introducing a variational parameter in a way that capture some of the information about the electron positions[1]. Allowing the value for nuclear charge to vary, that is replacing it by a variational parameter q, gives a more flexible model. This new parameter q represents the effective nuclear charge seen by an electron, or in other words the charge that an  electron feels from the nucleus as the result of shielding from the second electron. This allows for a sort of indirect capture of electron-electron interaction, although it is still lacking angular dependence. Optimizing the new trial wave function for minimum energy expectation value yields the ground state energy.

```math
\psi_{trial} = {8 \over \pi a^3} e^{-q(r_1+r_2) \over a}
```



## Getting Started

### Dependencies
* matplotlib 

* numpy

### Executing program

* Run main.py


## Authors

David Petrie (petrie.david.james@gmail.com)

## Acknowledgments
[1] Griffiths, David J. Introduction to quantum mechanics. Cambridge University Press, 2016.

[2] Julien Toulouse, Roland Assaraf, C. J. Umrigar. Introduction to the variational and diffusion
Monte Carlo methods. Advances in Quantum Chemistry, 2016, Electron Correlation in Molecules
ab initio Beyond Gaussian Quantum Chemistry, 73, pp.285.

[3] Nightingale, M. Peter, and Cyrus J. Umrigar, eds. Quantum Monte Carlo methods in physics
and chemistry. No. 525. Springer Science & Business Media, 1998.

