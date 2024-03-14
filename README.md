# Quantum Monte Carlo

A program which determines the ground state energy of the helium atom using a Quantum Monte Carlo approach.

## Description



In principle the ground state energy of a system can be detemined by solving Schrodinger's equation given the appropriate Hamiltonian. The eigenvalue equation yields the wavefunction as an eigenvector/eigenstate, and the energy as the corresponding eigenvalue. In practice anything more complicated than a hydrogen atom becomes an n-body problem, typically lacking a closed form solution. Even after constructing a good approximation for the wavefunction, the multidimensional integral necessary for determining the energy level remain problematic and computationally heavy. Instead, it is possible to find the energy level of a system with a given trial wavefunction stochastically, and apply a variational method to correct the wavefunction. This method sometimes known as a Variational Monte Carlo (VMC) allows us to find the ground state as the solution to an optimization problem, namely finding the minimum energy for the system given a wavefunction with a variational parameter.

The ground state energy will be the energy expectation value associated with the trial wavefunction for the system, and as such can be calculated by sampling the configuration space of the system according to its probability density function. This is achieved by taking an initially random homogeneous spherical distribution of electron pair positions, and updating the points according to a Metropolis-Hastings algorithm until they approach an energy equilibrium. Continuing to update the points this way after initially reaching equilibrium results in a Markov Chain of position data values which match the necessary distribution. This process, known as a Markov Chain Monte Carlo (MCMC), is particularly useful for the evaulation of multidimensional integrals. The average energy value for the system as determined by the electron positions in the Markov Chain gives the energy expectation value, or ground level energy.

```math
H = {-\hbar \over 2m} (\nabla_1^2 + \nabla_2^2) - {e^2 \over 4 \pi \epsilon_0}({2 \over r_1} + {2 \over r_2} - {1 \over |\vec{r_1} - \vec{r_2}|})
```

### The Trial Wavefunction

In place of an analytic solution for the wavefunction, we can canstruct a trial wavefunction as an approximation. A typical construction involves neglecting the inter-electron Coulomb interaction, and gives us a solution in the form of a product between two hydrogen-like wavefunctions with a nuclear charge of 2e. 

```math
\psi_{trial} = {8 \over \pi a^3} e^{-2(r_1+r_2) \over a}
```
A clever trick to improve the trial wavefunction involves introducing a variational parameter in a way that capture some of the information about the electron positions. Allowing the value for nuclear charge to vary, that is replacing it by a variational parameter q, gives a more flexible model. This new parameter q represents the effective nuclear charge seen by an electron, or in other words the charge that an  electron feels from the nucleus as the result of shielding from the second electron. This allows for a sort of indirect capture of electron-electron interaction, although it is still lacking angular dependence. Optimizing the new trial wavefunction for minimum energy expectation value yields the ground state energy.

```math
\psi_{trial} = {8 \over \pi a^3} e^{-q(r_1+r_2) \over a}
```



## Getting Started

### Dependencies



### Executing program

* Run main.py


## Authors

David Petrie (petrie.david.james@gmail.com)

## Acknowledgments
