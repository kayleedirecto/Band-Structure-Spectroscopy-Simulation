# Band-Structure-Spectroscopy-Simulation

This was my research during my summer 2024 NSERC USRA! I had the opportunity to learn about the field of photonics and hone my programming skills to simulate a fibre optic experiment. 

## The Background 

This work builds off of the research published in 2019 by Avik Dutt et al. 

### Synthetic Dimensions 

In the field of photonics, there has been increasing research in synthetic dimensions, which mimic additional spatial dimensions. Combining these with geometric space allows researchers to study higher dimensional physics, with lower dimensional setups. The concept has interesting applications in topological physics, quantum information processing, and more (Dutt et al., 2019). 

To create a synthetic dimension, a non-spatial degree of freedom of a system is manipulated in such a way that it forms a lattice. For example, one can use the frequency of light to form synthetic frequency dimensions, which was the focus of this work. A common approach uses a ring resonator integrated with a phase modulator. The ring resonator naturally supports resonant frequency modes, and by modulating the resonator, these frequencies can be coupled, creating artificial connectivity between them to form a synthetic lattice. Varying the modulations leads to different coupling, producing different dimensionalities and lattice geometries (Yuan et al., 2018). 

Synthetic lattices are characterized by a band structure, similar to structures in solid-state physics. By determining the band structure of the system through rigorous arithmetic, researchers can determine the lattice geometries that should be created by the modulations they have introduced. Research performed by Dutt et al. (2019) has determined a method to experimentally readout the band structure of a synthetic lattice, allowing them to verify their experimental results with theory. By dynamically modulating a ring resonator, they created various lattice geometries, and showed that time-resolved transmission measurements of the system provided a direct readout of the band structure which matched theoretical results. The experimental setup can be seen below. Details on the exact method can be found in their paper.

![Screenshot 2025-06-05 at 1 58 32 PM](https://github.com/user-attachments/assets/431c3c51-ea18-47e4-8e89-b67f22732d1d)

The math to determine band structure becomes increasingly difficult for more complicated lattice shapes, so we wanted a simpler way to verify experimental results. The goal of my work was to create an exact simulation of the experimental setup for band structure readout. Therefore, before using the experiment to create new geometries, one could verify and explore theories using this simulation.

## The Algorithm  

We use the nonlinear Schrödinger equation (NLSE) to model the propagation of light though the optical cable, and numerically solve it using the Split Step Fourier Method (SSFM). 

### SSFM 

The nonlinear Schrödinger equation can be written as 

$\frac{\delta A}{\delta z} = - \frac{i \beta_2}{2} \frac{\delta^2 A}{\delta t^2} + i \gamma|A|^2A = [\hat{D} + \hat{N}]A$ 


## The Results 

## The References 

Dutt, A., Minkov, M., Lin, Q. et al. Experimental band structure spectroscopy along a synthetic dimension. Nat Commun 10, 3122 (2019). https://doi.org/10.1038/s41467-019-11117-9

Yuan, Luqi & Lin, Qian & Xiao, Meng & Fan, Shanhui. (2018). Synthetic dimension in photonics. Optica. 5. 1396-1405. 10.1364/OPTICA.5.001396. 





