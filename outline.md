# Project Outline and plan


## Science questions:

- Quantify the resonances in the halo?
- Does the m=1 mode is the disk motion with respect to the halo?
- What does it term in the MW's expansion means? Can we find the wake potential?**


### Time evolution of the Gravitational Potential

### Study case: Orbits of satellite galacies and GC.

Additional ideas of what to do once this project is done are found [here](https://github.com/jngaravitoc/MW-LMC-SCF/blob/master/ideas.md)


## Methods:

### Implement the BFE expansion:

 - Noise in the coefficients.
    - Check computation, maybe there is a bug in the code.
    - Look at Wainberg variance methodology.
    
 - How to choose the number of coefficients in each expansion.
    - Again think about the variance?
    
 - What particles from the LMC use to compute the potential from the LMC.
    - First use all of them and see how many terms are needed.
    - Use the bonund particles only.

 - Technical details:
    - Compute coefficients in parallel. 

## Products:

Release Python library to integrate orbits of Halo tracers.
The project lives [here](https://github.com/jngaravitoc/BFE_integrator)

## Timeline: 
  - Week 1 (02.25-03.01): 
    - Check Coefficients computtion if there is a bug in the code.
    - Compute coefficients in parallel. Run the code done in Heidelberg with Adrian in ocelote and desktop. Decide if something else needs to be done in thise end.
    - Read Weimberg's variance theory.
   



