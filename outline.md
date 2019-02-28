# Project outline and plan:

The paper is [here](https://github.com/jngaravitoc/MW-LMC-SCF/blob/master/doc/BFE_MWLMC.pdf)


## Science questions and goals:

- What is the shape of the density profile and the potential of the MW and LMC. 
  - Not oblate, prolate or triaxial.  (need an expansion) 
  - Use the inertia tensor fit the best ellipsoid and argue the need of the expansion. 
  - Use the intertia tensor to define a coordenate sistem to apply the BFEs.


- With BFEs:
  - Compare the density and potential with those fomr the Sims.
  - If the coefficients are changing if the coordinate system is rotating.
  - What each terms mean in the expansion.
  - Can we identify Local and Global wake (maybe m=1 mode)? At t=0.
  - Can we identify Local and Global wake (in general)? At all times.
  - Quantify the resonances in the halo? 
  - Integrate orbits of particles and compare with the orbits from Gadget.


### Study case: Orbits of satellite galaxies and GC vs static MW and static MW + LMC.

Additional ideas of what to do once this project is done are found [here](https://github.com/jngaravitoc/MW-LMC-SCF/blob/master/ideas.md)


## Methods:

### Implement the BFE expansion:

 - Why high order l coefficients are so large ([notebook](https://github.com/jngaravitoc/SCF_tools/blob/master/code/computing_coefficients.ipynb)).
    - Check computation, maybe there is a bug in the code.
    - Does these coefficients matter for the potential.
    - Look at Weinberg variance method.
    - Understand the Noise in the expansions.
    
 - How to choose the number of coefficients in each expansion.
    - Again, think about the variance.
    
 - What particles should we use from the LMC use to compute the potential from the LMC.
    - First use all of them and see how many terms are needed.
    - Use the bound particles only.

 - Technical details:
    - Compute coefficients in parallel. 

## Products:

Release Python library to integrate orbits of Halo tracers.
The project lives [here](https://github.com/jngaravitoc/BFE_integrator)

## Timeline: 
  - Week 1 (02.25-03.01): 
    - Check Coefficients computation if there is a bug in the code.
    - Compute coefficients in parallel. Run the code done in Heidelberg with Adrian on ocelote and desktop. Decide if something else needs to be done on this end.
    - Read Weinberg's variance theory.
   



