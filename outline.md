# Project outline and plan:


## Science questions:

- Shape and density profile of the MW and LMC potential. Slices as a function of time.
  - Not oblate, prolate or triaxial.  (need an expansion) 
  - Inertia tensor try to fit the best ellipsoid. 

- With BFEs:
  - Can we identify Local and Global wake (maybe m=1 mode)? At t=0.
  - Can we identify Local and Global wake (maybe m=1 mode)? At all times.
  - Quantify the resonances in the halo? 


### Study case: Orbits of satellite galacies and GC vs static MW and static MW + LMC.

Additional ideas of what to do once this project is done are found [here](https://github.com/jngaravitoc/MW-LMC-SCF/blob/master/ideas.md)


## Methods:

### Implement the BFE expansion:

 - Noise in the coefficients.
    - Check computation, maybe there is a bug in the code.
    - Look at Wainberg variance method.
    
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
   



