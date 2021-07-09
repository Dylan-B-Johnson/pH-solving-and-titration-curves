# pH-solving-and-titration-curves
Generates approximate pH-titration curves, based on initial conditions.
- plot_titration() plots the pH-titration curve, as well as returns the generated data.
- get_pH() returns the pH for a titration at a specific point (of volume added).
- react() solves Stoichiometric neutralizations.
- All three functions have comments explaining their parameters. 
- If you get weird squiggles, decreasing the increment may help, but it is an approximation and the areas where there is a transition between methods of calculation will result in larger errors.
- This module uses the calculation method taught in IB-HL Chemistry and therefore may not be best method available.
