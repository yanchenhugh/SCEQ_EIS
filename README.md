# SCEQ_EIS
The code package for the project of SCEQ-EIS algorithm

# files:
## Code files:
*IRRf_can_error.gms: simulation and error checking for the model w/o EIS, or the Canonical Information Set (CIS) \n
*IRRf_eis_error.gms: simulation and error checking for the model w EIS
*IRR_params_simp.gms: the parameter file

## Simulation Results files:
*model_sols_can_ghq7.csv: results for the model w/o EIS, 1000 simulations.
**Column rgm1: 1 if the investment constraint is binding, 0 if nonbinding.
**Column succ: 1 if the simulation solution is successful, 0 if fails.
**Column succ err: 1 if solutions for all GHQ nodes in error checking are successful, 0 OW.
*model_sols_eis_ghq7_x.csv: results for the model w EIS, totally 1481 simulations. 

## Result retrieving files:
*get_Euler_error: .py and jupyter notebook file, the file to extract Euler errors from the simulation results files, written in python
*euler_err.txt: Euler error table, same w/ Table 1 in the paper.

## SCEQ_EIS.pdf: the paper draft.
