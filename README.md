# STSe-01D-Solver
Universal State-to-state modeling code. Allows to simulate 1D and 0D chemical-electronic-vabrational relaxation. It's importans fiture is accounting of electronic excitation. Only for atoms and diatomic molecules.
### todo
- add NO parameters (FHO parameters too)
- add NO exhange reactions
- finish NO Streicher RefSW experiments test cases
- add universal exchange reactions
- finish the Hubner 2012 experiment DC test case
- rewrite SW conditions recalculation before and after a SW (sw_cond_ar_f_eq, sw_cond_ar_f2, in_con_Ar, in_con_exp_p1v1T2, in_con_O2), add this recalculation to Shatalov's test case
- rewrite or exclude par_shatalov_f
- add a warning function to check for possible errors in the kinetic scheme
- add the universal output format
- add the universal plot builder
- add a comfy Arrhenius law switcher
- add Billing's VT, VV models
- add P4E VT model
- add Annusova's VT, VV model
- make separate k_diss rate coefficients functions (and rewrite k_diss_Aliat, k_diss_old, K_diss_Savelev)
- rewrite or exclude all R_diss functions (R_diss, R_diss_old, R_diss_Savelev(leave only k_diss_Savelev))
- transfer previous work for O2 isothermal discharge
- add the calculation of vibrational levels number below dissociation energy
- add an universal VE exchange function (k_VE and R_VE)
- check Oblapenko's FHO parameters
- add a separate functions for lambda