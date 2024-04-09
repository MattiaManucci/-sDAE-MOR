# sDAE-MOR
The folder contains the Matlab codes to perform Model Order Reduction (MOR), based on balancing, for control systems of switched Differential Algebraic Equation (sDAE). The method is described in **[1]**.


Main scripts:

* **MOR_sDAE**: the main script to run, it constructs the projection spaces used to build the reduced system. Two test problem are included:
  * Constrained mass-sping system
  * Instationary Stokes control system

* **Compute_RS**: run after MOR_sDAE to compute and visualize the reduced solution and the decay of the reduction error with respect to the size of the reduced problem.

Auxiliary scripts:

* **Solve\_LS_GLE**: solve large-scale Generalized Lyapunov Equations---->AX+XA^T+\sum_{j=1}^{M}(D\_jXD\_j^T+B\_jB\_j^T)=0
* **solve_KS**: solve large-scale Lyapunov equation using standard Krylov method: AX+XA^T+BB^T=0
* **solve\_KS_t**: solve large-scale Lyapunov equation using standard Krylov method: A^TX+XA+C^TC=0


References

**[1]** To be inserted...



