[![arXiv][arxiv-shield]][arxiv-url]
[![DOI][doi-shield]][doi-url]

# [sDAE-MOR][arxiv-url]
The folder contains the Matlab codes to perform Model Order Reduction (MOR), based on balancing, for control systems of switched Differential-Algebraic Equations (sDAE). The method is described in the paper [Balancing-based model reduction for switched descriptor systems][arxiv-url].


##Code info:

* **MOR_sDAE**: the main script to run, it constructs the projection spaces used to build the reduced system. Two test problems are included:
  * Constrained mass-sping system
  * Instationary Stokes control system

* **Compute_RS**: run after MOR_sDAE to compute and visualize the reduced solution and the decay of the reduction error with respect to the size of the reduced problem.

Functions:

* **Solve\_LS_GLE**: solve large-scale Generalized Lyapunov Equations----> $AX+X A^T +\Sigma_{j=1}^{M} (D_j X D_{j}^{T}+B_j B_{j}^{T})=0$ 
* **solve_KS**: solve large-scale Lyapunov equation of type $AX+XA^T+BB^T=0$ using standard Krylov method; 
* **solve\_KS_t**: solve large-scale Lyapunov equation of type $A^T X+XA+C^T C=0$ using standard Krylov method.


## Citing
If you use this project for academic work, please consider citing our
[publication][arxiv-url]:

    M. Manucci and B. Unger
    Balancing-based model reduction for switched descriptor systems
    ArXiv e-print 2404.10511, 2024.
    
## License
Distributed under the MIT License. See `LICENSE` for more information.


## Contacts

* Mattia Manucci - [mattia.manucci@simtech.uni-stuttgart.de](mattia.manucci@simtech.uni-stuttgart.de)
* Benjamin Unger - [benjamin.unger@simtech.uni-stuttgart.de](benjamin.unger@simtech.uni-stuttgart.de)



[doi-shield]: https://img.shields.io/badge/DOI-10.5281%20%2F%20zenodo.8335231-blue.svg?style=for-the-badge
[doi-url]: https://zenodo.org/records/10948132
[arxiv-shield]: https://img.shields.io/badge/arXiv-2204.13474-b31b1b.svg?style=for-the-badge
[arxiv-url]:http://arxiv.org/abs/2404.10511







