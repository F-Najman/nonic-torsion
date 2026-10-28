# Least CM Degree

This repository holds Magma code and lists for the paper *The least degree of a CM point on a modular curve* by Pete L. Clark, Tyler Genao, Paul Pollack, and Frederick Saia. 

See **Required Lists** descriptions for dependencies. 

## Sporadic Checks

- `sporadic_checks_X0_X1.m`: Knowing that there are finitely many N for which the modular curves X_1(N) and X_0(N) may not have a sporadic CM point, this code is for checking such N further to see if we can guarantee a sporadic point. Specifically, we use the theorem of Faltings-Frey combined with lower bounds on the gonality of these curves and upper bounds on the least degree of a CM point on them. 

- `sporadic_checks_XMN.m`: Knowing that there are finitely many pairs (M,N) with M|N for which the modular curve X(M,N) may not have a sporadic CM point, this code is for checking such pairs further to see if we can guarantee a sporadic point. Specifically, we use the theorem of Faltings-Frey combined with a lower bound on the gonality of X(M,N) and an upper bound on the least degree of a CM point on X(M,N). 

- `further_bads_X0_X1.m`: The sequence of N for which techniques from 'sporadic_checks_X0_X1' do not guarantee a sporadic CM point on X_0(N) and X_1(N). 

- `further_bads_XMN.m`: The sequence of pairs (M,N) for which techniques from 'sporadic_checks_XMN' do not guarantee a sporadic CM point on X(M,N). 

## Least Degrees

### X1

- `least_degreesX1.m`: The aim of this code is to compute, for an integer N >= 2, the least degree over Q of a CM point on the modular curve X_1(N) (via methods of Bourdon-Clark '19). These computations are then used to try to guarantee the existence of a sporadic CM point on X_1(N) via Frey-Faltings and lower bounds on the gonality of X_1(N). We also prove there can be no sporadic CM point on X_1(N) for some N using the least degree computation combined with upper bounds on the gonality of X_1(N) from Derickx-van Hoeij '13. 

- `hyper_bads_X1.m`: Sequence of 297 naturals N in `further_bads_X0_X1.m` such that X_1(N) is guaranteed a sporadic CM point via methods in `least_degreesX1.m`. This corresponds to the set F_1 in our paper. 

- `no_sporadic_CM_X1.m`: Sequence of the 67 values of N for which we prove X_1(N) has no sporadic CM point. 

- `unknown_X1.m`: Sequence of the 227 values of N for which we do not know whether X_1(N) has a sporadic CM point. 

- `dcm_list_X1_1mil.m`: List of least CM degree values for X_1(N) for N up to 1,000,000, each including a minimizing order. 

- `dcm_list_all_min_orders_X1_10k.m`: List of least CM degree values for X_1(N) for N up to 10,000, with count of and complete list of all the minimizing orders.

### X0 

- `least_degreesX0.m`: The aim of this code is to compute, for an integer N >= 2, the least degree over Q of a CM point on the modular curve X_0(N). We directly use computations from `least_degreesX1.m` to do this. These computations are then used to try to guarantee the existence of a sporadic CM point on X_0(N) via Frey-Faltings and lower bounds on the gonality of X_0(N). We also prove there can be no sporadic CM point on X_0(N) for some N in cases where delta(X_0(N)) <= 2. 

- `hyper_bads_X0.m`: Sequence of 359 naturals N in `further_bads_X0_X1.m` such that X_0(N) is guaranteed a sporadic CM point via methods in `least_degreesX0.m`. This corresponds to the set F_0 in our paper.

- `no_sporadic_CM_X0.m`: Sequence of the 50 values of N for which we prove X_0(N) has no sporadic CM point. 

- `unknown_X0.m`: Sequence of the 106 values of N for which we do not know whether X_0(N) has a sporadic CM point. 

- `dcm_list_X0_1mil.m`: List of least CM degree values for X_0(N) for N up to 1,000,000, each including a minimizing order. 

- 'dcm_list_all_min_orders_X0_10k.m': List of least CM degree values for X_0(N) for N up to 10,000, with count of and complete list of all the minimizing orders.

### XMN

- `least_degreesXMN.m`: The aim of this code is to compute, for integers M,N with M|N and M>1, the least degree over Q of a CM point on the modular curve X(M,N) (via methods of Bourdon-Clark '19). These computations are then used to try to guarantee the existence of a sporadic CM point on X(M,N) via Frey-Faltings and lower bounds on the gonality of X(M,N). We also prove there can be no sporadic CM point on X(M,N) for some M,N using the least degree computation combined with upper bounds on the gonality of X(M,N) derived from those of X_1(N) Derickx-van Hoeij '13. 

- `hyper_bads_XMN.m`: Sequence of 480 pairs (M,N) in `further_bads_XMN.m` such that X(M,N) is not guaranteed a sporadic CM point via methods in `least_degreesXMN.m`.

- `no_sporadic_CM_XMN.m`: Sequence of the 37 pairs (M,N) with M|N and M>1 for which we prove X(M,N) has no sporadic CM point. 

- `unknown_XMN.m`: Sequence of the 146 pairs (M,N) with M|N and M>1 for which we do not know whether X(M,N) has a sporadic CM point. For pairs with M=1 corresponding to X(1,N) = X_1(N), see `unknown_X1.m`. 

- `dcm_bounds_list_XMN_100`: List of upper bounds on least CM degree values for X(M,N) for all M|N for N up to 100, including a list of all orders of class number up to 100 yielding this upper bound. For all N<53, these least CM degree values are exact. 

## Required Lists

- `fields.m`: List of quadratic number fields of discriminant -3 to -40,000 downloaded from the LMFDB on 23 July 2019. Required for  
  - `sporadic_checks_X0_X1.m`  
  - `sporadic_checks_XMN.m`  
  - `least_degreesX1.m`  
  - `least_degreesXMN.m`  
  
- `classnum_disc_list.m`: Sequence of absolute values of discriminants of imaginary quadratic fields of class number up to 100. Modified from list of M. Watkins. Required for  
  - `least_degreesX1.m` *Only* if using the dcm_over_orders function (currently commented out, superseded by dcm_over_nonmax_orders)
  
- `cond_disc_list_all_O.m`: List of all (not just maximal) imaginary quadratic orders of class number up to 100. Generated using list of maximal orders by M. Watkins. Required for  
  - `least_degreesX1.m`  
  - `least_degreesXMN.m` 
  - 'least_degreesX0.m' (only used for the dcm_X0_all_minimizing_orders function)

- `computing_nonmax_orders.m`: Used to generate list `cond_disc_list_all_0.m`.

- `further_bads_X0_X1.m`: Copy of file from **Sporadic_Checks**. Required for  
  - `least_degreesX1.m` if computing `hyper_bads_X1.m`  
  - `least_degreesX0.m` if computing `hyper_bads_X0.m`  

- `further_bads_XMN.m`: Copy of file from **Sporadic_Checks**. Required for  
  - `least_degreesXMN.m` if computing `hyper_bads_XMN.m`  

- `dcm_list_X1_1mil.m`: Copy of file from **Least Degrees : X1**. Required for   
  - `least_degreesX0.m` if computing `hyper_bads_X0` or using dcm_X0_all_minimizing_orders function

- `dcm_list_X0_1mil.m`: Copy of file from **Least Degrees : X0**. Required for   
  - `least_degreesX0.m` (only if using dcm_X0_all_minimizing_orders function)

- `dcm_list_all_min_orders_X1_10k.m`: Copy of file from **Least Degrees : X1**. Required for   
  - `least_degreesX0.m` (only if using dcm_X0_all_minimizing_orders function)

- `hyper_bads_X0.m`: Copy of file from **Least Degrees : X0**. Required for  
  - `least_degreesX0.m` if computing `no_sporadic_CM_X0.m` without first generating `hyper_bads_X0`  

- `hyper_bads_X1.m`: Copy of file from **Least Degrees : X1**. Required for  
  - `least_degreesX1.m` if computing `no_sporadic_CM_X1.m` without first generating `hyper_bads_X1`  

- `hyper_bads_XMN.m`: Copy of file from **Least Degrees : XMN**. Required for  
  - `least_degreesXMN.m` if computing `no_sporadic_CM_XMN.m` without first generating `hyper_bads_XMN`
