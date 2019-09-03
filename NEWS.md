# lg 0.4.0

   * Added bootstrap confidence intervals as option in dlg.
   * Added the estimation method "trivariate" for joint estimation and independence testing 
     with trivariate data.

# lg 0.3.0
 
   * Added contagion test and independence tests.
   * Added bw_simple for quick construction of bandwidth object.
   * Switches to old logspline algorithm if the new fails (very rare). 
   * Better NA handling, especially in the test for conditional independence

# lg 0.2.0 

   * Added function partial_cor() for calculating the local partial correlations
     from a clg-object.
   * Added function ci_test() for conducting a test for conditional 
     independence.
   * Added functionality for normalizing estimated (conditional) densities using
     a MC approximation.
   * Local (partial) correlations now comes with approximate confidence 
     intervals.

# lg 0.1.0 

   * Initial release
