# rdlocrand

The **rdlocrand** package provides **Stata** and **R** implementations of statistical inference and graphical procedures for Regression Discontinuity designs employing local randomization methods. It provides point estimators, confidence intervals estimators, windows selectors, automatic plots, sensitivity analysis and other related features. This work is supported by the National Science Foundation through grant `SES-1357561`.

For technical, methodological and implementation details see:

* Cattaneo, Frandsen and Titiunik (2015): ***Randomization Inference in the Regression Discontinuity Design: An Application to Party Advantages in the U.S. Senate***, Journal of Causal Inference 3(1): 1-24.
* Cattaneo, Titiunik and Vazquez-Bare (2016): ***Inference in Regression Discontinuity Designs under Local Randomization***, Stata Journal 16(2): 331-367.
* Cattaneo, Titiunik and Vazquez-Bare (2017): ***Comparing Inference Approaches for RD Designs: A Reexamination of the Effect of Head Start on Child Mortality***, Journal of Policy Analysis and Management 36(3): 643-681, Summer 2017. [Replication Files]


**Implementation in Stata**

* To install/update in Stata type:

    	net install rdlocrand, from(https://sites.google.com/site/rdpackages/rdlocrand/stata) replace  
  or  

		github install iphone7725/rdlocrand  

* Help files: [rdrandinf](), [rdwinselect](), [rdsensitivity](), [rdrbounds]() -- Replication files: [do-file](), [data-senate]()  
* Repository for manual installation: [https://sites.google.com/site/rdpackages/rdlocrand/stata]()


**Implementation in R:**

* To install/update in R type:

    	install.packages('rdlocrand')

* [Manual]() -- Replication files: [Illustration](), [R-script](), [data-senate]()  
* [CRAN repository]()


**Last update:** March 13, 2018.  
