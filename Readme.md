# CPROP - Converge properties

  Converge properties is a library for thermodynamic and transport properties of fluids to be used
  in the optimization and simulation of engineering design and process analysis.
  It is programed in implicitly statically typed python so it can be compiled to other languages.
  
## Usage

Just download and put the CPROP.py file in the folder ;)

## Examples

~~~~ python
  >>> from CPROP import R134a
  >>> r = R134a()
  >>> r.rho(300,100) # Compute Density [kg/m3] as function of T[K] and P [kPa]
  4.1739312767380143 
  >>> r.rho(300,1000) 
  1198.8373358177078
  >>> r.sig(300) # Compute Superficial Tension as function of T[K]
  7.8332491148274812
  >>> r.eta(300,100) # Compute viscosity [microPa/s] as function of T[K] and P [kPa]
  12.057383150000001
  >>> r.eta(300,1000)
  193.77767982471246
  >>> r.kappa(300,100) #Compute thermal conductivity [W/m2*K]
  0.014550624000000003
  >>> r.kappa(300,1000)
  0.080673664499999978
  >>> help(r) #Show help!
~~~~  
  
## Copyright

  Copyright (c) 19th March 2012 Fabio C. Canesin &lt;canesin@converge-es.com&gt;

  CPROP is free software: you can redistribute it and/or modify it under the
  terms of the GNU General Public License as published by the Free Software
  Foundation, version 3 of the License.  See the website [http://www.gnu.org/licenses/],
  for a description of the GNU General Public License terms under which you can copy the files.