# GDit: gravity darkening in python

GDit implements the mathematical model for gravity darkening in rotating stars presented by [Espinosa Lara & Rieutord (2011)](http://adsabs.harvard.edu/abs/2011A%26A...533A..43E).
![Spheroids](https://github.com/aarondotter/GDit/blob/master/plots/spheroids.png)
The above figure shows a range of spheroids corresponding to rotating stars at different fractions of the Keplerian angular velocity &omega; and inclination angle i. In the figure the colorscale indicates the behavior of the stellar surface: blue is fainter/cooler and yellow is brighter/hotter.

The code allows one to calculate the gravity darkening coefficients C<sub>L</sub> and C<sub>T</sub>, which tell how both rotation and inclination alter the observed, surface-averaged temperature and luminosity of the star. These are inspired by [Georgy et al. (2014)](http://adsabs.harvard.edu/abs/2014A%26A...566A..21G) but calculated in a manner consistent with [Espinosa Lara & Rieutord (2011)](http://adsabs.harvard.edu/abs/2011A%26A...533A..43E) for Roche equipotential surfaces [MESA V]().

C<sub>T</sub> and C<sub>L</sub> allow one to quickly determine the observational consquences of rotating stellar models with arbitrary orientation.
![H-R diagram](https://github.com/aarondotter/GDit/blob/master/plots/HRD.png)
