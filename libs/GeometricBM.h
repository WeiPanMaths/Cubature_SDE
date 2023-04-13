#pragma once

// dS = rSdt + sigma S dW
double GeometricBM(double s, double r, double sigma, double t, double bm);

// drift is (r - .5*vol*vol)(T-t)
// bm here is sigma W
double GeometricBM(double s, double drift, double bm);