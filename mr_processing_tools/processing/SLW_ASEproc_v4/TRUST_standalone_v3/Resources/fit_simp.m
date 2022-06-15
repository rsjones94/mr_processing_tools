function err = fit_simp(parameter,X,Y)
% This function is called by LSQNONLIN.

fit = parameter(2)*(exp(X*parameter(1)));
err = fit - Y; 