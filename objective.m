function [logpost]=objective(theta)
% example: estimation two paramters
% =====================================================0
% Y:      data, endogenous variables
% X:      predetermined observable variables
% prior:  function gives prior pdf
% measureq:  construct matrices for measurement eq.
% simsform:  construct matrices for transition eq.
% *************
% Alex, december 2017
% **************
global Y X

p      = priors(theta);
likeli = dsge_liki(Y,X,theta,@measure_eq_basic,@sims_basic_model);

logpost = -(likeli+p);

end