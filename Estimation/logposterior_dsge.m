function [logpost]=logposterior_dsge(Y,X,theta,prior,measureq,simsform)

% Y:      data, endogenous variables
% X:      predetermined observable variables
% prior:  function gives prior pdf
% measureq:  construct matrices for measurement eq.
% simsform:  construct matrices for transition eq.
% *************
% Alex, december 2017
% **************

p      = prior(theta);
likeli = dsge_liki(Y,X,theta,measureq,simsform);

logpost = -(likeli+p);

end