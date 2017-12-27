function ll=invgammapdf(x,a,b)

% **************
% Alex, december 2017
% *******************

ll=b^a/gamma(a).*(x.^(-a-1)).*(exp(-b./x));

end