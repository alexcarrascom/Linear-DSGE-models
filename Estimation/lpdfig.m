function ldens = lpdfig(x, a, b)
% log inverse gamma
ldens = log(2) - lngam(b/2) + (b/2)*log(b*a^2/2) ...
         - ( (b+1)/2 )*log(x.^2) - b*a^2./(2*x.^2);
     