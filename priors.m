function [logprior,P]=priors(x)
% Dynare
% prior for estimation of two parameters
% h     ,  BETA_PDF, 0.7, 0.1;         % Castillo et al 2011
% varphi,  INV_GAMMA_PDF, 1,0.3;       % Castillo et al 2011
% sigma,   GAMMA_PDF, 2, .5;           % Herbs Schorfeihe 2015
% delta,   UNIFORM_PDF, , , 0,1;       % unknown
% gamma,	 BETA_PDF, 0.5, 0.1;         % Castillo et al 2011
% lambda,  BETA_PDF, 0.66, 0.1;        % Castillo et al 2011
% thetass, NORMAL_PDF, 6, 1;           % mean based in castillo
% // taylor
% phi_pi, NORMAL_PDF, 1.5, .3;         % Castillo et al 200   
% phi_c,  BETA_PDF, .5, .2;            % OWN
% phi_Dc, BETA_PDF, .5, .2;
% phi_s,  NORMMAL_PDF, .5, .1;            % Casillo et al 
% 
% // rho's (common)
% rho_i     , BETA_PDF, .5,.2;
% rho_pi    , BETA_PDF, .5,.2;
% rho_i_f   , BETA_PDF, .5,.2;
% rho_omega , BETA_PDF, .5,.2;
% rho_theta , BETA_PDF, .5,.2;
% rho_phi   , BETA_PDF, .5,.2;
% rho_a     , BETA_PDF, .5,.2;
% rho_v     , BETA_PDF, .5,.2;
% 
% // desv. 
% sigma_pi   , INV_GAMMA_PDF, .4, .3;
% sigma_i_f  , INV_GAMMA_PDF, .4, .3;
% sigma_omega , INV_GAMMA_PDF, .4, .3;
% sigma_theta , INV_GAMMA_PDF, .4, .3;
% sigma_phi   , INV_GAMMA_PDF, .4, .3;
% sigma_a     , INV_GAMMA_PDF, .4, .3;
% sigma_v     , INV_GAMMA_PDF, .4, .3;
%
%   Own code
% h      = x(1); sigma  = x(2);  varphi = x(3);  delta  = x(4);
% gamma  = x(5);  theta  = x(6); lambda = x(7); rho_i  = x(8);
% phi_pi = x(9);  phi_s  = x(10); phi_c  = x(11);  phi_Dc  = x(12);
% 
% % exogenous process
% rho_pi    = x(13); rho_i_f   = x(14); rho_omega = x(15);
% rho_theta = x(16);  rho_phi   = x(17);  rho_a     = x(18);
% rho_v     = x(19);
%
% sigma_pi   = x(20); sigma_i_f   = x(21); sigma_omega = x(22);
% sigma_theta = x(23);  sigma_phi   = x(24); sigma_a     = x(25);
% sigma_v     = x(26);

P=nan(numel(x),1);

%% Uniform distribution
P(4)=1*(x(4)<1 && x(4)>0);   % a=0,b=1;

%% Beta distribution
param=[1 5 7 8 11 12 13:19];
med=[0.7 0.5 0.66 0.5*ones(1,10)];
sd = [0.1 0.1 0.1 0.2*ones(1,10)];

% alpha= mu^2*(1-mu)/sigma^2 - mu
% beta = (1/mu-1)*alpha

for ii=1:numel(param)
    alpha = med(ii)^2*(1-med(ii))/sd(ii)^2 - med(ii);
    bet  = (1/med(ii)-1)*alpha;
    P(param(ii)) = betapdf(x(param(ii)),alpha,bet);
end

%% Gamma distribution
b = .5^2/2;
a = 2/b;

P(2)=gampdf(x(2),a,b);

%% Normal pdf
param=[6 9 10];
med=[6 1.5 0.5 ];
sd = [1 0.3 0.1];
for ii=1:numel(param)
    P(param(ii))=normpdf(x(param(ii)),med(ii),sd(ii))*(x(param(ii))>0);
end

%% Inverse gamma pdf
param=[3 20:26];
med=[1 .4*ones(1,7)];
sd = [.3 0.3*ones(1,7)];
for ii=1:numel(param)
    if x(param(ii))>0    % bounds (natural)
        P(param(ii))= exp(lpdfig(x(param(ii)),med(ii),sd(ii)));        
    else
        P(param(ii))= 0;
    end
end

if x(23)>1          % a bound for sigma theta
    P(23)=0;
end

if x(3)>2          % a bound for varphi
    P(3)=0;
end
%% end
logprior=log(prod(P));

end