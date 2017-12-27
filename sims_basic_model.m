function [Gamma0,Gamma1,Const,Psi,Pi,Se]=sims_basic_model(x)
% There are 26 parameters
% Law of motion
beta   = 0.99;
h      = x(1);
sigma  = x(2);
varphi = x(3);
delta  = x(4);
gamma  = x(5);
theta  = x(6);
lambda = x(7);
chi    = 1e-3;
rho_i  = x(8);
phi_pi = x(9);
phi_s  = x(10);
phi_c  = x(11);
phi_Dc  = x(12);

% exogenous process
rho_pi    = x(13);
rho_i_f   = x(14);
rho_omega = x(15);
rho_theta = x(16);
rho_phi   = x(17);
rho_a     = x(18);
rho_v     = x(19);

sigma_pi   = x(20);
sigma_i_f   = x(21);
sigma_omega = x(22);
sigma_theta = x(23);
sigma_phi   = x(24);
sigma_a     = x(25);
sigma_v     = x(26);

% function of parameters
Zeta   = (varphi*(1-h)+sigma)*(1-delta)/((1-h)*(1+varphi*delta));
Lambda = lambda*(1-(1-lambda)*beta)/(1-lambda);
rho_y  = sigma*h/(varphi*(1-h)+sigma);
Theta  = Lambda/(theta-1);

Const = 0;
Se    = eye(7); 

% matrices
Gamma0 = zeros(21,21);
Gamma1 = zeros(21,21);
Psi    = zeros(21,7);
Pi     = zeros(21,4);

% standard deviations
Psi(7:13,:) = diag([sigma_pi,sigma_i_f,sigma_omega,sigma_theta,sigma_phi,sigma_a,sigma_v]);

% ******************************
%   Eq1: Eq Euler
% y pi i mc rer Ds Ey+1 Epi+1 y-1 y-2 y-3 y-4 i-1 Ds-1 pif if omega theta phi a v ----> 21 variables
% ******************************
Gamma0(1,1)=-(1+h); Gamma0(1,[3 17])=-(1-h)/sigma;
Gamma0(1,7)=1;    Gamma0(1,8)=(1-h)/sigma;
Gamma1(1,1)=-h;

% ******************************
%   Eq2: Costos marginales
% y pi i mc rer Ds Ey+1 Epi+1 y-1 y-2 y-3 y-4 i-1 Ds-1 pif if omega theta phi a v ----> 21 variables
% ******************************
Gamma0(2,1)=-Zeta; Gamma0(2,[4 5])=[1 -(1+varphi)*delta/(1+varphi*delta)];
Gamma0(2,20)= (1+varphi)/(1+varphi*delta);
Gamma1(2,1)=-Zeta*rho_y;

% ******************************
%   Eq3: Phillips
% y pi i mc rer Ds Ey+1 Epi+1 y-1 y-2 y-3 y-4 i-1 Ds-1 pif if omega theta phi a v ----> 21 variables
% ******************************
Gamma0(3,2) = (1+beta*gamma); Gamma0(3,4)= -Lambda;
Gamma0(3,8) = -beta; Gamma0(3,18)= Theta;
Gamma1(3,2) = gamma;

% ******************************
%   Eq4: UIP
% y pi i mc rer Ds Ey+1 Epi+1 y-1 y-2 y-3 y-4 i-1 Ds-1 pif if omega theta phi a v ----> 21 variables
% ******************************
Gamma0(4,6) = 1; 
Gamma1(4,[3 5 16 19]) = [1 chi -1 -1]; 
Pi(4,3) = 1;

% ******************************
%   Eq5: Taylor
% y pi i mc rer Ds Ey+1 Epi+1 y-1 y-2 y-3 y-4 i-1 Ds-1 pif if omega theta phi a v ----> 21 variables
% ******************************
Gamma0(5,8) = (1-rho_i)*phi_pi; 
Gamma1(5,[3 13]) = [1 -rho_i]; Gamma1(5,1) = -(1-rho_i)*(phi_c+phi_Dc); 
Gamma1(5,[2 8]) = -(1-rho_i)*phi_pi;
Gamma1(5,12) = (1-rho_i)*phi_Dc;  Gamma1(5,21) = -1; 
Gamma1(5,[6 14])=-(1-rho_i)*phi_s;
Pi(5,4) = (1-rho_i)*phi_pi;

% ******************************
%   Eq6: RER
%  y pi i mc rer Ds Ey+1 Epi+1 y-1 y-2 y-3 y-4 i-1 Ds-1 pif if omega theta phi a v ----> 21 variables
% ******************************
Gamma0(6,[5 6 15 2]) = [1 -1 -1 1]; 
Gamma1(6,5) = 1;  

% *****************
%   Exogenous (7)
% *****************
Gamma0(7:13,15:end) = eye(7);
Gamma1(7:13,15:end) = diag([rho_pi rho_i_f rho_omega rho_theta rho_phi rho_a rho_v]);

% ********************
%   Auxiliar variables (6)
% ********************
Gamma0(14,1)=1; Gamma0(15,2)=1;
Gamma0(16:21,9:14)=eye(6);
Gamma1(14,7) = 1; Gamma1(15,8) = 1;
Gamma1(16,1) = 1;  Gamma1([17 18 19],[9 10 11]) = eye(3); 
Gamma1(20,3) = 1; Gamma1(21,6)=1;

Pi(14,1)=1;
Pi(15,2)=1;




end