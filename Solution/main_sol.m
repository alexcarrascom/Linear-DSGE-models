%% Solutions for linear RE models (Undetermined coefficients and GSchur)
% Student: Alex Carrasco - PUC-rio, 2017
% Suggestion: Run by predetermined sections [Ctrl+Enter]

clc; clear all; close all;
addpath(genpath([pwd '\Solution']));

%% [I] RBC model
% Calibrating
beta = 0.99; sigma=1; alpha=1/4; varphi=1; delta=1;
rho_a = 0.75; sigma_a = 0.01; 

% Steady-state
r=1/beta-1+delta; k_n= (alpha/r)^(1/(1-alpha));
y_n=k_n^alpha;
c_n=y_n - delta*k_n;

% [I.1] Matrices, Coeficientes indeterminados
A=[-k_n/y_n; 0; 0; 0; 0]; B=[k_n/y_n*(1-delta); alpha; 0; -r; 0]; 
C = [1,-c_n/y_n, 0,0,0; -1,0,1-alpha,0,0; 0,-sigma, -varphi,1,0; r,0,0,0,-1; 1,0,-1,-1,0];
D=[0;1;0;0;0];
J=[0 -sigma 0 0 beta]; K=[0 sigma 0 0 0];
%
% [I.2] Finding policy functions and irfs (Uhlig)
[P,Pz,R,Rz]=REM_ucm(A,B,C,D,0,0,0,J,K,0,0,rho_a);
vars={'output', 'consumpion','employment','real wage','real rate','capital','shock'};
Ay = [zeros(5,5) R; zeros(1,5) P]; Bz = [Rz;Pz];
Cz = rho_a;
IR1=IRF(Ay,Bz,Cz,'noprint','periods',16,'var_shocks',sigma_a,'names',vars,'transform',@(x) 100*x);

% print -depsc2 'fig_rbc1.eps';
% eps2pdf('fig_rbc1.eps');

%% [I.3] Matrices, Schur generalizado
A=[k_n/y_n 0 0 0 0 0; zeros(4,6); 0 0 sigma 0 0 -beta]; 
B=[[k_n/y_n*(1-delta); alpha;0;-r;0], C; 0 0 sigma 0 0 0 ];
C=[0; 1; 0; 0;0 ;0];

% [I.4] Finding policy functions and irfs (Uhlig)
[P,Pz,R,Rz]=REM_gschur(A,B,C, sigma_a,1);
vars={'output', 'consumpion','employment','real wage','real rate','capital','shock'};
Ay = [zeros(5,5) R; zeros(1,5) P]; Bz = [Rz;Pz];
Cz = rho_a;
IR2=IRF(Ay,Bz,Cz,'noprint','periods',16,'var_shocks',sigma_a,'names',vars,'transform',@(x) 100*x);

% print -depsc2 'fig_rbc2.eps';
% eps2pdf('fig_rbc2.eps');


%% [II] New Keynesian model
% Calibration
beta = 0.99; theta=3/4; sigma=1; alpha=1/4; varphi=5;
eps = 9; eta=4; phi_pi = 1.5; phi_y=0.5/4; 

rho_a = 0.9; rho_z = 0.5; rho_v = 0.5;
sigma_a = 0.01; sigma_z = -0.005; sigma_v = 0.0025; % Here will be the size of the shock

% Function of paramters
zeta= (varphi + alpha + sigma*(1-alpha))/(1-alpha+alpha*eps);
psi_ya = (1+varphi)/((1-alpha)*sigma+varphi+alpha);
kappa = (1-theta)*(1-theta*beta)*zeta/theta;

% [II] Stationary technology progress

%  [II.1] Matrices
%   A*E[x(t+1)',y(t+1)'|t]'=B*[x(t+1)',y(t+1)']' + C*z(t)
%   z(t+1)= N*z(t) + e(t+1)
%  Here, y=(inflation, output gap)'; z=(a,z,v)'; x=0;

A = [beta 0; 1/sigma 1];
B = [1 -kappa; phi_pi/sigma  1+phi_y/sigma];
C = [0 0 0; psi_ya*((1-rho_a)+phi_y/sigma), -(1-rho_z)/sigma, 1/sigma];
N = diag([rho_a,rho_z,rho_v]); D=diag([sigma_a,sigma_z,sigma_v]);

%  [II.2] Finding law of motion
[~,~,~,Rz]=REM_gschur(A,B,C,N,0); % Since m=0, the relevant matrix is Rz

%  [II.3] Impulse-response functions 
Ay = zeros(2);
Az = Rz;
Cz = N;
[IR,~]=IRF(Ay,Az,Cz,'nograph','periods',16,'var_shocks',D);
EDpi = IR(2:end,1,:); IR=IR(1:end-1,:,:); cIR=cumsum(IR);

%  [II.4] Replicating graphs
vars={'output gap', 'inflation', 'output', 'employment', 'real wage', ...
       'price level', 'nominal rate', 'real rate', 'money supply','shock'};
exo ={'a' 'z' 'v'};

Y=[];
Y(:,1,:)=IR(:,2,:);                                     % output gap
Y(:,2,:)=4*IR(:,1,:);                                   % inflation
Y(:,3,:)=Y(:,1,:)+psi_ya*IR(:,3,:);                       % output    
Y(:,4,:)=(1-alpha)^(-1)*(Y(:,3,:)-IR(:,3,:));           % employment
Y(:,5,:)=sigma*Y(:,3,:)+varphi*Y(:,4,:);                % real wage
Y(:,6,:)=cIR(:,1,:);                                    % price level
Y(:,7,:)=4*(phi_pi*IR(:,1,:)+phi_y*Y(:,3,:)+IR(:,5,:)); % nominal rate
Y(:,8,:)=Y(:,7,:) -4*EDpi;
Y(:,9,:)=Y(:,6,:)+Y(:,3,:)-eta/4*Y(:,7,:);

k=size(Y,3);
m=size(Y,2);

for jj=1:k
   figure(2+jj)
   for ii=1:m
       subplot(5,2,ii)
       plot(100*Y(:,ii,jj),'marker','o','markersize',6,'linewidth',1.5);
       title(vars{ii}); axis tight; 
   end
   subplot(5,2,10)
   plot(100*IR(:,2+jj,jj),'marker','o','markersize',6,'linewidth',1.5);
   title(exo{jj}); axis tight;
%    eval(['print -depsc2 fig_nk1_' num2str(jj) '.eps;']);
%    eps2pdf(['fig_nk1_' num2str(jj) '.eps']);
end




