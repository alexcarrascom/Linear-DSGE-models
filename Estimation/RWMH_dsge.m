function [theta,thetamode,Sigma,theta1,aratio,LIK]=RWMH_dsge(Y,X,scal,theta0,Sigma0,prior,measureq,simsform,varargin)
% Sampler that uses Random Walk Metropolis-Hastings algorithm 
% Y:       data
% X:       predetermined observable variable
% scal:    scaling parameter for covariance of the proposal distribution.
% theta0:  serves to initialize csminwel function
% prior:  function gives prior pdf
% measureq:  construct matrices for measurement eq.
% simsform:  construct matrices for transition eq.
% ====================================================================
% Obs: Uses Sims' codes to solve the linear model (gensys.m) and for
% numeric optimization (csminwel.m).
% *************************************************
%     By Alex Carrasco, december 2017
% *************************************************

%% [I] set up
burnin = 1/2;
N      = 5000;
scal0  = 2*scal; % following dynare
obj_name = 'objective'; 
compute_mode = 1;
for i=1:numel(varargin)
    if strcmp('n_drop',varargin{i}), burnin=varargin{i+1}; end
    if strcmp('n_replic',varargin{i}), N=varargin{i+1}; end
    if strcmp('init_scal',varargin{i}), scal0=varargin{i+1}; end
    if strcmp('fcn_name',varargin{i}), obj_name=varargin{i+1}; end
    if strcmp('mode_file',varargin{i}), mode_file=varargin{i+1}; compute_mode=0; end
end
nparam = numel(theta0);

%% [II] Finding posterior mode and hessian (% init hessian similar to dynare)
if compute_mode    
    [fh,thetamode,~,H,itct,fcount,retcodeh] = csminwel(obj_name,theta0,Sigma0,[],10e-6,200);
    Sigma = nhess(obj_name,thetamode);
    Sigma = inv(Sigma);     % new covariance for initialization
else
    load(mode_file);
end
%% [III] RWMH
aratio = nan(N,1);
theta = nan(nparam,N);
% scal0=scal0.*ones(1,nparam);
% scal=scal.*ones(1,nparam);

% initializing
% error  = abs()
theta1 = mvnrnd(thetamode,scal0*Sigma);
lik    = dsge_liki(Y,X,theta1,measureq,simsform);
p      = prior(theta1);
lpost0 = p + lik; 
theta(:,1) = theta1;
aratio(1)=1;
LIK(1)=lpost0;

h=waitbar(0,'acceptance ratio = 0.0000');
for i=2:N
    % Random walk MH sampler
    theta1 = mvnrnd(theta(:,i-1),scal*Sigma)'; 
    % solving
    lik    = dsge_liki(Y,X,theta1,measureq,simsform);
    p      = prior(theta1);
    lpost1 = p+lik;    

    if lpost1<=-1e12
        theta(:,i)  = theta(:,i-1);
        aratio(i) = 0;
        LIK(i)=LIK(i-1);
    else        % rejection or acceptance
        alpha = min([1,exp(lpost1-lpost0)]);
        u = rand(1);
        if u<= alpha
            aratio(i) =  1; 
            theta(:,i)= theta1; 
            lpost0=lpost1;
            LIK(i)=lpost1;
        else
            theta(:,i)  = theta(:,i-1);
            aratio(i) = 0;
            LIK(i)=LIK(i-1);            
        end
    end
    waitbar(i/N,h,['acceptance ratio = ', num2str(mean(aratio(~isnan(aratio))))]);
end
close(h);   

% Burn in
burnin = ceil(burnin*N);
theta(:,1:burnin)=[];
aratio(1:burnin)=[];

theta=theta';
% save mode_comp thetamode Sigma;

end