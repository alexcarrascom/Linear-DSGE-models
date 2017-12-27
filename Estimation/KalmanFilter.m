function [KF,LF,Lt]=KalmanFilter(Y,A1,Ae,Psi2,Se,Su,varargin)
%% KalmanFilter: 
% computes the forecast and smoother for state vector and the
% likelihood function assuming the following system
%       z_t = A1*z_{t-1} + Ae*e_t [transition equation]
%       y_t = Psi1*x_t + Psi2*z_t + u_t [measurement equation]
% where E[e_te'_t]=Se, E[u_tu'_t]=Su, and are serially uncorrelated and
% uncorrelated among them.
% **************************************
%   By Alex Carrasco, november 2017
% **************************************

%% [I] Default and options
smooth_flag = 1;
nz          = size(Ae,1);
[ny,T]      = size(Y);
X           = zeros(1,T);
Psi1        = zeros(ny,1);

% Initialization 
z0 = zeros(nz,1);
if max(abs(eig(A1)))<1
    P0 = reshape((eye(nz^2,nz^2)-kron(A1,A1))\vec(Ae*Se*Ae'),nz,nz);
end

zt_t1 = nan(nz,T);
zt_t = nan(nz,T+1);
zt_T = nan(nz,T);
Pt_t1 = nan(nz,nz,T);
Pt_t  = nan(nz,nz,T+1);
Pt_T  = nan(nz,nz,T);
yt_t1 = nan(ny,T);
Ft_t1 = nan(ny,ny,T);
Lt    = nan(1,T);

for ii=1:numel(varargin)
    if strcmp(varargin{ii},'init_z'), z0=varargin{ii+1}; end
    if strcmp(varargin{ii},'init_P'), P0=varargin{ii+1}; end
    if strcmp(varargin{ii},'exog'),   X =varargin{ii+1}; 
                                      Psi1=varargin{ii+2}; end
    if strcmp(varargin{ii},'nosmooth'), smooth_flag=0; end
end


%% [II] Kalman filter
zt_t(:,1) = z0; Pt_t(:,:,1)=P0; 
% Assuming normality
fy=@(y,Ey,F) -ny/2*log(2*pi) -0.5*log(det(F)) -0.5*(y-Ey)'*(F\(y-Ey)) ;

for t=1:T
 % forecasiting t based on t-1
   zt_t1(:,t)   = A1*zt_t(:,t);        
   Pt_t1(:,:,t) = A1*Pt_t(:,:,t)*A1' + Ae*Se*Ae';
   yt_t1(:,t)   = Psi1*X(:,t) + Psi2*zt_t1(:,t);
   Ft_t1(:,:,t) = Psi2*Pt_t1(:,:,t)*Psi2' + Su;
 % updating inference
   zt_t(:,t+1)    = zt_t1(:,t)+Pt_t1(:,:,t)*Psi2'*(Ft_t1(:,:,t)\(Y(:,t)-yt_t1(:,t))); 
   Pt_t(:,:,t+1)  = Pt_t1(:,:,t)-Pt_t1(:,:,t)*Psi2'*(Ft_t1(:,:,t)\(Psi2*Pt_t1(:,:,t)));
 % log-likelihood function
   Lt(t) = fy(Y(:,t),yt_t1(:,t),Ft_t1(:,:,t));
% Lt(t) = log(mvnpdf(Y(:,t)',yt_t1(:,t)',Ft_t1(:,:,t)));
end

KF.forecastZ.z = zt_t1;
KF.forecastZ.P = Pt_t1;
KF.forecastY.y = yt_t1;
KF.forecastY.F = Ft_t1;
KF.filterZ.z   = zt_t;
KF.filterZ.P   = Pt_t;
LF = sum(Lt);      

%% [III] Smoother
if smooth_flag
    zt_T(:,T) = zt_t(:,T+1); Pt_T(:,:,T) = Pt_t(:,:,T+1);
    for t=T-1:-1:1
        J           = Pt_t(:,:,t+1)*A1'/Pt_t1(:,:,t+1);
        zt_T(:,t)   = zt_t(:,t+1) + J*(zt_T(:,t+1)-zt_t1(:,t+1));
        Pt_T(:,:,t) = Pt_t(:,:,t+1) + J*(Pt_T(:,:,t+1)-Pt_t1(:,:,t+1))*J';
    end
    KF.smoothZ.z   = zt_T;
    KF.smoothZ.P   = Pt_T;
end

end