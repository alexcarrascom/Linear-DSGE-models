function [P,Pz,R,Rz,lambda]=REM_gschur(A,B,C,N,m)
%% REM_gschur: 
% Compute matrices which define the equilibrium law of motion of a
% rational expectations model following Klein's method  (Klein, 1997).
%
% The system is represented for three main blocks:
%   A*E[x(t+1)',y(t+1)'|t]'=B*[x(t+1)',y(t+1)']' + C*z(t)
%   z(t+1)= N*z(t) + e(t+1)
%
% where x(t) contains 'm' endogenous state variable. 
% The equilibrium law of motion are given by:
%   x(t)=P*x(t-1) + Pz*z(t)
%   y(t)=R*x(t-1) + Rz*z(t)
%
% ====================================================================
% Syntaxis: [P,Pz,R,Rz,lambda]=REM_gschur(A,B,C,N,m)
% ====================================================================
% lambda: Generalized eigenvalues with modulus less than one.
% ====================================================================
% By Alex Carrasco - PUCrio, september 2017
% ======================================================================

%% [I] Setting up
[mn,~]=size(A);
[~,k] = size(C);
n = mn-m;

[S,T,Q,Z]=qz(A,B,'real');
lambda=diag(T)./diag(S);
stable=(abs(lambda)<1);
if sum(stable)~=m
    error('Blanchard-Khan conditions are not satisfied');
end
[S,T,Q,Z]=ordqz(S,T,Q,Z,'udo'); % Ordering GSD (stable eigenvalues come first)

%% [II] Computing law of motion

Z_11 = Z(1:m,1:m); Z_12 = Z(1:m,m+1:end); Z_21 = Z(m+1:end,1:m); Z_22 = Z(m+1:end,m+1:end);
S_11 = S(1:m,1:m); S_12 = S(1:m,m+1:end); S_22 = S(m+1:end,m+1:end);
T_11 = T(1:m,1:m); T_12 = T(1:m,m+1:end); T_22 = T(m+1:end,m+1:end);
Q_1  = Q(1:m,:); Q_2 = Q(m+1:end,:);
M=reshape((kron(N',S_22) - kron(eye(k),T_22))\vec(Q_2*C),n,k);

P  = (Z_11/S_11)*(T_11/Z_11);
Pz = -P*Z_12*M + (Z_11/S_11)*(T_12*M-S_12*M*N+Q_1*C)+Z_12*M*N;
R  = Z_21/Z_11;
Rz = (Z_22-(Z_21/Z_11)*Z_12)*M;


end