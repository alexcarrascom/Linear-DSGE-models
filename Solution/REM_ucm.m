function [P,Pz,R,Rz,lambda,Omega]=REM_ucm(A,B,C,D,F,G,H,J,K,L,M,N)
%% REM_ucm: 
% Compute matrices which define the equilibrium law of motion of a
% rational expectations model following the undetermined coefficients 
% method (Uhlig, 1997).
%
% The system is represented for three main blocks:
%   0 = A*x(t) + B*x(t-1) + C*y(t)+D*z(t)
%   0=E[F*x(t+1) + G*x(t) + H*x(t-1) + J*y(t+1) +K*y(t) +L*z(t+1) + M*z(t) | t]
%   z(t+1)= N*z(t) + e(t+1)
%
% The equilibrium law of motion are given by:
%   x(t)=P*x(t-1) + Pz*z(t)
%   y(t)=R*x(t-1) + Rz*z(t)
%
% ====================================================================
% Syntaxis: [P,Pz,R,Rz,EIG_,V_1]=REM_ucm(A,B,C,D,F,G,H,J,K,L,M,N)
% ====================================================================
% lambda    : Generalized eigenvalues with modulus less than one.
% Omega     : Corresponding eigenvectors.
% ===================================================================================
% By Alex Carrasco - PUCrio, september 2017
% ===================================================================================

%% [I] Setting up 
[~,mendog]=size(A);
[r_c,nendog]=size(C);
[kexo,~]=size(N);

% I am supossing that matrix C has complete column rank
flag_pseudo=1;
if r_c==nendog  
    [~,ll]=eig(C);
    if all(diag(ll))
        flag_pseudo=0;
    end
elseif r_c < nendog
    error('Not possible to get a solution! Try to redefine variables or to check the calibration.');
end

%% [II] Finding law of motions

if flag_pseudo
    % Here the compution needs to calculate a Moore-Penrose pseudo inverse.
    % So, this block will be called only in strict necesary cases.
    
    % Getting P
    C_inv=pinv(C);
    C_cero=(null(C'))';
    Cero=zeros(r_c-nendog,mendog);

    aux1=[[C_cero*A; (J*C_inv*B)-G+(K*C_inv*A)] ...
        [C_cero*B; (K*C_inv*B)-H]; eye(mendog) zeros(mendog)];
    aux2=[[Cero;F-(J*C_inv*A)] zeros(mendog); zeros(mendog) eye(mendog)];        
       
    [V0,lambda0]=eig(aux1, aux2);
    lambda0=diag(lambda0);
        
    aux3=(abs(lambda0)<1);
    lambda = lambda0(aux3); nstable = numel(lambda);
        
    aux4 = V0*diag(aux3);
    Omega    = aux4(:,~all(aux4==0));
    Omega    = diag(lambda)\Omega(1:nstable,1:nstable);
                
    if nstable==mendog
        P=(Omega*diag(lambda))/Omega;   
        
        % Getting R
        R = -C_inv*(A*P+B);

        % Getting Pz y Rz
        V = [kron(eye(mendog),A) kron(eye(mendog),C); ...
            kron(N',F)+kron(eye(mendog),(F*P+J*R+G)) kron(N',J) + kron(eye(mendog),K)];
           
        Pz_Rz=V\[-vec(D);-vec(L*N+M)];
            
        Pz = reshape(Pz_Rz(1:mendog*kexo),mendog,kexo); 
        Rz = reshape(Pz_Rz(mendog*kexo+1:end),nendog,kexo);

    else
        error('Blanchard-Khan conditions are not satisfied');
    end
    
else
    % Here the computation is faster!
    % Getting P
    aux1=[(J/C)*B-G+(K/C)*A (K/C)*B-H; eye(mendog) zeros(mendog)];
    aux2=[F-(J/C)*A zeros(mendog); zeros(mendog) eye(mendog)];

    [V0,lambda0]=eig(aux1, aux2);
    lambda0=diag(lambda0);
        
    aux3=(abs(lambda0)<1);
    lambda = lambda0(aux3); nstable = numel(lambda);
        
    aux4 = V0*diag(aux3);
    Omega    = aux4(:,~all(aux4==0));
    Omega    = diag(lambda)\Omega(1:nstable,1:nstable);
        
    if nstable==mendog        
        P=Omega*diag(lambda)/Omega;   

        % Getting R
        R=-C\(A*P+B);
            
        % Getting Pz
        aux5 = kron(N',F-(J/C)*A)+kron(eye(kexo),F*P+G+J*R-(K/C)*A);
        aux6 = vec(((J/C)*D-L)*N + (K/C)*D-M);
        Pz =aux5\aux6;            
        Pz = reshape(Pz,mendog,kexo);
            
        % Getting  Rz
        Rz=-C\(A*Pz+D);
    else
        error('Blanchard-Khan conditions are not satisfied')
    end
end
    
end