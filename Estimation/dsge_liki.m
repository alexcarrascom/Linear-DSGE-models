function [ll,eu]=dsge_liki(Y,X,theta,measureq,simsform)
% Y:    data
% X:    predetermined observable process
% theta:     parameters
% measureq:  construct matrices for measurement eq.
% simsform:  construct matrices for transition eq.
% *********
% option,  to especify exogenous process
% ********

%% solving and obtaining likelihood
% measurement equation
[Psi2,Su,Psi1] = measureq(theta);   
% transition equation
[Gamma0,Gamma1,Const,Psi,Pi,Se]  = simsform(theta); 
[A1,~,Ae,~,~,~,~,eu] = gensys(Gamma0,Gamma1,Const,Psi,Pi);


% likelihood
if ~(sum(eu)<2)    % kalman filter is initialized at eye(statenumber)
    [~,ll]=KalmanFilter(Y,A1,Ae,Psi2,Se,Su,'nosmooth','exog',X,Psi1,'init_P',10*eye(size(A1,1)));
%        ll=kalman(zeros(size(Y,1),1),Psi2,Su,Ae,eye(size(Ae,2)),A1,Y');
%         ll=sum(ll);
else
    ll=-1e12;   % probability zero
end

end