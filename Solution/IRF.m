function [IR,cIR]=IRF(A,B,C,varargin)
%% IRF:
% Compute impulse-response functions of the varaibles of the system
%   Y(t)   = A*Y(t-1) + B*Z(t)
%   Z(t+1) = C*Z(t)   + D*u(t+1)
%
% D=0.01*I by default. Introducing different variance matrix for
% structural shocks is possible with the option 'var_shocks'.
% =======================================================================
% By Alex Carrasco, September 2017
% =======================================================================

%% [I] Setting up
% Default
flag_graph = 1; flag_print=1;
[nendo,~]  = size(A);
[kexo,~]   = size(C);
nsg        = nendo+kexo;        % Number of subplots per shock
D          = eye(kexo);
shocks     = (1:kexo)';
names      = cellstr([repmat('var ',nsg,1),num2str((1:nsg)')]);   
T          = 20;
f          = @(x) x;

% Options
for i=1:numel(varargin)
   if strcmp(varargin{i},'var_shocks'),  D = varargin{i+1}; shocks=(1:size(D,2))'; end
   if strcmp(varargin{i},'shocks'),      shocks=varargin{i+1}; end
   if strcmp(varargin{i},'names'),       names =varargin{i+1}; end
   if strcmp(varargin{i},'nograph'),     flag_graph=0; end
   if strcmp(varargin{i},'noprint'),     flag_print=0; end
   if strcmp(varargin{i},'periods'),     T=varargin{i+1}; end
   if strcmp(varargin{i},'transform'),   f=varargin{i+1}; end
end
% ... Getting ready
nsh = numel(shocks); 
Z = zeros(kexo,T+1,nsh);
Y = zeros(nendo,T+1,nsh);
U = zeros(nsh,T+1,nsh);
I = eye(nsh);
%% [II] IRF computation
for k=1:nsh
    U(:,2,k)=I(:,shocks(k));
    for t=2:T+1;
        Z(:,t,k) = C*Z(:,t-1,k) + D*U(:,t,k);
        Y(:,t,k) = A*Y(:,t-1,k) + B*Z(:,t,k);
    end
end

IR = permute([Y;Z],[2 1 3]); IR=IR(2:end,:,:);
cIR = cumsum(IR);

%% [III] Graphs
[nfigs,nr,nc,ls] = plotorg(nsg);
if flag_graph
for k=1:nsh
    for fig = 1:nfigs;
        if fig==nfigs
            figure(fig+nfigs*(k-1))
            for i=1:ls
                subplot(nr,nc,i)
                h=plot([f(IR(:,(nfigs-1)*9+i,k)),zeros(T,1)]);
                set(h(1),'color','b','linewidth',1.5,'marker','o','markersize',5);
                set(h(2),'color','k','linewidth',1.25);
                title(names{(nfigs-1)*9+i});
                axis tight;
            end
         else
            figure(fig+nfigs*(k-1))
            for i=1:9
                subplot(3,3,i)
                h=plot([f(IR(:,(fig-1)*9+i,k)),zeros(T,1)]);
                set(h(1),'color','b','linewidth',1.5,'marker','o','markersize',5);
                set(h(2),'color','k','linewidth',1.25);
                title(names{(fig-1)*9+i});
                axis tight;
            end
        end
        if flag_print
            eval(['print -depsc2 graph' num2str(fig) 'shock_' num2str(k) '.eps;']);
        end
    end
end
end

end

function [nfigs,nr,nc,ls]=plotorg(nsg)
    
nfigs  = ceil(nsg/9); ls = rem(nsg,9);

if     ls == 1
    nr=1; nc=1;
elseif ls == 2;
    nr=2; nc=1;
elseif ls == 3
    nr=3; nc=1;
elseif ls == 4
    nr=2; nc=2;
elseif (ls == 5 || ls == 6)
    nr=3; nc=2;
elseif (ls == 7 || ls==8 || ls==0)
    nr=3; nc=3;    
end

end