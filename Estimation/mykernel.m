function [f,F,Xgrid]=mykernel(X,varargin)

n_grids  = 100;
kernel = 1;
[N,k]  = size(X);
h      = (4/((k+2)*N))^(1/(k+4))*std(X); % Scott's rule of thumb
minX   = min(X)-3*h; maxX=max(X)+3*h;

for w=1:numel(varargin)
    if strcmp(varargin{w},'kernel'),    kernel=varargin{w+1}; end
    if strcmp(varargin{w},'n_grid'),    n_grids=varargin{w+1}; end
    if strcmp(varargin{w},'bandwidth'), h=varargin{w+1}; end
end

Xgrid  = ones(n_grids,1)*minX + diag(0:(n_grids-1))*ones(n_grids,1)*(maxX-minX)/(n_grids-1); 

dom='X1';  gg  = 'Xgrid(:,1)';
for m=2:k
    dom  = [dom, [' X' num2str(m) ] ];
    gg   = [gg, [', Xgrid(:,' num2str(m) ')'] ];
end
eval( ['[' dom ']=ndgrid(' gg ');'] );

switch kernel
   case 1 %'Uniform'
        K  = @(psi) (abs(psi)<=1)*1/2;
        IK = @(psi) (abs(psi)<=1).*(psi+1)/2 + (psi>1);
    case 2 %'Gaussian'
        K  = @normpdf;
        IK = @normcdf;
    case 3 %'Epanechnikov'
        K  = @(psi) (abs(psi)<=1).*(3/4*(1-psi.^2));
        IK = @(psi) (abs(psi)<=1).*(0.5+3/4*psi -1/4*psi.^3) + (psi>1); 
    case 4 %'Triangular'
        K  = @(psi) (abs(psi)<=1).*(1-abs(psi));
        IK = @(psi) (psi>=-1 && psi<0).*((1-abs(psi))*(psi+1)/2) ...
                    +(psi>=0 && psi<1).*(1-(1-abs(psi))*(-psi+1)/2) ...
                    + (psi>=1);
end

gridsize=size(X1);
f = zeros(gridsize);
F = zeros(gridsize);

for i=1:N
    aux1=1;
    aux2=1;
    for j=1:k
        eval(['aux3 = (X(i,j)-X' num2str(j) ')/h(j);']);
        aux1 = aux1.*K(aux3);     
        aux2 = aux2.*IK(-aux3);
    end
    f=f+aux1;
    F=F+aux2;
end

f=f/(N*prod(h));
F=F/N;

end