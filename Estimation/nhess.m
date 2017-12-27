function [hssn,varargout]=nhess(fcn,para,varargin);

if nargout>2
	error('maximum number of output arguments is 2.');
end

% printing format
fmt='%5.0f ';		

npara = length(para);

% Compute Hessian, element by element, fine tune with dxscale
ndx = 6;
h0  = exp(-(6:2:(6+2*(ndx-1))));	% step size

hssn     = zeros(npara,npara);
hessdiag = zeros(ndx,1);

dxscale  = ones(npara,1);			% specify different scales across parameters
dxscale  = sparse(1:npara,1:npara,dxscale,npara,npara);

fx    = feval(fcn,para,varargin{:});	% evaluate function 

disp(' ');
disp(' ');
disp('Computing Hessian..');
disp('--------------------------------------------------');
disp(sprintf(['  diagonal elements     :   ' fmt ],npara));


% Compute diagonal elements first
for seli=1:npara

	h = dxscale(:,seli)*h0;
	
	for i=1:ndx

		% forward point
		paradx = para + h(:,i);

		% backward point
		parady = para - h(:,i);
		
		% evaluate function at forward and backward points
		fdx   = feval(fcn,paradx,varargin{:});
		fdy   = feval(fcn,parady,varargin{:});

		% Hessian
		hessdiag(i) = -(2*fx-fdx-fdy)/(h(seli,i))^2; 

    end

    hssn(seli,seli) = 0.5*(hessdiag(3)+hessdiag(4));
	disp(sprintf(['                            ' fmt ],seli));
end


% Now compute off-diagonal elements
% Make sure that correlations are between -1 and 1
% errorij contains the index of elements that are invalid
disp(' ');
disp(sprintf(['  off-diagonal elements :   ' fmt ],npara*(npara-1)/2));

errorij = [];
k=1;

for seli=1:npara

	hi = dxscale(:,seli)*h0;
	
	for selj=(seli+1):npara

		hj = dxscale(:,selj)*h0;

		for i=1:ndx

			% forward to seli-th direction
			paradx = para + hi(:,i);

			% backward to selj-th direction
			parady = para - hj(:,i);

			% forward to seli-th direction and backward to selj-th direction
			paradxdy = paradx - hj(:,i);

			% evalutate functions
			fdx   = feval(fcn,paradx,varargin{:});
			fdy   = feval(fcn,parady,varargin{:});
			fdxdy = feval(fcn,paradxdy,varargin{:});

			hessdiag(i) = -(fx-fdx-fdy+fdxdy)/(hi(seli,i)*hj(selj,i));
		end
		
		hssn(seli,selj) = 0.5*(hessdiag(3)+hessdiag(4));
     
		% calculate correlation
		if (hssn(seli,seli)==0)|(hssn(selj,selj)==0)	% 0 when the variances are zero
			corrij = 0;
		else											
			corrij = hssn(seli,selj)/sqrt(hssn(seli,seli)*hssn(selj,selj));
		end

		if (abs(corrij)>0.98)		% abs(corr) is too big, error
			hssn(seli,selj)=0.9*sqrt(hssn(seli,seli)*hssn(selj,selj));
			errorij = [errorij; seli selj corrij ];
		elseif (abs(corrij)<0.005)	% abs(corr) is too small, make it 0
			hssn(seli,selj)=0;
		end

		hssn(selj,seli) = hssn(seli,selj);

		if mod(k,5)==0;
			disp(sprintf(['                            ' fmt ],k));
		end
		k=k+1;

	end
end

hssn = real(hssn);

varargout = {errorij};
