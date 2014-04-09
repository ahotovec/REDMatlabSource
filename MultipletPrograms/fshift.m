function y = fshift(x,s)
% FSHIFT Fractional circular shift
%
% USAGE
% y = fshift(x,s)
%
% INPUT
% x: Input time series to be shifted. If X is a matrix, shifting is done by
%    columns.
% s: The subsample (or sample) temporal shift
%
% OUTPUT
% y: The shifted time series
%
% FSHIFT circularly shifts the elements of vector x by a (possibly
% non-integer) number of elements s. FSHIFT works by applying a linear
% phase in the spectrum domain and is equivalent to CIRCSHIFT for integer
% values of argument s (to machine precision).
%
% (c) 2005 Francois Bouffard
%          fbouffar@gel.ulaval.ca

needtr = 0; 

if size(x,1) == 1; 
    x = x(:); 
    needtr = 1; 
end;

N = size(x,1); 
r = floor(N/2)+1; 
f = ((1:N)-r)/(N/2); 

%make 1i positive to ensure it's consistent with circshift.
p = exp(1i*s*pi*f)'; 

%if x is a matrix
p = repmat(p,1,size(x,2));

y = ifft(fft(x).*ifftshift(p)); 

if isreal(x); 
    y = real(y); 
end;

if needtr; y = y.'; end;
return;