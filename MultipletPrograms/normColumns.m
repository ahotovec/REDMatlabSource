function  [Wout,varargout]=normColumns(W)
% Treats input matrix columns as time series and normalizes each column,
% without loops. Outputs norm as well.
%
% USAGES: 
% [Wout] = normColumns(W);
% [W,n]  = normColumns(W);
%
% INPUT
% W: a matrix of data, arranged column-wise.
%
% OUTPUT:
% Wout: The normalized data
% n:    The norm of each column
%
% Example;
% X = peaks(48,48);
% W=detrend(X);
% [W]=normColumns(W);
% plotXmatrix(W); title('Normalized Data');
%
%-----------------------------------------------------------------------
% Latest Edit: 08.April.2010
% Joshua D Carmichael
% josh.carmichael@gmail.com
%-----------------------------------------------------------------------

%initialize output norm matrix
Wout            = W;
Wout(isnan(W))  = 0;

n               = sqrt( sum( W.*conj(W),1) );

if(mean(n) < eps)
    
    Wout    = W;
    if(nargout==2)
        varargout{1}=n;
    end;
    return;
    
end;

m               = n(n>eps);
%temp            = diag(m.^(-1));
temp            = repmat(m.^(-1),size(W,1),1);
Wout(:,n>eps)   = W(:,n>eps).*temp;
%Wout(:,n>eps)   = W(:,n>eps)*temp;

if(nargout==2)
    varargout{1}=n;
end;
return;