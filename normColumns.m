function  [Wout,varargout]=normColumns(W,varargin)
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
isnum   = cellfun(@isscalar,varargin);

if(nnz(isnum)>1/2),
    neps = varargin{isnum};
else
    neps            = 1*eps;
end;

%initialize output norm matrix
Wout            = W;
Wout(isnan(W))  = 0;

n               = sqrt( sum( W.*conj(W),1) );

if(mean(n) < neps)
    
    Wout    = W;
    if(nargout==2)
        varargout{1}=n;
    end;
    return;
    
end;

m               = n(n>neps);
%temp            = diag(m.^(-1));
temp            = repmat(m.^(-1),size(W,1),1);
Wout(:,n>neps)  = W(:,n>neps).*temp;
%Wout(:,n>eps)   = W(:,n>eps)*temp;

if(nargout==2)
    varargout{1}=n;
end;
return;