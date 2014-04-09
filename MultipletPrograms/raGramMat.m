function [C] = raGramMat(S,varargin)
% This function computes the Gram matrix for coral structures by
% treating each column of an MxN coral structure as a 'vector' in a
% Hilbert space whose elements are array seismograms
%
% Put another way, the coral structure under consideration
% has the form:
%
%  [S(1,1)  S(1,2)  ...  ...    S(N,1)]
%  [S(2,1)  S(2,2)  ...  ...    S(N,2)
%    ...     ...    ...  ...     ...
%  [S(M,1)                      S(M,N)];
%
% Each column gives M seismogram records given by an M receiver array.
% There are N separate array observations.
%
% USAGES
% 
% [C]   = raGramMat(S);
% [C]   = raGramMat(S,chan);
%   
% INPUT:
% S:    An MxN coral structure array
% chan: Optional input to project array observations only over the
%       specified channels.
%
% OUTPUT:
% C:    A sparse matrix giving Gram matrix. That is, element i,j is the
%       Frobenius inner product between S(:,i) and S(:,j). Only the upper
%       triangular part is retained, as C is symmetric.
%-----------------------------------------------------------------------
% Latest Edit: 28.Feb.2010
% Joshua D Carmichael
% josh.carmichael@gmail.com
%-----------------------------------------------------------------------

for k = 1:length(varargin)
    if(ischar(varargin{k}))
        S = coralExtract(S,'staChannel',varargin{k});
    end;
end;

C	= sparse(size(S,2),1);
S	= arrayfun(@coralDetrend,S);
S	= arrayfun(@coralTaper,S);

for m = 1:size(S,2)
    
    for n = m+1:size(S,2)
        
        X       = [S(:,m).data];
        Y       = [S(:,n).data];
        C(m,n)	= (trace(X'*Y))/(norm(X,'fro')*norm(Y,'fro'));
        
    end;
    
end;

C   = triu(C,1);

return;