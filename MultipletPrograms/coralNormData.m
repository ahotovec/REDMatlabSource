function  [Snorm,varargout]=coralNormData(S)
% coralNormData: normalizes each seismogram, returns norm if requested. If
% S is an Mx1 coral structure, the vector in each of the M 'data' fields is
% normalized using a dot product.
% If S is an MxN coral structure, the normalization is computed on each of
% the columns of the coral structure. The vectors in each of the M 'data'
% fields in a column is used to form a data matrix.  The size of the  data
% matrix from column k is is S(:,k).recNumData x M, where each row must
% have the same recNumData value. The norm of the entire column is then
% computed using the Frobenius inner product on the data matrix.
%
% USAGES
% [Snorm]   = coralNormData(S);
% [Snorm,N] = coralNormData(S);
%
% INPUT
% S: coral structure array.  Required fields: data, recNumData
%
% OUTPUT:
% Snorm:    A coral structure whose data is normalized
% N:        A matrix giving the norm of the data in each column
%
%-----------------------------------------------------------------------
% EXAMPLE
% %Suppose you have an MxN coral structure array, S
% S = arrayfun( @coralDetrend, S); S = arrayfun( @coralPad, S );
% [Snorm,N]=coralNormData(S);
% plot(N); title('Energy rms in seismogram');
%
% %Suppose you would like to normalize the data to unit norm for each
% %receiver, irrespective of size of S
%
% S = coralDetrend(S);
% S = arrayfun( @coralNormData, S);
%-----------------------------------------------------------------------
% Latest Edit: 09.Dec.2009
% Joshua D Carmichael
% josh.carmichael@gmail.com
%-----------------------------------------------------------------------

%initialize output norm matrix

S = coralPad(S);
N       = zeros(size(S,1),size(S,2));
Snorm   = S;
L       = S(1).recNumData;

M       = [S.data];

M(isnan(M)) = 0;
M(isinf(M)) = 0;

% If S is an Mx1 coral structure, normalize each seismogram using a
% Euclidean inner product, or dot product
if(min(size(S))==1)
    
    for k = 1:size(M,2)
        
        [temp1,temp2]   = normColumns(M(:,k));
        M(:,k)          = temp1;
        N(k)            = temp2';
        Snorm(k).data   = M(:,k);
        
    end;
end;

% If S is an MxN coral structure, treat the N-columns as elements of a
% product space, and normalize the data in each Mx1 structure array using a
% Frobenius inner product
if(size(S,2)>1)
    M       = reshape(M,L,size(S,1),size(S,2));
    
    N       = zeros(size(S,2),1);
    
    for k = 1:size(M,3)
        
        N(k)            = norm( (M(:,:,k)), 'fro');
        
        if(N(k)<eps),
            
            M = 0;
            
        else
            
            M(:,:,k)        = M(:,:,k)./N(k);
            
        end;
        
        for n = 1:size(S,1)
            
            Snorm(n,k).data = M(:,n,k);
            
        end;
        
    end;
end;

if(nargout==2)
    varargout{1}=N;
end;

return;