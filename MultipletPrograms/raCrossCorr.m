function [C,varargout] = raCrossCorr(S,varargin)
% This function computes the cross-correlation matrix for coral structures
% by treating each column of an MxN coral structure as a 'vector' in a
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
% Each column gives M record sections recorded by an M receiver array.
% There are N separate array observations.
%
% The translation done via cross correlation is (1) computed in the
% frequency domain to optimize efficiency and (2) done only along columns,
% so that there is no changing relative times between stations in a single
% record section. Data are tapered and detrended internally.
%
% USAGES
% [C]           = raCrossCorr(S);
% [C]           = raCrossCorr(S,chan);
% [C]           = raCrossCorr(S,chan,ind);
%              ...
%
% [C,I,dt,G,pt,lags]    = raCrossCorr(...);
%
%
% INPUT:
% S:    An MxN coral structure array
% chan: Optional input to cross-correlate the array observations only over
%       the specificed channel in field staChannel.
% ind:  A vector of indices pointing to which record section to
%       cross-correlate of the form ind = 1:size(S,2), or some subset.
%       Default value for ind is all, with auto-correlations omitted.
%
% OUTPUT:
% C:    A sparse matrix giving the maximum cross-correlation coefficents
%       for each pair of columns.
% I:    The shift index providing the maximum correlation. I(k,n) is amount
%       record section k is shifted with respect to record section n to
%       maximize C between record section k and n.
% dt:   Differential time lags giving maximum correlation, corresponding to
%       index I. May be used in double-difference operations.
% G:    A matrix of correlograms with M*N columns.
% pt:   A (M*N)x2 set of pointers. Row "k" indicates which record sections
%       in G(:,k) were correlated.
% lags: The lag times for the correlograms. The convention is positive lags
%       mean right-shifts, to be consistent with circshift.
%-----------------------------------------------------------------------
% %EXAMPLE 1
% %Suppose you have S and size(S) = [M,N]
%
% %compute array cross correlation and get max indices, called ind
%
% [cc,ind,dt] = raCrossCorr(S);
%
% %shift the coral structure relative to S(:,end) so that it is aligned
% with the array observation in S(:,end)
%
% S        = coralCircShift(S,ind(:,end));
%
% %EXAMPLE 2
% %Suppose you have Sc and size(Sc) = [1,1].
%
% %Make a set of duplicates of Sc, and forward model shifts with indShift
% %by using coralSubShift to shift the signals:
%
% numRep  = 5;
% S       = repmat(Sc(1,1),1,numRep);
% N       = [(S(:).recNumData)];
% N       = median(N);
% nRec    = size(S,1);
% nMult   = size(S,2);
% indShift    = 1 + (N-1)*rand(nMult,1);
% dtShift     = (S(1).recSampInt)*diff(indShift);
% S           = coralSubShift(S,indShift);
%
% %S now contains a set of identical seismograms that are translated by
% %fractional indices indShift and times dtShift.
% %Now take seismogram pairs, and find the time difference between the two.
%
% [C,ind]     = raCrossCorr(S);
%
% %To test your shifts, plot all possible seismograms with their inverted
% %shifts:
% for k = 1:numRep,
%     for n = k+1:numRep,
%         Skn = coralCircShift(S(k),ind(k,n));
%         temp = raGramMat([S(n),Skn]);
%         B(k,n) = temp(1,2);
%         figure; coralPlot(S(n));
%         h=coralPlot(Skn);set(h,'color','r');
%     end;
% end;
%
% %Each plot should show to nearly coincident seismograms. Note that the
% %input shift indices are all positive, but inverted shift indices may be
% %negative. This prevents wrap around effects.
%-----------------------------------------------------------------------
% Joshua D Carmichael
% josh.carmichael@gmail.com
%
% Edit Log
% 18.Nov.2009
% Made function
%
% 12.April.2010
% insert coralExtract in favor of similar code to extract the requested
% substructure when a channel option is input.
%
% 13.July.2010
% Added output of differential travel times, and changed output indices to
% remove wrap around effects. Previously, all output indices were positive.
% Now they are negative if the shift is greater than 1/2 the number of data
% samples. The changes appear below if(nargout>1) and if(nargout>2). I
% confirmed the code with EXAMPLE 2 above as a test.
%
% 04.August.2010
% I considered removing the coralTaper function. I think it is necessary to
% retain for the finiteness of the time series, to reduce spectral leakage.
% It is likely better to retain detrending and tapering for the sake of
% being careful.
% I also removed some redundant operations in the correlation computation.
% I replaced mean with sum, and included output option of correlogram for
% record sections, with lags, and index pointers.
%-----------------------------------------------------------------------

ind1    = cellfun(@ischar,varargin);
ind2    = cellfun(@isnumeric,varargin);

%Choose which channels to check
if(nnz(ind1)>0)
    
    S = coralExtract(S,'staChannel',varargin{ind1});
    
end;

%if given indices were input
if(nnz(ind2)>0)
    
    ind2 = varargin{ind2};
    
else
    
    ind2 = [];
    
end;


% for k = 1:length(varargin)
%
%     if(ischar(varargin{k}))
%         S = coralExtract(S,'staChannel',varargin{k});
%     end;
%
% end;
S=coralPad(S);
N       = [S(:).recNumData];
sInt    = [S(:).recSampInt];

%Check for sample length consistency
if(std(N)>1e3*eps),
    
    error('MATLAB:raCrossCorr','Data columns not same length. Use coralPad.');
    
end;

%Check for sample rate consistency
if(std(sInt)>1e3*eps),
    
    error('MATLAB:raCrossCorr','Data needs to have uniform sample rate. Up/Down sample.');
    
end;

N   = N(1);
sInt= sInt(1);
C	= sparse(size(S,2),1);
I	= C;
S	= arrayfun(@coralDetrend,S);
S	= arrayfun(@coralTaper,S);
G   = [];
pt  = [];

if(isempty(ind2))
    
    for m = 1:size(S,2)
        
        for n = m+1:size(S,2)
            
            Sx      = ifft( conj(fft([S(:,m).data])).*(fft([S(:,n).data])) );
            Cgrm    = sum(Sx,2)./(norm([S(:,m).data],'fro')*norm([S(:,n).data],'fro'));
            [tmp,k]	= max(abs(sum(Cgrm,2)));
            C(m,n)  = tmp;
            I(m,n) 	= k-1; %Shift Index
            
            if(nargout>3)
                
                G       = cat(2,G,Cgrm(:));
                pt      = cat(1,pt,[m,n]);
                
            end;
            
        end;
        
        % Keep old code: To illustrate what direction shifting is done.
        %     X       = circshift([S(:,m).data],k-1); %Shift to the right
        %     Y       = [S(:,n).data];
        %     C(m,n)	= abs(trace(X'*Y))/(norm(X,'fro')*norm(Y,'fro'));
        %     I(m,n) 	= k-1;
        
    end;
    
else
       
    
    for m = 1:size(S,2),
        
        for n = ind2
            
            Sx      = ifft( conj(fft([S(:,m).data])).*(fft([S(:,n).data])) );
            Cgrm    = sum(Sx,2)./(norm([S(:,m).data],'fro')*norm([S(:,n).data],'fro'));
            [tmp,k]	= max(abs(sum(Cgrm,2)));
            C(m,n)  = tmp;
            I(m,n) 	= k-1; %Shift Index
            
            if(nargout>3)
                
                G       = cat(2,G,Cgrm(:));
                pt      = cat(1,pt,[m,n]);
                
            end;
            
        end;
        
    end;
        
end;

%C   = triu(C,1);

if(nargout>1)
    
    %shift indices
    diffSign            = +1*ones(size(I));
    diffSign(I > N/2)   = -1;
    diffSign            = triu(diffSign,+1);
    I(I > N/2)          = N - I(I > N/2);
    I                   = I.*diffSign;
    varargout{1}        = I;
    
end;

if(nargout>2)
    
    %differential shift times, in seconds
    dt              = I.*sInt;
    varargout{2}	= dt;
    
end;

if(nargout>3)
    
    %correlograms
    varargout{3}    = G;
    
end;

if(nargout>4)
    
    %record section pointers
    varargout{4}    = pt;
    
end;

if(nargout>5)
    
    %lag indices for circular shifting. Right shift means positive.
    lags            = 0:size([S(:,m).data],1)-1;
    varargout{5}    = lags(:);
    
end;

return;