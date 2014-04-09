function [C,varargout] = subRaCrossCorr(S,varargin)
% This function computes the correlograms between all columns of the MxN
% input coral structure, including subsample shifts and stores each
% correlogram in an output cell. This function is intended to be used with
% subsample shifts. To get a sample-resolution shift first, see EXAMPLE
% contained herein.
% As with raCrossCorr each column of the MxN coral structure S is treated 
% as a 'vector' in a Hilbert space whose elements are array seismograms. A
% default set of subsample shifts (-0.95:+0.95)/10 are tested if no inputs
% are provided for test indices.
%
% As with raCrossCorr.m, the coral structure under consideration has form:
%
%  [S(1,1)  S(1,2)  ...  ...    S(N,1)]
%  [S(2,1)  S(2,2)  ...  ...    S(N,2)
%    ...     ...    ...  ...     ...
%  [S(M,1)                      S(M,N)];
%
% Each column gives M seismogram records given by an M receiver array.
% There are N separate array observations.
%
% The translation done via cross correlation is (1) computed in the frequency
% domain to optimize efficiency and (2) done only along columns, so that
% there is no changing relative times between stations in a single record
% section.
%
% USAGES
% 
% [C]       = subRaCrossCorr(S);
% [C]       = subRaCrossCorr(S,ind);
% [C]       = subRaCrossCorr(S,chan);
% [C]       = subRaCrossCorr(S,chan,ind);
% [C,I,dt]  = subRaCrossCorr(S);
%
%               ....
%
% [C,I] = subRaCrossCorr(S,chan,ind);
%   
% INPUT:
% S:    An MxN coral structure array
% ind:  A vector of shift indices to test; fractional samples accepted.
%       Default: ind = setdiff(linspace(-0.95,0.95,10),0);
% chan: Optional input to cross-correlate the array observations only over
%       the specificed channel in field staChannel.
%
% OUTPUT:
% C:    An MxN cell array containing length(ind) correlograms..
% I:    MxN array providing the shifts for maximum correlation. Unless the
%       user inputs integer shifts, these will usually be fractional shifts
% ind:  The shift in seconds corresponding to I.
%-----------------------------------------------------------------------
% EXAMPLE
% %Suppose you have S and size(S) = [M,N]
% 
% %compute array cross correlation and get max indices ind
% 
% [cc,ind] = raCrossCorr(S);
% 
% %shift the coral structure relative to S(:,end) so that it is aligned with
% %the array observation in S(:,end)
% 
% S        = coralCircShift(S,ind(:,end));
% 
% %create a list of sub-sample shift indices
% 
% subind   = -1:0.1:1;
% 
% %now compute NxN cell array of length(subind) correlograms
% 
% [C,I,ind]   = subRaCrossCorr(S,subind);
% 
% %plot the correlograms between S(:,1) and S(:,4) in blue and between S(:,2) 
% %and S(:,3) in red
% 
% figure; plot(subind,C{1,4}); hold on; plot(subind,C{2,3},'r');
% 
% %plot the aligned seismograms
% 
% Sstack = coralSubShift(S,ind(I(:,end)));
% 
% figure;
% 
% for k = 1:size(S,2)
%  h=coralPlot(Sstack(:,k)); set(h,'color','r'); hold on;
% end;
% 
% h=coralPlot(S(:,1)); set(h,'LineStyle','--');
%----------------------------------------------------------------------- 
% Joshua D Carmichael
% josh.carmichael@gmail.com
% 
% Edit Log
% 12.April.2010
% Made function
%
% 15.July.2010
% Added differential time output, and transposed C. I removed the redundant 
% index output. The first variabled output is now the sub-sample index. If
% S = [S(1),S(2)] then I(1,2) is the subsample shift on S(1) that aligns it
% with S(2). The second variable output is dt. dt(1,2) is t(2) - t(1)
% physically. It is the sample rate times I(1,2).
%
% 05.Aug.2010
% Added an absolute value operation to correlation matrix, C.
% Previous sub-sample shifts did not correctly use the absolute maximum.
%
% stackRaClusters currently outputs a non-optimum stack of traces. Needs
% revision.
%-----------------------------------------------------------------------

ind     = setdiff(linspace(-0.95,0.95,20),0);
N       = [S(:).recNumData];
sInt    = [S(:).recSampInt];

if(std(N) > 1e2*eps)
    
    error('MATLAB:subRaCrossCorr','Record Sections must have data of samel length. Zero pad data');
    
end;

N   = N(1);

if(std(sInt) > 1e2*eps)
    
    error('MATLAB:subRaCrossCorr','Record Sections must have uniform sample rate. Resample.');
    
end;
    
sInt    = sInt(1);

for k = 1:length(varargin)
    
    if(ischar(varargin{k}))
        S = coralExtract(S,'staChannel',varargin{k});
    end;
    
    if(isnumeric(varargin{k}))
        ind = varargin{k};
    end;
    
end;

S       = arrayfun(@coralDetrend,S);
S       = arrayfun(@coralTaper,S);

C       = cell(size(S,2));

% Compute correlograms between input S and it's shifted counterpart
for k = 1:length(ind)
    
    Ssub    = coralSubShift(S,ind(k).*ones(size(S,2),1));
    Sbig    = [S,Ssub];
    G       = raGramMat(Sbig);
    %only keep the upper right quadrant of the big matrix. These values
    %correspond to the inner product between the shifted and non-shifted
    %traces. The rows are the non-shifted traces, the columns are the
    %shifted traces, i.e. G(i,j) = < S(i), S'(j) > where the prime denotes
    %shifting.
    
    G       = G(1:size(S,2),size(S,2)+1:end);    
    %G       = triu(G,1) + tril(G,-1);
    
    %Cell stuffing loop: store correlograms.
    %Each element of C(i,j) contains a correlogram between S(:,i) and the
    %subsample shifted version S(:,j).
    for m = 1:size(S,2)
        
        for n = 1:size(S,2)
            
            %C{m,n} is < S(m), S'(n) > : the second trace/record section is
            %the shifted one. The columns of C{m,n} give the shifted trace
            C{m,n} = cat(1,C{m,n},abs(G(m,n)));
            
        end;
        
    end;
    
end;

%max is computed along columns for each element, meaning the rows of I
%give the max shift indcies for each correlogram C{m,n}. Transpose C so
%that the columns correspond to the non-shifted signals. So now:
%
% C{m,n} is < S(m), S'(n) >
%

C        = cellfun(@full,C,'UniformOutput',false);
[dum,I]    = cellfun(@max,C');

if(nargout>1)
 
    varargout{1} = ind(I);
    
end;

if(nargout>2)
 
    varargout{2} = sInt*ind(I);
    
end;

return;