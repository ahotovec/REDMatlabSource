function [D] = coralCircShift(S,I,varargin)
% This function circularly shifts the data in an MxN coral structure array.
% coralCircShift does not internally use additional coral functions.  Uses
% timeadd; it does require certain fields of the coral structure.
% Required fields: data, staChannel, recNumData, recSampInt, recStartTime.
%
% USAGES
%
% [D] = coralCircShift(S,I);
% [D] = coralCircShift(S,I,'ch');
%
% INPUT:
% S:    An MxN coral structure array containing data to be shifted.
% I:    An Nx1 shift vector.  Accepts sparse values that are output from
%       coralCrossCorr. Change from sparse to full value done internally.
% ch:   Input specifies to shift only the data corresponding to field 
%       staChannel = ch. Optional.
%
% OUTPUT:
% D:    The coral structure with it's data temporally shifted, for all
%       data corresponding to staChannel = ch, if applicable.  The field
%       recStartTime is changed to a value corresponding to the circular
%       shift, corrected for wrap around effects, using timeadd.
%
%-----------------------------------------------------------------------
% EXAMPLE
% %Suppose you have S and size(S) = [M,N]
% 
% %compute array cross correlation and get max indices ind
% 
% [cc,ind] = raCrossCorr(S); 
% C1 = raGramMat(S);
%
% %shift the coral structure relative to S(:,end) so that it is aligned with
% %the array observation in S(:,end)
% 
% %And compute correlation for comparison
% S  = coralCircShift(S,ind(:,end));
% C2 = raGramMat(S);
% 
% %now show that the shifted and stacked record sections are clean....
% 
% Sstack = coralRaStack(S);
% coralPlot(Sstack); hold on;
% 
% %compared to an original record section...
% 
% h = coralPlot(S(:,1)); set(h,'color','r');
%
%-----------------------------------------------------------------------
% Latest Edit: 27.Nov.2009
% Joshua D Carmichael
% josh.carmichael@gmail.com
%-----------------------------------------------------------------------

D       = S;
rows    = [1:size(S,1)]';
ind     = true(size(S,1),1);

%if I comes directly from coralCrossCorr or raCrossCorr, extract nonzero
%index of interest
if(issparse(I))
    I   = full(I);
end;

%remove wrap around effect 
temp     = (I > (S(1).recNumData/2));
I(temp)  = I(temp) - (S(1).recNumData);

%Get channel input if option is input and extract sub-structure channel
for k = 1:length(varargin)
    if(ischar(varargin{k}))
        chan            = varargin{k};
        ind             = regexp({S(:,1).staChannel},chan);
        ind             = logical(cellfun(@length,ind));
        S               = S(ind,:);
    end;
end;

% Extract data and locate time series to be shifted
sint    = mean([S.recSampInt]);
% sint    = [S.recSampInt];
% sint    = reshape(sint,size(S));
% sint    = mean(sint,1);
S=coralPad(S);
M       = [S.data];
M       = reshape(M, S(1).recNumData, size(S,1), size(S,2));
tadd    = sint(:).*I;
rows    = rows(ind);

%Shift each array response function by the amount indicated in I
for n = 1:size(M,3)
    
    M(:,:,n) = circshift(M(:,:,n),I(n));
    
    for k = 1:size(M,2)

        D(rows(k),n).data           = M(:,k,n);
        D(rows(k),n).recStartTime   = timeadd( S(k,n).recStartTime ,tadd(n));        
        
    end;
    
    %Improve this cPick stuffing! Inefficient with loops!
%     if(isfield(D(rows(k),n),'cPick'))
%         
%         cPick           	= D(rows(k),n).cPick;
%         
%         for m = 1:size(cPick,1)
%             
%             [cPick(m).pickTime]	= timeadd([cPick(m).pickTime],tadd(n));
%             
%         end;
%         
%         [D(rows(k),n).cPick] = cPick;
%        
%    end;
    
end;

return;