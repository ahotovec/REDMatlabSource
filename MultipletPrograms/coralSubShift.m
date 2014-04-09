function [D] = coralSubShift(S,I,varargin)
% This function circularly shifts the data in an MxN coral structure by
% a subsample time interval, using the Fourier shift theorem in the 
% frequency domain. It uses a worker function fshift.m
%
% coralSubShift does not internally use additional coral functions; uses
% timeadd; it does require certain fields of the coral structure.
% Required fields: data, staChannel, recNumData, recSampInt, recStartTime.
%
% USAGES
% [D] = coralSubShift(S,I);
% [D] = coralSubShift(S,I,'ch');
%
% INPUT
% S:    An MxN coral structure array containing data to be shifted.
% I:    An Nx1 shift vector, in fractions of a sample; Accepts sparse 
%       values. Change from sparse to full value done internally.
% ch:   Input specifies to shift only the data corresponding to field 
%       staChannel = ch. Optional.
%
% OUTPUT
% D:    The coral structure with it's data temporally shifted, for all
%       data corresponding to staChannel = ch, if applicable.  The field
%       recStartTime is changed to a value corresponding to the circular
%       shift, corrected for wrap around effects, using timeadd.
%
% NOTE: This function needs fshift.m. This is not available from MATLAB,
%       and comes from the MATLAB central file exchange. *If* it is not
%       included in the package your received, you can download it.
%       However, I generalized it to accept arrays. If you do not have the
%       corrected version, add this line in the function after p is
%       defined: p = repmat(p,1,size(x,2));
%-----------------------------------------------------------------------
% Latest Edit: 22.March.2010
% Joshua D Carmichael
% josh.carmichael@gmail.com
%
% Edit Log
% 22.March.2010
% Designed function
%
% 12.April.2010
% Corrected fshift.m to work with matrices as well, so shifting is done
% column wise. Contact me if you do not have this version of fshift. To
% correct yourself, add the following line in fshift.m after defining p:
%
% p = repmat(p,1,size(x,2)); 
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

%Get channel input if option there, and extract sub-structure channel
for k = 1:length(varargin)
    if(ischar(varargin{k}))
        S = coralExtract(S,'staChannel',varargin{k});
    end;
end;

% Extract data and locate time series to be shifted
sint    = mean([S.recSampInt]);
M       = [S.data];
M       = reshape(M, S(1).recNumData, size(S,1), size(S,2));
tadd    = sint(:).*I;
rows    = rows(ind);

%Shift each array response function by the amount indicated in I
for n = 1:size(M,3)
    
    M(:,:,n) = fshift(M(:,:,n),I(n));
    
    for k = 1:size(M,2)

        D(rows(k),n).data           = M(:,k,n);
        D(rows(k),n).recStartTime   = timeadd(S(k,n).recStartTime ,tadd(n));
        
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
%     end;
    
end;

return;