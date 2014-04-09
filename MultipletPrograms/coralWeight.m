function [S] = coralWeight(S,varargin)
% Weights data in MxN coral structure array. The user may input either a
% set of reference field names with weights, in the case where data weights
% are complicated and station name dependent, or simply an Mx1 set of
% weights, in the simpler case. See examples below.
%
% USAGES
% [S] = coralWeight(S, varargin);
%
% INPUT
% S:	An Mx1 or MxN coral structure
% W:	(Option 1) An Mx1 vector of weights to multiply the data in S.
% ref:  (Option 2) The reference field(s) to extract data for weighting.
% str:  A list of names for field ref to weight data.
% W:    A vector of weights the same size as str that selectively weights
%       or operates on the data.
% oper: An operation handle for scaling. Default: oper = @times.
%
% OUTPUT
% S:    The input coral structure with appropriately weighted data.
%
%-----------------------------------------------------------------------
% EXAMPLES: Suppose S is a 6x4 coral structure array. Assume the data in
% field 'data' have length 2001.
%
% %multiply all the data by a windowing function g:
% %g = gausswin(2001,5);
% %S = coralWeight(S,g);
%
% %zero-out data for stations CINDY and JAN
% [S] = coralWeight(S, 'staCode', {'CINDY','JAN'}, 0);
%
% %change polarity for stations with channel EPZ
% [S] = coralWeight(S,'staChannel',{'EPZ'}, -1);
%
% %window all data for receiver 'JAN': assume time series has uniform
% %length in structure for this example:
% x             = zeros(2001,1);
% x(900:1100)   = 1;
% [S] = coralWeight(S, 'staCode', {'CINDY','JAN'}, x);
%
% %Add normally distributed noise to station BOBBY
% [Scut] = coralWeight(S, 'staCode',{'BOBBY'},@plus,randn(2001,1));
%
%-----------------------------------------------------------------------
% Latest Edit: 10.Mar.2011
% Joshua D Carmichael
% josh.carmichael@gmail.com
%-----------------------------------------------------------------------

docut   = false;
doflip  = false;
indnum  = cellfun(@isnumeric,varargin);
indhan  = cellfun(@isa,varargin,repmat({'function_handle'},size(varargin)));
indchar = cellfun(@ischar,varargin);
indcell = cellfun(@iscell,varargin);

if(size(S,1)==1)
    
    doflip  = true;
    S       = S(:);
    
end;

%Get input data weights
if(nnz(indnum) > 0)
    
    W      = varargin{indnum};
    
end;

%Get input operation function
if(nnz(indhan) > 0)
    
    oper = varargin{indhan};
    
else
    
    oper = @times;
    
end;

%Get input reference field
if(nnz(indchar) > 0)
    
    docut   = true;
    ref     = varargin{indchar};
    
end;

%Get input cells
if(nnz(indcell) > 0)
    
    incell = varargin{indcell};
    
    %if data are being extracted for specified fields
    if(nnz(cellfun(@ischar,incell)) > 0)
        
        %ind = cellfun(@ischar,incell);
        str = incell;
    
    %the numerical weights for the extracted fields    
    elseif(nnz(cellfun(@isnumeric,incell)) > 0)
        
        ind = cellfun(@isnumeric,incell);
        W   = incell{ind};
        
    end;
    
end;


if(docut)
    
    ind = false(size(S,1),1);
    Li  = 1:size(S,1);
    
    %loop over string inputs to get extraction index
    for k = 1:length(str),
        
        temp    = regexp( {S(:,1).(ref)}, str{k} );
        ind     = ind + logical( cellfun(@length,temp(:)) );
        
    end;
    
    ind     = logical(ind);
    Lw      = Li(ind);
    
else
    
    Lw = 1:size(S,1);
    
end;

%Sp = arithmStruct(Sk,'data',@times,W);
a0      = num2cell(S(Lw,:));
a1      = repmat({'data'},size(S(Lw,:)));
a2      = repmat({oper},size(S(Lw,:)));
a3      = repmat({W},size(S(Lw,:)));

Sw      = cellfun(@arithmStruct,a0,a1,a2,a3);
S(Lw,:) = Sw;

if(doflip)
    
    S   = S(:)';
    
end;

return;