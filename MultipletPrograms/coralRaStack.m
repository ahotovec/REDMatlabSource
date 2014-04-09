function [D] = coralRaStack(S,varargin)
% This function stacks the data in an MxN coral struture. This function
% uses the median, not the mean to average the data. In order to create
% coherently aligned and stacked data, use coralCircShift first, or
% coralSubShift and subRaCrossCorr together for sub-sample shifted data
% first.
%
% USAGES
%
% [D] = coralRaStack(S);
%
% INPUT:
% S:    An MxN coral structure array containing data to be stacked to form
%       an Mx1 coral stucture.
%
% OUTPUT:
% D:    The data stack, output as an Mx1 coral structure. The data field 
%       contains the data stack, and the recStartTime fields each contain
%       a 6xN date vector for each observation in the stack. NOTE: An
%       additional field, 'multiplets' is added to D. D.multiplets contains
%       the date vectors corresponding to the stacked data, without date
%       repeats.
% 
%
% See Also: stackRaClusters.m, getRaClusterStack.m
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
% %now show that the shifted and stacked record sections are clean....
% 
% Sstack = coralRaStack(S);
% coralPlot(Sstack); hold on;
% 
% %compared to an original record section...
% 
% h = coralPlot(S(:,1)); set(h,'color','r');
%-----------------------------------------------------------------------
% Latest Edit: 13.Jan.2010
% Joshua D Carmichael
% josh.carmichael@gmail.com
%
% Edit Log
% 13.Jan.2010
% Made function
%
% 14.April.2010
% (1) Now no shift indices are computed internally if one is not provided.
%     The function does no shifting if no indices are provided.
% (2) coralSubShift has replaced coralCircShift in case fractional shifts
%     are provided via optional input ind.
%
% 17.May.2010
% Now adds a field called multiplets. Multiplets is an array of datevectors
% giving the observation time of each multiplet.
% Field recStartTime holds the first pick date corresponding to the first 
% multiplet.
%
% 08.June.2010
% If field multiplet already exists, the multiplet dates are catenated.
%
% 08.June.2010
% If field multiplet already exists, the multiplet dates are catenated.
%-----------------------------------------------------------------------

%get any column of input S in order to retain the organization of the coral
%structure
D       = S(:,1);

% Extract variable inputs
for k = 1:length(varargin),
    
    if(isnumeric(varargin{k}))
        
%add stuff in here if desired....like a weight matrix?
        
    end;
    
end;

% Extract data, stuff into a matrix, and then stack data (take median)

M       = [S.data];
M       = reshape(M,D(1).recNumData,size(S,1),size(S,2));
M       = median(M,3);

%initialize D's multiplet fields if they do not exist
if(not(isfield(D,'multiplets')))
    
    [D(1:size(S,1)).multiplets] = deal(0);
    
end;

%structure stuffing loop: loop over receivers
for k = 1:size(S,1)
    
    D(k).data           = M(:,k);
    
    %if there is no multiplet field, add one to the structure
    if(not( isfield(S(:,1), 'multiplets') ) )
        
        D(k).recStartTime   = [S(k,1).recStartTime];
        D(k).multiplets     = unique( [S(k,:).recStartTime]', 'rows')';
        comment             = sprintf('%s',' stacked record sections');
        D(k).recLog         = strcat(D(k).recLog,comment);
        
    else
        
        %if a multiplet field exists, catenate the other multiplets on,
        %taking care not to repeat S(k,1).multiplet which is already in D
        temp                = unique( [S(k,2:end).multiplets]', 'rows')';
        D(k).multiplets     = cat(2,[D(k).multiplets],temp);
        
    end;
    
%      if(isfield(S(:,1), 'cPick'))
%          
%          temp               = [S(k,2:end).cPick];
%          D(k).cPick         = cat(2,[D(k).cPick],temp);
%          
%      end;
    
end;

return;