function [varargout] = timediffMat(date1,date2,rsec)
% [varargout] = timediffMat(date1,date2,rsec)
%
% This function computes the difference in times between all pairs of date
% vectors given by 6xN matrices date1 and date2.  date1 and date2 do not
% have to be of the same size.
%
% USAGE
% [tmatch]      = timediffMat(date1,date2,rsec);
% [tmatch,I,J]  = timediffMat(date1,date2,rsec);
% [tmatch,I,J,D]= timediffMat(date1,date2,rsec);
%
% INPUT
% date1:    A 6xN set of date vectors
% date2:    A 6xM set of date vectors. M and N do not have to be equal.
% rsec:     The time window, in seconds, for date cluster cut-off.
%
% OUTPUT
% tmatch:   If date1 has more columns than date2, tmatch gives the columns 
%           of date1 that match within rsec of date2.
% I:        If date1 is longer than date2, 'I' gives the indices for date1 
%           corresponding to tmatch, so that tmatch = date1(:,I).
% J:        The indices for date2 corresponding to tmatch.
% D:        The MxN matrix giving time differences between each column 
%           of date1 and date2.
%-----------------------------------------------------------------------
% Latest Edit: 06.Nov.2009
% Joshua D Carmichael
% josh.carmichael@gmail.com
%-----------------------------------------------------------------------

if( size(date1,2) > size(date2,2) )
    M       = size(date1,2);
    N       = size(date2,2);
    dateLong = date1;
    dateShrt = date2;
else
    M       = size(date2,2);
    N       = size(date1,2);
    dateLong = date2;
    dateShrt = date1;
end;

D      = zeros( M,N );

for k = 1:N %loop over columns
    
    temp    = [abs(timediff( [dateShrt(:,k), dateLong] )) ]';
    D(:,k)  = temp(2:end);
    
end;

[long,short]    = find(D <= rsec);
tmatch          = [dateLong(:,long)];

if(nargout>=1), varargout{1} = tmatch;  end;
if(nargout>=2), varargout{2} = long;    end;
if(nargout>=3), varargout{3} = short;   end;
if(nargout>=4), varargout{4} = D;       end;