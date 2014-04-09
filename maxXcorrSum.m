function [ind,Csum] = maxXcorrSum(C)
%find record section that is maximizes the sum of correlations between all
%other record sections.
% 
% USAGE
% [ind,Csum] = maxXcorrSum(C,I);
%
% INPUT
% C:    The cross correlation matrix
% I:    The shift index matrix giving the maximum correlation for the
%       row i, column j of C.
%
% OUTPUT
% ind:  The element in the cross correlation matrix C that maximizes the
%       sum of correlation between all other elements.
% Csum: The maximum mean correlation.
%-----------------------------------------------------------------------
% Joshua D Carmichael
% josh.carmichael@gmail.com
%
% Edit Log
% 21.Feb.2012
% Made function
%-----------------------------------------------------------------------

if(length(C)<3/2)
    
    warning( 'MATLAB:NullCorrelationMatrix', 'Trivial matrix input');
    ind     = NaN;
    Csum    = NaN;
    
else
    
    Csum        = (sum(triu(C,1),2) + sum(triu(C,1),1)')/(length(C)-1);
    [Csum,ind]  = max(full(Csum));
    
end;

return;