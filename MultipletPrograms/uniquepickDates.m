function  C = uniquepickDates(C)
% Removes any date repeat within the cell, so that cell listing gives
% observations which observed on a fixed number of stations.
% That is, the observation dates in Cout{i,j,i} are observed on i stations,
% but none of them are observed on i-1 stations.
%
% USAGES
% C = uniquepickDates(Cin);
%
% INPUT
% Cin:  A multidimensional cell containing 6xM date vectors.
%
% OUTPUT
% C:    The cell with any date repeats removed. Unique dates are left
%       toward the bottom, if redundancy is between a top and bottom cell
%       entry.
%
% EXAMPLE
% This example uses a set of observations by 4 different observers and
% finds every combination of observations within 0.49 second of each other,
% and returns every combination, then finds the observation dates which are
% unique, and yield the largest number of simultaneous observations.
% 
% level     = 3;
% rsec      = 0.49;
% 
% %assign datecell values
%
% datecell          = cell(1,1,4);
% datecell{:,:,1}   = datevec(now)';
% datecell{:,:,2}   = timeadd(datecell{:,:,1},[-2:0.5:2]);
% datecell{:,:,3}   = timeadd(datecell{:,:,1},[-3:0.5:3]);
% datecell{:,:,4}   = timeadd(datecell{:,:,1},[-1:0.5:1]);
% 
% %cluster dates and find unique dates
% C = clusterdates(datecell,rsec,level);
% C = uniquepickDates(C)
%
% %count total number of unique dates 
% L = cellfun('size', C, 2);
% numUnique = sum(L)
%-----------------------------------------------------------------------
% Latest Edit: 19.Marc.2010
% Joshua D Carmichael
% josh.carmichael@gmail.com
%
% Edit Log
% 24.Nov.2009
% Authored function
%
% Edit Log
% 19.March.2010
% Removed datenum operation from set-difference, since datenum is
% imprecise.
%
% 24.June.2010
% Got rid of removal of empty cells because it lost the structure of the
% cell array
%-----------------------------------------------------------------------

%C = Cin;
%ind     = cellfun(@isempty,Cin);
%ind     = not(ind);
%C       = Cin(ind);
N       = numel(C);
L       = 1:N;

%Need to remove usage of datenum...
for k = N:-1:2
    
    Cend        = C(L>=k);
    Cend        = [Cend{:}];
    Ctest       = C(k-1);
    Ctest       = [Ctest{:}];
    %temp        = setdiff( datenum( [Ctest]' ), datenum( [Cend]' ) );
    temp        = setdiff( Ctest', Cend', 'rows')';
    C{k-1}      = temp;
    disp('');
    
end; 
       
return;
%--------------------------------------------------
% Old Code Below
%---------------------------------------------------
% Cout=C;
% for k = 1:numel(C)
%     
%     if(isempty(C{k})), continue; end;
%     
%     for j=k+1:numel(C)
%         
%         if(isempty(C{j})), continue, end;
%         
%         Cout{k}=[setdiff(Cout{k}',C{j}','rows')]';
%         
%     end;
% end;