function C = clusterdates(datecell,varargin)
% This function takes a cell whose (1,1,k) entry is a 6xN matrix of date
% vectors and clusters each entry of the output cell into shared observation
% dates within a given time window, rsec.  Default time window rsec is
% 1sec.  A cell C is returned.  If the input datecell is L entries long, then
% the {1,1,2}(:,:) entry in the return cell C contains the set of observation
% dates shared within rsec between observation set 1 and % observation set
% 2. The {1,1,3}(:,:) entry gives the set of observation dates shared within
% rsec between observation set 1 and observation set 3, and so on for the
% {1,:,:} entries.  The entries {k,m,n} contain the observation dates
% shared within rsec between k observation sets.  An illustrative example
% is given below.
%
% USAGES
% C = clusterdates(datecell)
% C = clusterdates(datecell,rsec);
% C = clusterdates(datecell,rsec,level);
%
% INPUT
% datecell: A (1,1,N) cell of 6xM date vectors, where M is the number
%           of independent observers, and the number of columns need not be
%           equal between cell entries
% rsec:     A time window around which shared observations cluster
% level:    The 'page' of the Cell to be returned if only 'level' number of
%           shared observations is desired as output.
%
% OUTPUT
% C:        An NxNxN-1 cluster cell containing common observations between
%           all combinations of observers
%
%-----------------------------------------------------------------------
% EXAMPLE 1
% This example uses a set of observations by 4 different observers and
% finds every combination of observations within 0.49 second of each other,
% and returns every combination. This example illustrates reducing cluster
% populations with more restrictive timing constraints.
%
% level     = 3;
% rsec      = 1;
% 
% datecell          = cell(1,1,4);
% datecell{:,:,1}   = datevec(now)';
% datecell{:,:,2}   = timeadd(datecell{:,:,1},[-2:0.5:2]);
% datecell{:,:,3}   = timeadd(datecell{:,:,1},[-3:0.5:3]);
% datecell{:,:,4}   = timeadd(datecell{:,:,1},[-1:0.5:1]);
% 
% C = clusterdates(datecell,rsec,level)
% 
% Using rsec = 0.505 as a time difference constraint instead yields more
% clusters per cell:
%
% rsec  = 0.505;
% C     = clusterdates(datecell,rsec,level)
%
%-----------------------------------------------------------------------
% EXAMPLE 2
% This example uses a set of observations by 3 different observers and
% finds every combination of observations within 0.49 second of each other,
% and checks the cluster returns against known cluster returns.
%
% level     = 3;
% rsec      = 1;
% N         = 10;
% 
% staCell   = nchoose({'A','B','C'});
% ind       = cellfun('size',staCell,2);
% staCell   = staCell(ind>1);
% 
% dA1       = datevec(now)';
% dA2       = timeadd(dA1, rsec/2);
% dA3       = timeadd(dA1, N*rsec);
% dA4       = timeadd(dA1, -N*rsec);
% A         = { [dA1,dA2,dA3,dA4] };
% 
% dB1       = dA1;
% dB1(1)    = dA1(1)-1; %NOT in A
% dB2       = dA2; %in A, and within rsec of dA1
% dB3       = timeadd(dB1,N*rsec);
% B         = { [dB1,dB2,dB3] };
% 
% dC1       = dA2; %in A, within rsec of dA2
% dC2       = dA3; %in A
% dC3       = timeadd(dB1,(N+1.05)*rsec);
% dC4       = timeadd(dC3,(N+2.05)*rsec);
% C         = { [dC1,dC2,dC3,dC4] };
% 
% datecell  = [A,B,C];
% D         = clusterdates( datecell, rsec, level);
% 
% Predicted answers:
% 
% D{1} = {combinations of {A,B} }   =>      [dA1, dA2];
% D{2} = {combinations of {A,C} }   =>      [dA1, dA2, dA3];
% D{3} = {combinations of {B,C} }   =>      [dA1 (or) dA2];
% D{4} = {combinations of {A,B,C} } =>      [dA1 (or) dA2];
%
% The (or) ambiguity arises because dA1 and dA2 are within rsec of each
% other, and either one may be retained in the clustering.
%-----------------------------------------------------------------------
% Latest Edit: 24.Nov.2009
% Joshua D Carmichael
% josh.carmichael@gmail.com
%
% Edit Log:
% 24.Nov.2009
% Authored 'OLD CODE' portion of function
%
% 15.March.2010
% Discovered error in combinatorics, found nchoose.m function online
%
% 19.March.2010
% Modified code by including nchoose.m function. Tested with Example 2
% above. Seems to work with that test, in that the computed answers equate
% to the predicted answers.
%-----------------------------------------------------------------------

LS      = length(datecell);
%C       = cell(LS,LS,LS-1); %initialize pTimecell
level   = LS; %page of cluster cell to return
rsec    = 1; %time window (sec) to check for dates around

if nargin>=2
    rsec = varargin{1};
    if(nargin==3)
        level = varargin{2};
    end;
end;

if(level > LS)
    disp('Not that many levels; returning last level')
end;

W   = nchoose(datecell); %returns all the possible combinations of one or more elements of the set datecell
ind = cellfun('size',W,2);

if level>1

[tmp,I] = sort(ind,'ascend');
ind     = ind(I);
W       = W(I);

W   = W(ind>1); %cut out the cells that are just one station (this is why running only one station doesn't work)
ind = ind(ind>1); 

W   = W(ind<=level);

C   = cell(size(W));

for k = 1:size(W,1)
    
    %re-sort D so that the smallest set is intersected first, in order to
    %reduce the size of the input matrices into timediffMat and speed the
    %process
    D       = W{k};
    L       = cellfun('size',D,2); %size of inside of W (how many events at each)
    [L,ind] = sort(L,'ascend');
    D       = D(ind);
    
    for n = 1:size(D,2)        
                
        if(n==1), Dn = D{n};
        
        else           
            
            tmatch = timediffMat(D{n},Dn,rsec);
            Dn      = (unique(tmatch','rows'))';
            
        end;
        
    end;
    
    C{k} = Dn;
    
end;

else
    
    
    
    
end

return;

%-----------------------------------------------------------------------
% OLD CODE
%-----------------------------------------------------------------------
% for k=1:LS;
%     C(k,1,1)=datecell(1,1,k);
%     for n=[1,k+1:LS]
%         if(n>1)
%             
%             temp    = timediffMat(datecell{1,1,k},datecell{1,1,n},rsec);
%             tempN   = unique(datenum(temp'));
%             temp    = datevec(tempN)';
%             
%             C{k,n,1}=(temp);
%         end;
%         if(n==1)
%             C{k,n,1}=[];
%         end;
%     end;
% end;
% 
% for i=2:level; %i=pages in cell-book
%     
%     for k=1:LS;
%         
%         for n=[1,k+1:LS]
%             
%             if(n>1)
%                 temp    = timediffMat(C{k,n-1,i-1},C{k,n,i-1},rsec) ;
%                 tempN   = unique(datenum(temp'));
%                 temp    = datevec(tempN)';
%                 C{k,n-1,i}=( temp );
%                 if(numel(C{k,n-1,i})==0)
%                     C{k,n-1,i}=[];
%                 end;
%             end;
%             
%         end;
%         
%     end;
% end;
% 
% % ind = cellfun(@isempty,C);
% % C   = C(not(ind));
%
%return;