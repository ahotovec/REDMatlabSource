function [varargout] = xcorMutClust(S,varargin)
% [varargout] = xcorMutClust(S,varargin)
% This function uses a partition based, agglomerative clustering method to
% group data from input coral structure S into mutually exclusive sets.
% Data in output coral structure is not normalized. See discussion of
% methodology below. 
% Required fields:recSampInt, data.
%
% USAGES:
% [Scluster,Sorphan]    = xcorMutClust(S,opt);
% [Scluster,Sorphan]    = xcorMutClust(S,opt,maxC);
% [Scluster,Sorphan,I]  = xcorMutClust(S,...);
%
% INPUT:
% S:    an MxN coral structure containing > 1 signals to cluster. S needs
%       two fields, 'data' and 'recSampInt'.
% opt:  an option structure.  Needs field 'thresh' which gives the value
%       for lower bound in the coefficient test.  Must be between 0 and 1.
% I:    The cluster index from which Scluster is obtained from S: e.g., 
%       S(I{k}) = Scluster{k}, k <= length(Scluster)
%
% OUTPUT
% Scluster: a structure containing the waveform clusters.
% Sorphan:  a strucuture containing the waveform orphans (optional).
%
% NOTE:
% (1) make sure the data are detrended prior to usage,
%     i.e., S = coralDetrend(S); OR
%     i.e., S = detrend(S); if S is a matrix
% (2) This function needs maxalign and xcorAlign
% (3) Required Fields: recSampInt, data
%--------------------------------------------------------------------------
% METHODLOGY FOR CLUSTERING
% This function's clustering algorithm operates as follows: it finds all 
% signals that cross-correlate above a level opt.thresh and returns 1 or 2 
% structures giving the stacked cluster waveforms and the 'orphans' 
% respectively.  For any number of signals, clusters are populated by those
% waveforms that don't cross-populate other clusters.
%
% In short, xcorMutClust populates the sets C_k defined by:
%
% C_k = { x : xcor( x1, x2 ) > opt.thresh }
%
% where x1 and x2 are signals stored as column matrices in the structure S.
% The function works by taking the cross-correlation matrix maxC and
% finding the signal pairs such that each signal can be uniquely
% paired with one cluster, which may contain any number of waveforms.
% example, the following index vectors [i,j] give the signal pair indices
% obtained from setdiff form two clusters:
%
% [i j] =
%
% 1     8
% 2     8
% 3     8
% 4     8
% 7     8
% 10    12
% 14    12
%
% There are no common indices in i and j.  Cluster one contains elements:
%
% (1,2,3,4,7,8)
%
% Cluster two contains:
%
% (10,12,14)
%
% Further, the cross correlation between the members of cluster 1 and those
% of cluster 2 correlate BELOW opt.thresh, so that C_1 and C_2 do not
% overlap.
%--------------------------------------------------------------------------
% A cluster is determined as follows:
% suppose (x,y) denotes the correlation coefficient between x and y.
%
% suppose (x,y) > opt.thresh
% suppose (y,z) > opt.thresh
% if (x,z) > opt.thresh, then x, y and z belong to a cluster.
% if (x,z) < opt.thresh, then x, y and z DO NOT belong to a cluster.
%--------------------------------------------------------------------------
%%EXAMPLE: Given coral structure; plot occurrences of a particular cluster
% opt.thresh = 0.80;
% opt.chan   = 'EPZ';
% Sin = coralDetrend(Sin);
% [Sc, So] = coralXcorClust(Sin,opt);
% coralPlot(Sc{1}); %plot the waveforms in cluster 1
% figure; hold on;
% for k=1:size(Sc,2)
%   plot(datenum( [Sc(k).recStartTime]' ), 1-k, 'd', 'markersize',12 );
% end;
% datetick('x','keeplimits');
%
%%EXAMPLE:  A matrix of data containing signals to be clustered
% %Suppose W is a data matrix, so that each column is a time series.
% %Also suppose that the data was recorded using a sampling rate of srate.
% %Then you must first construct a data structure, using the following:
% %Just uncomment and paste the below into matlab and run it.  You need to
% %because this function calls a coral function which needs the field
% %recSampInt.  You must first construct W to run this.  W is your data
% %matrix, described above.
%
% opt.thresh = 0.80;
% datastruc.data = W(:,1);
% datastruc.recSampInt = 0.02; %srate = 0.02 in my example. 
%                              %Replace as needed.
% datastruc     = repmat(datastruc,1,size(W,2));
%
% for k = 1:size(W,2);
%       datastruc(k).data = W(:,k);
% end;
%
% [Wclust,Worphans,clusterIndex] =xcorMutClust(datastruc,opt);
%--------------------------------------------------------------------------
% Joshua D Carmichael
% josh.carmichael@gmail.com
%
% Edit log
%
% 16.Feb.2009
% Authored
% 
% 30.Aug.2009
% Modified to include doublets as part of the cluster output.
%
% 30.Jan.2010
% Modified to ensure signals received equal weighting, by normalizing each
% trace by the euclidean inner product, prior to clustering. Original data
% is not normalized and output seismograms are not normalized.
%
% 06.Aug.2010
% Commented out coralNormData operation. Results in both higher and lower
% correlation values, depending on particular application. Only
% normalization now done in raCrossCorr.m. This methodology is more
% consistent with product space formulation (then kept it).
%
% 25.July.2011
% Removed coralNormData operation; should be in function calling
% xcorMutClust.m. coralWeight was added between this entry and the previous
% 06.Aug.2010 entry for cases in which some kind of cluster weighting is
% performed.
%--------------------------------------------------------------------------

%Sorphan  = S;
Sn       = S;

for k = 1:nargout
    varargout{k} = [];
end;

if(size(S,2) <= 1),
    
    return;
    
end;

default_cluster_thresh = 0.65;
default_maxNum_clust   = 200;
opt.thresh             = default_cluster_thresh;
opt.maxNumClust        = default_maxNum_clust;


for k = 1:length(varargin)
    
    if(isstruct(varargin{k}))
        
        opt = varargin{k};
        
    end;
    
    if(isnumeric(varargin{k}))
        
        maxC = varargin{k};
        
    end;
    
end;

if(~isfield(opt,'thresh'))
    
    opt.thresh = default_cluster_thresh;
    
end;

%--------------------------------------------------------------------------
% Check to see if the cross-correlation matrix maxC was input.  If it was
% not input, then obtain maxC using one of the cross correlation functions
%--------------------------------------------------------------------------
xcheck = exist('maxC','var');

if(~xcheck==1)
    
    if(size(S,1)==1) %time series as data 
        
        if(isfield(opt,'sta2Weight')),
            
            [Sn] = coralWeight(Sn,'staCode',opt.sta2Weight,opt.scalWeight);
            
        end;
        
        [~, ~, ~, ~, maxC] = coralCrossCorr(Sn,opt);
        
    else %array of time series as data
        
        %Sn     = mat2cell(Sn,size(Sn,1),ones(1,size(Sn,2)));
        %Sn     = cellfun(@coralNormData,Sn,'UniformOutput',false);
        %Sn     = [Sn{:}];
        %Sn     = arrayfun(@coralNormData,Sn);
        
        if(isfield(opt,'sta2Weight')),
            
            [Sn] = coralWeight(Sn,'staCode',opt.sta2Weight,opt.scalWeight);

        end;       
            
        if(isfield(opt,'chan'))
            
            [maxC] = raCrossCorr(Sn,opt.chan);
            %[maxC] = raCrossCorr(S,opt.chan);
        
        else
            
            [maxC] = raCrossCorr(Sn);
            %[maxC] = raCrossCorr(S);
            
        end;
        
    end;
    
end;

maxC    = full(triu(maxC,1));
I       = find(maxC >= opt.thresh);

%--------------------------------------------------------------------------
% If there are no cross-correlation matches, return the input structure
%--------------------------------------------------------------------------
if(isempty(I)),
    
    if(nargout>=1), varargout{1} = {};  end; %no above-correlation level clusters
    if(nargout>=2), varargout{2} = S;   end;
    
    return;
    
end;

%--------------------------------------------------------------------------
% Find the signals that cross-correlate with each other above the given
% threshold supplied as opt.thresh.  A cluster is determined as follows:
% suppose (x,y) denotes the correlation coefficient between x and y.
%
% suppose (x,y) > opt.thresh
% suppose (y,z) > opt.thresh
% if (x,z) > opt.thresh, then x, y and z belong to a cluster.
% if (x,z) < opt.thresh, then x, y and z DO NOT belong to a cluster.
%--------------------------------------------------------------------------
[i,j] 	= ind2sub(size(maxC),I);
[dum,Ic]	= getchunks(j);

%If there are only doublets, return the doublets as clusters, and return
%script.
if(isempty(Ic))
    
    clust    = mat2cell( [i,j], ones(size(i,1),1), 2 );
    clust    = cellfun(@transpose,clust,'UniformOutput',false);
    Scluster = clust; %intialize Scluster
    
    for n = 1:length(i)
        Scluster{n} = S(:,[i(n),j(n)]);
    end;
    
    if(nargout>0),
        if(nargout>=1),     varargout{1} = Scluster;                                        end;
        if(nargout>=2),     varargout{2} = S(:,setdiff([1:size(S,2)],vertcat(clust{:})));	end;
        if(nargout>=3),     varargout{3} = clust;                                           end;
    end;
    
    return;
    
end;

jclust          = j(Ic);
clust           = cell(length(jclust),1);

%--------------------------------------------------------------------------
%Go through the set of correlation matches, and identify which groups of
%signals cross correlate with each other above a given value: that is,
%ensure each signal in every cluster correlates above the given level with
%ALL signals in it's cluster.
%--------------------------------------------------------------------------
for n = 1:length(jclust)
    
    iclust      = i(j==jclust(n));
    clust{n}	= [iclust;jclust(n)];
    pairs       = nchoosek( clust{n}, 2 ); %test all pair-wise combinations
    
    for k = 1:length(pairs),
        
        Ck = maxC( pairs(k,1), pairs(k,2) );
        
        if(Ck < opt.thresh )
            
            %if a signal does not correlate highly with others in cluster,
            %remove it
            clust{n} = setdiff(clust{n},pairs(k,1));
            
        end;
        
    end;
    
end;

%--------------------------------------------------------------------------
%Get all possible multiplet pairs, and remove them from the original index
%vector I, to obtain indicies for potential doublets:
%Code block also removes any cluster repeats in doublets. Thus, a subset
%of another larger cluster is not counted. There is no cross-population 
%between clusters.
%--------------------------------------------------------------------------
pairs   = nchoosek(unique(cell2mat(clust)), 2 );
linpair = sub2ind( [size(S,2),size(S,2)], pairs(:,1), pairs(:,2) );
[it,jt] = ind2sub( [size(S,2),size(S,2)], setdiff(I,linpair) );

idoub   = [it,jt];
doublets= mat2cell(idoub', [2],  ones(size(idoub,1),1) )';
L       = cellfun(@length,clust);
[dum,it]  = sort(L,'descend');
clust   = clust(it);
clust   = [clust;doublets];

for k = 0:length(clust)-1,
    
    clust{end-k} = setdiff(clust{end-k}, cell2mat(clust(end-k-1:-1:1)) );
    
end;

L               = cellfun( @length, clust);
clust           = clust(L>1);
[dum,it]          = sort(L(L>1),'descend'); 
clust           = clust(it);
Scluster        = cell(length(clust),1);

for n = 1:length(Scluster)
    
    Scluster{n}	= S(:,clust{n});
    
end;

Sorphan = S(:,setdiff( [1:size(S,2)]',cell2mat(clust) ));

if(nargout>0),
    if(nargout>=1),     varargout{1} = Scluster;    end;
    if(nargout>=2),     varargout{2} = Sorphan;     end;
    if(nargout>=3),     varargout{3} = clust;       end;
end;

% OLD CODE. REMOVED SINCE IT NEGLECTED MATCHES
%--------------------------------------------------------------------------
%Removed singletons, and order the cluster from most to least populous
%--------------------------------------------------------------------------
% L               = cellfun( @length, clust);
% 
% if( length(L(L==1)) > 1),
%     
%     pairs           = nchoosek( cell2mat(clust(L==1)), 2 );
%     singletons      = cell(size(pairs,1),1);
%     
%     %--------------------------------------------------------------------------
%     %Add doublets that may be formed from singleton pairs that correlate above
%     %opt.thresh. Singletons are formed from the setdiff operation but may still
%     %correlate with another singleton from another cluster. Each singleton does
%     %not already belong to another cluster by construction
%     %--------------------------------------------------------------------------
%     if(numel(pairs)>=2)
%         
%         for k = 1:size(pairs,1),
%             
%             if(~(isempty(intersect(cell2mat(singletons),[pairs(k,1);pairs(k,2)]))))
%                 continue;
%             end;
%             
%             Ck = maxC( pairs(k,1), pairs(k,2) );
%             
%             if(Ck >= opt.thresh )
%                 
%                 singletons{k} = [pairs(k,1);pairs(k,2)];
%                 
%             end;
%             
%         end;
%         
%     end;
%     
% else
%     
%     singletons = {};
%     
% end;
% 
% clust           = [clust;singletons];
%-------------------------------------------------------------------------