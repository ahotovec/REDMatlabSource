function [varargout] = coralCluster(Scut,varargin)
% Clusters coral structures using xcorMutClust. This partitions data into
% mutually exclusive array response vectors and outputs clusters as a cell
% array and orphans as a structure. The data are filtered according to opt.
%
% USAGES
% [Sc]          = coralCluster(S)
% [Sc]          = coralCluster(S,opt);
% [Sc,So]       = coralCluster(S)
% [Sc,So]       = coralCluster(S,opt);
% [Sc,So,ind]   = coralCluster(S,opt);
%
% INPUT
% S:        An MxN coral structure.
%
% opt:      An optional input structure with the following fields:
% thresh        The correlation threshold above which defines a cluster
% maxNumClust   The maximum number of observations to include in a cluster.
%               Used to avoid running out of memory
% cutTime       The amount of time to cut a seismogram out
% fcut          The percentage of Nyquist to include in pass band of filter
% chan          The channel of the station to cluster over
%
% OUTPUT
% Sc: 	The cell array of coral structure clusters. Date are low pass 
%     	filtered data (options for filtering described by field fcut in 
%     	opt.)
% So: 	The coral structure of orphans (non-clusters)
% ind:	A length N vector of column indices for S such that Sc{k} =
%       S(:,ind{k}).
%-----------------------------------------------------------------------
% Latest Edit: 28.Jan.2010
% Joshua D Carmichael
% josh.carmichael@gmail.com
%
% Edit Log
%
% 28.Jan.2010
% Added sorting at end of function, so that Scluster is returned, with the
% cells sorted by the size of the clusters
%
% 28.Jan.2010
% Scut = coralNormData(S); is now commented out. This operation is done
% internally by raCrossCorr.m. Doing it here is redundant.
%
% 28.March.2011
% Filtering now uses minimum phase to avoid attenuating arrival times of
% the waveforms, and is only done if fcut is a field in opt.
%
% 19.May.2011
% A loop that filters the data now is done via cellfun rather than using
% the (previous) loop. Added index output. Over-sized structure condition
% removed (opt.maxNumClust) since it should be added to raPicks2raCluster.m
% function instead. Index output added.
%
% 04.Nov.2011
% Removal of duplicate columns has been placed in raPicks2raCluster now.
%-----------------------------------------------------------------------

%check input
ind = cellfun(@isstruct,varargin);

if(nnz(ind)>0),
    
    opt = varargin{ind};
    
else
    
    %intialize inputs and optional opt structure used in different places    
    opt.thresh      = 0.7;
    opt.maxNumClust = 200;
    opt.cutTime     = 5;
    opt.fcut        = 0.35;
    opt.chan        = '';
    
end;

%cut out the station structure in Scut, so that the copy is processed.
%Scut    = S;

%Low pass filter the data for correlating and clustering
if(isfield(opt,'fcut')),
    
    sintr       = median([Scut.recSampInt],2);
    Fnyq        = 0.5*1/(sintr);
    Fnyq        = repmat({(opt.fcut)*Fnyq},1,size(Scut,2));
    filtType    = repmat({'low'},1,size(Scut,2));
    filtOrder   = num2cell(8*ones(1,size(Scut,2)));
    filtPhase   = repmat({'minimum'},1,size(Scut,2));
    Scut        = cellfun(@coralFilter,num2cell(Scut,1),Fnyq,filtType,...
                            filtOrder,filtPhase,'UniformOutput',false);
    Scut        = cell2mat(Scut);
    
end;

%Perform mutually exclusive partition based clustering
%if( size(Scut,2) < opt.maxNumClust )
    
[Scluster,Sorphan,ind]  = xcorMutClust(Scut,opt);
%ind                     = indOrg(ind);

if(not(isempty(Scluster)))
    
    L           = cellfun('size', Scluster, 2);
    [~,ij]      = sort(L,'descend');
    Scluster    = Scluster(ij);
    
end;

if(nargout>=1)
    
    varargout{1} = Scluster;
    
end;

if(nargout>=2)
    
    varargout{2} = Sorphan;
    
end;

if(nargout>=3)
    
    varargout{3} = ind;
    
end;

return;

%-----------------------------------------------------------------------
%OLD CODE: USED FOR OVER-SIZED CORAL STRUCTURES DUE TO MEMORY PROBLEMS
%else
    
%     count = 0;
%     %get number of columns of Scut in initial input before clustering
%     fullN = min(opt.maxNumClust,size(Scut,2)));
%     
%     while(size(Scut,2)>=1)
%         
%         [Sc,So,I]  = xcorMutClust(Scut(:,1:min(opt.maxNumClust,size(Scut,2))),opt);
%         
%         if(count==0),
%             
%             Scluster    = Sc; 
%             Sorphan     = So; 
%             ind         = I;
%             
%         end;
%         
%         if(count > 0)
%             
%             Scluster	= cat(1,Scluster,Sc);
%             Sorphan  	= cat(2,Sorphan,So);
%             %temp        = max(cellfun(@max,ind));
%             I           = cellfun(@plus,I,num2cell(fullN*ones(size(I))),'UniformOutput',false);
%             %I           = cellfun(@plus,
%             ind         = cat(1,ind,I);
%             
%         end;
%         
%         count = count + 1;
%         Scut(:,1:min(opt.maxNumClust,size(Scut,2))) = [];
%         
%         if(isempty(Scut)), 
%             
%             break, 
%             
%         end;
%         
%     end;
%     
% end;