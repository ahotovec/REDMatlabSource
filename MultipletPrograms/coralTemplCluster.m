function [varargout]=coralTemplCluster(S,W,varargin)
% Clusters data using a template record section. If S is a set of record
% sections and W is a template record section, then S(:,k) is belongs to a
% multiplet containing W if:
%
% (S(:,k),W) >= rho
%
% Where S(:,k) is the k_th record section stored in S. coralTemplCluster
% also works for single receiver correlation.
%
% USAGES
% [Sc]      = coralTemplCluster(S,W)
% [Sc]      = coralTemplCluster(S,W,opt);
% [Sc,So]   = coralTemplCluster(S,W)
% [Sc,So]   = coralTemplCluster(S,W,opt);
% [Sc,So,I] = coralTemplCluster(...);
% [Sc,So,I,xcorval] = coralTemplCluster(...);
%
% INPUT
% S:        An MxN coral structure of N record sections.
% W:        An Mx1 template record section used to cluster around.
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
% Sc:     The cell array of coral structure clusters. Date are low pass
%         filtered data (options for filtering described by field fcut in
%         opt.)
% So:     The coral structure of orphans (non-clusters)
% I:      The index vector from which Sc is obtained from S: that is,
%         S(:,I(k)) = Sc{k}, k <= length(I)
% xcorval: A vector of the xcorr values for all events above the threshold
%-----------------------------------------------------------------------
% Latest Edit: 28.March.2011
% Joshua D Carmichael
% josh.carmichael@gmail.com
%
% Edit Log
%-----------------------------------------------------------------------

%check variable input
ind = cellfun(@isstruct,varargin);

if(nnz(ind) > 0)
    
    opt = varargin{ind};
    
end;    

%intialize inputs and optional opt structure used in different places
if(not(exist('opt','var')))
    
    opt.thresh      = 0.65;
    opt.maxNumClust = 200;
    opt.cutTime     = 5;
    opt.fcut        = 0.35;
    opt.chan        = 'EPZ';
    
end;

%Take array input and identify unique start times
if(size(S,2)>=2)
    
    sepTimes    = diff(timediff([S(1,:).recStartTime,datevec(date)']));
    S           = S(:,abs(sepTimes)>(opt.cutTimeAF));
    
end;

%cut out the station structure in Scut, so that the copy is processed.
Scut        = S;

%Low pass filter the data for correlating and clustering
if(isfield(opt,'fcut'))
    
    sintr       = median([Scut.recSampInt],2);
    Fnyq        = 0.5*1/(sintr);
    Fcut        = repmat({(opt.fcut)*Fnyq},1,size(Scut,2));
    filtType    = repmat({'low'},1,size(Scut,2));
    filtOrder   = num2cell(8*ones(1,size(Scut,2)));
    filtPhase   = repmat({'minimum'},1,size(Scut,2));
    Scut        = cellfun(@coralFilter,num2cell(Scut,1),Fcut,filtType,...
                            filtOrder,filtPhase,'UniformOutput',false);
    Scut        = cell2mat(Scut);      
    W           = coralFilter(W,Fcut{1},'low',8,'minimum');
    
end;

if(isfield(opt,'chan'))
    
    [C]  = raCrossCorr([W,Scut],opt.chan,1);
    
else
    
    [C]  = raCrossCorr([W,Scut],1);
    
end;

tmp = 1:(length(C)-1);
ind = C(2:end)>=opt.thresh;
ind = tmp(ind);
xcorval=C(ind); %added by Kate 1-17-11
Sc  = S(:,ind);
tmp = setdiff(1:(length(C)-1),ind);
So  = S(:,tmp);

if(nargout>=1)
    
    varargout{1} = Sc;
    
end;

if(nargout>=2)
    
    varargout{2} = So;
    
end;


if(nargout>=3)
    
    varargout{3} = ind;
    
end;

%added by Kate 1-17-11
if(nargout>=4)
   
    varargout{4}=xcorval;
    
end