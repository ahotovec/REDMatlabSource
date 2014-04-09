function  [clusterstack] = getRaClusterStack(opt, varargin)
% Stacks array CLUSTERHISTORY .mat files; the raw data corresponding to the
% filtered data stored in the cells of CLUSTERHISTORY is what is stacked.
%
% USAGES
% [clusterstack] = getRaClusterStack(opt, varargin)
%
% INPUT
% opt:      A required input structure with the following fields:
% chan      The channel as it appears in the pick date save file/the coral
%           structure for the array.
%
%                   ...
%
% sampRate: Additional file specifiers for reading in files. This example
%           is sampRate.
%
% OUTPUT
% clusterstack: A structure array with field coralstruc containing the
%               coral structure of stacks.
%
% A file named CLUSTERSTACK.*mat is saved containing clusterstack
%
% See Also: stackRaClusters.m, coralRaStack.m
%-----------------------------------------------------------------------
% Joshua D Carmichael
% josh.carmichael@gmail.com
%
% Edit Log
% 28.Feb.2010
% Made function
%
% 17.April.2010
% Function now uses sub-sample optimization of correlation function prior
% to stacking. Uses subRaCrossCorr.m
%
% 08.June.2010
% Function now saves 'clusterstack' as a coral structure with extra field
% 'multiplets' rather than a cell: every stack has only column
%-----------------------------------------------------------------------

saveFlag        = true;
clusterstack    = [];

if(isempty(varargin))
    
    fthresh     = sprintf('%i%i',0,100*(opt.thresh));
    flist       = getFiles('CLUSTERHISTORY','mat',(opt.chan),fthresh);
    
else
    
    if(length(varargin)>=1)
        
        flist = getFiles(varargin{:},'mat',(opt.chan));
        
    end;
    
    
end;

for k = 1:length(flist)
    
    tic
    
    fcell       = flist{k};
    output      = load(flist{k});
    
    if(not(isfield(output,'Scluster')))
        
        continue;
        
    else
        
        Sc          = output.Scluster;
        
    end   
    
    %output D is an MxN coral structure array
    [D] = stackRaClusters(Sc,opt);
    
    if( saveFlag )
        
        %Sc is a cell, So is a structure array
        [Sc,So] = coralCluster(D,opt);
        
        %Sc is now a structure array
        Sc      = stackRaClusters(Sc,opt);
        
        D       = [Sc,So];
        L       = zeros(size(D,2),1);
        
        for kk = 1:size(D,2)
            
            L(kk) = size( [D(1,kk).multiplets], 2);
            
        end;
        
        [L,ind]= sort(L,'descend');
        
        D      = D(:,ind);
        
        %Output coral structure
        clusterstack.coralstruc = D;
        clear D;
        
        C      = textscan(fcell,'%s','Delimiter','.');
        C      = C{1};
        this   = strcat(C(2:end-1),'.');
        this   = strcat('CLUSTERSTACK','.',[this{:}],'mat');
        
        disp(sprintf('Saving Stacked Waveforms for %s',fcell));
        
        saveFile = sprintf('%s',this);
        save(saveFile,'clusterstack');
        
    end;
    
end;

return;