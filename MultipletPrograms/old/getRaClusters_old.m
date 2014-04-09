function  [multtempl] = getRaClusters(opt, varargin)
% Returns a cell array of clusters obtained from combining all CLUSTERPICK
% files in the current directory and clustering them together.  A .mat file
% is saved in the current directory containing the cluster cell array, with
% the multiplet template and multiplet structures
%
% USAGES
% [Scut] = getRaClusters(opt, varargin)
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
% multtempl:    A structure containing the following fields and contents:
% 
% multiplets    A Px1 cell array of N x M_k coral structures that were
%               obtained by clustering the database with Nx1 coral
%               structure templates. M_k is the number of multiplets found
%               using coral structure template "k", and P is the total
%               number of multiplets found. The data is unnormalized.
% template      A cell array of Nx1 coral structures that were used as
%               templates to identify the multiplets above.
% clustIndex    The indices of the coral structure columns that correspond
%               to multiplets. The coral structure column, in this case, is
%               the very large coral structure obtained by catenating all
%               CLUSTERPICK*.mat files.
% orphIndex     The orphan events for the coral structure obtained by
%               catenating all CLUSTERPICK*.mat files.
%-----------------------------------------------------------------------
% Joshua D Carmichael
% josh.carmichael@gmail.com
%
% Edit Log
% 31.Mar.2010
% Had to alter code so that files without Scluster fields are passed over
% and not including in the picking process; otherwise dayStart, yearInput,
% and dayEnd fields might not be correctly saved.
%
% 25.Oct.2011
% multtempl structure saved with template and clusters saved independently
% 
% 20 Jan 2011
% Kate edited to make it save the xcorrelation values and to save the
% unnormalized data instead of normalized.
%-----------------------------------------------------------------------

if(isempty(varargin))
    
    flist   = getFiles(['CLUSTERPICKS.ARRAY.0',num2str(opt.thresh*100)],'mat');
    
else
    
    if(length(varargin)>=1)
        
        flist = getFiles(['CLUSTERPICKS.ARRAY.0',num2str(opt.thresh*100)],'mat',varargin{:});
        
    end;    
    
end;

Scat        = [];
yearStart   = [];
dayStart    = [];

for k = 1:length(flist)
    
    tic;
    
    output  = load(flist{k});
    
    if(not(isfield(output,'Scluster'))),
        
        continue;
        
    end;
    
    Sc      = output.Scluster;
    
    temp        = textscan(flist{k},'%s','Delimiter','.');
    yearStart   = cat(1,yearStart,str2double(temp{1}(end-2)));
    dayStart    = cat(1,dayStart,str2double(temp{1}(end-1)));        
    Scat        = cat(1,Scat,Sc);
    
    toc
    
end;

day1        = min(dayStart(yearStart<=min(yearStart)));
day2        = max(dayStart(yearStart>=max(yearStart)));
year1       = min(yearStart);
year2       = max(yearStart);

%make room in memory for the intensive clustering process
clear flist;
clear output;
clear Sc;
clear Scluster;

%Sort the catenated cluster arrays in descending order in multiplet size
Lclus       = cellfun('size',Scat,2);
[~,ind]     = sort(Lclus,'descend');
Scat        = Scat(ind);

%KEEP FOR SYNTAX REFERENCE
%Scat  	= cellfun(@arrayfun,repmat({@coralNormData},size(Scat)),Scat,'UniformOutput',false);

%Make a cell array of stacked waveforms, aligned to sub-sample precision.
%Detrending, tapering, and normalizing are internally performed in
%stackRaClusters. This pre-processing determines optimal shifts. However,
%the stacked data is not not processed and contains un-normalized data. St
%contains template record sections.

disp(sprintf('Stacking %i multiplet record sections',length(Scat))); tic;
[St] = stackRaClusters(Scat,opt);
disp(sprintf('Done Stacking the %i multiplet record sections',length(Scat))); toc

%Detrend, taper, and normalize a copy of the data for template clustering
Sk   = cellfun(@coralDetrend,Scat,'UniformOutput',false);
Sk   = cellfun(@coralTaper,Sk,'UniformOutput',false);
Skorig= Sk; %save a copy of unnormalized data
Sk   = cellfun(@arrayfun,repmat({@coralNormData},size(Sk)),Sk,'UniformOutput',false);

St   = cellfun(@coralDetrend,St,'UniformOutput',false);
St   = cellfun(@coralTaper,St,'UniformOutput',false);
St   = cellfun(@arrayfun,repmat({@coralNormData},size(St)),St,'UniformOutput',false);
Storig=St;
%---------------------------------

%Replace current row of template cell array St with matches.
%clust = cell(size(Sk,2),1);
if isempty(Scat)
    fprintf('No clusters in month, ending getRaClusters\n');
    multtempl=[];
return
end
disp(sprintf('Clustering %i multiplet record sections',length(Scat))); tic;


%Scat is the original record section data. Is there a reason to save
%original Scat, pre-processing? I don't think so....
%Scat        = [Scat{:}];

%Sk is the processed record section data, Skorig is the unnormalized
%processed record section data
Sk          = [Sk{:}];
Skorig       = [Skorig{:}];

%Scat is a copy of the processed record section data, Scatorig is copy of
%unnormlaized data
Scat        = Sk;
Scatorig     = Skorig;

%Sm is a copy of the template record sections
Sm          = St;


%Get the index vectors for the cluster templates, and record sections
clustInd    = cell(size(St,1),1);
allInd      = 1:size(Sk,2);
remInd      = allInd;

%Consider the possibility that Sk and St have different fields for each
%element, as Sk or St may have been processed by other functions that may
%have changed their field elements or numbers.
f1 = fieldnames(Sk(:,1));

for k = 1:length(St),
    
    f2        = fieldnames(St{k});
    df        = setxor(f1,f2);
    
    if(length(f1)>length(f2)),
        
        [~,~,ind,xcorval] = coralTemplCluster(rmfield(Sk,df),St{k},opt);
        
    elseif(length(f2)>length(f1)),
        
        [~,~,ind,xcorval] = coralTemplCluster(Sk,rmfield(St{k},df),opt);
        
    else
        
        [~,~,ind,xcorval] = coralTemplCluster(Sk,St{k},opt);
        
    end;
    
    %pickInd are the indices from Scat that cluster
    pickInd     = remInd(ind);
    
    %remInd are the indices of Scat that do not cluster
    remInd      = setdiff(remInd,pickInd);
    
    %ind are the indices of Sk that DO cluster. I re-stuff St with the
    %cluster elements of Sk. Do same with original data
    St{k}       = Sk(:,ind);
    Storig{k}   = Skorig(:,ind);
    
    %place the xcorrelation value between the template and each event with its corresponding event
    %added by Kate 1-17-2012
    %St{k}=arrayfun(@setfield,St{k},'xcorval',xcorval'); 
    for q=1:size(St{k},2)
       St{k}(1,q).xcorval=xcorval(q);
       Storig{k}(1,q).xcorval=xcorval(q);
    end
    
    %reinitialize Sk to be the remainder of Scat that does not contain
    %previously found clusters. same with Sorig
    Sk          = Scat(:,remInd);    
    Skorig       = Scatorig(:,remInd);
    
    %get the cluster pick indices
    clustInd{k} = pickInd;
    
    if(min(size(Sk))<1), break; end;
    
end

clear Sk;
clear Skorig;
ind = setdiff(allInd,[clustInd{:}]);
St  = cat(1,St,{Scat(:,remInd)});
Storig = cat(1,Storig,{Scatorig(:,remInd)});

%----------------

%disp(sprintf('Stacking the %i multiplet record sections',length(St))); tic;
%[St] = stackRaClusters(St,opt);
%disp(sprintf('Done Stacking the %i multiplet record sections',length(St))); toc

saveFile = sprintf('%s.%s.%.3i.%i.%.3i.%i.%.3i.mat','MULTTEMPL',...
    (opt.chan),round(100*(opt.thresh)),year1,day1,year2,day2);
disp(saveFile)
indFull                = not(cellfun(@isempty,Sm));
%indFull                = not(indFull);
St                     = St(indFull);
Storig                 = Storig(indFull);
clustInd               = clustInd(indFull);
Sm                     = Sm(indFull);
indSize                = cellfun('size',St,2);
%get rid of sets with 2 or fewer sets
St                     = St(indSize>2);
Sm                     = Sm(indSize>2);
clustInd               = clustInd(indSize>2);
Storig                 = Storig(indSize>2);
multtempl.multiplets   = Storig;
%multtempl.multnorm     = St;
multtempl.template     = Sm;
multtempl.clustIndex   = clustInd;
multtempl.orphIndex    = ind;
save(saveFile,'multtempl');

disp(sprintf('%i cluster sets found',size(St,1)))
%disp(sprintf('%i orphans left',size(Sorphan,2)))

return;