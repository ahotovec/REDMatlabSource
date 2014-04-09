function [Sc]=quickStackTop(Sc,opt) 

%Same at quickStack but uses only the first station to get the lag times
%instead of cross correlating the entire station matrix - Kate A wrote this

%(copied and slightly modified from stackSubRaCLusters)
% Stacks coral structure record sections contained in Mx1 cell arrays. Each
% row of cell array contains a set of record sections, possibly of
% different sizes. That is, row j of Sc may contain a K x N_j coral
% structure array, where K is the number of receivers in the array, and N_j
% is the number of observations/record sections. Output data is detrended,
% tapered, and traces are individually normalized.
% NOTE:
% Alternative function for stackRaClusters.m. This function permits usage
% of a sub-array for time-lag alignment, but all record sections are
% stacked. See fields staWeight and numWeight.
%
% USAGES
% [D] = stackSubRaClusters(Sc,opt);
%
% INPUT
% Sc:       An Mx1 cell array, with each row containing a coral structure
%           array or record section, arranged by columns.
% opt:      A required input structure with the following fields:
%
% chan      The channel as it appears in the pick date save file/the coral
%           structure for the array.
% staWeight (optional) A set of station codes to weight in order to reject
%           from the correlating
% numWeight (optional) A set of weights (such as 0) to suppress or amplify
%           the effects of stacking. For example, the true amplitudes of
%           the ground motion could be used as well.
% lpfilt    lowpass filter of data
%
% OUTPUT
% Sc:       A strucutre array containing the cluster stacks, with an added
%           field multiplets, containing the date vectors indicating the
%           stacked data.

%compute the array correlation for each cell element
ind     = cellfun(@isstruct,Sc);
Sc      = Sc(ind);
Sc      = cellfun(@coralResample,Sc,'UniformOutput',false);
Sc      = cellfun(@coralPad,Sc,'UniformOutput',false);
Sc      = cellfun(@coralDemean,Sc,'UniformOutput',false);

for k = 1:size(Sc,1)
%Filter first if specified in opt file
    if length(opt.bpfilt)==2
        Sc{k}=coralDemean(Sc{k});
        Sc{k}=coralTaper(Sc{k});
        Sc{k}=coralFilter(Sc{k},opt.bpfilt,'bandpass',4,'minimum');
    elseif opt.lpfilt>0 && opt.bpfilt==0
        Sc{k}=coralDemean(Sc{k});
        Sc{k}=coralTaper(Sc{k});
        Sc{k}=coralFilter(Sc{k},opt.lpfilt,'low',4,'minimum');
    end
end

for i=1:size(Sc,1) %Sk is just the data for the first station
Sk{i,1}      = Sc{i}(1,:);

end
%optc 	= repmat(num2cell(opt),1,size(Sc,1));

%Detrend, taper and normalize the data. Tapering and detrending are also
%done internally by raCrossCorr. It is done here in order to ensure the
%traces in each record section receive equal weighting when
%cross-correlating in raCrossCorr.m
Sk      = cellfun(@coralDetrend, Sk,'UniformOutput',false);
Sk      = cellfun(@coralTaper, Sk,'UniformOutput',false);
Sk      = cellfun(@arrayfun,repmat({@coralNormData},size(Sk)),Sk,'UniformOutput',false);

%Check optional weighting of stations with fields staWeight and numWeight
if(and(isfield(opt,'staWeight'),isfield(opt,'numWeight'))),
        
    [~,Ic]  = cellfun(@raCrossCorr,cellfun(@coralWeight,Sk,repmat({'staCode'},size(Sk)),...
                       repmat({opt.staWeight},size(Sk)),repmat({opt.numWeight},size(Sk)),...
                       'UniformOutput',false),repmat({opt.chan},size(Sk)),'UniformOutput',false);
                   
else
    
    [~,Ic]  = cellfun(@raCrossCorr,Sk,repmat({opt.chan},size(Sk)),'UniformOutput',false);
    
end;

for k = 1:size(Sc,1)
    
    %shift each record section to get maximum correlation with the last
    %record section in Sc using correlation values from Sk
       
    Sc{k}   = coralCircShift(Sc{k},Ic{k}(:,end));       
    
    %Sk{k}   = coralRaStack(Sk{k});
    Sc{k}   = coralRaStackNoZero(Sc{k});   %stacking that doesn't include flatlined data in median calculation  
    
end;