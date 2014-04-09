function  [Sc] = stackRaClusters(Sc,opt)
% Stacks coral structure record sections contained in Mx1 cell arrays. Each
% row of cell array contains a set of record sections, possibly of
% different sizes. That is, row j of Sc may contain a K x N_j coral
% structure array, where K is the number of receivers in the array, and N_j
% is the number of observations/record sections. Output data is detrended,
% tapered, and traces are individually normalized.
%
% USAGES
% [D] = stackRaClusters(Sc,opt);
%
% INPUT
% Sc:       An Mx1 cell array, with each row containing a coral structure
%           array or record section, arranged by columns.
% opt:      A required input structure with the following fields:
% chan      The channel as it appears in the pick date save file/the coral
%           structure for the array.
%
% OUTPUT
% Sc:       A strucutre array containing the cluster stacks, with an added
%           field multiplets, containing the date vectors indicating the
%           stacked data.
%
% See Also: getRaClusterStack.m
%-----------------------------------------------------------------------
% Latest Edit: 15.May.2010
% Joshua D Carmichael
% josh.carmichael@gmail.com
%
% Edit Log
% 08.June.2010
% Changing output of function to be a structure, not a cell.
%
% 11.July.2010
% Made alteration to subRaCrossCorr. Changed index picking to ind(I(1,:)).
% This is the correct optimization.
%
% 19.July.2011
% Body of code now heavily uses cellfun now. Used script_testSubShifting
% code in place of old stackRaClusters.
%
% 25.July.2011
% Data returned is now not normalized, but left in original form.
% Processed copy of record sections are used to determine shift indices and
% 
%-----------------------------------------------------------------------

%compute the array correlation for each cell element
ind     = cellfun(@isstruct,Sc);
Sc      = Sc(ind);
Sc      = cellfun(@coralResample,Sc,'UniformOutput',false);
Sc      = cellfun(@coralPad,Sc,'UniformOutput',false);
Sk      = Sc;
%optc 	= repmat(num2cell(opt),1,size(Sc,1));

%Detrend, taper and normalize the data. Tapering and detrending are also
%done internally by raCrossCorr. It is done here in order to ensure the
%traces in each record section receive equal weighting when
%cross-correlating in raCrossCorr.m
Sk      = cellfun(@coralDetrend, Sk,'UniformOutput',false);
Sk      = cellfun(@coralTaper, Sk,'UniformOutput',false);
Sk      = cellfun(@arrayfun,repmat({@coralNormData},size(Sk)),Sk,'UniformOutput',false);
[~,Ic]  = cellfun(@raCrossCorr,Sk,repmat({opt.chan},size(Sk)),'UniformOutput',false);

%compute the Gram Matrix between array elements
% [M]     = cellfun(@raGramMat,Sc,repmat({'EPZ'},size(Sc)),'UniformOutput',false); 
% Mk      = M;
% Mj      = M;

%loop over each cluster and circularly shift it to match the last record
%section. 
%Ensure subind has the zero shift as an index.
subind  = linspace(-1.05,+1.05,24);
subind  = unique([0;subind(:)]);

for k = 1:size(Sk,1)
    
    %shift each record section to get maximum correlation with the last
    %record section in Sc
    Sk{k}   = coralCircShift(Sk{k},Ic{k}(:,end));   
    
    %compare sample-width aligned correlations
    %Mk{k}   = raGramMat(Sk{k},'EPZ');
    %Mk{k}   = Mk{k} - M{k};
        
    %now sub-sample shift:
    [dum,tmp2]    = subRaCrossCorr(Sk{k},opt.chan,subind);
    Sk{k}       = coralSubShift(Sk{k},tmp2(:,end));
    
    M           = raGramMat(Sk{k},opt.chan);
    
    %test Negative correlation values. The negative polarity signals are
    %given by logical true ind values.
    isNeg       = sign(M) + eye(size(M));    
    isNeg       = isNeg(:,end);
    ind         = isNeg<0;
    
    %Now shift original data, and check for needed polarity correction
    Sc{k}       = coralSubShift(Sc{k},tmp2(:,end)+Ic{k}(:,end));
    
    if(nnz(ind)>0)
        
        %Sk{k}(:,ind)= coralWeight(Sk{k}(:,ind),-1);
        Sc{k}(:,ind)= coralWeight(Sk{k}(:,ind),-1);
        
    end;
    
    %Sk{k}   = coralRaStack(Sk{k});
    Sc{k}   = coralRaStack(Sc{k});     
    
end;

return;