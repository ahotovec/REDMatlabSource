function h_plot = plotRaClusterStacks(Sc,opt,varargin)
% Returns a set of plots of record sections for coral structures in cell
% array Sc. Required fields in Sc include cPicks.
%
% USAGES
% [h] = plotRaClusterStacks(Sc,opt);
% [h] = plotRaClusterStacks(Sc,opt,phLabel);
% [h] = plotRaClusterStacks(Sc,opt,...,Nclus);
%
% INPUT
% Sc:       An Mx1 cell array containing the coral structures for each
%           cluster. Coral structures can be arrays.
% opt:      A structure with field chan.
% phLabel:  (optional) String. Input 'phase' specifies labeling each text
%           arrow on channel trace opt.chan.
% t_unc:    (optional). 2-vector specifying time window buffer around
%           pick times which to plot record sections. Default: [0.1,0.1];
%
% OUTPUT
% h:        An output of plot handles
%-----------------------------------------------------------------------
% Joshua D Carmichael
% josh.carmichael@gmail.com
%
% 28.Feb.2010
% Made function
%
% 03.Aug.2011
% Changed function to plot record sections.
%-----------------------------------------------------------------------
t_unc   = 0.1*ones(1,2);
phLabel = 'none';
Nclus   = length(Sc);

charInd = cellfun(@ischar,varargin);
vecInd  = cellfun(@isnumeric,varargin);
scalInd = cellfun(@isscalar,varargin);
vecInd  = logical(vecInd.*not(scalInd));

if(nnz(charInd)>0)
    
    phLabel = varargin{charInd};
    
end;

if(nnz(vecInd)>0)
    
    t_unc  = varargin{vecInd};
    
end;

if(nnz(scalInd)>0)
    
    Nclus   = varargin{scalInd};
    
end;

LogcPick = cellfun(@isfield,Sc,repmat({'cPick'},size(Sc)));

if(Nclus > nnz(LogcPick))
    
    disp(sprinf('Only %i clusters input, only %i record sections plotted',Nclus));
    
end;

%Initial cell arrays for pick times. The pick time array rows give (in
%order) for each station (column):
%(1) record view origin time
%(2) p-wave arrival time
%(3) s-wave arrival time
%(4) c-wave arrival time
%(5) record view ending time
%This cell array is given by tLoc.
tLoc    = cell(Nclus,1);
%chanLoc = tLoc;

for k = 1:Nclus
    
    S = Sc{k};
    
    %S           = coralExtract(Sc{k},'staChannel',opt.chan);
    %pickInd     = strncmpi({Sc{k}.staChannel},opt.chan,3);
    t0          = [S.recStartTime];
    ti          = t0;
    %[~,ind]     = min(timediff(t0));
    %t0          = t0(:,ind);
    cPick       = [S.cPick];
    phase       = {cPick.phase};
    %polarity    = {cPick.polarity};      
    indp        = strncmpi(phase, 'P', 1);
    inds        = strncmpi(phase, 'S', 1);
    indc        = strncmpi(phase, 'C', 1);
    pickDate    = [cPick.pickTime];
    pwave       = pickDate(:,indp);
    %pPol        = polarity(:,indp);
    swave       = pickDate(:,inds);
    %sPol        = polarity(:,inds);
    cwave       = pickDate(:,indc);
    ind         = repmat(logical([1,0]),1,size(S,1));
    cwave       = cwave(:,ind);
    
    % *pDiff gives the times after the record start time the phase is
    % observed. These give the time the pick arrow will be plotted.
    ppDiff      = timediff(pwave,ti);
    spDiff      = timediff(swave,ti);
    cpDiff      = timediff(cwave,ti);
    
    %ppDiff      = timediff(pwave);
    %start the plotting at this time, less the uncertainty
    %[~,ind]     = min(ppDiff);
    %dt          = (timediff([t0,pwave]));
    %t0          = min(dt(2:end));
    %tend        = timeadd(cwave,t0);
    %tend        = t0 + max(timediff(tend));   
    %tend        = max(cpDiff);    
    %tend        = timeadd(ti,max(cpDiff));
    t0          = repmat(min(ppDiff) - t_unc(1), 1, size(ppDiff,2));
    tend        = repmat(max(cpDiff) + t_unc(2), 1, size(cpDiff,2));
    tLoc{k}     = [t0; ppDiff; spDiff; cpDiff; tend];
    
end;

%Deconvolve instrument response and correct any errors in polarity
%Sc      = cellfun(@coralPolarityCheck,Sc,'UniformOutput',false);
%Sc      = cellfun(@coral2DeconDisp,Sc,repmat({opt},size(Sc)),'UniformOutput',false);
%Sc      = cellfun(@coralDetrend,Sc,'UniformOutput',false);

h_plot  = [];
ymarks  = 0:-1:(-size(Sc{1},1)+1);
chLoop  = 1:size(Sc{k},1);
%Plot all record sections. Plot arrows only the requested channel.

for k = 1:length(tLoc),
    
    figure;  
    h       = coralPlotLabel(Sc{k});
    ind     = strncmpi({Sc{k}.staChannel},opt.chan,3);
    
    for n = chLoop(ind),
                        
        if(strncmpi('phase',phLabel,4)),
            
            text(tLoc{k}(2,n),ymarks(n),'P','HorizontalAlignment','right',...
                'VerticalAlignment','bottom',...
                'color','k','fontsize',16,'fontweight','bold');
            text(tLoc{k}(3,n),ymarks(n),'S','HorizontalAlignment','right',...
                'VerticalAlignment','bottom',...
                'color','k','fontsize',16,'fontweight','bold');
            text(tLoc{k}(4,n),ymarks(n),'C','HorizontalAlignment','right',...
                'VerticalAlignment','bottom',...
                'color','k','fontsize',16,'fontweight','bold');
            
        end;
        
        %Plot pick arrow
        text(tLoc{k}(2,n),ymarks(n),'\downarrow',...
            'VerticalAlignment','bottom',...
            'color','k','fontsize',18,'fontweight','bold');
        text(tLoc{k}(3,n),ymarks(n),'\downarrow',...
            'VerticalAlignment','bottom',...
            'color','r','fontsize',18,'fontweight','bold');
        text(tLoc{k}(4,n),ymarks(n),'\downarrow',...
            'VerticalAlignment','bottom',...
            'color','b','fontsize',18,'fontweight','bold');
        h_plot  = cat(1,h_plot,h);
        
    end;
    
    xlim([median(tLoc{k}(1,:)),median(tLoc{k}(end,:))]);
    
end;

return;
