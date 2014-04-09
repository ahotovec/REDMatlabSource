function [varargout]=getDailyClusters(year,day,varargin)
%Takes ARRAYPICKS .mat files according to input year, day, channel, etc...
% and clusters array vector observation structures, according to a user
% specified threshold.
% (Kate A edited raPicks2raCluster to run two arraypicks files if there are two for
% any given day, renamed for simplicity)
%
% USAGE
% [varargout]=getClusters(year,day,opt)
%
% INPUT
% year:     The four digit year corresponding to the pick dates.
% day:      The three digit Julian day corresponding to the pick dates.
%
% opt:      An optional input structure with the following fields:
% thresh        The correlation threshold above which defines a cluster
% maxNumClust   The maximum number of observations to include in a cluster.
%               Used to avoid running out of memory
% cutTime       One-half the amount of time to cut a seismogram out
% fcut          The percentage of Nyquist to include pass band of filter
% chan          The channel of the station to cluster over
%
% OUTPUT
% Scluster:     The cell of coral structure clusters
% Sorphan:      The coral structure of orphans (non-clusters)
% ind:          An index vector giving the columns of the ARRAYPICK*.mat
%               file loaded internally.
% saveFile:     The name of the saved file
%-----------------------------------------------------------------------
% Latest Edit: 25.Nov.2009
% Joshua D Carmichael
% josh.carmichael@gmail.com
%
% Edit Log
% 01.June.2011
% Modified handling of too-large of ARRAYPICKS, changed input handling to
% be more efficient
%
% 07.Nov.2011
% Duplicate columns are removed from coral structure prior to running
% coralCluser.m. Outputs are corrected. Remember the saved Scluster and
% Sorphan may correlate below opt.thresh. This is only because coralCluster
% does filtering operations internally, and the correlation is only
% required to exceed rho on the correlated data. To remove this, choose
% opt.fcut = the Nyquist frequency (or 90% of it). Example, if opt.fcut =
% 0.35, and opt.recSampInt = 200, set opt.fcut = 0.9 instead.
%-----------------------------------------------------------------------

%intialize output
if(nargout > 1)
    
    varargout = cell(nargout,1);        
    
end;


%check input for opt input
ind = cellfun(@isstruct,varargin);

if(nnz(ind)>0),
    
    opt = varargin{ind};

else

    opt = defaultOptStruc;
    
end;

%Use getFiles to find ARRAYPICK files in current directory matching the
%date and year inputs
dayNum      = day;
day         = sprintf('%s.%.3i%s.','\',day,'\');

ind         = isfield(opt,'chan');

if(not(ind))
    
    opt.chan = '';
    
end;

%look for first file
flist       = getFiles('ARRAYPICKS',(opt.chan),num2str(day),num2str(year),'.1.mat');

if(length(flist)>1)
    
    ind     = regexp(flist,'ARRAYPICKS','start');
    ind     = cell2mat(ind);
    flist   = flist(ind(ind==1));
    
    if(length(flist)>1)
        
        error('Too many arraypick returns');
        
    end;
    
end;

arraypicks  = load(flist{1});
arraypicks  = arraypicks.arraypicks;
S           = arraypicks.coralstruc;
S=coralPad(S);

%Take array input and identify unique start times
if(size(S,2)>=2)
    
    %get the original indices for the columns of the input structure
%    indOrg      = 1:size(S,2);
    sepTimes    = diff(timediff([S(1,:).recStartTime,datevec(date)']));
    S           = S(:,abs(sepTimes)>(opt.cutTimeAF));
%    indOrg      = indOrg(abs(sepTimes)>(opt.cutTime));
    
%     if( min(abs(sepTimes)) < (opt.cutTimeAF) )
%         
%         disp(sprintf('Warning: Input structure contains duplicate columns.'));
%         disp(sprintf('Duplicate columns are now removed.'));
%         disp(sprintf('Outputs have been corrected.'));
%         
%     end;
    
end;

Ix          = 1:size(S,2);
S           = arrayfun(@coralDetrend,S);
S           = arrayfun(@coralTaper,S);
%S           = arrayfun(@coralPad,S);
Sorig     = S; %save copy of original data to save later

%filter before clustering if specified in defaultOptStruc
if length(opt.bpfilt)==2
    S=coralFilter(S,opt.bpfilt,'bandpass',2,'minimum');
elseif opt.lpfilt>0 && opt.bpfilt==0
    S=coralFilter(S,opt.lpfilt,'low',4,'minimum');
end

[~,~,ind1] 	= coralCluster(arrayfun(@coralNormData,S),opt);

Scluster1    = cell(size(ind1));

for k = 1:length(Scluster1),
    
    Scluster1{k} = Sorig(:,ind1{k});
    
end;

Sorphan1     = Sorig(:,setdiff(Ix,cell2mat(ind1)));
% Old Code Marker -------------------------------------------------------
% Code in coralClusters was previously in here
%-------------------------------------------------------------------------



%% Second file if it exists
clear S Sorig flist

%if size(dir(['ARRAYPICKS.',(opt.chan),'.',num2str(year),'.',num2str(day),'.2.mat']),1)~=0
if size(dir(['ARRAYPICKS.',(opt.chan),'.',num2str(year),'.*',num2str(dayNum,'%.3i'),'.2.mat']),1)~=0
fprintf('There are two arraypicks files, picking second one\n')
flist       = getFiles('ARRAYPICKS',(opt.chan),day,num2str(year),'.2.mat');

if(length(flist)>1)
    
    ind     = regexp(flist,'ARRAYPICKS','start');
    ind     = cell2mat(ind);
    flist   = flist(ind(ind==1));
    
    if(length(flist)>1)
        
        error('Too many arraypick returns');
        
    end;
    
end;

arraypicks  = load(flist{1});
arraypicks  = arraypicks.arraypicks;
S           = arraypicks.coralstruc;
S=coralPad(S);

%Take array input and identify unique start times
if(size(S,2)>=2)
    
    %get the original indices for the columns of the input structure
%    indOrg      = 1:size(S,2);
    sepTimes    = diff(timediff([S(1,:).recStartTime,datevec(date)']));
    S           = S(:,abs(sepTimes)>(opt.cutTimeAF));
%    indOrg      = indOrg(abs(sepTimes)>(opt.cutTime));
    
    if( min(abs(sepTimes)) < (opt.cutTimeAF) )
        
        disp(sprintf('Warning: Input structure contains duplicate columns.'));
        disp(sprintf('Duplicate columns are now removed.'));
        disp(sprintf('Outputs have been corrected.'));
        
    end;
    
end;

Ix          = 1:size(S,2);
S           = arrayfun(@coralDetrend,S);
S           = arrayfun(@coralTaper,S);
%S           = arrayfun(@coralPad,S);
Sorig     = S; %save copy of original data to save later

%filter before clustering if specified in defaultOptStruc
if length(opt.bpfilt)==2
    S=coralFilter(S,opt.bpfilt,'bandpass',2,'minimum');
elseif opt.lpfilt>0 && opt.bpfilt==0
    S=coralFilter(S,opt.lpfilt,'low',4,'minimum');
end

[~,~,ind2] 	= coralCluster(arrayfun(@coralNormData,S),opt);

Scluster2   = cell(size(ind2));

for k = 1:length(Scluster2),
    
    Scluster2{k} = Sorig(:,ind2{k});
    
end;

Sorphan2     = Sorig(:,setdiff(Ix,cell2mat(ind2)));
% Old Code Marker -------------------------------------------------------
% Code in coralClusters was previously in here
%-------------------------------------------------------------------------

%Combine first and second file outputs
Scluster=[Scluster1; Scluster2];
Sorphan=[Sorphan1 Sorphan2];
ind=[ind1; ind2];
else
    fprintf('just one file, fool\n')
    Scluster=Scluster1;
    Sorphan=Sorphan1;
    ind=ind1;
end

%Create a file-to-save name for the cluster picks
saveFile    = sprintf('%s.%s.%.3i.%s.%i.%.3i.mat','CLUSTERPICKS',...
                      'ARRAY',round(100*(opt.thresh)),(opt.chan),year,dayNum);
disp(saveFile)                 
save(saveFile,'Scluster','Sorphan','ind');

if(nargout>=1)
    
    varargout{1}=Scluster;
    
end;

if(nargout>=2)
    
    varargout{2}=Sorphan;
    
end;

if(nargout>=3)
    
    varargout{3}=ind;
    
end;

if(nargout>=4)
    
    varargout{4}=saveFile;
    
end