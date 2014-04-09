function [varargout]=raPicks2raCluster(year,day,varargin)
% Takes ARRAYPICKS .mat files according to input year, day, channel, etc...
% and clusters array vector observation structures, according to a user
% specified threshold.
%
% USAGE
% [varargout]=raPicks2raCluster(year,day,opt)
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

flist       = getFiles('ARRAYPICKS',(opt.chan),day,num2str(year),'.mat');

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

[~,~,ind] 	= coralCluster(arrayfun(@coralNormData,S),opt);

Scluster    = cell(size(ind));

for k = 1:length(Scluster),
    
    Scluster{k} = S(:,ind{k});
    
end;

Sorphan     = S(:,setdiff(Ix,cell2mat(ind)));
%% Old Code Marker -------------------------------------------------------
% Code in coralClusters was previously in here
%-------------------------------------------------------------------------

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
    
end;

%% Old Code Material
% Code below has been incorporated into coralCluster
%-------------------------------------------------------------------------
% opt.thresh      = 0.65;
% opt.maxNumClust = 500;
% opt.cutTime     = 5;
% opt.fcut        = 0.35;
% opt.chan        = 'EPZ';
% Scluster        = [];
% Sorphan         = [];
% 
% %intialize output
% if(nargout > 1)
%     for k = 1:nargout-1;
%         varargout{k} = [];
%     end;
% end;
% 
% %check input
% for k = 1:length(varargin)
%     if(isstruct(varargin{k}))
%         opt = varargin{k};
%     end;
% end;
% 
% %Take array input and identify unique start times
% if(size(S,2)>=2)
%     sepTimes    = diff(timediff([S(1,:).recStartTime,datevec(date)']));
%     S           = S(:,sepTimes>(opt.cutTime));
% end;
% 
% %cut out the station structure in Scut, and normalize the data for
% %processing and clustering
% Scut        = S;
% Scut        = arrayfun(@coralNormData,Scut);
% sintr       = median([Scut.recSampInt],2);
% Fnyq        = 0.5*1/(mean(sintr));
% 
% %Low pass filter the data for correlating and clustering
% for m = 1:size(Scut,2)
%     Scut(:,m)   = coralFilter(Scut(:,m),(opt.fcut)*Fnyq,'low',8,'zero');
% end;
% 
% %Perform mutually exclusive partition based clustering
% if( size(Scut,2) < opt.maxNumClust )
%     
%     [Scluster,Sorphan] = xcorMutClust(Scut,opt);
%     
% else
%     
%     count=0;
%     
%     while(size(Scut,2)>=1)
%         
%         [Sc,So]  = xcorMutClust(Scut(:,1:min(opt.maxNumClust,size(Scut,2))),opt);
%         
%         if(count==0), Scluster = Sc; Sorphan = So; end;
%         
%         if(count > 0)
%             Scluster           = cat(1,Scluster,Sc);
%             Sorphan            = cat(2,Sorphan,So);
%         end;
%         
%         count=count+1;
%         Scut(:,1:min(opt.maxNumClust,size(Scut,2))) = [];
%         if(isempty(Scut)), break, end;
%         
%     end;
%     
% end;
% 
% if(nargout>=1)
%     varargout{1}=Scluster;
% end;
% if(nargout>=2)
%     varargout{2}=Sorphan;
% end;