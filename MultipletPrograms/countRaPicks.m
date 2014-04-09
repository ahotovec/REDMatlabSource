function [varargout] = countRaPicks(opt,varargin)
% Returns number of arraywide seismic picks. Does not return picks from
% single stations.
%
% USAGES
% [S] = countRaPicks(opt);
% [S] = countRaPicks(opt,year);
% [S] = countRaPicks(opt,day,year);
%
% INPUT
% tbin:     The time (sec) to bin observations
% optional: Any set of file numerics that further specify ARRAYPICKS.*.mat
%           files, like year, or day.
%
% OUTPUT
% S:        The output structure that contains observation metrics as 
%           determined by ARRAYPICKS.*.mat. Output S contains the
%           following fields:
%          
% day       A vector of observation days over which the
%           number of array picks were counted.
% tbin      A vector of bin locations over which the total
%           number of array pick events are binned.
% nPicksDay A vector of the number different of array picks
%           observed each day corresponding to vector day
% nPicksBin The binned set of array picks corresponding
%           to tbin
%-----------------------------------------------------------------------
% Latest Edit: 14.March.2010
% Joshua D Carmichael
% josh.carmichael@gmail.com
%
% Edit Log
% 14.March.2010
% Designed function
%
% 04.April.2010
% Changed input arguments so that numerical day or year inputs are
% accepted, in preference to year or day string inputs.
%-----------------------------------------------------------------------

if(nargout<1)
    
    plotopt = 'plot';
    
else
    
    plotopt = 'none';
    
end;

for k=1:nargout
    varargout{k} = [];
end;

L = length(varargin);

if(L==0)
    flist       = getFiles('ARRAYPICKS','.mat',(opt.chan));
end;

if(L==1)
    flist       = getFiles('ARRAYPICKS','.mat',(opt.chan),num2str(varargin{1}));
end;
if(L==2)
    flist       = getFiles('ARRAYPICKS','.mat',(opt.chan),num2str(varargin{1}),num2str(varargin{2}));
end;

%Initialize Output Structure
Tall        = [];
Tdayb       = [];
Tdaye       = [];
Nc          = [];        

for k = 1:length(flist)    
    
    %Load files, extra cluster and orphans
    loadOut     = load(flist{k});
    
    if(not(isfield(loadOut,'arraypicks'))), continue; end;
    
    raPick      = loadOut.arraypicks;
    raPick      = raPick.pTime;
    
    if(numel(raPick)==0), continue; end;
    
    Nc          = cat(1,Nc,size(raPick,2));
    
    temp        = textscan(flist{k},'%s','Delimiter','.');
    day         = str2double(temp{1}(end-1));
    year        = str2double(temp{1}(end-2));
    daystart    = [year,1,1,0,0,0]';
    daystart    = timeadd(daystart,(day-1)*24*3600 );
    dayend      = timeadd(daystart,24*3600);
    Tdayb       = cat(1,Tdayb,datenum(daystart'));
    Tdaye       = cat(1,Tdaye,datenum(dayend'));
    
    Tall        = cat(1,Tall,datenum(raPick'));
    
end;

tmin    = datevec( min(Tdayb) )';
tmax    = datevec( max(Tdaye) )';
Tobs    = timediff(tmax,tmin);
tbin    = [ 0.5/(opt.tbin):(opt.tbin):Tobs-0.5/(opt.tbin)]';
tbin    = datenum(timeadd(tmin,tbin')');

if(strcmp(plotopt,'plot'))

    hist(Tall,tbin); hold on;
    plot(Tdayb,Nc,'-r');
    datetick('x');
    
end;

if(nargout>0)
    
    count.day       = Tdayb;
    count.tbin      = tbin;
    count.nPicksDay = Nc;
    count.nPicksBin = Tall;
    varargout{1}    = count;
    
end;
return;