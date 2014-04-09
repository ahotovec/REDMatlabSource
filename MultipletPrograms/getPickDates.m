function [varargout]=getPickDates(day,year,varargin)
% Returns observed pick dates from current directory, given Julian day,
% year, and another extension option, such as channel.
%
% USAGES
% [datecell] = getPickDates(day,year);
% [datecell] = getPickDates(day,year,ext);
% [datecell,files] = getPickDates(day, year);
% [datecell,files] = getPickDates(day, year, ext);
% [datecell,files,names] = getPickDates(day, year);
% [datecell,files,names] = getPickDates(day, year, ext);
%
% INPUT
% day:  Julian day that appears in PICKDATES file name
% year: Year that appears in PICKDATES file name
% ext:  Optional string that specifies which files to select for PICKDATES
%
% OUTPUT
% datecell: A cell of pick times for each station match for the input day 
%           and year
% files:    The files found for the requested pick dates
% names:    The station names corresponding to the datcells
%------------------------------------------------------------------------
% Latest Edit: 25.Nov.2009
% Joshua D Carmichael
% josh.carmichael@gmail.com
%------------------------------------------------------------------------

ext   = '.mat';

for k=1:nargout
    varargout{k} = [];
end;

for k = 1:length(varargin)
    
    if(ischar(varargin{k}))
        ext = varargin{k};
    end;

end;

day         = sprintf('%s.%.3i%s.','\',day,'\');
flist       = getFiles('PICKDATES',day,num2str(year),'.mat',ext);
fn          = length(flist);

% Determine the number of stations in the array:  Loop through flist
% files, find which files have a unique station name, and catenate the list
% of station names
for k = 1:fn
    
    if(k==1),
        temp    = textscan(flist{k}, '%s', 'Delimiter','.');
        names   = temp{1}(2);
        continue;
    else
        temp    = textscan(flist{k}, '%s', 'Delimiter','.');
        temp    = temp{1}(2);
        names   = unique( vertcat(temp, names) );
    end;
    
end;

if(fn==0),  return; end;
LS          = length(names);
datecell    = cell(1,1,LS);

for n = 1:fn
    
    load(flist{n});
    
    tcell   = pTime;
    
    L       = cellfun(@length,tcell);
    tcell   = tcell(L>1);
    t1      = [tcell{:}];
    
    datecell{1,1,n} = t1;
    
end;

if(nargout>=1), varargout{1} = datecell;	end;
if(nargout>=2), varargout{2} = flist;       end;
if(nargout>=3), varargout{3} = names;       end;

return;