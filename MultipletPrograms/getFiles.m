function [varargout] = getFiles(varargin)
% [Sacs] = getFiles(varargin);
% This function returns a cell of file names from the current directory,
% given any number of inputs of regular expressions. It is analogous to the
% LINUX 'ls' command, but the output is given to the user as a cell, when 
% using getFiles.m
% 
% USAGE
% [Sacs] = getFiles('REGEXP1','REGEXP2',...,'REGEXPN');
%
% INPUT
%
% REGEXP1:  A character regular expression, such as 'get' or 'mat'.
%
% OUTPUT
% 
% Sacs:     A cell array containing the file names that match the
%           requested character matches.
%
% EXAMPLES
% 
% %Suppose you have a directory where JPEG images are stored under names
% %that include 'November' somewhere in the title. This returns those file
% %names:
%
% [Sacs] = getFiles('November','.jpg');
%
% %Suppose you have directory containing data of the form:
% %STATION.CHAN.YEAR.DAY.HR.SAC
% %The following command returns all .SAC files for station 'JAN', channel
% %'EPZ', Julian day 340, and year 2009:
%
% [Sacs] = getFiles('SAC','JAN',340,2009);
%
%-----------------------------------------------------------------------
% Latest Edit: 14.Nov.2009
% Joshua D Carmichael
% josh.carmichael@gmail.com
%-----------------------------------------------------------------------

% Start function
D = dir;
I = not([D.isdir]);
D = D(I);
C = struct2cell(D);
C = {C{1,:}};
Sacs    = C(:); 

if(nargin==0)
    
    return;

else

    for k = 1:length(varargin)
        
        rec     = varargin{k};
        
        if(~isempty(rec))
            
            I       = regexp(C,rec);
            full    = cellfun(@length,I);
            I       = logical(full);
            C       = C(I);
            
            Sacs    = intersect(Sacs,C(:));
            
        end;
        
    end;
    
end;
    
Sacs = unique(Sacs);
L    = size(Sacs,1);

if(L==0), 
    
    disp(sprintf('\nNo files found in current directory')), 
    
    for k=1:nargout
        
        varargout{k} = [];
        
    end;
    
end;

if(nargout>=1), varargout{1} = Sacs;    end;
if(nargout>=2), varargout{2} = L;       end;