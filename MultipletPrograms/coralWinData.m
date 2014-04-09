function D1 = coralWinData(staList, startTime, endTime, server, port)
% coralWinData      Get data from winston wave server.
%   D1 = coralWinData(staList, startTime, endTime, server, port);
%
%   coralWinData queries a winston wave server and returns data in a
%   coral structure.  All code is self contained except for usgs.gov, which
%   is located in several places.  If the code runs, its already in your
%   path, if it doesn't then there are directions on where to find it in 
%   the error messages that result.
% Input:
%   staList: cell array of station,channel and network code separated by
%            periods.  Will also accept location codes, though it isn't
%            required.  You may also indicate one of the following strings
%            as the first cell in staList to get presets:
%               msh
%               rainier
%               sisters
%               hood
%   startTime: row vector of start time (matlab format)
%   endTime: row vector of end time (matlab format)
%   server: name of computer with winston. The following servers are
%	    preconfigured: 
%			tremito
%			pele 
%			eruption
%			jeannie
%           mazama
%           If you want something else, make the server
%           the full network name or ip and add the port in the next variable.
%   port(optional): port of winston wave server. If one of the preconfigured
%                   machines are not specified, then this variable must 
%                   be assigned or it will be set the the default value of 16017.
% Output:
%   D1: coral structure of waveforms with station information (if available)
%
% IMPORTANT: IF USING A MAC, THE JAVA DEFAULTS ARE SCREWED UP.  ADD THE
% FILE JAVA.OPTS (IN /U0/IRIS/MATLAB/TEST) TO THE DIRECTORY:
% /APPLICATIONS/MATLAB_<RELEASE>.APP 
% WHERE RELEASE IS THE MATLAB VERSION. RIGHT CLICK THE APPLICATION AND
% CLICK "SHOW CONTENTS" TO ACCESS THE DIRECTORY. SEE DETAILS AT MATHWORKS
% WEBSITE.  http://www.mathworks.com/support/solutions/en/data/1-18I2C/ 
%
% See also lookWinData.m
%
% Written by Weston Thelen 11/14/09
% thelenwes@gmail.com
% Based on code from the waveform and correlation package (see
% getwinstondata.m)
% Modified 1/25/2012 by Ken Creager to accept new servers
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% Default parameters-------------------------------------------------------
if nargin<5;    % if not specified, set default value for port number
  port = 16017;
end

if strcmp( server, 'tremito' )
    server = 'tremito.ess.washington.edu';
    port = 16017;
elseif strcmp( server, 'eruption' )
    server = 'eruption.ess.washington.edu';
    port = 16022;
elseif strcmp( server, 'pele' )
    server = 'pele.ess.washington.edu';
    port = 16017;
elseif strcmp( server, 'jeannie' )
    server = 'jeannie.ess.washington.edu';
    port = 16017;
elseif strcmp( server, 'importnq' )
    server = 'importnq.ess.washington.edu';
    port = 16017;
elseif strcmp( server, 'avo' )
    server = 'pubavo1.wr.usgs.gov';
    port = 16022;
elseif strcmp( server, 'mazama')
    server= 'mazama.ess.washington.edu';
    port = 16017;
end

verbose=0;

% Check for column vectors
if size(startTime,1) ~= 1
    startTime = startTime';
end
if size(endTime,1) ~= 1
    endTime = endTime';
end

t1 = mep2dep(datenum(startTime));

t2 = mep2dep(datenum(endTime));

%keyboard;
if t1 > t2, error('Start Time > End Time');  end

%--------------------------------------------------------------------------
%Preset datasets ----------------------------------------------------------
if strcmp(staList{1}, 'msh')
    staList = {'SEP.EHZ.CC', 'NED.EHZ.CC', 'VALT.BHZ.CC', 'SUG.EHZ.UW',...
        'STD.BHZ.CC', 'B202.EH1.PB', 'HSR.EHZ.UW', 'SHW.EHZ.UW',...
        'EDM.EHZ.UW', 'JUN.EHZ.UW', 'B204.EH1.PB', 'B203.EH1.PB', 'CDF.EHZ.UW', 'FL2.EHZ.UW',...
        'B201.EH1.PB', 'ELK.EHZ.UW', 'TDL.EHZ.UW', 'LVP.EHZ.UW'};
elseif strcmp(staList{1}, 'rainier')
    staList = {'RCM.EHZ.UW', 'RCS.EHZ.UW', 'STAR.EHZ.UW', 'PANH.BHZ.CC', 'FMW.EHZ.UW',...
        'RER.EHZ.UW', 'RVC.EHZ.UW', 'LON.BHZ.UW', 'WPW.EHZ.UW',...
        'GLK.EHZ.UW'};
elseif strcmp(staList{1}, 'sisters')
    staList = {'HUO.EHZ.UO', 'MOON.EHZ.UW', 'WIFE.BHZ.CC',...
        'PRLK.BHZ.CC', 'BKC.EHZ.UW', 'FRIS.EHZ.U0'};
elseif strcmp(staList{1}, 'hood')
    staList = {'TIMB.EHZ.CC', 'HOOD.BHZ.UW', 'TDH.EHZ.UW',...
        'VFP.EHZ.UW', 'VLL.EHZ.UW'};
elseif strcmp(staList{1}, 'redoubt')
    staList = {'RSO.EHZ.AV', 'REF.EHZ.AV', 'RDN.EHZ.AV'};
end

%-----------------------------------------------------------------------
% Open winston -----------------------------------------------------------
try
    WWS=gov.usgs.winston.server.WWSClient(server, port);
catch
    disp('Getting usgs.jar into your path');
    getLib();
    try
        WWS=gov.usgs.winston.server.WWSClient(server, port);
    catch
        disp('usgs.jar not found.  Go find it. ');
    end
end
rr= WWS.connect;
if rr == 0
    error('\nDid not connect to winston.  Check server and port settings.');
end

%-----------------------------------------------------------------------
% Get station information
p = WWS.getChannels.toArray;
for i = 1 : length(p)
    codes{i} = char(p(i).getCode);
    try
        stTime(i) = datenum(sec2cal(p(i).getMinTime));
    catch
        stTime(i) = NaN;
    end
    try
        enTime(i) = datenum(sec2cal(p(i).getMaxTime));
    catch
        enTime(i) = NaN;
    end
end

%-----------------------------------------------------------------------
% Cycle through stations in station list
cntr = 1;
for i = 1 : length(staList)
    t1 = mep2dep(datenum(startTime));
    %---------------------------------------------------------------------
    % Setup coral stuff --------------------------------------------------
    D(cntr) = getEmptyCoralStruct;
    [sta, chan, net, loc] = strread(staList{i}, '%s%s%s%s', 'delimiter', '.');
    if isempty(net)
        net = {'UW'};
    end
    if isempty(loc)
        loc = {'--'};
    end
    D(cntr).staCode = sta{1};
    D(cntr).staChannel = chan{1};
    D(cntr).staNetworkCode = net{1};
    if verbose; disp(sprintf('\nGetting %s.%s.%s.%s from %s', sta{1}, chan{1},...
        net{1}, loc{1}, server)); end
    %---------------------------------------------------------------------
    % Start Winston Stuff ------------------------------------------------
    %keyboard;
    % Check that starttime is during data availability
    [sta, chan, net] = strread(staList{i}, '%s%s%s', 'delimiter', '.');
    staInfo = strcat(sta{1}, '$', chan{1}, '$', net{1});
    R = find(strcmp(staInfo, codes)==1);
    if datenum(startTime) > stTime(R)
        pad = 0;
        if datenum(endTime) > enTime(R)
             warning('End time for %s.%s.%s exceeds time in winston', sta{1}, chan{1}, net{1}); 
        end
    elseif datenum(startTime) < stTime(R)
        tpadS = timediff(datevec(stTime(R))', startTime');
        t1 = mep2dep(datenum(stTime(R)));
        if t1 > t2
            warning(sprintf('Check start time for %s.%s.%s', sta{1}, chan{1}, net{1}));
        else
            warning(sprintf('Zero padding %s.%s.%s because startTime is out of range!',...
                sta{1}, chan{1}, net{1}));
        end
        pad = 1;
    end
    
    tr = WWS.getRawData(sta{1}, chan{1}, net{1}, loc{1}, t1, t2);
    %keyboard;
    if ~isempty(tr)
        tr.buffer(tr.buffer == intmin('int32')) = 0; %fix odd spikes
        D(cntr).data = double(tr.buffer); % d.buffer is int32, and must be converted
        if verbose; 
          disp(sprintf('Start Time: %s', datestr(startTime, 'mm/dd/yy HH:MM:SS')));
          disp(sprintf('End Time: %s', datestr(endTime, 'mm/dd/yy HH:MM:SS')));
        end
        D(cntr).recSampInt = 1/tr.getSamplingRate;
        D(cntr).recNumData = tr.samples;
        D(cntr).recStartTime = datevec(dep2mep(tr.getStartTime))';
        D(cntr).recComment = sprintf('coralWinData:%s:%d', server, port);
        if pad == 1
            %keyboard;
            nzero = round((1/D(cntr).recSampInt)*tpadS+1); % The plus one is for cutting purposes (good to have a tail)
            padVec = zeros(nzero,1);
            D(cntr).recStartTime = timeadd(D(cntr).recStartTime, -nzero*D(cntr).recSampInt);
            D(cntr).data = [padVec; D(cntr).data];
            D(cntr).recNumData = length(D(cntr).data);
        end
    else 
        disp('-----------------------------------------------------------');
        disp('-----------------------------------------------------------');
        disp(sprintf('%s.%s.%s.%s may not exist or data may not be available', sta{1}, chan{1},...
            net{1}, loc{1}));
        disp('-----------------------------------------------------------');
        disp('-----------------------------------------------------------');
    end
    %---------------------------------------------------------------------
    %---------------------------------------------------------------------
    % Put in station location, if in database-----------------------------
    if strcmp(loc{1}, '--')
        staInfo = strcat(sta{1},'$', chan{1}, '$', net{1});
    else
        staInfo = strcat(sta{1},'$', chan{1}, '$', net{1}, '$', loc{1});
    end
    R = find(strcmp(staInfo, codes)==1);
%     if ~isempty(R)
%         lat = p(R).getLat;
%         if ~isnan(lat)
%             D(cntr).staLat = lat;
%             D(cntr).staLon = p(R).getLon;
%         end
%     end
    %---------------------------------------------------------------------
    cntr = cntr + 1;
end
WWS.close;
%keyboard;
opt.cutType = 'absTime';
opt.absStartTime = startTime';
opt.absEndTime = endTime';
for i = 1 : length(D)
    if ~isempty(D(i).data)
        [D1(i), ierr] = coralCut(D(i), opt);
    else
        D1(i) = D(i);
    end
end


function D = getEmptyCoralStruct
%-------------------------------------------------------------------------
% getEmptyCoralStruct
% Function to return an empty coral structure.
% Usage:
%   D = getEmptyCoralStruct
% Output:
%   D: Empty coral structure
%-------------------------------------------------------------------------


D.staCode = '';
D.staChannel = '';
D.staNetworkCode = '';
D.staLocationCode = '';
D.staQualityCode = '';
D.staType = '';
D.staLat = [];
D.staLon = [];
D.staElev = [];
D.staRespType = '';
D.staGain = [];
D.staGainUnit = '';
D.staGainFreq = [];
D.staNormalization = [];
D.staNormalizationFreq = [];
D.staPoles = [];
D.staZeros = [];
D.eqLat = [];
D.eqLon = [];
D.eqDepth = [];
D.eqOriginTime = [];
D.eqMagnitude = [];
D.eqMagnitudeType = {};
D.eqMomentTensor = [];
D.eqComment = '';
D.recSampInt = [];
D.recStartTime = [];
D.recDip = [];
D.recLog = '';
D.recAzimuth = [];
D.recNumData = 0;
D.data = [];
D.recComment = '';

function matTime = mep2dep(matepoch)
matTime = (matepoch - 719529) * 24 * 3600;
%matTime = (matepoch - 730486.5) * 24 * 3600;

function depTime = dep2mep(dep)
depTime = dep / 86400 + 719529;
%depTime = dep / 86400 + 730486.5;

function getLib()
warning off backtrace
jcp = javaclasspath('-all');
RequiredFiles = 'usgs.jar';
% Go through a number of possibilities for finding usgs.jar
if exist(RequiredFiles, 'file') == 2
    javaaddpath(which('usgs.jar'));
end
if exist(RequiredFiles, 'file') ~= 2
    javaaddpath('http://grasso.ess.washington.edu/LOCAL/archive/wethelen/usgs.jar');
end
if exist(RequiredFiles, 'file') ~= 2
    javaaddpath('./usgs.jar');
end
if exist(RequiredFiles, 'file') ~= 2
    javaaddpath('/opt/local/swarm2/usgs.jar');
end
if exist(RequiredFiles, 'file') ~= 2
    javaaddpath('/u0/iris/MATLAB/TEST/usgs.jar');
end
if exist(RequiredFiles, 'file') ~= 2
    javaaddpath('http://www.avo.alaska.edu/Input/celso/swarmstuff/usgs.jar');
end
if exist(RequiredFiles, 'file') ~= 2
    disp('You must have a copy of usgs.jar in your path.');
    disp('It is in several locations, try:');
    disp('/u0/iris/MATLAB/TEST');
    disp('/opt/local/swarm2 on a pnsn machine');
    disp('http://www.avo.alaska.edu/Input/celso/swarmstuff/usgs.jar');
end
warning on backtrace

function cal=sec2cal(sec)
%SEC2CAL  Returns calendar date given seconds before 2000/01/01 12:00.
%     cal=sec2cal(sec)

%This version is slightly slower but more numerically stable than
%the original: cal=datevec(sec/86400+730486.5)

%Make sure sec is a column vector

    sec=sec(:);
    
%Calculate calendar day

    fsec=floor(sec);
    epsilon=sec-fsec;
    
    fsec=fsec+63114033600;
    rs=mod(fsec,86400);
    rh=mod(rs,3600);
    rm=mod(rh,60);
    
    cal=datevec((fsec-rs)/86400);
    cal(:,4)=(rs-rh)/3600;
    cal(:,5)=(rh-rm)/60;
    cal(:,6)=rm+round(epsilon*1e8)/1e8;

function sec=cal2sec(cal)
%CAL2SEC   sec=cal2sec(cal)  
%
%Returns seconds before 2000/01/01 12:00 given a calendar date.
%Input 'cal' should be: [year month day hour minute second]
%(i.e., the dates are stored in the rows of 'cal'.)


%cal(:,1)=yy2yyyy(cal(:,1));

n=size(cal,2);

switch size(cal,2)

    case 1

        cal(:,2)=1;
        cal(:,3)=1;

    case 2
    
        cal(:,3)=1;


end

cal(:,7)=0;

sec=(datenum(cal(:,1),cal(:,2),cal(:,3),0,0,0)-730486.5)*86400+cal(:,4)*3600+cal(:,5)*60+cal(:,6);
%sec=(datenum(cal(:,1),cal(:,2),cal(:,3),0,0,0)-719529)*86400+cal(:,4)*3600+cal(:,5)*60+cal(:,6);
