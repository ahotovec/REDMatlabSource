function [ DATA ] = getIRISdata(staList,startTime,endTime)
%Download data directly from the IRIS data center to a coral structure
%Any data that is in their database can be accessed. 

%INPUTS

%staList = list of stations, channels and network codes to access 
%    must be in this form:{'RCS.EHZ.UW','RCS.EHZ.UW','PANH.BHZ.CC'}
%startTime = start time in any date form you wish that matlab can 
%    interpret, vector is easiest 
%endTime = end time

%OUTPUT

%Data is a structure array containing coral structures for all requested
%traces. Any data that is not available will be filled with a coral
%structure of zeros

%example of usage
%DATA=getIRISdata({'RCS.EHZ.UW','STAR.EHZ.UW','RCM.EHZ.UW'},[2012 02 01 00 00 00],[2012 02 01 01 00 00]);

%set options for merging chopped up data
opt.time_tol=0.01; %time at the start of each input seismogram that
%      is merged must match time in merged seismogram to this 
%      tolerance (s) 
opt.fill_max=3000000; %if there is a gap between seismograms, merge them
%      anyway if the gap is less than opt.fill_max samples (keep big so
%      empty sections just get filled with zeros)
opt.fill_type=0; %fill with zeros

%go through list of stations
for i=1:length(staList)
emptydata=0; %marker to tell
%try to fetch trace from IRIS, if something goes wrong it will fail
try
m = irisFetch.Traces(staList{i}(end-1:end),staList{i}(1:end-7),'*',...
   staList{i}(end-5:end-3),datestr(startTime,'yyyy-mm-dd HH:MM:SS'),...
   datestr(endTime,'yyyy-mm-dd HH:MM:SS'));
catch %if data doesn't exist at IRIS, fill it in with zeros so your programs don't get messed up
    fprintf(['NO data at IRIS for ',staList{i},'! created placeholder structure\n'])
temp.staCode=staList{i}(1:end-7);
temp.staChannel=staList{i}(end-5:end-3);
temp.staNetworkCode=staList{i}(end-1:end);
temp.staLocationcode='';temp.staQualityCode=[];temp.staType=[];
temp.staLat=[];temp.staLon=[];temp.staElev=[];
temp.recSampInt = 0.01;
temp.recStartTime=startTime';
temp.recDip=[];temp.recAzimuth=[];temp.recNumData=[];
temp.data=zeros((datenum(endTime)-datenum(startTime))*24*60*60/temp.recSampInt,1);
%fill in empty stuff
temp.staRespType=[]; temp.staGain=[]; temp.staGainUnit=[]; temp.staGainFreq=[];
temp.staNormalization=[];temp.staNormalizationFrequency=[]; temp.staPoles=[];
temp.staZeros=[]; temp.staZeros=[]; temp.eqLat=[];temp.eqLon=[];temp.eqDepth=[];
temp.eqOriginTime=[];temp.eqMagnitude=[];temp.eqMagnitudeType=[];temp.eqMomentTensor=[];
temp.eqComment='';temp.recLog=[];temp.recComment='data accessed using irisFetch.Traces';
DATA(1,i)=temp;
clear m temp  
emptydata=1;
    continue %go to next station
end

if emptydata==0;
for j=1:length(m)
%change from IRIS format to coral structure
temp(j).staCode=m(j).station;
temp(j).staChannel=m(j).channel;
temp(j).staNetworkCode=m(j).network;
temp(j).staLocationcode=m(j).location;
temp(j).staQualityCode=m(j).quality;
temp(j).staType=[];
temp(j).staLat=m(j).latitude;
temp(j).staLon=m(j).longitude;
temp(j).staElev=m(j).elevation;    
temp(j).recSampInt=1./m(j).sampleRate;

temp(j).recStartTime=datevec(m(j).startTime)';
temp(j).recDip=m(j).dip;
temp(j).recAzimuth=m(j).azimuth;
temp(j).recNumData=m(j).sampleCount;
temp(j).data=double(m(j).data);

%fill in empty stuff
temp(j).staRespType=[]; temp(j).staGain=[]; temp(j).staGainUnit=[]; temp(j).staGainFreq=[];
temp(j).staNormalization=[];temp(j).staNormalizationFrequency=[]; temp(j).staPoles=[];
temp(j).staZeros=[]; temp(j).staZeros=[]; temp(j).eqLat=[];temp(j).eqLon=[];temp(j).eqDepth=[];
temp(j).eqOriginTime=[];temp(j).eqMagnitude=[];temp(j).eqMagnitudeType=[];temp(j).eqMomentTensor=[];
temp(j).eqComment='';temp(j).recLog=[];temp(j).recComment='data accessed using irisFetch.Traces';

end
end


DATA(1,i)=coralMerge(temp,opt);

clear m temp

end

