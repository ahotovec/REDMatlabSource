% Script to run multiplet analysis using predefined templates
clear all
close all
clc

%define some things
threshold=0.5; %threshold correlation coefficient to put an event in a cluster
maxtime=60; %maximum number of days between day being searched and templates being compared to the data

%add path to all the scripts required (edit to match their location on your
%own computer
addpath C:\Users\geo-user\Documents\KateA\MultCodes
addpath C:\Users\geo-user\Documents\KateA\MultCodes\MultipletPrograms
%javaaddpath /Users/Kate/MATLAB1/Kates/usgs.jar
javaaddpath C:\Users\geo-user\Documents\KateA\MultCodes\IRIS-WS-1.3.jar

cd C:\Users\geo-user\Documents\KateA\MultCodes

%load in template
load C:\Users\geo-user\Documents\KateA\MultCodes\MasterStacks_RCS.mat


%enter name of the station to use in the analysis in the form
%stations={'STAR.EHZ.UW'}; only one station can be used, should be same
%station as the templates are
station={'RCS.EHZ.UW'};

%enter starting month and year
startYear=2008;
startMonth=1; %number
endYear=2013;
endMonth=1;

%make vector of starts of each month in datenum format
monthvec=[];
for j=startYear:endYear
    if j==startYear && j~=endYear
        addl=[];
        for k=startMonth:12;
            addl=[addl datenum([j k 1 0 0 0])];
        end
    elseif j==endYear && j~=startYear
        addl=[];
        for k=1:endMonth
            addl=[addl datenum([j k 1 0 0 0])];
        end
    elseif j==endYear && j==startYear
        addl=[];
        for k=startMonth:endMonth
            addl=[addl datenum([j k 1 0 0 0])];
        end
    else
        addl=[];
        for k=1:12
            addl=[addl datenum([j k 1 0 0 0])];
        end
    end
    monthvec=[monthvec addl];
end
%% p=1;
% if p==1
%matlabpool(4)

for i=1:length(monthvec)-1
    %if month has already been done, skip to next one
    if isempty(dir([datestr(monthvec(i),'yyyy-mm_'),MasterStack(1).staCode,'.mat']));
        
        
        %initiate empty cell array
        D=cell(size(MasterStack));
        
        %change back to home directory if not already there
        cd C:\Users\geo-user\Documents\KateA\MultCodes\RCStemplsearch
        
        %enter starttime and endtime for the dates you want to search over, use
        %at least half days
        startTime=datevec(monthvec(i));
        endTime=datevec(monthvec(i+1));
        
        %Run a template search
        
        %first define which templates can be compared for this time period (within
        %maxtime days)
        I=find(startdates>datenum(startTime)-maxtime & enddates<datenum(endTime)+maxtime);
        templates=MasterStack(I);
        %compare them and pull out data
        tic;
        %[allEv]=coralTemplateSearchMany(templates,datenum(startTime),datenum(startTime)+0.01,threshold);%outputs a cell array of same number of columns as MasterStack containing
        [allEv]=coralTemplateSearch(templates,startTime,endTime,threshold);%outputs a cell array of same number of columns as MasterStack containing
        toc;
        
        %place these events in a master cell array in the right locations for their
        %corresponding templates
        for m=1:length(I)
            D{I(m)}=[D{I(m)} allEv{m}];
        end
        
        %remove duplicate events (keep the event with the higher xcorval)
        [D]=coralRemoveDuplicates(D,5);
        
        saveTemplSearch(D,monthvec(i),templates(1).staCode)
        %save([datestr(monthvec(i),'yyyy-mm_'),templates(1).staCode,'.mat'],'D');
    else
        fprintf([datestr(monthvec(i),'yyyy-mm_'),' already done\n']);
    end
end
%matlabpool close