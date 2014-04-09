%script_obliqProjCompare.m
% HEY KATE!
disp(sprintf('HEY KATE! Read the script is commented if you care.\n'));
disp(sprintf('It shows the power of Oblique Projectors at keeping\n'));
disp(sprintf('Good data, and rejected Bad data, even if there is\n'));
disp(sprintf('significant correlation between the good and the bad.\n'));
%This script compares simliar waveforms/signals using an oblique projector
%matrix. This demo script shows how one waveform type can be rejected from
%a detection statistic, while a similar, but different "enough" waveform
%that is "real" is kept in full, undistorted form. This script SHOULD be
%generalized so that it can be implemented as a more efficient detector,
%and used ARRAY based correlation. There are some time/computational
%savings that I can 
clc; close all;
%(1) Build 2 electronic spikes: both are electronic 'blips' that we want to
%ignore.
%

%(2) Build 2 real signals: both of these waveforms are real and we want to
%keep them, while ignoring/rejecting the electronic spikes in (1).

%change into the directory where the following multiplet records are stored
load MULTTEMPL.EHZ.070.2006.033.2006.059.mat

%Now unpack the multiplet templates and construct a set of template signals
%where there is "garbage" detections and "real" detections. These garbage
%and real detections will respectively be the NULL space for an OBLIQUE
%projector matrix and the RANGE of the OBLIQUE projector matrix.

%Indices for the good signals: save these
rangeInd    = [1;2;5;14;76;109;174;(201:205)';(209:211)'];

%Indices for the garbage signals: reject these
nullInd     = setdiff(1:211,rangeInd);

D           = [multtempl.template];
Dr          = [D{rangeInd}];
Dn          = [D{nullInd}];
D           = [Dr,Dn];

%NOW consider the WORSE case scenario: when all the signals (GOOD AND BAD)
%are aligned to maximize total correlation.
[C,I]   = raCrossCorr(D);
ind     = maxXcorrSum(C);
I       = triu(I,1) - triu(I,1)';
D       = coralCircShift(D,I(:,ind));

%Cut out the Good and Bad comparison part of the correlation matrix.
C       = C(1:length(rangeInd),length(rangeInd)+1:end);

%Make a plot of the signals to reject, and which ones to keep.
scrsz = get(0,'ScreenSize');
figure('Position',[1 scrsz(4) scrsz(3)/2 scrsz(4)]);
t = 0:D(1).recSampInt:(D(1).recSampInt)*(D(1).recNumData-1);
h = plotColumns(t,[D.data]); hold on;
set(gca,'fontsize',18,'fontweight','b','linewidth',2);
set(gcf,'color','w');
set(h(1:length(rangeInd)),'color','k');
set(h(length(rangeInd)+1:length(rangeInd)+length(nullInd)),'color','r');
title('Signals to KEEP are BLACK/REJECT are RED; BEFORE PROJECTION');

%Now use Josh Carmichael function to compute oblique projector that keeps
%good signals, and rejects crappy ones (so bad ones are in null space of
%projector and return effectively zero results).

%Make the projectors. Note the comments on your screen when this happens.
H           = [Dr.data];
S           = [Dn.data];

%Verify that H and S are respectively full column rank
nr	    = rank(H);
nn      = rank(S);

if(and( nr>(length(rangeInd)-1/2), nn>(length(nullInd)-1/2)))
    
    disp(sprintf('Both H and S have full column rank: Oblique projector condition number OK.'));
    disp(sprintf('Condition numbers are: H ~ %i, while S ~ %i', cond(H), cond(S)));
    
else
    
    disp(sprintf('H and S have linearly dependent columns. Projector degenerate'));
    disp(sprintf('HUGE Condition numbers are: H ~ %i, while S ~ %i', cond(H), cond(S)));
    
end;

%NOW adjust the column space and null space, using the SVD. Take only the
%columns corresponding to the largest singular vectors for each.

%RANGE: ONLY USE SOME OF THE DATA....if the projector is GOOD, it will
%keep the good data NOT included, too. Cut the data off when the
%condition number of the data matrix is still low. Say, order 20.
[Ur,Sr,Vr]  = svd(normColumns(H),0);
Ur          = Ur(:,diag(Sr) > 1/20);

%NULL: ONLY USED SOME OF THE DATA...if the projector is GOOD, it will
%reject the crappy data NOT included, too. Cut the data off when the
%condition number of the data matrix is still low. Say, order 20. The less
%columns you use, the faster, but more approximate, too.
[Un,Sn,Vn]  = svd(normColumns(S),0);
Un          = Un(:,diag(Sn) > 1/20);

[P,Ehs,Esh] = proj2Obliq(Ur,Un);

%NOW RUN DATA THROUGH PROJECTOR!
Hp          = Ehs*H;
Sp          = Ehs*S;

scrsz = get(0,'ScreenSize');
figure('Position',[scrsz(3)/2 scrsz(4) scrsz(3)/2 scrsz(4)]);
h = plotColumns(t,[Hp,Sp]); hold on;
set(gca,'fontsize',18,'fontweight','b','linewidth',2);
set(gcf,'color','w');
set(h(1:length(rangeInd)),'color','k');
set(h(length(rangeInd)+1:length(rangeInd)+length(nullInd)),'color','r');
title('Signals AFTER PROJECTION');

scrsz = get(0,'ScreenSize');
figure('Position',[1 scrsz(4) scrsz(3)/2 scrsz(4)]);
h = imagesc(C); hold on;
colorbar;
set(gca,'fontsize',18,'fontweight','b','linewidth',2);
set(gcf,'color','w');
title('Correlation between Good and Bad data ONLY');

scrsz = get(0,'ScreenSize');
figure('Position',[scrsz(3)/2 scrsz(4) scrsz(3)/2 scrsz(4)]);
[~,nb] = normColumns(S);
[~,np] = normColumns(Sp);
h = zeros(1,2);
h(1) = plot(1:length(nb),nb,'+r'); hold on;
h(2)= plot(1:length(np),np,'ok'); hold on;
set(h(1),'color','r','markersize',12);
set(h(2),'color','k','markersize',12);
set(gca,'fontsize',18,'fontweight','b','linewidth',2);
set(gcf,'color','w');
title('Norm Before Projection, and After Projection of BAD');