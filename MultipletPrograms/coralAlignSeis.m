function [data, timeshift, C, beam]=coralAlignSeis(data,dataMat,iflag, index); 
%   coralAlignSeis    align coral seismograms by cross correlation or trace extrema
% USAGE: [data, timeshift, C, beam]=coralAlignSeis(data,dataMat,iflag, index); 
%
%  Coral function to align all of the traces on the maximum value of the traces
%    (iflag==1), or on the minimum value of the traces (iflag==2),
%    or by cross correlation (iflag==3).  For the later case the traces are
%    stacked (without normalizing) and each trace is cross correlated with
%    the stack to determine the time shift.
%    C is the maximum correlation coefficient if flag==3

t=[0:length(dataMat(:,1))-1] * data(1).recSampInt;  % time vector starting at 0
C=min(abs(dataMat));
index=find(C~=0);                        % find non-zero traces
beam = zeros(size(dataMat(:,1)));
if     iflag==3,                         % use cross correlation

  [A,B]=size(dataMat);
  if nargin<5; index=1:B; end
  beam=zeros(A,1);  % beam is normalized sum of non-zero traces in the index list
  for i=1:B
    if C(i)>0 & any(i==index)
      beam=beam+dataMat(:,i)/(max(dataMat(:,i))-min(dataMat(:,i))) ;
    end
  end
    
  C=zeros(1,B);
  timeshift=zeros(1,B);
  for i=1:B
    if max(abs(dataMat(:,i)))==0
      C(i)=0; timeshift(i)=0;
    else
      correl = xcorr(beam,dataMat(:,i),'coeff');
      [C(i),timeshift(i)] = max(correl);
    end
  end

elseif iflag==2,                          % use minumum
  [a,timeshift]=min(dataMat);
else,                                     % use maximum
  [a,timeshift]=max(dataMat);
end

% if there are no data for a trace (all values are zero)
% then C is zero for that trace.  In that case do not shift 
% that trace and force the mean of the rest of the time
% shifts to be zero

timeshift=timeshift.*[data.recSampInt];       % convert from samples to seconds
index=find(C~=0);                        % find non-zero traces
mean_shift=mean(timeshift(index));       % find mean shift of non-zero traces
del_t=timeshift';                    
del_t(index)=del_t(index)-mean_shift;    % remove mean from non-zero traces
index=find(C==0);                        % find indices of zero traces
if length(C)>0,
  del_t(index)=zeros(size(index));       % set their time shifts to zero
end
timeshift=del_t;
%Tdur=header1(2,:);
%window=[del_t';Tdur;(1:length(Tdur))];
%obs1(1,:)=header1(1,:)+del_t';
%obs1(1,:)=obs1(1,:) - mean(obs1(1,:)-header1(9,:)); % differential travel time (s)
%obs1(2,:)=obs1(1,:) - header1(9,:);                 % differential travel time residual (s)
%[data1,header1,label1,obs1]=apply_window(data1,header1,label1,obs1,window);
%disp(label1');
