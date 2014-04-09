function PlotSTACKSauto(stack)

%will automatically create plot of each of the stacks in input cell array
%of stacks. Will label stations and name the event based on its first
%appearance month day and hour ex.J01-1
freqlim=[0.5 10];

for i=1:size(stack,1)
C=stack{i,1};
%get name of stack
earliest=datevec(min(datenum(C(1,1).multiplets')));
figure(i)


for j=1:size(C,1)
 %plot seismogram
 subplot(size(C,1),3,j*3-2:j*3-1)
 tvec=0:C(j,1).recSampInt:(length(C(j,1).data)-1)*C(j,1).recSampInt;
 plot(tvec,C(j,1).data)
 ylabel(C(j,1).staCode);
 A=get(gca,'Position');
 set(gca,'Position',A+[0 -0.01 0 0.01])
 
 if j==size(C,1)
 xlabel('Seconds')    
 elseif j==1     
 title(['multiplet ',datestr(earliest,'mdd-HH'),' -> ',num2str(length(C(1,1).multiplets)),' total events'])
 else
 set(gca,'XTickLabel',[])
 end
 
 %plot spectrum
 subplot(size(C,1),3,j*3)
 [f,a]=amplitudespectrum(tvec,C(j,1).data);
 plot(f,a)
 xlim(freqlim)
 A=get(gca,'Position');
 set(gca,'Position',A+[0 -0.01 0 0.01])
 
 if j==size(C,1)
 xlabel('Freq (Hz)')
 else 
 set(gca,'XTickLabel',[])    
 end
    
end
    
 
   
end

