function plotSeisClusters(year, days, station, chan)
%
% Need to make a function that plots the cluster activity as a funciton of
% time.  If there are N clusters, the plot will give a plot of cluster
% activity over time, with little markers at for y = 1:N. Simple, fast.

%h(1) = subplot(3,3,[1:2,4:5,7:8]); hold on;
%h(2) = subplot(3,3,[3:3:9]); hold on;
%opt.thresh=0.85;



loadFile=sprintf('%s.%s.%s.%i.%i.%i.mat','CLUSTERPICKS',...
    char(station),char(chan),year,days(1),days(2));

load(loadFile);
%Scluster=coralClusterSeis(Scluster,opt);

D       = cell(1,size(Scluster,2));
O       = cell(1,size(Sorphan,2));
L       = zeros(size(Scluster,2),1);

t1      = zeros(size(Scluster,2),1);
t2      = t1;
t3      = zeros(size(Scluster,2),1);
t4      = t3;

for k=1:size(Scluster,2)

    D(k)    = {Scluster(k).recStartTime};
    L(k)    = median([Scluster(k).recNumData]);
    t1(k)   = min( datenum( [D{k}(:,:)]' ) );
    t2(k)   = max( datenum( [D{k}(:,:)]' ) );

end;

for k = 1:size(Sorphan,2)

    O(k)    = {Sorphan(k).recStartTime};
    t3(k)   = min( datenum( [O{k}(:,:)]' ) );
    t4(k)   = max( datenum( [O{k}(:,:)]' ) );

end;

t   = [linspace(min(min(t1),min(t3)), max(max(t2),max(t4)), max(L))]';


thecolor = {'r','g','b'};


for k=1:size(D,2)

    s=plot(datenum([D{k}(:,:)]'),(k-1),'d'); hold on;
    set(s,'markersize',8); set(s,'linewidth',2);
    set(s,'color',char(thecolor(mod(k,3)+1)));

end;

for k = 1:size(O,2)

    s=plot(datenum([O{k}(:,:)]'),size(D,2),'.'); hold on;
    set(s,'markersize',12);
    set(s,'color','k');

end;

grey = 0.5*[1,1,1];

X = [Scluster.data];
X = maxalign(X);

%[w]=signalwidth(X,0.99);

w = [1,size(X,1)];
w = repmat(w,size(X,2),1);

%plot a seismogram without overlap
for k = 1:size(Scluster,2)
    I = w(k,1):w(k,2);
    s = plot(t(I),X(I,k) + (k-1),'color',grey); hold on;
    %set(s,'linewidth',1.25);
end;

set(gca,'fontsize',16,'fontweight','b');
%set(gca,'fontweight','bold');

set(gca,'xtick',linspace(t(1),t(end),12));
datetick('x',2,'keepticks');
ylim([-1,size(D,2)+0.25]);
ylabel('Cluster Waveform Stack','fontsize',16,'fontweight','bold');
%xlabel('Observation Period','fontsize',18,'fontweight','bold');
set(gcf,'color','w');
s=title(sprintf('Seismic Waveform Clusters for Station %s.%s',char(station),char(chan)));
set(s,'fontsize',18,'fontweight','b');

disp('ta-duh');