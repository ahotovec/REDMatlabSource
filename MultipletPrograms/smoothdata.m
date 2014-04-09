function [dataout]=smoothdata(datain,varargin)
% Smooths data with a boxcar and plots it.
%------------------------------------------------------------------------------------------%
% INPUT
%datain: A time series vector.  If a row vector is input, a column vector
%        will be output though.
%Npts:   (Optional) As an option, you may input the number of points you would like to smooth over.
%        If nothing is input other than datain, the default of 1/10th the data length is used.
%plotopt: (Optional) A string, 'plot'.  If some other string is input, you get no plot.
%
%OUTPUT
%dataout: The smoothed data.
%NOTE: If plotopt is input, a plot pops up too.
%
%Example: 
%t=linspace(0,1,100); t=t';
%x=cos(t*pi)+rand(100,1);
%y=smoothdata(x,10,'plot');
% 
% NOTE: Default values are used if no '10' is input. 'plot' is optional.
%------------------------------------------------------------------------------------------%
[Mx,Nx]=size(datain);

if(Nx~=1)
    datain=datain';
    [Mx,Nx]=size(datain);
end;

plotopt='none';
Npts=floor(1/10*Mx);

% Check plotting, Npts options
if(nargin==2&&isnumeric(varargin{1}))
    Npts=varargin{1};
end;
if(nargin==2&&strcmp('plot',varargin{1}))
    plotopt=varargin{1};
end;
if(nargin==3&&strcmp('plot',varargin{1}))
    plotopt=varargin{1};
end;
if(nargin==3&&strcmp('plot',varargin{2}))
    plotopt=varargin{2};
end;
if(nargin==3&&isnumeric(varargin{1}))
    Npts=varargin{1};
end;
if(nargin==3&&isnumeric(varargin{2}))
    Npts=varargin{2};
end;

b=ones(Npts,1);
datatemp=zeros(3*Mx,1);
datatemp(1:Mx)=flipud(datain);
datatemp(Mx+1:2*Mx)=datain;
datatemp(2*Mx+1:end)=flipud(datain);
datatemp=conv(datatemp,b);
n=round(Npts/2);
dataout=1/(Npts)*datatemp(Mx+n:end-Mx-n);

if(strcmp(plotopt,'plot'));
    plot(datain,'k','linewidth',1.7); hold on;
    plot(dataout,'r','linewidth',1.7);
    legend('Input data','Smoothed Data');
    s=title(sprintf('Data Smoothed via Convolution with a %i pt. Boxcar',Npts));
    set(s,'fontsize',14);
    set(gca,'Color', [0.85,0.9,0.95]);
    set(gcf,'Color', [1,1,1]);
end;