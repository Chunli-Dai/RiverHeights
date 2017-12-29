% Filtering and fitting of river height profiles.
% Plot Fig.2b in Dai et al., GRL, 2018.
% Modified based on code by Mike Durand, https://github.com/mikedurand/SmoothRiverElevations
% Modified by Chunli Dai, 2017

clear all
%  macdir=[''];
addpath([macdir,'//home/chunli/scripts'])
addpath([macdir,'//home/chunli/scripts/riverfilter/'])

%%
c=shaperead('tan_cl_Close7.shp');
c=c(1);

%get date list
filename='gagefo.txt';
fid = fopen(filename);
n = linecount(fid);
fid = fopen(filename);
ymdg=zeros(n);
for i=1:n
   ymdg(i)=fscanf(fid, '%f', [1, 1])';ifile=fgetl(fid);
end

xmin=[];xmax=[];
xminb=[];xmaxb=[];
xming=[];xmaxg=[];
slope=[];

fid2 = fopen('gageft.txt', 'a');
fid3 = fopen('gageslope.txt', 'a');
for i=1:n
    ymd=ymdg(i);
    
    ymds=num2str(ymd);
    mon=str2num(ymds(5:6));
    if 0 & ( mon <=4 || mon>=11)  %winter
    continue
    else
    display(['Run i=',num2str(i),';ymd=',num2str(ymd)])
    end
    
    display(['i=',num2str(i),';ymd ',num2str(ymd)])
    rivprof1=['rivprof',num2str(ymd),'s'];
    rivprof2=['rivprof',num2str(ymd),'b'];
    
    load([rivprof1,'.dat']);
    load([rivprof2,'.dat']);

    % example of Fig. 2b.
    rivprof1='rivp20110804b';rivprof2='rivp20121011a';
    load('rivp20121011a.mat');
    load('rivp20110804b.mat');

%assign the loaded profiles to the cell s.
s{1}=eval(rivprof1);
s{2}=eval(rivprof2);
 
%parameters
reAttach=true; %best to keep set to true
p.distmax=500; %max distance from centerline to include pixels
p.ncut=20;     %min # of pixels required for each node
p.StdMax=2.5;  %max pixel standard deviation to still include
p.pct=25;      %percentile of pixels chosen to represent that node
p.dx=.1;       %node spacing in km
p.x=1:p.dx:120;
p.N=length(p.x);
p.HrefCut=1;   %cutoff for maximum elevation difference allowed from "reference" (i.e. first pass) 

%% first pass
Href{1}=nan(1,p.N); Href{2}=nan(1,p.N);
[Est]=ProcessData(c,s,reAttach,p,Href,'Fairbanks','LP');

%% second pass
reAttach=true;
Href{1}=Est{1}.Hc; Href{2}=Est{2}.Hc;
[Est,Data]=ProcessData(c,s,reAttach,p,Href,'Fairbanks','SLM'); % 5 km might be too short for this

%%
xobs=93; % from CL #7 from Rui
% using Vdatum, 404.93' above NAVD88 = 400.67' above EGM08
Hobs=([21.9 18.15]+400.67).*.3048; %above EGM08
dHxform=-11.23; % using Vdatum, 0 m WGS84 = -11.23 m EGM08

if 0
%% Plotting 
xmin=min([Data{1}.lon; Data{2}.lon;]);
xmax=max([Data{1}.lon; Data{2}.lon;]);
ymin=min([Data{1}.lat; Data{2}.lat;]);
ymax=max([Data{1}.lat; Data{2}.lat;]);

%% Plot the shoreline location.
figure(1)
subplot(311)
mapshow(c); hold on;
plot(Data{1}.lon,Data{1}.lat,'.',Data{1}.lon(Data{1}.iuse),Data{1}.lat(Data{1}.iuse),'.'); hold off;
axis([xmin xmax ymin ymax]); hold off;
legend('All','Use')
title('rivp20110804')
subplot(312)
mapshow(c); hold on;
plot(Data{1}.lon,Data{1}.lat,'.',Data{2}.lon(Data{2}.iuse),Data{2}.lat(Data{2}.iuse),'.'); hold off;
axis([xmin xmax ymin ymax]); hold off;
title('rivp20121011')
subplot(313)
mapshow(c); hold on;
plot(Data{1}.lon(Data{1}.iuse),Data{1}.lat(Data{1}.iuse),'.',Data{2}.lon(Data{2}.iuse),Data{2}.lat(Data{2}.iuse),'.'); hold off;
hold on;plot( - (147+50/60+20/3600), 64+47/60+34/3600, 'bs','Markersize',12) %gauge
axis([xmin xmax ymin ymax]); hold off

%% Plot the river elevation profiles
figure(2)
subplot(411)
h311=plot(Data{1}.FDh(Data{1}.iuse),Data{1}.h(Data{1}.iuse),'.',p.x,Est{1}.Hhat,'o');
set(h311(1),'Color',[.7 .7 .7])
title('August 2011')
subplot(412)
h312=plot(Data{2}.FDh(Data{2}.iuse),Data{2}.h(Data{2}.iuse),'.',p.x,Est{2}.Hhat,'o');
set(h312(1),'Color',[.7 .7 .7])
title('October 2012')
subplot(413)
plot(p.x,Est{1}.Hhat,'o',p.x,Est{2}.Hhat,'o','LineWidth',2)
set(gca,'FontSize',14)
grid on; legend('August 2011','October 2012','Location','Best')
subplot(414)
plot(p.x,Est{1}.HhatStd,'o',p.x,Est{2}.HhatStd,'o','LineWidth',2)
grid on
set(gca,'FontSize',14)
end
%% Plot the fitting curve. Fig. 2b
figure(3);
set(gcf,'Color','white')
hold all
han1=plot(p.x(Est{1}.Use),Est{1}.Hhat(Est{1}.Use)+dHxform,'o',p.x(Est{2}.Use),Est{2}.Hhat(Est{2}.Use)+dHxform,'o'); hold on;
han2=plot(p.x(Est{2}.iSolUse),Est{2}.Hc(Est{2}.iSolUse)+dHxform,'g',p.x(Est{1}.iSolUse),Est{1}.Hc(Est{1}.iSolUse)+dHxform,'r','LineWidth',2); 
han3=plot(xobs,flip(Hobs)+1.23,'rs','LineWidth',2,'MarkerSize',8);% hold off;
set(gca,'FontSize',16)
% set(han1(1:2),'Color',[.7 .7 .7])
set(han2(1),'Color',[0.93 0.69 0.13])
set(han2(2),'Color',[0.49 0.18 0.56])
set(han3(1),'Color',get(han2(1),'Color'))
set(han3(2),'Color',get(han2(2),'Color'))
set(han1(2),'Color',get(han2(1),'Color'))
set(han1(1),'Color',get(han2(2),'Color'))
legend([han2; han3;],'ArcticDEM 11 Oct 2012','ArcticDEM 4 Aug 2011','USGS 11 Oct 2012','USGS 4 Aug 2011','Location','Best')
xlabel('Main channel flow distance, km')
ylabel('Water surface elevation (EGM08), m')
set(gca,'XLim',[85 94])
box on

xming(i)=min(p.x(Est{1}.iSolUse));xmaxg(i)=max(p.x(Est{1}.iSolUse));

xmin=max([xmin min(p.x(Est{1}.iSolUse))]);
xmax=min([xmax max(p.x(Est{1}.iSolUse))]);

xminb=max([xminb min(p.x(Est{2}.iSolUse))]);
xmaxb=min([xmaxb max(p.x(Est{2}.iSolUse))]);

 %use fitting line for average
id=find(abs(Est{1}.xp-xobs)<=500e-3 );nn=length(id);
idb=find(abs(Est{2}.xp-xobs)<=500e-3);nnb=length(idb);
hest=Est{1}.yp(id);heststd=Est{1}.stdres*ones(size(hest));
hestb=Est{2}.yp(idb);heststdb=Est{2}.stdres*ones(size(hestb));
hgage=mean(hest);hstd=Est{1}.stdres; %1./nn*sqrt(sum(heststd.^2));
hgageb=mean(hestb);hstdb=Est{2}.stdres;  %1./nnb*sqrt(sum(heststdb.^2));

% 92.9000  100.1000; 7.2km
xl=92.9;xr=100.1;
id=find(Est{1}.xp>=xl&Est{1}.xp<=xr);
idb=find(Est{2}.xp>=xl&Est{2}.xp<=xr);
slope(1:2,i)=[mean(Est{1}.Slope(id)) mean(Est{2}.Slope(idb))];

fprintf(fid2,'%d %12.6f %12.6f  %12.6f %12.6f\n',ymd,hgage,hstd,hgageb,hstdb);
fprintf(fid3,'%d %12.6f %12.6f \n',ymd,slope(1:2,i));

end %for i=1:n

