%Plot the time series of river heights at the USGS gage in Fairbanks, Alaska.
%Plot Fig.3a in Dai et al., GRL, 2018.
% Chunli Dai, July 2017

%Gage station: https://nwis.waterdata.usgs.gov/nwis/inventory/?site_no=15485500&agency_cd=USGS
ghf=load('gagefo.txt'); %gage height data
ghfw=load('gageft.txt'); % Winter data
pL5=[1.8855   0 0]; %Vertical offset between ArcticDEM DEM (2011/10/8) with LiDAR.
%Control surface:  FAirport2.gmt FAirport3.gmt Pondsh.gmt road.gmt
gh2=importdata('usgsgage.txt');
ft2m=0.3048;%feet to meter

% Vertical datum: WGS84
w2d=-11.23; %meter from WGS84 TO egm08 BY Mike Durand
%-11.23 from web https://geographiclib.sourceforge.io/cgi-bin/GeoidEval?input=64%B047%2734%22+-147%B050%2720%22&option=Submit
pL5(1)=pL5(1)+w2d;

epoch=datenum(num2str(ghf(:,1)),'yyyymmdd');
epochw=datenum(num2str(ghfw(:,1)),'yyyymmdd');

%replace ghf winter times with ghfw;
idnn=1:length(epoch);
idw=idnn(ismember(epoch,epochw));
ids=idnn(~ismember(epoch,epochw));
ghf(idw,:)=ghfw;
T6f=ghf(:,2)+pL5(1);T6fstd=ghf(:,3);
T6bf=ghf(:,4)+pL5(1);T6bfstd=ghf(:,5);
hs=[epoch(ids),T6f(ids),T6fstd(ids)];

% Datum of gage: 404.93 feet above   NAVD88. %https://nwis.waterdata.usgs.gov/nwis/inventory/?site_no=15485500&agency_cd=USGS
% using Vdatum, 404.93' above NAVD88 = 400.67' above EGM08 by Mike Durand
datum=400.67 *ft2m ; % EGM08 

%Load USGS gage time series
for j=1:length(gh2.textdata)
datet{j}=[gh2.textdata{j,1} gh2.textdata{j,2}];
end
epoch2=datenum(datet,'yyyy-mm-ddHH:MM');%2010-06-01 03:15
gh2=gh2.data*ft2m+datum; % to meter

id=find((epoch2(2:end)-epoch2(1:end-1))==0);
epoch2(id)=[];gh2(id)=[]; %Delete repeat data
% Interpolate USGS gage measurements to the same epoch of the ArcticDEM measurements
gh2i = interp1( epoch2 , gh2 , epoch,'linear') ;
dh=T6f-gh2i;
% if there is no gage, the interpolation is not accurate.
%exclude the epochs that has no gage data, for more accurate bias calculation
id=[2,3,4,5,6,8, 10,11,12, 16,19]; % 8 20110804; 19 20130413
dh(id)=[];
RMSE=std(dh)
Bias=mean(dh) %
Bias=2.3; %winter2/gageft.txt 0.3261 rms
gh2=gh2+Bias; %Adjust the USGS gage measurements to the vertical datum of ArcticDEM;

RMSEdt=std(dh(2:end)-dh(1:(end-1)))
Biasdt=mean(dh(2:end)-dh(1:(end-1)))

id=find(epoch==datenum('2013-03-29'));epoch(id)=[];T6f(id)=[];T6fstd(id)=[];
id=find(epochw==datenum('2013-03-29'));epochw(id)=[];ghfw(id,:)=[];
stdm=1.;

%Plot Fig 3a.
figure  
set(gcf,'Color','white')
set(gca,'FontSize', 12);
set(gcf, 'PaperPosition', [0 0 8 3]);
set(gcf, 'PaperSize', [8 3]);
hold all;
plot(epoch2,gh2,'b.')
hold on
errorbar(epoch,T6f,T6fstd,'r.','linewidth',1,'markersize',12)
errorbar(epochw,ghfw(:,2)+pL5(1),ghfw(:,3),'k.','linewidth',1,'markersize',12)
plot(epochw,ghfw(:,2)+pL5(1),'ko','markersize',4)
legend('USGS gage','ArcticDEM','ArcticDEM Winter')
box on
ylabel('Elevation (m)')
datetick('x','yyyy')
ofile=['gagefawinp'];
axis([datenum('2010-01-01','yyyy-mm-dd') datenum('2016-01-01','yyyy-mm-dd') 128 133.5])
saveas(gcf,ofile,'fig')
print('-dpng','-r400',ofile)
print('-dpdf','-r400',ofile)

