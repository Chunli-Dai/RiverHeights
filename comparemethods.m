%Plotting the river height time series by two different methods.

clear
close all
ghf=load('gagefo.txt');
ghfw=load('winter2/gageft.txt'); % correction of the odd lows in 2013
epoch=datenum(num2str(ghf(:,1)),'yyyymmdd');
epochw=datenum(num2str(ghfw(:,1)),'yyyymmdd');
%replace ghf winter times with ghfw;
idnn=1:length(epoch);
idw=idnn(ismember(epoch,epochw));
ids=idnn(~ismember(epoch,epochw));
ghf(idw,:)=ghfw;

id=[8,16,19];%[20110804] [20130329] [20130413]
id=[2,3,4,5,6,8, 10,11,12, 16,19]; %No gage data % 8 20110804 19 20130413 
ghf(id,:)=[];
pL5=[1.8855   0 0]; %  2011/10/8 with LiDAR FAirport2.gmt FAirport3.gmt Pondsh.gmt road.gmt
gh2=importdata('usgsgage.txt');
ft2m=0.3048;%feet to meter
w2d=-11.23; %meter from WGS84 TO egm08 BY Mike Durand
pL5(1)=pL5(1)+w2d;
epoch=datenum(num2str(ghf(:,1)),'yyyymmdd');
epocht=epoch;
idn=1:length(epocht);
mon=str2num(datestr(epocht,'mm'));
idw=find( mon <=4 | mon>=11);
ids=find( mon >4 & mon<11);

datum=400.67 *ft2m ; % EGM08
for j=1:length(gh2.textdata)
datet{j}=[gh2.textdata{j,1} gh2.textdata{j,2}];
end
epoch2=datenum(datet,'yyyy-mm-ddHH:MM');%2010-06-01 03:15
gh2=gh2.data*ft2m+datum; % to meter
id=find((epoch2(2:end)-epoch2(1:end-1))==0);
epoch2(id)=[];gh2(id)=[];
gh2i = interp1( epoch2 , gh2 , epoch,'linear') ;
Bias=2.097; %winter/gageft.txt add 20130413
Bias=2.3; %winter2/gageft.txt 0.3261 rms
gh2=gh2+Bias; %align with T6f ArcticDEM;
T6f=ghf(:,2)+pL5(1);
dh=T6f-(gh2i+Bias);

if 0
    %RMS of dh
RMSEdt=std(dh(idw)) %sqrt(mean((dh(idw)-Bias).^2))
Biasdt=mean(dh(idw))
RMSEdt=std(dh(ids)) %sqrt(mean((dh(ids)-Bias).^2))
Biasdt=mean(dh(ids))
RMSEdt=std(dh) %sqrt(mean((dh-Bias).^2)) %std 0.3261m rms 0.3189 m
Biasdt=mean(dh)
T6f=zeros(size(gh2i));
T6f(idw)=[ghf(idw,2)]+pL5(1);
T6f(ids)=[ghf(ids,4)]+pL5(1);
dh=T6f-(gh2i+Bias);

end
epoch(idw)=[];
ghf(idw,:)=[];gh2i(idw,:)=[];
T6f=ghf(:,2)+pL5(1);T6fstd=ghf(:,3);
T6bf=ghf(:,4)+pL5(1);T6bfstd=ghf(:,5);

dh=T6f-(gh2i+Bias);
mean(dh) % 0.0717
std(dh) %0.355
sqrt(mean(dh.^2)) %0.3495 %0.30 m
dh=T6bf-(gh2i+Bias);
mean(dh) %-0.4
std(dh) %0.28
sqrt(mean(dh.^2)) %0.48 %0.33 m
T6bf=T6bf+0.3587;
dh=T6bf-T6f;
std(dh) % 0.26 sqrt(mean(dh.^2)) %0.21 m
mean(dh) %-0.1141

%
figure
set(gcf,'Color','white')
set(gca,'FontSize', 12);
set(gcf, 'PaperPosition', [0 0 8 3]);
set(gcf, 'PaperSize', [8 3]);
hold all;
plot(epoch2,gh2,'b.')
errorbar(epoch,T6f,T6fstd,'r.','linewidth',1,'markersize',12)
errorbar(epoch,T6bf,T6bfstd,'g.','linewidth',1,'markersize',12)
legend('USGS Gage','Direct Method','Imagery-Altimetry Method')
box on
ylabel('Elevation (m)')
% datetick('x','mm/yyyy')
datetick('x','yyyy')
% axis([datenum('2010-01-01','yyyy-mm-dd') datenum('2016-01-01','yyyy-mm-dd') 128 133])
axis([datenum('2010-01-01','yyyy-mm-dd') datenum('2016-01-01','yyyy-mm-dd') 128 133.5])
ofile=['gagefbp'];
saveas(gcf,ofile,'fig')
print('-dpng','-r400',ofile)
print('-dpdf','-r400',ofile)
