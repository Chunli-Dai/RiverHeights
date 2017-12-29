%Plot Fig. 3b in Dai et al., GRL, 2018.

close all

ft2m=0.3048; %feet to meter
% macdir='/Users/chunli/surge/';
filename='uv10to16.txt';
fid = fopen(filename);
nheader=30;
n = linecount(fid)-nheader;
fid = fopen(filename);
for i=1:nheader
ifile=fgetl(fid);
end
% range=fscanf(fid, '%f', [4, n))';
data=zeros(n,2);
epoch=zeros(n,1);
for i=1:n
   line1 = textscan(fid, '%5s%15s%20s%s%6s%14s%10s%14s%10s\r\n %*[^\n]');
   if ~isempty(line1{9}{1}) %there is discharge
       data(i,1)= str2double(line1{8}{1}); %stage
       data(i,2)= str2double(line1{6}{1}); %discharge
       epochs=[line1{3}{1} line1{4}{1}];
   else
       data(i,1)= str2double(line1{6}{1}); %stage
       epochs=[line1{3}{1} line1{4}{1}];
   end
   epoch(i)=datenum(epochs,'yyyy-mm-ddHH:MM');
end
id=find(data(:,2)~=0);

%compare with https://waterwatch.usgs.gov/index.php
figure;
plot(data(id,2),data(id,1),'b.-') %stage
xlabel('Discharge (ft^3/s)') 
ylabel('Stage (ft)')

datum=400.67 *ft2m ;
Bias=2.097; % gage height to EGM08 
Bias=2.3; %winter2/gageft.txt 0.3261 rms
figure;
plot(data(id,2)*ft2m^3,data(id,1)*ft2m+datum+Bias,'b.-') %stage
xlabel('Discharge (m^3/s)') 
ylabel('Stage (m)')


% fitting 
nmax=2;% nmax=3 too flat for the lowest stage
t=data(id,1); %stage
% tl=14; tr=30;%water surface height in feet
tl=13; tr=26;%water surface height in feet
tstar=scale4legs(t,[tl,tr]);
n=length(t);
m=nmax;
AM=legs(m,tstar);
condA=cond(AM)
yobs=data(id,2);
est=inv(AM'*AM)*AM'*yobs;
etilde=yobs-AM*est;
sigma02hat=etilde'*etilde/(length(yobs)-1); % sigma=1.5978e+03 ft^3/s = 45.2456 m^3/s
var=inv(AM'*AM);
var=var*sigma02hat;

%get ArcticDEM water height data;
ghf=load('gagefo.txt');
idd=[8,16,19];%[20110804] [20130329] [20130413]
ghf(idd,:)=[];

pL5=[1.8855   0 0];w2d=-11.23; %meter from WGS84 TO egm08 BY Mike Durand
pL5(1)=pL5(1)+w2d;
epocho=datenum(num2str(ghf(:,1)),'yyyymmdd');
idn=1:length(epocho);
mon=str2num(datestr(epocho,'mm'));
ids=find( mon >4 & mon<11);
T6f=ghf(ids,2)+pL5(1);T6fstd=ghf(ids,3);% in meter
xfit=[tl:0.1:tr]';
xfit=(T6f-datum-Bias)/ft2m; %in ft
tstar=scale4legs(xfit,[tl,tr]);
Af=legs(m,tstar);
yfit=Af*est; %discharge  in ft^3/s
%yfitstd
[~,Afd]=legsd(m,tstar);
G=[Afd*2/(tr-tl)*est,Af];
yfitstd=zeros(size(xfit));
for i=1:length(xfit)
    Din=diag([(T6fstd(i)/ft2m)^2,zeros(1,m+1)]);
%     Din=zeros(5,5);
    Din(2:end,2:end)=var;
    Gi=G(i,:);
    yfitstd(i)=sqrt(Gi*Din*Gi');
end

hold on;
plot(yfit*ft2m^3,xfit*ft2m+datum+Bias,'k.-')
plot(yfit*ft2m^3,xfit*ft2m+datum+Bias,'r>')

gh2i = interp1( epoch(id) , data(id,2) , epocho(ids),'linear') ;
dh=(yfit-gh2i)*ft2m^3;
RMSE=std(dh) %194 m^3/s 206
RMSE=sqrt(mean(dh.^2)) %220 m^3/s; formal std is about 100.
%For updated stage (gagefsv4.txt from fitting curve) RMSE 234 (std 232); formal std
%around 220.
mean(dh) % 71.5

%Plot Fig. 3b.
figure;
set(gcf,'Color','white')
set(gca,'FontSize', 12);
set(gcf, 'PaperPosition', [0 0 8 2.5]); 
set(gcf, 'PaperSize', [8 2.5]); 
set(gcf, 'PaperPosition', [0 0 8 2.8]);
set(gcf, 'PaperSize', [8 2.8]);
hold all;box on
plot(epoch(id),data(id,2)*ft2m^3,'b.','markersize',9);
errorbar(epocho(ids),yfit*ft2m^3,yfitstd*ft2m^3,'r.','linewidth',1,'markersize',12);
legend('USGS gage','ArcticDEM')
% axis([datenum('2010-01-01','yyyy-mm-dd') datenum('2016-01-01','yyyy-mm-dd') 0 2500])
axis([datenum('2010-01-01','yyyy-mm-dd') datenum('2016-01-01','yyyy-mm-dd') 0 2700])
datetick('x','yyyy','keeplimits')
ylabel('Discharge (m^3/s)') 
ofile='disch';
saveas(gcf,ofile,'fig')
print('-dpng','-r400',ofile)
print('-dpdf','-r400',ofile)



