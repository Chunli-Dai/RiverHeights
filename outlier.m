function [idkp]=outlier(demp,dempmt,epochorg,ymd)
y0s=630998;%y0;
t=epochorg;
ni=length(demp);
idnn=find(abs(demp--9999)>100);idout=[];idoutpre=[];
idkp= idnn(~ismember(idnn,idout));
mp=2;
minnpt=mp;
P=1;
yobs=demp(idkp);
a=2/(max(epochorg)-min(epochorg));
b=-1-a*min(epochorg);
AM=[];
AM=zeros(length(idkp),mp);  %
AMa=ones(ni,1);
AM(:,1)=AMa(idkp,:);
AM(:,2)=epochorg(idkp)*a+b;
est=inv(AM'*P*AM)*AM'*P*yobs;
fit=AM*est;
etilde=yobs-AM*est;
sigma02hat=etilde'*P*etilde/(length(yobs)-1);
meani=mean(fit);  %est(1); bad value for less keeped points and a bad slant angle.
stdi=sqrt(sigma02hat);
multi=3;%6;%3;
id=find(abs(etilde)>=multi*stdi);idout=[idoutpre;idkp(id)];idsign=id;
tkp=t(idkp);[ts,idsort]=sort(tkp);

figure % (1)
set(gcf,'Color','white')
set(gca,'FontSize', 18);
set(gcf, 'PaperPosition', [0.25 2.5 4 3]); %Default [0.25 2.5 8.0 6.0];Prefer 6 by 4.5 or 4 by 3
hold all
plot(t(idkp(idsort)),demp(idkp(idsort)),'b>-','MarkerSize',12,'linewidth',4)
%                     plot(t(idkp),demp(idkp),'b>-','MarkerSize',12,'linewidth',4)
plot(t(dempmt==0),demp(dempmt==0),'co','MarkerSize',18,'linewidth',4)
plot(t(idout),demp(idout),'ks','MarkerSize',12,'linewidth',4)
plot(t(idkp(idsort)),fit(idsort),'g*-','MarkerSize',12,'linewidth',4)
plot(t(idkp(idsort)),fit(idsort)-multi*stdi,'r-','linewidth',4)
plot(t(idkp(idsort)),fit(idsort)+multi*stdi,'r-','linewidth',4)
datetick('x','mm/yy')
box on
ylabel('DEM time series')

idkp= idkp(~ismember(idkp,idout));
idout2=find(dempmt==0);idkp= idkp(~ismember(idkp,idout2));
tkp=t(idkp);[ts,idsort]=sort(tkp);
figure;plot(t(idkp(idsort)),demp(idkp(idsort)),'b>-','MarkerSize',12,'linewidth',4)

%get the elevation time series without outlier
zp=demp(idkp);yp=t(idkp);
y0=min(yp);y1=max(yp); 
% ty=linspace(y0,y1,200);
ty=y0:100:y1;
tz=zeros(size(ty(1:end-1)));tzs=zeros(size(ty(1:end-1)));
for j=1:(length(ty)-1);
    tzt=zp(yp>=ty(j)&yp<ty(j+1));
    tz(j)=mean(tzt,'omitnan');tzs(j)=std(tzt,'omitnan');
end
figure;hold on ; plot(ty(1:end-1),tz,'r.-')

epoch=ty(1:end-1);T6=tz;T6std=tzs;
figure  
set(gcf,'Color','white')
set(gca,'FontSize', 12);
set(gcf, 'PaperPosition', [0.25 2.5 6 2.5]); 
hold all;
plot((yp-y0s)*1e-3,zp,'g.','MarkerSize',1,'linewidth',1)
plot((epoch-y0s)*1e-3,T6(:),'b.-','linewidth',1,'markersize',1)
hold on
shadedErrorBar((epoch-y0s)*1e-3,T6(:),T6std(:),'b',1)
box on
ylabel('Elevation (m)')
xlabel('x (km)')
ofile=['rivprof',num2str(ymd)];
% axis([0 18 115 140])
% axis([0 18 120 145])
% [dir,str,ext] =fileparts(infile);
title([num2str(ymd)])
saveas(gcf,ofile,'fig')

end
