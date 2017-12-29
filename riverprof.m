% Elevation extraction.
%
% Chunli Dai, July 2017
% Chunli Dai, December 2017
macdir='/Users/chunlidai/';
%Gage at Fairbanks
filename='filelist';
iref=1; %20111008 the lowest stage. 

fprintf ('\n Step 0: geting the boundary for all files in the region.')
%READ INPUT PARAMETERS; getting the boundaries for all files
% filename='boundaries_reg31.dat';
fid = fopen(filename);
n = linecount(fid);
fid = fopen(filename);
% range=fscanf(fid, '%f', [4, n))';
range=zeros(n,4);
datadir='/data3/ArcticDEM/region_34_alaska_north/tif_results/2m/';
for i=1:n
   range(i,1:4)=fscanf(fid, '%f', [4, 1])';ifile=fgetl(fid);
   [demdir,name,ext] =fileparts([macdir,datadir,strtrim(ifile)]);
   f{i}=[name,ext];
   fdir{i}=[demdir,'/'];%working on two regions 
end
display(['demdir=',demdir])
demdir='';
id=1:n;

%control/stable surfaces for coregistration
riv=load('FAirport2.gmt');
[xap,yap]=polarstereo_fwd(riv(:,2),riv(:,1),[],[],70,-45);
park=load('park.gmt');park2=load('park2.gmt');road=load('road.gmt');
[xpark,ypark]=polarstereo_fwd(park(:,2),park(:,1),[],[],70,-45);
[xpark2,ypark2]=polarstereo_fwd(park2(:,2),park2(:,1),[],[],70,-45);
[xroad,yroad]=polarstereo_fwd(road(:,2),road(:,1),[],[],70,-45);

%for j=1:length(id)%0%length(id)
for j=[1:(length(id)-4),(length(id)-2):length(id)]%0%length(id)
i=id(j);
demdir=fdir{i};
infile= strrep([demdir,'/',f{i}],'meta.txt','dem.tif');
ymd=str2num(f{i}(6:13));

mon=str2num(f{i}(10:11));
if mon <=4 || mon>=11 || j==1; %winter
display(['j=',num2str(j),';ymd=',num2str(ymd)])
else
% continue
end

%coregistration
resr=40;
data=readGeotiff(infile);
res=data.info.map_info.dx;
nsr=resr/res; dsr=res/resr;

%filtering out bad edges, e.g. 20110804
Me=~(imdilate((data.z==-9999),ones(round(30*2*8/res)))); 
data.z(~Me)=-9999;

mtFile = strrep(infile,'dem.tif','matchtag.tif');
mt=readGeotiff(mtFile);

% Reduce resolution of DEM for coregistration.
ranget=[[min(data.x) max(data.x) min(data.y) max(data.y)]/resr];
ranget=[ceil(ranget(1)) floor(ranget(2)) ceil(ranget(3)) floor(ranget(4))]*resr;
%ranget=[-2706760,-2700000,614480,621160];
tx=ranget(1):resr:ranget(2);ty=ranget(4):-resr:ranget(3);
idrxs=find(abs(data.x-ranget(1))<1e-3);idrxe=find(abs(data.x-ranget(2))<1e-3);
idrys=find(abs(data.y-ranget(4))<1e-3);idrye=find(abs(data.y-ranget(3))<1e-3);
idrx=idrxs:nsr:idrxe;
idry=idrys:nsr:idrye;
dd=[idrx(2:end)-idrx(1:end-1),idry(2:end)-idry(1:end-1)];
if isempty(dd)||any(abs(dd-nsr)>1e-3);warning('Resized grid is not constant spacing.');end
tz=data.z(idry,idrx); %0.09s
datar= struct();
datar.x=tx;datar.y=ty;  datar.z=tz;

if iref==i
data0=data;mt0=mt;
data0.z(data0.z==-9999)=NaN;

data0r=datar;%datar=datarsv;
data0r.z(data0r.z== -9999) = NaN;

p=zeros(3,1); %Translation parameters are zeros for the reference DEM.

%Get the mask of control surfaces.
idx=round((xap-datar.x(1))/resr)+1;
idy=round(-(yap-datar.y(1))/resr)+1;
mp2 = poly2mask(idx,idy, length(datar.y),length(datar.x)); % build edge mask
mpairport=mp2;

idx=round((xroad-datar.x(1))/resr)+1;
idy=round(-(yroad-datar.y(1))/resr)+1;
% mproad = poly2mask(idx,idy, length(datar.y),length(datar.x)); % build edge mask
k = boundary(idx,idy,1); % find outer data boundary - this function
mproad = poly2mask(idx(k),idy(k), length(datar.y),length(datar.x));

idx=round((xpark-datar.x(1))/resr)+1;
idy=round(-(ypark-datar.y(1))/resr)+1;
mppark = poly2mask(idx,idy, length(datar.y),length(datar.x)); % build edge mask
idx=round((xpark2-datar.x(1))/resr)+1;
idy=round(-(ypark2-datar.y(1))/resr)+1;
mppark2 = poly2mask(idx,idy, length(datar.y),length(datar.x)); % build edge mask
mppark=mppark|mppark2;
mp1=mpairport|mproad|mppark;

else
datar.z(datar.z== -9999) = NaN; %For accurate interpolation, use NaNs for none data.
% [z2out,p,d0] = coregisterdems(data0r.x,data0r.y,double(data0r.z),datar.x,datar.y,double(datar.z));

%Get the mask of control surfaces.
idx=round((xap-datar.x(1))/resr)+1;
idy=round(-(yap-datar.y(1))/resr)+1;
mp2 = poly2mask(idx,idy, length(datar.y),length(datar.x)); % build edge mask
mpairport=mp2;

idx=round((xroad-datar.x(1))/resr)+1;
idy=round(-(yroad-datar.y(1))/resr)+1;
k = boundary(idx,idy,1); % find outer data boundary - this function
mproad = poly2mask(idx(k),idy(k), length(datar.y),length(datar.x));

idx=round((xpark-datar.x(1))/resr)+1;
idy=round(-(ypark-datar.y(1))/resr)+1;
mppark = poly2mask(idx,idy, length(datar.y),length(datar.x)); % build edge mask
idx=round((xpark2-datar.x(1))/resr)+1;
idy=round(-(ypark2-datar.y(1))/resr)+1;
mppark2 = poly2mask(idx,idy, length(datar.y),length(datar.x)); % build edge mask
mppark=mppark|mppark2;
mp2=mpairport|mproad|mppark;

mp1r = interp2(data0r.x' ,data0r.y,mp1 ,datar.x',datar.y,'*nearest',0);
mpcom=logical(mp1r)&logical(mp2);

% Ratio of control points to total points, larger the better.
RC=sum(mpcom(:))/sum(sum(~isnan(data0r.z)))*100;
display(['RC=',num2str(RC)])
if RC > 0.3 %enough control points
    %Coregistration algorithm in Nuth and Kaab, 2010.
    [z2out,p,d0] = coregisterdems(data0r.x,data0r.y,double(data0r.z),datar.x,datar.y,double(datar.z),mp1,mp2);
else % given control points too less, do not use them.
    [z2out,p,d0] = coregisterdems(data0r.x,data0r.y,double(data0r.z),datar.x,datar.y,double(datar.z));

end

% The case when coregistration is not applied.
z2n = interp2(datar.x' ,datar.y,datar.z ,data0r.x',data0r.y,'*linear');
        
[X,Y]=meshgrid(data0r.x,data0r.y);
[LATa,LONa]=polarstereo_inv(X,Y,[],[],70,-45);

figure;set(gcf,'Color','white');set(gca,'FontSize', 12);set(gcf, 'PaperPosition', [0.25 2.5 6 2.5]);
% imagesc(data0r.x*1e-3,data0r.y*1e-3,z2n-data0r.z);caxis([-5 5]);colorbar
surf(LONa,LATa,z2n-data0r.z); shading interp;
colorbar;colormap jet;view(0,90)
xlabel('{\itx} (km)')
ylabel('y (km)')
caxis([-5 5])
title('2011/10/08-2013/05/26 DEM (m); No coregistration')
title('2012/10/11 - 2011/8/4 DEM (m); No coregistration')

figure;set(gcf,'Color','white');set(gca,'FontSize', 12);set(gcf, 'PaperPosition', [0.25 2.5 6 2.5]);
surf(LONa,LATa,z2out-data0r.z); shading interp;
% imagesc(data0r.x*1e-3,data0r.y*1e-3,z2out-data0r.z);caxis([-5 5]);colorbar
% title('2011/10/08-2013/05/26 DEM (m); After coregistration')
title('2012/10/11 - 2011/8/4 DEM (m); After coregistration')
colorbar;colormap jet;view(0,90)
hl=xlabel('Longitude ($^{\circ}$)');
set(hl, 'Interpreter', 'latex');
hl=ylabel('Latitude ($^{\circ}$)');
set(hl, 'Interpreter', 'latex');
caxis([-8 8])
hold on;plot(riv(:,1),riv(:,2),'r-','linewidth',4)
hold on;plot(park(:,1),park(:,2),'r-','linewidth',4)
hold on;plot(park2(:,1),park2(:,2),'r-','linewidth',4)
hold on;plot(road(:,1),road(:,2),'r.','linewidth',4)
plot3( - (147+50/60+20/3600), 64+47/60+34/3600,1e9, 'ks','Markersize',12,'linewidth',6) %gauge

end
            
M=maskentropy(infile);
M=M.z;
y0s=630998;%y0;
ymd=str2num(f{i}(6:13));
        
% Align the DEM coordinates to reference DEM using the translation parameters from coregistration.
data.x=data.x- p(2);data.y=data.y- p(3);data.z=data.z- p(1);
[X,Y]=meshgrid(data.x,data.y);
[LAT,LON]=polarstereo_inv(X(M),Y(M),[],[],70,-45);

%Output the extracted elevation along shorelines
output=[LAT(:),LON(:),double(data.z(M)),X(M),Y(M)];
ofile=['rivprof',num2str(ymd),'raw.dat'];
save(ofile,'output','-ascii') 

[X,Y]=meshgrid(data.x,data.y);
% data.z=z;
demp=data.z(M);
dempmt=mt.z(M);
epochorg=Y(M);
[idkp ]=outlier(demp,dempmt,epochorg,ymd);
yp=epochorg(idkp);

output=[LAT(idkp),LON(idkp),double(demp(idkp)),(yp(:)-y0s)];
ofile=['rivprof',num2str(ymd),'s.dat']; %remove outlier
save(ofile,'output','-ascii') 
%output=[epoch(:)-y0s,T6(:),T6std(:)];
%ofile=['rivprof',num2str(i),'ms.dat']; %mean and std
%save(ofile,'output','-ascii') 

% %second method
% Retrieve the shoreline location in the lowest-stage DEM coordinates system.
y0s=630998;%y0;
xriv=X(M); yriv=Y(M); %data.x- p(2), in reference DEM coordinates
p=zeros(3,1); % % lowest-stage DEM - p =  reference DEM.
xriv=xriv+ p(2);yriv=yriv+ p(3);%coregistered to the lowest-stage dem. 2012/10/11
% Interpolate the shoreline heights from the lowest-stage DEM.
zriv= interp2(data0.x,data0.y,double(data0s.z),xriv,yriv,'*linear',NaN);
dempmt= interp2(mt0.x,mt0.y,double(mt0.z),xriv,yriv,'*nearest',0);
zriv(isnan(zriv))=-9999;
zriv=zriv-p(1); %back to the reference dem

%Detect the outliers and mean and std
[idkpb]=outlier(zriv,dempmt,epochorg,ymd);
yp=epochorg(idkpb);

%Write the shoreline heights
output=[LAT(idkpb),LON(idkpb),double(zriv(idkpb)),(yp(:)-y0s)];
ofile=['rivprof',num2str(ymd),'b.dat']; %remove outlier
save(ofile,'output','-ascii') 

%plot the orthoimage and shorelines
orFile = strrep(infile,'dem.tif','ortho.tif');
or=readGeotiff(orFile);
resr=40;
res=or.info.map_info.dx;
nsr=resr/res; dsr=res/resr;

or.z = imresize(or.z,dsr);
or.x = imresize(or.x,dsr);
or.y = imresize(or.y,dsr);
[X,Y]=meshgrid(or.x,or.y);
[LATa,LONa]=polarstereo_inv(X- p(2),Y- p(3),[],[],70,-45);

figure;
set(gcf,'Color','white')
set(gca,'FontSize', 12);
set(gcf, 'PaperPosition', [0.25 2.5 6 2.5]);
hold all;
surf(LONa,LATa,or.z); % too slow
shading interp;
colorbar;colormap gray;view(0,90)
caxis([0 10000])
hold on;plot3(LON(idkp),LAT(idkp),1e9*ones(size(LAT(idkp))),'r.','Markersize',8)
 view(0,90)
box on
hl=xlabel('Longitude ($^{\circ}$)');
set(hl, 'Interpreter', 'latex');
hl=ylabel('Latitude ($^{\circ}$)');
set(hl, 'Interpreter', 'latex');
ofile=['rivloc',num2str(ymd)]; %mean and std
hold on;plot( - (147+50/60+20/3600), 64+47/60+34/3600, 'bs','Markersize',12) %gauge
hold on;plot(riv(:,1),riv(:,2),'b-','linewidth',4) %airport
hold on;plot3(LON(idkp(in)),LAT(idkp(in)),1e9*ones(size(LAT(idkp(in)))),'b>','Markersize',8)
hold on;plot3(LON(idkpb(inb)),LAT(idkpb(inb)),1e9*ones(size(LAT(idkpb(inb)))),'co','Markersize',8)
title(num2str(ymd))
saveas(gcf,ofile,'fig')

end %for j


