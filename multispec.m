% get the shoreline from multispectral images using NDWI index

macdir='/Users/chunlidai/surge/';
addpath(genpath([macdir,'/home/chunli/scripts/']));
addpath(genpath([macdir,'/home/chunli/coastline']));

%
%River Tanana; Dai et al., 2018 (GRL), Fig.S1
infile=[macdir,'/home/chunli/scripts/TananaRiver/workNDWI/GE01_20130526211927_10504100029EAE00_13MAY26211927-M1BS-053537405030_01_P001t2.tif'];

data=readGeotiff(infile);

%GEOEYE-1: 1, 2, 3, 4 BLUE GREEN RED NIR1

Green=double(data.z(:,:,2));NIR1=double(data.z(:,:,4)); 
NDWI=double(Green-NIR1)./double(Green+NIR1);%McFeeters 1996

figure;
set(gcf,'Color','white');set(gca,'FontSize', 18);set(gcf, 'PaperPosition', [0 0 6 4]);set(gcf, 'PaperSize', [ 6 4]);
imagesc(data.x,data.y,NDWI);
title('NDWI')
colorbar
caxis([-1 1])
caxis([-0.3 0.1])

M=NDWI>0; %threshold

figure;imagesc(data.x,data.y,M);colorbar

Modj=M;
Modj= bwareaopen(Modj, 1000*500); %remove small clusters
Modfil = bwareaopen(~Modj, 1000*5); %remove small holes
Modfil=~Modfil;

Md1 = imdilate(Modfil, ones(3));
M=logical(Md1-Modfil);
% M=M&~(imdilate(or == 0,ones(round(9*2*8/resr)))); 

figure;imagesc(data.x,data.y,double(Modfil))
hold on;plot(X(M),Y(M),'g.')

[X,Y]=meshgrid(data.x,data.y);

%river shoreline data for Fig. S1(b).
[LAT,LON]=polarstereo_inv(X,Y,[],[],70,-45);
output=[LAT(M),LON(M)];
save coastndwi.dat output -ascii 

Ms=LAT>=64.697&LAT<=64.76&LON>-148.175&LON<-148.03;
output=[LAT(Ms),LON(Ms),double(NDWI(Ms))];
save ndwi.dat output -ascii 

