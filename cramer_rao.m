close all hidden;
clear all;

Matrixbinning=0.02;
Msize=2*1;

for n=1:Msize/Matrixbinning+1
for i=1:Msize/Matrixbinning+1

isc=(i-Msize/2/Matrixbinning-1)*Matrixbinning;
nsc=(n-Msize/2/Matrixbinning-1)*Matrixbinning;

% tsk transit time average
binning=0.5e-12;
LY=4577; % 
s=66e-12; % sigma trans
t_d=30.3e-9*10^nsc;
t_r=70e-12*10^isc;

R=0;
t_d2=t_d;
t_r2=t_r;


if s<4*t_r
    FS=@(t,tsk) (1-R)/(t_d-t_r)*(1/2*exp((s^2-2*t*t_d+2*t_d*tsk)/(2*t_d^2)).*(1-erf((s^2+t_d*(tsk-t))/(sqrt(2)*s*t_d)))- ...
        1/2*exp((s^2-2*t*t_r+2*t_r*tsk)/(2*t_r^2)).*(1-erf((s^2+t_r*(tsk-t))/(sqrt(2)*s*t_r)))) + ...
        R/(t_d2-t_r2)*(1/2*exp((s^2-2*t*t_d2+2*t_d2*tsk)/(2*t_d2^2)).*(1-erf((s^2+t_d2*(tsk-t))/(sqrt(2)*s*t_d2)))- ...  
        1/2*exp((s^2-2*t*t_r2+2*t_r2*tsk)/(2*t_r2^2)).*(1-erf((s^2+t_r2*(tsk-t))/(sqrt(2)*s*t_r2))));  
else
    FS=@(t,tsk) (1-R)/(t_d-t_r)*(1/2*exp((s^2-2*t*t_d+2*t_d*tsk)/(2*t_d^2)).*(1-erf((s^2+t_d*(tsk-t))/(sqrt(2)*s*t_d))));
end


skew=0e-9;    
x=-10*s:binning:100e-9;
 figure;
 hold on;
 plot(x,FS(x,skew)); 
 plot(x,gradient(FS(x,skew)),'red');

FSnum=FS(x,skew);
dFSnum=-gradient(FSnum)/binning;

% figure;
% hold on;
% plot(x,FS(x,skew)); 
% plot(x,dFSnum,'red');

gpd=find(FSnum(1:length(FSnum))~=0);
I=sum(1./FSnum(gpd).*dFSnum(gpd).^2*binning); %Fisher information
CTRscan(i,n)=sqrt(1/I*1/LY)*3.33*1e12; %CTR
trscan(i)=10^(isc-Matrixbinning/2);
trscancon(i)=10^(isc);
sscan(n)=10^(nsc-Matrixbinning/2);
sscancon(n)=10^(nsc);

end
CTRscan
end

i=length(trscan)+1;
isc=(i-Msize/2/Matrixbinning-1)*Matrixbinning;
trscan(i)=10^(isc-Matrixbinning/2);
n=length(sscan)+1;
nsc=(n-Msize/2/Matrixbinning-1)*Matrixbinning;
sscan(n)=10^(nsc-Matrixbinning/2);

     ErgName='CTRMatrix.mat';
     FErg=fopen(ErgName,'w');
     for iiE=1:1:length(trscan)-1
         for nnE=1:1:length(sscan)-1
             fprintf(FErg,'%g \t %g \t %g \n',trscan(iiE),sscan(nnE),CTRscan(iiE,nnE));
             fprintf(FErg,'%g \t %g \t %g \n',trscan(iiE),sscan(nnE+1),CTRscan(iiE,nnE));
         end
         fprintf(FErg,'\n');
         for nnE=1:1:length(sscan)-1
             fprintf(FErg,'%g \t %g \t %g \n',trscan(iiE+1),sscan(nnE),CTRscan(iiE,nnE));
             fprintf(FErg,'%g \t %g \t %g \n',trscan(iiE+1),sscan(nnE+1),CTRscan(iiE,nnE));
         end
         fprintf(FErg,'\n');
     end
     
     size(CTRscan)
     size(sscan)
     size(trscan)
     
dlmwrite('CTRscanMatrix.mat', CTRscan, 'delimiter', '\t');
     
clear v;
z=1;
for i=-2:0.1:5
v(z)=floor(10^i);
z=z+1;
end
figure;
C=contour(log10(trscancon),log10(sscancon),CTRscan',v,'ShowText','on','LabelSpacing',150,'linewidth',1.5);
% set(gca,'XScale','log')
% set(gca,'YScale','log')

% ylabel('SPTR [log(multiplication factor)]');
% xlabel('scintillation rise time [log(multiplication factor)]');
ylabel('log(fall time normalized to 30.3ns)');
xlabel('log(rise time normalized to 70ps)');

%axis tight
set(gca,'FontSize',16)
h = get(gca,'ylabel');
set(h,'FontSize',16)
h = get(gca,'xlabel');
set(h,'FontSize',16)
h = get(gca,'title');
set(h,'FontSize',16)
print('-djpeg','-r300','contour.jpg');
