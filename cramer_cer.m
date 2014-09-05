close all hidden;
clear all;

Matrixbinning=0.02;
Msize=2*1;

for n=1:Msize/Matrixbinning+1
for i=1:Msize/Matrixbinning+1

isc=(i-Msize/2/Matrixbinning-1)*Matrixbinning;
nsc=(n-Msize/2/Matrixbinning-1)*Matrixbinning;

% tsk transit time average
binning = 0.5e-12;
LY = 4577; % Light Yield 
CY = 30; % Cerenkov yield
a = LY / CY * (1 / (1 + LY/CY)); % normalization light yield
b = 1 / (1 + LY / CY); % normalization cerenkov yield
s = 66e-12; % sigma trans
t_d = 30.3e-9*10^nsc; % decay time
t_r = 70e-12*10^isc; % rise time
l = 10.0e-12; % cerenkov sigma
theta = 200e-12; % gamma time, parameter CR LB
tsk = 2e-9; % media trans

cer_min = theta - 2 * l; % cer min

%norm_s = 1.0 / (t_d-t_r); % normalization shao
norm_s = a * 1.0 / (t_d-t_r); % normalization shao

%norm_irf = 2 / sqrt(pi) / sqrt(2) / s / erfc(- tsk / s / sqrt(2)); % normalization gaussian
norm_irf = 2. / s / sqrt(2) / sqrt(pi) / erfc(- tsk / s / sqrt(2));

%norm_c =  1. / l / sqrt(2. * pi) * (1. + l / 2. * sqrt(2. * pi) * erfc(1. / sqrt(2) / l * (l + 4.*tsk_c))); % cerenkov normalization 
norm_c = b * 2 / (l * sqrt(2 * pi) * erfc(-sqrt(2))); 

FS=@(t) sqrt(pi / 2) * s * norm_irf * norm_s * ...
        (exp((s*s-2*t*t_d+2*t_d*tsk+2*t_d*theta)/(2*t_d*t_d)) .* ...
          (erf((t_d*(tsk-t+theta)+s*s)/(sqrt(2)*s*t_d))-erf((t_d*tsk+s*s)/(sqrt(2)*s*t_d))) + ...
         exp((s*s-2*t*t_r+2*t_r*tsk+2*t_r*theta)/(2*t_r*t_r)) .* ...
          (erf((t_r*(tsk-t+theta)+s*s)/(sqrt(2)*s*t_r))-erf((t_r*tsk+s*s)/(sqrt(2)*s*t_r)))) + ...
        1.0 / (sqrt(l*l+s*s)) * sqrt(pi / 2) * l * s * norm_c * norm_irf * ...
        (exp(-(tsk+theta-t)*(tsk+theta-t)/2/(l*l+s*s)) .* ...
          erf((-tsk*l*l+theta*s*s+l*l*(t-cer_min)-s*s*t+s*s*(t-cer_min))/(sqrt(2)*s*l*sqrt(l*l+s*s))) - ...
         exp(-(tsk+theta-t)*(tsk+theta-t)/2/(l*l+s*s)) .* ...
          erf((-tsk*l*l+theta*s*s-s*s*t)/(sqrt(2)*s*l*sqrt(l*l+s*s))));


%sqrt(pi / 2) * s * a * norm_irf * norm_s * ...
%        (exp((-s^2-2*t*t_d+2*t_d*tsk)/(2*t_d^2)).*(erf((tsk+s^2/t_d)/(sqrt(2)*s))+erf((t_d*(t-tsk)-s^2)/(sqrt(2)*s*t_d)))- ...
%        exp((-s^2-2*t*t_r+2*t_r*tsk)/(2*t_r^2)).*(erf((tsk+s^2/t_r)/(sqrt(2)*s))+erf((t_r*(t-tsk)-s^2)/(sqrt(2)*s*t_r))))+ ...
%        1.0 / (sqrt(l*l+s*s)) * sqrt(pi / 2) * l * s * b * norm_c *...
%        (exp(-(tsk+tsk_c-t).*(tsk+tsk_c-t)/2/(l*l+s*s)).* ...
%        erf((-tsk*l*l+tsk_c*s*s+l*l*(t-tsk_c+4*l)-s*s*t+s*s*(t-tsk_c+4*l))/(sqrt(2)*s*l*sqrt(l*l+s*s))) - ...
%        exp(-(tsk+tsk_c-t).*(tsk+tsk_c-t)/2/(l*l+s*s)).* ...
%        erf((-tsk*l*l+tsk_c*s*s-s*s*t)/(sqrt(2)*s*l*sqrt(l*l+s*s))));
 
x=0*s:binning:100e-9;
%figure;
%hold on;
%plot(x,FS(x)); 
%plot(x,gradient(FS(x)),'red');

FSnum=FS(x);
dFSnum=-gradient(FSnum)/binning;

% figure;
% hold on;
% plot(x,FS(x,skew)); 
% plot(x,dFSnum,'red');

gpd=find(FSnum(1:length(FSnum))~=0);
I=sum(1./FSnum(gpd).*dFSnum(gpd).^2*binning); %Fisher information
CTRscan(i,n)=sqrt(1/I*1/(LY+CY))*3.33*1e12; %CTR
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
