close all hidden;
clear all;

t_d_min = 1e-9;
t_d_max = 1001e-9;
t_d_bin = 100;
t_d_pitch = (t_d_max - t_d_min)/t_d_bin; % 100e-9

t_r_min = 1e-12;
t_r_max = 301e-12;
t_r_bin = 30;
t_r_pitch = (t_r_max - t_r_min)/t_r_bin; % 10e-12

s_min = 1e-12;
s_max = 501e-12;
s_bin = 50;
s_pitch = (s_max - s_min)/s_bin; % 10e-12;

LY_min = 100;
LY_max = 10100;
LY_bin = 20;
LY_pitch = (LY_max - LY_min)/LY_bin; % 500;
CY = 30;

l = 10.0e-12; % cerenkov sigma
theta = 200e-12; % gamma time, parameter CR LB
tsk = 2e-9; % media trans
cer_min = theta; % cer min
cer_mean = theta + 2*l; % cer min

CTRscan = zeros(t_r_bin, t_d_bin, s_bin, LY_bin);

for i=1:t_r_bin
disp('Evento i : ');
disp(i);
for n=1:t_d_bin
%disp('Evento n : ');
%disp(n);
for l=1:s_bin
%disp('Evento l : ');
%disp(l);
for m=1:LY_bin
%disp('Evento m : ');
%disp(m);
binning = 0.1e-12;
LY = LY_min + m*LY_pitch; % Light Yield 
t_d = t_d_min + n*t_d_pitch; % decay time
t_r = t_r_min + i*t_r_pitch; % rise time
s = s_min + l*s_pitch; % sigma trans

a = LY / CY * (1 / (1 + LY/CY)); % normalization light yield
b = 1 / (1 + LY / CY); % normalization cerenkov yield

%norm_s = 1.0 / (t_d-t_r); % normalization shao
norm_s = 1.0 / (t_d-t_r); % normalization shao

%norm_irf = 2 / sqrt(pi) / sqrt(2) / s / erfc(- tsk / s / sqrt(2)); % normalization gaussian
norm_irf = 1. / (s * sqrt(pi / 2) * (1 + erf(tsk / sqrt(2) / s)));

%norm_c =  1. / l / sqrt(2. * pi) * (1. + l / 2. * sqrt(2. * pi) * erfc(1. / sqrt(2) / l * (l + 4.*tsk_c))); % cerenkov normalization 
norm_c = 1. / (l * sqrt(pi / 2) * (1 + erf((cer_mean - cer_min) / sqrt(2) / l)));

FS=@(t) (1.0/(sqrt(l*l+s*s)))*sqrt(pi/2)*l*s*norm_c*b*norm_irf*...
        (exp(-(tsk+cer_mean-t).*(tsk+cer_mean-t)/2/(l*l+s*s)).*...
          (erf((-tsk*l*l+cer_mean*s*s+l*l*(t-cer_min)-s*s*t+s*s*(t-cer_min))/(sqrt(2)*l*s*sqrt(l*l+s*s)))-...
          erf((-tsk*l*l+cer_mean*s*s-s*s*t)/(sqrt(2)*l*s*sqrt(l*l+s*s)))))+...
        sqrt(pi/2)*s*norm_irf*norm_s*a*...
        (exp((s*s-2*t*t_d+2*t_d*tsk+2*t_d*theta)/(2*t_d*t_d)).*...
	  (erf((t_d*(t-theta-tsk)-s*s)/(sqrt(2)*s*t_d))+...
          erf((t_d*tsk +s*s)/(sqrt(2)*s*t_d)))-...
        exp((s*s-2*t*t_r+2*t_r*tsk+2*t_r*theta)/(2*t_r*t_r)).*...
          (erf((t_r*(t-theta-tsk)-s*s)/(sqrt(2)*s*t_r))+...
	  erf((t_r*tsk +s*s)/(sqrt(2)*s*t_r))));

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
CTRscan(i,n,l,m)=sqrt(1/I*1/(LY+CY))*3.33*1e12; %CTR

end
%CTRscan
end
end
end
save('cramer_tot.mat', 'CTRSCAN');