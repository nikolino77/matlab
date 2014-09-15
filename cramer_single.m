close all hidden;
clear all;

CY = 3;

l = 10.0e-12; % cerenkov sigma
theta = 2000e-12; % gamma time, parameter CR LB
tsk = 2e-9; % media trans
cer_min = theta; % cer min
cer_mean = theta + 2*l; % cer min

binning = 0.1e-12;
LY = 4000; % Light Yield 
t_d = 40e-09; % decay time
t_r = 1e-12; % rise time
s = 70e-12; % sigma trans

a = LY / CY * (1 / (1 + LY/CY)); % normalization light yield
b = 1 / (1 + LY / CY); % normalization cerenkov yield

%a = 1.;
%b = 0.;

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

figure;
hold on;
plot(x,FS(x)); 
plot(x,gradient(FS(x)),'red');

FSnum=FS(x);
dFSnum=-gradient(FSnum)/binning;
gpd=find(FSnum(1:length(FSnum))~=0);
I=sum(1./FSnum(gpd).*dFSnum(gpd).^2*binning); %Fisher information
CTR=sqrt(1/I*1/(LY+CY))*3.33*1e12; %CTR

 SS=@(t) exp((s*s-2*t*t_r+2*t_r*tsk+2*t_r*theta)/(2*t_r*t_r));
  disp(exp((0.5*s*s/t_r/t_r-(1e-9-theta-tsk)/(t_r))));
  root = 6.75722e+18;
  disp((s*s-2*1e-9*t_r+2*t_r*tsk+2*t_r*theta)/(2*t_r*t_r));
x=0*s:binning:1e-9;
%disp(SS(x))
figure;
hold on;
plot(x,SS(x)); 
plot(x,gradient(SS(x)),'red');

CS=@(t) (1.0/(sqrt(l*l+s*s)))*sqrt(pi/2)*l*s*norm_c*b*norm_irf*...
        (exp(-(tsk+cer_mean-t).*(tsk+cer_mean-t)/2/(l*l+s*s)).*...
          (erf((-tsk*l*l+cer_mean*s*s+l*l*(t-cer_min)-s*s*t+s*s*(t-cer_min))/(sqrt(2)*l*s*sqrt(l*l+s*s)))-...
          erf((-tsk*l*l+cer_mean*s*s-s*s*t)/(sqrt(2)*l*s*sqrt(l*l+s*s)))));

x=0*s:binning:100e-9;

figure;
hold on;
plot(x,CS(x)); 
plot(x,gradient(CS(x)),'red');

S=@(t)  norm_s*(exp(-(t-theta)/t_d) - exp(-(t-theta)/t_r));

x=theta:binning:100e-9;

figure;
hold on;
plot(x,S(x)); 
plot(x,gradient(S(x)),'red');

IRF=@(t)  norm_irf*exp(-(t-tsk).*(t-tsk)/s/s/2.);

x=theta:binning:100e-9;

figure;
hold on;
plot(x,IRF(x)); 
plot(x,gradient(IRF(x)),'red');


disp(CTR);
