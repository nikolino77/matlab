close all hidden;
clear all;

Matrixbinning=0.02;
Msize=2*1;
ErgName='test.txt';
FErg=fopen(ErgName,'w');

for n=1:Msize/Matrixbinning+1
for i=1:Msize/Matrixbinning+1

isc=(i-Msize/2/Matrixbinning-1)*Matrixbinning;
nsc=(n-Msize/2/Matrixbinning-1)*Matrixbinning;

t_d=30.3e-9*10^nsc;
t_r=70e-12*10^isc;
fprintf(FErg,'%g\n',isc);

end
%fprintf(FErg,'%g\n',t_d);
end
