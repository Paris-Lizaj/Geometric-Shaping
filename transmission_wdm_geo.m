clear all;
close all;
clc;
format long;

Power_dBm=[-13:1:-6];   % This is the input power that you wan't to use (dBm)
%Power_dBm=[-8 -9]
Power_Watt=10.^((Power_dBm-30)/10)/10;     % This is the input power that you wan't to use (W)
Power_Watt2=1*Power_Watt;   % for wdm channels . Se to 0*Power_Watt for single channel transmission   
Spans=[10:1:20]; %that means that we take into account Spans*L total transmission length
%Spans=[20]
Nt=2^19; %%%%%%%%%% number of samples  (You can change this for more or less transmitted symbols)

num_of_total_channels =9;  % Please do not change this for now  
M=16; % M-QAM modulation format for the main channel
M2=16; % M-QAM modulation format for the wdm channels

downDBP=4;   %% Samples per symbol for Digital Backing Propogation 
dt=4e-11/32; %%%%% dt: sampling time  na ginei 32
dtDBP=dt*downDBP;

Tw=dt*Nt;  %%%%%%%%%% time window
NtDBP=Nt/downDBP;


Amp=4;
GreyCode={};

% 16=8+8
theta1=[0:(1/8):1-1/8];
realconst=[sin(theta1*2*pi)*Amp sin(theta1*2*pi+pi/8)*Amp/2]; 
imagconst=[cos(theta1*2*pi)*Amp cos(theta1*2*pi+pi/8)*Amp/2]; 
const=realconst+ 1i*imagconst; const=const.';
mygrey=[11 3 16 8 10 2 9 1 12 4 15 7 13 5 14 6];
coordinates=[];
for k=1:M
    temp=const(mygrey(k));
    coordinates=[coordinates temp];
end
%coordinates=coordinates.';
GreyCode{1}=coordinates;

%16=4+4+4+4   (-8)
theta1=[0:(1/4):1-1/4];
realconst=[sin(theta1*2*pi+pi/4)*Amp sin(theta1*2*pi)*Amp*3/4 sin(theta1*2*pi+pi/4)*Amp*2.2/4 sin(theta1*2*pi)*Amp*1.1/4]; 
imagconst=[cos(theta1*2*pi+pi/4)*Amp cos(theta1*2*pi)*Amp*3/4 cos(theta1*2*pi+pi/4)*Amp*2.2/4 cos(theta1*2*pi)*Amp*1.1/4]; 
const=realconst+ 1i*imagconst; const=const.';
mygrey=[13 16 14 15 9 12 10 11 5 8 6 7 1 4 2 3];
coordinates=[];
for k=1:M
    temp=const(mygrey(k));
    coordinates=[coordinates temp];
end
%coordinates=coordinates.';
GreyCode{2}=coordinates;

%16=1+6+9
theta1=[0:(1/6):1-1/6];
realconst1=sin(theta1*2*pi+3.1*pi/6)*Amp*0.9/2; 
imagconst1=cos(theta1*2*pi+3.1*pi/6)*Amp*0.9/2; 
theta2=[0:(1/9):1-1/9];
realconst2=sin(theta2*2*pi)*Amp; 
imagconst2=cos(theta2*2*pi)*Amp; 
realconst=[0 realconst1 realconst2]; 
imagconst=[0 imagconst1 imagconst2]; 
const=realconst+ 1i*imagconst; const=const.';
mygrey=[1 2 6 7 8 16 14 15 4 3 5 11 9 10 13 12];
coordinates=[];
for k=1:M
    temp=const(mygrey(k));
    coordinates=[coordinates temp];
end
%coordinates=coordinates.';
GreyCode{3}=coordinates;


% %16=4+8+4
theta1=[0:(1/4):1-1/4];
realconst1=sin(theta1*2*pi+pi/4)*Amp*1/4; 
imagconst1=cos(theta1*2*pi+pi/4)*Amp*1/4; 
theta2=[0:(1/8):1-1/8];
realconst2=sin(theta2*2*pi+pi/8)*Amp*2.4/4; 
imagconst2=cos(theta2*2*pi+pi/8)*Amp*2.4/4; 
%theta1
realconst3=sin(theta1*2*pi+pi/4)*Amp; 
imagconst3=cos(theta1*2*pi+pi/4)*Amp; 
realconst=[realconst1 realconst2 realconst3]; %realconst=realconst.';
imagconst=[imagconst1 imagconst2 imagconst3]; %imagconst=imagconst.';
const=realconst+ 1i*imagconst; const=const.';
mygrey=[4 11 1 6 12 16 5 13 3 10 2 7 9 15 8 14];
coordinates=[];
for k=1:M
    temp=const(mygrey(k));
    coordinates=[coordinates temp];
end
GreyCode{4}=coordinates;theta1=[0:(1/4):1-1/4];

clear temp coordinates const realconst imagconst mygrey ;
%%
%Wavs=4; % number of wavelengths
df=50e9; % wavelength grid = 50 GHz = 0.4 nmout (channels space (GHz))  
Symbol_rate=25e9; %%% symbol rate
samples_symbol=1/Symbol_rate/dt; %%%% samples per symbol
symbols=ceil(Nt/samples_symbol*2/3); %%% number of symbols

%Fourier properties
dw=1/dt/Nt*2*pi;

dwDBP=1/dt/Nt*2*pi;

wmax=1/dt*pi-dw;
wmin=-1/dt*pi;
w=wmin:dw:wmax;

wmaxDBP=1/dtDBP*pi-dwDBP;
wminDBP=-1/dtDBP*pi;
wDBP=wminDBP:dwDBP:wmaxDBP;

%%%%%%%%%%%%% transmission data%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

L=50e3;  %%%%%%%%%% transmission distance of each span in m
b2=-20e-27; %%%%%%%%%%% second order dispersion in s^2/m 
b3=0; %%%%%%%%%% third order dispersion in s^3/m
g=1.3e-3; %%%%%%%%%%%%%%%% nonlinear parameter in 1/W/m, you can put it to zero in the beginning
a=0.2/4.343/1000; %%%%%%%%%%%% losses in m^-1 0.34 dB for O-band
G=exp(a*L); %%%%%%%%%%%%% amplifier gain: usually equal to overall losses
dz=0.1e3; %%%%%%%%%%%% simulation step for the NLSE method

Nz=L/dz; %%%%%%%%%%% number of steps

%%%%%%%%%%%%%%%%%%%%%% generation of the OFDM signal %%%%%%%%%%%%%%%%

t=dt:dt:Tw; %%%%%%%%%%% time matrix
tDBP=dtDBP:dtDBP:Tw;

%% LOOP Input Power soans and geoshape
for igrey=1:length(GreyCode)
for iii=1:length(Power_dBm) % From -5 dBm to 3 dBm  with step 1 dB 
for ispan=1:length(Spans)


rng('shuffle');

for chan=1:2*num_of_total_channels

    data{chan}=randi([0,M-1],1,symbols,'uint8');    
    %GreyCode{igrey}'
    mod_data{chan}=genqammod(data{chan},GreyCode{igrey}); %modulated data
    
    meanP{chan}=mean(abs(mod_data{chan}).^2);
   
    % Initially Ax = sqrt(Pow)
    A{chan}=ones(1, Nt)*sqrt(Power_Watt(iii));
    Aloc{chan}=A{chan};
   
    if (chan==1) && (iii==1) && (ispan==1)
        index=[0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15];
        scatterplot(mod_data{1})
        hold on
        for k = 1:M
            text(real(GreyCode{igrey}(k))-0.4,imag(GreyCode{igrey}(k))+0.4,num2str(index(k)));
            %num2str(data(k))
        end
        hold off
    end
%oversampling of modulation to analog sampling
for kk=1:length(cell2mat(mod_data(chan)))
    Resampled{chan}(samples_symbol*(kk-1)+1:samples_symbol*kk)=mod_data{chan}(kk);
end

Diff{chan}=length(A{chan})-length(Resampled{chan});

%final modulation
A{chan}(Diff{chan}/2:length(Resampled{chan})+Diff{chan}/2-1)=sqrt(Power_Watt(iii))*Resampled{chan};

end

Power1=mean(abs(A{1}(Diff{1}/2:length(Resampled{1})+Diff{1}/2-1)).^2); % print the actual mean power after modulation, this is the power you should include in your x axis. 
PdBm1=10*log10(1000*Power1);
% 


h=6.626070040e-34;
fo=228.85e12;
nsp=1.5;
Geq=10; % equivalent gain of an EDFA which determines the OSNR:  You change this parameter for noise tests

LPMD=dz;
Dp=0.0e-12/sqrt(1000); % PMD coefficient 0.1 ps/km^1/2 typical value for optical fibers
dPMD=Dp/sqrt(LPMD)*sqrt(3*pi/8)/2;

Leff=(1-exp(-a*L))/a;

NzDBP=2; % number of spans of the digital back-propagation
dzDBP=L/NzDBP;

Linear=exp(-a/2*dz+1i*(w.^2*b2*dz/2.0+w.^3*b3*dz/6.0));

LinearDBP=exp(a/2*dzDBP-1i*(wDBP.^2*b2*dzDBP/2.0+wDBP.^3*b3*dzDBP/6.0));

Linear_eq=exp(-1i*(w.^2*b2*Spans(ispan)*L/2.0+w.^3*b3*Spans(ispan)*L/6));


Ain=A;
 
Ax=A{1}+...
    +1*A{3}.*exp(-1i*df*2*pi*t)...
    +1*A{5}.*exp(1i*df*2*pi*t)...
    +1*A{7}.*exp(-1i*2*df*2*pi*t)...
    +1*A{9}.*exp(1i*2*df*2*pi*t)...
    +1*A{11}.*exp(-1i*3*df*2*pi*t)...
    +1*A{13}.*exp(1i*3*df*2*pi*t)...
    +1*A{15}.*exp(-1i*4*df*2*pi*t)...
    +1*A{17}.*exp(1i*4*df*2*pi*t) ;

Ay=A{2}...
    +1*A{4}.*exp(-1i*df*2*pi*t)...
    +1*A{6}.*exp(1i*df*2*pi*t)...
    +1*A{8}.*exp(-1i*2*df*2*pi*t)...
    +1*A{10}.*exp(1i*2*df*2*pi*t)...
    +1*A{12}.*exp(-1i*3*df*2*pi*t)...
    +1*A{14}.*exp(1i*3*df*2*pi*t)...
    +1*A{16}.*exp(-1i*4*df*2*pi*t)...
    +1*A{18}.*exp(1i*4*df*2*pi*t) ;

AllMean=10*log10(mean(abs(Ax).^2+abs(Ay).^2));

%% Transmission --- 
tic       
[Ax, Ay]= span(Spans(ispan), Nz, Ax, g, dz, Ay, Linear, G, h, fo, nsp, dt);
toc

% optical spectrum representation
%figure
%Axf=fft(Ax);
%Axf=fftshift(Axf);
%semilogy(w/2/pi, abs(Axf).^2, 'r');



Ax=butterworth(Ax, Nt, Tw, 10, Symbol_rate, 0); 
Ay=butterworth(Ay, Nt, Tw, 10, Symbol_rate, 0); 


% optical spectrum representation
%figure
%Axf=fft(Ax);
%Axf=fftshift(Axf);
%semilogy(w/2/pi, abs(Axf).^2, 'r');


AxDBP=downsample(Ax,downDBP);
AyDBP=downsample(Ay,downDBP);

AxDBP=AxDBP*sqrt(1/G);
AyDBP=AyDBP*sqrt(1/G);


%% dbp
tic
[AxDBP, AyDBP]= span(Spans(ispan), NzDBP, AxDBP, -g*0.0, dzDBP, AyDBP, LinearDBP, 1/G, h, fo, 0*nsp, dtDBP);
toc

AxDBP=butterworth(AxDBP, Nt/downDBP, Tw, 10, Symbol_rate, 0); 
AyDBP=butterworth(AyDBP, Nt/downDBP, Tw, 10, Symbol_rate, 0); 


Axout=AxDBP;
Ayout=AyDBP;

%% coherent receiver
realx=real(Axout);
imagx=imag(Axout);

realy=real(Ayout);
imagy=imag(Ayout);


out_complex_x=(realx+1i*imagx);
out_complex_y=(realy+1i*imagy);


%% xcorrelation for lag 
[xC,lagx]=xcorr(Ain{1},Axout);
xC=xC/max(xC);
[Mcx,Icx]=max(xC);
t_cor_x=lagx(Icx);

[yC,lagy]=xcorr(Ain{2},Ayout);
yC=yC/max(yC);
[Mcy,Icy]=max(yC);
t_cor_y=lagy(Icy);


Diff_x=Diff{1}/downDBP;
Diff_y=Diff{2}/downDBP;
samples_symbol2=samples_symbol/downDBP;

for kk=1:length(mod_data{1})
out_data_x(kk)=out_complex_x(Diff_x/2+samples_symbol2*(kk-1)+samples_symbol2/2-0*t_cor_x);
out_data_y(kk)=out_complex_y(Diff_y/2+samples_symbol2*(kk-1)+samples_symbol2/2-0*t_cor_y);
end
%%

out_data_x=sqrt(meanP{1})*out_data_x/sqrt(mean(abs(out_data_x).^2));
out_data_y=sqrt(meanP{2})*out_data_y/sqrt(mean(abs(out_data_y).^2));


%% constallation rotation for find best BER
for k=1:360
out_data_shift_x=out_data_x*exp(-k*1i*2*pi/360);
out_data_shift_y=out_data_y*exp(-k*1i*2*pi/360);


Rec_data_x=genqamdemod(out_data_shift_x, GreyCode{igrey});
Rec_data_y=genqamdemod(out_data_shift_y, GreyCode{igrey});


[NNx, RRx]=biterr(data{1}, Rec_data_x);
[NNy, RRy]=biterr(data{2}, Rec_data_y);


BER_x(k)=RRx; % teliko bit error rate
BER_y(k)=RRy; % teliko bit error rate

end

[minBERx,Ix] = min(BER_x);
minBERx
[minBERy,Iy] = min(BER_y);
minBERy

sprintf('minBERx= %f', minBERx);
sprintf('minBERy= %f', minBERy);
%scatterplot(out_data_x*exp(-Ix*1i*2*pi/360))
%scatterplot(out_data_y*exp(-Iy*1i*2*pi/360))



%% Export data
%3 for geo qam
%igrey for each grey 1 2 3 4
writtenData=[3 symbols igrey Spans(ispan) Power_dBm(iii) AllMean minBERx minBERy];
dlmwrite('PowerBERgeo.txt',writtenData, '-append');

X_x_real=real(out_data_x*exp(-Ix*1i*2*pi/360))';
X_x_imag=imag(out_data_x*exp(-Ix*1i*2*pi/360))';

X_y_real=real(out_data_y*exp(-Iy*1i*2*pi/360))';
X_y_imag=imag(out_data_y*exp(-Iy*1i*2*pi/360))';

X_x=[X_x_real, X_x_imag];
X_y=[X_y_real, X_y_imag];

str=sprintf("_%d_%d_%d_%d.txt",M,Spans(ispan),Power_dBm(iii),igrey);
% EXPORTS central
filenameXx = "X_x"+str;
csvwrite(filenameXx,X_x)

filenameYx ="Y_x"+str;
csvwrite(filenameYx,Rec_data_x')

filenameXy ="X_y"+str;
csvwrite(filenameXy,X_y)

filenameYy = "Y_y"+str;
csvwrite(filenameYy,Rec_data_y')


    
clear data;
clear mod_data;
clear mod_data_rrc;
clear meanP;
clear A;
clear Aloc;
clear Resampled;
clear Diff;

clear Ain;
clear Aout;
clear imagn;
clear out_complex;
clear realn;
clear lag;
clear imagpart;
clear realpart;

clear Ax AxDBP Axout Ay AyDBP Ayout;
clear BER_x BER_y chan Diff_x Diff_y Dp dPMD;
clear fo h Icx Icy imagximagy Ix Iy k kk lagx Leff Linear LinearDBP LPMD;
clear Mcx Mcy minBERx minBERy minimum_x minimum_y NNx NNy nsp;
clear out_complex_x out_complex_y out_data_shift_x out_data_shift_y out_data_x out_data_y out_data_y_rcc out_data_y_rcc;
clear Power_dBm_central realx realy Rec_data_x Rec_data_y RRx RRy str;
clear t_cor_x t_cor_y writtenData X_x X_x_imag X_x_real X_y X_y_imag X_y_real xC yC;
    
    
end
end
end
% clear xC;
% clear A;
% clear Ain;
% clear Aout;
% clear imagn;
% clear outcomplex;
% clear realn;
% clear lag;
% clear Ax;
% clear Axf;
% clear Ay;
% clear Ayf;
% clear Ain;
% clear Linear;
% clear Resampled;
% clear Aloc;
% clear imagpart;
% clear realpart;
% clear mod_data;
% clear out_data;
% clear out_data_shift;
% clear X;
% clear data;
% clear Rec_data;
% clear X_imag;
% clear X_real;
% clear BER;
% clear Mcx;
% clear Diff;
% clear I;
% clear Ic;
% clear meanP;
% clear meanPhi;
% clear minBER;
% clear NN;
% clear RR;
% clear t_cor;
% close all;
