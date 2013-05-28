% Pulsar timing array data simulator
% Yan Wang, March 8, 2013

% Note:
% currently white Gaussian noise only for this function

function [alpha,delta,omega,phi0,phiI,pname] =PTAsimulator()

% ======  Useful constants  ======
pc2ly=3.261563777;  % 1 pc=3.26 ly (Julian)
dy2yr=1.0/365.25;  % 1 day=365.25 yr (Julian)  ??
kilo=1.0*10^3;  % kilo 1000

% set global variables for functions need them
global Np alphaP deltaP kp N timingResiduals sd yr


% ===============================
% Constructing a pulsar timing array using Np pulsars
Np=8;  % number of pulsars in the timing array

% sky location of the pulsars in the equatorial coordinate
% we need to transfer from hr angle and degree to radian
alphaP=zeros(Np,1);  % right ascension, in radian
deltaP=zeros(Np,1);  % declination, in radian
distP=zeros(Np,1);  % (parallax) distance from SSB to pulsars, from mas/pc to ly
kp=zeros(Np,3);  % unit vector pointing from SSB to pulsars,


% -------------------------------
% transfer (hr,min,sec) and (degree,min,sec) to radian; mas/kpc to ly
tmp1='J0030+0451';
alphaP(1)=(0*15+30*15/60)*pi/180;
deltaP(1)=(4+51/60)*pi/180;
distP(1)=1.376*kilo*pc2ly;  % in ly

tmp2='J0613-0200';
alphaP(2)=(6*15+13*15/60)*pi/180;
deltaP(2)=-(2+0/60)*pi/180;
distP(2)=6.318*kilo*pc2ly;

tmp3='J1713+0747';
alphaP(3)=(17*15+13*15/60)*pi/180;
deltaP(3)=(7+47/60)*pi/180;
distP(3)=7.524*kilo*pc2ly;

tmp4='J1909-3744';
alphaP(4)=(19*15+9*15/60)*pi/180;
deltaP(4)=-(37+44/60)*pi/180;
distP(4)=3.532*kilo*pc2ly;

%-----------
tmp5='J1012+5307';
alphaP(5)=(10*15+12*15/60)*pi/180;
deltaP(5)=(53+7/60)*pi/180;
distP(5)=1.045*kilo*pc2ly;

tmp6='J1455-3330';
alphaP(6)=(14*15+55*15/60)*pi/180;
deltaP(6)=-(33+30/60)*pi/180;
distP(6)=6.593*kilo*pc2ly;

tmp7='J1600-3053';
alphaP(7)=(16*15+0*15/60)*pi/180;
deltaP(7)=-(30+53/60)*pi/180;
distP(7)=13.532*kilo*pc2ly;

tmp8='J1640+2224';
alphaP(8)=(16*15+40*15/60)*pi/180;
deltaP(8)=(22+24/60)*pi/180;
distP(8)=3.675*kilo*pc2ly;

pname={tmp1 tmp2 tmp3 tmp4 tmp5 tmp6 tmp7 tmp8};

% sky location of pulsars in Cartesian coordinate
for i=1:1:Np
    kp(i,1)=cos(deltaP(i))*cos(alphaP(i));
    kp(i,2)=cos(deltaP(i))*sin(alphaP(i));
    kp(i,3)=sin(deltaP(i));
end


% ===============================
% Simulate the post-fit residuals, composing of the GW contributions and 
% other noises (white, red, Gaussian, non-Gaussian, systematic ...)

% sky location of the single source in the equatorial coordinate
alpha=pi/4+1.2;  % right ascension [0, 2*pi]
delta=pi/4-0.16;  % declination [-pi/2, pi/2]
%dist=3.12*10^6;  % distance from SSB to source, in pc

% sky location in Cartesian coordinate
k=zeros(1,3);  % unit vector pointing from SSB to source
k(1)=cos(delta)*cos(alpha);
k(2)=cos(delta)*sin(alpha);
k(3)=sin(delta);

% other parameters for the source
period=0.3925;  % period of GW measured in the SSB, in unit of yr
omega=2*pi/period;  % angular frequency of GW, two times of orbital freq.
iota=pi/4+0.6;  % inclination between orbital plane and plane of the sky
thetaN=pi/4+0.3;  % angle to the line of nodes
phi0=1.6;  % initial orbital phase
% mu=1.42*10^8;  % reduced mass for the binary, in solar mass


% -------------------------------
% starting epoch of the observations
start=53187;  % Modified Julian Day, 'July 1, 2004'
%finish=start+ceil(365.25*5);  % approximately, set 5 yrs
deltaT=14;  % observation cadence, in days, set biweekly
N=128;  % number of biweekly observations for all pulsars, fft prefer 2^n
dy=zeros(N,1);  % observation epoch, in day
yr=zeros(N,1);  % observation epoch, in year
for i=1:1:N
    dy(i)=start+(i-1)*deltaT;  % dates conducting observations, in MJD
    yr(i)=2004.5+(i-1)*deltaT*dy2yr;
end


% -------------------------------
% calculate timing residuals induced by GW for each pulsar
Amp=1.3*10^(-7);  % overall amplitude of timing residuals, the same for all pulsars, sec
timingResiduals=zeros(Np,N);  % signal, i.e. GW induced timing residuals, Np pulsars, N observations
phiI=zeros(Np,1);  % arbitrary phase for each pulsar, relative distance

for i=1:1:Np

    theta=acos(k*kp(i,:)');
    %sprintf('%d pulsar theta=%g',i,theta)
    phiI(i)=mod(phi0-omega*distP(i)*(1-cos(theta)), 2*pi);  % modulus after division
    
    timingResiduals(i,:)=FullResiduals(alpha,delta,omega,phi0,phiI(i),alphaP(i),deltaP(i),...
                         Amp,iota,thetaN,theta,yr);

end


disp('True Parameters')
disp([alpha,delta,omega,phi0,phiI']');

% -------------------------------
% plot the timing residuals for each pulsar
figure
for i=1:1:Np
    subplot(4,2,i)
    plot(yr,timingResiduals(i,:),'.-');
    grid on;
    %hold on
    %plot(dy,re(i,:),'r.-');
    %xlabel('Modified Juliant Day');
    xlabel('years');
    ylabel('noise free (sec)');
    title(pname(i));
end

% -------------------------------
% noise models for each pulsar, make it a function, more realistic noise
% in the future...

%noise=zeros(Np,N);  % noise
sd=0.2*10^(-7);  % standard deviation of the normal distribtion (sec)

% we fix the noise realization here, #*10^(-8)
%tmp=load('fixNoise1.mat');  % sd = 10^(-8) sec
tmp=load('fixNoise2.mat');
%tmp=load('fixNoise5.mat');
%tmp=load('fixNoise10.mat');

noise=tmp.fixNoise;
for i=1:1:Np
    %noise(i,:)=sd*randn(1,N);  % Gaussian noise
    timingResiduals(i,:)=timingResiduals(i,:)+noise(i,:);
end

% -------------------------------
% plot the timing residuals for each pulsar
figure
for i=1:1:Np
    subplot(4,2,i)
    plot(yr,timingResiduals(i,:),'.-');
    grid on;
    %hold on
    %plot(dy,re(i,:),'r.-');
    %xlabel('Modified Juliant Day');
    xlabel('years');
    ylabel('residuals with noise (sec)');
    title(pname(i));
end

% -------------------------------
% calculate the signal-noise ratio of some kind... sqrt(Ps/Pn);
snr=zeros(Np,1);
temporal1=0.0;
temporal2=0.0;

for i=1:1:Np
    for j=1:1:N
        temporal1=temporal1+timingResiduals(i,j)^2;
        temporal2=temporal2+noise(i,j)^2;
    end
    
    snr(i,1)=sqrt(temporal1/temporal2);
    temporal1=0.0;
    temporal2=0.0;
    
end

overallSNR=0.0;
for i=1:1:Np
    overallSNR=overallSNR+snr(i,1);
end
overallSNR=overallSNR/Np;

disp('signal to noise ratio for each pulsar')
disp(snr');
disp('overall SNR: ') 
disp(overallSNR)

%{

% -------------------------------
% histogram for the noise
figure
for i=1:1:Np
    subplot(2,2,i)
    hist(noise(i,:));
end


%}


% END of function
