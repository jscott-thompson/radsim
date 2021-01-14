%%
%   REVISED RADAR MODEL -- NO CLUTTER
%%

clear all;  close all;
tic
c=3.0e8;
d2r=180/pi;

%%  DEFINE WAVEFORM

f=1.0e9;
lambda=c/f;
k=2*pi/lambda;

PRF=50e3;
PRI=1/PRF;
DUTY=0.1;
PW=DUTY*PRI;
NPULSES=64;
TCPI=NPULSES*PRI;

%   fast-time sample space
NSAMPPRI=1000;
NSAMPPULSE=ceil(NSAMPPRI*DUTY);

%   misc derived quantities
RUNAMB=0.5*c*PRI;
RMIN=c*PW;
FDUNAMB=0.5*PRF;


%%  SET UP BASEBAND WAVEFORM

t=linspace(0,PRI,NSAMPPRI);

%   simple pulse envelope

s=zeros(1,NSAMPPRI);    s(1:NSAMPPULSE+1)=1;

%   LFM
BW=5e6; PHASESTART=0;   CHIRPDIR=+1;
s=s.*exp(1j*CHIRPDIR*BW/PW*(t-PHASESTART).^2);

%   CPI
S=ones(NPULSES,1)*s;

%%  SET UP TARGET

%   target velocity and range
FDTGT=[-0.5;0.33;0.5]*FDUNAMB;  
RTGT=[0.7;2.3;0.8]*RUNAMB;    
RCSTGT=[10;10;10];                         %   m^2

VTGT=0.5*lambda*FDTGT;
R4TGT=RTGT.^-4;

tgtPhase=zeros(NPULSES,1);
sdel=zeros(1,NSAMPPRI);
Y=zeros(NPULSES,NSAMPPRI);

for i=1:numel(FDTGT)
    %   phase progression of target across pulses
    tgtPhase=exp(1j*2*pi*FDTGT(i)*((0:NPULSES-1)'*PRI))*exp(-1j*2*k*RTGT(i));

    %   apply time delay to pulse
    TTGT=mod(RTGT(i),RUNAMB)/c;
    sdel=circshift(s,ceil(TTGT*NSAMPPRI/PRI));
    
    %   received signal reflected from a point target
    Y=Y+tgtPhase*sdel;
end

%   blank received signal during the pulse
Y(:,1:NSAMPPULSE)=0;  

%%  PULSE DOPPLER PROCESSING

X=zeros(NPULSES,2*NSAMPPRI-1);
X_test=zeros(NPULSES,NSAMPPRI);

%   match filter
for i=1:NPULSES
    X(i,:)=xcorr(Y(i,:),conj(s));
 
    %   test the transform method
    fy=fftshift(fft(Y(i,:)));
    fsp=fftshift(fft(conj(s)));
    X_test(i,:)=ifft(fy.*fsp);
 end

X(:,1:NSAMPPRI-1)=[];

%   Doppler process
for i=1:NSAMPPRI
    RD(:,i)=fftshift(fft(X(:,i)));
    RD_test(:,i)=fftshift(fft(X_test(:,i)));

end

%%  NORMALIZE AND ADD NOISE
RD=RD./max(abs(RD(:)));
RD_test=RD_test./max(abs(RD_test(:)));
SNRPK=10;                           %   dB
SNRPKL=10^(0.2*SNRPK);
N=sqrt(1/SNRPKL)*(randn(NPULSES,NSAMPPRI)+1j*randn(NPULSES,NSAMPPRI));
RD=RD+N;
RD=abs(RD);
RD_test=RD_test+N;
RD_test=abs(RD_test);

%%  ID DETECTION CANDIDATES USING CFAR ALGORITHM


[Thresh,NumCells]=getCFARLevel(RD,6,3);

% PD=0.5; PFA=0.11;   
% SINRCA=((PD/PFA).^(1./NumCells)-1)./(1-PD.^(1./NumCells));  %   the SINR required to achieve PFA/PD spec  

SINRCA=10^(0.1*6);
PFA=exp(-1*10*log10(SINRCA)/SNRPKL);
DetCands=RD./Thresh>SINRCA;
DetCell=find(DetCands==1);
[CT,VEL]=meshgrid(c*t,0.5*lambda*PRF*linspace(-0.5,0.5,NPULSES));
DetRange=CT(DetCell);
DetVel=VEL(DetCell);

    %%
    subplot(2,4,[1,2,5,6]);
    imagesc(CT(1,:),VEL(:,1),10*log10(RD)); title('Range-Doppler'); hold on;
    plot(DetRange,DetVel,'rx')
    ylabel('Velocity [m/s]');
    xlabel('Range [m]');
%     caxis([10 40])
    subplot(2,4,3);
    imagesc(c*t,0.5*lambda*PRF*linspace(-0.5,0.5,NPULSES),10*log10(SINRCA)); title('SINR Requirement');
    ylabel('Velocity [m/s]');
    xlabel('Range [m]');
%     caxis([10 40])
    subplot(2,4,4);
    imagesc(c*t,0.5*lambda*PRF*linspace(-0.5,0.5,NPULSES),Thresh); title('Cell Average');
    ylabel('Velocity [m/s]');
    xlabel('Range [m]');
    subplot(2,4,7);
    imagesc(c*t,0.5*lambda*PRF*linspace(-0.5,0.5,NPULSES),DetCands); title('Candidate Detections');
    ylabel('Velocity [m/s]');
    xlabel('Range [m]');
    %%

%%  CLUSTER THRESHOLD EXCEEDANCES TO FINALIZE DETECTIONS

EUCLDIST=zeros(size(DetRange));

%   normalize the range/velocity coordinates
DetRange0=DetRange./max(DetRange);
DetVel0=DetVel./max(DetVel);

%   compute the euclidian distance between each exceedance
for i=1:numel(DetRange0)
    for j=1:numel(DetVel0)        
        EUCLDIST(i,j)=sqrt((DetRange0(i)-DetRange0(j))^2+(DetVel0(i)-DetVel0(j))^2);
    end
end

%   cluster 
distThresh=0.5;
EUCLDISTBIN=EUCLDIST<distThresh;
DETECTIONS=[];
i=1;    j=1;

while i<numel(DetRange0)
    DETTEMP=find(EUCLDISTBIN(i,:)>0);
    DETECTIONS(j).inds=DETTEMP;
    DETECTIONS(j).r=DetRange(DETTEMP);
    DETECTIONS(j).d=DetVel(DETTEMP);
    DETECTIONS(j).centroid=[mean(DETECTIONS(j).r) mean(DETECTIONS(j).d)];
    j=1+j;
    i=DETTEMP(end)+1;
end

NUMDET=numel(DETECTIONS);
cmap=jet(NUMDET);
        
    %%
    subplot(2,4,[1,2,5,6]);
    for i=1:NUMDET
        plot(DETECTIONS(i).r,DETECTIONS(i).d,'o','markeredgecolor',cmap(i,:))
    end
    
    ylabel('Velocity [m/s]');
    xlabel('Range [m]');    
    
    subplot(2,4,8)
    imagesc(CT(1,:),VEL(:,1),10*log10(RD)); title('Finalized Detections'); hold on;
    for i=1:NUMDET
        plot(DETECTIONS(i).centroid(1),DETECTIONS(i).centroid(2),'v','markeredgecolor',cmap(i,:))
    end
    ylabel('Velocity [m/s]');
    xlabel('Range [m]');
    %%
    toc