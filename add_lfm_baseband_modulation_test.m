c = get_c();

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


%  SET UP BASEBAND WAVEFORM

t=linspace(0,PRI,NSAMPPRI);

%   simple pulse envelope

s=zeros(1,NSAMPPRI);    s(1:NSAMPPULSE+1)=1;




%   LFM
BW=5e6; PHASESTART=0;   CHIRPDIR=+1;

radar_signal.t = t;
radar_signal.voltage = s;
radar_signal = add_lfm_baseband_modulation(radar_signal,PW,BW,PHASESTART,CHIRPDIR);

s=s.*exp(1j*CHIRPDIR*BW/PW*(t-PHASESTART).^2);


assert(all(radar_signal.voltage == s))