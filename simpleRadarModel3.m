
total_run_timer = tic;

% Load packages and constants
% Uncomment the next 2 lines to run with Octave
% pkg load dataframe;
% pkg load signal;
rad_clear_vars;
rad_load_constants;

% Set random seed
if shuffle_seed
    rng_settings = rng('shuffle');
else
    rng_settings = rng(0);
end

% Set output file
logfile = "rad_output.csv";

% Load run matrix file
run_matrix_file = "final_run_matrix.csv";
params = rad_load_run_matrix(run_matrix_file);
N_runs = size(params,1);
results = zeros(N_runs,6);

max_targets_det = 0;
run_time = 0;
f_waitbar = waitbar(0,['Case 0 of ' num2str(N_runs)]);

if ~exist('result','var')
    result = cell(N_runs,3);
end

runs = 1:N_runs;
runs(2197) = [];

% Loop over runs
for i_run = runs
%     if params(i_run,N_TARGETS) == 0
%         continue;
%     end
    close all;
    waitbar(i_run/N_runs,f_waitbar,['Case ' num2str(i_run) ' of ' num2str(N_runs) ': ' num2str(run_time) ' s; ETC: ' num2str(N_runs*toc(total_run_timer)/i_run/60) ' min']);
    this_run_timer = tic;

    vr=params(i_run,RADAR_VEL)*Nmi2m*3600^-1;   %   velocity of the platform, m/sec
    vvec=[1;0;0];           %   velocity vector of the platform
    hr=params(i_run,RADAR_ALT)*f2m;             %   altitude of the radar
    x0=[0;0;hr];            %   position vector of the radar
    fr=params(i_run,FREQUENCY);                 %   frequency of the radar 
    lambda=c*fr^-1;         %   wavelength of the radar 
    k=2*pi*lambda^-1;       %   wavenumber

    Ny=10;  Nz=10;          %   number of antenna elements in y and z directions
    Nel=Ny*Nz;              %   total number of antenna elements 
    dyz=0.5*lambda;         %   spacing between row/column in rectangular array

    [Array]=getArrayElements(Ny,Nz,dyz);  %   Nel x 3 matrix of element positions

    % Waveform
    PRF=params(i_run,I_PRF);               %   pulse repetition frequency 
    PRI=PRF^-1;             %   pulse repetition interval
    tau=params(i_run,DUTY_FACTOR);                %   duty factor
    pw=tau/PRF;             %   pulse width
    Npulses= params(i_run,N_PULSES);             %   number of pulses in a CPI
    NFFT=Npulses;               %   number of points in the FFT **** Needs to match number of pulses?
    dr=c*pw;                %   uncompressed range resolution
    Nsamp=1000;             %   number of samples per PRI **** Needs to scale with frequency?


    %   set up a single pulse
    Tpri=linspace(0,PRI,Nsamp); 
    dt=Tpri(2)-Tpri(1);
    a0=zeros(1,numel(Tpri));    %   baseband waveform -- single pulse
    a0(Tpri<=pw)=1;
    s0=a0.*exp(1j*2*pi*fr*Tpri);    %   waveform at RF -- single pulse 

    %   expand into the CPI
    Tcpi=linspace(0,PRI*Npulses,Nsamp*Npulses);
    A0=ones(Npulses,1)*a0;  
    S0=ones(Npulses,1)*s0;

    %%  scene geometry
    Runamb=0.5*c*PRI;       %   unambiguous range of radar 
    Rlosmax=sqrt((Re+hr)^2-Re^2);   %   maximum line of sight range 
    RmaxSim=0.25*Rlosmax;   %   maximum range of the simulation

    %%  radar electronic scan angles 
    elsq=-10*d2r;   %   elevation (theta) 
    azsq=0*d2r;     %   azimith (phi)

    rsq=[cos(azsq)*cos(elsq);sin(azsq)*cos(elsq);sin(elsq)];    %   unit vector for the scanning
    phsq=exp(-1j*k*Array*rsq);   %   Nel x 1 vector of phase shifts for the scan direction

    %%  begin the analysis

    %%  1.  Set up the clutter
    if simulate_clutter_response
        %   ClutterMap parameters 
        % clut_dr=0.25*dr;    %   ClutterMap array spacing, m
        clut_dr=200;
        clut_ring=20*Nmi2m;  %   ClutterMap ring width about target range   
        rcs_std=params(i_run,CLUTTER_RCS_SD);

        %   get the ClutterMap  
        ClutterMap=getClutterArray(-1*x0(3)/sin(elsq),clut_dr,clut_ring);  %   Mx3 matrix of ClutterMap array points

        rcv=ClutterMap-ones(numel(ClutterMap(:,1)),1)*x0';   %   vector from radar to each clutter map point 
        rc=sqrt(rcv(:,1).^2+rcv(:,2).^2+rcv(:,3).^2);   %   range form radar to each clutter point
        rch=rcv./rc; %   unit vector from the radar phase center to each ClutterMap point

        %   calculate the clutter normalized RCS
        theta_ClutterMap=atan(-1*rch(:,3)./rch(:,1));   %   incident angle to each clutter point
        rcsc=abs(rcs_std*randn(numel(rc),1).*sin(theta_ClutterMap));                     %   normalized clutter rcs

        %%  2.  Get the array factor toward the clutter points
        AFc=phsq'*exp(-1j*k*Array*rch');    %   array gain toward the ClutterMap points
        AFc=AFc./max(abs(AFc(:)));         %   normalize the array gain

        AFc(sign(rch(:,1))<0)=1e-4; %   set the back-lobe gain to -40 dB

        %%  3.  Get the clutter power and reduce the computation size
        powc=abs(AFc').^2.*rcsc.*rc.^-4;
        %powc=abs(AFc').^2;
        powc=powc./max(powc(:));    %   normalize the clutter power

        %   plot
    % 	figure; scatter3(ClutterMap(:,1),ClutterMap(:,2),ClutterMap(:,3),10,10*log10(powc),'filled')
    % 	title('Pre-threshold');

        %   set the threshold for eliminating extra points, then do it
        thresh=1e-4;
        clutterReduce=powc<thresh;

        ClutterMap(clutterReduce,:)=[];
        AFc(clutterReduce)=[];
        powc(clutterReduce)=[];
        rcv(clutterReduce,:)=[];
        rc(clutterReduce)=[];
        rch(clutterReduce,:)=[];

        nptsc=numel(rc);    %   number of clutter points remaining

        %   plot
        % ##figure; scatter3(ClutterMap(:,1),ClutterMap(:,2),ClutterMap(:,3),10,10*log10(powc),'filled')
        % ##title('Post-threshold');

        %%  4.  get clutter range Doppler response 

        jitterVal=1e-4; 

        while numel(rc)~=numel(unique(rc))
    % 	    disp('jittering');
            rcv=rcv+randn(nptsc,3)*jitterVal; %   jitter the points to ensure that no two have the same slant range
            rc=sqrt(rcv(:,1).^2+rcv(:,2).^2+rcv(:,3).^2);   %   new slant range
            rch=rcv./rc; %   still need unit vectors
        end

        %   need to handle range ambiguous clutter
        if any(rc>0.5*c*PRI)
            rc_unamb=rc;    %   store the unambiguous clutter range for debugging purposes
            rc=mod(rc,0.5*c*PRI);   %   reset the ranges to be ambiguous
        end

        %   need to create a strictly monotonically increasing range grid, so sort
        [rc,order]=sort(rc);    
        rch=rch(order,:); 
        rcv=rcv(order,:);

        fdc=2*vr*lambda^-1*rch*vvec;    %   Doppler shift of each clutter patch
        phc=exp(1j*2*pi*lambda^-1*(fdc*ones(1,Npulses)).*((PRI*ones(numel(rc),1))...
            *(0:Npulses-1)+rc*ones(1,Npulses)));    %   complex phase of patches across slow time

        %   set up an interpolation to get back to the sample grid
        [Pulses,Rc]=meshgrid((0:Npulses-1),rc); %   make a grid of fast time, slow time cells
        [Pulses_q,Rc_q]=meshgrid((0:Npulses-1),0.5*c*Tpri); %   the query grid is the true sample space

        phc_grid_Re=interp2(Pulses,Rc,powc.*real(phc),Pulses_q,Rc_q); %   interpolate real part of phase to the true sample space
        phc_grid_Im=interp2(Pulses,Rc,powc.*imag(phc),Pulses_q,Rc_q); %   interpolate imag part of phase to the true sample space

        phc_grid=phc_grid_Re+1j*phc_grid_Im;
        phc_grid(isnan(phc_grid))=0;    %   remove any NaNs that may be generated 

        %   tidy up
        phic=phc_grid;
        Rc=Rc_q;

        clear phic_grid_Re phic_grid_Im phic_grid Rc_q Pulses_q

        %%  5. generate the phase history of the clutter scene by convolving with the waveform

        %   pre-allocate arrays for speed
        yc=zeros(numel(Tpri)+numel(s0)-1,Npulses);  %   this is the convolution array
        xc_est=zeros(2*(numel(Tpri)+numel(s0)-1)-1,Npulses);    %   this is the xcorr (matched filter) array
        lag=zeros(2*(numel(Tpri)+numel(s0)-1)-1,1); %   this is the lag of the xcorr operation

        parfor i=1:Npulses
            yc(:,i)=conv(phic(:,i),s0);

    % 	    if mod(i,20)==0
    % 	        disp(['Convolution ' num2str(i/Npulses*100),'%']);
    % 	    end
        end

        for i=1:Npulses
            [xc_est(:,i),lag]=xcorr(yc(:,i),s0);

    % 	    if mod(i,20)==0
    % 	        disp(['Matched filter ' num2str(i/Npulses*100),'%']);
    % 	    end
        end

        %   remove non-causal lags and those beyond the unambiguous range 
        xc_est(lag<0,:)=[]; 
        lag(lag<0)=[];  
        xc_est(lag>Nsamp,:)=[];
        lag(lag>Nsamp)=[];

        %   Doppler process across slow time to build the range-Doppler map
        RDc=zeros(numel(lag),NFFT);

        for i=1:Nsamp
            RDc(i,:)=fftshift(fft(xc_est(i,:),NFFT));
        end

        RDc=RDc./max(abs(RDc(:)));  %   normalize the clutter range-Doppler map
        fax=linspace(-1,1,NFFT)*PRF/2;  %   frequency axis 

        %   plot the range-Doppler map of the clutter
        figure;
        imagesc(fax*1e-3,0.5*c*Tpri,10*log10(abs(RDc)));
        xlabel('Doppler shift [kHz]');
        ylabel('Ambiguous range [m]');
        title('Clutter');
        set(gca,'ydir','normal'); caxis([-40 0]);
    end

    %%  6. introduce targets into the scene
    if simulate_target_response
        rcst=[params(i_run,TARGET_RCS);params(i_run,TARGET_RCS2)];
        if params(i_run,N_TARGETS) == 0
            rcst = zeros(2,1);
        elseif params(i_run,N_TARGETS) == 1
            rcst(2) = 0;
        end

        elt=[elsq*r2d;elsq*r2d]*d2r;                 %   elevation of the target
        azt=[azsq*r2d;azsq*r2d]*d2r;               %   azimuth of the target
        rtgh=[cos(azt).*cos(elt) sin(azt).*cos(elt) sin(elt)]; %  unit vector in the direction of the target

        rt_unamb=[params(i_run,TARGET_RANGE);params(i_run,TARGET_RANGE2)]*Nmi2m;;
        rt=mod(rt_unamb,0.5*c*PRI);             %   ambiguous range 
        fdt=2*lambda^-1*[params(i_run,TARGET_VEL);params(i_run,TARGET_VEL2)];          %   Doppler for the target

        At=phsq'*exp(-1j*k*Array*rtgh');    %   array factor toward the targets
        powt=abs(At').^2.*rcst.*rt_unamb.^-4;     %   scattered power from the target
        powt=powt./max(powt(:));            %   normalize the scattered power

        pht=exp(1j*2*pi*lambda^-1*(fdt*ones(1,Npulses)).*((PRI*ones(numel(rt),1))...
            *(0:Npulses-1)+rt*ones(1,Npulses)));    %   complex phase of targets across slow time

        [Pulses,rt_grid]=meshgrid(0:Npulses-1,rt);    %   make a grid to insert the targets
        [Pulses_q,rt_grid_q]=meshgrid(0:Npulses-1,0.5*c*Tpri);    %   make the query grid

        pht_grid_Re=interp2(Pulses,rt_grid,powt.*real(pht),Pulses_q,rt_grid_q);
        pht_grid_Im=interp2(Pulses,rt_grid,powt.*imag(pht),Pulses_q,rt_grid_q);
        pht_grid=pht_grid_Re+1j*pht_grid_Im;

        pht_grid(isnan(pht_grid))=0; 

        %%  7. doppler processing 

        %   pre-allocate arrays for speed
        yt=zeros(numel(Tpri)+numel(s0)-1,Npulses);  %   this is the convolution array
        xt_est=zeros(2*(numel(Tpri)+numel(s0)-1)-1,Npulses);    %   this is the xcorr (matched filter) array
        lag=zeros(2*(numel(pht)+numel(s0)-1)-1,1); %   this is the lag of the xcorr operation

        for i=1:Npulses

            temp=conv(pht_grid(:,i),s0);
            yt(:,i)=temp;

    % 	    if mod(i,20)==0
    % 	        disp(['Convolution ' num2str(i/Npulses*100),'%']);
    % 	    end
        end

        for i=1:Npulses
            [xt_est(:,i),lag]=xcorr(yt(:,i),s0);

    % 	    if mod(i,20)==0
    % 	        disp(['Matched filter ' num2str(i/Npulses*100),'%']);
    % 	    end
        end

        %   remove non-causal lags and those beyond the unambiguous range 
        xt_est(lag<0,:)=[];
        lag(lag<0)=[];
        xt_est(lag>Nsamp,:)=[];
        lag(lag>Nsamp)=[];

        %   Doppler process across slow time to build the range-Doppler map
        RDt=zeros(numel(lag),NFFT);

        for i=1:Nsamp
            RDt(i,:)=fftshift(fft(xt_est(i,:),NFFT));
        end

        RDt=RDt./max(abs(RDt(:)));  %   normalize the target range-Doppler map

        %   plot the range-Doppler map of the target
        figure;
        imagesc(fax*1e-3,c*Tpri,10*log10(abs(RDt)));
        xlabel('Doppler shift [kHz]');
        ylabel('Ambiguous range [m]');
        title('Target(s)');
        set(gca,'ydir','normal'); caxis([-40 0]);
        saveas(gcf,['img/tgt_rd_' num2str(i_run,'%04d') '.png']);

    end	
	%   plot the sum of the target and the clutter
% 	figure;
% 	imagesc(fax*1e-3,c*Tpri,10*log10(abs(RDt+RDc)));
% 	xlabel('Doppler shift [kHz]');
% 	ylabel('Ambiguous range [m]');
% 	title('Clutter + target(s)');
% 	set(gca,'ydir','normal'); caxis([-40 0]);
	
	%% 8. Set the clutter, noise, and target signal levels relative to each other
	if ~simulate_clutter_response && simulate_target_response
        RDc = zeros(size(RDt));
    end
	RDn=randn(size(RDc))+1j*randn(size(RDc));
	
	SNR=20; %   max signal to noise level, dB
	CNR=15; %   max clutter to noise level, dB
	SCR=SNR-CNR;
	
	maxPowT=max(abs(RDt(:)));
	maxPowC=max(abs(RDc(:)));
	
	%   normalize to the target, adjust the noise and clutter relative to that
	%   using CNR and SNR
	RDt_final=RDt/maxPowT;
	RDn_final=10^(-0.1*SNR)*RDn;
	RDc_final=10^(-0.1*SCR)*RDc;
	RD=RDt_final+RDc_final+RDn_final;
	
	%   plot 
	figure;
	imagesc(fax*1e-3,c*Tpri,10*log10(abs(RD)));
	xlabel('Doppler shift [kHz]');
	ylabel('Ambiguous range [m]');
	title('Clutter + target(s) + noise');
	set(gca,'ydir','normal'); caxis([-40 0]);
    saveas(gcf,['img/tgt_ct_ns_rd_' num2str(i_run,'%04d') '.png']);
	
	%% 9. CFAR processing
	
	Ntest=10;    %   number of test cells in the range and Doppler dimension (half)
	Nguard=5;   %   number of guard cells in the range and Doppler dimension (half)
	
	Thresh=getCFARLevel(abs(RD),Ntest,Nguard);
	TargetsRegions=abs(RD)>5*Thresh;
	
	%   plot detected target regions 
	figure; hold on;
	imagesc(fax*1e-3,c*Tpri,TargetsRegions);
	xlabel('Doppler shift [kHz]');
	ylabel('Ambiguous range [m]');
	title('Targets detected');
	set(gca,'ydir','normal'); axis('tight');
    saveas(gcf,['img/tgt_det_' num2str(i_run,'%04d') '.png']);
    imwrite(TargetsRegions,colormap(gcf),['img/tgt_det_' num2str(i_run,'%04d') '.jpg'],'jpg');
	
	%   tally number of detections and false detections
	Ndet=sum(TargetsRegions(:));
	
	%   determine the nearest range-Doppler cell to each target
	for i=1:numel(rt)
	    rtcell(i)=min(find(rt(i)<c*Tpri));
	    
	    if fdt(i)>=0
	        fdtcell(i)=min(find(fdt(i)<fax));
	    else
	        fdtcell(i)=max(find(fdt(i)>fax));
	    end
	end
	
	Nfalsedet=Ndet-sum(TargetsRegions(rtcell,fdtcell));
	
% 	plot(fax(fdtcell)*1e-3,0.5*c*Tpri(rtcell),'xr');
% 	legend('True target position'); axis tight;

    %%	Find a single range value for each target
    fig = gcf;
    fig_data = fig.Children.Children;
    fig_x_data = fig_data.XData;
    fig_y_data = fig_data.YData;
    [target_rows,target_cols] = find(TargetsRegions);
    unique_cols = unique(target_cols);
    target_range = zeros(numel(unique_cols),1);
    target_doppler = target_range;
    n_dets_to_keep = 0;
    for i_det = 1:size(unique_cols,1)
        rows_this_target = find(TargetsRegions(:,unique_cols(i_det)));
        if numel(rows_this_target) > 50
            n_dets_to_keep = n_dets_to_keep + 1;
            target_range(n_dets_to_keep) = fig_y_data(floor(median(rows_this_target)));
            target_doppler(n_dets_to_keep) = fig_x_data(unique_cols(i_det));
        else
            TargetsRegions(:,unique_cols(i_det)) = 0;
        end
        
    end
    if n_dets_to_keep < numel(unique_cols)
        target_range(n_dets_to_keep+1:end) = [];
        target_doppler(n_dets_to_keep+1:end) = [];
    end
    
    
	%%
    max_targets_det = max(max_targets_det,n_dets_to_keep);
    
%     results(i_run,N_TARGETS_DETECTED) = 2;
%     results(i_run,RANGE:DOPPLER2) = 0;
%     results(i_run,RUN_TIME) = toc(this_run_timer);
    run_time = toc(this_run_timer);
    result{i_run,1} = target_range;
    result{i_run,2} = target_doppler;
    result{i_run,3} = run_time;
%     return
end

% rad_write_result(logfile,params,results);
disp(['Max targets found: ' num2str(max_targets_det)]);
toc(total_run_timer);
