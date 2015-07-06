%% dPAC: A method for debiasing phase-amplitude cross-frequency coupling
% Joram van Driel, Roy Cox & Mike X Cohen
% 2014/2015
% --
% This code accompanies the paper titled "Phase clustering bias in
% phase-amplitude cross-frequency coupling and its removal". Below, 
% simulations are run that show how non-sinusoidal properties can produce
% phase clustering bias in different coupling measures.
% In addition, band-pass filtering is used to approximate regular time
% series analysis approaches to real data.
% Using the code without following the paper may lead to confusion, 
% incorrect data analyses, and misinterpretations of results. 
% The authors assume no responsibility for inappropriate or incorrect use 
% of this code. 


clear all

cd('Z:\PhD\dPAC'); % -- change directory if you want to save output and/or plots

%% parameters

srate = 1000;                  % -- sampling rate
t = 1/srate:1/srate:12;        % -- time points: 10000 1ms steps; 12 seconds (10 seconds plus padded 1000 ms on each side for edge artifacts)
ntimepoints = length(t)-2000;  % -- get number of timepoints minus padded buffer zones
centTimes = 0:0.2:12;          % -- center times for 'oscillatory' peaks; corresponds to 5 Hz
fpow = 30;                     % -- frequency for power: gamma
nbins = 18;                    % -- number of bins used for Tort's MI

gausWidth = linspace(0.001,0.08,50); % -- gaussian width of the 'oscillatory' peaks

[pac, dpac, mi, pc, plv] = deal(zeros(2,50)); % -- initialize output matrices


%%
q=0;

for in_anti = 1:2 % -- simulate two scenarios: phase clustering bias as in vs. anti phase with coupling angle
    for g = 1:50 % -- loop over 50 widths of the gaussian 'cycles'
        
        gW = gausWidth(g); % -- set the width of the gaussian
        
        % -- create time series
        % -- this codes concatenates a series of gaussians that are
        % -- then detrended
        cmd = 'ts = detrend(';
        for i=1:length(centTimes)
            cmd = [cmd 'exp(-((t-' num2str(centTimes(i)) ').^2)/(2*' num2str(gW) '^2)) + '];
        end
        cmd = [cmd '0);'];
        eval(cmd)
        
        % -- gamma is a complex sine wave that is phase-modulated by the
        % -- above time series
        gamma = ((ts+0.5)  .* exp( 1i * 2 * pi * fpow * t ));
        
        if in_anti==2
            ts = ts([101:end 1:100]); % -- shift the time series with pi (100 ms)
        end
        
        thetagamma = real(gamma) + ts;
        
        %% basic filter settings
        
        % -- band-pass filter
        thetaband           = [3 7];
        gammaband           = [25 35];
        
        theta_filt_order    = round(3*(srate/thetaband(1)));
        theta_filterweights = fir1(theta_filt_order,[mean(thetaband)-(thetaband(2)-thetaband(1)) mean(thetaband)+(thetaband(2)-thetaband(1))]/(srate/2));
        
        gamma_filt_order    = round(3*(srate/gammaband(1)));
        gamma_filterweights = fir1(gamma_filt_order,[mean(gammaband)-(gammaband(2)-gammaband(1)) mean(gammaband)+(gammaband(2)-gammaband(1))]/(srate/2));
        
        thetafilt = filtfilt(theta_filterweights,1,thetagamma); thetafilt = thetafilt(1001:end-1000);
        gammafilt = filtfilt(gamma_filterweights,1,thetagamma); gammafilt = gammafilt(1001:end-1000);

        %% PAC, dPAC, MI, PLV (note: dPAC not shown in the paper)
        
        % -- extract phase and power
        thetaphase = angle(hilbert(thetafilt));
        gammapow = abs(hilbert(gammafilt));
        
        % -- equations to calculate PAC, dPAC and PC (phase clustering vector length)
        pac(in_anti,g)  = abs(mean(exp(1i*thetaphase) .* gammapow));
        debias_term = mean(exp(1i*thetaphase)); % -- this is the phase clustering bias
        dpac(in_anti,g) = abs(mean( (exp(1i*thetaphase) - debias_term) .* gammapow)); % -- which is subtracted here
        pc(in_anti,g) = abs(mean(exp(1i*thetaphase)));
        
        % -- Tort's Modulation Index
        
        thetaphase_bin = ceil( tiedrank( thetaphase ) / (ntimepoints / nbins) ); % -- bin the theta phase angles into nbins
        gammapow_bin = zeros(1,nbins);
        for k=1:nbins
            gammapow_bin(k) = squeeze(mean(gammapow(thetaphase_bin==k))); % -- compute mean gamma power in each bin
        end
        gammapow_bin = gammapow_bin ./ sum(gammapow_bin); % -- normalize
        
        mi(in_anti,g) = (log(nbins) + sum(gammapow_bin.*log(gammapow_bin)) ) ./ log(nbins); % -- compute MI
        
        % -- Phase-locking value (Cohen, 2008; Colgin et al 2009)        

        plv(in_anti,g) = abs(mean(exp(1i*(thetaphase- angle(hilbert(detrend(gammapow)))))));
        
                
    end
end


%% plot of different CFC measures

figure('position',[400 100 400 600])

subplot(321)
plot(gausWidth,pc(1,:),'b','linewidth',1); hold on
plot(gausWidth,pc(2,:),':b','linewidth',2);
set(gca,'xscale','log','xtick',[0.001 0.004 0.016 0.064],'ylim',[0 0.8]) %
title('PC')
box off

% -- small inlay plot to zoom in on small bump of PC increase for PC being
% -- anti-phase with CFC angle (not visible because of strong PC increase
% -- for low-g and PC-CFC in-phase
ax1 = gca;
ax1_pos = get(ax1,'position');
ax2 = axes('Position',[0.27 0.78 0.15 0.1],...
    'XAxisLocation','bottom',...
    'YAxisLocation','left',...
    'Color','none');
plot(gausWidth,pc(1,:),'parent',ax2,'color','b','linewidth',1); hold on
set(gca,'xscale','log','xtick',[0.001 0.064],'ylim',[0 0.006])
box off

subplot(322)
plot(gausWidth,pac(1,:),'r','linewidth',1); hold on
plot(gausWidth,pac(2,:),':r','linewidth',2); 
set(gca,'xscale','log','xtick',[0.001 0.004 0.016 0.064])
title('(d)PAC')
box off
legend('anti-phase','in-phase')

subplot(323)
plot(gausWidth,mi(1,:),'r','linewidth',1); hold on
plot(gausWidth,mi(2,:),':r','linewidth',2);
set(gca,'xscale','log','xtick',[0.001 0.004 0.016 0.064],'ylim',[0 0.1])
title('MI')
box off

subplot(324)
plot(gausWidth,plv(1,:),'r','linewidth',1); hold on
plot(gausWidth,plv(2,:),':r','linewidth',2);
set(gca,'xscale','log','xtick',[0.001 0.004 0.016 0.064])
title('PLV')
box off

