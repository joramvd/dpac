%% dPAC: A method for debiasing phase-amplitude cross-frequency coupling
% Joram van Driel, Roy Cox & Mike X Cohen
% 2014
% --
% This code accompanies the paper titled "dPAC: A method for debiasing 
% phase-amplitude cross-frequency coupling". Below, simulations are run
% that show how non-sinusoidal oscillations produce phase clustering which
% introduces a bias in Canolty's PAC method.
% Using the code without following the paper may lead to confusion, 
% incorrect data analyses, and misinterpretations of results. 
% The authors assume no responsibility for inappropriate or incorrect use 
% of this code. 


clear all

cd('path\to\folder'); % -- change directory if you want to save output and/or plots

%% parameters

srate = 1000;              % -- sampling rate
t = 1/srate:1/srate:10;    % -- time points: 10000 1ms steps; 10 seconds
ntimepoints = length(t);   % -- get number of timepoints
centTimes = 0:0.2:10;      % -- center times for 'oscillatory' peaks; corresponds to 5 Hz
fpow = 30;                 % -- frequency for power: gamma
nbins = 18;                % -- number of bins used for Tort's MI

gausWidth = linspace(0.001,0.08,50); % -- gaussian width of the 'oscillatory' peaks

[pac dpac mi itpc] = deal(zeros(2,50)); % -- initialize output matrices


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
        
        %% PAC, dPAC, MI (note: dPAC not shown in the paper)
        
        % -- extract phase and power
        thetaphase = angle(hilbert(ts));
        gammapow = abs(gamma);
        
        % -- equations to calculate PAC, dPAC and ITPC (phase clustering vector length)
        pac(in_anti,g)  = abs(mean(exp(1i*thetaphase) .* gammapow));
        debias_term = mean(exp(1i*thetaphase)); % -- this is the phase clustering bias
        dpac(in_anti,g) = abs(mean( (exp(1i*thetaphase) - debias_term) .* gammapow)); % -- which is subtracted here
        itpc(in_anti,g) = abs(mean(exp(1i*thetaphase)));
        
        % Tort's Modulation Index
        
        thetaphase_bin = ceil( tiedrank( thetaphase ) / (ntimepoints / nbins) ); % -- bin the theta phase angles into nbins
        gammapow_bin = zeros(1,nbins);
        for k=1:nbins
            gammapow_bin(k) = squeeze(mean(gammapow(thetaphase_bin==k))); % -- compute mean gamma power in each bin
        end
        gammapow_bin = gammapow_bin ./ sum(gammapow_bin); % -- normalize
        
        mi(in_anti,g) = (log(nbins) + sum(gammapow_bin.*log(gammapow_bin)) ) ./ log(nbins); % -- compute MI
        
        
    end
end

%% figure of three gW examples (Figure 1)

% some annoying code to make nice subplot figure...

plotidx(1,1,:)=[1 2];
plotidx(2,1,:)=[7 8];
plotidx(3,1,:)=[13 14];

plotidx(:,2,:)=plotidx(:,1,:)+3;

plotidx(1,3,:)=[3 6];
plotidx(2,3,:)=[9 12];
plotidx(3,3,:)=[15 18];

figure('position',[500 200 600 400])

gausexamp = [0.01 0.03 0.05]; % -- plot three examples of gaussian 'oscillations'

for q=1:3 % -- loop over the three examples
    
    gW = gausexamp(q);

    % create time series

    cmd = 'ts = detrend(';
    for i=1:length(centTimes)
        cmd = [cmd 'exp(-((t-' num2str(centTimes(i)) ').^2)/(2*' num2str(gW) '^2)) + '];
    end
    cmd = [cmd '0);'];
    eval(cmd)

    %% 

    thetaphase = angle(hilbert(ts));
    
    % -- plot of 'theta' time series
    subplot(6,3,squeeze(plotidx(q,1,:)))
    plot(t,ts,'k','linewidth',1)
    set(gca,'xlim',[0 1.6],'ylim',[-1.5 1.5],'xticklabel',{})
    ylabel('Ampl.')
    box off
    
    % -- polar plot of phase angles, with Phase clustering (PC)
    subplot(6,3,squeeze(plotidx(q,3,:)))
    randsel = randperm(ntimepoints);
    randsel = randsel(1:100);
    h1=polar([thetaphase(randsel); thetaphase(randsel)],[zeros(1,100); ones(1,100)]); hold on
    for hh=1:length(h1)
        hasbehavior(h1(hh),'legend',false);
    end
    set(h1,'color',[0.8 0.8 0.8]);
    txt = findall(gca,'type','text'); delete(txt);
    
    meanvect = mean(exp(1i*thetaphase));
    h=polar([0 angle(meanvect)],[0 abs(meanvect)],'k'); hold on
    set(h,'linewidth',3);
    legend(['PC = ' num2str(abs(meanvect))]);

    % -- phase angle time series
    subplot(6,3,squeeze(plotidx(q,2,:)))
    plot(t,thetaphase,'k','linewidth',1);
    set(gca,'xlim',[0 1.6]);
    ylabel('Phase (rad.)')
    if q==3, xlabel('Time (s)'); end
    box off
    
end


%% figure of two gW examples with coupling (Figure 2)

% some annoying code to make nice subplot figure...
clear plotidx
plotidx(1,1,:)=[1 2];
plotidx(2,1,:)=[4 5];
plotidx(1,2,:)=3;
plotidx(2,2,:)=6;

figure('position',[500 200 600 500])

gausexamp = [0.01 0.05]; % -- plot two examples of gaussian 'oscillations'

for q=1:2 % -- loop over the two examples
    
    gW = gausexamp(q);

    % create time series

    cmd = 'ts = detrend(';
    for i=1:length(centTimes)
        cmd = [cmd 'exp(-((t-' num2str(centTimes(i)) ').^2)/(2*' num2str(gW) '^2)) + '];
    end
    cmd = [cmd '0);'];
    eval(cmd)

    gamma = ((ts+0.5)  .* exp( 1i * 2 * pi * fpow * t ));

    %% PAC

    thetaphase = angle(hilbert(ts));
    gammapow = abs(gamma);
    
    % -- plot of 'theta' time series plus phase-modulated gamma
    subplot(3,3,plotidx(q,1,:))
    plot(t,ts,'k','linewidth',1)
    set(gca,'xlim',[0 0.8],'ylim',[-1.5 1.5])
    hold on
    plot(t,real(gamma),'r','linewidth',0.5);
    ylabel('Amplitude')
    xlabel('Time (s)');
    legend('theta','gamma')
    box off
    
    % -- polar plot of phase angles + PAC
    subplot(3,3,plotidx(q,2,:))
    randsel = randperm(ntimepoints);
    randsel = randsel(1:100);
    h1=polar([thetaphase(randsel); thetaphase(randsel)],[zeros(1,100); ones(1,100)]); hold on
    for hh=1:length(h1)
        hasbehavior(h1(hh),'legend',false);
    end
    set(h1,'color',[0.8 0.8 0.8]);
    txt = findall(gca,'type','text'); delete(txt);
    
    tmppac  = mean(exp(1i*thetaphase) .* gammapow);
    h=polar([0 angle(tmppac)],[0 3*abs(tmppac)],'r'); hold on
    set(h,'linewidth',3);
    legend(['PAC = ' num2str(abs(tmppac))]);
    
end
%%
subplot(337)
plot(gausWidth,itpc,'k') % -- plot PC as a function of the gaussian widths
set(gca,'xlim',[gausWidth(1) gausWidth(end)])
ylabel('PC')
xlabel('Gaussian width')
title('Phase clustering bias')
box off

subplot(338) % -- plot the coupling measures; zscored to have them in the same scale
plot(gausWidth,zscore(pac(1,:)),'b'); hold on
plot(gausWidth,zscore(mi(1,:)),'g');
% plot(gausWidth,zscore(dpac(1,:)),'r'); % -- uncomment to plot dPAC as well (not shown in paper)
ylabel('zscored coupling')
xlabel('Gaussian width')
title('Coupling-clustering anti-phase')
set(gca,'xlim',[gausWidth(1) gausWidth(end)],'ylim',[-2 2.5])
box off

subplot(339) % -- the same but now with clustering and coupling in-phase
plot(gausWidth,zscore(pac(2,:)),'b'); hold on
plot(gausWidth,zscore(mi(2,:)),'g');
% plot(gausWidth,zscore(dpac(2,:)),'r'); % -- uncomment to plot dPAC as well (not shown in paper)
set(gca,'xlim',[gausWidth(1) gausWidth(end)],'ylim',[-2 2.5])
ylabel('zscored coupling')
xlabel('Gaussian width')
title('Coupling-clustering in-phase')
legend('pac','mi')
% legend('pac','mi','dpac')
legend boxoff
box off




