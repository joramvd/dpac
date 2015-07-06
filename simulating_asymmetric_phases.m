%% dPAC: A method for debiasing phase-amplitude cross-frequency coupling
% Joram van Driel, Roy Cox & Mike X Cohen
% 2014/2015
% --
% This code accompanies the paper titled "Phase clustering bias in
% phase-amplitude cross-frequency coupling and its removal". Below, 
% simulations are run that show how non-sinusoidal properties can produce
% phase clustering bias in different coupling measures.
% Using the code without following the paper may lead to confusion, 
% incorrect data analyses, and misinterpretations of results. 
% The authors assume no responsibility for inappropriate or incorrect use 
% of this code. 


clear all

cd('Z:\PhD\dPAC'); % -- change directory if you want to save output and/or plots

%% parameters

srate = 1000;              % -- sampling rate
t = 1/srate:1/srate:10;    % -- time points: 10000 1ms steps; 10 seconds
ntimepoints = length(t);   % -- get number of timepoints
centTimes = 0:0.2:10;      % -- center times for 'oscillatory' peaks; corresponds to 5 Hz
fpow = 30;                 % -- frequency for power: gamma
nbins = 18;                % -- number of bins used for Tort's MI

gausWidth = linspace(0.001,0.08,50); % -- gaussian width of the 'oscillatory' peaks

[pac, dpac, mi, pc, plv] = deal(zeros(2,50)); % -- initialize output matrices


%%
q=0;

for in_anti = 1:2 % -- simulate two scenarios: coupling angle in/anti-phase with theta peak
    for g = 1:50 % -- loop over 50 widths of the gaussian 'cycles'
        
        gW = gausWidth(g); % -- set the width of the gaussian
        
        % -- create time series
        % -- this codes concatenates a series of gaussians that are
        % -- then detrended to 'oscillate' around 0
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
        
        %% PAC, dPAC, MI, PLV (note: dPAC not shown in the paper)
        
        % -- extract phase and power
        thetaphase = angle(hilbert(ts));
        gammapow = abs(gamma);
        
        % -- equations to calculate PAC, dPAC and PC (phase clustering vector length)
        pac(in_anti,g)  = abs(mean(exp(1i*thetaphase) .* gammapow));
        debias_term = mean(exp(1i*thetaphase)); % -- this is the phase clustering bias
        dpac(in_anti,g) = abs(mean( (exp(1i*thetaphase) - debias_term) .* gammapow)); % -- which is subtracted here
        pc(in_anti,g) = abs(mean(exp(1i*thetaphase)));
        
        % Tort's Modulation Index (Tort et al., 2010)
        
        thetaphase_bin = ceil( tiedrank( thetaphase ) / (ntimepoints / nbins) ); % -- bin the theta phase angles into nbins
        gammapow_bin = zeros(1,nbins);
        for k=1:nbins
            gammapow_bin(k) = squeeze(mean(gammapow(thetaphase_bin==k))); % -- compute mean gamma power in each bin
        end
        gammapow_bin = gammapow_bin ./ sum(gammapow_bin); % -- normalize
        
        mi(in_anti,g) = (log(nbins) + sum(gammapow_bin.*log(gammapow_bin)) ) ./ log(nbins); % -- compute MI
       
        
        % -- Phase Locking Value (Cohen, 2008; Colgin et al 2009)
        plv(in_anti,g) = abs(mean(exp(1i*( thetaphase - angle(hilbert(detrend(gammapow))) ))));
        
        
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
plotidx(1,1,:)=1:2;
plotidx(1,2,:)=[5 5];
plotidx(1,3,:)=[6 6];
plotidx(2,1,:)=3:4;
plotidx(2,2,:)=[7 7];
plotidx(2,3,:)=[8 8];

figure('position',[500 200 900 300])

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
    subplot(2,4,plotidx(q,1,:))
    plot(t,ts,'k','linewidth',1)
    set(gca,'xlim',[0 0.8],'ylim',[-1.5 1.5])
    hold on
    plot(t,real(gamma),'r','linewidth',0.5);
    ylabel('Amplitude')
    xlabel('Time (s)');
    legend('theta','gamma')
    box off
    
    % -- polar plot of phase angles + PAC
    subplot(2,4,plotidx(q,2,:))
    randsel = randperm(ntimepoints);
    randsel = randsel(1:100);
    h1=polar([thetaphase(randsel); thetaphase(randsel)],[zeros(1,100); ones(1,100)]); hold on
    for hh=1:length(h1)
        hasbehavior(h1(hh),'legend',false);
    end
    set(h1,'color',[0.8 0.8 0.8]);
    txt = findall(gca,'type','text'); delete(txt);
    
    % compute PAC and plot as average vector in polar space
    tmppac  = mean(exp(1i*thetaphase) .* gammapow);
    h=polar([0 angle(tmppac)],[0 3*abs(tmppac)],'r'); hold on
    set(h,'linewidth',3);
    legend(['PAC = ' num2str(abs(tmppac))]);

    % compute MI
    thetaphase_bin = ceil( tiedrank( thetaphase ) / (ntimepoints / nbins) ); % -- bin the theta phase angles into nbins
    gammapow_bin = zeros(1,nbins);
    for k=1:nbins
        gammapow_bin(k) = squeeze(mean(gammapow(thetaphase_bin==k))); % -- compute mean gamma power in each bin
    end
    gammapow_bin = gammapow_bin ./ sum(gammapow_bin); % -- normalize
    tmpmi = (log(nbins) + sum(gammapow_bin.*log(gammapow_bin)) ) ./ log(nbins);
    
    % plot the binned gammapower as a function of bin #, as a histogram
    subplot(2,4,plotidx(q,3,:)); 
    bar(1:18,gammapow_bin,'style','hist')
    set(gca,'ylim',[0 0.2])
    set(gca,'xlim',[0 19])
    text(2,0.18,['MI = ' num2str(tmpmi)]);
    box off
    
end


%% plot of different CFC measures (figure 3A)

figure('position',[400 100 400 400])

subplot(221)
plot(gausWidth,pc(1,:),'b'); hold on
plot(gausWidth,pc(2,:),'--b','linewidth',1.5);
set(gca,'xscale','log','xtick',[0.001 0.004 0.016 0.064])
title('PC')
box off
legend('anti-phase','in-phase')

subplot(222)
plot(gausWidth,pac(1,:),'r'); hold on
plot(gausWidth,pac(2,:),'--r','linewidth',1.5); 
set(gca,'xscale','log','xtick',[0.001 0.004 0.016 0.064])
title('(d)PAC')
box off

subplot(223)
plot(gausWidth,mi(1,:),'r'); hold on
plot(gausWidth,mi(2,:),'--r','linewidth',1.5);
set(gca,'xscale','log','xtick',[0.001 0.004 0.016 0.064])
title('MI')
box off

subplot(224)
plot(gausWidth,plv(1,:),'r'); hold on
plot(gausWidth,plv(2,:),'--r','linewidth',1.5);
set(gca,'xscale','log','xtick',[0.001 0.004 0.016 0.064])
title('Coh.')
box off














%% Additional plot to illustrate PLV behavior to asymmetric oscillations
%  Polar plots of phase angle differences between theta phase and phase of
%  gamma power envelope; note how a narrow guass-width of the asymmetric
%  'cycles' with clustering in-phase with angle of CFC results in reduced
%  PLV (top-right subplot)
%  This plot is not shown in the paper, but results are discussed

figure('position',[500 200 400 400])

gausexamp = [0.01 0.05]; % -- plot two examples of gaussian 'oscillations'
k=1;
for q=1:2 % -- loop over the two examples
    
    for in_anti=1:2
        gW = gausexamp(q);
        
        % create time series
        
        cmd = 'ts = detrend(';
        for i=1:length(centTimes)
            cmd = [cmd 'exp(-((t-' num2str(centTimes(i)) ').^2)/(2*' num2str(gW) '^2)) + '];
        end
        cmd = [cmd '0);'];
        eval(cmd)
        
        gamma = ((ts+0.5)  .* exp( 1i * 2 * pi * fpow * t ));
        
        if in_anti==2
            ts = ts([101:end 1:100]); % -- shift the time series with pi (100 ms)
        end
        
        % -- extract phase and power
        thetaphase = angle(hilbert(ts));
        gammapow = abs(gamma);
        gammapowphase = angle(hilbert(detrend(gammapow)));
        
        % -- Phase Locking Value (Cohen, 2008; Colgin et al 2009)
        
        phasediffs = (thetaphase-mean(thetaphase)) - (gammapowphase - mean(gammapowphase));
        plvtmp = abs(mean(exp(1i*( phasediffs ))));
        
        subplot(2,2,k)        
        randsel = randperm(ntimepoints);
        randsel = randsel(1:100);
        h1=polar([phasediffs(randsel); phasediffs(randsel)],[zeros(1,100); ones(1,100)]); hold on
        for hh=1:length(h1)
            hasbehavior(h1(hh),'legend',false);
        end
        set(h1,'color','r');
        if in_anti==1, inantitle='anti-phase'; else inantitle='in-phase'; end;
        title(['g = ' num2str(gW) ' - ' inantitle])

        k=k+1;
    end
    
end

           
