%% dPAC: A method for debiasing phase-amplitude cross-frequency coupling
% Joram van Driel, Roy Cox & Mike X Cohen
% 2014/2015
% --
% This code accompanies the paper titled "Phase clustering bias in
% phase-amplitude cross-frequency coupling and its removal". Below, 
% simulations are run to test four phase-amplitude cross-frequency 
% copuling measures (PAC, dPAC MI and PLV) as a function of phase 
% clustering bias and coupling angle.
% Using the code without following the paper may lead to confusion, 
% incorrect data analyses, and misinterpretations of results. 
% The authors assume no responsibility for inappropriate or incorrect use 
% of this code. 

clear, close all

cd('Z:\PhD\dPAC'); % -- change directory if you want to save output and/or plots

%% Ingredients

save_results = true;
coupling = true;    % -- set to false to run simulation without coupling (to test for false positives)
permutation = true; % -- run a permutation test? slows down simulation!
nperm = 1000;       % -- number of permutations; set to lower number to speed up analysis
            
srate = 1000;              % -- sampling rate
t = 1/srate:1/srate:10;    % -- time points: 10 seconds

fpow = 30;                 % -- frequency for power: gamma
fphase = 5;                % -- frequency for phase: theta
ntimepoints = length(t);   % -- get number of timepoints

biasincrease = logspace(log10(5),log10(0.5)); % -- parameter used for different widths of the theta phase distribution, to simulate clustering

%% Initialize variables

num_time_shifts = 50; % -- angle-shift loop-resolution; set this to a lower value to speed up the analysis
time_shifts = round(linspace(1,51,num_time_shifts));

num_bias_levels = 51; % -- bias-level loop-resolution; set this to a lower value to speed up the analysis
bias_levels = round(linspace(1,50,num_bias_levels));

 % -- if the simulation is run with non-coupled signals, temporally shifting theta will have no effect, so we can skip this part
if ~coupling, 
    num_time_shifts = 1;
    time_shifts = 0; 
end

[PAC,dPAC,MI,PLV, PACz,dPACz,MIz,PLVz] = deal( zeros(num_bias_levels,num_time_shifts) );  % -- initialize output matrices: 50 clustering levels and 51 PAC angles  

%% Simulation

% -- create complex sine waves
theta = 5.*(exp( 1i * 2 * pi * fphase * t ));
if coupling
    gamma = ((real(theta)+6)  .* exp( 1i * 2 * pi * fpow * t )); % -- gamma is phase-modulated by theta; i.e. pure cross-frequency coupling
else
    gamma = exp( 1i * 2 * pi * fpow * t ); % -- pure gamma sine wave; i.e. no coupling
end

% -- compute theta phase angles and gamma power
thetaphase = angle(theta);
gammapower = abs(gamma);

% -- compute PAC and dPAC without bias
PAC_nobias = abs(mean(exp(1i*thetaphase) .* gammapower)); % -- this is the regular PAC equation (Canolty et al., 2006)

debias_term = mean(exp(1i*thetaphase)); % -- this is the phase clustering bias
dPAC_nobias = abs(mean( (exp(1i*thetaphase) - debias_term) .* gammapower)); % -- which is removed from every phase angle: dPAC

% -- note that both show the same value (because the debias_term is zero)
% -- this value will also not change as a function of coupling angle

%% -- compute MI (Tort et al., 2010)

nbins = 18;
thetaphase_bin = ceil( tiedrank( thetaphase ) / (ntimepoints / nbins) );

gammapower_bin = zeros(1,nbins);
for k=1:nbins
    gammapower_bin(k) = squeeze(mean(gammapower(thetaphase_bin==k)));
end
gammapower_bin = gammapower_bin ./ sum(gammapower_bin);

MI_nobias = (log(nbins) + sum(gammapower_bin.*log(gammapower_bin)) ) ./ log(nbins);

%% -- Phase-locking value (Cohen, 2008; Colgin et al 2009)

gammapower_phase = angle(hilbert(detrend(gammapower))); % -- note: hilbert transform of power envelope of complex signal
PLV_nobias = abs(mean(exp(1i*(thetaphase-gammapower_phase)))); % -- gives coherence of 1.0


%%
% -- empty strings to report progress of loop
msg='';
reverseStr = '';

% -- first loop: over angle difference between PAC and theta phase bias

tii = 0; % -- initialize loop-counter variables
for ti = time_shifts % -- ti is the index of the temporal shift, to change the angle of maximal coupling
    tii = tii+1; % -- increment counter
    
    if ti>1
        theta_shift = theta([ 2*(ti-1):end 1:2*(ti-1)-1]); % -- time-shift theta; step size depends on num_time_shift variable
        thetaphase = angle(theta_shift); % -- recompute theta phase angles
    end
    
    % -- second loop: over different levels of phase clustering bias
    
    nii = 0; % -- initialize loop-counter variables
    for ni = bias_levels % -- ni is the index of degree of phase clustering bias (bias will increase with increasing ni)
        nii = nii+1; % -- increment counter
        
        % -- display progress
        msg = sprintf('Running gauss-level %i/%i anglediff %i/%i...',  nii,num_bias_levels,tii,num_time_shifts);
        fprintf([reverseStr, msg]);
        reverseStr = repmat(sprintf('\b'), 1, length(msg));
        
        % -- create G-shaped phase distribution
        gaussiandist = randn(1,ntimepoints*2)*biasincrease(ni); % -- gaussian distribution with a particular width; width decreases (clustering increases) over loop iterations
        pidist = gaussiandist;
        pidist(pidist<-pi | pidist>pi) = []; % -- cut the distribution off at -pi and pi
        
        % -- use gauss-based angle-distribution to sample from theta phase
        % -- distribution
        clear idx
        for k=1:length(pidist)
            [~,idx(k)]=min(abs(pidist(k)-thetaphase));
        end
        idx=sort(idx); % -- note that this gives an extreme case of 
        % non-uniform phase angles; re-plotting the time series shows 
        % highly distorted time series; however, using the 'idx' variable,
        % the corresponding gammapower values are retained along with the
        % (unnatural) selection of thetaphase values, so (d)PAC can still
        % be calculated

        % -- in below variables, G stands for 'gaussian'
        thetaphaseG = thetaphase(idx); % -- theta phase angles now have clustering around 0pi
        gammapowerG = gammapower(idx); % -- save the corresponding gamma power values
        ntimepointsG = length(idx);
        
        % -- compute PAC and dPAC, now with phase clustering bias; note
        % -- that the debias term is now incorporated in dPAC into one line of
        % -- code
        PAC(nii,tii)  = abs(mean(exp(1i*thetaphaseG) .* gammapowerG));
        dPAC(nii,tii) = abs(mean( (exp(1i*thetaphaseG) - mean(exp(1i*thetaphaseG))) .* gammapowerG));
                
        %% -- compute MI (Tort et al., 2010)
        
        thetaphaseG_bin = ceil( tiedrank( thetaphaseG ) / (ntimepointsG / nbins) );
        
        gammapowerG_bin = zeros(1,nbins);
        for k=1:nbins
            gammapowerG_bin(k) = squeeze(mean(gammapowerG(thetaphaseG_bin==k)));
        end
        gammapowerG_bin = gammapowerG_bin ./ sum(gammapowerG_bin);
        
        MI(nii,tii) = (log(nbins) + sum(gammapowerG_bin.*log(gammapowerG_bin)) ) ./ log(nbins);

        %% -- compute PLV
        
        gammapowerG_phase = angle(hilbert(detrend(gammapowerG))); % -- note: hilbert transform of power envelope of complex signal
        PLV(nii,tii) = abs(mean(exp(1i*(thetaphaseG-gammapowerG_phase)))); 
        
        %% now with permutation testing
        
        if permutation

            [fake_PAC,fake_dPAC,fake_MI,fake_PLV,fake_PAC_diff,fake_dPAC_diff,fake_MI_diff,fake_PLV_diff] = deal(zeros(1,nperm)); % -- initialize matrices that will store the null-disbribution of surrogate data

            for permi = 1:nperm

                % -- cut-and-paste a random portion of the data; this preserves
                % -- temporal autocorrelation while removing the coupling
                cutLoc = 5 + randperm(ntimepointsG-10); % -- 5 and 10 prevent the first and last time points from being selected
                cutLoc = cutLoc(1);
                thetaphaseG_shuffled = thetaphaseG([cutLoc:end 1:cutLoc-1]);
                
                fake_PAC(permi)  = abs(mean(exp(1i*thetaphaseG_shuffled) .* gammapowerG)); % -- compute surrogate PAC
                fake_dPAC(permi) = abs(mean((exp(1i*thetaphaseG_shuffled) - mean(exp(1i*thetaphaseG_shuffled))) .* gammapowerG)); % -- compute surrogate dPAC
                fake_PLV(permi)  = abs(mean(exp(1i*(thetaphaseG_shuffled-gammapowerG_phase)))); % -- compute surrogate PLV
                
                % -- compute MI (Tort et al., 2010)
                thetaphaseG_bin_shuffled = ceil( tiedrank( thetaphaseG_shuffled ) / (ntimepointsG / nbins) );
                gammapowerG_bin = zeros(1,nbins);
                for k=1:nbins
                    gammapowerG_bin(k) = squeeze(mean(gammapowerG(thetaphaseG_bin_shuffled==k)));
                end
                gammapowerG_bin = gammapowerG_bin ./ sum(gammapowerG_bin);
                
                fake_MI(permi) = (log(nbins) + sum(gammapowerG_bin.*log(gammapowerG_bin)) ) ./ log(nbins);
            end

            % -- below, the zscore is defined by observed value minus mean of
            % -- null-distribution, divided by the standard deviation of
            % -- null-distribution
            PACz(nii,tii)  = (squeeze(PAC(nii,tii))' - squeeze(mean(fake_PAC))) ./ squeeze(std(fake_PAC)); 
            dPACz(nii,tii) = (squeeze(dPAC(nii,tii))' - squeeze(mean(fake_dPAC))) ./ squeeze(std(fake_dPAC)); 
            MIz(nii,tii)   = (squeeze(MI(nii,tii))' - squeeze(mean(fake_MI))) ./ squeeze(std(fake_MI));
            PLVz(nii,tii)  = (squeeze(PLV(nii,tii))' - squeeze(mean(fake_PLV))) ./ squeeze(std(fake_PLV));
            
        end
    end
    
end
fprintf(' done\n');

%%
if save_results
    % -- save the results if you want, in an appropriate directory
    if coupling
        filename = 'dPAC_MI_PLV_simresults_biasincrease.mat';
    else
        filename = 'dPAC_MI_PLV_simresults_biasincrease_NOcoupling.mat';
    end
    save(filename,....
        'PAC_nobias', 'dPAC_nobias', 'MI_nobias','PLV_nobias',...
        'PAC', 'dPAC', 'MI', 'PLV', 'PACz', 'dPACz', 'MIz', 'PLVz',...
        'num_bias_levels', 'num_time_shifts');
end

%% Plot the result (Figure 4) -- only when coupling set to true

if coupling,
    
    % -- in the paper, a particular colormap and colorscaling is used that
    % -- necessitates some third-party functions that are provided with the
    % -- code;
    % -- set to false to use default jet map and symmetric color scaling
    
    plot_like_paper = true;
    
    % -- when set to true, the othercolor function is used, which can be downloaded here: 
    % -- [http://nl.mathworks.com/matlabcentral/fileexchange/30564-othercolor]
    % -- the diverging_map function can be downloaded here:
    % -- [http://www.sandia.gov/~kmorel/documents/ColorMaps/diverging_map.m]
    
    %%
    angleaxis = 0:pi/(num_time_shifts-1):pi; % -- the x-axis will show coupling angles in fractions of pi
    
    figure('position',[600 300 1200 400]); % -- these values may need to be changed depending on screen settings
    
    % -- plot PAC; from every angle-by-phase-clustering point, the "pure" PAC
    % -- (no bias and 0pi coupling) is subtracted, and %change is computed
    subplot(241);
    PAC_perc = 100.*((PAC - PAC_nobias)./PAC_nobias);
    contourf(angleaxis,1:num_bias_levels,PAC_perc,40,'linestyle','none'); % -- here, contourf is used for smooth contours; alternatively, you can use imagesc (which requires flipping the y-axis direction)
    cl = max(abs(get(gca,'clim')));
    set(gca,'clim',[-cl cl],...
        'ytick',[1 num_bias_levels],'yticklabel',{'min','max'},...
        'xtick',0:pi/2:pi,'xticklabel',{'0', '1/2', '1'}); colorbar
    ylabel('PC bias')
    xlabel('Coupling angle (\pi)')
    title('PAC')
    
    PAC_thresh = zeros(size(PAC_perc));
    PAC_thresh(PAC_perc>2*std(PAC_perc(:)) | PAC_perc<-2*std(PAC_perc(:)))=1;
    hold on
    contour(angleaxis,1:num_bias_levels,PAC_thresh,1,'k','LineWidth',1) % -- plot 2SD as a qualitative threshold of strength of deviation from true coupling
    
    % -- the same now for dPAC
    subplot(242);
    dPAC_perc = 100.*((dPAC - dPAC_nobias)./dPAC_nobias);
    contourf(angleaxis,1:num_bias_levels,dPAC_perc,40,'linestyle','none');
    cl = max(abs(get(gca,'clim')));
    set(gca,'clim',[-cl cl],...
        'ytick',[1 num_bias_levels],'yticklabel',{},...
        'xtick',0:pi/2:pi,'xticklabel',{'0', '1/2', '1'}); colorbar
    title('dPAC')
    
    dPAC_thresh = zeros(size(dPAC_perc));
    dPAC_thresh(dPAC_perc>2*std(dPAC_perc(:)) | dPAC_perc<-2*std(dPAC_perc(:)))=1;
    hold on
    contour(angleaxis,1:num_bias_levels,dPAC_thresh,1,'k','LineWidth',1) % -- plot p<0.001 clusters as overlaid black line
    
    % -- the same now for MI
    MI_perc = 100.*((MI - MI_nobias)./MI_nobias);
    subplot(243);
    contourf(angleaxis,1:num_bias_levels,MI_perc,40,'linestyle','none');
    cl = max(abs(get(gca,'clim')));
    set(gca,'clim',[-cl cl],...
        'ytick',[1 num_bias_levels],'yticklabel',{ },...
        'xtick',0:pi/2:pi,'xticklabel',{'0', '1/2', '1'}); colorbar
    title('MI')
    
    MI_thresh = zeros(size(MI_perc));
    MI_thresh(MI_perc>2*std(MI_perc(:)) | MI_perc<-2*std(MI_perc(:)))=1;
    hold on
    contour(angleaxis,1:num_bias_levels,MI_thresh,1,'k','LineWidth',1) % -- plot p<0.001 clusters as overlaid black line

    % -- the same now for PLV
    PLV_perc = 100.*((PLV - PLV_nobias)./PLV_nobias);
    subplot(244);
    contourf(angleaxis,1:num_bias_levels,PLV_perc,40,'linestyle','none');
    cl = max(abs(get(gca,'clim')));
    set(gca,'clim',[-cl cl],...
        'ytick',[1 num_bias_levels],'yticklabel',{ },...
        'xtick',0:pi/2:pi,'xticklabel',{'0', '1/2', '1'}); colorbar
    title('PLV')
    
    PLV_thresh = zeros(size(PLV_perc));
    PLV_thresh(PLV_perc>2*std(PLV_perc(:)) | PLV_perc<-2*std(PLV_perc(:)))=1;
    hold on
    contour(angleaxis,1:num_bias_levels,PLV_thresh,1,'k','LineWidth',1) % -- plot p<0.001 clusters as overlaid black line

    if plot_like_paper,
        addpath(genpath('Z:\PhD\dPAC\plot_like_paper')); % -- change path accordingly
        colormap(othercolor('BuDRd_18'))
    end
    
%%
    % -- plot PACz
    if plot_like_paper, figure('position',[600 300 1200 400]); end % -- new figure needed to have different colormaps
    subplot(245); % -- PACz
    contourf(angleaxis,1:num_bias_levels,squeeze(PACz),40,'linestyle','none');
    cl = max(abs(get(gca,'clim')));
    set(gca,'clim',[-cl cl],...
        'ytick',[1 num_bias_levels],'yticklabel',{'min','max'},...
        'xtick',0:pi/2:pi,'xticklabel',{'0', '1/2', '1'}); colorbar
    ylabel('PC bias')
    title('PACz')
    
    threshold_001 = icdf('normal',1-0.001/2,0,1); % -- z-value corresponding to 99.9% or p = 0.001 of z-distribution
    threshold_05 = icdf('normal',1-0.05/2,0,1); % -- z-value corresponding to 95% or p = 0.05 of z-distribution

    PACz_thres_001 = squeeze(PACz);
    PACz_thres_001(PACz_thres_001<threshold_001 & PACz_thres_001>(-1*threshold_001))=0; % -- binary threshold: set all values with p>0.001 to zero
    PACz_thres_001(PACz_thres_001~=0)=1; % -- everything else to 1
    PACz_thres_05 = squeeze(PACz);
    PACz_thres_05(PACz_thres_05<threshold_05 & PACz_thres_05>(-1*threshold_05))=0; % -- binary threshold: set all values with p>0.05 to zero
    PACz_thres_05(PACz_thres_05~=0)=1; % -- everything else to 1

    % -- cluster size thresholding (image processing toolbox required)
    if exist('bwlabel','file')
        [r,c]=bwlabel(PACz_thres_001,4); % -- function that searches for 2D clusters
        for ci=1:c
            if sum(any(r==ci,1))<5 || sum(any(r==ci,2))<5 % -- arbitrary constraint: at least 5 contiguous points
                PACz_thres_001(r==ci)=0;
            end
        end
        hold on
        contour(angleaxis,1:num_bias_levels,PACz_thres_001,1,'k','LineWidth',1) % -- plot p<0.001 clusters as overlaid black line
        [r,c]=bwlabel(PACz_thres_05,4); % -- function that searches for 2D clusters
        for ci=1:c
            if sum(any(r==ci,1))<5 || sum(any(r==ci,2))<5 % -- arbitrary constraint: at least 5 contiguous points
                PACz_thres_05(r==ci)=0;
            end
        end
        hold on
        contour(angleaxis,1:num_bias_levels,PACz_thres_05,1,'r','LineWidth',1) % -- plot p<0.05 clusters as overlaid red line
    end
    
    % -- plot dPACz
    subplot(246); 
    contourf(angleaxis,1:num_bias_levels,dPACz,40,'linestyle','none');
    cl = max(abs(get(gca,'clim')));
    set(gca,'clim',[-cl cl],...
        'ytick',[1 num_bias_levels],'yticklabel',{},...
        'xtick',0:pi/2:pi,'xticklabel',{'0', '1/2', '1'}); colorbar
    xlabel('Clustering angle (\pi)')
    title('dPACz')
    
    % -- plot MIz
    subplot(247); 
    contourf(angleaxis,1:num_bias_levels,MIz,40,'linestyle','none');
    cl = max(abs(get(gca,'clim')));
    set(gca,'clim',[-cl cl],...
        'ytick',[1 num_bias_levels],'yticklabel',{},...
        'xtick',0:pi/2:pi,'xticklabel',{'0', '1/2', '1'}); colorbar
    title('MIz')
    
    % -- plot PLVz
    subplot(248); 
    contourf(angleaxis,1:num_bias_levels,PLVz,40,'linestyle','none');
    cl = max(abs(get(gca,'clim')));
    set(gca,'clim',[-cl cl],...
        'ytick',[1 num_bias_levels],'yticklabel',{},...
        'xtick',0:pi/2:pi,'xticklabel',{'0', '1/2', '1'}); colorbar
    title('PLVz')

    if plot_like_paper,
        oldmap = colormap;
        newmap = diverging_map(0:1/length(oldmap):1,[0 0 1],[0 1 0]);
        colormap(newmap)
    end
    
    %% Example plot of different angle distributions + coupling angles (Figure 3)

    figure('position',[400 100 500 500])
    bias2plot = [5 2.5 1.25 0.625];
    for N=1:4

        % -- create G-shaped phase distribution
        gaussiandist = randn(1,ntimepoints*2)*bias2plot(N);
        pidist = gaussiandist;
        pidist(pidist>pi)=[];
        pidist(pidist<-1*pi)=[];

        clear idx
        for k=1:length(pidist)
            [~,idx(k)]=min(abs(pidist(k)-thetaphase));
        end
        idx=sort(idx);

        thetaphaseG = thetaphase(idx);

        subplot(3,4,N)
        hist(thetaphaseG,50)
        h = findobj(gca,'Type','patch');
        set(h,'edgecolor',[0.5 0.5 0.5],'facecolor',[0.5 0.5 0.5])
        set(gca,'ytick',[],'ylim',[0 1500],'xtick',-1*pi:pi:pi,'xticklabel',{'-pi', '0', 'pi'});
        if N==1, xlabel('Phase angle (rad.)'); end
        title(['g = ' num2str(bias2plot(N))]);

        ax=subplot(3,4,N+4);
        angles2plot = randperm(length(idx)); angles2plot=angles2plot(1:100);
        h1=polar([zeros(1,100);thetaphaseG(angles2plot)],repmat([0 1],100,1)','-k'); hold on
        for hh=1:length(h1)
            hasbehavior(h1(hh),'legend',false);
        end
        set(h1,'color',[0.8 0.8 0.8]);
        txt = findall(ax,'type','text'); delete(txt);
        meanvect = mean(exp(1i*thetaphaseG));
        h2=polar([0 angle(meanvect)],[0 abs(meanvect)],'r'); hold on
        set(h2,'linewidth',3);
        text(-1,1.2,['PC = ' num2str(abs(meanvect))])
        
    end
    %
    subplot(3,4,9:12)
    plot(t(1:1000),real(gamma(1:1000)),'r'); hold on
    plot(t(1:1000),real(theta(1:1000)),'k'); 

    theta_shift = theta([ 2*(25):end 1:2*(25)-1]); % -- time-shift theta; step size depends on num_time_shift variable
    plot(t(1:1000),real(theta_shift(1:1000)),'--k');hold on

    theta_shift = theta([ 2*(50):end 1:2*(50)-1]); % -- time-shift theta; step size depends on num_time_shift variable
    plot(t(1:1000),real(theta_shift(1:1000)),':k');hold on

    legend('gamma','theta - 0\pi','theta - \pi/2','theta - \pi');
    box off

end

%% Line plots in case of no coupling (Not in paper)

if ~coupling,
    figure('position',[500 200 800 400])
    
    % -- plot contrast between bias > no bias
    subplot(241)
    plot(1:num_bias_levels,(100.*((PAC - PAC_nobias)./PAC_nobias)),'k','linewidth',1);
    yl = max(get(gca,'ylim'));
    set(gca,'ylim',[-yl/5 yl],'xlim',[1 num_bias_levels]);
    ylabel('Coupling value')
    title('PAC')
    box off
    
    subplot(242)
    plot(1:num_bias_levels,(100.*((dPAC - dPAC_nobias)./dPAC_nobias)),'k','linewidth',1);
    yl = max(get(gca,'ylim'));
    set(gca,'ylim',[-yl/5 yl],'xlim',[1 num_bias_levels]);
    title('dPAC')
    xlabel('Phase bias level');
    box off
    
    subplot(243)
    plot(1:num_bias_levels,(100.*((MI - MI_nobias)./MI_nobias)),'k','linewidth',1); 
    yl = max(get(gca,'ylim'));
    set(gca,'ylim',[-yl/5 yl],'xlim',[1 num_bias_levels]);
    set(gca,'xlim',[1 num_bias_levels]);
    title('MI')
    box off
    
    subplot(244)
    plot(1:num_bias_levels,(100.*((PLV - PLV_nobias)./PLV_nobias)),'k','linewidth',1); 
    yl = max(get(gca,'ylim'));
    set(gca,'ylim',[-yl/5 yl],'xlim',[1 num_bias_levels]);
    set(gca,'xlim',[1 num_bias_levels]);
    title('PLV')
    box off

    % -- plot z-values of permutation test
    threshold_001 = icdf('normal',1-0.001/2,0,1); % -- z-value corresponding to 99.9% or p = 0.001 of z-distribution
    threshold_05 = icdf('normal',1-0.05/2,0,1); % -- z-value corresponding to 95% or p = 0.05 of z-distribution

    subplot(245)
    plot(1:num_bias_levels,PACz,'k','linewidth',1); hold on
    plot([1 num_bias_levels],[threshold_001 threshold_001],'--k'); % -- 3.3 is the z-value threshold corresponding to p=0.001
    plot([1 num_bias_levels],[-threshold_001 -threshold_001],'--k');
    plot([1 num_bias_levels],[threshold_05 threshold_05],'--r'); % -- 1.96 is the z-value threshold corresponding to p=0.05
    plot([1 num_bias_levels],[-threshold_05 -threshold_05],'--r');
    set(gca,'xlim',[1 num_bias_levels]);
    ylabel('Z-value')
    title('PACz')
    box off
    
    subplot(246)
    plot(1:num_bias_levels,dPACz,'k','linewidth',1); hold on
    plot([1 num_bias_levels],[threshold_001 threshold_001],'--k'); % -- 3.3 is the z-value threshold corresponding to p=0.001
    plot([1 num_bias_levels],[-threshold_001 -threshold_001],'--k');
    plot([1 num_bias_levels],[threshold_05 threshold_05],'--r'); % -- 1.96 is the z-value threshold corresponding to p=0.05
    plot([1 num_bias_levels],[-threshold_05 -threshold_05],'--r');
    set(gca,'xlim',[1 num_bias_levels]);
    xlabel('Phase bias level')
    title('dPACz')
    box off
    
    subplot(247)
    h1=plot(1:num_bias_levels,MIz,'k','linewidth',1); hold on
    h2=plot([1 num_bias_levels],[threshold_001 threshold_001],'--k'); % -- 3.3 is the z-value threshold corresponding to p=0.001
    h3=plot([1 num_bias_levels],[threshold_05 threshold_05],'--r'); % -- 1.96 is the z-value threshold corresponding to p=0.05
    hasbehavior(h1,'legend',false); hasbehavior(h2,'legend',false);  hasbehavior(h3,'legend',false);
    plot([1 num_bias_levels],[-threshold_001 -threshold_001],'--k');
    plot([1 num_bias_levels],[-threshold_05 -threshold_05],'--r');
    set(gca,'xlim',[1 num_bias_levels]);
    %legend('p = 0.001','p = 0.05')
    title('MIz')
    box off
       
    subplot(248)
    h1=plot(1:num_bias_levels,MIz,'k','linewidth',1); hold on
    h2=plot([1 num_bias_levels],[threshold_001 threshold_001],'--k'); % -- 3.3 is the z-value threshold corresponding to p=0.001
    h3=plot([1 num_bias_levels],[threshold_05 threshold_05],'--r'); % -- 1.96 is the z-value threshold corresponding to p=0.05
    hasbehavior(h1,'legend',false); hasbehavior(h2,'legend',false);  hasbehavior(h3,'legend',false);
    plot([1 num_bias_levels],[-threshold_001 -threshold_001],'--k');
    plot([1 num_bias_levels],[-threshold_05 -threshold_05],'--r');
    set(gca,'xlim',[1 num_bias_levels]);
    legend('p = 0.001','p = 0.05')
    title('PLVz')
    box off
    
end
