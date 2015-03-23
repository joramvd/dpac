%% dPAC: A method for debiasing phase-amplitude cross-frequency coupling
% Joram van Driel, Roy Cox & Mike X Cohen
% 2014
% --
% This code accompanies the paper titled "dPAC: A method for debiasing 
% phase-amplitude cross-frequency coupling". Below, simulations are run to
% test three phase-amplitude cross-frequency copuling measures (PAC, dPAC
% and MI) as a function of phase clustering bias, coupling angle, and noise
% Using the code without following the paper may lead to confusion, 
% incorrect data analyses, and misinterpretations of results. 
% The authors assume no responsibility for inappropriate or incorrect use 
% of this code. 


clear, close all

cd('path\to\folder'); % -- change directory if you want to save output and/or plots

%% Ingredients

save_results = true;
coupling = true;    % -- set to false to run simulation without coupling (to test for false positives)
permutation = true; % -- run a permutation test? slows down simulation!
nperm = 1000;       % -- number of permutations; set to lower number to speed up analysis
            
srate = 1000;              % -- sampling rate
t = 1/srate:1/srate:10;    % -- time points: 10000 1ms steps; 10 seconds

fpow = 30;                 % -- frequency for power: gamma
fphase = 5;                % -- frequency for phase: theta
ntimepoints = length(t);   % -- get number of timepoints

noiseincrease = linspace(0,100,50);  % -- parameter used for gaussian variance, to simulate white noise

%% Initialize variables

num_time_shifts = 51; % -- angle-shift loop-resolution; set this to a lower value to speed up the analysis
time_shifts = round(linspace(1,51,num_time_shifts));

num_noise_levels = 50; % -- noise-level loop-resolution; set this to a lower value to speed up the analysis
noise_levels = round(linspace(1,50,num_noise_levels));

 % -- if the simulation is run with non-coupled signals, temporally shifting theta will have no effect, so we can skip this part
if ~coupling, 
    num_time_shifts = 1;
    time_shifts = 0; 
end

[PAC,dPAC,MI, PACz,dPACz,MIz] = deal( zeros(num_noise_levels,num_time_shifts) );  % -- initialize output matrices: 50 clustering levels and 51 PAC angles  
[PAC_nonoise,dPAC_nonoise,MI_nonoise] = deal( zeros(1,num_time_shifts) ); % -- PAC/dPAC/MI without noise, with bias; only change over coupling angle
[PAC_nobias,dPAC_nobias,MI_nobias] = deal( zeros(1,num_noise_levels) ); % -- PAC/dPAC/MI without bias, with noise; only change over noise levels


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

% -- compute PAC, dPAC and MI in the "pure" coupling case 
% -- note that the debias term [mean(exp(1i*thetaphase))] is now 
% -- incorporated in dPAC into one line of code
PAC_nobias_nonoise  = abs(mean(exp(1i*thetaphase) .* gammapower));
dPAC_nobias_nonoise = abs(mean( (exp(1i*thetaphase) - mean(exp(1i*thetaphase))) .* gammapower));

nbins = 18;
thetaphase_bin = ceil( tiedrank( thetaphase ) / (ntimepoints / nbins) );
gammapower_bin = zeros(1,nbins);
for k=1:nbins
    gammapower_bin(k) = squeeze(mean(gammapower(thetaphase_bin==k)));
end
gammapower_bin = gammapower_bin ./ sum(gammapower_bin);
MI_nobias_nonoise = (log(nbins) + sum(gammapower_bin.*log(gammapower_bin)) ) ./ log(nbins);

% -- in this simulation, we'll set the phase clustering bias at a fixed
% -- level
% -- create G-shaped phase distribution
gaussiandist = randn(1,ntimepoints*2)*1.5; % -- gaussian distribution with a fixed width
pidist = gaussiandist;
pidist(pidist<-pi | pidist>pi) = []; % -- cut the distribution off at -pi and pi
    

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
        thetaphase = angle(theta_shift);
    end
    
    % -- use gauss-based angle-distribution to sample from theta phase
    % -- distribution
    clear idx
    for k=1:length(pidist)
        [~,idx(k)]=min(abs(pidist(k)-thetaphase));
    end
    idx=sort(idx);
    
    % -- in below variables, G stands for 'gaussian'
    thetaphaseG = thetaphase(idx); % -- theta phase angles now have clustering around 0pi
    gammapowerG = gammapower(idx); % -- save the corresponding gamma power values
    ntimepointsG = length(idx);
    
    % -- compute PAC, dPAC and MI, now with phase clustering bias; no noise yet
    PAC_nonoise(tii)  = abs(mean(exp(1i*thetaphaseG) .* gammapowerG));
    dPAC_nonoise(tii) = abs(mean( (exp(1i*thetaphaseG) - mean(exp(1i*thetaphaseG))) .* gammapowerG));
        
    thetaphaseG_bin = ceil( tiedrank( thetaphaseG ) / (ntimepointsG / nbins) );
    gammapowerG_bin = zeros(1,nbins);
    for k=1:nbins
        gammapowerG_bin(k) = squeeze(mean(gammapowerG(thetaphaseG_bin==k)));
    end
    gammapowerG_bin = gammapowerG_bin ./ sum(gammapowerG_bin);
    MI_nonoise(tii) = (log(nbins) + sum(gammapowerG_bin.*log(gammapowerG_bin)) ) ./ log(nbins);

    %%
    % -- second loop: over different levels of noise added to gamma

    nii = 0; % -- initialize loop-counter variables
    for ni = noise_levels % -- ni is the index of noise strength
        nii = nii+1; % -- increment counter
        
        % -- display progress
        msg = sprintf('Running noise-level %i/%i anglediff %i/%i...',  nii,num_noise_levels,tii,num_time_shifts);
        fprintf([reverseStr, msg]);
        reverseStr = repmat(sprintf('\b'), 1, length(msg));
        
        noise = noiseincrease(ni) *randn(1,ntimepoints);
        
        % -- band-pass filter broadband noise
        band           = [25 35]; % -- we filter the noise to be in a "gamma" band of 25-35 Hz
        filt_order     = round(3*(srate/band(1)));
        filterweights  = fir1(filt_order,[mean(band)-(band(2)-band(1)) mean(band)+(band(2)-band(1))]/(srate/2));
        noisefilt = filtfilt(filterweights,1,noise); % -- uses the zero-phase-lag filtfilt function of Matlab
        
        if coupling
            gamma = ((real(theta)+6)  .* exp( 1i * 2 * pi * fpow * t )) + hilbert(noisefilt); % -- complex gamma is phase-modulated by theta and combined with hilbert transform of noise
        else
            gamma = exp( 1i * 2 * pi * fpow * t )  + hilbert(noisefilt); % -- only gamma with noise; i.e. no coupling
        end

        gammapower_noise = abs(gamma); 
        gammapowerG_noise = gammapower_noise(idx); % -- gamma power values from the clustering distribution
        
        if tii==1
            % -- compute PAC, dPAC and MI, now without phase clustering bias, with noise
            PAC_nobias(nii)  = abs(mean(exp(1i*thetaphase) .* gammapower_noise));
            dPAC_nobias(nii) = abs(mean( (exp(1i*thetaphase) - mean(exp(1i*thetaphase))) .* gammapower_noise));
            
            gammapower_noise_bin = zeros(1,nbins);
            for k=1:nbins
                gammapower_noise_bin(k) = squeeze(mean(gammapower_noise(thetaphase_bin==k)));
            end
            gammapower_noise_bin = gammapower_noise_bin ./ sum(gammapower_noise_bin);
            MI_nobias(nii) = (log(nbins) + sum(gammapower_noise_bin.*log(gammapower_noise_bin)) ) ./ log(nbins);
        end
        
        % -- compute PAC, dPAC and MI, now with phase clustering bias
        PAC(nii,tii)  = abs(mean(exp(1i*thetaphaseG) .* gammapowerG_noise));
        dPAC(nii,tii) = abs(mean( (exp(1i*thetaphaseG) - mean(exp(1i*thetaphaseG))) .* gammapowerG_noise));
                        
        thetaphaseG_bin = ceil( tiedrank( thetaphaseG ) / (ntimepointsG / nbins) );
        gammapowerG_noise_bin = zeros(1,nbins);
        for k=1:nbins
            gammapowerG_noise_bin(k) = squeeze(mean(gammapowerG_noise(thetaphaseG_bin==k)));
        end
        gammapowerG_noise_bin = gammapowerG_noise_bin ./ sum(gammapowerG_noise_bin);
        MI(nii,tii) = (log(nbins) + sum(gammapowerG_noise_bin.*log(gammapowerG_noise_bin)) ) ./ log(nbins);

        %% now with permutation testing
        
        if permutation
            
            [fake_PAC,fake_dPAC,fake_MI] = deal(zeros(1,nperm));
            
            for permi = 1:nperm
                
                % -- cut-and-paste a random portion of the data; this preserves
                % -- temporal autocorrelation while removing the coupling
                cutLoc = 5 + randperm(ntimepointsG-10); % -- 5 and 10 prevent the first and last time points from being selected
                cutLoc = cutLoc(1);
                thetaphaseG_shuffled = thetaphaseG([cutLoc:end 1:cutLoc-1]);
                
                fake_PAC(permi)  = abs(mean(exp(1i*thetaphaseG_shuffled) .* gammapowerG_noise)); % -- compute surrogate PAC
                fake_dPAC(permi) = abs(mean( (exp(1i*thetaphaseG_shuffled) - mean(exp(1i*thetaphaseG_shuffled))) .* gammapowerG_noise)); % -- compute surrogate dPAC
                
                % -- compute MI (Tort et al., 2010)
                thetaphaseG_bin_shuffled = ceil( tiedrank( thetaphaseG_shuffled ) / (ntimepointsG / nbins) );
                gammapowerG_noise_bin = zeros(1,nbins);
                for k=1:nbins
                    gammapowerG_noise_bin(k) = squeeze(mean(gammapowerG_noise(thetaphaseG_bin_shuffled==k)));
                end
                gammapowerG_noise_bin = gammapowerG_noise_bin ./ sum(gammapowerG_noise_bin);
                fake_MI(permi) = (log(nbins) + sum(gammapowerG_noise_bin.*log(gammapowerG_noise_bin)) ) ./ log(nbins);
                
            end
            
            % -- below, the zscore is defined by observed value minus mean of
            % -- null-distribution, divided by the standard deviation of
            % -- null-distribution
            PACz(nii,tii)  = (squeeze(PAC(nii,tii))' - squeeze(mean(fake_PAC))) ./ squeeze(std(fake_PAC)); 
            dPACz(nii,tii) = (squeeze(dPAC(nii,tii))' - squeeze(mean(fake_dPAC))) ./ squeeze(std(fake_dPAC)); 
            MIz(nii,tii)   = (squeeze(MI(nii,tii))' - squeeze(mean(fake_MI))) ./ squeeze(std(fake_MI));
        end
    end
    
end
fprintf(' done\n');

%% save results

if save_results
    if coupling
        filename = 'dPAC_MI_simresults_noiseincrease.mat';
    else
        filename = 'dPAC_MI_simresults_noiseincrease_NOcoupling.mat';
    end
    save(filename,...
        'PAC_nobias', 'PAC_nonoise', 'PAC_nobias_nonoise',...
        'dPAC_nobias', 'dPAC_nonoise', 'dPAC_nobias_nonoise',...
        'MI_nobias', 'MI_nonoise', 'MI_nobias_nonoise',...
        'PAC', 'dPAC', 'MI', 'PACz', 'dPACz', 'MIz',...
        'num_noise_levels', 'num_time_shifts');
    
end

%% Plot the result (Figure 6) -- only when coupling set to true

if coupling,


    % -- in the paper, a particular colormap and colorscaling is used that 
    % -- necessitates some third-party functions that are provided with the 
    % -- code;
    % -- set to false to use default jet map

    plot_like_paper = true;

    % -- when set to true, the othercolor function is used, which can be downloaded here: 
    % -- [http://nl.mathworks.com/matlabcentral/fileexchange/30564-othercolor]
    % -- the diverging_map function can be downloaded here:
    % -- [http://www.sandia.gov/~kmorel/documents/ColorMaps/diverging_map.m]

%%
    %%%%%%
    % -- below, noise levels of PAC/dPAC/MI are compared with situation without
    % -- noise nor bias; to compare with situation without bias but 
    % -- with noise, or without noise but with bias, comment-out the desired line below:

    % comparePAC = repmat(PAC_nonoise,[num_noise_levels 1]); comparedPAC = repmat(dPAC_nonoise,[num_noise_levels 1]); compareMI = repmat(MI_nonoise,[num_noise_levels 1]);
    % comparePAC = repmat(PAC_nobias, [num_time_shifts  1])'; comparedPAC = repmat(dPAC_nobias, [num_time_shifts  1])'; compareMI = repmat(MI_nobias, [num_time_shifts  1])';
    comparePAC = PAC_nobias_nonoise; comparedPAC = dPAC_nobias_nonoise; compareMI = MI_nobias_nonoise; 

    angleaxis = 0:pi/(num_time_shifts-1):pi; % -- the x-axis will show coupling angles in radian

    figure('position',[600 300 1000 400]); % -- these values may need to be changed depending on screen settings

    % -- plot PAC; from every angle-by-phase-clustering point, PAC without
    % -- noise, but with fixed bias is subtracted
    subplot(231) 
    contourf(angleaxis,1:num_noise_levels,(PAC - comparePAC),40,'linestyle','none'); % -- here, contourf is used for smooth contours; alternatively, you can use imagesc (which requires flipping the y-axis direction)
    cl = max(abs(get(gca,'clim')));
    set(gca,'clim',[-cl cl],...
        'ytick',[1 num_noise_levels],'yticklabel',{'min','max'},...
        'xtick',0:pi/2:pi,'xticklabel',{'0', 'pi/2', 'pi'}); colorbar
    ylabel('Gamma noise')
    title('PAC')

     % -- the same now for dPAC
    subplot(232)
    contourf(angleaxis,1:num_noise_levels,(dPAC - comparedPAC),40,'linestyle','none');
    set(gca,'clim',[-cl cl],...
        'ytick',[1 num_noise_levels],'yticklabel',{},...
        'xtick',0:pi/2:pi,'xticklabel',{'0', 'pi/2', 'pi'}); colorbar
    xlabel('Clustering angle (rad.)')
    title('dPAC')

     % -- the same now for MI
    subplot(233)
    contourf(angleaxis,1:num_noise_levels,(MI - compareMI),40,'linestyle','none');
    cl = max(abs(get(gca,'clim')));
    set(gca,'clim',[-cl cl],...
        'ytick',[1 num_noise_levels],'yticklabel',{},...
        'xtick',0:pi/2:pi,'xticklabel',{'0', 'pi/2', 'pi'}); colorbar
    title('MI')

    if plot_like_paper, 
        addpath(genpath('path\to\dPAC\scripts\plot_like_paper\')); % -- change path accordingly
        colormap(othercolor('BuDRd_18'))
    end
    
%%
    % -- plot PACz
    if plot_like_paper, figure('position',[600 300 1000 400]); end % -- new figure needed to have different colormaps

    subplot(234) % -- PACz
    contourf(angleaxis,1:num_noise_levels,PACz,40,'linestyle','none');
    set(gca,'clim',[-4 4],...
        'ytick',[1 num_noise_levels],'yticklabel',{'min','max'},...
        'xtick',0:pi/2:pi,'xticklabel',{'0', 'pi/2', 'pi'}); colorbar
    ylabel('Gamma noise')
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
        contour(angleaxis,1:num_noise_levels,PACz_thres_001,1,'k','LineWidth',1) % -- plot p<0.001 clusters as overlaid black line
        [r,c]=bwlabel(PACz_thres_05,4); % -- function that searches for 2D clusters
        for ci=1:c
            if sum(any(r==ci,1))<5 || sum(any(r==ci,2))<5 % -- arbitrary constraint: at least 5 contiguous points
                PACz_thres_05(r==ci)=0;
            end
        end
        hold on
        contour(angleaxis,1:num_noise_levels,PACz_thres_05,1,'r','LineWidth',1) % -- plot p<0.05 clusters as overlaid red line
    end
    
    subplot(235) % -- dPACz
    contourf(angleaxis,1:num_noise_levels,dPACz,40,'linestyle','none');
    set(gca,'clim',[-4 4],...
        'ytick',[1 num_noise_levels],'yticklabel',{},...
        'xtick',0:pi/2:pi,'xticklabel',{'0', 'pi/2', 'pi'}); colorbar
    xlabel('Clustering angle (rad.)')
    title('dPACz')

    dPACz_thres_001 = squeeze(dPACz);
    dPACz_thres_001(dPACz_thres_001<threshold_001 & dPACz_thres_001>(-1*threshold_001))=0; % -- binary threshold: set all values with p>0.001 to zero
    dPACz_thres_001(dPACz_thres_001~=0)=1; % -- everything else to 1
    dPACz_thres_05 = squeeze(dPACz);
    dPACz_thres_05(dPACz_thres_05<threshold_05 & dPACz_thres_05>(-1*threshold_05))=0; % -- binary threshold: set all values with p>0.05 to zero
    dPACz_thres_05(dPACz_thres_05~=0)=1; % -- everything else to 1

    % -- cluster size thresholding (image processing toolbox required)
    if exist('bwlabel','file')
        [r,c]=bwlabel(dPACz_thres_001,4); % -- function that searches for 2D clusters
        for ci=1:c
            if sum(any(r==ci,1))<5 || sum(any(r==ci,2))<5 % -- arbitrary constraint: at least 5 contiguous points
                dPACz_thres_001(r==ci)=0;
            end
        end
        hold on
        contour(angleaxis,1:num_noise_levels,dPACz_thres_001,1,'k','LineWidth',1) % -- plot p<0.001 clusters as overlaid black line
        [r,c]=bwlabel(dPACz_thres_05,4); % -- function that searches for 2D clusters
        for ci=1:c
            if sum(any(r==ci,1))<5 || sum(any(r==ci,2))<5 % -- arbitrary constraint: at least 5 contiguous points
                dPACz_thres_05(r==ci)=0;
            end
        end
        hold on
        contour(angleaxis,1:num_noise_levels,dPACz_thres_05,1,'r','LineWidth',1) % -- plot p<0.05 clusters as overlaid red line
    end
    
    subplot(236) % -- MIz
    contourf(angleaxis,1:num_noise_levels,MIz,40,'linestyle','none');
    set(gca,'clim',[-4 4],...
        'ytick',[1 num_noise_levels],'yticklabel',{},...
        'xtick',0:pi/2:pi,'xticklabel',{'0', 'pi/2', 'pi'}); colorbar
    title('MIz')

    MIz_thres_001 = squeeze(MIz);
    MIz_thres_001(MIz_thres_001<threshold_001 & MIz_thres_001>(-1*threshold_001))=0; % -- binary threshold: set all values with p>0.001 to zero
    MIz_thres_001(MIz_thres_001~=0)=1; % -- everything else to 1
    MIz_thres_05 = squeeze(MIz);
    MIz_thres_05(MIz_thres_05<threshold_05 & MIz_thres_05>(-1*threshold_05))=0; % -- binary threshold: set all values with p>0.05 to zero
    MIz_thres_05(MIz_thres_05~=0)=1; % -- everything else to 1

    % -- cluster size thresholding (image processing toolbox required)
    if exist('bwlabel','file')
        [r,c]=bwlabel(MIz_thres_001,4); % -- function that searches for 2D clusters
        for ci=1:c
            if sum(any(r==ci,1))<5 || sum(any(r==ci,2))<5 % -- arbitrary constraint: at least 5 contiguous points
                MIz_thres_001(r==ci)=0;
            end
        end
        hold on
        contour(angleaxis,1:num_noise_levels,MIz_thres_001,1,'k','LineWidth',1) % -- plot p<0.001 clusters as overlaid black line
        [r,c]=bwlabel(MIz_thres_05,4); % -- function that searches for 2D clusters
        for ci=1:c
            if sum(any(r==ci,1))<5 || sum(any(r==ci,2))<5 % -- arbitrary constraint: at least 5 contiguous points
                MIz_thres_05(r==ci)=0;
            end
        end
        hold on
        contour(angleaxis,1:num_noise_levels,MIz_thres_05,1,'r','LineWidth',1) % -- plot p<0.05 clusters as overlaid red line
    end
    
    
    if plot_like_paper, 
        oldmap = colormap;
        newmap = diverging_map(0:1/length(oldmap):1,[0 0 1],[0 1 0]);
        colormap(newmap)
    end

    %% Plot noise manipulation illustration (Figure 5)

    figure('position',[700 200 400 600]);

    gamma = ((real(theta)+6)  .* exp( 1i * 2 * pi * fpow * t )); % -- gamma is phase-modulated by theta
    subplot(411)
    plot(t,real(gamma),'r'); hold on
    plot(t,abs(gamma),'k'); box off
    set(gca,'xlim',[t(1) t(1000)],'xticklabel',{},'ylim',[-30 30])
    title('Pure theta-gamma coupling');

    noise = noiseincrease(25) *randn(1,ntimepoints);
    subplot(412);
    plot(t,noise,'k'); box off
    set(gca,'xlim',[t(1) t(1000)],'xticklabel',{})
    title('Broadband noise')

    % -- band-pass filter broadband noise
    band           = [25 35]; % -- we filter the noise to be in a "gamma" band of 25-35 Hz
    filt_order     = round(3*(srate/band(1)));
    filterweights  = fir1(filt_order,[mean(band)-(band(2)-band(1)) mean(band)+(band(2)-band(1))]/(srate/2));
    noisefilt = filtfilt(filterweights,1,noise); % -- uses the zero-phase-lag filtfilt function of Matlab
    
    subplot(413);
    plot(t,noisefilt,'k'); box off
    set(gca,'xlim',[t(1) t(1000)],'xticklabel',{},'ylim',[-30 30])
    title('Gamma filtered noise')

    gamma = ((real(theta)+6)  .* exp( 1i * 2 * pi * fpow * t ));
    subplot(414)
    plot(t,real(gamma) + noisefilt,'r'); hold on
    plot(t,abs(gamma + hilbert(noisefilt)),'k'); box off
    set(gca,'xlim',[t(1) t(1000)],'ylim',[-30 30])
    ylabel('Amplitude')
    xlabel('Time (s)')
    title('Theta modulated gamma with noise')

end

%% Line plots in case of no coupling (Not in paper) 

if ~coupling,
    figure('position',[500 200 800 400])
    
    % -- plot contrast between bias > no bias
    subplot(231)
    plot(1:num_noise_levels,(PAC - PAC_nobias_nonoise),'k','linewidth',1);
    yl = max(get(gca,'ylim'));
    set(gca,'ylim',[-yl/5 yl],'xlim',[1 num_noise_levels]);
    ylabel('Coupling value (diff.)')
    title('PAC')
    box off
    
    subplot(232)
    plot(1:num_noise_levels,(dPAC - dPAC_nobias_nonoise),'k','linewidth',1);
    set(gca,'ylim',[-yl/5 yl],'xlim',[1 num_noise_levels]);
    title('dPAC')
    xlabel('Noise level');
    box off
    
    subplot(233)
    plot(1:num_noise_levels,(MI - MI_nobias_nonoise),'k','linewidth',1);
    yl = max(get(gca,'ylim'));
    set(gca,'ylim',[-yl/5 yl],'xlim',[1 num_noise_levels]);
    title('MI')
    box off
    
    % -- plot z-values of permutation test
    threshold_001 = icdf('normal',1-0.001/2,0,1); % -- z-value corresponding to 99.9% or p = 0.001 of z-distribution
    threshold_05 = icdf('normal',1-0.05/2,0,1); % -- z-value corresponding to 95% or p = 0.05 of z-distribution
    
    subplot(234)
    plot(1:num_noise_levels,PACz,'k','linewidth',1); hold on
    plot([1 num_noise_levels],[threshold_001 threshold_001],'--k'); % -- 3.3 is the z-value threshold corresponding to p=0.001
    plot([1 num_noise_levels],[-threshold_001 -threshold_001],'--k');
    plot([1 num_noise_levels],[threshold_05 threshold_05],'--r'); % -- 1.96 is the z-value threshold corresponding to p=0.05
    plot([1 num_noise_levels],[-threshold_05 -threshold_05],'--r');
    set(gca,'xlim',[1 num_noise_levels]);
    ylabel('Z-value')
    title('PACz')
    box off
    
    subplot(235)
    plot(1:num_noise_levels,dPACz,'k','linewidth',1); hold on
    plot([1 num_noise_levels],[threshold_001 threshold_001],'--k'); % -- 3.3 is the z-value threshold corresponding to p=0.001
    plot([1 num_noise_levels],[-threshold_001 -threshold_001],'--k');
    plot([1 num_noise_levels],[threshold_05 threshold_05],'--r'); % -- 1.96 is the z-value threshold corresponding to p=0.05
    plot([1 num_noise_levels],[-threshold_05 -threshold_05],'--r');
    set(gca,'xlim',[1 num_noise_levels]);
    xlabel('Noise level');
    title('dPACz')
    box off
    
    subplot(236)
    h1=plot(1:num_noise_levels,MIz,'k','linewidth',1); hold on
    h2=plot([1 num_noise_levels],[threshold_001 threshold_001],'--k'); % -- 3.3 is the z-value threshold corresponding to p=0.001
    h3=plot([1 num_noise_levels],[threshold_05 threshold_05],'--r'); % -- 1.96 is the z-value threshold corresponding to p=0.05
    hasbehavior(h1,'legend',false); hasbehavior(h2,'legend',false);  hasbehavior(h3,'legend',false);
    plot([1 num_noise_levels],[-threshold_001 -threshold_001],'--k');
    plot([1 num_noise_levels],[-threshold_05 -threshold_05],'--r');
    set(gca,'xlim',[1 num_noise_levels]);
    legend('p = 0.001','p = 0.05')
    title('MIz')
    box off
    
end
