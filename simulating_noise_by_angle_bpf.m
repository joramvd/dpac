%% dPAC: A method for debiasing phase-amplitude cross-frequency coupling
% Joram van Driel, Roy Cox & Mike X Cohen
% 2014/2015
% --
% This code accompanies the paper titled "Phase clustering bias in
% phase-amplitude cross-frequency coupling and its removal". Below, 
% simulations are run to test four phase-amplitude cross-frequency 
% copuling measures (PAC, dPAC MI and PLV) as a function of pink or white
% noise and coupling angle.
% In addition, band-pass filtering is used to approximate regular time
% series analysis approaches to real data.
% Using the code without following the paper may lead to confusion, 
% incorrect data analyses, and misinterpretations of results. 
% The authors assume no responsibility for inappropriate or incorrect use 
% of this code. 


clear, close all

cd('Z:\PhD\dPAC\'); % -- change directory if you want to save output and/or plots

%% Ingredients

save_results = false;
coupling = true;    % -- set to false to run simulation without coupling (to test for false positives)
permutation = true; % -- run a permutation test? slows down simulation!
nperm = 1000;       % -- number of permutations; set to lower number to speed up analysis
            
srate = 1000;              % -- sampling rate
t = 1/srate:1/srate:12;    % -- time points: 10000 1ms steps; 12 seconds (10 seconds plus padded 1000 ms on each side for edge artifacts)

fpow = 30;                 % -- frequency for power: gamma
fphase = 5;                % -- frequency for phase: theta
ntimepoints = length(t)-2000;   % -- get number of timepoints (subtract buffer zone)

type_of_noise = 'pink'; % -- change this to 'white' to do the simulation with white noise; the results are not very different, but pink noise is a slightly more similar to noise as seen in EEG
% -- note that the pink noise simulation calls a third-party Matlab
% -- function called pinknoise.m; this function is provided along with the
% -- simulation code, and was downloaded from 
% -- [http://www.mathworks.com/matlabcentral/fileexchange/42919]

noiseincrease = linspace(0,10,50);  % -- parameter used for gaussian variance, to simulate white noise

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

[PAC,dPAC,MI,PLV, PACz,dPACz,MIz,PLVz]            = deal( zeros(num_noise_levels,num_time_shifts) );  % -- initialize output matrices: 50 clustering levels and 51 PAC angles  
[PAC_nonoise,dPAC_nonoise,MI_nonoise,PLV_nonoise] = deal( zeros(1,num_time_shifts) ); % -- PAC/dPAC/MI/PLV without noise, with bias; only change over coupling angle
[PAC_nobias,dPAC_nobias,MI_nobias,PLV_nobias]     = deal( zeros(1,num_noise_levels) ); % -- PAC/dPAC/MI/PLV without bias, with noise; only change over noise levels


%% Simulation

% -- create complex sine waves
theta = 5.*(exp( 1i * 2 * pi * fphase * t ));
if coupling
    gamma = ((real(theta)+6)  .* exp( 1i * 2 * pi * fpow * t )); % -- gamma is phase-modulated by theta; i.e. pure cross-frequency coupling
else
    gamma = exp( 1i * 2 * pi * fpow * t ); % -- pure gamma sine wave; i.e. no coupling
end

thetagamma = theta+gamma;

%% basic filter settings

% -- band-pass filter
thetaband           = [3 7];
gammaband           = [25 35];

theta_filt_order    = round(3*(srate/thetaband(1)));
theta_filterweights = fir1(theta_filt_order,[mean(thetaband)-(thetaband(2)-thetaband(1)) mean(thetaband)+(thetaband(2)-thetaband(1))]/(srate/2));

gamma_filt_order    = round(3*(srate/gammaband(1)));
gamma_filterweights = fir1(gamma_filt_order,[mean(gammaband)-(gammaband(2)-gammaband(1)) mean(gammaband)+(gammaband(2)-gammaband(1))]/(srate/2));

thetafilt = filtfilt(theta_filterweights,1,real(thetagamma));
gammafilt = filtfilt(gamma_filterweights,1,real(thetagamma));


%%

% -- compute theta phase angles and gamma power
% -- first and last 1000 pnts are removed to account for edge artifacts
thetaphase = angle(hilbert(thetafilt(1001:end-1000)));
gammapower = abs(hilbert(gammafilt(1001:end-1000)));

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

% -- PLV (Cohen, 2008; Colgin et al 2009)
gammapower_phase = angle(hilbert(detrend(gammapower))); % -- note: hilbert transform of power envelope of complex signal
PLV_nobias_nonoise = abs(mean(exp(1i*(thetaphase-gammapower_phase)))); % -- gives PLV of 1.0


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
        thetagamma = theta_shift+gamma;
        thetafilt = filtfilt(theta_filterweights,1,real(thetagamma));
        gammafilt = filtfilt(gamma_filterweights,1,real(thetagamma));

        thetaphase = angle(hilbert(thetafilt(1001:end-1000)));
        gammapower = abs(hilbert(gammafilt(1001:end-1000)));
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

    % -- PLV value (Cohen, 2008; Colgin et al 2009)
    gammapowerG_phase = angle(hilbert(detrend(gammapowerG))); % -- note: hilbert transform of power envelope of complex signal
    PLV_nonoise(tii) = abs(mean(exp(1i*(thetaphaseG-gammapowerG_phase))));

    %%
    % -- second loop: over different levels of noise added to gamma

    nii = 0; % -- initialize loop-counter variables
    for ni = noise_levels % -- ni is the index of noise strength
        nii = nii+1; % -- increment counter
        
        % -- display progress
        msg = sprintf('Running noise-level %i/%i anglediff %i/%i...',  nii,num_noise_levels,tii,num_time_shifts);
        fprintf([reverseStr, msg]);
        reverseStr = repmat(sprintf('\b'), 1, length(msg));
        
        if strcmp(type_of_noise,'white');
            noise = noiseincrease(ni) *randn(1,length(t));
        elseif strcmp(type_of_noise,'pink');
            noise = noiseincrease(ni) *pinknoise(length(t)); % -- calls the pinknoise.m function, see above
        end
        
        % -- first, all signals are added together: theta, theta-modulated
        % -- gamma, and pink/white noise
        thetagammanoise = real(thetagamma)+noise;
        
        % -- next, theta and gamma are recreated through band-pass
        % -- filtering
        thetafilt = filtfilt(theta_filterweights,1,thetagammanoise);
        gammafilt = filtfilt(gamma_filterweights,1,thetagammanoise);

        % -- phase and power is extracted, first and last 1000 pnts are
        % -- removed to account for edge artifacts
        thetaphase = angle(hilbert(thetafilt(1001:end-1000)));
        gammapower = abs(hilbert(gammafilt(1001:end-1000)));

        thetaphaseG = thetaphase(idx);
        gammapowerG = gammapower(idx); % -- gamma power values from the clustering distribution
        
        if tii==1
            % -- compute PAC, dPAC and MI, now without phase clustering bias, with noise
            PAC_nobias(nii)  = abs(mean(exp(1i*thetaphase) .* gammapower));
            dPAC_nobias(nii) = abs(mean( (exp(1i*thetaphase) - mean(exp(1i*thetaphase))) .* gammapower));
            
            gammapower_bin = zeros(1,nbins);
            for k=1:nbins
                gammapower_bin(k) = squeeze(mean(gammapower(thetaphase_bin==k)));
            end
            gammapower_bin = gammapower_bin ./ sum(gammapower_bin);
            MI_nobias(nii) = (log(nbins) + sum(gammapower_bin.*log(gammapower_bin)) ) ./ log(nbins);
            
            % -- PLV (Cohen, 2008; Colgin et al 2009)
            gammapower_phase = angle(hilbert(detrend(gammapower))); % -- note: hilbert transform of power envelope of complex signal
            PLV_nobias(nii) = abs(mean(exp(1i*(thetaphase-gammapower_phase)))); % -- gives PLV of 1.0

        end
        
        % -- compute PAC, dPAC and MI, now with phase clustering bias
        PAC(nii,tii)  = abs(mean(exp(1i*thetaphaseG) .* gammapowerG));
        dPAC(nii,tii) = abs(mean( (exp(1i*thetaphaseG) - mean(exp(1i*thetaphaseG))) .* gammapowerG));
                        
        thetaphaseG_bin = ceil( tiedrank( thetaphaseG ) / (ntimepointsG / nbins) );
        gammapowerG_bin = zeros(1,nbins);
        for k=1:nbins
            gammapowerG_bin(k) = squeeze(mean(gammapowerG(thetaphaseG_bin==k)));
        end
        gammapowerG_bin = gammapowerG_bin ./ sum(gammapowerG_bin);
        MI(nii,tii) = (log(nbins) + sum(gammapowerG_bin.*log(gammapowerG_bin)) ) ./ log(nbins);

        % -- PLV (Cohen, 2008; Colgin et al 2009)
        gammapowerG_phase = angle(hilbert(detrend(gammapowerG))); % -- note: hilbert transform of power envelope of complex signal
        PLV(nii,tii) = abs(mean(exp(1i*(thetaphaseG-gammapowerG_phase))));

        %% now with permutation testing
        
        if permutation
            
            [fake_PAC,fake_dPAC,fake_MI,fake_PLV] = deal(zeros(1,nperm));
            
            for permi = 1:nperm
                
                % -- cut-and-paste a random portion of the data; this preserves
                % -- temporal autocorrelation while removing the coupling
                cutLoc = 5 + randperm(ntimepointsG-10); % -- 5 and 10 prevent the first and last time points from being selected
                cutLoc = cutLoc(1);
                thetaphaseG_shuffled = thetaphaseG([cutLoc:end 1:cutLoc-1]);
                
                fake_PAC(permi)  = abs(mean(exp(1i*thetaphaseG_shuffled) .* gammapowerG)); % -- compute surrogate PAC
                fake_dPAC(permi) = abs(mean( (exp(1i*thetaphaseG_shuffled) - mean(exp(1i*thetaphaseG_shuffled))) .* gammapowerG)); % -- compute surrogate dPAC
                
                % -- compute MI (Tort et al., 2010)
                thetaphaseG_bin_shuffled = ceil( tiedrank( thetaphaseG_shuffled ) / (ntimepointsG / nbins) );
                gammapowerG_bin = zeros(1,nbins);
                for k=1:nbins
                    gammapowerG_bin(k) = squeeze(mean(gammapowerG(thetaphaseG_bin_shuffled==k)));
                end
                gammapowerG_bin = gammapowerG_bin ./ sum(gammapowerG_bin);
                fake_MI(permi) = (log(nbins) + sum(gammapowerG_bin.*log(gammapowerG_bin)) ) ./ log(nbins);
                
                % -- PLV value (Cohen, 2008; Colgin et al 2009)
                fake_PLV(permi) = abs(mean(exp(1i*(thetaphaseG_shuffled-gammapowerG_phase))));
                
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

%% save results

if save_results
    if coupling
        filename = ['dPAC_MI_simresults_' type_of_noise 'noise_filtered.mat'];
    else
        filename = ['dPAC_MI_simresults_' type_of_noise 'noise_filtered_NOcoupling.mat'];
    end
    save(filename,...
        'PAC_nobias', 'PAC_nonoise', 'PAC_nobias_nonoise',...
        'dPAC_nobias', 'dPAC_nonoise', 'dPAC_nobias_nonoise',...
        'MI_nobias', 'MI_nonoise', 'MI_nobias_nonoise',...
        'PLV_nobias', 'PLV_nonoise', 'PLV_nobias_nonoise',...
        'PAC', 'dPAC', 'MI', 'PLV', 'PACz', 'dPACz', 'MIz', 'PLVz',...
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

    % comparePAC = repmat(PAC_nonoise,[num_noise_levels 1]); comparedPAC = repmat(dPAC_nonoise,[num_noise_levels 1]); compareMI = repmat(MI_nonoise,[num_noise_levels 1]); comparePLV = repmat(PLV_nonoise,[num_noise_levels 1]);
    % comparePAC = repmat(PAC_nobias, [num_time_shifts  1])'; comparedPAC = repmat(dPAC_nobias, [num_time_shifts  1])'; compareMI = repmat(MI_nobias, [num_time_shifts  1])'; comparePLV = repmat(PLV_nobias, [num_time_shifts  1])';
    comparePAC = PAC_nobias_nonoise; comparedPAC = dPAC_nobias_nonoise; compareMI = MI_nobias_nonoise; comparePLV = PLV_nobias_nonoise;

    angleaxis = 0:pi/(num_time_shifts-1):pi; % -- the x-axis will show coupling angles in radian

    figure('position',[600 300 1200 400]); % -- these values may need to be changed depending on screen settings

    % -- plot PAC; from every angle-by-phase-clustering point, PAC without
    % -- noise, but with fixed bias is subtracted and perc. signal change
    % -- computed
    subplot(241) 
    PAC_perc = 100.*((PAC - comparePAC)./comparePAC);
    contourf(angleaxis,1:num_noise_levels,PAC_perc,40,'linestyle','none'); % -- here, contourf is used for smooth contours; alternatively, you can use imagesc (which requires flipping the y-axis direction)
    cl = max(abs(get(gca,'clim')));
    set(gca,'clim',[-cl cl],...
        'ytick',[1 num_noise_levels],'yticklabel',{'min','max'},...
        'xtick',0:pi/2:pi,'xticklabel',{'0', 'pi/2', 'pi'}); colorbar
    ylabel('Gamma noise')
    title('PAC')

    PAC_thresh = zeros(size(PAC_perc));
    PAC_thresh(PAC_perc>2*std(PAC_perc(:)) | PAC_perc<-2*std(PAC_perc(:)))=1;
    hold on
    contour(angleaxis,1:num_noise_levels,PAC_thresh,1,'k','LineWidth',1) % -- plot p<0.001 clusters as overlaid black line

    % -- the same now for dPAC
    subplot(242)
    dPAC_perc = 100.*((dPAC - comparedPAC)./comparedPAC);
    contourf(angleaxis,1:num_noise_levels,dPAC_perc,40,'linestyle','none');
    cl = max(abs(get(gca,'clim')));
    set(gca,'clim',[-cl cl],...
        'ytick',[1 num_noise_levels],'yticklabel',{},...
        'xtick',0:pi/2:pi,'xticklabel',{'0', 'pi/2', 'pi'}); colorbar
    xlabel('Clustering angle (rad.)')
    title('dPAC')

    dPAC_thresh = zeros(size(dPAC_perc));
    dPAC_thresh(dPAC_perc>2*std(dPAC_perc(:)) | dPAC_perc<-2*std(dPAC_perc(:)))=1;
    hold on
    contour(angleaxis,1:num_noise_levels,dPAC_thresh,1,'k','LineWidth',1) % -- plot p<0.001 clusters as overlaid black line

    % -- the same now for MI
    subplot(243)
    MI_perc = 100.*((MI - compareMI)./compareMI);
    contourf(angleaxis,1:num_noise_levels,MI_perc,40,'linestyle','none');
    cl = max(abs(get(gca,'clim')));
    set(gca,'clim',[-cl cl],...
        'ytick',[1 num_noise_levels],'yticklabel',{},...
        'xtick',0:pi/2:pi,'xticklabel',{'0', 'pi/2', 'pi'}); colorbar
    title('MI')

    MI_thresh = zeros(size(MI_perc));
    MI_thresh(MI_perc>2*std(MI_perc(:)) | MI_perc<-2*std(MI_perc(:)))=1;
    hold on
    contour(angleaxis,1:num_noise_levels,MI_thresh,1,'k','LineWidth',1) % -- plot p<0.001 clusters as overlaid black line

    % -- the same now for PLV
    subplot(244)
    PLV_perc = 100.*((PLV - comparePLV)./comparePLV);
    contourf(angleaxis,1:num_noise_levels,PLV_perc,40,'linestyle','none');
    cl = max(abs(get(gca,'clim')));
    set(gca,'clim',[-cl cl],...
        'ytick',[1 num_noise_levels],'yticklabel',{},...
        'xtick',0:pi/2:pi,'xticklabel',{'0', 'pi/2', 'pi'}); colorbar
    title('PLV')
    
    PLV_thresh = zeros(size(PLV_perc));
    PLV_thresh(PLV_perc>2*std(PLV_perc(:)) | PLV_perc<-2*std(PLV_perc(:)))=1;
    hold on
    contour(angleaxis,1:num_noise_levels,PLV_thresh,1,'k','LineWidth',1) % -- plot p<0.001 clusters as overlaid black line

    if plot_like_paper, 
        addpath(genpath('Z:\PhD\dPAC\plot_like_paper')); % -- change path accordingly
        colormap(othercolor('BuDRd_18'))
    end
    
%%
    % -- plot PACz
    if plot_like_paper, figure('position',[600 300 1200 400]); end % -- new figure needed to have different colormaps

    subplot(245) % -- PACz
    contourf(angleaxis,1:num_noise_levels,PACz,40,'linestyle','none');
    set(gca,'clim',[-4 4],...
        'ytick',[1 num_noise_levels],'yticklabel',{'min','max'},...
        'xtick',0:pi/2:pi,'xticklabel',{'0', 'pi/2', 'pi'}); colorbar
    ylabel('White noise')
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
    
    subplot(246) % -- dPACz
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
    
    subplot(247) % -- MIz
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
    
    
    subplot(248) % -- PLVz
    contourf(angleaxis,1:num_noise_levels,PLVz,40,'linestyle','none');
    set(gca,'clim',[-4 4],...
        'ytick',[1 num_noise_levels],'yticklabel',{},...
        'xtick',0:pi/2:pi,'xticklabel',{'0', 'pi/2', 'pi'}); colorbar
    title('PLVz')

    PLVz_thres_001 = squeeze(PLVz);
    PLVz_thres_001(PLVz_thres_001<threshold_001 & PLVz_thres_001>(-1*threshold_001))=0; % -- binary threshold: set all values with p>0.001 to zero
    PLVz_thres_001(PLVz_thres_001~=0)=1; % -- everything else to 1
    PLVz_thres_05 = squeeze(PLVz);
    PLVz_thres_05(PLVz_thres_05<threshold_05 & PLVz_thres_05>(-1*threshold_05))=0; % -- binary threshold: set all values with p>0.05 to zero
    PLVz_thres_05(PLVz_thres_05~=0)=1; % -- everything else to 1

    % -- cluster size thresholding (image processing toolbox required)
    if exist('bwlabel','file')
        [r,c]=bwlabel(PLVz_thres_001,4); % -- function that searches for 2D clusters
        for ci=1:c
            if sum(any(r==ci,1))<5 || sum(any(r==ci,2))<5 % -- arbitrary constraint: at least 5 contiguous points
                PLVz_thres_001(r==ci)=0;
            end
        end
        hold on
        contour(angleaxis,1:num_noise_levels,PLVz_thres_001,1,'k','LineWidth',1) % -- plot p<0.001 clusters as overlaid black line
        [r,c]=bwlabel(PLVz_thres_05,4); % -- function that searches for 2D clusters
        for ci=1:c
            if sum(any(r==ci,1))<5 || sum(any(r==ci,2))<5 % -- arbitrary constraint: at least 5 contiguous points
                PLVz_thres_05(r==ci)=0;
            end
        end
        hold on
        contour(angleaxis,1:num_noise_levels,PLVz_thres_05,1,'r','LineWidth',1) % -- plot p<0.05 clusters as overlaid red line
    end
    
    if plot_like_paper,
        oldmap = colormap;
        newmap = diverging_map(0:1/length(oldmap):1,[0 0 1],[0 1 0]);
        colormap(newmap)
    end

    %% Plot noise manipulation illustration (Figure 5)

    figure('position',[700 100 400 600]);

    gamma = ((real(theta)+6)  .* exp( 1i * 2 * pi * fpow * t )); % -- gamma is phase-modulated by theta
    % -- construct noise
    if strcmp(type_of_noise,'white')
        noise = 20 *randn(1,length(t));
    elseif strcmp(type_of_noise,'pink')
        noise = 20 *pinknoise(length(t));
    end
    % -- add all signals together
    thetagamma = theta+gamma;
    thetagammanoise = real(thetagamma)+noise;
    % -- filter in gamma and theta band
    thetafilt = filtfilt(theta_filterweights,1,thetagammanoise);
    gammafilt = filtfilt(gamma_filterweights,1,thetagammanoise);
    
    % -- extract power and phase
    thetaphase = angle(hilbert(thetafilt(1001:end-1000)));
    gammapower = abs(hilbert(gammafilt(1001:end-1000)));

    subplot(411)
    plot(t,real(gamma),'r'); hold on
    plot(t,abs(gamma),'r'); 
    plot(t,real(theta),'k'); box off
    set(gca,'xlim',[t(1) t(1000)],'xticklabel',{},'ylim',[-30 30])
    title('Pure theta-gamma coupling');
    
    subplot(412);
    plot(t,thetagammanoise,'k'); box off
    set(gca,'xlim',[t(3000) t(4000)],'xticklabel',{})
    title('Noise plus theta plus gamma')

    subplot(413);
    plot(t(1001:end-1000),(real(hilbert(thetafilt(1001:end-1000)))),'k'); hold on
    plot(t(1001:end-1000),(abs(hilbert(gammafilt(1001:end-1000)))),'r');
    plot(t(1001:end-1000),(real(hilbert(gammafilt(1001:end-1000)))),'r'); box off
    set(gca,'xlim',[t(3000) t(4000)],'ylim',[-40 40])
    title('Filtered signals')
    xlabel('Time (s)')
    ylabel('Amplitude')

end

%% Line plots in case of no coupling (Not in paper) 

if ~coupling,
    figure('position',[500 200 800 400])
    
    % -- plot contrast between bias > no bias (% sign. change)
    subplot(241)
    plot(1:num_noise_levels,(100.*((PAC - PAC_nobias_nonoise)./PAC_nobias_nonoise)),'k','linewidth',1);
    yl = max(get(gca,'ylim'));
    set(gca,'ylim',[-yl/5 yl],'xlim',[1 num_noise_levels]);
    ylabel('Coupling value (diff.)')
    title('PAC')
    box off
    
    subplot(242)
    plot(1:num_noise_levels,(100.*((dPAC - dPAC_nobias_nonoise)./dPAC_nobias_nonoise)),'k','linewidth',1);
    set(gca,'ylim',[-yl/5 yl],'xlim',[1 num_noise_levels]);
    title('dPAC')
    xlabel('Noise level');
    box off
    
    subplot(243)
    plot(1:num_noise_levels,(100.*((MI - MI_nobias_nonoise)./MI_nobias_nonoise)),'k','linewidth',1);
    yl = max(get(gca,'ylim'));
    set(gca,'ylim',[-yl/5 yl],'xlim',[1 num_noise_levels]);
    title('MI')
    box off
    
    subplot(244)
    plot(1:num_noise_levels,(100.*((PLV - PLV_nobias_nonoise)./PLV_nobias_nonoise)),'k','linewidth',1);
    yl = max(get(gca,'ylim'));
    set(gca,'ylim',[-yl/5 yl],'xlim',[1 num_noise_levels]);
    title('PLV')
    box off
    
    % -- plot z-values of permutation test
    threshold_001 = icdf('normal',1-0.001/2,0,1); % -- z-value corresponding to 99.9% or p = 0.001 of z-distribution
    threshold_05 = icdf('normal',1-0.05/2,0,1); % -- z-value corresponding to 95% or p = 0.05 of z-distribution
    
    subplot(245)
    plot(1:num_noise_levels,PACz,'k','linewidth',1); hold on
    plot([1 num_noise_levels],[threshold_001 threshold_001],'--k'); % -- 3.3 is the z-value threshold corresponding to p=0.001
    plot([1 num_noise_levels],[-threshold_001 -threshold_001],'--k');
    plot([1 num_noise_levels],[threshold_05 threshold_05],'--r'); % -- 1.96 is the z-value threshold corresponding to p=0.05
    plot([1 num_noise_levels],[-threshold_05 -threshold_05],'--r');
    set(gca,'xlim',[1 num_noise_levels]);
    ylabel('Z-value')
    title('PACz')
    box off
    
    subplot(246)
    plot(1:num_noise_levels,dPACz,'k','linewidth',1); hold on
    plot([1 num_noise_levels],[threshold_001 threshold_001],'--k'); % -- 3.3 is the z-value threshold corresponding to p=0.001
    plot([1 num_noise_levels],[-threshold_001 -threshold_001],'--k');
    plot([1 num_noise_levels],[threshold_05 threshold_05],'--r'); % -- 1.96 is the z-value threshold corresponding to p=0.05
    plot([1 num_noise_levels],[-threshold_05 -threshold_05],'--r');
    set(gca,'xlim',[1 num_noise_levels]);
    xlabel('Noise level');
    title('dPACz')
    box off
    
    subplot(247)
    h1=plot(1:num_noise_levels,MIz,'k','linewidth',1); hold on
    h2=plot([1 num_noise_levels],[threshold_001 threshold_001],'--k'); % -- 3.3 is the z-value threshold corresponding to p=0.001
    h3=plot([1 num_noise_levels],[threshold_05 threshold_05],'--r'); % -- 1.96 is the z-value threshold corresponding to p=0.05
    hasbehavior(h1,'legend',false); hasbehavior(h2,'legend',false);  hasbehavior(h3,'legend',false);
    plot([1 num_noise_levels],[-threshold_001 -threshold_001],'--k');
    plot([1 num_noise_levels],[-threshold_05 -threshold_05],'--r');
    set(gca,'xlim',[1 num_noise_levels]);
    %legend('p = 0.001','p = 0.05')
    title('MIz')
    box off
    
    subplot(248)
    h1=plot(1:num_noise_levels,PLVz,'k','linewidth',1); hold on
    h2=plot([1 num_noise_levels],[threshold_001 threshold_001],'--k'); % -- 3.3 is the z-value threshold corresponding to p=0.001
    h3=plot([1 num_noise_levels],[threshold_05 threshold_05],'--r'); % -- 1.96 is the z-value threshold corresponding to p=0.05
    hasbehavior(h1,'legend',false); hasbehavior(h2,'legend',false);  hasbehavior(h3,'legend',false);
    plot([1 num_noise_levels],[-threshold_001 -threshold_001],'--k');
    plot([1 num_noise_levels],[-threshold_05 -threshold_05],'--r');
    set(gca,'xlim',[1 num_noise_levels]);
    legend('p = 0.001','p = 0.05')
    title('PLVz')
    box off
    
end
