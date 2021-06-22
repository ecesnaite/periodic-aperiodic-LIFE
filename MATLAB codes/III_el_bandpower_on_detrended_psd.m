% (C) Elena Cesnaite, email: e.cesnaite@gmail.com, page: https://www.researchgate.net/profile/Elena-Cesnaite

% This code was created to analyze data described in 'Alterations in rhythmic and non-rhythmic resting-state
% EEG activity and their link to cognition in older age' paper.
% The code uses PSD as an input and aperiodic information (derived from FOOOF) is used to detrend
% it. Peak parameters are then estimated on detrended PSD.
% Last updated 22.06.2021

%% Store aperiodic info for those participants whose r^2 of FOOOF model fit is higher or equal to 0.8
dataDir = '' % FOOOF output directory

files = dir([dataDir,'L*'])
match = 0;
load('') % chan_all variable. Load channel information for each participant in case of missing channels
saveDir = ''
ap = [] % aperiodic component

fullch = chan_all(1).chan_names; %full channel setup

for p = 1:length(files)
    ap(p).ID = files(p).name(1:end-4);
    data = load(files(p).name, '-mat'); %1 - background params (offset and exponent (slope)); 2 - peak params (central freq, ampl, bandwi); 3 - r_squared; 5 - gaussian_params
    ch = 1 %number of channels might differ
    offset = [];
    
    while ch <= length(data.results)
        aparams = data.results{1,ch}{1};%
        l =length(offset)+1;
        
        if data.results{1,ch}{3} >= 0.8 %use data if r^2 is higher or equal to 0.8
            offset(l) = aparams(1);
            slope(l) = aparams(2);
        else
            offset(l) = NaN;
            slope(l) = NaN;
        end
        ch = ch + 2; %every second input is detrended psd;
        clearvars aparams
    end
    
    if l < 31 % if there are less than 31 channels
        indx = ismember({chan_all.ID}, files(p).name(1:end-4));
        missch =chan_all(indx).chan_names;
        
        [name,pl] = setdiff(fullch, missch) % find a missing channel
        
        if length(pl) == 1
            offset = [offset(1:pl-1), NaN, offset(pl:end)];
            slope = [slope(1:pl-1), NaN, slope(pl:end)];
        elseif length(pl) == 2
            offset = [offset(1:pl(1)-1), NaN, offset(pl(1):end)];
            offset = [offset(1:pl(2)-1), NaN, offset(pl(2):end)];
            slope = [slope(1:pl(1)-1), NaN, slope(pl(1):end)];
            slope = [slope(1:pl(2)-1), NaN, slope(pl(2):end)];
        elseif length(pl)> 2
            error("too many missing channels") % no dataset had more than 2 missing channels
        end
    end
    
    offset_mat(p,:) = offset;
    slope_mat(p,:) = slope;
    
    ap(p).offset = offset;
    ap(p).slope = slope;
    
    indx = ismember({chan_all.ID}, files(p).name(1:end-4));
    ap(p).chan = chan_all(indx).chan_names;
    clearvars -except ap p dataDir chan_all files saveDir saveDirSpec offset_mat slope_mat num_low fullch
    
end

%% Calculate alpha band power on detrended PSD

psd_all = load('') % load original (not detrended) PSD
psd_all.psd = psd_all.psd(all(~cellfun(@isempty,struct2cell(psd_all.psd)))); % find empty rows if any
saveDir = ''
f_orig = [0:0.25:45]; %original frequency sampling
fs = [2:0.25:40]; %fs after FOOOF
offset_id = {ap.ID}';

for p = 1:length(psd_all.psd)
    id = psd_all.psd(p).ID;
    original = psd_all.psd(strcmp({psd_all.psd.ID},id));
    offset1 = offset_mat(strcmp(offset_id, id),:);
    
    if isempty(offset1)
        clearvars offset1 original id
        continue
    end
    
    slope1 = slope_mat(strcmp(offset_id, id),:);
    
    for ch = 1:size(psd_all.psd(p).spect,1)
        if ~isnan(offset1(ch)) % if it's nan, it's below 0.8
            logL = offset1(ch) - slope1(ch)*log10(fs);
            L = 10.^logL;
            
            % subtract the slope from the original psd
            detrend = original.spect(ch,f_orig>=2 & f_orig<=40) - L;
            
            % search for peaks on detrended PSD
            [pow,ifreq,width,~,bounds] = findpeaks_adjusted_new(detrend,fs,'MinPeakProminence',0.05,'Annotate','extents');% min peak for power threshold
            any_alp = (ifreq >= 7 & ifreq <= 13); % look for alpha peak
            
            if any(any_alp)
                if sum(any_alp)==1 %if one peak si detected
                    freq = ifreq(any_alp);
                    range = bounds(any_alp,:);
                    
                    if diff(range) > 6      %whenever the width of the peak is more than 6Hz, set it to 3
                        range = [freq-3, freq+3];
                    end
                    
                    alpha_broad_power = sum(detrend(fs>=round(range(1),2) & fs<=round(range(2),2)))*0.25;
                    theta_range = [range(1) - 3.25 , range(1)-0.25];% theta is 3 Hz prior to alpha begining
                    
                    if theta_range(1)<=1 % theta should not be estimated prior to 1Hz
                        theta_range = [1 theta_range(2)]
                    end
                    
                    theta_broad_power = sum(original.spect(ch,f_orig>=round(theta_range(1),2) & f_orig<=round(theta_range(2),2)))*0.25; % on the original PSD
                    alpha_pow(ch) = double(alpha_broad_power);
                    theta_pow(ch) = double(theta_broad_power);
                    iaf(ch) = freq;
                    
                elseif sum(any_alp)>1 % if there are more than a single alpha peak
                    
                    [val alp_loc] = max(pow(any_alp)); %out of two alpha peaks which one is higher in amplitude
                    loc = find(any_alp); % find two alpha peak locations out of every peak
                    alp_loc = loc(alp_loc); %higher peak location out of every peak
                    any_alp(:) = 0; % set everything to 0 as if no alpha peak was found
                    any_alp(alp_loc) = 1; %set the maxima location of alpha peak to 1
                    
                    freq = ifreq(any_alp);
                    range = bounds(any_alp,:);
                    
                    if diff(range) > 6      %whenever the width of the peak is more than 6Hz, set it to 3
                        range = [freq-3, freq+3];
                    end
                    
                    alpha_broad_power = sum(detrend(fs>=round(range(1),2) & fs<=round(range(2),2)))*0.25;
                    theta_range = [range(1) - 3.25 , range(1)-0.25];% theta is 3 Hz prior to alpha begin
                    
                    if theta_range(1)<=1 % theta should not be estimated prior to 1Hz
                        theta_range = [1 theta_range(2)]
                    end
                    
                    theta_broad_power = sum(original.spect(ch,f_orig>=round(theta_range(1),2) & f_orig<=round(theta_range(2),2)))*0.25;
                    alpha_pow(ch) = double(alpha_broad_power);
                    theta_pow(ch) = double(theta_broad_power);
                    iaf(ch) = freq;
                end
            else
                
                alpha_pow(ch) =NaN;
                iaf(ch) = NaN;
                theta_pow(ch) = NaN;
                clearvars val alp_loc loc any_alp freq range alpha_broad_power pow ifreq width detrend L logL bounds theta_broad_power
            end
        else
            
            alpha_pow(ch) =NaN;
            iaf(ch) = NaN;
            theta_pow(ch) = NaN;
            clearvars val alp_loc loc any_alp freq range alpha_broad_power pow ifreq width detrend L logL theta_broad_power
        end
        clearvars val alp_loc loc any_alp freq range alpha_broad_power pow ifreq width detrend L logL bounds theta_broad_power
    end
    
    alpha_params(p).ID = id
    alpha_params(p).power = alpha_pow
    alpha_params(p).freq = iaf
    alpha_params(p).theta = theta_pow
    
    indx = ismember({chan_all.ID}, id);
    alpha_params(p).chan = chan_all(indx).chan_names
    
    clearvars id logi pl offset1 slope1 original alpha_pow iaf theta_pow
    
end

alpha_params = alpha_params(all(~cellfun(@isempty,struct2cell(alpha_params))));

for i = 1:length(alpha_params)
    
    chind = ismember(all_chan, alpha_params(i).chan)%find missing channels
    
    if sum(~chind) == 1
        alpha_broad.broadpow(i,:)=[alpha_params(i).power(1:find(~chind)-1) NaN alpha_params(i).power(find(~chind):end)];
    elseif sum(~chind) == 2
        nval = find(~chind)
        new_mat =[alpha_params(i).power(1:nval(1)-1) NaN alpha_params(i).power(nval(1):end)]
        new_mat = [new_mat(1:nval(2)-1) NaN new_mat(nval(2):end)]
        alpha_broad.broadpow(i,:)=new_mat;
        
    else
        error('too many missing channels')
    end
    
    clearvars new_mat chind nval
end
save([saveDir, 'alpha_params_findpeaks_detrend_2_40_7_13_threshold_05_a_stages_5min'], 'alpha_params')

% Remove negative values that appeared due to a small alpha peak and its estimation under the curve

[row channel] = find(alpha_broad.broadpow<0)
for n = 1:length(row)
    if alpha_broad.broadpow(row(n), channel(n))>0
        error('its more than 0')
    end
    alpha_broad.broadpow(row(n), channel(n))=NaN %check if it's really -
    alpha_params(row(n)).power(channel(n)) = NaN 
    alpha_params(row(n)).freq(channel(n)) = NaN
    alpha_params(row(n)).theta(channel(n)) = NaN
    
end

save([saveDir, 'alpha_params_findpeaks_detrend_2_40_7_13_threshold_05_a_stages_5min'], 'alpha_params')
save([saveDir, 'alpha_broad_pow_mtx_a_stages_5min'], 'alpha_broad')



