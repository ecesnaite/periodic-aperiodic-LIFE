
% (C) Elena Cesnaite, email: e.cesnaite@gmail.com, page: https://www.researchgate.net/profile/Elena-Cesnaite

% This code was created to analyze data described in 'Alterations in rhythmic and non-rhythmic resting-state 
% EEG activity and their link to cognition in older age' paper. The code rejects bad data segments and calculates power 
% spectral density for each channel. 
% Last updated 22.06.2021

clc; clear all; close all

%call eeglab
eeglab

aftDir ='' %A stages have been cleaned from artifacts (trimOut + ICA) and cut to 5 min
saveDir = ''

files = dir([aftDir, '*.set']);

for isb = 1:length(files)
    
    EEG = pop_loadset(files(isb).name, aftDir);
    name = extractBefore(files(isb).name, '_prep')
    
    %find events associated with bad segments   
    events = EEG.event(:,(ismember({EEG.event.type}, [{'auto_start'},{'manual_start'}])),:);
    rejData(:,1) = [events.latency];
    rejData(:,2) = [events.latency] + [events.duration];
    rejData = sortrows(rejData);
    EEG.rejData = rejData;
    
    % remove bad data segments if there were any
    if ~isempty(rejData)
        EEG_rej = pop_select(EEG, 'nopoint', [rejData(:,1) rejData(:,2)]); %reject data for ICA.
    else
        EEG_rej = EEG
    end
    
    %make sure data is 5 min long
    if length(EEG_rej.data)<EEG_rej.srate*5*60
        clearvars EEG_rej EEG rejData events name
        continue
    end
    
    %calculate PSD
        [spect freq] = el_plot_spec_freq_adjust(EEG.data',EEG.srate,45)
        spect = spect';
        psd(isb).spect = spect
        psd(isb).freq = freq
        psd(isb).ID = name
    
    clear rejData EEG EEG_rej EEG_rej_pas act pas sim_gco new_passive spect freq name%UPDATE
end

save([saveDir, 'psd_for_FOOOF_a_stages_5min'], 'psd')








