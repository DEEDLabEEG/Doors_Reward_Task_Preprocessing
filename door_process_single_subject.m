% ********************************************************************** %
% Preprocessing Script for Doors Task EEG Data [Script 2]
% Authors: Armen Bagdasarov & Kenneth Roberts
% Institution: Duke University
% Year created: 2021
% Year last modified: 2021
% ********************************************************************** %

% Run looping script (door_loop_over_subjects.m) to run this script, which is a function

function summary_info = door_process_single_subject(varargin)

global proj % Declare variables as global (variables that you can access in other functions)

% ********************************************************************** %

%% Reset the random number generator to try to make results replicable 
% This produces somewhat consistent results for clean_rawdata functions and ICA
rng('default');

% ********************************************************************** %

%% Import data
% Raw data is in .mff format
mff_filename = fullfile(proj.data_location, proj.mff_filenames{proj.currentSub}); % Get file name
EEG = pop_mffimport({mff_filename},{'code'}); % Import .mff
summary_info.currentId = {proj.currentId}; % Save subject ID in summary info

% ********************************************************************** %

%% Remove outer ring of electrodes
% The outer ring is often noisy in high-density nets, so let's remove these channels

% Save variable with all channel locations (129 channels)
full_chan_locs = EEG.chanlocs; 

% Outer layer of channels to be removed (24 of them)
outer_chans = {'E17' 'E38' 'E43' 'E44' 'E48' 'E49' 'E113' 'E114' 'E119' 'E120' 'E121'...
    'E125' 'E126' 'E127' 'E128' 'E56' 'E63' 'E68' 'E73' 'E81' 'E88' 'E94' 'E99' 'E107'};

% Remove outerlayer_channels
EEG = pop_select(EEG, 'nochannel', outer_chans);

% ********************************************************************** %

%% Downsample from 1000 to 250 Hz
EEG = pop_resample(EEG, 250);

% ********************************************************************** %

%% Remove segments without events (i.e., breaks)

% Remove the first 20 seconds of data first (contains SESS, CELL, etc.)
EEG = pop_select(EEG, 'notime',[0 20]);

% Delete segments of data between event codes if the length of the segment is greater than 4 seconds (4000 ms)
% Save 1000 ms of time before the first event code and 1000 ms after the end of the last event code
% This keeps the practice events, which is fine
EEG  = pop_erplabDeleteTimeSegments(EEG , 'displayEEG',  0, 'endEventcodeBufferMS',  ...
    1000, 'ignoreUseType', 'Ignore', 'startEventcodeBufferMS',  1000, 'timeThresholdMS',  4000 );

% ********************************************************************** %

%% Filter & remove line noise

% 0.1 Hz high-pass and 30 Hz low-pass Butterworth filter
EEG  = pop_basicfilter(EEG, 1:EEG.nbchan, 'Boundary', 'boundary', 'Cutoff', [0.1 30], ...
    'Design', 'butter', 'Filter', 'bandpass', 'Order',  4, 'RemoveDC', 'off');

% CleanLine to remove 60 Hz line noise 
% Because even though we low pass filtered at 30 Hz, there is still some line noise creeping in
% This will help ICA later, so let's do it
EEG = pop_cleanline(EEG, 'bandwidth',2,'chanlist',1:EEG.nbchan ,'computepower',1,...
    'linefreqs',60,'newversion',0,'normSpectrum',0,'p',0.01,'pad',2,...
    'plotfigures',0,'scanforlines',1,'sigtype','Channels','taperbandwidth',2,...
    'tau',50,'verb',1,'winsize',4,'winstep',1);

% ********************************************************************** %

%% Re-reference to the mastoids
%  E57 = Left mastoid
%  E100 = Right mastoid
chans_to_ref = find(ismember({EEG.chanlocs(:).labels}, {'E57', 'E100'}));
EEG = pop_reref(EEG,chans_to_ref, 'keepref','on'); % Keep reference (Cz/E129)

% Save variable with reduced (104) channel locations after re-referencing but 
% before rejecting bad channels (will be needed later when interpolating bad channels)
reduced_chan_locs = EEG.chanlocs;

% ********************************************************************** %

%% Reject bad channels 

% These are all default settings
% Only rejecting bad channels, so everything else (e.g., ASR) is turned off
all_eeg_chanlabels = {EEG.chanlocs(:).labels}; % Store list before rejection
EEG = pop_clean_rawdata(EEG, 'FlatlineCriterion',5,'ChannelCriterion',0.8,...
    'LineNoiseCriterion',4,'Highpass','off','BurstCriterion','off',...
    'WindowCriterion','off','BurstRejection','off','Distance','Euclidian',...
    'MaxMem', 48); % MaxMem set to 48gb for reproducibility 
kept_eeg_chanlabels = {EEG.chanlocs(:).labels}; % Store list after rejection 
bad_chans = setdiff(all_eeg_chanlabels, kept_eeg_chanlabels);

% Save which channels were bad in summary info
if isempty(bad_chans) % If no bad chans...
    summary_info.bad_chans = {[]}; % ...then leave blank
else
    summary_info.bad_chans = {strjoin(bad_chans)};

    % Plot bad channels to identify whether there are clusters of bad channels
    bad_chan_ind = find(ismember({reduced_chan_locs(:).labels}, bad_chans));
    figure; topoplot(bad_chan_ind, reduced_chan_locs, 'style', 'blank', ...
        'emarker', {'.','k',[],10}, 'electrodes', 'ptslabels');

    % Save bad channels plot
    set(gcf, 'Units', 'Inches', 'Position', [0, 0, 10, 10], 'PaperUnits', 'Inches', 'PaperSize', [10, 10])
    bad_chan_plot_path = 'INSERT_PATH'; 
    bad_chan_plot_name = [proj.currentId '_door_bad_chans_plot'];
    saveas(gca, fullfile(bad_chan_plot_path, bad_chan_plot_name), 'png');
    close(gcf);
end

% Save the number of bad channels in summary info
summary_info.n_bad_chans = length(bad_chans);

% ********************************************************************** %

%% Save interim file

% Save file before ASR and ICA
interim_path = 'INSERT_PATH'; 
interim_name = [proj.currentId '_door_interim_1'];
pop_saveset(EEG, fullfile(interim_path, interim_name));

% ********************************************************************** %

%% Remove large artifacts with ASR

% Artifact Subspace Reconstruction (ASR) + additional removal of bad data periods

% First, let's save our data before ASR
% We will run ICA later on the EEG data after ASR
% But apply the ICA fields to the full EEG data
EEG_no_rej = EEG;

% ASR and ICA work best with data filtered at 1 Hz high pass, so let's high pass filter again temporarily
EEG  = pop_basicfilter(EEG, 1:EEG.nbchan, 'Boundary', 'boundary', 'Cutoff', 1, ...
    'Design', 'butter', 'Filter', 'highpass', 'Order',  4, 'RemoveDC', 'off');

% ASR
% These are all default settings
% Most importantly the burst criterion is set to 20 and burst rejection is set to on 
EEG = pop_clean_rawdata(EEG, 'FlatlineCriterion','off','ChannelCriterion','off',...
    'LineNoiseCriterion','off','Highpass','off','BurstCriterion',20,...
    'WindowCriterion',0.25,'BurstRejection','on','Distance','Euclidian',...
    'WindowCriterionTolerances',[-Inf 7], 'MaxMem', 48); % MaxMem set to 48gb for reproducibility 

% Save how many seconds of data is left after ASR in summary info
% This will be important later for excluding participants 
% For example, if ICA was run only on 30 seconds of data because ASR cut
% out the rest, we probably should get rid of this participant (i.e., their
% data was probably very noisy)
summary_info.post_ASR_data_length = EEG.xmax;

% ********************************************************************** %

%% ICA

% Extended infomax ICA with PCA dimension reduction
% PCA dimension reduction is necessary because we have a lot of channels and a relatively short amount of data
EEG = pop_runica(EEG, 'icatype', 'runica', 'extended', 1, 'pca', 50);

% Save ICA plot
ica_plot_path = 'INSERT_PATH'; 
ica_plot_name = [proj.currentId '_door_ica_plot'];
pop_topoplot(EEG, 0, [1:50], 'Independent Components', 0, 'electrodes','off');
set(gcf, 'Units', 'Inches', 'Position', [0, 0, 10, 10], 'PaperUnits', 'Inches', 'PaperSize', [10, 10])
saveas(gca, fullfile(ica_plot_path, ica_plot_name), 'png');
close(gcf);

% ********************************************************************** %

%% Save interim file

% Save file after ASR and ICA decomposition is complete but before ICs are actually removed
interim_path = 'INSERT_PATH'; 
interim_name = [proj.currentId '_door_interim_2'];
pop_saveset(EEG, fullfile(interim_path, interim_name));

% ********************************************************************** %

%% Select IC components related to eye artifact only

% Automatic classification with ICLabel
EEG = pop_iclabel(EEG, 'default');

% Flag components with >= 70% of being eye
EEG = pop_icflag(EEG, [NaN NaN; NaN NaN; 0.7 1; NaN NaN; NaN NaN; NaN NaN; NaN NaN]);

% Select components with >= 70% of being eye
eye_prob = EEG.etc.ic_classification.ICLabel.classifications(:,3);
eye_rej = find(eye_prob >= .70);
eye_rej = [eye_rej];
eye_rej = eye_rej';

% Save retained variance post-ICA in summary info
[projected, pvar] = compvar(EEG.data, {EEG.icasphere, EEG.icaweights}, EEG.icawinv, eye_rej);
summary_info.var_retained = 100-pvar;

% Plot only the removed components
ica_rej_plot_path = 'INSERT_PATH'; 
ica_rej_plot_name = [proj.currentId '_door_removed_ics'];
figure % This line is necessary for if there is only 1 component to plot
pop_topoplot(EEG, 0, eye_rej, 'Independent Components', 0, 'electrodes','off');
set(gcf, 'Units', 'Inches', 'Position', [0, 0, 10, 10], 'PaperUnits', 'Inches', 'PaperSize', [10, 10])
saveas(gca, fullfile(ica_rej_plot_path, ica_rej_plot_name), 'png');
close(gcf); close(gcf);

% ********************************************************************** %

%% Copy over EEG ICA fields to EEG_no_rej and remove components with >= 70% of being eye
% Basically, we are back-projecting the ICA information from the ASR-reduced data to the full data

EEG_no_rej.icawinv = EEG.icawinv;
EEG_no_rej.icasphere = EEG.icasphere;
EEG_no_rej.icaweights = EEG.icaweights;
EEG_no_rej.icachansind = EEG.icachansind;

EEG = EEG_no_rej; % Set EEG to the one with full data length, so it's not longer the ASR-reduced one

% Remove components with >= 70% of being eye
EEG = pop_subcomp(EEG, eye_rej, 0);

% Save which components were removed in summary info
if isempty(eye_rej) % If no ICs removed...
    summary_info.ics_removed = {[]}; % ...then leave blank
else
    % Save which components were removed in summary info
    summary_info.ics_removed = {num2str(eye_rej)};
end

% Save the number of components removed in summary info
summary_info.n_ics_removed = length(eye_rej);

% ********************************************************************** %

%% Interpolate removed bad channels
EEG = pop_interp(EEG, reduced_chan_locs, 'spherical');

% ********************************************************************** %

%% Save file after IC removal & chan interpolation but before anything else
interim_path = 'INSERT_PATH'; 
interim_name = [proj.currentId '_door_interim_3'];
pop_saveset(EEG, fullfile(interim_path, interim_name));

% ********************************************************************** %

%% Plot channel spectra

% NOTE: This is just to check for extreme line noise after all of the steps
% above. The spectra will not look fantastic because artifact rejection was
% not done yet. However, plotting spectra before epoching is preferred. So,
% let's do it now.

% Plot channel spectra 
spectra_time = EEG.xmax * 1000;
figure; pop_spectopo(EEG, 1, [0  spectra_time], 'EEG' , 'freqrange',[2 80],'electrodes','off');

% Save channel spectra plot
spectra_plot_path = 'INSERT_PATH'; 
spectra_plot_name = [proj.currentId '_door_spectra_plot'];
saveas(gca, fullfile(spectra_plot_path, spectra_plot_name), 'png');
close(gcf);

% ********************************************************************** %

%% Segmentation

% Select only the events that matter
EEG = pop_selectevent( EEG, 'type',{'Rew','Pun','Neu'},'deleteevents','on');

% Load event list
el_filename = ['INSERT_PATH' proj.currentId '_door_event_list.txt'];
EEG  = pop_editeventlist( EEG , 'BoundaryNumeric', { -99}, 'BoundaryString', { 'boundary' }, ...
    'ExportEL', el_filename, 'List', 'INSERT_PATH\door_equation_list.txt', ...
    'SendEL2', 'EEG&Text', 'UpdateEEG', 'code', 'Warning', 'off' );

% Upload bins file
el_import_filename =  fullfile('', el_filename);
EEG  = pop_binlister( EEG , 'BDF', 'INSERT_PATH\door_bins.txt', 'ImportEL', ...
  el_import_filename, 'IndexEL',  1, 'SendEL2', 'EEG', 'Voutput', 'EEG' ); 

% Extract bin-based epochs and apply baseline correction
EEG = pop_epochbin(EEG , [-500  800],  [-200 0]); 

% ********************************************************************** %

%% Artifact rejection using TBT plugin

EEG_orig_erp = EEG; % Save the EEG struct

% Max-min amplitude difference
% Basically, the same as ERPLAB's moving window peak-to-peak threshold
% Peak-to-peak amplitude exceeding 100 uV within 200 ms windows sliding by 20 ms
EEG = pop_eegmaxmin(EEG, 1:EEG.nbchan,[-500  796], 100, 200, 20, 0);

% Abnormal values - simple voltage thresholding
% -150/+150 uV
EEG = pop_eegthresh(EEG, 1, 1:EEG.nbchan, -150 , 150 , -0.5 , 0.796, 1, 0); 

% Improbable data - based on joint probability
% SD = 3 for both local and global thresholds 
EEG = pop_jointprob(EEG, 1, 1:EEG.nbchan, 3, 3, 1, 0, 0);

% Reject based on the epochs selected above
[EEG badlist ]= pop_TBT(EEG, EEG.reject.rejmaxminE | EEG.reject.rejthreshE | EEG.reject.rejjpE, 10,1,0); % Do chan/epoch interpolation on ALL types of artifact rejection at once
    % Criteria must be met in at least 10 channels for the epoch to be rejected
    % We don't want to remove channels, so criteria must be met in all channels for the channel to be removed (unlikely to happen so we set at 100%)

% Now for a bit of a hack: We want to copy the interpolated epochs from pop_TBT back into the original struct
    
    % Use eeg_epochformat to get epoch fields to timelocking event
    [epoch_arr, fn] = eeg_epochformat(EEG.epoch, 'array');
    fn_ind = find(strcmp(fn, 'eventbepoch'));
    eventbepoch_tbt = [ epoch_arr{:, fn_ind}];
    
    % Push good events back into the struct
    EEG_orig_erp.data(:,:,eventbepoch_tbt) = EEG.data;
    
    % Set the EEGLAB flags in EEG.reject to reflect epochs flagged by pop_TBT
    EEG_orig_erp.reject.rejmanual = true(1, size(EEG_orig_erp.data, 3));
    EEG_orig_erp.reject.rejmanual(eventbepoch_tbt) = false;
    
    % Sync from EEGLAB to ERPLAB - Mark the bad epochs in an ERPLAB-compatible way so ERPLAB knows which epochs are bad
    EEG_orig_erp = pop_syncroartifacts(EEG_orig_erp, 'Direction', 'eeglab2erplab');
    
    % Copy back into EEG
    EEG = EEG_orig_erp;

% Save how many epochs are kept 
summary_info.n_epochs_kept = sum(~EEG.reject.rejmanual);

%% Compute averaged ERPs

ERP = pop_averager(EEG , 'Criterion', 'good', 'DQ_flag', 1, 'ExcludeBoundary', 'on', 'SEM', 'on' );

% And save them
erp_name = [proj.currentId '_door_averaged_erps'];
erp_file_name = [erp_name '.erp'];
erp_path = 'INSERT_PATH';
ERP = pop_savemyerp(ERP, 'erpname', erp_name, 'filename', erp_file_name, 'filepath', erp_path, 'Warning', 'off');

%% Save trials rejected & accepted per bin

summary_info.trials_rej_b1_rew = ERP.ntrials.rejected(1);
summary_info.trials_rej_b2_pun = ERP.ntrials.rejected(2);
summary_info.trials_rej_b3_neu = ERP.ntrials.rejected(3);

summary_info.trials_accept_b1_rew = ERP.ntrials.accepted(1);
summary_info.trials_accept_b2_pun = ERP.ntrials.accepted(2);
summary_info.trials_accept_b3_neu = ERP.ntrials.accepted(3);

%% Save final preprocessed files
% .set format
set_path = 'INSERT_PATH'; 
set_name = [proj.currentId '_door_complete'];
pop_saveset(EEG, fullfile(set_path, set_name));

% ****************************** THE END ******************************* %
