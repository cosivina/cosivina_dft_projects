%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% COSIVINA IOWA Task Simulator
%
% John Spencer
% building on in-line simulator from Sebastian Schneegans
% and IOWA-C expansion from 2022 paper
%
% origin: stripped down version of the biased competition model (spatial pathway only)
% model components:
% field_a - spatial attention field
% field_s - saccade motor field
% neuron_r - saccade reset neuron
% neuron_x - fixation neuron, excites foveal position in field_a
% neuron_g - gaze-change neuron, inhibits foveal position in field_a
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% constant identifiers %%
%%%%%%%%%%%%%%%%%%%%%%%%%%

close all; clear all; 

mode = 1; % 0 for batch mode, 1 for auto mode with visualization, 2 for multicore (also need to switch to 'parfor' in line 193 below)

%% use HPC?
HPC=0;

if HPC
    run('../../cosivina/setpath.m') % add cosivina and jsoblab files to the matlab path
    parpool('SlurmProfile1',96)
end

randn('state',sum(100*clock));
rand('state',sum(100*clock));

%% Task identifiers
DOUBLE_CUE = 1;
INVALID_CUE = 2;
NO_CUE = 3;
TONE_CUE = 4;
VALID_CUE = 5;
CDOUBLE_CUE = 6;
CINVALID_CUE = 7;
CNO_CUE = 8;
CTONE_CUE = 9;
CVALID_CUE = 10;
N_CONDITIONS = 10; %first 5 IOWA, second 5 IOWA-C

% Age identifiers
INFANT_5MO = 1;
INFANT_7MO = 2;
INFANT_10MO = 3;
N_AGES = 3;


%% experiment settings %%
%%%%%%%%%%%%%%%%%%%%%%%%%

% stimulus sizes (one field unit = 0.1 degrees of visual angle)
targetSize = 46; %30; %46; % matches mean diameter of target in experiment
targetEccentricity = 110;

% stimulus presentation timing (in milliseconds)
tCueStart = 400; %in original, tFixEnd = 300 then cue at 400;
tCueEnd = 500; %tCueStart + 100;
tTargetStart = tCueEnd + 100;

minSaccLatency = 100;
maxSaccLatency = 2600;

tMax = tTargetStart + maxSaccLatency;
tTargetEnd = tMax;
tToneBoostStart = tCueStart;
tToneBoostEnd = tCueEnd;

breakOnFirstSaccade = true; % current method for storing saccade metrics requires this to be true


%% model parameters %%
%%%%%%%%%%%%%%%%%%%%%%

% Parameters are all defined in a struct p for easier management and
% storage. For each parameter, either a single value should be specified,
% or a vector of three values for ages [5MO, 7MO, 10MO]. Additional
% modifications on these parameters (for trying out different settings in a
% batch) can be specified below.

p.cueSize = 17; % 10 matches diameter of cue in experiment
p.fixationSize = 27; % 27 matches mean diameter of fixation cue in experiment

p.spatialFieldSize = 301; % FIXED; 1 deg equals 10 field units
spatialHalfSize = (p.spatialFieldSize-1)/2;

p.tau_a = 60;
p.tau_s = 60;
p.h_a = -5; p.beta_a = 1; p.q_a = [0.55*2, 0.5*2, 0.45*2]; %was [.55, .5, .45]
p.h_s = -5; p.beta_s = 4; p.q_s = [0.6, 0.55, 0.5]; %[0.6, 0.55, 0.5];
p.h_r = -5; p.beta_r = 4; p.q_r = 0.05;
p.h_x = -5; p.beta_x = 1; p.q_x = 0.05;
p.h_g = -5; p.beta_g = 1; %p.q_g = 0.05;

p.sigma_q = 2;

% smoothing of visual input
p.c_v_exc = [8, 8, 8]; p.sigma_v_exc = 2.5;
p.c_v_inh = 0; p.sigma_v_inh = 5;
p.c_v_gi = 0;

% spatial attention field
p.c_aa_exc = [22, 27, 27]; %was [20,27,27]
p.sigma_aa_exc = 8; 
p.c_aa_inh = [20, 24, 30]; p.sigma_aa_inh = 20; %[20, 24, 30]
p.c_aa_gi = [0.045*1.95, 0.1*1.95, 0.1*1.95]; %[0.045, 0.1, 0.13] *1.9

% saccade motor field
p.c_ss_exc = [30, 40, 50]; p.sigma_ss_exc = 8;
p.c_ss_inh = 0; p.sigma_ss_inh = 8;
p.c_ss_gi = [0.75, 1.0, 1.25]; %[0.75, 1.0, 1.25]

p.c_sa = [6.75, 7.4, 9.75]; p.sigma_sa = 10; % [6.75, 7.4, 9.75]
p.c_as = [0, 0, 0]; p.sigma_as = 10;

p.c_rr = [2.5, 2.5, 2.5];
p.c_rs = [0.9, 1.2, 1.5];
p.c_sr_gi = [18, 24, 30];
p.c_ar_gi = [18, 24, 30];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% NEW -- some details needed for 2022 paper
p.c_gg = 2.5; %0.15; %0.1; %was 0
p.c_xx = 5.0; %0.15; %[3, 3, 3]; %was 0 1.5
p.c_gx_inh = 1.0; %0.4; % inhibition from g+x to both of them--keep set to 0.

p.c_ax = [5, 5, 5]; %was 5 %strength of fix node influence
p.sigma_ax_exc = [8, 8, 8]; %width of fixation node influence

%gaze change node has inverted Gaussian kernel with c_ag_exc constant and
%inhib parameters define Gaussian (c_ag_inh = neg strength; sigma_ag_exc is
%width
p.c_ag_inh = 10; %says 0 in paper!! go with 10 since that's value in final simulator
p.c_ag_exc = 2; %says 10 in paper!! go with 2 since that's value in final simulator
p.sigma_ag_exc = 8;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p.c_xr_inh = 5; % was 5
p.c_gr_inh = 5;

p.theta_saccStart = 0.95; % threshold of saccade reset neuron for saccade begin
p.theta_saccEnd = 0.05; % threshold of saccade reset neuron for saccade end
p.c_sacc = 0.0015; % size of eye movement (deg) per saccade field output unit

p.stimStrength_g_tone = 4; %5; %was4
p.stimStrength_g_base = 4; %0; %was4

p.stimStrength_x = [7, 7, 7]; %[0.1, 0.1, 0.1]; %was 10 %new--John had set to .1

p.saccadeLatencyOffset = [75, 75, 75]; % must be 3-element vector for evaluation to work properly

p.kernelWidthMultiplier = 4;


%% simulation settings %%
%%%%%%%%%%%%%%%%%%%%%%%%%

if mode == 1
    plotLatencies = false;
    plotErrorRates = false;
    printTrials = true;
    keepFullResults = false;%
    saveResults =false;%
    storeHistory = true;
    %conditions = [DOUBLE_CUE, INVALID_CUE, NO_CUE, TONE_CUE, VALID_CUE, CDOUBLE_CUE, CINVALID_CUE, CNO_CUE, CTONE_CUE, CVALID_CUE];
    conditions = [VALID_CUE, INVALID_CUE, DOUBLE_CUE, TONE_CUE, NO_CUE, CVALID_CUE, CINVALID_CUE, CDOUBLE_CUE, CTONE_CUE, CNO_CUE];
    ages = [INFANT_10MO]; 
    nRepeats = 1; %400;
    noise = true;
else
    plotLatencies = false;
    plotErrorRates = false;
    printTrials = false;
    keepFullResults = true;%
    saveResults =true;%
    resultBaseFilename = 'Sims2_';
    storeHistory = false;
    %conditions --> first 5 IOWA, second 5 IOWA-C
    conditions = [DOUBLE_CUE, INVALID_CUE, NO_CUE, TONE_CUE, VALID_CUE, CDOUBLE_CUE, CINVALID_CUE, CNO_CUE, CTONE_CUE, CVALID_CUE];

% %     %repeating conditions --> first 5 IOWA, second 5 IOWA-C, then 4
% %     %IOWA-VWM, then 4 IOWAC-VWM
% %     conditions = [DOUBLE_CUE, INVALID_CUE, NO_CUE, TONE_CUE, VALID_CUE, CDOUBLE_CUE, CINVALID_CUE, CNO_CUE, CTONE_CUE, CVALID_CUE, ...
% %         DOUBLE_CUE, INVALID_CUE, NO_CUE, VALID_CUE, CDOUBLE_CUE, CINVALID_CUE, CNO_CUE, CVALID_CUE];
    ages = [INFANT_5MO, INFANT_7MO, INFANT_10MO];
    nRepeats = 400;
    noise = true;
end


%% parameter modifications for batch trials %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

paramMods = struct('param', {}, 'ages', {}, 'factors', {});

%Sims optimized 2019
% % paramMods(end + 1) = struct('param', 'stimStrength_x', 'ages', [INFANT_5MO, INFANT_7MO, INFANT_10MO], 'factors', [50:50:250]);
% % paramMods(end + 1) = struct('param', 'c_xx', 'ages', [INFANT_5MO, INFANT_7MO, INFANT_10MO], 'factors', [0.5:0.5:1.5]);
% % paramMods(end + 1) = struct('param', 'c_ax', 'ages', [INFANT_5MO, INFANT_7MO, INFANT_10MO], 'factors', [6:0.5:8]);
% % paramMods(end + 1) = struct('param', 'cueSize', 'ages', [INFANT_5MO, INFANT_7MO, INFANT_10MO], 'factors', [19:0.5:21]);

%Sims optimized 2020
% % paramMods(end + 1) = struct('param', 'stimStrength_x', 'ages', [INFANT_5MO, INFANT_7MO, INFANT_10MO], 'factors', [1:.3:1.3]); 
% % paramMods(end + 1) = struct('param', 'c_gg', 'ages', [INFANT_5MO, INFANT_7MO, INFANT_10MO], 'factors', [1:.075:1.15]);
% % paramMods(end + 1) = struct('param', 'c_xx', 'ages', [INFANT_5MO, INFANT_7MO, INFANT_10MO], 'factors', [1:.075:1.15]);
% % paramMods(end + 1) = struct('param', 'c_gx_inh', 'ages', [INFANT_5MO, INFANT_7MO, INFANT_10MO], 'factors', [1:0.2:1.2]);
% % paramMods(end + 1) = struct('param', 'c_aa_gi', 'ages', [INFANT_5MO, INFANT_7MO, INFANT_10MO], 'factors', [1.7:.1:1.9]);
% % paramMods(end + 1) = struct('param', 'cueSize', 'ages', [INFANT_5MO, INFANT_7MO, INFANT_10MO], 'factors', [1:.2:1.6]);
% % paramMods(end + 1) = struct('param', 'q_x', 'ages', [INFANT_5MO, INFANT_7MO, INFANT_10MO], 'factors', [1:4:9]);

% % paramMods(end + 1) = struct('param', 'fixationSize', 'ages', [INFANT_5MO, INFANT_7MO, INFANT_10MO], 'factors', [0.5:0.25:1]);
% % paramMods(end + 1) = struct('param', 'stimStrength_g_tone', 'ages', [INFANT_5MO, INFANT_7MO, INFANT_10MO], 'factors', [0.75:.25:1.25]); 


%% prepare batches of trials %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nParamMods = numel(paramMods);
nFactors = zeros(1, nParamMods);
for m = 1 : nParamMods
    nFactors(m) = numel(paramMods(m).factors);
end
nRuns = prod(nFactors);
nFactorsPre = [1, cumprod(nFactors(1:end-1))];
nFactorsPost = fliplr([1, cumprod(fliplr(nFactors(2:end)))]);

modAges = cell(nRuns, nParamMods);
for m = 1 : nParamMods
    paramMods(m).factors = reshape(repmat(reshape(paramMods(m).factors, 1, []), ...
        [nFactorsPost(m), nFactorsPre(m)]), [], 1);
    modAges(:, m) = reshape(repmat(cat(2, {ages}, repmat({paramMods(m).ages}, [1, nFactors(m) - 1])), ...
        [nFactorsPost(m), nFactorsPre(m)]), [], 1);
end

agesInRuns = repmat({ages}, [nRuns, 1]);
for r = 1 : nRuns
    for m = 1 : nParamMods
        agesInRuns{r} = intersect(agesInRuns{r}, modAges{r, m});
        if isempty(agesInRuns{r})
            agesInRuns{r} = []; % necessary to ensure that loop over runs does not produce errors
        end
    end
end

ageRunPointers = nan(nRuns, N_AGES);
for r = 1 : nRuns
    for a = ages
        for j = r : -1 : 1
            if any(agesInRuns{j} == a)
                ageRunFound = true;
                for m = 1 : nParamMods
                    if any(paramMods(m).ages == a) && paramMods(m).factors(j) ~= paramMods(m).factors(r)
                        ageRunFound = false;
                        break;
                    end
                end
                if ageRunFound
                    ageRunPointers(r, a) = j;
                    break;
                end
            end
        end
    end
end


%% prepare results struct and file %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nConditions = length(conditions);
nAges = length(ages);

result = struct('conditions', {}, 'ages', {}, 'parameters', {}, 'saccStartTimes', {}, 'saccEndTimes', {}, ...
    'saccMetrics', {}, 'latenciesCorrect', {}, 'ratesCorrect', {}, 'ageRunPointers', {}, ...
    'noise', {}, 'nRepeats', {}, 'errorRatesEstimated', {});

if saveResults
    resultFilename = [resultBaseFilename, datestr(now, 'yyyy-mm-dd-THHMMSS') '.mat'];
end


%% simulation %%
%%%%%%%%%%%%%%%%
latenciesCorrect = nan(N_CONDITIONS, N_AGES);
ratesCorrect = nan(N_CONDITIONS, N_AGES);

for r = 1 : nRuns % loop over runs (different paramter modifications)
    disp(['Run ' num2str(r) ' out of ' num2str(nRuns) ', started ' datestr(now)])
    
    saccStartTimes = cell(N_AGES, 1);
    saccEndTimes = cell(N_AGES, 1);
    saccMetrics = cell(N_AGES, 1);
    
    % apply parameter modifications for this run and age group
    pr = p;
    for j = 1 : nParamMods
        if isscalar(pr.(paramMods(j).param))
            pr.(paramMods(j).param) = p.(paramMods(j).param) * paramMods(j).factors(r);
        else
            for k = 1 : numel(paramMods(j).ages)
                a = paramMods(j).ages(k);
                pr.(paramMods(j).param)(a) = p.(paramMods(j).param)(a) * paramMods(j).factors(r);
            end
        end
    end
    
    for a = agesInRuns{r} % loop over ages for this run
        
        disp(['Age ' num2str(a)])
        
        % set parameters for this age group
        paramNames = fieldnames(p);
        pa = pr;
        for j = 1 : length(paramNames)
            if ~isscalar(pr.(paramNames{j}))
                pa.(paramNames{j}) = pr.(paramNames{j})(a);
            end
        end               

        %% initialize COSIVINA
        createIOWASim_2022;
        sim.init();

        gui_speed = 3; %display every X timesteps if gui engaged
        createIOWAGUI_2022;

        if mode == 1
            gui.init();
        end

        %kernel used to create visual input
        kSize_v = min(round(pa.kernelWidthMultiplier * max(pa.sigma_v_exc, pa.sigma_v_inh)), pa.spatialFieldSize);
        g1 = exp(-0.5 * ((-kSize_v:kSize_v)-0).^2 / pa.sigma_v_exc^2);
        g1 = g1 / sum(g1);
        g2 = exp(-0.5 * ((-kSize_v:kSize_v)-0).^2 / pa.sigma_v_inh^2);
        g2 = g2 / sum(g2);
        kernel_v = pa.c_v_exc * g1 - pa.c_v_inh * g2;

        saccStartTimes{a} = nan(N_CONDITIONS, nRepeats);
        saccEndTimes{a} = nan(N_CONDITIONS, nRepeats);
        saccMetrics{a} = nan(N_CONDITIONS, nRepeats);
        
        
        %% simulation %%
        %%%%%%%%%%%%%%%%
        
        for j = 1 : nConditions
            c = conditions(j);
            
            disp(['Condition ' num2str(c)])
            
            %set up IOWA vs IOWA-C contrast
            if (c < 6) 
                tFixEnd = 300;
            else 
                tFixEnd = 2600; %300;
            end
            
            % all stimuli have a size of 3*p.spatialFieldSize (to account for gaze shifts)
            spatialOrigin = pa.spatialFieldSize + spatialHalfSize + 1;
            
            input_x = zeros(tMax, 1);
            input_g = zeros(tMax, 1);
            
            input_g(:) = pa.stimStrength_g_base;
            input_x(1:tFixEnd) = pa.stimStrength_x;
            switch c
                case VALID_CUE
                    stimPos_v = [0, targetEccentricity, targetEccentricity];
                    stimSize_v = [pa.fixationSize, pa.cueSize, targetSize]; % should be odd
                    stimTimes_v = [1, tFixEnd; tCueStart, tCueEnd; tTargetStart, tTargetEnd];
                    input_g(tToneBoostStart:tToneBoostEnd) = pa.stimStrength_g_base + pa.stimStrength_g_tone;
                case INVALID_CUE
                    stimPos_v = [0, -targetEccentricity, targetEccentricity];
                    stimSize_v = [pa.fixationSize, pa.cueSize, targetSize]; % should be odd
                    stimTimes_v = [1, tFixEnd; tCueStart, tCueEnd; tTargetStart, tTargetEnd];
                    input_g(tToneBoostStart:tToneBoostEnd) = pa.stimStrength_g_base + pa.stimStrength_g_tone;
                case DOUBLE_CUE
                    stimPos_v = [0, -targetEccentricity, targetEccentricity, targetEccentricity];
                    stimSize_v = [pa.fixationSize, pa.cueSize, pa.cueSize, targetSize]; % should be odd
                    stimTimes_v = [1, tFixEnd; tCueStart, tCueEnd; tCueStart, tCueEnd; tTargetStart, tTargetEnd];
                    input_g(tToneBoostStart:tToneBoostEnd) = pa.stimStrength_g_base + pa.stimStrength_g_tone;
                case TONE_CUE
                    stimPos_v = [0, targetEccentricity];
                    stimSize_v = [pa.fixationSize, targetSize];
                    stimTimes_v = [1, tFixEnd; tTargetStart, tTargetEnd];
                    input_g(tToneBoostStart:tToneBoostEnd) = pa.stimStrength_g_base + pa.stimStrength_g_tone;
                case NO_CUE
                    stimPos_v = [0, targetEccentricity];
                    stimSize_v = [pa.fixationSize, targetSize];
                    stimTimes_v = [1, tFixEnd; tTargetStart, tTargetEnd];
                case CVALID_CUE
                    stimPos_v = [0, targetEccentricity, targetEccentricity];
                    stimSize_v = [pa.fixationSize, pa.cueSize, targetSize]; % should be odd
                    stimTimes_v = [1, tFixEnd; tCueStart, tCueEnd; tTargetStart, tTargetEnd];
                    input_g(tToneBoostStart:tToneBoostEnd) = pa.stimStrength_g_base + pa.stimStrength_g_tone;
                case CINVALID_CUE
                    stimPos_v = [0, -targetEccentricity, targetEccentricity];
                    stimSize_v = [pa.fixationSize, pa.cueSize, targetSize]; % should be odd
                    stimTimes_v = [1, tFixEnd; tCueStart, tCueEnd; tTargetStart, tTargetEnd];
                    input_g(tToneBoostStart:tToneBoostEnd) = pa.stimStrength_g_base + pa.stimStrength_g_tone;
                case CDOUBLE_CUE
                    stimPos_v = [0, -targetEccentricity, targetEccentricity, targetEccentricity];
                    stimSize_v = [pa.fixationSize, pa.cueSize, pa.cueSize, targetSize]; % should be odd
                    stimTimes_v = [1, tFixEnd; tCueStart, tCueEnd; tCueStart, tCueEnd; tTargetStart, tTargetEnd];
                    input_g(tToneBoostStart:tToneBoostEnd) = pa.stimStrength_g_base + pa.stimStrength_g_tone;
                case CTONE_CUE
                    stimPos_v = [0, targetEccentricity];
                    stimSize_v = [pa.fixationSize, targetSize];
                    stimTimes_v = [1, tFixEnd; tTargetStart, tTargetEnd];
                    input_g(tToneBoostStart:tToneBoostEnd) = pa.stimStrength_g_base + pa.stimStrength_g_tone;
                case CNO_CUE
                    stimPos_v = [0, targetEccentricity];
                    stimSize_v = [pa.fixationSize, targetSize];
                    stimTimes_v = [1, tFixEnd; tTargetStart, tTargetEnd];                    
            end
            
            % generate actual stimuli
            nStimuli_v = length(stimPos_v);
            
            stimuli_v = zeros(nStimuli_v, 3*pa.spatialFieldSize);
            stimWeights_v = zeros(tMax, nStimuli_v);
            
            for k = 1 : nStimuli_v
                if stimPos_v(k) >= 0
                    stimStart = stimPos_v(k) - ceil((stimSize_v(k)-1)/2) + spatialOrigin;
                    stimEnd = stimPos_v(k) + floor((stimSize_v(k)-1)/2) + spatialOrigin;
                else
                    stimStart = stimPos_v(k) - floor((stimSize_v(k)-1)/2) + spatialOrigin;
                    stimEnd = stimPos_v(k) + ceil((stimSize_v(k)-1)/2) + spatialOrigin;
                end
                spatialInput = zeros(1, 3*pa.spatialFieldSize);
                spatialInput(stimStart:stimEnd) = 1;
                stimuli_v(k, :) = spatialInput;
                stimWeights_v(stimTimes_v(k, 1):stimTimes_v(k, 2)-1, k) = 1;
            end
            
            saccStartTimesTemp = zeros(1,nRepeats);
            saccEndTimesTemp = zeros(1,nRepeats);
            saccMetricsTemp = zeros(1,nRepeats);


            %% loop over repetitions

            %need to track over ages = a and condition = c; repack after parfor
            for k = 1 : nRepeats
            %parfor k = 1 : nRepeats
                tic
                
                %% COSIVINA
                sim2 = 5;
                if (mode == 0) || (mode == 1)
                    sim2 = sim;
                elseif mode == 2
                    sim2 = sim.copy();
                end

                gazeDirection = spatialOrigin; % initial center of gaze within stimuli_v
                saccadeInProgress = false;
                saccMetricsIntegrator = 0;
                gazeChangeCounter = 0;
                
                shiftedStimuli_v = stimuli_v(:, gazeDirection-spatialHalfSize:gazeDirection+spatialHalfSize);
                convStimuli_v = conv2(1, kernel_v, shiftedStimuli_v, 'same');
                
                %% update COSIVINA noise parameters -- scale by tau
                sim2.setElementParameters({'noise kernel a'}, {'amplitude'}, [pa.q_a*pa.tau_a]); %multiplying by tau to move from in-line params to cosivina params
                sim2.setElementParameters({'noise kernel s'}, {'amplitude'}, [pa.q_s*pa.tau_s]); %multiplying by tau to move from in-line params to cosivina params
                sim2.setElementParameters('noise r', 'amplitude', pa.q_r*pa.tau_s); %multiplying by tau to move from in-line params to cosivina params
                sim2.setElementParameters('noise x', 'amplitude', pa.q_x*pa.tau_a); %multiplying by tau to move from in-line params to cosivina params
                sim2.setElementParameters('noise g', 'amplitude', pa.q_x*pa.tau_a); %multiplying by tau to move from in-line params to cosivina params               

                %% set up handles/structures to store history 
                if storeHistory
                    history_a = zeros(tMax, pa.spatialFieldSize);
                    history_s = zeros(tMax, pa.spatialFieldSize);
                    history_r = zeros(tMax, 1);
                    history_x = zeros(tMax, 1);
                    history_g = zeros(tMax, 1);                   
                end
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % actual trial starting here %
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                %% initialize sim and field values
                sim2.init();
                sim2.t=1;

                %% initialize fields
                handle_a = sim2.getElement('field a');
                handle_a.activation = pa.h_a;
                handle_s = sim2.getElement('field s');
                handle_s.activation = pa.h_s;
                handle_r = sim2.getElement('r');
                handle_r.activation = pa.h_r;
                handle_x = sim2.getElement('x');
                handle_x.activation = pa.h_x;
                handle_g = sim2.getElement('g');
                handle_g.activation = pa.h_g;
                handle_input_v = sim2.getElement('input_v');
                handle_input_x = sim2.getElement('input_x');
                handle_input_g = sim2.getElement('input_g');

                while sim2.t <= tMax

                    t = sim2.t;
                    
                    % determine start of saccade and suppress input
                    if (~saccadeInProgress && handle_r.output >= pa.theta_saccStart)
                        gazeChangeCounter = gazeChangeCounter + 1;
                        saccStartTimesTemp(k) = t;
                        saccadeInProgress = true;
                    end
                    
                    % determine end of saccade and update input
                    if (saccadeInProgress && handle_r.output < pa.theta_saccEnd)
                        saccEndTimesTemp(k) = t;
                        saccMetricsTemp(k) = saccMetricsIntegrator;
                        saccMetricsIntegrator = 0;
                        if breakOnFirstSaccade
                            break;
                        end
                        gazeDirection = gazeDirection + round(saccMetricsTemp(k));
                        shiftedStimuli_v = stimuli_v(:, gazeDirection-spatialHalfSize:gazeDirection+spatialHalfSize);
                        convStimuli_v = conv2(1, kernel_v, shiftedStimuli_v, 'same');
                        saccadeInProgress = false;
                    end
                    
                    % accumulate eye movement command
                    saccMetricsIntegrator = saccMetricsIntegrator + pa.c_sacc * sum(handle_s.output .* (-spatialHalfSize:spatialHalfSize));
                    
                    % update inputs
                    input_vTmp = ~saccadeInProgress * sum(convStimuli_v .* repmat(stimWeights_v(t, :)', [1, pa.spatialFieldSize]), 1);
                    handle_input_v.output = input_vTmp;
                    handle_input_x.output = input_x(t);
                    handle_input_g.output = input_g(t);
                    
                    %% store activation history
                    if storeHistory
                        history_a(t, :) = handle_a.activation;
                        history_s(t, :) = handle_s.activation;
                        history_r(t) = handle_r.activation;
                        history_x(t) = handle_x.activation;
                        history_g(t) = handle_g.activation;
                    end

                    %% COSIVINA step and manage GUI
                    if mode == 1 && (mod(t,gui_speed)==0)
                        gui.updateVisualizations();
                        gui.checkAndUpdateControls();
                    end
                    if ~gui.pauseSimulation
                        sim2.step();
                    end
                    if gui.loadParameters
                        if gui.loadParametersFromFile()
                            sim2.init();
                            sim2.t = 1; % done manually in this file, would be better to define tZero = 1 for simulator object
                        end
                        gui.loadParameters = false;
                    end
                    if gui.saveParameters
                        gui.saveParametersToFile();
                        gui.saveParameters = false;
                    end
                    if gui.quitSimulation
                        gui.close();
                        break;
                    end

                end %time
                
                if printTrials
                    disp(['Run ' num2str(k) ', Condition: ' num2str(c) ', Time: ' num2str(toc) ' sec, RT: ' ...
                        num2str(saccEndTimesTemp(k) - tTargetStart + pa.saccadeLatencyOffset) ' ms, Metrics: ' ...
                        num2str(saccMetricsTemp(k))])
                end

            end %repeats (parfor)
            
            saccStartTimes{a}(c, :) = saccStartTimesTemp(:);
            saccEndTimes{a}(c, :) = saccEndTimesTemp(:);
            saccMetrics{a}(c, :) = saccMetricsTemp(:);
            
        end %conditions
                
    end % end loop over ages
    
    
    %% analysis of results %%
    %%%%%%%%%%%%%%%%%%%%%%%%%
    
    for a = ages
        if ageRunPointers(r, a) == r
            for c = conditions
                ratesCorrect(c, a) = sum(saccMetrics{a}(c, :) > 0)/nRepeats;
                %         validTrials = saccMetrics{a}(c, :) > 0 & ...
                %           saccStartTimes{a}(c, :) >= minSaccLatency & saccStartTimes{a}(c, :) <= maxSaccLatency;
                RTs = saccEndTimes{a}(c, :) - tTargetStart + pr.saccadeLatencyOffset(a);
                validTrials = saccMetrics{a}(c, :) > 0 & RTs >= minSaccLatency & RTs <= maxSaccLatency;
                latenciesCorrect(c, a) = mean(RTs(validTrials));
            end
        else
            ratesCorrect(:, a) = result(ageRunPointers(r, a)).ratesCorrect(:, a);
            latenciesCorrect(:, a) = result(ageRunPointers(r, a)).latenciesCorrect(:, a);
        end
    end
    
    plotColors = cell(4, 1);
    plotColors{INFANT_5MO} = [0, 0, 0.8];
    plotColors{INFANT_7MO} = [0.8, 0, 0];
    plotColors{INFANT_10MO} = [0, 0.8, 0];
    
    conditionLabels = cell(N_CONDITIONS, 1);
    conditionLabels{DOUBLE_CUE} = 'double';
    conditionLabels{INVALID_CUE} = 'invalid';
    conditionLabels{NO_CUE} = 'none';
    conditionLabels{TONE_CUE} = 'tone';
    conditionLabels{VALID_CUE} = 'valid';
    conditionLabels{CDOUBLE_CUE} = 'doubleC';
    conditionLabels{CINVALID_CUE} = 'invalidC';
    conditionLabels{CNO_CUE} = 'noneC';
    conditionLabels{CTONE_CUE} = 'toneC';
    conditionLabels{CVALID_CUE} = 'validC';
    
    result(r) = struct('conditions', conditions, 'ages', ages, 'parameters', pr, ...
        'saccStartTimes', [], 'saccEndTimes', [], 'saccMetrics', [], 'latenciesCorrect', latenciesCorrect, ...
        'ratesCorrect', ratesCorrect, 'ageRunPointers', ageRunPointers(r, :), ...
        'noise', noise, 'nRepeats', nRepeats, 'errorRatesEstimated', false);
    if keepFullResults
        result(r).saccStartTimes = saccStartTimes;
        result(r).saccEndTimes = saccEndTimes;
        result(r).saccMetrics = saccMetrics;
    end
    if saveResults
        save(resultFilename, 'result');
    end
       
    if plotLatencies
        figure
        hLatencyAxis = axes('nextPlot', 'add', 'XLim', [1, nConditions], 'XTick', 1:nConditions, ...
            'XTickLabel', conditionLabels(conditions), 'YLim', [100, 550]);
        title(['Latencies - Run ' num2str(r)]);
        for a = ages
            plot(hLatencyAxis, latenciesCorrect(:, a), 'Color', plotColors{a}, 'Marker', 'square');
        end
    end
    
    if plotErrorRates
        figure
        hErrorAxis = axes('nextPlot', 'add', 'XLim', [1, nConditions], 'XTick', 1:nConditions, ...
            'XTickLabel', conditionLabels(conditions), 'YLim', [0.4, 1.0]);
        title(['Proportion of correct saccades - Run ' num2str(r)]);
        for a = ages
            plot(hErrorAxis, ratesCorrect(:, a), 'Color', plotColors{a}, 'Marker', 'square');
        end
    end
    
end % loop over modifications



