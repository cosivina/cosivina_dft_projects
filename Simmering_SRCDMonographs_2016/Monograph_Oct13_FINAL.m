%   Simulator from capacity paper (February 2012)
%   Entirely self-contained--no other model file (just need gauss and sigmoid)
%   June 2012
%   Starting point for new fits of thesis data for SRCD monograph
%   adding dynamic long-term memory
%   December 2012, adding fixation system from Dev Science paper
%   May 2013, adding a few new measures of looking dynamics

clear;
close all;

randn('state',sum(100*clock));
rand('state',sum(100*clock));

if 0  %repeated trials
    STORAGE=1;
    NOISE=1;        %use spatially correlated noise; if=0, use ODE45
    n_trials=6;     %12 for CD (trials per set size per "subject"), 6 for PL
    no_changeTR=n_trials/2;
    total_runs=40;  %number of "subjects" 20 for CD, 40 for PL
    STORE1=1;       %flag used to tag the first time an output file is created
    SimName='PL_retest'; %CD_AD_grant_SS7act'; %dev_FINAL';
    %for CD only
    startSS=1;
    stopSS=7;
    nSS=stopSS-startSS+1;
else
    STORAGE=0;
    NOISE=1; 		%use spatially correlated noise
    n_trials=1;
    total_runs=1;
    %for CD only
    startSS=3;
    stopSS=3;
    nSS=stopSS-startSS+1;
end

h_noise = 1;
track_act_U=0;
track_act_dec=0;
track_act_W=0;

nD=1;  %Specify the number of ages to be tested
lD=[1]; %List the age indices that will be tested; 1=3yo, 2=4yo, 3=5yo, 4=adult
dTrack=0; %Keep track of which index of the lD you are on, initialize at 0

%%%%% Parameter values

n_field=361;           %must be odd;
C=fix(n_field/2)+1;
D2U=1;
U2D=1;
MS2TS=.5;
TS2MS=2;

task=2;
%1=change detection (CD)
%2=preferential looking (PL)
%DON'T FORGET TO COMMENT OUT BLOCKS AND CHANGE TRIALS FOR PL

%if running full set
if task == 1
    total_iterations=total_runs*nSS*n_trials*nD;  %add dev manually
else
    total_iterations=total_runs*n_trials*nD;  %add dev manually
end
current_iteration=0;

disp('starting simulation');
tic_total=tic;
running_time=0;

% timing parameters: 1 time step = 2ms
%ENTER VALUES IN MILLISECONDS FROM STUDY, FORMLUAE TRANSLATE FOR YOU
init_time=floor(100*MS2TS);
update_time=floor(0*MS2TS);
if task == 1
    sample_time=floor(2000*MS2TS);
    delay=floor(900*MS2TS);
    test_time=floor(2000*MS2TS);    %maximum (i.e., without response - code added below to stop trial after response is given
    n_time =init_time+sample_time+delay+test_time+update_time;
else
    sample_time=floor(500*MS2TS);   %all "on" periods
    delay=floor(250*MS2TS);         %all "off" periods
    test_time=floor(0*MS2TS);
    n_stream_events=27;             %27 for 20s trial; 13 for 10s trial
    att_getter=(400*MS2TS);         
    pre_dur=(150*MS2TS);            %before stimuli come on, after attention getter
    stream_on=n_stream_events*(sample_time+delay)-delay;  %time from first to last presentation of stimuli (20s); have to subtract final delay
    n_time=init_time+att_getter+pre_dur+stream_on;  %this was wrong before, still not sure
end
t_init=0.; t_final=n_time-1;


%SimName='Testing';
DataName='DATA_testing';

if task == 1
    gateTrack=zeros(nSS,n_trials,n_time);
    idTrack=zeros(nSS,n_trials,n_time);
    isTrack=zeros(nSS,n_trials,n_time);
end

if track_act_U == 1
    SS1_U_time=zeros(n_time,n_trials);
    SS2_U_time=zeros(n_time,n_trials);
    SS3_U_time=zeros(n_time,n_trials);
    SS4_U_time=zeros(n_time,n_trials);
    SS5_U_time=zeros(n_time,n_trials);
    SS6_U_time=zeros(n_time,n_trials);
    
    SS1_peakU_time=zeros(n_time,n_trials);
    SS2_peakU_time=zeros(n_time,n_trials);
    SS3_peakU_time=zeros(n_time,n_trials);
    SS4_peakU_time=zeros(n_time,n_trials);
    SS5_peakU_time=zeros(n_time,n_trials);
    SS6_peakU_time=zeros(n_time,n_trials);
end

if track_act_W == 1
    SS1_W_time=zeros(n_time,n_trials);
    SS2_W_time=zeros(n_time,n_trials);
    SS3_W_time=zeros(n_time,n_trials);
    SS4_W_time=zeros(n_time,n_trials);
    SS5_W_time=zeros(n_time,n_trials);
    SS6_W_time=zeros(n_time,n_trials);
    
    SS1_peakW_time=zeros(n_time,n_trials);
    SS2_peakW_time=zeros(n_time,n_trials);
    SS3_peakW_time=zeros(n_time,n_trials);
    SS4_peakW_time=zeros(n_time,n_trials);
    SS5_peakW_time=zeros(n_time,n_trials);
    SS6_peakW_time=zeros(n_time,n_trials);
end

if track_act_dec == 1 && task == 1
    SS1_S_time=zeros(n_time,n_trials);
    SS2_S_time=zeros(n_time,n_trials);
    SS3_S_time=zeros(n_time,n_trials);
    SS4_S_time=zeros(n_time,n_trials);
    SS5_S_time=zeros(n_time,n_trials);
    SS6_S_time=zeros(n_time,n_trials);

    SS1_D_time=zeros(n_time,n_trials);
    SS2_D_time=zeros(n_time,n_trials);
    SS3_D_time=zeros(n_time,n_trials);
    SS4_D_time=zeros(n_time,n_trials);
    SS5_D_time=zeros(n_time,n_trials);
    SS6_D_time=zeros(n_time,n_trials);
end

peakStore=zeros(nSS,n_trials);

for tdct = 1:nD
    dct = lD(tdct);
    dev_sigma_input=[ 1.5,  1.5,  1.5,  1, 1];
    dev_winput=     [ 0.3,0.325,0.355,  1, 1]; 
    dev_uv=         [0.35, 0.35, 0.35,  1, 1];
    dev_wv=         [0.55, 0.58, 0.62,  1, 1]; %CHANGED FROM 3y TO 5y IN THESIS  %same as dev_uv in thesis
    dev_global=     [0.75, 0.75, 0.75,  1, 1];
    dev_uu=         [0.85, 0.85, 0.85,  1, 1];
    dev_ww=         [ 0.4, 0.45, 0.53,  1, 1]; %CHANGED FROM 3y TO 5y IN THESIS
    %dev_winput_w=   [ 1.0,  1.0,  1.0,  1, 1]; %CHANGED FROM 3y TO 5y IN THESIS
    dev_noise=      [1.25, 1.25, 1.25,  1, 1]; %for fields %%NOISE WAS DIFFERENT IN THESIS
    dev_pre_build=  [0.25,  0.5, 0.75,  1, 1]; %not in capacity params; ON TAU, NOT STRENGTH
    dev_pre_decay=  [ 1.5,  1.5,  1.5,  1, 1];
    if task == 1
        %decision nodes
        dev_h_is=       [0.25, 0.25, 0.25,  0, .25]; %ADDED, not multiplied
        dev_noise_node= [1.25, 1.25, 1.25,  1, 1];
        dev_wi_node=    [0.75, 0.75, 0.75,  1, 1];
        dev_we_node=    [0.75, 0.75, 0.75,  1, 1];
        dev_diff=       [0.5357,0.75, .89,  1, 1.25];
        dev_same=       [0.66, 0.66, 0.66,  1, 1]; %CHANGED FROM 3y TO 5y IN THESIS  %0.3 in thesis
        dev_w_gate=     [ 5/3,  5/3,  5/3,  1, 1]; %not in thesis params
    else
        %dev_noise_fix=  [   1,   1,   1,   1,    1];  %6.67 for thesis value
        %dev_we_fix=     [   1,   1,   1,   1,1.425];  %not in thesis
        dev_wi_fix=     [  .9,  .9,  .9,   1,    2];   %8 for thesis value
        dev_const=      [ .93, .93, .93,   1,  1.1];   %2.89 for thesis value   %to L/R nodes (when stimulus is on)
        dev_w_in=       [1.14,1.14,1.14,   1,0.475];   %3.4 for thesis value    %from U to L/R nodes
        %dev_slope=      [   1,   1,   1,   1,    1];  %1.2 for thesis value    %multiplied by h_slope
        dev_h=          [ .82,.825, .84,   1,1.325];   %4.85 for thesis value   %multiplied by h_down
        %making this stronger decreased total looking and increased switching
    end
    
    for rct=1:total_runs
        if task == 1
            changeTrack=zeros(nSS,n_trials);
        end
        %COMMENT OUT FOR PL
%         for SS=startSS:stopSS
%             maxTrialCt=[n_trials/2 n_trials/2];
%             trialTypeCt=[0 0];
            
            %COMMENT OUT FOR CD
                    trials=[2 2 4 4 6 6];
                    order=perms(trials);
                    SS_order=order(round(rand*720),:);  %1+rand*720
                    end_time = n_time;
            
            for trial_ct=1:n_trials
                tic_iteration=tic;
                
                if trial_ct == 1
                    upre_field=zeros(1,n_field);
                    wpre_field=zeros(1,n_field);
                end
                if task == 1
                    TrialOver=0;
                    RT=0;
                    end_time = n_time; %setting this to be reset from RT of prev trial
                    if trialTypeCt(1)==maxTrialCt(1)
                        change = 1;
                        changeTrack(SS,trial_ct)=1;
                        trialTypeCt(2)=trialTypeCt(2)+1;
                    elseif trialTypeCt(2)==maxTrialCt(2)
                        change=0;
                        changeTrack(SS,trial_ct)=0;
                        trialTypeCt(1)=trialTypeCt(1)+1;
                    else
                        change=round(rand);
                        changeTrack(SS,trial_ct)=change;
                        trialTypeCt(change+1)=trialTypeCt(change+1)+1;
                    end
                else
                    switches=0;  %reset from switches of prev trial
                    n_looks=0;
                    last_L=0;
                    last_R=0;
                    blank_ends=zeros(1,27);
                    blank_ends(1,1)=700;
                    ct_PFact_NC=0;
                    ct_PFact_CH=0;
                    peak_ct_PL=zeros(1,27);
                end
                
                for SS = SS_order(trial_ct)            %COMMENT OUT FOR CD
                
                stim_str=(30*dev_winput(dct)); 
                stim_width=3*dev_sigma_input(dct);
                
                stim_values=[-160,-120,-80,-40,0,40,80,120,160];
                stim_order=perms(stim_values);
                stim1_pos=stim_order(round(1+rand*362879),:);
                
                stim2_pos =stim1_pos;
                if task == 1
                    Gstim_str=30; 
                    Gstim_width=3;
                    Gstim1_strength=[Gstim_str,Gstim_str,Gstim_str,Gstim_str,Gstim_str,Gstim_str,Gstim_str,Gstim_str];
                    Gstim1_width  =[Gstim_width,Gstim_width,Gstim_width,Gstim_width,Gstim_width,Gstim_width,Gstim_width,Gstim_width];
                    Gstim2_strength=Gstim1_strength;
                    Gstim2_width  =Gstim1_width;
                    
                    if change == 1
                        stim2_pos(1,1) =stim1_pos(1,9);
                    end
                end
                
                % stim1
                n_stim1_on  =init_time;
                n_stim1_off  =n_stim1_on+sample_time;
                ct_stim1=SS;
                
                stim1_strength=[stim_str,stim_str,stim_str,stim_str,stim_str,stim_str,stim_str,stim_str];
                stim1_width  =[stim_width,stim_width,stim_width,stim_width,stim_width,stim_width,stim_width,stim_width];
                
                % stim2
                n_stim2_on =n_stim1_off+delay;
                n_stim2_off =n_stim2_on+test_time;
                ct_stim2=ct_stim1;
                
                stim2_strength=stim1_strength;
                stim2_width  =stim1_width;
                
                % h  parameters
                h_u=-6.75;  %thesis = -7 for PL
                h_v=-12;
                h_w=-4.5;
                
                if task == 1
                    %Decision nodes
                    %beta_g = 0.5;
                    tau_g = 80;
                    h_g = -4.8; 
                    w_igw = .03*dev_w_gate(dct);           %from w to gate node
                    w_gS = .01;                            %from stim input to gate node
                    transient_time = floor(30*MS2TS);      %duration of transient associated with each stim
                    w_trans = 18;                          %strength of transient input
                    w_gg = 4.;                             %self-excitation of gate node
                    w_ig = 4.25;                           %boost from gate node to decision nodes when gate node "on"
                    
                    beta_i = 5;
                    tau_i = 80;
                    h_is = -4.75 + (dev_h_is(dct)); %-4.35 for capacity paper
                    h_id = -5;
                    w_ii=1.85*dev_we_node(dct);
                    ci=14*dev_wi_node(dct);
                    
                    w_idu = 1.4*dev_diff(dct); 
                    w_isw = .0275*dev_same(dct);
                    w_uid = 1;
                    w_wis = 1;
                else
                    %Fixation nodes - UPDATED WITH DEV SCIENCE PARAMETERS 1/3/2013
                    beta_i=5;
                    tau_i=80;
                    h_i=-5; %-5.5 in thesis;
                    w_ii=2.85;                                             %self-excitation
                    ci=5*dev_wi_fix(dct); %20 in thesis?                   %competitive inhibition
                    
                    c_s=6.16*dev_const(dct); %14.5 in thesis               %strength of constant input to L/R nodes (s_input in thesis)
                    %this was multiplied by a dev param in thesis, but it worked differently there than described in DevSci (was on when each stimulus appeared, which is trans now)
                    c_sa=5.6; %2.5 in thesis;                               %strength of constant input for A node  (a_input in thesis)
                    %ADD DEV CHANGE?
                    
                    c_ag=5.5; %10 in thesis;                               %strength of att-getter To center node (center_boost in thesis)
                    c_b=1; %2 in thesis;                                   %strength of boost to L/R nodes to begin each trial
                    %descrition in thesis, didn't match DevSci: Strength of Transient Attention Boost To Left and Right Nodes
                    %(boost in thesis) - in thesis it was on while stimuli were on
                    n_steps_boost=125; %75 ts in thesis                    %duration of attention boost to L/R nodes at trial beginning
                    
                    c_t=0.4;                                               %strength of transient boost to L/R nodes for first half of each stimulus
                    %wasn't in thesis, boost/constant(spatial)/trans didn't match with description in DevSci
                    
                    noisyinput=3; 
                    %*dev_noise_fix(dct); %20 in thesis         %Strength of Noisy Input
                    %uses noisy inputs to nodes, rather than noise on nodes
                    %NOT REPORTED IN DEV SCI, FOUND IN SIMULATOR
                    w_iu=1.4*dev_w_in(dct); %1 in thesis                   %Strength of U to L/R Nodes (c_fu in DevSci)
                    w_ui=1;                                                %Strength of L/R Nodes to U (c_uf in DevSci)
                    
                    % h dynamics for nodes
                    tau_h=80;
                    h_slope=.75; %.9*dev_slope(dct) in thesis              %(not reported in DevSci, found in simulator file from Sammy)
                    h_rest=5; %5.5 in thesis;                              %(reported as negative in DevSci)
                    h_down=4.1075*dev_h(dct); %15 in thesis;               %(reported as negative in DevSci)
                end
                
                noise_sigma=1;
                noise_strength=0.04*dev_noise(dct);
                if task == 1
                    noise_strength_nodes=0.065*dev_noise_node(dct);
                    noise_strength_gate=.025;
                else
                    noise_strength_fix=0.065; %*dev_noise_fix(dct) in thesis
                end
                if h_noise == 1
                    h_noisestr_u=6;
                    h_noisestr_v=6;
                    h_noisestr_w=6;
                else
                    h_noisestr_u=0;
                    h_noisestr_v=0;
                    h_noisestr_w=0;
                end
                tau_noise=50;
                
                % preshape dynamics
                upre_sigma=6.4;
                w_upre=0.2;     %into Upre
                wpre_sigma=6.4;
                w_wpre=0.2;     %into Wpre
                % Lipinski PBR paper: stregths = 0.2, sigmas = 6.4 (AD)
                % Samuelson PLoS One: stregths = 0.8, sigmas = 1 (1.5yo)
                % Perone Dev Sciece: strenghts = 0.25, sigmas = 5(u), 15(w)
                tau_pre_u_build=3000*dev_pre_build(dct);
                tau_pre_w_build=3000*dev_pre_build(dct);
                tau_pre_u_decay=100000*dev_pre_decay(dct);
                tau_pre_w_decay=100000*dev_pre_decay(dct);
                % Lipinski PBR paper: build = 3000, decay = 100000 (AD)
                % Samuelson PLoS One: build = 4500, decay = 50000 (1.5yo)
                % Perone Dev Sciece: build = 5000, decay = 25000 (6-10m)

                % Amari dynamics parameters
                tau_u=80.;
                tau_v=10.;
                tau_w=80.;
                beta=5;
                
                % kernels: width and strength parameters
                uu_sigma=3;
                w_uu=2.0*dev_uu(dct);
                
                uv_sigma=26;
                w_uv=1.85*dev_uv(dct);
                w_uv_const=0.05*dev_global(dct);
                
                vu_sigma=10;
                w_vu=2.0;
                
                vw_sigma=5;
                w_vw=1.95;
                
                wv_sigma=42;
                w_wv=(0.325*dev_wv(dct));
                w_wv_const=0.08*dev_global(dct);
                
                wu_sigma=5;
                w_wu=1.6;
                
                ww_sigma=3;
                w_ww=(3.15*dev_ww(dct));
                w_wS=.2; %*dev_winput_w(dct) in thesis 
                
                % input kernels
                stim1_input=zeros(1,n_field);
                for ct=1:ct_stim1
                    stim1_input=stim1_input+gauss(n_field-1,(fix(stim1_pos(ct)*D2U)+C)-1,stim1_strength(ct),stim1_width(ct),0);
                end
                
                stim2_input=zeros(1,n_field);
                for ct=1:ct_stim2
                    stim2_input=stim2_input+gauss(n_field-1,(fix(stim2_pos(ct)*D2U)+C)-1,stim2_strength(ct),stim2_width(ct),0);
                end
                
                if task == 1
                    Gstim1_input=zeros(1,n_field);
                    for ct=1:ct_stim1
                        Gstim1_input=Gstim1_input+gauss(n_field-1,(fix(stim1_pos(ct)*D2U)+C)-1,Gstim1_strength(ct),Gstim1_width(ct),0);
                    end
                    
                    Gstim2_input=zeros(1,n_field);
                    for ct=1:ct_stim2
                        Gstim2_input=Gstim2_input+gauss(n_field-1,(fix(stim2_pos(ct)*D2U)+C)-1,Gstim2_strength(ct),Gstim2_width(ct),0);
                    end
                end
                
                %%%% Constructing the stimulus
                if task == 1
                    stimulus=zeros(n_time,n_field);
                    stimulus(n_stim1_on+1:n_stim1_off,:)=ones(n_stim1_off-n_stim1_on,1)*(stim1_input);
                    stimulus(n_stim2_on+1:n_stim2_off,:)=ones(n_stim2_off-n_stim2_on,1)*(stim2_input);
                    
                    Gstimulus=zeros(n_time,n_field);
                    Gstimulus(n_stim1_on+1:n_stim1_off,:)=ones(n_stim1_off-n_stim1_on,1)*(Gstim1_input);
                    Gstimulus(n_stim2_on+1:n_stim2_off,:)=ones(n_stim2_off-n_stim2_on,1)*(Gstim2_input);
                    
                    transient=zeros(n_time,1);
                    transient(n_stim1_on+1:n_stim1_on+transient_time,1)=ones(transient_time,1)*w_trans;
                    transient(n_stim2_on+1:n_stim2_on+transient_time,1)=ones(transient_time,1)*w_trans;
                else
                    PL_stimulus     %separate file that constructs the streams
                    trial_begin = init_time + att_getter + pre_dur;
                    
                    stimulus_pre_att=zeros(1,init_time);                   %before attention getter (relaxation time)
                    stimulus_pre=zeros(n_field,trial_begin);               %"blank" time before stimuli are presented (includes att getter and relaxation time)
                    
                    %STREAM STRUCTURE
                    cycle_on=ones(1,sample_time);           %500 ms
                    cycle_off=zeros(1,delay);               %250 ms
                    cycle_duration=[cycle_on cycle_off];    %750 ms per cycle
                    
                    %CENTER NODE ATTENTION GETTER
                    stimulus_center=zeros(1,n_time);
                    center_on=ones(1,att_getter);
                    stimulus_center(1,(1+init_time):(init_time+att_getter))=center_on;
                    %sets input to 1 during attention getter, 0 before and after
                    
                    %BOOST TO LEFT AND RIGHT NODES
                    stimulus_boost=zeros(1,n_time);
                    boost_on=ones(1,n_steps_boost);
                    stimulus_boost(1,(1+trial_begin:(trial_begin+n_steps_boost)))=boost_on;
                    %first 125ms of each trial, follows relaxation time, attention getter, and "pre_dur"
                    
                    %CONSTANT INPUT FOR LEFT/RIGHT/AWAY NODES
                    stimulus_const=ones(1,n_time);
                    const_off=zeros(1,(init_time + att_getter + pre_dur));
                    stimulus_const(1,1:(init_time + att_getter + pre_dur))=const_off;
                    %set to 1 for constant input once stimuli start; does not turn off during delays (multipliers for L/R vs A set above)
                    
                    %TRANSIENT FOR LEFT AND RIGHT NODES
                    trans_on=ones(1,(sample_time/2));
                    trans_off=zeros(1,(sample_time/2)+delay);
                    stimulus_trans_single=[trans_on trans_off];
                    %on for the first half (250ms) of each "on" cycle, constructed below
                    
                    %STIMULUS STREAM (on/off) FOR PREFERENTIAL LOOKING
                    for j=1:n_stream_events
                        if j==1
                            stimulus_NC=[stimulus_pre stim_NC_input(j,:)'*cycle_duration];
                            stimulus_CH=[stimulus_pre stim_CH_input(j,:)'*cycle_duration];
                            stimulus_trans=[const_off c_t.*stimulus_trans_single];
                        else
                            stimulus_NC=[stimulus_NC stim_NC_input(j,:)'*cycle_duration];
                            stimulus_CH=[stimulus_CH stim_CH_input(j,:)'*cycle_duration];
                            stimulus_trans=[stimulus_trans c_t.*stimulus_trans_single];
                        end  %j
                    end  %j stream events
                    stimulus_trans=stimulus_trans(1,1:n_time);
                    
                    %complete inputs for looking nodes
                    stimulus_C     =(c_ag.*stimulus_center);
                    stimulus_A     =(c_sa.*stimulus_const);
                    stimulus_LR    =(c_b.*stimulus_boost) + (c_s.*stimulus_const) + (c_t.*stimulus_trans);
                    %needed to make noisy inputs indepenent
                    stimulus_L=stimulus_LR;
                    stimulus_R=stimulus_LR;
                    
                    for bt=2:27
                        blank_ends(1,bt)=blank_ends(1,bt-1)+375;
                    end
                    blank=1;
                    
                    
                end
                
                %%%% preparing the kernels
                kernel_width_multiplier=3; %DO NOT SET TO SMALLER THAN 3
                
                uu_kernel_width=floor(min(kernel_width_multiplier*uu_sigma,(n_field-1)/2));
                uu_int_kernel=(w_uu)*(exp(-.5*(-uu_kernel_width:uu_kernel_width).^2/uu_sigma^2));
                uu_ext_index=[n_field-uu_kernel_width+1:n_field, 1:n_field, 1: uu_kernel_width];
                
                uv_kernel_width=floor(min(kernel_width_multiplier*uv_sigma,(n_field-1)/2));
                uv_int_kernel=(w_uv)*(exp(-.5*(-uv_kernel_width:uv_kernel_width).^2/uv_sigma^2));
                uv_ext_index=[n_field-uv_kernel_width+1:n_field, 1:n_field, 1: uv_kernel_width];
                
                vu_kernel_width=floor(min(kernel_width_multiplier*vu_sigma,(n_field-1)/2));
                vu_int_kernel=(w_vu)*(exp(-.5*(-vu_kernel_width:vu_kernel_width).^2/vu_sigma^2));
                vu_ext_index=[n_field-vu_kernel_width+1:n_field, 1:n_field, 1: vu_kernel_width];
                
                vw_kernel_width=floor(min(kernel_width_multiplier*vw_sigma,(n_field-1)/2));
                vw_int_kernel=(w_vw)*(exp(-.5*(-vw_kernel_width:vw_kernel_width).^2/vw_sigma^2));
                vw_ext_index=[n_field-vw_kernel_width+1:n_field, 1:n_field, 1: vw_kernel_width];
                
                wu_kernel_width=floor(min(kernel_width_multiplier*wu_sigma,(n_field-1)/2));
                wu_int_kernel=(w_wu)*(exp(-.5*(-wu_kernel_width:wu_kernel_width).^2/wu_sigma^2));
                wu_ext_index=[n_field-wu_kernel_width+1:n_field, 1:n_field, 1: wu_kernel_width];
                
                wv_kernel_width=floor(min(kernel_width_multiplier*wv_sigma,(n_field-1)/2));
                wv_int_kernel=(w_wv)*(exp(-.5*(-wv_kernel_width:wv_kernel_width).^2/wv_sigma^2));
                wv_ext_index=[n_field-wv_kernel_width+1:n_field, 1:n_field, 1: wv_kernel_width];
                
                ww_kernel_width=floor(min(kernel_width_multiplier*ww_sigma,(n_field-1)/2));
                ww_int_kernel=(w_ww)*(exp(-.5*(-ww_kernel_width:ww_kernel_width).^2/ww_sigma^2));
                ww_ext_index=[n_field-ww_kernel_width+1:n_field, 1:n_field, 1: ww_kernel_width];
                
                upre_kernel_width=floor(min(kernel_width_multiplier*upre_sigma,(n_field-1)/2));
                upre_int_kernel=(w_upre)*(exp(-.5*(-upre_kernel_width:upre_kernel_width).^2/upre_sigma^2));
                upre_ext_index=[n_field-upre_kernel_width+1:n_field, 1:n_field, 1: upre_kernel_width];
                
                wpre_kernel_width=floor(min(kernel_width_multiplier*wpre_sigma,(n_field-1)/2));
                wpre_int_kernel=(w_wpre)*(exp(-.5*(-wpre_kernel_width:wpre_kernel_width).^2/wpre_sigma^2));
                wpre_ext_index=[n_field-wpre_kernel_width+1:n_field, 1:n_field, 1: wpre_kernel_width];
                
                % stochastics
                noise_kernel_width=floor(min(kernel_width_multiplier*noise_sigma,(n_field-1)/2));
                noise_kernel=(noise_strength)*(exp(-.5*(-noise_kernel_width:noise_kernel_width).^2/noise_sigma^2));
                noise_size=n_field+2*noise_kernel_width;
                
                %%%%% initiatlization, initial conditions
                u_field =zeros(n_time,n_field) + h_u;
                v_field =zeros(n_time,n_field) + h_v;
                w_field =zeros(n_time,n_field) + h_w;
                upre_field =zeros(n_time,n_field);
                wpre_field =zeros(n_time,n_field);
                if h_noise == 1
                    h_noise_u=zeros(n_time,1);
                    h_noise_v=zeros(n_time,1);
                    h_noise_w=zeros(n_time,1);
                end
                if task == 1
                    gate_node=zeros(n_time,1) + h_g;
                    is_node =zeros(n_time,1) + h_is;
                    id_node =zeros(n_time,1) + h_id;
                else
                    il_node =zeros(n_time,1) + h_i;
                    ir_node =zeros(n_time,1) + h_i;
                    ic_node =zeros(n_time,1) + h_i;
                    ia_node =zeros(n_time,1) + h_i;
                    h_il =zeros(n_time,1) + h_i;
                    h_ir =zeros(n_time,1) + h_i;
                    h_ic =zeros(n_time,1) + h_i;
                    h_ia =zeros(n_time,1) + h_i;
                    input=zeros(n_time,n_field);  %to track different stimulus inputs over looking
                    look_L =zeros(n_time,1);
                    look_R =zeros(n_time,1);
                    switchLtoR =zeros(n_time,1);
                    switchRtoL =zeros(n_time,1);
                    stable_L =zeros(n_time,1);
                    stable_R =zeros(n_time,1);
                end
                
                for time=2:end_time
                    
                    f_u=sigmoid(u_field(time-1,:),beta,0);
                    f_v=sigmoid(v_field(time-1,:),beta,0);
                    f_w=sigmoid(w_field(time-1,:),beta,0);
                    uu_conv=conv2(1,uu_int_kernel,f_u(uu_ext_index),'valid');
                    uv_conv=conv2(1,uv_int_kernel,f_v(uv_ext_index),'valid');
                    vu_conv=conv2(1,vu_int_kernel,f_u(vu_ext_index),'valid');
                    vw_conv=conv2(1,vw_int_kernel,f_w(vw_ext_index),'valid');
                    wu_conv=conv2(1,wu_int_kernel,f_u(wu_ext_index),'valid');
                    wv_conv=conv2(1,wv_int_kernel,f_v(wv_ext_index),'valid');
                    ww_conv=conv2(1,ww_int_kernel,f_w(ww_ext_index),'valid');
                    upre_conv=conv2(1,upre_int_kernel,upre_field(time-1,upre_ext_index),'valid');
                    wpre_conv=conv2(1,wpre_int_kernel,wpre_field(time-1,wpre_ext_index),'valid');
                    if(NOISE)
                        u_noise=conv2(1,noise_kernel,randn(1,noise_size),'valid');
                        v_noise=conv2(1,noise_kernel,randn(1,noise_size),'valid');
                        w_noise=conv2(1,noise_kernel,randn(1,noise_size),'valid');
                    else
                        u_noise=0;
                        v_noise=0;
                        w_noise=0;
                    end %noise
                    if h_noise == 1
                        h_noise_u(time,:) = h_noise_u(time-1,:) + (-h_noise_u(time-1,:) + h_noisestr_u*randn)/tau_noise;
                        h_noise_v(time,:) = h_noise_v(time-1,:) + (-h_noise_v(time-1,:) + h_noisestr_v*randn)/tau_noise;
                        h_noise_w(time,:) = h_noise_w(time-1,:) + (-h_noise_w(time-1,:) + h_noisestr_w*randn)/tau_noise;
                    else
                        h_noise_u=0;
                        h_noise_v=0;
                        h_noise_w=0;
                    end
                    
                    upre_noise=0;
                    wpre_noise=0;
                    
                    if task == 1
                        %COUPLING DECISION NODES TO FIELDS
                        s_ig=sigmoid(gate_node(time-1,:),beta_i,0);
                        s_id=sigmoid(id_node(time-1,:),beta_i,0);
                        s_is=sigmoid(is_node(time-1,:),beta_i,0);
                        input_from_id=w_uid*s_id;
                        input_from_is=w_wis*s_is;
                        input_u_to_id=w_idu*sum(sigmoid(u_field(time-1,:),beta,0));
                        input_w_to_is=w_isw*sum(sigmoid(w_field(time-1,:),beta,0));
                        input_w_to_ig=w_igw*sum(sigmoid(w_field(time-1,:),beta,0));
                        if(NOISE)  %adding extra noise to nodes to try to reduce perfect performance
                            id_noise=noise_strength_nodes*randn(1,1);
                            is_noise=noise_strength_nodes*randn(1,1);
                            ig_noise=noise_strength_gate*randn(1,1);
                        else
                            id_noise=0;
                            is_noise=0;
                            ig_noise=0;
                        end %noise
                        
                        S = stimulus(time,:);
                        GS = Gstimulus(time,:);
                        stim_to_ig = w_gS*sum(sigmoid(GS,beta,0));
                        
                        % equations
                        gate_node(time,:) = gate_node(time-1,:) + (- gate_node(time-1,:) + h_g + w_gg*s_ig + input_w_to_ig + stim_to_ig + transient(time,1))/tau_g + ig_noise;
                        id_node(time,:) = id_node(time-1,:) + (- id_node(time-1,:) + h_id + w_ii*s_id + (input_u_to_id + w_ig)*s_ig - ci*s_is)/tau_i + id_noise;
                        is_node(time,:) = is_node(time-1,:) + (- is_node(time-1,:) + h_is + w_ii*s_is + (input_w_to_is + w_ig)*s_ig - ci*s_id)/tau_i + is_noise;
                        
                        gateTrack(SS,trial_ct,time)=gate_node(time);
                        idTrack(SS,trial_ct,time)=id_node(time);
                        isTrack(SS,trial_ct,time)=is_node(time);
                        
                        if h_noise == 1
                            u_field(time,:) = u_field(time-1,:) + (- u_field(time-1,:) + h_u + h_noise_u(time-1,:) + input_from_id + S      + uu_conv           - uv_conv - w_uv_const*sum(f_v)+ upre_conv)/tau_u + u_noise;
                            v_field(time,:) = v_field(time-1,:) + (- v_field(time-1,:) + h_v + h_noise_v(time-1,:)                          + vu_conv + vw_conv                                           )/tau_v + v_noise;
                            w_field(time,:) = w_field(time-1,:) + (- w_field(time-1,:) + h_w + h_noise_w(time-1,:) + input_from_is + S*w_wS + wu_conv + ww_conv - wv_conv - w_wv_const*sum(f_v)+ wpre_conv)/tau_w + w_noise;
                        else
                            u_field(time,:) = u_field(time-1,:) + (- u_field(time-1,:) + h_u + input_from_id + S      + uu_conv           - uv_conv - w_uv_const*sum(f_v)+ upre_conv)/tau_u + u_noise;
                            v_field(time,:) = v_field(time-1,:) + (- v_field(time-1,:) + h_v +                        + vu_conv + vw_conv                                           )/tau_v + v_noise;
                            w_field(time,:) = w_field(time-1,:) + (- w_field(time-1,:) + h_w + input_from_is + S*w_wS + wu_conv + ww_conv - wv_conv - w_wv_const*sum(f_v)+ wpre_conv)/tau_w + w_noise;
                        end
                        
                        theta_upre=u_field(time,:)>0;
                        theta_wpre=w_field(time,:)>0;
                        upre_field(time,:) = upre_field(time-1,:) + (theta_upre).*(- upre_field(time-1,:) + f_u)/tau_pre_u_build + (1-theta_upre).*(-upre_field(time-1,:))/tau_pre_u_decay + upre_noise;
                        wpre_field(time,:) = wpre_field(time-1,:) + (theta_wpre).*(- wpre_field(time-1,:) + f_w)/tau_pre_w_build + (1-theta_wpre).*(-wpre_field(time-1,:))/tau_pre_w_decay + wpre_noise;  
                        
                        %counting peaks in U at each time step
                        peakDistance = 1;
                        supraThreshold = find(f_u >= 0.5)';
                        if size(supraThreshold, 1) < 2  %updated 8/16/13
                            nPeaksU = 0;
                        else
                            % the actual clustering
                            clusterIndices = clusterdata(supraThreshold, 'distance', 'euclidean', 'linkage',...
                                'single', 'cutoff', peakDistance, 'criterion', 'inconsistent');
                            % determining the centers of the clusters
                            nPeaksU = max(clusterIndices);
                        end
                        
                        %counting peaks in W at each time step
                        peakDistance = 1;
                        supraThreshold = find(f_w >= 0.5)';
                        if size(supraThreshold, 1) < 2  %updated 8/16/13
                            nPeaksW = 0;
                        else
                            % the actual clustering
                            clusterIndices = clusterdata(supraThreshold, 'distance', 'euclidean', 'linkage',...
                                'single', 'cutoff', peakDistance, 'criterion', 'inconsistent');
                            % determining the centers of the clusters
                            nPeaksW = max(clusterIndices);
                        end
                        
                        %%% for reaction times
                        if TrialOver==0
                            if (id_node(time,:)>0 && is_node(time,:)<0) || (id_node(time,:)<0 && is_node(time,:)>0)
                                RT=time;
                                TrialOver=1;
                                end_time = RT + 250;
                                if end_time > n_time
                                    end_time = n_time;
                                end
                            end
                        end
                        
                        %%% counting the number of peaks at the end of delay in CD
                        if time == n_stim2_on-1
                            peak_ct_CD = nPeaksW;
                        end
                        
                        if track_act_dec == 1 && rct == total_runs
                            if SS == 1
                                SS1_S_time(time,trial_ct)=is_node(time,:);
                                SS1_D_time(time,trial_ct)=id_node(time,:);
                            elseif SS == 2
                                SS2_S_time(time,trial_ct)=is_node(time,:);
                                SS2_D_time(time,trial_ct)=id_node(time,:);
                            elseif SS == 3
                                SS3_S_time(time,trial_ct)=is_node(time,:);
                                SS3_D_time(time,trial_ct)=id_node(time,:);
                            elseif SS == 4
                                SS4_S_time(time,trial_ct)=is_node(time,:);
                                SS4_D_time(time,trial_ct)=id_node(time,:);
                            else
                                SS5_S_time(time,trial_ct)=is_node(time,:);
                                SS5_D_time(time,trial_ct)=id_node(time,:);
                            end
                        end
                        
                    else %task
                        s_il=sigmoid(il_node(time-1),beta_i,0);
                        s_ir=sigmoid(ir_node(time-1),beta_i,0);
                        s_ic=sigmoid(ic_node(time-1),beta_i,0);
                        s_ia=sigmoid(ia_node(time-1),beta_i,0);
                        input_u_to_i=w_iu*sum(sigmoid(u_field(time-1,:),beta,0));
                        if(NOISE)  %adding extra noise to nodes to try to reduce perfect performance
                            il_noise=noise_strength_fix*randn(1,1);
                            ir_noise=noise_strength_fix*randn(1,1);
                            ic_noise=noise_strength_fix*randn(1,1);
                            ia_noise=noise_strength_fix*randn(1,1);
                        else
                            il_noise=0;
                            ir_noise=0;
                            ic_noise=0;
                            ia_noise=0;
                        end %noise
                        
                        %TRYING THIS DIFFERENTLY FOR NOW (1/4/13)
                        %adding noisy input for nodes
                        if stimulus_A(time)>0
                            stimulus_A(time)=stimulus_A(time) + noisyinput*randn;
                        end
                        if stimulus_L(time)>0
                            stimulus_L(time)=stimulus_L(time) + noisyinput*randn;
                        end
                        if stimulus_R(time)>0
                            stimulus_R(time)=stimulus_R(time) + noisyinput*randn;
                        end
                        if stimulus_C(time)>0
                            stimulus_C(time)=stimulus_C(time) + noisyinput*randn;
                        end
                        
                        % dynamic h for nodes
                        h_il(time,:) = h_il(time-1,:) + (-h_slope*(h_il(time-1,:) + (h_rest + s_il*h_down)))/tau_h;
                        h_ir(time,:) = h_ir(time-1,:) + (-h_slope*(h_ir(time-1,:) + (h_rest + s_ir*h_down)))/tau_h;
                        h_ic(time,:) = h_ic(time-1,:) + (-h_slope*(h_ic(time-1,:) + (h_rest + s_ic*h_down)))/tau_h;
                        h_ia(time,:) = h_ia(time-1,:) + (-h_slope*(h_ia(time-1,:) + (h_rest + s_ia*h_down)))/tau_h;
                        
                        S_NC=stimulus_NC(:,time)';
                        S_CH=stimulus_CH(:,time)';
                        S=(S_NC*s_il + S_CH*s_ir);
                        input(time,:)=S;
                        
                        %equations
                        %noisy input rather than noise on nodes
                        il_node(time,:) = il_node(time-1,:) + (- il_node(time-1,:) + h_il(time-1,:) + w_ii*s_il + input_u_to_i*s_il + stimulus_L(1,time) - ci*(s_ir + s_ic + s_ia))/tau_i;
                        ir_node(time,:) = ir_node(time-1,:) + (- ir_node(time-1,:) + h_ir(time-1,:) + w_ii*s_ir + input_u_to_i*s_ir + stimulus_R(1,time) - ci*(s_il + s_ic + s_ia))/tau_i;
                        ic_node(time,:) = ic_node(time-1,:) + (- ic_node(time-1,:) + h_ic(time-1,:) + w_ii*s_ic                     + stimulus_C(1,time) - ci*(s_il + s_ir + s_ia))/tau_i;
                        ia_node(time,:) = ia_node(time-1,:) + (- ia_node(time-1,:) + h_ia(time-1,:) + w_ii*s_ia                     + stimulus_A(1,time) - ci*(s_il + s_ir + s_ic))/tau_i;
                        
                        if h_noise == 1
                            u_field(time,:) = u_field(time-1,:) + (- u_field(time-1,:) + h_u + h_noise_u(time-1,:) + S + w_ui*(s_il + s_ir + s_ic) + uu_conv           - uv_conv - w_uv_const*sum(f_v)+ upre_conv)/tau_u + u_noise;
                            v_field(time,:) = v_field(time-1,:) + (- v_field(time-1,:) + h_v + h_noise_v(time-1,:)                                 + vu_conv + vw_conv                                           )/tau_v + v_noise;
                            w_field(time,:) = w_field(time-1,:) + (- w_field(time-1,:) + h_w + h_noise_w(time-1,:) + w_wS*S                        + wu_conv + ww_conv - wv_conv - w_wv_const*sum(f_v)+ wpre_conv)/tau_w + w_noise;
                        else
                            u_field(time,:) = u_field(time-1,:) + (- u_field(time-1,:) + h_u + S + w_ui*(s_il + s_ir + s_ic) + uu_conv           - uv_conv - w_uv_const*sum(f_v)+ upre_conv)/tau_u + u_noise;
                            v_field(time,:) = v_field(time-1,:) + (- v_field(time-1,:) + h_v                                 + vu_conv + vw_conv                                           )/tau_v + v_noise;
                            w_field(time,:) = w_field(time-1,:) + (- w_field(time-1,:) + h_w + w_wS*S                        + wu_conv + ww_conv - wv_conv - w_wv_const*sum(f_v)+ wpre_conv)/tau_w + w_noise;
                        end
                        
                        theta_upre=u_field(time,:)>0;
                        theta_wpre=w_field(time,:)>0;
                        upre_field(time,:) = upre_field(time-1,:) + (theta_upre).*(- upre_field(time-1,:) + f_u)/tau_pre_u_build + (1-theta_upre).*(-upre_field(time-1,:))/tau_pre_u_decay + upre_noise;
                        wpre_field(time,:) = wpre_field(time-1,:) + (theta_wpre).*(- wpre_field(time-1,:) + f_w)/tau_pre_w_build + (1-theta_wpre).*(-wpre_field(time-1,:))/tau_pre_w_decay + wpre_noise;
                        
                        %counting peaks at each time step
                        peakDistance = 1;
                        supraThreshold = find(f_w >= 0.5)';
                        if size(supraThreshold, 1) < 2  %updated 8/16/13
                            nPeaksW = 0;
                        else
                            % the actual clustering
                            clusterIndices = clusterdata(supraThreshold, 'distance', 'euclidean', 'linkage',...
                                'single', 'cutoff', peakDistance, 'criterion', 'inconsistent');
                            % determining the centers of the clusters
                            nPeaksW = max(clusterIndices);
                        end
                  
                        if time>200
                            if il_node(time,:)>0 && ir_node(time,:)<0 && ia_node(time,:)<0 && ic_node(time,:)<0
                                look_L(time,:) = 1;
                            elseif ir_node(time,:)>0 && il_node(time,:)<0 && ia_node(time,:)<0 && ic_node(time,:)<0
                                look_R(time,:) = 1;
                            end
                            
                            if sum(look_L(time-20:time,:)) > 18
                                stable_L(time,:)=1;
                            end
                            if sum(look_R(time-20:time,:)) > 18
                                stable_R(time,:)=1;
                            end
                            
                            if stable_L(time,:)==0 && stable_L(time-1,:)==1
                                last_L=time;
                                n_looks=n_looks+1;
                            end
                            if stable_R(time,:)==0 && stable_R(time-1,:)==1
                                last_R=time;
                                n_looks=n_looks+1;
                            end
                            
                            if stable_L(time,:)==1 && stable_L(time-1,:)==0
                                if last_R > last_L
                                    switch_LtoR(time,:)=1;
                                    switches=switches+1;
                                end
                            end
                            if stable_R(time,:)==1 && stable_R(time-1,:)==0
                                if last_L > last_R
                                    switch_RtoL(time,:)=1;
                                    switches=switches+1;
                                end
                            end
                            
                            %%% summing activation in PF when L vs R looking node is above threshold
                            if il_node(time-1,:) > 0 && input_u_to_i > 0.5
                                ct_PFact_NC=ct_PFact_NC+1;
                            end
                            if ir_node(time-1,:) > 0 && input_u_to_i > 0.5
                                ct_PFact_CH=ct_PFact_CH+1;
                            end
                            
                            %%% counting the number of peaks at the end of each delay in PL
                            if time == blank_ends(blank)-1
                                peak_ct_PL(1,blank)=nPeaksW;
                                blank=blank+1;
                            end
                            
                        end %time>200
                    end % task
                     
                    if rct == total_runs
                        if track_act_U == 1
                            if SS == 1
                                SS1_U_time(time,trial_ct)=(sum((u_field(time,:)+abs(u_field(time,:)))./2));
                                SS1_peakU_time(time,trial_ct)=nPeaksU;
                            elseif SS == 2
                                SS2_U_time(time,trial_ct)=(sum((u_field(time,:)+abs(u_field(time,:)))./2));
                                SS2_peakU_time(time,trial_ct)=nPeaksU;
                            elseif SS == 3
                                SS3_U_time(time,trial_ct)=(sum((u_field(time,:)+abs(u_field(time,:)))./2));
                                SS3_peakU_time(time,trial_ct)=nPeaksU;
                            elseif SS == 4
                                SS4_U_time(time,trial_ct)=(sum((u_field(time,:)+abs(u_field(time,:)))./2));
                                SS4_peakU_time(time,trial_ct)=nPeaksU;
                            elseif SS == 5
                                SS5_U_time(time,trial_ct)=(sum((u_field(time,:)+abs(u_field(time,:)))./2));
                                SS5_peakU_time(time,trial_ct)=nPeaksU;
                            else
                                SS6_U_time(time,trial_ct)=(sum((u_field(time,:)+abs(u_field(time,:)))./2));
                                SS6_peakU_time(time,trial_ct)=nPeaksU;
                            end
                        end
                        
                        if track_act_W == 1
                            if SS == 1
                                SS1_W_time(time,trial_ct)=(sum((w_field(time,:)+abs(w_field(time,:)))./2));
                                SS1_peakW_time(time,trial_ct)=nPeaksW;
                            elseif SS == 2
                                SS2_W_time(time,trial_ct)=(sum((w_field(time,:)+abs(w_field(time,:)))./2));
                                SS2_peakW_time(time,trial_ct)=nPeaksW;
                            elseif SS == 3
                                SS3_W_time(time,trial_ct)=(sum((w_field(time,:)+abs(w_field(time,:)))./2));
                                SS3_peakW_time(time,trial_ct)=nPeaksW;
                            elseif SS == 4
                                SS4_W_time(time,trial_ct)=(sum((w_field(time,:)+abs(w_field(time,:)))./2));
                                SS4_peakW_time(time,trial_ct)=nPeaksW;
                            elseif SS == 5
                                SS5_W_time(time,trial_ct)=(sum((w_field(time,:)+abs(w_field(time,:)))./2));
                                SS5_peakW_time(time,trial_ct)=nPeaksW;
                            else
                                SS6_W_time(time,trial_ct)=(sum((w_field(time,:)+abs(w_field(time,:)))./2));
                                SS6_peakW_time(time,trial_ct)=nPeaksW;
                            end
                        end
                        
                    end %rct for track_act

                end %time loop
                
                n_timesteps=n_time;
                time=[1:n_time];
                
                upre_last=upre_field(n_timesteps,:);
                wpre_last=wpre_field(n_timesteps,:);
                
                stim1_pos(1,1:SS);
                stim2_pos(1,1:SS);
                
                if task == 1
                    steps_diff_above=find(id_node(1:end_time)>0);
                    steps_same_above=find(is_node(1:end_time)>0);
                    decision=[size(steps_diff_above) size(steps_same_above)];
                    
                    if decision(1) == 0 && decision(3) == 0
                        response=-2;
                        class = 0;
                        %CD_node_plot
                    elseif decision(1) > 50 && decision(3) > 50
                        %CD_node_plot
                        if decision(1) > decision(3)
                            response=-1;
                        else response=-0.1;
                        end
                        class = 0;
                    elseif decision(1) > decision(3)
                        response=1;
                        if change == 1;
                            class = 6; %HIT
                        else class = 8; %FA
                        end
                    else response=0;
                        if change == 0;
                            class = 5; %CR
                        else class = 7; %MISS
                        end
                    end
                else %task
                    steps_above_left=find(il_node>0);
                    steps_above_right=find(ir_node>0);
                    steps_above_away=find(ia_node>0);
                    steps_above_center=find(ic_node>0);
                    time_above_left=size(steps_above_left);
                    time_above_right=size(steps_above_right);
                    time_above_away=size(steps_above_away);
                    time_above_center=size(steps_above_center);
                    Total_Look=time_above_left+time_above_right;
                    %No change is left screen!
                    Change_Pref=time_above_right/Total_Look; %PRINTED
                    
                    %check for multiple nodes above threshold
                    left_look=0;
                    right_look=0;
                    away_look=0;
                    center_look=0;
                    no_look=0;
                    multi=0;
                    for lc=1:n_timesteps
                        if il_node(lc)>0 & ir_node(lc)<0 & ia_node(lc)<0 & ic_node(lc)<0
                            left_look=left_look+1;
                        elseif ir_node(lc)>0 & il_node(lc)<0 & ia_node(lc)<0 & ic_node(lc)<0
                            right_look=right_look+1;
                        elseif ia_node(lc)>0 & il_node(lc)<0 & ir_node(lc)<0 & ic_node(lc)<0
                            away_look=away_look+1;
                        elseif ic_node(lc)>0 & il_node(lc)<0 & ir_node(lc)<0 & ia_node(lc)<0
                            center_look=center_look+0;
                        elseif il_node(lc)<0 & ir_node(lc)<0 & ia_node(lc)<0 & ic_node(lc)<0
                            no_look=no_look+1;
                        else
                            multi=multi+1;
                        end
                    end
                    
                end  %task
                
                if STORAGE == 1
                    OutName = sprintf('%s.out',SimName);
                    OutFile = fopen(OutName,'a');
                    if STORE1 == 1
                        if task == 1
                            fprintf(OutFile,'AGE  E  SS  Trial  Peaks  CH  Resp  RT  S1  S2  S3  S4  S5  S6  S7  S8  S9  class\n');
                        else
                            fprintf(OutFile,'AGE  E  SS  Trial  TotTimeStep  ChPref  Left  Right  None  Multi  Switches  Looks  PFact_NC  PFact_CH  FirstPeaks LastPeaks  AvePeaks  MaxPeaks  MinPeaks\n');
                        end
                        STORE1=0;
                    end
                    if task == 1
                        fprintf(OutFile,'%i ',dct,rct,ct_stim1,trial_ct,peak_ct_CD,change,response,RT,stim1_pos,class);
                    else
                        fprintf(OutFile,'%i ',dct,rct,ct_stim1,trial_ct,Total_Look(1),Change_Pref,left_look,right_look,no_look,multi,switches,n_looks,ct_PFact_NC,ct_PFact_CH,peak_ct_PL(1),peak_ct_PL(26),mean(peak_ct_PL(1,1:26)),max(peak_ct_PL(1,1:26)),min(peak_ct_PL(1,1:26)));
                    end
                    fprintf(OutFile,'\n');
                    fclose(OutFile);
                else
                    if task == 1
                        fprintf('E %i SS %i Peaks %i CH %i Resp %i RT %i \n',rct,ct_stim1,nPeaksW,change,response,RT);
                    else
                        fprintf('E %i SS %i ChPref %i Switches %i \n',rct,ct_stim1,Change_Pref,switches);
                    end
                end %storage
                
                current_iteration=current_iteration+1;
                running_time=running_time+toc(tic_iteration)/60; % in minutes
                est_total_time=(total_iterations*running_time)/(current_iteration);
                fprintf('%i/%i      %0.0fh %0.0fm (completed) / %0.0fh %0.0fm (est. total)\n',current_iteration,total_iterations,...
                    floor(running_time/60),mod(running_time,60),floor(est_total_time/60),mod(est_total_time,60));
                
            end %trials in CD, SS in PL
        end  %SS in CD, trials in PL
    end %rct
end %dct
