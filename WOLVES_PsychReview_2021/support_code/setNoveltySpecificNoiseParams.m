  function  setNoveltySpecificNoiseParams(sim, noise_ior_s, noise_wm_f )   
    sim.setElementParameters({'noise ior_s'}, {'amplitude'}, noise_ior_s ); % 1.3   
     %sim.setElementParameters({'atn_sa -> atn_c'}, {'amplitude'}, 6.3);
     %sim.setElementParameters({'wm_s -> wm_c'}, {'amplitude'}, 2);
     %sim.setElementParameters({'word -> wf'}, {'amplitude'}, 3.5); %4.25 
     %sim.setElementParameters({'hword -> word'}, {'amplitude'}, 0); 
     %testString = [char(D1) 'H' num2str(H*10) 'M' num2str(M*10) 'W' num2str(W*10) 'A' num2str(A*10)];
     %simNamePar = ['Hab' num2str(hab_Val(hab)*10) 'Atn' num2str(atn_Val(atn)*10) 'atnc2atnsa' num2str(at2at_Val(a2a)*10)];
     %simNamePar = ['wfChanges_conwmf0_hwmc1_fix48_xxx_12h' num2str(parVal2)];
    nFeatures = 2;
      for i = 1 : nFeatures
          n = num2str(i);
% % % %wolvesPaperPR diff from wfChanges.. varied fixations, habitaution, stop topdown from atn, now conf
% % % contrast to wm ihibition
          %sim.setElementParameters({['hwm_c' n ' -> wm_c' n]}, {'amplitude'}, 2.5); %H 1.5 for looking 5
          %sim.setElementParameters({['hwm_f' n ' -> wm_f' n]}, {'amplitude'}, 4.0); % F 0
          % sim.setElementParameters({['hcon_f' n ' -> con_f' n]}, {'amplitude'}, 0.3); % F 0
          %sim.setElementParameters({['wm_f' n ' -> wm_f' n]}, {'amplitudeExc'}, 19.0); % F 0
%          sim.setElementParameters({['wf' n ' -> atn_f' n]}, {'amplitude'}, 0);%3.5our average %A
          %sim.setElementParameters({['wf' n ' -> con_f' n]}, {'amplitude'}, 20);
%          %sim.setElementParameters({['wf' n ' -> atn_c' n]}, {'amplitude'}, -0.25);%1.3our average
%          %sim.setElementParameters({['hwf' n ' -> wf' n]}, {'amplitude'}, 4);% M %4 
%          %sim.setElementParameters({['wf' n ' -> wf' n]}, {'amplitudeExc'}, 35);%   W  
%          %sim.setElementParameters({['vis_f' n ' -> vis_f' n ' (global)']}, {'amplitude'}, -0.00005);%
%          %sim.setElementParameters({['vis_f' n ' -> con_f' n]}, {'amplitude'}, 0.5); % C
%          %sim.setElementParameters({['con_f' n ' -> atn_f' n]}, {'amplitude'}, 5);%1.
%          sim.setElementParameters({['con_f' n ' -> wm_f' n]}, {'amplitude'}, 0);%1.
%          %sim.setElementParameters({['atn_f' n ' -> atn_c' n]}, {'amplitude'}, 4);%1.
%          %sim.setElementParameters({['atn_c' n ' -> wf' n]}, {'amplitude'}, 0.7);%1.
          %sim.setElementParameters({['wm_f' n ' -> con_f' n]}, {'amplitude'}, -22);
%          %sim.setElementParameters({['atn_f' n ' -> wf' n]}, {'amplitude'}, 0);%1. noise wm_f1
            sim.setElementParameters({['noise wm_f' n]}, {'amplitude'}, noise_wm_f); %2
            %sim.setElementParameters({['noise vis_f' n]}, {'amplitude'}, 3);
      end