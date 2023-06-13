train=[];test=[];
for subject=1:numSubjects
   try
       OutName1 = [simName,num2str(subject),'_train.mat'];              
       tempTrn=load(OutName1);
       [tempTrn(:).subject] = subject;
       train = [train; tempTrn];
       delete(OutName1);
    catch
       disp('Error on concatenating training file for subject number ');
       disp(subject);
       continue;
   end
end
for subject=1:numSubjects
    try
        OutName2 = [simName,num2str(subject),'_test.mat'];
        tempTst=load(OutName2);
        [tempTst(:).subject] = subject;
        test = [test; tempTst];
        delete(OutName2);
    catch
        disp('Error on concatenating test data for subject number ');
        disp(subject);
        continue;
    end
end
run_file = fileread('XSIT_Manual_run.m');sim_file = fileread('createComboSim.m');task_file = fileread([char(taskName),'.m']);
OutName = [simName,'results.mat'];
save(OutName,'train','test','sim','sim_seed','notes','run_file','sim_file','task_file');