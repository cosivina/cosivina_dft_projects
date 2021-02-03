%         train=[];test=[];
%         simName= 'wmc2_3000Kachergis_Yu_Shiffrin_2012_';
%         numSubjects = 35;
%         
%          LCond=3;ECond=4;
%     
%         
%         %% Save Data          
%         for subject=1:numSubjects
%             try
%                 OutName1 = [simName,num2str(subject),'_train.mat'];              
%                 OutName2 = [simName,num2str(subject),'_test.mat'];
%                 tempTrn=load(OutName1);
%                 tempTst=load(OutName2);
%                 iteration=mod(subject,LCond*ECond);       
%                 if ismember(iteration, [1,4,7,10])%[3,6,9,0])%% ) %%if 3 late condition %% 6 latee [2,5,8,11]) %% 9 late [3,6,9,0]
%                     if  ismember(iteration, [1,5,9]) %% if 0 early condition
%                         train1 = [train1; tempTrn];
%                         test1 = [test1; tempTst];
%                     elseif ismember(iteration, [2,6,10])
%                         train10 = [train10; tempTrn];
%                         test10 = [test10; tempTst];
%                     elseif ismember(iteration, [3,7,11])
%                         train7 = [train7; tempTrn];
%                         test7 = [test7; tempTst];
%                     elseif ismember(iteration, [4,8,0])
%                         train4 = [train4; tempTrn];
%                         test4 = [test4; tempTst];
%                     end
%                 elseif ismember(iteration, [2,5,8,11])%[3,6,9,0])%% [1,4,7,10])
%                     if  ismember(iteration, [1,5,9]) %% if 0 early condition
%                         train5 = [train5; tempTrn];
%                         test5 = [test5; tempTst];
%                     elseif ismember(iteration, [2,6,10])
%                         train2 = [train2; tempTrn];
%                         test2 = [test2; tempTst];
%                     elseif ismember(iteration, [3,7,11])
%                         train11 = [train11; tempTrn];
%                         test11 = [test11; tempTst];
%                     elseif ismember(iteration, [4,8,0])
%                         train8 = [train8; tempTrn];
%                         test8 = [test8; tempTst];
%                     end
%                 elseif ismember(iteration, [3,6,9,0])%[2,5,8,11])%%% [1,4,7,10])
%                     if  ismember(iteration, [1,5,9]) %% if 0 early condition
%                         train9 = [train9; tempTrn];
%                         test9 = [test9; tempTst];
%                     elseif ismember(iteration, [2,6,10])
%                         train2 = [train2; tempTrn];
%                         test2 = [test2; tempTst];
%                     elseif ismember(iteration, [3,7,11])
%                         train11 = [train11; tempTrn];
%                         test11 = [test11; tempTst];
%                     elseif ismember(iteration, [4,8,0])
%                         train8 = [train8; tempTrn];
%                         test8 = [test8; tempTst];
%                     end    
%                 end    
%                     
%                 delete(OutName1);
%                 delete(OutName2);
%             catch
%                 disp('Error on concatenating subject number ');
%                 disp(subject);
%                 continue;
%             end
%         end
%         OutName = [simName,'results.mat'];
%         save(OutName,'train','test');


        train=[];test=[];
        simName= 'wmc2_3000Kachergis_Yu_Shiffrin_2012_';
        %numSubjects = 35;
        %% Save Data          
        for subject=17:740%numSubjects
            try
                OutName1 = [simName,num2str(subject),'_train.mat'];              
                OutName2 = [simName,num2str(subject),'_test.mat'];
                tempTrn=load(OutName1);
                tempTst=load(OutName2);
                [tempTrn(:).subject] = subject;
                train = [train; tempTrn];
                test = [test; tempTst];
                delete(OutName1);
                delete(OutName2);
            catch
                disp('Error on concatenating subject number ');
                disp(subject);
                continue;
            end
        end
        OutName = [simName,'results.mat'];
        save(OutName,'train','test');















