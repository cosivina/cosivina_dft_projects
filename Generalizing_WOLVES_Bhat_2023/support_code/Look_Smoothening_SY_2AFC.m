% DATA SMOOTHENING for Right and Left Look data or 2-choice looking 2AFC
for subject=1:numSubjects
    for trt=1:nTestTrials
        for side=1:2      
         if side == 1
             ldata = round(xsit_result.test(subject).historyLt(trt,:));
         else
             ldata = round(xsit_result.test(subject).historyRt(trt,:));
         end
         
            prevlooking=0;
            gaplen=0;
            for time=1:length(ldata)
                if (round(ldata(time)) == 0)
                    if prevlooking==0 %%continue adding off looking gap length
                        gaplen=gaplen+1;
                    else
                        prevlooking=0;
                        gaplen=0;
                    end
                else
                    if gaplen <= MIN_LOOK_DURATION
                        ldata(time-gaplen:time)=1;
                    end
                    gaplen=0;
                end
            end
            if side == 1
                xsit_result.test(subject).historyLt(trt,:)=ldata;
            else
                xsit_result.test(subject).historyRt(trt,:)=ldata;
            end    
        end
    end
    for tr=1:nTrainTrials
        for side=1:2      
         if side == 1
             ldata = round(xsit_result.train(subject).historyL(tr,:));
         else
             ldata = round(xsit_result.train(subject).historyR(tr,:));
         end
            prevlooking=0;
            gaplen=0;
            for time=1:length(ldata)
                if (round(ldata(time)) == 0)
                    if prevlooking==0 %%continue adding off looking gap length
                        gaplen=gaplen+1;
                    else
                        prevlooking=0;
                        gaplen=0;
                    end
                else
                    if gaplen <= MIN_LOOK_DURATION
                        ldata(time-gaplen:time)=1;
                    end
                    gaplen=0;
                end
            end
            if side == 1
                xsit_result.train(subject).historyL(tr,:)= ldata;
            else
                xsit_result.train(subject).historyR(tr,:)= ldata;
            end
        end
    end
end