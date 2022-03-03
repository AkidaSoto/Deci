%Load the Participant translation sheet

Convert = tdfread('C:\Users\John\Documents\Data\1\participants.tsv');

mkdir(['C:\Users\John\Desktop\SourceData\1_new'])

%Loop Subj



Subjs = unique(cellfun(@(c) c(1:14),CleanDir('C:\Users\John\Desktop\SourceData\1'),'UniformOutput',false));

for Subj = 1:length(Subjs)
    %Load the Event Data
    hdr = ft_read_header(['C:\Users\John\Desktop\SourceData\1' filesep Subjs{Subj} '_task-ReinforcementLearning_eeg.set']);
    events = ft_read_event(['C:\Users\John\Desktop\SourceData\1' filesep Subjs{Subj} '_task-ReinforcementLearning_eeg.set']);
    events2 = events(1);

    Location = ismember(Convert.participant_id,{Subjs{Subj}(1:7)});
    Behv = Convert.Original_ID(Location);

    disp(['----------'])
    display(['Subject   ' num2str(Behv) '_sess' Subjs{Subj}(end)]);

    Behv= load(['C:\Users\John\Documents\Data\1\code\BEH' filesep num2str(Behv) '_sess' Subjs{Subj}(end) '_VVbeh']);

    trialcount =1;

    trialinfo = find(ismember({events.value},{'S  6','S  7','S 10','S 11'}));


    display(['Total # Trial is ' num2str(length(trialinfo))]);
    display(['Last Trial is at ' num2str(trialinfo(end))]);
    test = trialinfo(end);

    for event = 1:length(events)

        if trialcount < length(Behv.TRAIN) + 7


            if ismember(events(event).value,{'S  6','S  7','S 10','S 11'}) && trialcount < 7
                trialcount = trialcount +1;
                continue
                
            elseif trialcount >= 7

                if ismember(events(event).value,{'S  1','S  2'})

                    events2(end+1) = events(event);
                    events2(end).sample = events(event).sample - [500/2];
                    events2(end).value = 'S 20';

                    events2(end+1) = events(event);
                    events2(end).value = ['S ' num2str(Behv.TRAIN(trialcount-6,1) + 20)];
                end

                events2(end+1) = events(event);

                if ismember(events(event).value,{'S  4','S  5'})
                    events2(end+1) = events(event);
                    events2(end).value = 'S 30';

                    events2(end+1) = events(event);
                    events2(end).value = ['S ' num2str(~logical(Behv.TRAIN(trialcount-6,2)) + 31)];
                end

                if ismember(events(event).value,{'S  6','S  7','S 10','S 11'})
                    if ismember(events(event).value,{'S 10','S 11'})
                        events2(end+1) = events(event);

                        if ismember(events(event).value,{'S 10'}) ~= ~Behv.TRAIN(trialcount-6,4)
                            error('wrong')
                        elseif ismember(events(event).value,{'S 11'}) ~= Behv.TRAIN(trialcount-6,4)
                            error('wrong')
                        end


                        events2(end).value = ['S ' num2str(~logical(Behv.TRAIN(trialcount-6,4)) + 41)];
                        trialcount = trialcount +1;
                    else
                        events2(end+1) = events(event);
                        events2(end).value = ['S 43'];
                    end

                    events2(end+1) = events(event);
                    events2(end).value = 'S 40';
                    events2(end).sample = events(event).sample + [750/2];
                end

            end
        end
    end
    display(['End at ' num2str(events2(end).bvmknum)]);
    check = events2(end).bvmknum;

    if test ~= check
        warning('mismatch')
    end

    data = ft_read_data(['C:\Users\John\Desktop\SourceData\1' filesep Subjs{Subj} '_task-ReinforcementLearning_eeg.set']);
   

    hdr.orig.MarkerFile = ['C:\Users\John\Desktop\SourceData\1_new' filesep Subjs{Subj} '_task-ReinforcementLearning_eeg.set'];
    hdr.orig.DataFile = ['C:\Users\John\Desktop\SourceData\1_new' filesep Subjs{Subj} '_task-ReinforcementLearning_eeg.fdt'];
    
    ft_write_data(['C:\Users\John\Desktop\SourceData\1_new' filesep Subjs{Subj} '_task-ReinforcementLearning_eeg.fdt'],data,'header',hdr,'dataformat','brainvision_eeg','event',events2)
    

end

