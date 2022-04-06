%Load the Participant translation sheet

mkdir(['C:\Users\John\Desktop\SourceData\2_new'])

%Loop Subj

Convert = tdfread('C:\Users\John\Documents\Data\2\participants.tsv');
Convert = struct2table(Convert);

Convert2 = Convert(false,:);

Subjs = unique(cellfun(@(c) c(1:7),CleanDir('C:\Users\John\Desktop\SourceData\2'),'UniformOutput',false));

for Subj = 1:length(Subjs)
    %Load the Event Data
    hdr = ft_read_header(['C:\Users\John\Desktop\SourceData\2' filesep Subjs{Subj} '_task-ProbabilisticSelection_eeg.set']);
    events = ft_read_event(['C:\Users\John\Desktop\SourceData\2' filesep Subjs{Subj} '_task-ProbabilisticSelection_eeg.set']);

    if ~isempty(find(ismember({events.value},{'245'}),1))
        start = find(ismember({events.value},{'245'}),1);
    else
        start = 1;
    end

    events = events(start:find(ismember({events.value},{'104','94'}),1,'last'));
    events2 = events(1);

    trialinfo = find(ismember({events.value},{'104','94'}));

    display(['Subj is ' Subjs{Subj}]);
    if length(trialinfo) < 300
        continue
    end

    pos = find(ismember(table2cell(Convert(:,1)),Subjs{Subj}));
    Convert2(end+1,:) = Convert(pos,:);

    display(['Total # Trial is ' num2str(length(trialinfo))]);

    if Subj == length(Subjs)
        writetable(Convert2,'C:\Users\John\Documents\Data\2\participants2')

    end

    %     if exist(['C:\Users\John\Desktop\SourceData\2_new' filesep Subjs{Subj} '_task-ProbabilisticSelection_eeg.vmrk']) == 2
    %         continue
    %     end

    rspmade = 0;

    for event = 1:length(events)


        if ismember(events(event).value,arrayfun(@num2str,[10:21],'UniformOutput',false))

            if rspmade ~=0
                idx = find(ismember({events2.value},{'62'}),1,'last');
                events2 = events2(1:idx);
            end

            %Trial Start
            events2(end+1) = events(event);
            events2(end).sample = events(event).sample - [hdr.Fs*.25];
            events2(end).value = '61';

            %Stim Onset
            events2(end+1) = events(event);
            events2(end).sample = events(event).sample;
            events2(end).value = '30';

            condi = events(event);

            events2(end+1) = events(event);
            events2(end).sample = events(event).sample;
            events2(end).value = num2str(30+ceil([str2num(condi.value) - 9.9]/4));
            rspmade = 1;

        end

        if ismember(events(event).value,{'keypad1','keypad2'})

            if rspmade == 2
                continue
            elseif rspmade == 0
                idx = find(ismember({events2.value},{'62'}),1,'last');
                events2 = events2(1:idx);
                rspmade = 0;
                continue
            end

            rspmade = 2;
            events2(end+1) = events(event);
            events2(end).value = '40';

            events2(end+1) = events(event);
            events2(end).value = [num2str(find(ismember({'keypad2','keypad1'},events(event).value)) + 40)];

            events2(end+1) = events(event);
            events2(end).value = num2str(51 - [[rem([str2num(condi.value)-10],4) >= 2] == [find(ismember({'keypad2','keypad1'},events(event).value))-1]]);

        end

        if ismember(events(event).value,[{'104','94'}])

            if rspmade ~= 2
                idx = find(ismember({events2.value},{'62'}),1,'last');
                events2 = events2(1:idx);
                rspmade = 0;
                continue
            end
            
            rspmade = 0;
            events2(end+1) = events(event);

            if  ismember(events(event).value,{'104','94'})
                events2(end+1) = events(event);
                events2(end).value = '60';
                events2(end).sample = events(event).sample;
            end

            events2(end+1) = events(event);
            events2(end).value = '62';
            events2(end).sample = events(event).sample + [hdr.Fs*.25];
        end

    end


    data = ft_read_data(['C:\Users\John\Desktop\SourceData\2' filesep Subjs{Subj} '_task-ProbabilisticSelection_eeg.set']);
    hdr.orig.MarkerFile = ['C:\Users\John\Desktop\SourceData\2_new' filesep Subjs{Subj} '_task-ProbabilisticSelection_eeg.set'];
    hdr.orig.DataFile = ['C:\Users\John\Desktop\SourceData\2_new' filesep Subjs{Subj} '_task-ProbabilisticSelection_eeg.fdt'];

    ft_write_data(['C:\Users\John\Desktop\SourceData\2_new' filesep Subjs{Subj} '_task-ProbabilisticSelection_eeg.fdt'],data,'header',hdr,'dataformat','brainvision_eeg','event',events2)


end

