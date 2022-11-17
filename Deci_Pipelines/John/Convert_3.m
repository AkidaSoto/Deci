%Load the Participant translation sheet

mkdir(['C:\Users\John\Desktop\SourceData\3_new'])

%Loop Subj

Subjs = unique(cellfun(@(c) c(1:14),CleanDir('C:\Users\John\Desktop\SourceData\3'),'UniformOutput',false));

for Subj = 1:length(Subjs)
    %Load the Event Data
    hdr = ft_read_header(['C:\Users\John\Desktop\SourceData\3' filesep Subjs{Subj} '_task-SimonConflict_eeg.set']);
    events = ft_read_event(['C:\Users\John\Desktop\SourceData\3' filesep Subjs{Subj} '_task-SimonConflict_eeg.set']);
    events = StandardizeEventMarkers(events);
    events2 = events(1);

    trialinfo = length(find(arrayfun(@(c) length(c.value) == 3 && str2num(c.value) < 106,events)));

    display(['Total # Trial is ' num2str(trialinfo)]);

    rspmade = 0;

    for event = 1:length(events)

        if length(events(event).value) == 3 && str2num(events(event).value) > 105


            if rspmade ~=0
                idx = find(ismember({events2.value},{'71'}),1,'last');
                events2 = events2(1:idx);
            end

            rspmade = 1;

            %Trial Start
            events2(end+1) = events(event);
            events2(end).sample = events(event).sample - [hdr.Fs*.25];
            events2(end).value = '70';

            events2(end+1) = events(event);
            events2(end).sample = events(event).sample;
            events2(end).value =  '10';


            events2(end+1) = events(event);
            events2(end).sample = events(event).sample;
            events2(end).value =  num2str(10 + str2num(events2(end).value(1)) );

            events2(end+1) = events(event);
            events2(end).sample = events(event).sample;
            events2(end).value =  num2str(20 + str2num(events2(end).value(2)) );

            events2(end+1) = events(event);
            events2(end).sample = events(event).sample;
            events2(end).value =  num2str(30 + str2num(events2(end).value(3)) );

        elseif length(events(event).value) == 3 && str2num(events(event).value) < 106

            if rspmade ~= 1
                idx = find(ismember({events2.value},{'71'}),1,'last');
                events2 = events2(1:idx);
                rspmade = 0;
                continue
            end

            rspmade = 2;

            if ismember(events(event).value,{'101','102','103','104'})
                events2(end+1) = events(event);
                events2(end).sample = events(event).sample;
                events2(end).value =  '40';
            else
                events2(end+1) = events(event);
            end

            if ismember(events(event).value,{'101','102'})

                events2(end+1) = events(event);
                events2(end).sample = events(event).sample;
                events2(end).value =  '41';

            elseif ismember(events(event).value,{'103','104'})

                events2(end+1) = events(event);
                events2(end).sample = events(event).sample;
                events2(end).value =  '42';

            end

            if ismember(events(event).value,{'101','103'})

                events2(end+1) = events(event);
                events2(end).sample = events(event).sample;
                events2(end).value =  '43';

            elseif ismember(events(event).value,{'102','104'})

                events2(end+1) = events(event);
                events2(end).sample = events(event).sample ;
                events2(end).value =  '44';

            end

        elseif ismember(events(event).value,{'6','7','8','9'})


            if rspmade ~= 2
                idx = find(ismember({events2.value},{'71'}),1,'last');
                events2 = events2(1:idx);
                rspmade =0;
                continue
            end

            rspmade = 0;
            events2(end+1) = events(event);

            if ismember(events(event).value,{'6','8','9'})
                events2(end+1) = events(event);
                events2(end).sample = events(event).sample;
                events2(end).value =  '50';
            end

            if ismember(events(event).value,{'8','9'})
                events2(end+1) = events(event);
                events2(end).sample = events(event).sample;
                events2(end).value =  '51';
            end

            events2(end+1) = events(event);
            events2(end).sample = events(event).sample + [hdr.Fs*.25];
            events2(end).value =  '71';

        end


    end


    data = ft_read_data(['C:\Users\John\Desktop\SourceData\3' filesep Subjs{Subj} '_task-SimonConflict_eeg.set']);
    hdr.orig.MarkerFile = ['C:\Users\John\Desktop\SourceData\3_new' filesep Subjs{Subj} '_task-SimonConflict_eeg.set'];
    hdr.orig.DataFile = ['C:\Users\John\Desktop\SourceData\3_new' filesep Subjs{Subj} '_task-SimonConflict_eeg.fdt'];

    ft_write_data(['C:\Users\John\Desktop\SourceData\3_new' filesep Subjs{Subj} '_task-SimonConflict_eeg.fdt'],data,'header',hdr,'dataformat','brainvision_eeg','event',events2)


end

