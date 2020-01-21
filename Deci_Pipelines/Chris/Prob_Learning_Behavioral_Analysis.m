%Prob_Learning_Behavioral_Analysis

load('/Volumes/Seagate Portable Drive/ProcessedData/PreProc/Prob_EEG_001BU_info.mat');
load('');

cond_info = data.condinfo{1,2}(1:360,1:6);


Mrk =[events.value]';

Block_Starts = find(Mrk==100);
Block_Ends = find(Mrk==101);

ab_total = 0;
cd_total = 0;
ef_total = 0;

ab_correct = 0;
cd_correct = 0;
ef_correct = 0;

for i = 1:6
    
    if i == 1
        Block{1} = Mrk(1:Block_Starts(1));
    else
        Block{i} = Mrk(Block_Starts(i-1):Block_Starts(i));
    end
    
    temp=[Block{i}];
    AB_count{i} = find(temp == 14);
    ab_opt{i} = find(
    
    CD_count{i} = find(temp == 15);
    cd_opt{i} = 
    
    EF_count{i} = find(temp == 16);
    
    ab_total = ab_total + numel(AB_count{i});
    cd_total = cd_total + numel(CD_count{i});
    ef_total = ef_total + numel(EF_count{i});

    
end