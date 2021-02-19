
function chans = dc_getchans(type)

switch type
    
    case 'all'
        chans =  [{'AF3'  } {'AF4'  } {'AF7'  } ...
            {'AF8'  } {'AFz'  } {'C1'   } {'C2'   } {'C3'   } {'C4'   } {'C5'   } ...
            {'C6'   } {'CP1'  } {'CP2'  } {'CP3'  } {'CP4'  } {'CP5'  } {'CP6'  } ...
            {'CPz'  } {'Cz'   } {'F1'   } {'F2'   } {'F3'   } {'F4'   } {'F5'   } ...
            {'F6'   } {'F7'   } {'F8'   } {'FC1'  } {'FC2'  } {'FC3'  } {'FC4'  } ...
            {'FC5'  } {'FC6'  } {'FCz'  } {'FT7'  } {'FT8'  } {'Fz'   } {'O1'   } ...
            {'O2'   } {'Oz'   } {'P1'   } {'P2'   } {'P3'   } {'P4'   } {'P5'   } ...
            {'P6'   } {'P7'   } {'P8'   } {'PO3'  } {'PO4'  } {'PO7'  } {'PO8'  } ...
            {'POz'  } {'Pz'   } {'T7'   } {'T8'   } {'TP10' } {'TP7'  } {'TP8'  } ...
            {'TP9'  } {'RHEOG'} {'BVEOG'}] ;
    case 'noeyes'
        chans =  [{'AF3'  } {'AF4'  }...
            {'AFz'  } {'C1'   } {'C2'   } {'C3'   } {'C4'   } {'C5'   } ...
            {'C6'   } {'CP1'  } {'CP2'  } {'CP3'  } {'CP4'  } {'CP5'  } {'CP6'  } ...
            {'CPz'  } {'Cz'   } {'F1'   } {'F2'   } {'F3'   } {'F4'   } {'F5'   } ...
            {'F6'   } {'F7'   } {'F8'   } {'FC1'  } {'FC2'  } {'FC3'  } {'FC4'  } ...
            {'FC5'  } {'FC6'  } {'FCz'  } {'FT7'  } {'FT8'  } {'Fz'   } {'O1'   } ...
            {'O2'   } {'Oz'   } {'P1'   } {'P2'   } {'P3'   } {'P4'   } {'P5'   } ...
            {'P6'   } {'P7'   } {'P8'   } {'PO3'  } {'PO4'  } {'PO7'  } {'PO8'  } ...
            {'POz'  } {'Pz'   } {'T7'   } {'T8'   } {'TP10' } {'TP7'  } {'TP8'  } ...
            {'TP9'  } ] ;
    case 'topo'
        chans =   [{'AF7', 'AF3', 'AFz','AF4','AF8','F7','F5','F3','F1','Fz','F2','F4','F6', ....
            'F8','FT7','FC5','FC3','FC1','FCz','FC2','FC4','FC6','FT8','T7','C5','C3', ...
            'C1','Cz','C2','C4','C6','T8','TP7','CP5','CP3','CP1','CPz','CP2','CP4','CP6', ...
            'TP8','P7','P5', 'P3','P1','Pz','P2','P4','P6','P8','PO7','PO3','POz','PO4','PO8', ...
            'O1','Oz','O2'}];
        
    case 'even'
        chans =   [{'AFz','AF4','AF8','Fz','F2','F4','F6', ....
            'F8','FCz','FC2','FC4','FC6','FT8', ...
            'Cz','C2','C4','C6','T8','CPz','CP2','CP4','CP6', ...
            'TP8','Pz','P2','P4','P6','P8','POz','PO4','PO8', ...
            'Oz','O2'}];
    case 'odd'
        chans =   [{'AF7', 'AF3', 'AFz','F7','F5','F3','F1','Fz', ....
            'FT7','FC5','FC3','FC1','FCz','T7','C5','C3', ...
            'C1','Cz','TP7','CP5','CP3','CP1','CPz', ...
            'P7','P5', 'P3','P1','Pz','PO7','PO3','POz', ...
            'O1','Oz'}];
    case 'Reinhart-All'
        chans =  [{'AF3'  } {'AF4'  } {'AF7'  } ...
            {'AF8'  } {'AFz'  } {'C1'   } {'C2'   } {'C3'   } {'C4'   } {'C5'   } ...
            {'C6'   } {'CP1'  } {'CP2'  } {'CP3'  } {'CP4'  } {'CP5'  } {'CP6'  } ...
            {'CPz'  } {'Cz'   } {'F1'   } {'F2'   } {'F3'   } {'F4'   } {'F5'   } ...
            {'F6'   } {'F7'   } {'F8'   } {'FC1'  } {'FC2'  } {'FC3'  } {'FC4'  } ...
            {'FC5'  } {'FC6'  } {'FCz'  } {'FT7'  } {'FT8'  } {'Fz'   } {'O1'   } ...
            {'O2'   } {'Oz'   } {'P1'   } {'P2'   } {'P3'   } {'P4'   } {'P5'   } ...
            {'P6'   } {'P7'   } {'P8'   } {'PO3'  } {'PO4'  } {'PO7'  } {'PO8'  } ...
            {'POz'  } {'Pz'   } {'T7'   } {'T8'   } {'TP10' } {'TP7'  } {'TP8'  } ...
            {'TP9'  } {'RHEOG'} {'BVEOG'}] ;
    case 'RNET'
         chans = [ {'AF3'}    {'AF4'}    {'AF7'}    {'AF8'}    {'AFF1h'}    {'AFF2h'}    {'AFF3h'}    {'AFF4h'}    {'AFF5h'}    {'AFF6h'}    {'AFp1'}    {'AFp2'}    {'AFz'}    {'C1'}    {'C2'}    {'C3'}    {'C4'}    {'C5'}    {'C6'} ...
         {'CCP1h'}    {'CCP2h'}    {'CCP3h'}    {'CCP4h'}    {'CCP5h'}    {'CCP6h'}    {'CP1'}    {'CP2'}    {'CP3'}    {'CP4'}    {'CP5'}    {'CP6'}    {'CPP1h'}    {'CPP2h'}    {'CPP3h'}    {'CPP4h'}    {'CPP5h'}    {'CPP6h'} ...
         {'CPz'}    {'Cz'}    {'F1'}    {'F10'}    {'F2'}    {'F3'}    {'F4'}    {'F5'}    {'F6'}    {'F7'}    {'F8'}    {'F9'}    {'FC1'}    {'FC2'}    {'FC3'}    {'FC4'}    {'FC5'}    {'FC6'}    {'FCC1h'}    {'FCC2h'}    {'FCC3h'} ...
         {'FCC4h'}    {'FCC5h'}    {'FCC6h'}    {'FCz'}    {'FFC1h'}    {'FFC2h'}    {'FFC3h'}    {'FFC4h'}    {'FFC5h'}    {'FFC6h'}    {'FFT7h'}    {'FFT8h'}    {'FT10'}    {'FT7'}    {'FT8'}    {'FT9'}    {'FTT7h'}    {'FTT8h'} ...
         {'Fp1'}    {'Fp2'}    {'Fz'}    {'I1'}    {'I2'}    {'O1'}    {'O2'}    {'OI1h'}    {'OI2h'}    {'Oz'}    {'P1'}    {'P10'}    {'P2'}    {'P3'}    {'P4'}    {'P5'}    {'P6'}    {'P7'}    {'P8'}    {'P9'}    {'PO10'} ...
         {'PO3'}    {'PO4'}    {'PO7'}    {'PO8'}    {'PO9'}    {'POO1'}    {'POO10h'}    {'POO2'}    {'POO9h'}    {'POz'}    {'PPO10h'}    {'PPO1h'}    {'PPO2h'}    {'PPO3h'}    {'PPO4h'}    {'PPO5h'}    {'PPO6h'}    {'PPO9h'} ...
         {'Pz'}    {'T7'}    {'T8'}    {'TP10'}    {'TP7'}    {'TP8'}    {'TP9'}    {'TPP10h'}    {'TPP7h'}    {'TPP8h'}    {'TPP9h'}    {'TTP7h'}    {'TTP8h'}];

end