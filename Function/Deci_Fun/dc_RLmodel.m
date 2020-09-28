function [out] = dc_RLmodel(Deci,info,dat,params)
% dc_RLmodel is ran once for all block per state (ex:  All blocks for Reward) per subject.

% 7 types of models are currently available

%1.) Classic model with 2 parameters, alpha and beta
%2.) Classic model with 2 alphas for - and + Pe
%3.) Classic model with 2 alphas for - and + Pe, letting intial Q vary
%4.) Classic model with 2 alphas for - and + Q
%5.) Classic model with additional d parameters
%6.) Actor Critic model with 4 parameters, alphaA, AlphaC, d and beta
%7.) Hybrid AC/Q model with 6
    %Classic model with alphaQ
    %AC model with alphaA and alphaC
    %with outside parameters d,c and beta

% Model fit works by using fmincon to find the local minima (highest probability of state | action)
%from the randomized starting free parameter. The model is rurun multiple
%times to find the global minima.

%out - outgoing parameters of best fit, cell array 1x5

%QL - Free Parameters
%Ev - Expected value [q,w,v]
%Pe - Prediction Error
%PseudoR - Model Fitness
%AIC - Model comparison fitness

%Deci - Configuration File parameters
%.Folder.Analysis - Folder where all analysis data is saved
%.SubjectList - List of subject
%.Analysis.CondTitle - List of conditions

%info - Unique parameters for this instance of function run
%.subject_list - current subject index
%.Cond - current condition index

%dat - Output of step 3 from Deci, fieldtrip data file
%Artifact trial rejected
% Function is usually called from Analyzor line 172 so you can read
% how Analyzor loads an example Step 3 datafile.

%params - Free variables depending on model
%.State - Event code for all possible conditions
%.Actions - Event code for all possible Actions
%.Reward - Event code for all possible Rewards
%.Value - Value for each rewardxstate combo, cell format
%.Reps - number of iterations to run randomized starting free
%parameters
%.ModelNum - number array of models to run
%.Start - Starting Intial Expected Value for q,w,v
%.Seed - rng seed;

%Ex:
%params.States = [20 21 23 24];
%params.Value  = {[20 10] [10 0] [0 -10] [-10 -20]};
%params.Actions = [31 32];
%params.Reward  = [51 52];
%params.Start = {[0 0],[.01 .01],[0]};


trialtypes = dat.events;

%Basic Data Maintence for finding block numbers
% Deci segments blocks by negative numbers
if any(any(trialtypes < 0))
    trialtypes(:,find(trialtypes(1,:) < 1)) = ceil(trialtypes(:,find(trialtypes(1,:) < 1)));
    blocknumbers = sort(unique(trialtypes(:,find(trialtypes(1,:) < 1))),'descend');
else
    blocknumbers = -1;
end

Values = find(ismember(params.States,trialtypes));

%Sorting data into blocks
for blk = 1:length(blocknumbers)
    
    blkmrk = sum(ismember(trialtypes,[blocknumbers(blk)]),2) ==  1;
    
    Actmrk = logical(sum(ismember(trialtypes(blkmrk,:),[params.Actions]),1));
    Actions{blk} = trialtypes(blkmrk,Actmrk);
    
    Rewmrk = logical(sum(ismember(trialtypes(blkmrk,:),[params.Reward]),1));
    Rewards{blk} = trialtypes(blkmrk,Rewmrk);
    
    for rew = 1:length(params.Reward)
        Rewards{blk}(Rewards{blk} == params.Reward(rew)) = params.Value{Values}(rew);
    end
    
end

% Loop for running models N times with randomized starting conditinos

LLE = nan([7 params.Reps]);
LLE2 = nan([7 params.Reps]);


%1.) Classic model with 2 parameters, alpha and beta
%2.) Classic model with 2 alphas for - and + Pe
%3.) Classic model with 2 alphas for - and + Pe, letting intial Q vary
%4.) Classic model with 2 alphas for - and + Q
%5.) Classic model with additional d parameters
%6.) Actor Critic model with 4 parameters, alphaA, AlphaC, d and beta
%7.) Hybrid AC/Q model with 6
    %Classic model with alphaQ
    %AC model with alphaA and alphaC
    %with outside parameters d,c and beta

rng(params.Seed)

%try

A_rep = [];
R_rep = [];


A_rep{1,1} = permute([Actions(:)],[2 3 1]);
R_rep{1,1} = permute([Rewards(:)],[2 3 1]);


for Dim2 = 1:size(A_rep,2)
    for Dim1 = 1:size(A_rep,1)
        
        
        Actions = A_rep{Dim1,Dim2};
        Rewards = R_rep{Dim1,Dim2};
        
        % Loop for running models N times with randomized starting conditinos
        
        LLE = nan([7 Deci.Analysis.Extra.QL.Reps]);
        LLE2 = nan([7 Deci.Analysis.Extra.QL.Reps]);
        
        
        %1.) Classic model with 2 parameters, alpha and beta
        %2.) Classic model with 2 alphas for - and + Pe
        %3.) Classic model with 2 alphas for - and + Pe, letting intial Q vary
        %4.) Classic model with 2 alphas for - and + Q
        %5.) Classic model with additional d parameters
        %6.) Actor Critic model with 4 parameters, alphaA, AlphaC, d and beta
        %7.) Hybrid AC/Q model with 6
        %Classic model with alphaQ
        %AC model with alphaA and alphaC
        %with outside parameters d,c and beta
        
        rng(Deci.Analysis.Extra.QL.Seed)
        
        %try
        
        ColNames = {};
        SheetNames = {};
        
        
        for init = 1:Deci.Analysis.Extra.QL.Reps
            
            disp(['Dim1 #' num2str(Dim1) ' Dim2 #' num2str(Dim2) ' init #' num2str(init)]);
            
            %1.) Classic model with 2 parameters, alpha and beta
            
            if ismember(1, Deci.Analysis.Extra.QL.ModelNum)
                Model1.LB = [0 1e-6];
                Model1.UB = [1 30];
                Model1.init =rand(1,length(Model1.LB)).*(Model1.UB-Model1.LB)+Model1.LB;
                
                
                [Value{1,init}] = ...
                    fmincon(@(x) SimpleQ(Deci.Analysis.Extra.QL.Start{1},Deci.Analysis.Extra.QL.Actions,Actions,Rewards,x(1),x(2)),...
                    [Model1.init],[],[],[],[],[Model1.LB],[Model1.UB],[],...
                    optimset('TolX', 0.00001, 'TolFun', 0.00001, 'MaxFunEvals', 9e+9, 'Algorithm', 'interior-point','Display','off'));
                
                [LLE(1,init),out{1,init}] = SimpleQ(Deci.Analysis.Extra.QL.Start{1},Deci.Analysis.Extra.QL.Actions,Actions,Rewards, Value{1,init}(1),Value{1,init}(2));
                LLE(1,init) = -LLE(1,init); %change fmin to fmax     
                LLE2(1,init) = aicbic(LLE(1,init),[2]);
                
                if init == 1
                    ColNames{1} = {'1_Classic_Alp' '1_Classic_Beta'};
                    SheetNames{1} = ('1_Classic_Model');
                end
                
            end
            
            %2.)  Classic model with 2 alphas for - and + Pe
            if ismember(2, Deci.Analysis.Extra.QL.ModelNum)
                Model2.LB = [0 0 1e-6];
                Model2.UB = [1 1 30];
                Model2.init =rand(1,length(Model2.LB)).*(Model2.UB-Model2.LB)+Model2.LB;
                
                [Value{2,init}] = ...
                    fmincon(@(x) Model2Pes(Deci.Analysis.Extra.QL.Start{1},Deci.Analysis.Extra.QL.Actions,Actions,Rewards,x(1),x(2),x(3)),...
                    [Model2.init],[],[],[],[],[Model2.LB],[Model2.UB],[],...
                    optimset('TolX', 0.00001, 'TolFun', 0.00001, 'MaxFunEvals', 9e+9, 'Algorithm', 'interior-point','Display','off'));
                
                [LLE(2,init),out{2,init}] = Model2Pes(Deci.Analysis.Extra.QL.Start{1},Deci.Analysis.Extra.QL.Actions,Actions,Rewards, Value{2,init}(1),Value{2,init}(2),Value{2,init}(3));
                LLE(2,init) = -LLE(2,init); %change fmin to fmax    
                LLE2(2,init) = aicbic(LLE(2,init),[3]);
                
                if init == 1
                    ColNames{2} = {'2_2PEs_AlpR' '2_2PEs_AlpP' '2_2PEs_Beta'};
                    SheetNames{2} = ('2_2PEs');
                end
            end
            
            %3.)  Classic model with 2 alphas for - and + Pe, letting intial Q vary
            if ismember(3,Deci.Analysis.Extra.QL.ModelNum)
                
                %Define Lower and Upper bounds of free parameters for FminCon
                Model3.LB = [min([Deci.Analysis.Extra.QL.Value{:}]) min([Deci.Analysis.Extra.QL.Value{:}]) 0 0 1e-6];
                Model3.UB = [max([Deci.Analysis.Extra.QL.Value{:}]) max([Deci.Analysis.Extra.QL.Value{:}]) 1 1 30];
                Model3.init =rand(1,length(Model3.LB)).*(Model3.UB-Model3.LB)+Model3.LB;
                
                [Value{3,init}] = ...
                    fmincon(@(x)  Model2Pes([x(1) x(2)],Deci.Analysis.Extra.QL.Actions,Actions,Rewards,x(3),x(4),x(5)),...
                    [Model3.init],[],[],[],[],[Model3.LB],[Model3.UB],[],...
                    optimset('TolX', 0.00001, 'TolFun', 0.00001, 'MaxFunEvals', 9e+9, 'Algorithm', 'interior-point','Display','off'));
                
                [LLE(3,init),out{3,init}] = Model2Pes([Value{3,init}(1) Value{3,init}(2)],Deci.Analysis.Extra.QL.Actions,Actions,Rewards, Value{3,init}(3),Value{3,init}(4),Value{3,init}(5));
                LLE(3,init) = -LLE(3,init); %change fmin to fmax    
                LLE2(3,init) = aicbic(LLE(3,init),[5]);
                
                if init == 1
                    ColNames{3} = {'3_2PEsM_QOpt' '3_2PEsM_QWst' '3_2PEsM_AlpR' '3_2PEsM_AlpP' '3_2PesM_Beta'};
                    SheetNames{3} = ('3_2PEsM');
                end
            end
            
            %4.)  Classic model with 2 alphas for - and + Pe
            if ismember(4, Deci.Analysis.Extra.QL.ModelNum)
                Model4.LB = [0 0 1e-6];
                Model4.UB = [1 1 30];
                Model4.init =rand(1,length(Model4.LB)).*(Model4.UB-Model4.LB)+Model4.LB;
                
                [Value{4,init}] = ...
                    fmincon(@(x) Model2Q(Deci.Analysis.Extra.QL.Start{1},Deci.Analysis.Extra.QL.Actions,Actions,Rewards,x(1),x(2),x(3)),...
                    [Model4.init],[],[],[],[],[Model4.LB],[Model4.UB],[],...
                    optimset('TolX', 0.00001, 'TolFun', 0.00001, 'MaxFunEvals', 9e+9, 'Algorithm', 'interior-point','Display','off'));
                
                [LLE(4,init),out{4,init}] = Model2Q(Deci.Analysis.Extra.QL.Start{1},Deci.Analysis.Extra.QL.Actions,Actions,Rewards, Value{4,init}(1),Value{4,init}(2),Value{4,init}(3));
                LLE(4,init) = -LLE(4,init); %change fmin to fmax    
                LLE2(4,init) = aicbic(LLE(4,init),[3]);
                
                if init == 1
                    ColNames{4} = {'4_2Qs_AlpR' '4_2Qs_AlpP' '4_2Qs_Beta'};
                    SheetNames{4} = ('4_2Qs');
                end
            end
            
            %5.) Classic model with additional d parameters
            
            if ismember(5, Deci.Analysis.Extra.QL.ModelNum)
                Model5.LB = [0 1e-6 0];
                Model5.UB = [1 30 1];
                Model5.init =rand(1,length(Model5.LB)).*(Model5.UB-Model5.LB)+Model5.LB;
                
                
                [Value{5,init}] = ...
                    fmincon(@(x) NormalizedQ(Deci.Analysis.Extra.QL.Start{1},Deci.Analysis.Extra.QL.Actions,Actions,Rewards,x(1),x(2),x(3)),...
                    [Model5.init],[],[],[],[],[Model5.LB],[Model5.UB],[],...
                    optimset('TolX', 0.00001, 'TolFun', 0.00001, 'MaxFunEvals', 9e+9, 'Algorithm', 'interior-point','Display','off'));
                
                [LLE(5,init),out{5,init}] = NormalizedQ(Deci.Analysis.Extra.QL.Start{1},Deci.Analysis.Extra.QL.Actions,Actions,Rewards, Value{5,init}(1),Value{5,init}(2),Value{5,init}(3));
                LLE(5,init) = -LLE(5,init); %change fmin to fmax    
                LLE2(5,init) = aicbic(LLE(5,init),[3]);
                
                
                if init == 1
                    ColNames{5} = {'5_ClassicD_Alp' '5_ClassicD_Beta' '5_Classic_d'};
                    SheetNames{5} = ('5_ClassicD');
                end
            end
            
            %6.) Actor Critic model with 4 parameters, alphaA, AlphaC, d and beta
            
            if ismember(6,Deci.Analysis.Extra.QL.ModelNum)
                
                %Define Lower and Upper bounds of free parameters for FminCon
                Model6.LB = [0 0 1e-6 0];
                Model6.UB = [1 1 30 1];
                Model6.init =rand(1,length(Model6.LB)).*(Model6.UB-Model6.LB)+Model6.LB;
                
                [Value{6,init}] = ...
                    fmincon(@(x)  ActorCritic(Deci.Analysis.Extra.QL.Start{2},Deci.Analysis.Extra.QL.Start{3},Deci.Analysis.Extra.QL.Actions,Actions,Rewards,x(1),x(2),x(3),x(4)),...
                    [Model6.init],[],[],[],[],[Model6.LB],[Model6.UB],[],...
                    optimset('TolX', 0.00001, 'TolFun', 0.00001, 'MaxFunEvals', 9e+9, 'Algorithm', 'interior-point','Display','off'));
                
                [LLE(6,init),out{6,init}] = ActorCritic(Deci.Analysis.Extra.QL.Start{2},Deci.Analysis.Extra.QL.Start{3},Deci.Analysis.Extra.QL.Actions,Actions,Rewards, Value{6,init}(1),Value{6,init}(2),Value{6,init}(3),Value{6,init}(4));
                LLE(6,init) = -LLE(6,init); %change fmin to fmax    
                LLE2(6,init) = aicbic(LLE(6,init),[4]);
                
                if init == 1
                    ColNames{6} = {'6_ActorCritic_AlpC' '6_ActorCritic_AlpA' '6_ActorCritic_Beta' '6_ActorCritic_d'};
                    SheetNames{6} = ('6_ActorCritic');
                end
            end
            
            %7.) Hybrid AC/Q model with 6
            %Classic model with alphaQ
            %AC model with alphaA and alphaC
            %with outside parameters d,c and beta
            
            if ismember(7,Deci.Analysis.Extra.QL.ModelNum)
                
                %Define Lower and Upper bounds of free parameters for FminCon
                Model7.LB = [0 0 0 1e-6 0 0];
                Model7.UB = [1 1 1 30 1 1];
                Model7.init =rand(1,length(Model7.LB)).*(Model7.UB-Model7.LB)+Model7.LB;
                
                [Value{7,init}] = ...
                    fmincon(@(x)  HybridModel(Deci.Analysis.Extra.QL.Start{1},Deci.Analysis.Extra.QL.Start{2},Deci.Analysis.Extra.QL.Start{3},Deci.Analysis.Extra.QL.Actions,Actions,Rewards,x(1),x(2),x(3),x(4),x(5),x(6)),...
                    [Model7.init],[],[],[],[],[Model7.LB],[Model7.UB],[],...
                    optimset('TolX', 0.00001, 'TolFun', 0.00001, 'MaxFunEvals', 9e+9, 'Algorithm', 'interior-point','Display','off'));
                
                [LLE(7,init),out{7,init}] = HybridModel(Deci.Analysis.Extra.QL.Start{1},Deci.Analysis.Extra.QL.Start{2},Deci.Analysis.Extra.QL.Start{3},Deci.Analysis.Extra.QL.Actions,Actions,Rewards, Value{7,init}(1),Value{7,init}(2),Value{7,init}(3),Value{7,init}(4),Value{7,init}(5),Value{7,init}(6));
                LLE(7,init) = -LLE(7,init); %change fmin to fmax    
                LLE2(7,init) = aicbic(LLE(7,init),[6]);
                
                
                if init == 1
                    ColNames{7} = {'7_Hybrid_AlpQ' '7_Hybrid_AlpA' '7_Hybrid_AlpC' '7_Hybrid_Beta' '7_Hybrid_d' '7_Hybrid_c'};
                    SheetNames{7} = ('7_Hybrid');
                end
            end
            
                        %6.) Actor Critic model with 4 parameters, alphaA, AlphaC, d and beta
            
            if ismember(8,Deci.Analysis.Extra.QL.ModelNum)
                
                %Define Lower and Upper bounds of free parameters for FminCon
                Model8.LB = [0 0 1e-6];
                Model8.UB = [1 1 30];
                Model8.init =rand(1,length(Model8.LB)).*(Model8.UB-Model8.LB)+Model8.LB;
                
                [Value{8,init}] = ...
                    fmincon(@(x)  ActorCritic2(Deci.Analysis.Extra.QL.Start{2},Deci.Analysis.Extra.QL.Start{3},Deci.Analysis.Extra.QL.Actions,Actions,Rewards,x(1),x(2),x(3)),...
                    [Model8.init],[],[],[],[],[Model8.LB],[Model8.UB],[],...
                    optimset('TolX', 0.00001, 'TolFun', 0.00001, 'MaxFunEvals', 9e+9, 'Algorithm', 'interior-point','Display','off'));
                
                [LLE(8,init),out{8,init}] = ActorCritic2(Deci.Analysis.Extra.QL.Start{2},Deci.Analysis.Extra.QL.Start{3},Deci.Analysis.Extra.QL.Actions,Actions,Rewards, Value{8,init}(1),Value{8,init}(2),Value{8,init}(3));
                LLE(8,init) = -LLE(8,init); %change fmin to fmax    
                LLE2(8,init) = aicbic(LLE(8,init),[3]);
                
                if init == 1
                    ColNames{8} = {'8_ActorCritic_AlpC' '8_ActorCritic_AlpA' '8_ActorCritic_Beta'};
                    SheetNames{8} = ('8_ActorCritic');
                end
            end
            
            %7.) Hybrid AC/Q model with 6
            %Classic model with alphaQ
            %AC model with alphaA and alphaC
            %with outside parameters d,c and beta
            
            if ismember(9,Deci.Analysis.Extra.QL.ModelNum)
                
                %Define Lower and Upper bounds of free parameters for FminCon
                Model9.LB = [0 0 0 1e-6 0];
                Model9.UB = [1 1 1 30 1];
                Model9.init =rand(1,length(Model9.LB)).*(Model9.UB-Model9.LB)+Model9.LB;
                
                [Value{9,init}] = ...
                    fmincon(@(x)  HybridModel2(Deci.Analysis.Extra.QL.Start{1},Deci.Analysis.Extra.QL.Start{2},Deci.Analysis.Extra.QL.Start{3},Deci.Analysis.Extra.QL.Actions,Actions,Rewards,x(1),x(2),x(3),x(4),x(5)),...
                    [Model9.init],[],[],[],[],[Model9.LB],[Model9.UB],[],...
                    optimset('TolX', 0.00001, 'TolFun', 0.00001, 'MaxFunEvals', 9e+9, 'Algorithm', 'interior-point','Display','off'));
                
                [LLE(9,init),out{9,init}] = HybridModel2(Deci.Analysis.Extra.QL.Start{1},Deci.Analysis.Extra.QL.Start{2},Deci.Analysis.Extra.QL.Start{3},Deci.Analysis.Extra.QL.Actions,Actions,Rewards, Value{9,init}(1),Value{9,init}(2),Value{9,init}(3),Value{9,init}(4),Value{9,init}(5));
                LLE(9,init) = -LLE(9,init); %change fmin to fmax    
                LLE2(9,init) = aicbic(LLE(9,init),[5]);
                
                
                if init == 1
                    ColNames{9} = {'9_Hybrid_AlpQ' '9_Hybrid_AlpA' '9_Hybrid_AlpC' '9_Hybrid_Beta' '9_Hybrid_c'};
                    SheetNames{9} = ('9_Hybrid');
                end
            end
            
        end
        % catch
        %     %sometimes the model crashes from the given randomized starting free
        %     %parameers, so we try to catch it. Still unsure why, but I usually just
        %     %rerun the code again if it happens.
        %     error(['init ' num2str(Model3.init) ' ' num2str(Model2.init) ' ' num2str(Model1.init) ' ']);
        % end
        
        
        [Best,I] = nanmin(LLE,[],2);
        [Best2,I2] = nanmin(LLE2,[],2);
        
        
        for m = 1:length(I2)
            if ismember(m,Deci.Analysis.Extra.QL.ModelNum)
                BestMod(m) = Value(m,I2(m));
                BestOut(m) = out(m,I2(m));
                PseudoR(m) = 1 - [Best(m)/[length(cat(1,Actions{:}))*log(.05)]];
                
            end
        end
        
        warning('off', 'MATLAB:MKDIR:DirectoryExists');
        
        % Let's save some data
        
        if Deci.Analysis.Extra.QL.saveexcel
            
            mkdir([Deci.Folder.Version filesep 'Plot']);
            if sum([Dim1,Dim2]) == 2
                delete([Deci.Folder.Version filesep 'Plot' filesep 'ModelOutputs' '.xlsx']);
            end
            
            if sum([Dim1,Dim2]) == 2
                writecell(['Subject' 'Cond' SheetNames(:)'],[ Deci.Folder.Version filesep 'Plot' filesep 'ModelOutputs' '.xlsx'],'Sheet','Summary_LLE','Range','A1');
                writecell(['Subject' 'Cond' SheetNames(:)'],[ Deci.Folder.Version filesep 'Plot' filesep 'ModelOutputs' '.xlsx'],'Sheet','Summary_AIC','Range','A1');
                writecell(['Subject' 'Cond' SheetNames(:)'],[ Deci.Folder.Version filesep 'Plot' filesep 'ModelOutputs' '.xlsx'],'Sheet','Summary_PseudoR','Range','A1');
            end
            
            
            writecell([Deci.SubjectList(Dim1) Deci.Analysis.CondTitle(Dim2) mat2cell([Best'],[1],ones([length([Best]) 1]))],[ Deci.Folder.Version filesep 'Plot' filesep 'ModelOutputs' '.xlsx'],'Sheet','Summary_LLE','Range',['A' num2str(sub2ind(size(A_rep),Dim1,Dim2)+1)]);
            writecell([Deci.SubjectList(Dim1) Deci.Analysis.CondTitle(Dim2) mat2cell([Best2'],[1],ones([length([Best2]) 1]))],[ Deci.Folder.Version filesep 'Plot' filesep 'ModelOutputs' '.xlsx'],'Sheet','Summary_AIC','Range',['A' num2str(sub2ind(size(A_rep),Dim1,Dim2)+1)])
            writecell([Deci.SubjectList(Dim1) Deci.Analysis.CondTitle(Dim2) mat2cell([PseudoR],[1],ones([length([PseudoR]) 1]))],[ Deci.Folder.Version filesep 'Plot' filesep 'ModelOutputs' '.xlsx'],'Sheet','Summary_PseudoR','Range',['A' num2str(sub2ind(size(A_rep),Dim1,Dim2)+1)])
            
        end
        
        if Deci.Analysis.Extra.QL.saveparams
            for m = 1:length(I2)
                if ismember(m,Deci.Analysis.Extra.QL.ModelNum)
                    
                    ModelDeci.Analysis.Extra.QL = fields(BestOut{m});
                    
                    for MP = 1:length(ModelDeci.Analysis.Extra.QL)
                        
                        param = BestOut{m}.(ModelDeci.Analysis.Extra.QL{MP});
                        
                        mkdir([Deci.Folder.Analysis filesep 'Extra' filesep SheetNames{m} filesep ModelDeci.Analysis.Extra.QL{MP} filesep Deci.SubjectList{info.subject_list}])
                        save([Deci.Folder.Analysis filesep 'Extra' filesep SheetNames{m} filesep ModelDeci.Analysis.Extra.QL{MP} filesep Deci.SubjectList{info.subject_list}  filesep Deci.Analysis.CondTitle{info.Cond}],'param');
                        
                    end
                    
                    
                end
            end
        end
        
        
    end
end



function [LLE,out] = SimpleQ(starting_q,possible_actions,a,r,alp,beta)

% This is the Classic model with no change.

%What are the possible Actions?
possible_actions;

% Function can run model on multiple blocks at once.

for subj = 1:size(a,1)
    
    for cond = 1:size(a,2)
        
        for block = 1:size(a,3)
            
            %get the intial starting expected value
            qout{subj,cond,block}(1,:) = starting_q;
            q{subj,cond,block}(1) = starting_q(1);
            
            %Now loop through the actions and find reward, Pe, Q(t+1) and P
            for Act = 1:length(a{subj,cond,block})
                
                %Which action is it?
                which = find(a{subj,cond,block}(Act) == possible_actions);
                
                %What is the expected value of that action?
                qout{subj,cond,block}(Act,which);
                
                %What is the Pe?
                PE{subj,cond,block}(Act) = (r{subj,cond,block}(Act) - qout{subj,cond,block}(Act,which));
                
                %What is the Probability of that action?
                
                P{subj,cond,block}(Act) = exp(qout{subj,cond,block}(Act,which)/beta)/ ...
                    [exp(qout{subj,cond,block}(Act,which)/beta) + exp(qout{subj,cond,block}(Act,find(a{subj,cond,block}(Act) ~= possible_actions))/beta)];
                
                % sometimes the probabilty is 1, but the
                % formula outputs a nan and so we have to
                % compensate for that.
                if isnan(P{subj,cond,block}(end))
                    P{subj,cond,block}(end) = qout{subj,cond,block}(Act,which) > qout{subj,cond,block}(Act,find(a{subj,cond,block}(Act) ~= possible_actions));
                end
                
                %update the expected value of the chosen action
                qout{subj,cond,block}(Act+1,which) = qout{subj,cond,block}(Act,which) + alp*PE{subj,cond,block}(Act);
                
                %keep the unchosen action's expected value the same.
                qout{subj,cond,block}(Act+1,find(a{subj,cond,block}(Act) ~= possible_actions))  = qout{subj,cond,block}(Act,find(a{subj,cond,block}(Act) ~= possible_actions));
           
                % Save Value;
                q{block}(Act+1) = qout{subj,cond,block}(Act,which);
            end
            
            %remove last, because it's not needed.
            q_next{block} = q{block}(2:end);
            q{block} = q{block}(1:end-1);
        end
        
    end
end

%fmincon finds the the local minima, so we will sum up all the
%probabilities and get the negative value of that. In theory, the
%higher the P, the higher the model were able to predict that the
%probability of the choice.
LLE = -sum(log([P{:}]));

if ~isfinite(LLE)
    LLE = -sum(numel(a{:})*log(.01));
end

out.q = [q{:}];
out.q_next = [q_next{:}];
out.P = [P{:}];
out.Pe = [PE{:}];

end

function [LLE,out] = Model2Pes(starting_q,possible_actions,a,r,alpP,alpN,beta)

% This is the model with that has two different alphas, for
% negative and positive Pe


%What are the possible Actions?
possible_actions;

% Function can run model on multiple blocks at once.
for subj = 1:size(a,1)
    
    for cond = 1:size(a,2)
        
        for block = 1:size(a,3)
            
            %get the intial starting expected value
            qout{subj,cond,block}(1,:) = starting_q;
            
            %Now loop through the actions and find reward, Pe, Q(t+1) and P
            for Act = 1:length(a{subj,cond,block})
                
                %Which action is it?
                which = find(a{subj,cond,block}(Act) == possible_actions);
                
                %What is the expected value of that action?
                q{block}(Act) = qout{subj,cond,block}(Act,which);
                
                %What is the Pe?
                PE{subj,cond,block}(Act) = (r{subj,cond,block}(Act) - qout{subj,cond,block}(Act,which));
                
                %What is the Probability of that action?
                
                P{subj,cond,block}(Act) = exp(qout{subj,cond,block}(Act,which)/beta)/ ...
                    [exp(qout{subj,cond,block}(Act,which)/beta) + exp(qout{subj,cond,block}(Act,find(a{subj,cond,block}(Act) ~= possible_actions))/beta)];
                
                % sometimes the probabilty is 1, but the
                % formula outputs a nan and so we have to
                % compensate for that.
                if isnan(P{subj,cond,block}(end))
                    P{subj,cond,block}(end) = qout{subj,cond,block}(Act,which) > qout{subj,cond,block}(Act,find(a{subj,cond,block}(Act) ~= possible_actions));
                end
                
                %update the expected value of the chosen action
                % depending on the Pe's abs val
                if PE{subj,cond,block}(Act) > 0
                    qout{subj,cond,block}(Act+1,which) = qout{subj,cond,block}(Act,which) + alpP*PE{subj,cond,block}(Act);
                else
                    qout{subj,cond,block}(Act+1,which) = qout{subj,cond,block}(Act,which) + alpN*PE{subj,cond,block}(Act);
                end
                
                %keep the unchosen action's expected value the same.
                qout{subj,cond,block}(Act+1,find(a{subj,cond,block}(Act) ~= possible_actions))  = qout{subj,cond,block}(Act,find(a{subj,cond,block}(Act) ~= possible_actions));
            end
            
            %remove last, because it's not needed.
            q{block} = qout{block}(1:end-1);
        end
    end
end
%fmincon finds the the local minima, so we will sum up all the
%probabilities and get the negative value of that. In theory, the
%higher the P, the higher the model were able to predict that the
%probability of the choice.
LLE = -sum(log([P{:}]));

if ~isfinite(LLE)
    LLE = -sum(numel(a{:})*log(.01));
end

out.q = qout(:,1);
out.P = P;
out.Pe = PE;
end


function [LLE,out] = Model2Q(starting_q,possible_actions,a,r,alpP,alpN,beta)

% This is the model with that has two different alphas, for
% negative and positive Q (expected values)


%What are the possible Actions?
possible_actions;

% Function can run model on multiple blocks at once.
for subj = 1:size(a,1)
    
    for cond = 1:size(a,2)
        
        for block = 1:size(a,3)
            
            %get the intial starting expected value
            qout{subj,cond,block}(1,:) = starting_q;
            
            %Now loop through the actions and find reward, Pe, Q(t+1) and P
            for Act = 1:length(a{subj,cond,block})
                
                %Which action is it?
                which = find(a{subj,cond,block}(Act) == possible_actions);
                
                %What is the expected value of that action?
                qout{subj,cond,block}(Act,which);
                
                %What is the Pe?
                PE{subj,cond,block}(Act) = (r{subj,cond,block}(Act) - qout{subj,cond,block}(Act,which));
                
                %What is the Probability of that action?
                
                P{subj,cond,block}(Act) = exp(qout{subj,cond,block}(Act,which)/beta)/ ...
                    [exp(qout{subj,cond,block}(Act,which)/beta) + exp(qout{subj,cond,block}(Act,find(a{subj,cond,block}(Act) ~= possible_actions))/beta)];
                
                % sometimes the probabilty is 1, but the
                % formula outputs a nan and so we have to
                % compensate for that.
                if isnan(P{subj,cond,block}(end))
                    P{subj,cond,block}(end) = qout{subj,cond,block}(Act,which) > qout{subj,cond,block}(Act,find(a{subj,cond,block}(Act) ~= possible_actions));
                end
                
                %update the expected value of the chosen action
                % depending on the Pe's abs val
                if r{subj,cond,block}(Act) > 0
                    qout{subj,cond,block}(Act+1,which) = qout{subj,cond,block}(Act,which) + alpP*PE{subj,cond,block}(Act);
                else
                    qout{subj,cond,block}(Act+1,which) = qout{subj,cond,block}(Act,which) + alpN*PE{subj,cond,block}(Act);
                end
                
                %keep the unchosen action's expected value the same.
                qout{subj,cond,block}(Act+1,find(a{subj,cond,block}(Act) ~= possible_actions))  = qout{subj,cond,block}(Act,find(a{subj,cond,block}(Act) ~= possible_actions));
            end
            
            %remove last, because it's not needed.
            qout{subj,cond,block} = qout{subj,cond,block}(1:end-1,:);
        end
    end
end
%fmincon finds the the local minima, so we will sum up all the
%probabilities and get the negative value of that. In theory, the
%higher the P, the higher the model were able to predict that the
%probability of the choice.
LLE = -sum(log([P{:}]));

if ~isfinite(LLE)
    LLE = -sum(numel(a{:})*log(.01));
end

out.q = qout;
out.P = P;
out.Pe = PE;
end



function [LLE,out] = NormalizedQ(starting_q,possible_actions,a,r,alp,beta,d)

% This is the Classic model with no change.

%What are the possible Actions?
possible_actions;

% Function can run model on multiple blocks at once.
for subj = 1:size(a,1)
    
    for cond = 1:size(a,2)
        
        for block = 1:size(a,3)
            
            
            %get the intial starting expected value
            qout{subj,cond,block}(1,:) = starting_q;
            
            %Now loop through the actions and find reward, Pe, Q(t+1) and P
            for Act = 1:length(a{subj,cond,block})
                
                %Which action is it?
                which = find(a{subj,cond,block}(Act) == possible_actions);
                
                %What is the expected value of that action?
                qout{subj,cond,block}(Act,which);
                
                %determine the outcome
                if r{subj,cond,block}(Act) > 0
                    outcome = 1 - d;
                elseif r{subj,cond,block}(Act) < 0
                    outcome = -d;
                elseif r{subj,cond,block}(Act) == 0
                    outcome = 0;
                else
                    error('nan in reward structure?')
                end
                
                %What is the Pe?
                PE{subj,cond,block}(Act) = (outcome - qout{subj,cond,block}(Act,which));
                
                %What is the Probability of that action?
                
                P{subj,cond,block}(Act) = exp(qout{subj,cond,block}(Act,which)/beta)/ ...
                    [exp(qout{subj,cond,block}(Act,which)/beta) + exp(qout{subj,cond,block}(Act,find(a{subj,cond,block}(Act) ~= possible_actions))/beta)];
                
                % sometimes the probabilty is 1, but the
                % formula outputs a nan and so we have to
                % compensate for that.
                
                if isnan(P{subj,cond,block}(end))
                    P{subj,cond,block}(end) = qout{subj,cond,block}(Act,which) > qout{subj,cond,block}(Act,find(a{subj,cond,block}(Act) ~= possible_actions));
                end
                
                
                
                %update the expected value of the chosen action
                qout{subj,cond,block}(Act+1,which) = qout{subj,cond,block}(Act,which) + alp*PE{subj,cond,block}(Act);
                
                %keep the unchosen action's expected value the same.
                qout{subj,cond,block}(Act+1,find(a{subj,cond,block}(Act) ~= possible_actions))  = qout{subj,cond,block}(Act,find(a{subj,cond,block}(Act) ~= possible_actions));
            end
            
            %remove last, because it's not needed.
            qout{subj,cond,block} = qout{subj,cond,block}(1:end-1,:);
        end
    end
end

%fmincon finds the the local minima, so we will sum up all the
%probabilities and get the negative value of that. In theory, the
%higher the P, the higher the model were able to predict that the
%probability of the choice.
LLE = -sum(log([P{:}]));

if ~isfinite(LLE)
    LLE = -sum(numel(a{:})*log(.01));
end

out.q = qout;
out.P = P;
out.Pe = PE;
end


function [LLE,out] = ActorCritic(starting_w,starting_v,possible_actions,a,r,alpC,alpA,beta,d)

% This is the basic ActorCritic Model


%What are the possible Actions?
possible_actions;

% Function can run model on multiple blocks at once.
for subj = 1:size(a,1)
    
    for cond = 1:size(a,2)
        
        for block = 1:size(a,3)
            
            
            %get the intial starting expected value
            vout{subj,cond,block}(1) = starting_v;
            
            %get the intial starting actor values
            wout{subj,cond,block}(1,:) = starting_w;
            
            %Now loop through the actions and find reward, Pe, Q(t+1) and P
            for Act = 1:length(a{subj,cond,block})
                
                %Which action is it?
                which = find(a{subj,cond,block}(Act) == possible_actions);
                
                % Which action isn't it?
                isnt = find(a{subj,cond,block}(Act) ~= possible_actions);
                
                %What is the expected value?
                vout{subj,cond,block}(Act);
                
                %determine the outcome
                if r{subj,cond,block}(Act) > 0
                    outcome = 1 - d;
                elseif r{subj,cond,block}(Act) < 0
                    outcome = -d;
                elseif r{subj,cond,block}(Act) == 0
                    outcome = 0;
                else
                    error('nan in reward structure?')
                end
                
                %What is the Pe?
                PE{subj,cond,block}(Act) = (outcome - vout{subj,cond,block}(Act));
                
                
                %Update the critic
                vout{subj,cond,block}(Act+1) = vout{subj,cond,block}(Act) + alpC*PE{subj,cond,block}(Act);
                
                %Update the Actor
                wout{subj,cond,block}(Act+1,which) = wout{subj,cond,block}(Act,which) + alpA*PE{subj,cond,block}(Act);
                
                %keep the unchosen Actor the same.
                wout{subj,cond,block}(Act+1,isnt)  = wout{subj,cond,block}(Act,isnt);

                
                %Since P calcuation is outside the loop, we will save the
                %outs into a seperate variable delinated by chosen and not
                %action type
                
                w_chosen{subj,cond,block}(Act,1) = wout{subj,cond,block}(Act,which);
                w_chosen{subj,cond,block}(Act,2) = wout{subj,cond,block}(Act,isnt);
                
            end
            
            %remove last, because it's not needed.
            vout{subj,cond,block} = vout{subj,cond,block}(1:end-1,:);
            wout{subj,cond,block} = wout{subj,cond,block}(1:end-1,:);
            
            
            %Normalize the Actors
            w_ph = w_chosen{subj,cond,block};
            
            w_chosen{subj,cond,block}(:,1) = w_ph(:,1)./[abs(w_ph(:,1)) + abs(w_ph(:,2))];
            w_chosen{subj,cond,block}(:,2) = w_ph(:,2)./[abs(w_ph(:,1)) + abs(w_ph(:,2))];
            
            %What is the Probability of that action?
            
            P{subj,cond,block} = exp(w_chosen{subj,cond,block}(:,1)/beta)./ ...
                [exp(w_chosen{subj,cond,block}(:,1)/beta) + exp(w_chosen{subj,cond,block}(:,2)/beta)];
            
%             % sometimes the probabilty is 1, but the
%             % formula outputs a nan and so we have to
%             % compensate for that.
            if any(isnan(P{subj,cond,block}))
                P{subj,cond,block}(isnan(P{subj,cond,block})) = w_chosen{subj,cond,block}(isnan(P{subj,cond,block}),1) > w_chosen{subj,cond,block}(isnan(P{subj,cond,block}),2);
            end
            
            
            
        end
    end
end

%fmincon finds the the local minima, so we will sum up all the
%probabilities and get the negative value of that. In theory, the
%higher the P, the higher the model were able to predict that the
%probability of the choice.
LLE = -sum(log([cat(1,P{:})]));

if ~isfinite(LLE)
    LLE = -sum(numel(a{:})*log(.01));
end

out.v = vout;
out.w = wout;
out.w_c = w_chosen;
out.P = P;
out.Pe = PE;
end

function [LLE,out] = HybridModel(starting_q,starting_w,starting_v,possible_actions,a,r,alpQ,alpC,alpA,beta,d,c)

% This is the basic ActorCritic Model


%What are the possible Actions?
possible_actions;

% Function can run model on multiple blocks at once.
for subj = 1:size(a,1)
    
    for cond = 1:size(a,2)
        
        for block = 1:size(a,3)
            
            %get the intial starting expected value
            qout{subj,cond,block}(1,:) = starting_q;
            
            %get the intial starting critic value
            vout{subj,cond,block}(1) = starting_v;
            
            %get the intial starting actor values
            wout{subj,cond,block}(1,:) = starting_w;
            
            %Now loop through the actions and find reward, Pe, Q(t+1) and P
            for Act = 1:length(a{subj,cond,block})
                
                %Which action is it?
                which = find(a{subj,cond,block}(Act) == possible_actions);
                
                % Which action isn't it?
                isnt = find(a{subj,cond,block}(Act) ~= possible_actions);
                
                %What is the expected value of that action?
                qout{subj,cond,block}(Act,which);
                
                
                %determine the outcome
                if r{subj,cond,block}(Act) > 0
                    outcome = 1 - d;
                elseif r{subj,cond,block}(Act) < 0
                    outcome = -d;
                elseif r{subj,cond,block}(Act) == 0
                    outcome = 0;
                else
                    error('nan in reward structure?')
                end
                
                % ------------------- QL -------------------
                %What is the QL Pe?
                PE_QL{subj,cond,block}(Act) = (outcome - qout{subj,cond,block}(Act,which));
                
                
                %update the expected value of the chosen action
                qout{subj,cond,block}(Act+1,which) = qout{subj,cond,block}(Act,which) + alpQ*PE_QL{subj,cond,block}(Act);
                
                %keep the unchosen action's expected value the same.
                qout{subj,cond,block}(Act+1,isnt)  = qout{subj,cond,block}(Act,isnt);
                
                
                %Since P calcuation is outside the loop, we will save the
                %outs into a seperate variable delinated by chosen and not
                %action type
                
                q_chosen{subj,cond,block}(Act,1) = qout{subj,cond,block}(Act,which);
                q_chosen{subj,cond,block}(Act,2) = qout{subj,cond,block}(Act,isnt);
                
                
                
                % ------------------- AC -------------------
                
                %What is the AC Pe?
                PE_AC{subj,cond,block}(Act) = (outcome - vout{subj,cond,block}(Act));
                
                
                %Update the critic
                vout{subj,cond,block}(Act+1) = vout{subj,cond,block}(Act) + alpC*PE_AC{subj,cond,block}(Act);
                
                %Update the Actor
                wout{subj,cond,block}(Act+1,which) = wout{subj,cond,block}(Act,which) + alpA*PE_AC{subj,cond,block}(Act);
                
                %keep the unchosen Actor the same.
                wout{subj,cond,block}(Act+1,isnt)  = wout{subj,cond,block}(Act,isnt);
                
                
                w_chosen{subj,cond,block}(Act,1) = wout{subj,cond,block}(Act,which);
                w_chosen{subj,cond,block}(Act,2) = wout{subj,cond,block}(Act,isnt);

            end
            
            %remove last, because it's not needed.
            qout{subj,cond,block} = qout{subj,cond,block}(1:end-1,:);
            vout{subj,cond,block} = vout{subj,cond,block}(1:end-1,:);
            wout{subj,cond,block} = wout{subj,cond,block}(1:end-1,:);
            
            %Normalize the Actors
            w_ph = w_chosen{subj,cond,block};
            
            w_chosen{subj,cond,block}(:,1) = w_ph(:,1)./[abs(w_ph(:,1)) + abs(w_ph(:,2))];
            w_chosen{subj,cond,block}(:,2) = w_ph(:,2)./[abs(w_ph(:,1)) + abs(w_ph(:,2))];
            
            % ------------------- Hybrid -------------------
            
            %Determine the Model Contributions
            H{subj,cond,block}(:,1) = [1-c]*w_chosen{subj,cond,block}(:,1) + [c]*q_chosen{subj,cond,block}(:,1);
            H{subj,cond,block}(:,2) = [1-c]*w_chosen{subj,cond,block}(:,2) + [c]*q_chosen{subj,cond,block}(:,2);
            
            %What is the Probability of that action?
            
            P{subj,cond,block} = exp(H{subj,cond,block}(:,1)/beta)./ ...
                [exp(H{subj,cond,block}(:,1)/beta) + exp(H{subj,cond,block}(:,2)/beta)];
            
            % sometimes the probabilty is 1, but the
            % formula outputs a nan and so we have to
            % compensate for that.
            if any(isnan(P{subj,cond,block}))
                P{subj,cond,block}(isnan(P{subj,cond,block})) = H{subj,cond,block}(isnan(P{subj,cond,block}),1) > H{subj,cond,block}(isnan(P{subj,cond,block}),2);
            end
            
            
            
            
        end
    end
end

%fmincon finds the the local minima, so we will sum up all the
%probabilities and get the negative value of that. In theory, the
%higher the P, the higher the model were able to predict that the
%probability of the choice.
LLE = -sum(log([cat(1,P{:})]));


if ~isfinite(LLE)
    LLE = -sum(numel(a{:})*log(.01));
end

out.q = qout;
out.v = vout;
out.w = wout;
out.P = P;
out.QPe = PE_QL;
out.QAC = PE_AC;
end

function [LLE,out] = ActorCritic2(starting_w,starting_v,possible_actions,a,r,alpC,alpA,beta)

% This is the basic ActorCritic Model


%What are the possible Actions?
possible_actions;

% Function can run model on multiple blocks at once.
for subj = 1:size(a,1)
    
    for cond = 1:size(a,2)
        
        for block = 1:size(a,3)
            
            
            %get the intial starting expected value
            vout{subj,cond,block}(1) = starting_v;
            %get the intial starting actor values
            wout{subj,cond,block}(1,:) = starting_w;
            
            v{subj,cond,block}(1) = starting_v(1);
            %get the intial starting actor values
            w{subj,cond,block}(1) = starting_w(1);
            
            
            %Now loop through the actions and find reward, Pe, Q(t+1) and P
            for Act = 1:length(a{subj,cond,block})
                
                %Which action is it?
                which = find(a{subj,cond,block}(Act) == possible_actions);
                
                % Which action isn't it?
                isnt = find(a{subj,cond,block}(Act) ~= possible_actions);
                
                %What is the expected value?
                vout{subj,cond,block}(Act);
                
                %determine the outcome
%                 if r{subj,cond,block}(Act) > 0
%                     outcome = 1 - d;
%                 elseif r{subj,cond,block}(Act) < 0
%                     outcome = -d;
%                 elseif r{subj,cond,block}(Act) == 0
%                     outcome = 0;
%                 else
%                     error('nan in reward structure?')
%                 end
                
                %What is the Pe?
                PE{subj,cond,block}(Act) = (r{subj,cond,block}(Act) - vout{subj,cond,block}(Act));
                
                
                %Update the critic
                vout{subj,cond,block}(Act+1) = vout{subj,cond,block}(Act) + alpC*PE{subj,cond,block}(Act);
                
                %Update the Actor
                wout{subj,cond,block}(Act+1,which) = wout{subj,cond,block}(Act,which) + alpA*PE{subj,cond,block}(Act);
                
                %keep the unchosen Actor the same.
                wout{subj,cond,block}(Act+1,isnt)  = wout{subj,cond,block}(Act,isnt);

                
                %Since P calcuation is outside the loop, we will save the
                %outs into a seperate variable delinated by chosen and not
                %action type
                
                w_chosen{subj,cond,block}(Act,1) = wout{subj,cond,block}(Act,which);
                w_chosen{subj,cond,block}(Act,2) = wout{subj,cond,block}(Act,isnt);
                
                w{subj,cond,block}(Act+1) = wout{subj,cond,block}(Act,which);
                v{subj,cond,block}(Act+1) = vout{subj,cond,block}(Act);
            end
            
            %remove last, because it's not needed.
            v_next{subj,cond,block} = v{subj,cond,block}(2:end);
            w_next{subj,cond,block} = w{subj,cond,block}(2:end);
            
            v{subj,cond,block} = v{subj,cond,block}(1:end-1);
            w{subj,cond,block} = w{subj,cond,block}(1:end-1);
            
            %Normalize the Actors
            w_ph = w_chosen{subj,cond,block};
            
            w_chosen{subj,cond,block}(:,1) = w_ph(:,1)./[abs(w_ph(:,1)) + abs(w_ph(:,2))];
            w_chosen{subj,cond,block}(:,2) = w_ph(:,2)./[abs(w_ph(:,1)) + abs(w_ph(:,2))];
            
            %What is the Probability of that action?
            
            P{subj,cond,block} = exp(w_chosen{subj,cond,block}(:,1)/beta)./ ...
                [exp(w_chosen{subj,cond,block}(:,1)/beta) + exp(w_chosen{subj,cond,block}(:,2)/beta)];
            
%             % sometimes the probabilty is 1, but the
%             % formula outputs a nan and so we have to
%             % compensate for that.
            if any(isnan(P{subj,cond,block}))
                P{subj,cond,block}(isnan(P{subj,cond,block})) = w_chosen{subj,cond,block}(isnan(P{subj,cond,block}),1) > w_chosen{subj,cond,block}(isnan(P{subj,cond,block}),2);
            end
            
            
            
        end
    end
end

%fmincon finds the the local minima, so we will sum up all the
%probabilities and get the negative value of that. In theory, the
%higher the P, the higher the model were able to predict that the
%probability of the choice.
LLE = -sum(log([cat(1,P{:})]));

if ~isfinite(LLE)
    LLE = -sum(numel(a{:})*log(.01));
end

out.v = [v{:}];
out.v_next = [v_next{:}];
out.w = [w{:}];
out.w_next = [w_next{:}];

out.P = cat(1,P{:});
out.Pe = [PE{:}];
end

function [LLE,out] = HybridModel2(starting_q,starting_w,starting_v,possible_actions,a,r,alpQ,alpC,alpA,beta,c)

% This is the basic ActorCritic Model


%What are the possible Actions?
possible_actions;

% Function can run model on multiple blocks at once.
for subj = 1:size(a,1)
    
    for cond = 1:size(a,2)
        
        for block = 1:size(a,3)
            
            %get the intial starting expected value
            qout{subj,cond,block}(1,:) = starting_q;
            
            %get the intial starting critic value
            vout{subj,cond,block}(1) = starting_v;
            
            %get the intial starting actor values
            wout{subj,cond,block}(1,:) = starting_w;
            
            
            %get the intial starting expected value
            q{subj,cond,block}(1) = starting_q(1);
            
            %get the intial starting critic value
            v{subj,cond,block}(1) = starting_v(1);
            
            %get the intial starting actor values
            w{subj,cond,block}(1) = starting_w(1);
            
            %Now loop through the actions and find reward, Pe, Q(t+1) and P
            for Act = 1:length(a{subj,cond,block})
                
                %Which action is it?
                which = find(a{subj,cond,block}(Act) == possible_actions);
                
                % Which action isn't it?
                isnt = find(a{subj,cond,block}(Act) ~= possible_actions);
                
                %What is the expected value of that action?
                qout{subj,cond,block}(Act,which);
                
                
                %determine the outcome
%                 if r{subj,cond,block}(Act) > 0
%                     outcome = 1 - d;
%                 elseif r{subj,cond,block}(Act) < 0
%                     outcome = -d;
%                 elseif r{subj,cond,block}(Act) == 0
%                     outcome = 0;
%                 else
%                     error('nan in reward structure?')
%                 end
                
                % ------------------- QL -------------------
                %What is the QL Pe?
                PE_QL{subj,cond,block}(Act) = (r{subj,cond,block}(Act) - qout{subj,cond,block}(Act,which));
                
                
                %update the expected value of the chosen action
                qout{subj,cond,block}(Act+1,which) = qout{subj,cond,block}(Act,which) + alpQ*PE_QL{subj,cond,block}(Act);
                
                %keep the unchosen action's expected value the same.
                qout{subj,cond,block}(Act+1,isnt)  = qout{subj,cond,block}(Act,isnt);
                
                
                %Since P calcuation is outside the loop, we will save the
                %outs into a seperate variable delinated by chosen and not
                %action type
                
                q_chosen{subj,cond,block}(Act,1) = qout{subj,cond,block}(Act,which);
                q_chosen{subj,cond,block}(Act,2) = qout{subj,cond,block}(Act,isnt);
                
                
                
                % ------------------- AC -------------------
                
                %What is the AC Pe?
                PE_AC{subj,cond,block}(Act) = (r{subj,cond,block}(Act) - vout{subj,cond,block}(Act));
                
                
                %Update the critic
                vout{subj,cond,block}(Act+1) = vout{subj,cond,block}(Act) + alpC*PE_AC{subj,cond,block}(Act);
                
                %Update the Actor
                wout{subj,cond,block}(Act+1,which) = wout{subj,cond,block}(Act,which) + alpA*PE_AC{subj,cond,block}(Act);
                
                %keep the unchosen Actor the same.
                wout{subj,cond,block}(Act+1,isnt)  = wout{subj,cond,block}(Act,isnt);
                
                
                w_chosen{subj,cond,block}(Act,1) = wout{subj,cond,block}(Act,which);
                w_chosen{subj,cond,block}(Act,2) = wout{subj,cond,block}(Act,isnt);

                q{subj,cond,block}(Act+1) =  qout{subj,cond,block}(Act,which);
                w{subj,cond,block}(Act+1) =  wout{subj,cond,block}(Act,which);
                v{subj,cond,block}(Act+1) =  vout{subj,cond,block}(Act);
            end
            
            %remove last, because it's not needed.
            
            q_next{subj,cond,block} = q{subj,cond,block}(2:end);
            v_next{subj,cond,block} = v{subj,cond,block}(2:end);
            w_next{subj,cond,block} = w{subj,cond,block}(2:end);
            
            q{subj,cond,block} = q{subj,cond,block}(1:end-1);
            v{subj,cond,block} = v{subj,cond,block}(1:end-1);
            w{subj,cond,block} = w{subj,cond,block}(1:end-1);
            
            
            %Normalize the Actors
            w_ph = w_chosen{subj,cond,block};
            
            w_chosen{subj,cond,block}(:,1) = w_ph(:,1)./[abs(w_ph(:,1)) + abs(w_ph(:,2))];
            w_chosen{subj,cond,block}(:,2) = w_ph(:,2)./[abs(w_ph(:,1)) + abs(w_ph(:,2))];
            
            % ------------------- Hybrid -------------------
            
            %Determine the Model Contributions
            H{subj,cond,block}(:,1) = [1-c]*w_chosen{subj,cond,block}(:,1) + [c]*q_chosen{subj,cond,block}(:,1);
            H{subj,cond,block}(:,2) = [1-c]*w_chosen{subj,cond,block}(:,2) + [c]*q_chosen{subj,cond,block}(:,2);
            
            %What is the Probability of that action?
            
            P{subj,cond,block} = exp(H{subj,cond,block}(:,1)/beta)./ ...
                [exp(H{subj,cond,block}(:,1)/beta) + exp(H{subj,cond,block}(:,2)/beta)];
            
            % sometimes the probabilty is 1, but the
            % formula outputs a nan and so we have to
            % compensate for that.
            if any(isnan(P{subj,cond,block}))
                P{subj,cond,block}(isnan(P{subj,cond,block})) = H{subj,cond,block}(isnan(P{subj,cond,block}),1) > H{subj,cond,block}(isnan(P{subj,cond,block}),2);
            end
            
            
            Hout{subj,cond,block} = H{subj,cond,block}(:,1);
            
        end
    end
end

%fmincon finds the the local minima, so we will sum up all the
%probabilities and get the negative value of that. In theory, the
%higher the P, the higher the model were able to predict that the
%probability of the choice.
LLE = -sum(log([cat(1,P{:})]));


if ~isfinite(LLE)
    LLE = -sum(numel(a{:})*log(.01));
end

out.q = [q{:}];
out.v = [v{:}];
out.w = [w{:}];

out.q_next = [q_next{:}];
out.v_next = [v_next{:}];
out.w_next = [w_next{:}];

out.P = cat(1,P{:});
out.QPe = [PE_QL{:}];
out.ACPe = [PE_AC{:}];
out.H = cat(1,Hout{:});
end

end