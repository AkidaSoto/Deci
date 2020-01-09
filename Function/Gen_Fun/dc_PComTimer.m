function dc_PComTimer(Deci)

global PCom

Deci.PComP = Exist(Deci.PComP,'Visual', true);
Deci.PComP = Exist(Deci.PComP,'Audio', true);

PComTime = timer;

PComTime.ExecutionMode = 'fixedSpacing';
PComTime.Period = Deci.PComP.Period;

time = clock;
Sounds.Ding = audioread('Ding.mp3');
Sounds.Error = audioread('ErrorSound.wav');
PComTime.TimerFcn = {@babyTime,time,Sounds,Deci.PComP};
PComTime.StartFcn = {@babyTime,time,Sounds,Deci.PComP};

start(PComTime)

    function babyTime(obj,event,time,Sounds,PComP)
        
        persistent initFinished

        switch event.Type
            case 'StartFcn'
                if PComP.Visual
                    disp('Started PComTimer');
                    disp(['Found ' num2str(length(find(ismember({PCom.State},{'running'})))) ' running jobs']);
                end
                
                if any(~arrayfun(@(c) c.Parent.Parent.Connected,PCom))
                    disp('PCom has old Jobs from deleted pools, will clean them up now')
                    PCom = PCom(arrayfun(@(c) c.Parent.Parent.Connected,PCom));
                end
                
                initFinished = ismember({PCom.State},{'finished'});
                
            case 'TimerFcn'
                
                if ~isempty([PCom.Error])
                    
                    errmessage = arrayfun(@(c) c.message, [PCom.Error],'un',0);
                    
                    if any(ismember(errmessage,{'No workers are available for FevalQueue execution.'}))
                       disp('PCom has old Error Jobs, will clean them up now')
                        PCom = PCom(~ismember(errmessage,{'No workers are available for FevalQueue execution.'}));
                        initFinished = ismember({PCom.State},{'finished'});  
                    end
                    
                    if any(~arrayfun(@(c) c.Parent.Parent.Connected,PCom))
                        disp('PCom has old Jobs from deleted pools, will clean them up now')
                        PCom = PCom(arrayfun(@(c) c.Parent.Parent.Connected,PCom));
                    end
                    
                    if ~isempty([PCom.Error])
                        if PComP.Visual
                            disp('Error!!!');
                            arrayfun(@(c) disp(c.message), [PCom.Error])
                        end
                        
                        if PComP.Audio
                            sound(Sounds.Error,48000);
                            pause(1)
                            sound(Sounds.Error,48000);
                        end
                    end
                    
                    
                elseif ismember({'running'},{PCom.State})
                    time2 = clock;
                    NumFinished = ismember({PCom.State},{'finished'});
                
                    if PComP.Visual
                        disp([num2str(round(etime(time2,time)/60)) ' min has elapsed']);
                        
                        if length(find(NumFinished)) > length(find(initFinished))
                            new = NumFinished & ~initFinished;
                            still = ~NumFinished & ~initFinished;
                            
                            disp(['Finished ' num2str(length(find(new))) ' more job(s)!'])
                            disp(['Still Waiting on ' num2str(length(find(still))) ' job(s).'])
                            
                            initFinished = NumFinished;
                            meantime = arrayfun(@(c) etime(datevec(PCom(c).FinishDateTime),datevec(PCom(c).StartDateTime))/60,find(initFinished),'UniformOutput',false);
                            disp(['Mean Time of Jobs Execution is ' num2str(round(mean([meantime{:}]))) ' mins'])
                            
                        end
                    end   
                else
                    time2 = clock;

                    if PComP.Visual
                        disp([num2str(round(etime(time2,time)/60)) '  min has elapsed']);
                        
                        NumFinished = ismember({PCom.State},{'finished'});
                        initFinished = NumFinished;
                        meantime = arrayfun(@(c) etime(datevec(PCom(c).FinishDateTime),datevec(PCom(c).StartDateTime))/60,find(initFinished),'UniformOutput',false);
                        disp(['Mean Time of Jobs Execution is ' num2str(round(mean([meantime{:}]))) ' mins'])
                        disp([num2str(length(find(initFinished))) ' jobs finished!'])
                    end
                    
                    if PComP.Audio
                        sound(Sounds.Ding,48000); %Test it
                        pause(1)
                        sound(Sounds.Ding,48000); %Test it
                    end
                    
                    delete(timerfindall)
                end
        end
    end

end