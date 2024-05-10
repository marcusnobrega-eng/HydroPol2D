function [state,options,optchanged]=ga_save_each_gen(options,state,flag)
        Score_gen=state.Score;
        Population_gen=state.Population;
        Generation_gen=state.Generation;
        optchanged=[];
        save(['gen_' num2str(Generation_gen,'%.4d') '.mat'],'Score_gen','Population_gen','Generation_gen')
end