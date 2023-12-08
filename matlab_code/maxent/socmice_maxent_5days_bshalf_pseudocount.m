function socmice_maxent_5days_bshalf_pseudocount
% Need to combine 5 days and learn from 5-day segments
% Learn maxent model from bootstrapped "random" halves of the data

%%
addpath(genpath('../../public_code'))
%%

% nstrain = 'male_before_bsa_200713'; nNodes = 10; % 
% Notice: for male before BSA, there are 12 mice, but 
% after BSA there are only 10 mice. Mouse 3 and mouse 11
% disappeared. So we will only build models for the 10 mice.
% 12;
nstrain = 'male_after_bsa_200727'; nNodes = 10; 
% nstrain = 'female_after_timp_180430'; nNodes = 14;
% nstrain = 'female_before_timp_180413'; nNodes = 14;
% nstrain = 'male_before_timp_180511'; nNodes = 15;
% nstrain = 'male_after_timp_180526'; nNodes = 15;

% nstrain = 'male_c57_181012'; nNodes = 13;
% nstrain = 'male_c57_200701'; nNodes = 10;

cd(['../output/' nstrain])
load bootstrap_index bsIdx2
%
plotFlag = 0;

nDays = 10; 
nStates = 4;
nbs = 10; % alternatively, set nbs = 1 if we want no bootstrapping
lambda = 8;

%% 
for ibs = 1:nbs 

    %%
    s2_bs = [];

%     for dayId  = 0:4  % uncomment for first 5 days of experiment
%         days = 'first5days'; % uncomment for first 5 days of experiment
    for dayId = 5:9 % uncomment for last 5 days of experiment
        days = 'last5days'; % uncomment for last 5 days of experiment


        fn = ['etho_' nstrain '_12hr_dt2s_d' num2str(dayId) '.mat'];
    
        load(fn,'s');
        s = s';
    
%         if strcmp(nstrain, 'male_before_bsa_200713')
%             s = s([1:2 4:10 12],:); % For male before BSA only
%         end

        if strcmp(nstrain, 'male_before_bsa_200713')
            s = s([1:2 4:8 10 12],:); % For male before BSA only. 
%             now also take out the inactive mice number 9 (save as _nim_)
        end

  
        if strcmp(nstrain, 'male_after_bsa_200727')
            s = s([1:7 9:10],:); % nim (= no immobilized mouse)
        end 



        if strcmp(nstrain, 'female_before_timp_180413')
            s = s([1:12 14],:); % For female, exclude mouse # 13
        end
    
        if strcmp(nstrain, 'female_after_timp_180430')
            s = s([1:12 14],:); % For female, exclude mouse # 13
        end
    
        s = int64(s);
        s = remove_tunnel(s);
        s = s(:,1801:1800*7);
        nNodes = size(s,1);

        if nbs > 1
            s_bs = s(:, bsIdx2(ibs + (rem(dayId,5) * 10) ,:));
            s2_bs = [s2_bs s_bs];
        else
            s2_bs = [s2_bs s];
        end

    end

    %% now do the learning
    nSamples = size(s2_bs,2);
    [edgeStruct, nodeMap, edgeMap] = model1_escp_notunnel(nNodes, nStates);
    suffStat = UGM_MRF_computeSuffStat(s2_bs',nodeMap,edgeMap,edgeStruct);
    
    cijExp_eachday = (suffStat((nStates-1)*nNodes+1:end) + lambda/nStates)...
        /(lambda+nSamples) ;
    mirExp_eachday = (suffStat(1:(nStates-1)*nNodes) + lambda/nStates)...
        /(lambda+nSamples) ;
    
    ccijExp_eachday = cm_to_cc2(cijExp_eachday, mirExp_eachday, nNodes, nStates);

    %
    [jij_gd_eachday, hir_gd_eachday, cijExp_eachday, cijExp_bsstd, ...
        mirExp_eachday, mirExp_bsstd, cijMod_eachday, mirMod_eachday,...
        meanj2_gd, stdj2_gd] = ...
        learn_maxent_gd_pseudocount(s2_bs, lambda, nStates, nodeMap, edgeMap, ...
        edgeStruct, plotFlag);




%%
    if nbs > 1

        if strcmp(nstrain, 'female_before_timp_180413') || ...
                strcmp(nstrain, 'female_after_timp_180430') || ...
                strcmp(nstrain, 'male_before_bsa_200713') || ...
                strcmp(nstrain, 'male_after_bsa_200727') 

            fp = ['sme_6h_nbs10_pseudocount_' nstrain ...
                '_nim_' days '_ibs' num2str(ibs) '.mat'];
        else
            fp = ['sme_6h_nbs10_pseudocount_' nstrain ...
                '_' days '_ibs' num2str(ibs) '.mat'];
        end

        save(fp, 's2_bs', 'nNodes', 'nSamples', 'nStates', 'lambda', ...
            'cijExp_eachday','mirExp_eachday','ccijExp_eachday',...
            'jij_gd_eachday', 'hir_gd_eachday', ...
            'meanj2_gd', 'stdj2_gd');  
    else
        
        s2 = s2_bs;


        if strcmp(nstrain, 'female_before_timp_180413') || ...
                strcmp(nstrain, 'female_after_timp_180430') || ...
                strcmp(nstrain, 'male_before_bsa_200713') || ...
                strcmp(nstrain, 'male_after_bsa_200727') 

            fp = ['sme_6h_nbs10_pseudocount_' nstrain ...
                '_nim_' days '_full6hr.mat'];
        else
            fp = ['sme_6h_nbs10_pseudocount_' nstrain ...
                '_' days '_full6hr.mat'];
        end

        save(fp, 's2', 'nNodes', 'nSamples', 'nStates', 'lambda', ...
            'cijExp_eachday','mirExp_eachday','ccijExp_eachday',...
            'jij_gd_eachday', 'hir_gd_eachday', ...
            'meanj2_gd', 'stdj2_gd');

    end

end


end









