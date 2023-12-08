function socmice_maxent_bshalf_pseudocount
% Learn maxent model from bootstrapped random halves of the data

%%
addpath(genpath('../../public_code'))

%%

% nstrain = 'male_before_bsa_200713'; nNodes = 10; % 
% Notice: for male before BSA, there are 12 mice, but 
% after BSA there are only 10 mice. Mouse 3 and mouse 11
% disappeared. So we will only build models for the 10 mice.

% nstrain = 'male_after_bsa_200727'; nNodes = 10; 
% nstrain = 'female_after_timp_180430'; nNodes = 14;
% nstrain = 'female_before_timp_180413'; nNodes = 14;
% nstrain = 'male_before_timp_180511'; nNodes = 15;
% nstrain = 'male_after_timp_180526'; nNodes = 15;

nstrain = 'male_c57_181012'; nNodes = 13;
% nstrain = 'male_c57_200701'; nNodes = 10;
cd(['../output/' nstrain])

load bootstrap_index bsIdx2
%
plotFlag = 0;

nDays = 10; 
nStates = 4;
[edgeStruct, nodeMap, edgeMap] = model1_escp_notunnel(nNodes, nStates);
nbs = 10;
lambda = 8;

%% First load the data, and check if there are some "non-engaging" mice
% "Non-engaging" mice are the ones that did not go to all 4 boxes
for dayId = 0:nDays-1

    %%
    fn = ['etho_' nstrain '_12hr_dt2s_d' num2str(dayId) '.mat'];

    load(fn,'s');
    s = s';

%     if strcmp(nstrain, 'male_before_bsa_200713')
%         s = s([1:2 4:8 10 12],:); % For male before BSA only
%     end
% 
%     if strcmp(nstrain, 'male_after_bsa_200727')
%         s = s([1:7 9:10],:); % For male after BSA only
%     end
% 
%     if strcmp(nstrain, 'female_after_timp_180430')
%         s = s([1:12 14],:); % For female, exclude mouse # 13
%     end
% 
%     if strcmp(nstrain, 'female_before_timp_180413')
%         s = s([1:12 14],:); % For female, exclude mouse # 13
%     end

    s = int64(s);
    s = remove_tunnel(s);
    s = s(:,1801:1800*7);

%     imagesc(s)

end

%%

for dayId = 0:nDays-1 

    dayId = dayId
    
    fn = ['etho_' nstrain '_12hr_dt2s_d' num2str(dayId) '.mat'];

    load(fn,'s');
    s = s';
    if strcmp(nstrain, 'male_before_bsa_200713')
        s = s([1:2 4:8 10 12],:); % no dead or immobilized mouse
    end

    if strcmp(nstrain, 'male_after_bsa_200727')
            s = s([1:7 9:10],:); % no dead or immobilized mouse
    end

    if strcmp(nstrain, 'female_before_timp_180413')
        s = s([1:12 14],:); % no dead or immobilized mouse
    end

    if strcmp(nstrain, 'female_after_timp_180430')
        s = s([1:12 14],:); % no dead or immobilized mouse
    end

    nNodes = size(s,1);
    s = int64(s);
    s = remove_tunnel(s);
    s = s(:,1801:1800*7);


    %% bootstrap nbs = 10 (start small)

    for ibs = 1:nbs % and just don't do the bootstrap :nbs %
        disp(['dayId = ' num2str(dayId) ', ibs = ' num2str(ibs)]);
        s_bs = s(:, bsIdx2(ibs,:));


        %%
        nSamples = size(s_bs,2);
        [edgeStruct, nodeMap, edgeMap] = model1_escp_notunnel(nNodes, nStates);
        suffStat = UGM_MRF_computeSuffStat(s_bs',nodeMap,edgeMap,edgeStruct);
        
        cijExp_eachday = (suffStat((nStates-1)*nNodes+1:end) + lambda/nStates)...
            /(lambda+nSamples) ;
        mirExp_eachday = (suffStat(1:(nStates-1)*nNodes) + lambda/nStates)...
            /(lambda+nSamples) ;

        
        ccijExp_eachday = cm_to_cc2(cijExp_eachday, mirExp_eachday,nNodes,nStates);

        %%
        [jij_gd_eachday, hir_gd_eachday, cijExp_eachday, cijExp_bsstd, ...
            mirExp_eachday, mirExp_bsstd, cijMod_eachday, mirMod_eachday, ...
            meanj2_gd, stdj2_gd] = ...
            learn_maxent_gd_pseudocount(s_bs, lambda, nStates, nodeMap, edgeMap, ...
            edgeStruct, plotFlag);

        %% 

        fp = ['sme_6h_nbs10_pseudocount_' nstrain ...
            '_iday' num2str(dayId+1) '_ibs' num2str(ibs) '.mat'];
         save(fp, 's_bs', 'nNodes', 'nSamples', 'nStates', 'lambda', ...
            'cijExp_eachday','mirExp_eachday','ccijExp_eachday',...
            'jij_gd_eachday', 'hir_gd_eachday', ...
            'meanj2_gd','stdj2_gd'); 
           
    end
end


end









