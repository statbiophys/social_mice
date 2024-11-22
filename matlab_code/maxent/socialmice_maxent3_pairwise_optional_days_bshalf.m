function socialmice_maxent3_pairwise_optional_days_bshalf...
    (di, df, nbs)
% Learn pairwise maximum entropy model on 5 day aggregate of data
% we want to generate results with errorbar 
% by bootstrapping random halves of the data
%%
addpath(genpath('../public_code'))

%%
% nstrain = 'male_after_timp_180526'; nNodes = 15;
nstrain = 'male_before_timp_180511'; nNodes = 15;
% nstrain = 'male_c57_181012'; nNodes = 13;
% nstrain = 'male_c57_200701'; nNodes = 10;
% nstrain = 'male_before_bsa_200713'; nNodes = 9; %12;
% nstrain = 'male_after_bsa_200727'; nNodes = 9; 
% there are 10 in the data, but one very slow
% nstrain = 'female_before_timp_180413'; nNodes = 13; %14;
% nstrain = 'female_after_timp_180430'; nNodes = 13; %14;

nstrain_save = nstrain;

load bootstrap_index bsIdx2
%%

% log10_beta_G2 = -1:0.5:1;
% beta_G2 = 10.^log10_beta_G2;

% beta_G2 = [0];


%%
plotFlag = 0;

% nDays = 10; 
nStates = 4;

lambda = 8;

[edgeStruct, nodeMap, edgeMap] = ...
    model1_escp_notunnel(nNodes, nStates);
edgeStruct.beta_G = 0;


jijFinal3_nbs = zeros(nNodes*(nNodes-1)/2, nbs);
hirFinal3_nbs = zeros(nNodes*(nStates-1), nbs);

for ibs = 1:nbs
    disp(['ibs = ' num2str(ibs)]);
    s2_bs = [];
    
    for dayId = di:df
            
        fn = ['etho_' nstrain '_12hr_dt2s_d' num2str(dayId) '.mat'];
    
        
        
        load(fn,'s');
        s = s';
        s = int64(s);
        s = remove_tunnel(s);
        
        s = s(:,1801:12600);
    
        if strcmp(nstrain, 'male_after_bsa_200727')
            s = s([1:7 9:10],:); % no dead or immobilized mouse
            nstrain_save = [nstrain '_sub9mice'];
        end
    
        if strcmp(nstrain, 'male_before_bsa_200713')
            s = s([1:2 4:8 10 12],:); % no dead or immobilized mouse
            nstrain_save = [nstrain '_sub9mice'];
        end
    
        if strcmp(nstrain, 'female_before_timp_180413')
            s = s([1:12 14],:); % no dead or immobilized mouse
            nstrain_save = [nstrain '_sub13mice'];
        end
    
        if strcmp(nstrain, 'female_after_timp_180430')
            s = s([1:12 14],:); % no dead or immobilized mouse
            nstrain_save = [nstrain '_sub13mice'];
        end
    
        if nbs > 1
            s_bs = s(:, bsIdx2(ibs + (rem(dayId,5) * 10) ,:));
            s2_bs = [s2_bs s_bs];
        else
            s2_bs = [s2_bs s];
        end
    end
    
   
    [jijFinal, hirFinal] = ...
        learn_maxent_pairwise_pseudocount_regularized_l2 ...
        (s2_bs, lambda, nStates, nodeMap, ...
        edgeMap, edgeStruct, plotFlag);

    jijFinal3_nbs(:, ibs) = jijFinal;
    hirFinal3_nbs(:, ibs) = hirFinal;


end


fp = ['sme_pairwise_l2_' nstrain_save '_di' num2str(di) ...
    '_df' num2str(df) '_bshalf.mat'];

save(fp,'jijFinal3_nbs','hirFinal3_nbs');

end


