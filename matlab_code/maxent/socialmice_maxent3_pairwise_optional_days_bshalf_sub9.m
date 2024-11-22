function socialmice_maxent3_pairwise_optional_days_bshalf_sub9...
    (istrain, di, df)
% select 9 mice among the ones for male timp and female timp,
% such that it has the same number of mice as the one in bsa
%
% we want to generate results with errorbar 
% by bootstrapping random halves of the data
%%
addpath(genpath('../public_code'))

%%
if istrain == 1
    nstrain = 'male_before_timp_180511';
elseif istrain == 2
    nstrain = 'male_after_timp_180526';
elseif istrain == 3
    nstrain = 'female_before_timp_180413';
elseif istrain == 4
    nstrain = 'female_after_timp_180430';
else 
    error('istrain need to be 1,2,3,4');
end
% nstrain = 'male_before_timp_180511'; nNodes = 15;
% nstrain = 'male_after_timp_180526'; nNodes = 15;
% nstrain = 'female_before_timp_180413'; nNodes = 13; %14;
% nstrain = 'female_after_timp_180430'; % nNodes = 13; %14;

nstrain_save = nstrain;
nbs = 10;

load bootstrap_index bsIdx2
load sub9_index f1_sub9_idx m1_sub9_idx



%%
plotFlag = 0;

% nDays = 10; 
nStates = 4;

lambda = 8;

nNodes = 9;
[edgeStruct, nodeMap, edgeMap] = ...
    model1_escp_notunnel(nNodes, nStates);
edgeStruct.beta_G = 0;




nsubsample = 10;

jijFinal3_nbs_sub9 = zeros(nNodes*(nNodes-1)/2, nbs, nsubsample);
hirFinal3_nbs_sub9 = zeros(nNodes*(nStates-1), nbs, nsubsample);

for isubsample = 1:nsubsample
    for ibs = 1:nbs
        disp(['istrain = ' num2str(istrain) ...
              ',di = ' num2str(di) ...
              ',df = ' num2str(df) ...
              ',isubsample = ' num2str(isubsample) ...
              ',ibs = ' num2str(ibs)]);
        s2_bs = [];
        
        for dayId = di:df
                
            fn = ['etho_' nstrain '_12hr_dt2s_d' num2str(dayId) '.mat'];
        
            
            
            load(fn,'s');
            s = s';
            s = int64(s);
            s = remove_tunnel(s);
            
            s = s(:,1801:12600);
        
        
            if strcmp(nstrain, 'female_before_timp_180413')
                s = s([1:12 14],:); % no dead or immobilized mouse
                s = s(f1_sub9_idx(isubsample, :),:);
                nstrain_save = [nstrain '_sub9mice'];
            end
        
            if strcmp(nstrain, 'female_after_timp_180430')
                s = s([1:12 14],:); % no dead or immobilized mouse
                s = s(f1_sub9_idx(isubsample, :),:);
                nstrain_save = [nstrain '_sub9mice'];
            end
    
            if strcmp(nstrain, 'male_before_timp_180511')
                s = s(m1_sub9_idx(isubsample, :),:);
                nstrain_save = [nstrain '_sub9mice'];
            end
        
            if strcmp(nstrain, 'male_after_timp_180526')
                s = s(m1_sub9_idx(isubsample, :),:);
                nstrain_save = [nstrain '_sub9mice'];
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
    
        jijFinal3_nbs_sub9(:, ibs, isubsample) = jijFinal;
        hirFinal3_nbs_sub9(:, ibs, isubsample) = hirFinal;
    end
end


fp = ['sme_pairwise_l2_' nstrain_save '_di' num2str(di) ...
    '_df' num2str(df) '_bshalf_sub9.mat'];

save(fp,'jijFinal3_nbs_sub9','hirFinal3_nbs_sub9');


end


