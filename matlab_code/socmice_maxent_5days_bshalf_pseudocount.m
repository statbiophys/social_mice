function socmice_maxent_5days_bshalf_pseudocount
% XC 2023-1-18
% XC 2023-1-27
% XC 2023-3-22 changed from socmice_maxent_bshalf_pseudocount
% Need to combine 5 days and learn from 5-day segments
% Learn maxent model from bootstrapped "random" halves of the data

%%
addpath(genpath('~/MyDocuments/biophysics/projects/maxent/UGM'))
addpath(genpath('~/MyDocuments/biophysics/projects/maxent/miscTools'))
addpath(genpath('~/MyDocuments/biophysics/projects/maxent_fdf/matlab_code'))
addpath(genpath('~/MyDocuments/biophysics/projects/maxent/energy_landscape'))
addpath(genpath('~/MyDocuments/biophysics/projects/social_mice'))

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
cd(['~/MyDocuments/biophysics/projects/social_mice/output/' nstrain])
load bootstrap_index bsIdx2
%
plotFlag = 0;

nDays = 10; 
nStates = 4;
% [edgeStruct, nodeMap, edgeMap] = model1_escp_notunnel(nNodes, nStates);


%
% s2 = [];
% s3 = zeros(nNodes, 21600/2, nDays);
nbs = 10;
% cijExp_flat2_bs = zeros(nNodes*(nNodes-1)/2, nDays, nbs);
% ccijExp_flat2_bs = zeros(nNodes*(nNodes-1)/2, nDays, nbs);
% mirExp_flat2_bs = zeros(nNodes*(nStates-1), nDays, nbs);

% cijMod_flat2 = cijExp_flat2;
% mirMod_flat2 = mirExp_flat2;
% ccijMod_flat2 = ccijExp_flat2;


% jij2_gd_bs = cijExp_flat2_bs;
% jij2_l1_bs = cijExp_flat2_bs;
% hir2_gd_bs = mirExp_flat2_bs;
% hir2_l1_bs = mirExp_flat2_bs;

lambda = 8;

%% 
for ibs = 1 %:nbs 

    %%
    s2_bs = [];

%     for dayId  = 0:4 %0:nDays-1
    for dayId = 5:9 %0:nDays-1

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

%         if strcmp(nstrain, 'male_before_bsa_200713')
%             s = s([1 4:8 10 12],:); % For male before BSA only. 
%             % now also take out the inactive mice number 2 and 9 (save as _nim2_)
%         end
%     
        if strcmp(nstrain, 'male_after_bsa_200727')
            s = s([1:7 9:10],:); % nim
%             s = s([1 3:7 9:10],:);  % nim2
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
%         s_bs = s(:, bsIdx2(ibs + (rem(dayId,5) * 10) ,:));
%         s2_bs = [s2_bs s_bs];

        s2_bs = [s2_bs s];
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

    figure(1)
    plot(meanj2_gd)
    hold on
    plot(stdj2_gd)
    hold off

    %
%     [jij_l1_eachday, hir_l1_eachday] = ...
%         learn_maxent_pseudocount(s2_bs, lambda, nStates, nodeMap, edgeMap, ...
%         edgeStruct, plotFlag);

    % save to file (for each random half)
%     fp = ['sme_6h_nbs10_pseudocount_' nstrain ...
%         '_nim_first5days_full6hr.mat'];
     fp = ['sme_6h_nbs10_pseudocount_' nstrain ...
            '_nim_last5days_full6hr.mat'];

%     fp = ['sme_6h_nbs10_pseudocount_' nstrain ...
%         '_first5days_full6hr.mat'];
% %      fp = ['sme_6h_nbs10_pseudocount_' nstrain ...
% %             '_last5days_full6hr.mat'];

%     fp = ['sme_6h_nbs10_pseudocount_' nstrain ...
%         '_first5days_ibs' num2str(ibs) '.mat'];
%     fp = ['sme_6h_nbs10_pseudocount_' nstrain ...
%         '_last5days_ibs' num2str(ibs) '.mat'];

%     fp = ['sme_6h_nbs10_pseudocount_' nstrain ...
%         '_nim_first5days_ibs' num2str(ibs) '.mat'];

%     fp = ['sme_6h_nbs10_pseudocount_' nstrain ...
%         '_nim2_first5days_ibs' num2str(ibs) '.mat'];

%     fp = ['sme_6h_nbs10_pseudocount_' nstrain ...
%         '_nim2_last5days_ibs' num2str(ibs) '.mat'];


%%
    s2 = s2_bs;
    save(fp, 's2', 'nNodes', 'nSamples', 'nStates', 'lambda', ...
        'cijExp_eachday','mirExp_eachday','ccijExp_eachday',...
        'jij_gd_eachday', 'hir_gd_eachday', ...
        'meanj2_gd', 'stdj2_gd');


%     save(fp, 's2_bs', 'nNodes', 'nSamples', 'nStates', 'lambda', ...
%         'cijExp_eachday','mirExp_eachday','ccijExp_eachday',...
%         'jij_gd_eachday', 'hir_gd_eachday', ...
%         'meanj2_gd', 'stdj2_gd');
%         'jij_l1_eachday', 'hir_l1_eachday');
end


end









