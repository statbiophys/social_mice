function socmice_maxent_bshalf_pseudocount
% XC 2023-1-18
% XC 2023-1-27
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
% nstrain = 'male_after_bsa_200727'; nNodes = 10; 
% nstrain = 'female_after_timp_180430'; nNodes = 14;
% nstrain = 'female_before_timp_180413'; nNodes = 14;
% nstrain = 'male_before_timp_180511'; nNodes = 15;
% nstrain = 'male_after_timp_180526'; nNodes = 15;

nstrain = 'male_c57_181012'; nNodes = 13;
% nstrain = 'male_c57_200701'; nNodes = 10;
cd(['~/MyDocuments/biophysics/projects/social_mice/output/' nstrain])
load bootstrap_index bsIdx2
%
plotFlag = 0;

nDays = 10; 
nStates = 4;
[edgeStruct, nodeMap, edgeMap] = model1_escp_notunnel(nNodes, nStates);


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

%% First load the data, and check if there are some "non-engaging" mice
% "Non-engaging" mice are the ones that did not go to all 4 boxes
%%
% % % dayId = 3 5 6
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

    imagesc(s)

    pause
end

%%
plotFlag = 0;
lambda = 8;
for dayId = 0:nDays-1 %0:nDays-1

    dayId = dayId
    
    fn = ['etho_' nstrain '_12hr_dt2s_d' num2str(dayId) '.mat'];

    load(fn,'s');
    s = s';
    if strcmp(nstrain, 'male_before_bsa_200713')
%         s = s([1:2 4:10 12],:); % no dead. For male before BSA only
        s = s([1:2 4:8 10 12],:); % nim
%         s = s([1 4:8 10 12],:); % nim2
    end

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

    nNodes = size(s,1);
    s = int64(s);
    s = remove_tunnel(s);
    s = s(:,1801:1800*7);
%     s3(:,:,dayId+1) = s;

%     s2 = [s2 s];

    %% bootstrap nbs = 10 (start small)

    for ibs = 1:10 % and just don't do the bootstrap :nbs %
        disp(['dayId = ' num2str(dayId) ', ibs = ' num2str(ibs)]);
        s_bs = s(:, bsIdx2(ibs,:));
% % %         s_bs(13,end-3:end) = 1:4;
%         s_bs = s;

        %%
        nSamples = size(s_bs,2);
        [edgeStruct, nodeMap, edgeMap] = model1_escp_notunnel(nNodes, nStates);
        suffStat = UGM_MRF_computeSuffStat(s_bs',nodeMap,edgeMap,edgeStruct);
%         cijExp_eachday = suffStat((nStates-1)*nNodes+1:end)/nSamples;
%         mirExp_eachday = suffStat(1:(nStates-1)*nNodes)/nSamples;
        
        cijExp_eachday = (suffStat((nStates-1)*nNodes+1:end) + lambda/nStates)...
            /(lambda+nSamples) ;
        mirExp_eachday = (suffStat(1:(nStates-1)*nNodes) + lambda/nStates)...
            /(lambda+nSamples) ;

%         cijExp_flat2_bs(:,dayId+1,ibs) = cijExp_eachday;
%         mirExp_flat2_bs(:,dayId+1,ibs) = mirExp_eachday;
        
        ccijExp_eachday = cm_to_cc2(cijExp_eachday, mirExp_eachday,nNodes,nStates);
%         ccijExp_flat2_bs(:,dayId+1,ibs) = ccijExp_eachday;

        %%

%         fp_full = ['sme_6h_nbs10_pseudocount_' nstrain ...
%             '_iday' num2str(dayId+1) '_ibs' num2str(ibs) '.mat'];
%         load(fp_full, 'jij_gd_eachday','hir_gd_eachday');
%         jij_temp_mat = recoverUpt(jij_gd_eachday, 10);
%         jijInitial = jij_temp_mat([1:7 9:10],[1:7 9:10]);
%         jijInitial = makeUpt(jijInitial);
%         hir_temp_mat = reshape(hir_gd_eachday, 10, 3);
%         hirInitial = hir_temp_mat([1:7 9:10],:);
%         hirInitial = reshape(hirInitial, 27, 1);


%         [jij_gd_eachday, hir_gd_eachday, cijExp_eachday, cijExp_bsstd, ...
%             mirExp_eachday, mirExp_bsstd, cijMod_eachday, mirMod_eachday, ...
%             meanj2_gd, stdj2_gd] = ...
%             learn_maxent_gd_pseudocount(s_bs, lambda, nStates, nodeMap, edgeMap, ...
%             edgeStruct, plotFlag, jijInitial, hirInitial);

        [jij_gd_eachday, hir_gd_eachday, cijExp_eachday, cijExp_bsstd, ...
            mirExp_eachday, mirExp_bsstd, cijMod_eachday, mirMod_eachday, ...
            meanj2_gd, stdj2_gd] = ...
            learn_maxent_gd_pseudocount(s_bs, lambda, nStates, nodeMap, edgeMap, ...
            edgeStruct, plotFlag);

%         figure(1)
%         plot(meanj2_gd)
%         hold on
%         plot(stdj2_gd)
%         hold off
        %%
%         [jij_l1_eachday, hir_l1_eachday, ~, ~, ...
%             ~, ~, ~, ~, ...
%             meanj2_l1, stdj2_l1] = ...
%             learn_maxent_pseudocount(s_bs, lambda, nStates, nodeMap, edgeMap, ...
%             edgeStruct, plotFlag);

        %% if the first time gradient descent takes too long, now use 
% %         jij_l1 and hir_l1 as initial points
%         [jij_gd_eachday, hir_gd_eachday, cijExp_eachday, cijExp_bsstd, ...
%             mirExp_eachday, mirExp_bsstd, cijMod_eachday, mirMod_eachday] = ...
%             learn_maxent_gd_pseudocount(s_bs, lambda, nStates, nodeMap, edgeMap, ...
%             edgeStruct, plotFlag, jij_l1_eachday, hir_l1_eachday);
        
        %%
%         jij2_gd_bs(:,dayId+1,ibs) = jij_gd_eachday;
%         hir2_gd_bs(:,dayId+1,ibs) = hir_gd_eachday;
%         jij2_l1_bs(:,dayId+1,ibs) = jij_l1_eachday;
%         hir2_l1_bs(:,dayId+1,ibs) = hir_l1_eachday;
    
        %% save to file (for each day, each bootstrapped half)
        % nim = no inactive mouse
%         fp = ['sme_6h_nbs10_pseudocount_' nstrain ...
%             '_iday' num2str(dayId+1) '_full6hr.mat'];
%         fp = ['sme_6h_nbs10_pseudocount_' nstrain ...
%             '_nim_iday' num2str(dayId+1) '_full6hr.mat'];
%         save(fp, 's', 'nNodes', 'nSamples', 'nStates', 'lambda', ...
%             'cijExp_eachday','mirExp_eachday','ccijExp_eachday',...
%             'jij_gd_eachday', 'hir_gd_eachday', ...
%             'meanj2_gd','stdj2_gd'); 


%           fp = ['sme_6h_nbs10_pseudocount_' nstrain ...
%             '_nim_iday' num2str(dayId+1) '_ibs' num2str(ibs) '.mat'];
        fp = ['sme_6h_nbs10_pseudocount_' nstrain ...
            '_iday' num2str(dayId+1) '_ibs' num2str(ibs) '.mat'];
         save(fp, 's_bs', 'nNodes', 'nSamples', 'nStates', 'lambda', ...
            'cijExp_eachday','mirExp_eachday','ccijExp_eachday',...
            'jij_gd_eachday', 'hir_gd_eachday', ...
            'meanj2_gd','stdj2_gd'); 
         
         
         %, 'meanj2_l1','stdj2_l1');
%             'jij_l1_eachday', 'hir_l1_eachday', ...
            
%         save(fp, 's_bs', 'nNodes', 'nSamples', 'nStates', 'lambda', ...
%             'cijExp_eachday','mirExp_eachday','ccijExp_eachday',...
%             'jij_gd_eachday', 'hir_gd_eachday', ...
%             'jij_l1_eachday', 'hir_l1_eachday');
    
    end
end
%
% fp = ['sme_6h_nbs10_' nstrain '.mat'];
% fp = ['sme_6h_nbs10_pseudocount_' nstrain '_n10.mat'];

% fp = ['sme_6h_nbs10_pseudocount_' nstrain '_nmc20.mat'];
% save(fp);


end









