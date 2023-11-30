function socialmice_maxent3_pseudocount
% XC 2022-1-18
% Learn static maxent model for different strains of mice
% A cleaner version of socialmice_maxent2
% XC 2023-1-26: add pseudocount to computation of experimental values of
% cij and mir
%%
addpath(genpath('~/MyDocuments/biophysics/projects/maxent/UGM'))
addpath(genpath('~/MyDocuments/biophysics/projects/maxent/miscTools'))
addpath(genpath('~/MyDocuments/biophysics/projects/maxent_fdf/matlab_code'))
addpath(genpath('~/MyDocuments/biophysics/projects/maxent/energy_landscape'))
addpath(genpath('~/MyDocuments/biophysics/projects/social_mice'))
%%
% nstrain = 'btbr_190212'; nNodes = 10;
% nstrain = 'female_after_timp_180430'; nNodes = 14;
% nstrain = 'female_before_timp_180413'; nNodes = 14;
nstrain = 'male_before_timp_180511'; nNodes = 15;
% nstrain = 'male_after_timp_180526'; nNodes = 15;

% nstrain = 'male_before_bsa_200713'; nNodes = 12;
% nstrain = 'male_after_bsa_200727'; nNodes = 10; 
% Warning! The number of mice are different


% nstrain = 'male_c57_181012'; nNodes = 13;
% nstrain = 'male_c57_200701'; nNodes = 10;


%%
plotFlag = 0;

nDays = 10; 
nStates = 4;
[edgeStruct, nodeMap, edgeMap] = ...
    model1_escp_notunnel(nNodes, nStates);
% s3 = zeros(nNodes, 21600, nDays);
s3 = zeros(nNodes, 21600/2, nDays);


cijExp_flat2 = zeros(nNodes*(nNodes-1)/2, nDays);
ccijExp_flat2 = zeros(nNodes*(nNodes-1)/2, nDays);
mirExp_flat2 = zeros(nNodes*(nStates-1), nDays);

% cijMod_flat2 = cijExp_flat2;
% mirMod_flat2 = mirExp_flat2;
% ccijMod_flat2 = ccijExp_flat2;


jij2_gd = cijExp_flat2;
jij2_l1 = cijExp_flat2;
hir2_gd = mirExp_flat2;
hir2_l1 = mirExp_flat2;
 

s2 = [];
%%
for dayId = 0:nDays-1 %nDays-1
    %%
    %%
    dayId = dayId
    
    fn = ['etho_' nstrain '_12hr_dt2s_d' num2str(dayId) '.mat'];
%     fn = ['etho_' nstrain '_12hr_dt1s_d' num2str(dayId) '.mat'];


    load(fn,'s');
    s = s';
    s = int64(s);
    s = remove_tunnel(s);
%     s = s(:,3601:3600*7);

    s = s(:,1801:12600);
    s3(:, :, dayId+1) = s;
        
    %%
    nSamples = size(s,2);
    [edgeStruct, nodeMap, edgeMap] = model1_escp_notunnel(nNodes, nStates);
    suffStat = UGM_MRF_computeSuffStat(s',nodeMap,edgeMap,edgeStruct);
%     cijExp_eachday = suffStat((nStates-1)*nNodes+1:end)/nSamples;
%     mirExp_eachday = suffStat(1:(nStates-1)*nNodes)/nSamples;
    lambda = 1;
    cijExp_eachday = (suffStat((nStates-1)*nNodes+1:end) + lambda/nStates)...
        /(lambda+nSamples) ;
    mirExp_eachday = (suffStat(1:(nStates-1)*nNodes) + lambda/nStates)...
    /(lambda+nSamples) ;

    cijExp_flat2(:,dayId+1) = cijExp_eachday;
    mirExp_flat2(:,dayId+1) = mirExp_eachday;
    
    ccijExp_eachday = cm_to_cc2(cijExp_eachday, mirExp_eachday,nNodes,nStates);
    ccijExp_flat2(:,dayId+1) = ccijExp_eachday;

    
    %%
%     [edgeStruct, nodeMap, edgeMap] = ...
%         model0p8_uniquej_boxsym(nNodes, nStates);
%     [jij_gd_eachday, hir_gd_eachday, cijExp_eachday, cijExp_bsstd, ...
%         mirExp_eachday, mirExp_bsstd, cijMod_eachday, mirMod_eachday] = ...
%         learn_maxent_uniquej_boxsym(s, nStates, nodeMap, edgeMap, ...
%         edgeStruct, plotFlag);

%     %%
    [jij_gd_eachday, hir_gd_eachday, cijExp_eachday, cijExp_bsstd, ...
        mirExp_eachday, mirExp_bsstd, cijMod_eachday, mirMod_eachday] = ...
        learn_maxent_gd_pseudocount(s, nStates, nodeMap, edgeMap, ...
        edgeStruct, plotFlag);
    %%
    [jij_l1_eachday, hir_l1_eachday] = ...
        learn_maxent_pseudocount(s, nStates, nodeMap, edgeMap, ...
        edgeStruct, plotFlag);
    
%%
    jij2_gd(:,dayId+1) = jij_gd_eachday;
    hir2_gd(:,dayId+1) = hir_gd_eachday;
    jij2_l1(:,dayId+1) = jij_l1_eachday;
    hir2_l1(:,dayId+1) = hir_l1_eachday;
    
%     s = random_cyclic_shuffle(s);

% %% same j for all pairs, unique h, for each day
%     [edgeStruct, nodeMap, edgeMap] = ...
%         model0p2_equalj_uniqueh_boxsym(nNodes, nStates);
% %         model0p5_equalj_uniqueh(nNodes, nStates);
% %%
%     [jij_mfj_eachday, hir_mfj_eachday] = ...
%         learn_maxent_meanfieldj_boxsym(s, nStates, nodeMap, ...
%         edgeMap, edgeStruct, plotFlag);
% 
%     %%
%     %
% %     [jij_mfj_eachday, hir_mfj_eachday] = ...
% %         learn_maxent_meanfieldj(s, nStates, nodeMap, edgeMap, edgeStruct, plotFlag);
%     jij2_mfj(:,dayId+1) = jij_mfj_eachday;
%     hir2_mfj(:,dayId+1) = hir_mfj_eachday;
    %%
    cijExp_flat2(:,dayId+1) = cijExp_eachday;
    mirExp_flat2(:,dayId+1) = mirExp_eachday;
    ccijExp_eachday = cm_to_cc2(cijExp_eachday, mirExp_eachday,nNodes,nStates);
    ccijExp_flat2(:,dayId+1) = ccijExp_eachday;
% %     
%     cijMod_flat2(:,dayId+1) = cijMod_eachday;
%     mirMod_flat2(:,dayId+1) = mirMod_eachday;
%     ccijMod_eachday = cm_to_cc2(cijMod_eachday, mirMod_eachday,nNodes,nStates);
%     ccijMod_flat2(:,dayId+1) = ccijMod_eachday;
%     
    %%
    s2 = [s2 s];
end

%%

nSamples = size(s2,2);
suffStat = UGM_MRF_computeSuffStat(s2',nodeMap,edgeMap,edgeStruct);
% cijExp = suffStat((nStates-1)*nNodes+1:end)/nSamples;
% mirExp = suffStat(1:(nStates-1)*nNodes)/nSamples;

cijExp = (suffStat((nStates-1)*nNodes+1:end) + lambda/nStates)...
    /(lambda+nSamples) ;
mirExp = (suffStat(1:(nStates-1)*nNodes) + lambda/nStates)...
/(lambda+nSamples) ;

mirExpMat = zeros(nNodes,nStates);
mirExpMat(:,1:nStates-1) = reshape(mirExp,nNodes,nStates-1);
mirExpMat(:,end) = 1-sum(mirExpMat(:,1:nStates-1),2);

% then learn for all s2
[jij_gd_allday, hir_gd_allday, cijExp, cijExp_bsstd, ...
    mirExp, mirExp_bsstd, cijMod, mirMod] = ...
    learn_maxent_gd_pseudocount(s2, nStates, nodeMap, edgeMap, ...
    edgeStruct, plotFlag);

[jij_l1_allday, hir_l1_allday] = ...
    learn_maxent_pseudocount(s2, nStates, nodeMap, edgeMap, ...
    edgeStruct, plotFlag);
%% Alternative model 1: use the same j for all pairs, unique h
% [edgeStruct, nodeMap, edgeMap] = model0p5_equalj_uniqueh(nNodes, nStates);
% [jijFinal_mfj, hirFinal_mfj] = ...
%     learn_maxent_meanfieldj(s2, nStates, nodeMap, edgeMap, edgeStruct, plotFlag);

%% Alternative model 2: same j for all pairs, symmetric h for boxes with food, same preference for all mice
% [edgeStruct, nodeMap, edgeMap] = model0p1_equalj_boxsym(nNodes, nStates);
% [jijFinal_mfj_mfhsym, hirFinal_mfj_mfhsym] = ...
%     learn_maxent_meanfieldj_mfboxsym(s2, nStates, nodeMap, edgeMap, edgeStruct, plotFlag);
% 
% 
% %% Alternative model 3: same j for all pairs, symmetric h for boxes with food, but different h for each mouse
% [edgeStruct, nodeMap, edgeMap] = model0p2_equalj_uniqueh_boxsym(nNodes, nStates);
% [jijFinal_mfj_hsym, hirFinal_mfj_hsym] = ...
%     learn_maxent_meanfieldj_boxsym(s2, nStates, nodeMap, edgeMap, edgeStruct, plotFlag);

%%
% f5 = first 5 days. l5 = last 5 days
% fp = ['sme_sym_6h_' nstrain '.mat'];
fp = ['sme_6h_pseudocount_' nstrain '.mat'];
save(fp);
% save(fp,'cijExp_flat2', 'mirExp_flat2',...
%     'ccijExp_flat2','-append');

% save(fp,'cijExp','mirExp','-append');
% save(fp, 's2', 's3', 'jij_gd_allday', 'hir_gd_allday', ...
%     'jij2_gd', 'hir2_gd', ...
%     'jijFinal_mfj', 'hirFinal_mfj'); %, ...
%     'jijFinal_mfj_mfhsym', 'hirFinal_mfj_mfhsym', ...
%     'jijFinal_mfj_hsym', 'hirFinal_mfj_hsym');

%     'jij_l1_allday', 'hir_l1_allday', ...

%     'jij2_l1', 'hir2_l1', ...

end


