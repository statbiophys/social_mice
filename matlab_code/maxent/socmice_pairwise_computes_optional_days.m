function socmice_pairwise_computes_optional_days(di, df, ttscheme_number)
% Compute partition function Z, entropy S, and log(EZ) for a given pairwise
% maxent model
%%
addpath(genpath('../public_code/'))

  %%
nNodes = 15; %
nStates = 4;

[edgeStruct, nodeMap, edgeMap] = ...
    model1_escp_notunnel(nNodes, nStates);

log10_beta_G2 = -1:0.5:1;
beta_G2 = 10.^log10_beta_G2;

% use this if want to compute for no regularization
% beta_G2 = [0];

%%
nstrain = 'male_before_timp_180511'; 

if ttscheme_number == 1
    ttscheme = '1h';
elseif ttscheme_number == 2
    ttscheme = '1day';
end

% fp = ['sme_pairwise_l2_' nstrain '_di' num2str(di) ...
%     '_df' num2str(df) '_' ttscheme '.mat'];

fp = ['sme_pairwise_l2_' nstrain '_di' num2str(di) ...
    '_df' num2str(df) '_' ttscheme '_beta0.mat'];


load(fp,'jijFinal3','hirFinal3');

%%
 
switch ttscheme
    case '1h'
        ntest_max = 6;
    case '1day'
        ntest_max = df-di+1;
    otherwise
        error('ttscheme has to be 1h or 1day')
end

%%
logZ3 = zeros(ntest_max, length(beta_G2));
logEZ3 = zeros(ntest_max, length(beta_G2));
S3 = zeros(ntest_max, length(beta_G2));

%%
for itest = 1:ntest_max
    for ibeta = 1:length(beta_G2)
        disp(['itest = ' num2str(itest) ', ibeta = ' num2str(ibeta)]);

        %%
        jijFinal = jijFinal3(:, itest, ibeta);
        hirFinal = hirFinal3(:, itest, ibeta);

        %%
        [nodePot,edgePot] = UGM_MRF_makePotentials...
            ([hirFinal; jijFinal], ...
            nodeMap,edgeMap,edgeStruct);

        %%
        edgeStruct.useMex = 1;
        tic
        [~, ~, logZ, logEZ, S] = ...
            xcUGM_ComputeS_Exact ...
            (nodePot, edgePot,edgeStruct);
        disp(['N = 15. Pairwise. logZ = ' num2str(logZ) ...
            ', logEZ = ' num2str(logEZ) ...
            ', S = ' num2str(S) ])
        toc

        %%
        logZ3(itest, ibeta) = logZ;
        logEZ3(itest, ibeta) = logEZ;
        S3(itest, ibeta) = S;
    end
end

%%

% fp_new = ['computes_pairwise_l2_' nstrain '_di' num2str(di) ...
%     '_df' num2str(df) '_' ttscheme '.mat'];

fp_new = ['computes_pairwise_l2_' nstrain '_di' num2str(di) ...
    '_df' num2str(df) '_' ttscheme '_beta0.mat'];

save(fp_new, 'logZ3','logEZ3','S3');

% save(['computes_pairwise_iday' num2str(iday) '.mat'],...
%     'iday','logZ','logEZ','S');



end







