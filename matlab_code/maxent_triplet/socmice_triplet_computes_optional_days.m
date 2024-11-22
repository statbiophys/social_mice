function socmice_triplet_computes_optional_days(di, df, ttscheme_number, ibeta)
% Compute partition function Z, entropy S, and log(EZ) for a given triplet
% maxent model
%%
addpath(genpath('../public_code/'))

  %%
nNodes = 15; %
nStates = 4;

[tripletStruct, nodeMap, edgeMap, tripletMap] = ...
    model5_escp_notunnel_triplet(nNodes, nStates);

log10_beta_G2 = -1:0.5:1;
beta_G2 = 10.^log10_beta_G2;


%%
nstrain = 'male_before_timp_180511'; 

if ttscheme_number == 1
    ttscheme = '1h';
elseif ttscheme_number == 2
    ttscheme = '1day';
end

fp = ['sme_triplet_l2_' nstrain '_di' num2str(di) ...
    '_df' num2str(df) '_' ttscheme '.mat'];

load(fp,'jijFinal3','hirFinal3','gijkFinal3');

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
logZ2 = zeros(ntest_max, 1)  ;%length(beta_G2));
logEZ2 = zeros(ntest_max, 1) ; % length(beta_G2));
S2 = zeros(ntest_max, 1); %length(beta_G2));

%%
for itest = 1:ntest_max
%     for ibeta = 1:length(beta_G2)
        disp(['itest = ' num2str(itest) ', ibeta = ' num2str(ibeta)]);

        %%
        jijFinal = jijFinal3(:, itest, ibeta);
        hirFinal = hirFinal3(:, itest, ibeta);
        gijkFinal = gijkFinal3(:, itest, ibeta);

        %%
        [nodePot,edgePot,tripletPot] = xcUGM_MRF_makePotentials_triplet...
            ([hirFinal; jijFinal; gijkFinal], ...
            nodeMap,edgeMap,tripletMap,tripletStruct);

        %%
        edgeStruct.useMex = 1;
        tic
        [~, ~, logZ, logEZ, S] = ...
            xcUGM_ComputeS_Exact_Triplet ...
            (nodePot, edgePot, tripletPot, tripletStruct);

        disp(['N = 15. Pairwise. logZ = ' num2str(logZ) ...
            ', logEZ = ' num2str(logEZ) ...
            ', S = ' num2str(S) ])
        toc

%         fp_ind = ['computes_triplet_l2_' nstrain '_di' num2str(di) ...
%             '_df' num2str(df) '_' ttscheme ...
%             '_itest' num2str(itest) '_ibeta' num2str(ibeta) ...
%             '.mat'];
%         save(fp_ind, 'logZ','logEZ','S');

        %%
        logZ2(itest) = logZ;
        logEZ2(itest) = logEZ;
        S2(itest) = S;
%     end
end

%%
fp_new = ['computes_triplet_l2_' nstrain '_di' num2str(di) ...
    '_df' num2str(df) '_' ttscheme ...
    '_ibeta' num2str(ibeta) '.mat'];
save(fp_new, 'logZ2','logEZ2','S2');

% save(['computes_pairwise_iday' num2str(iday) '.mat'],...
%     'iday','logZ','logEZ','S');



end







