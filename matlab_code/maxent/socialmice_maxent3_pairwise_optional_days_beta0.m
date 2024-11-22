function socialmice_maxent3_pairwise_optional_days_beta0...
    (di, df, ttscheme_number)
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
%%

% log10_beta_G2 = -1:0.5:1;
% beta_G2 = 10.^log10_beta_G2;

beta_G2 = [0];

%%
if ttscheme_number == 1
    ttscheme = '1h';
elseif ttscheme_number == 2
    ttscheme = '1day';
end

%%
plotFlag = 0;

nDays = 10; 
nStates = 4;

lambda = 8;

[edgeStruct, nodeMap, edgeMap] = ...
    model1_escp_notunnel(nNodes, nStates);


switch ttscheme
    case '1h'

        jijFinal3 = zeros(nNodes*(nNodes-1)/2, 6, length(beta_G2));
        hirFinal3 = zeros(nNodes*(nStates-1), 6, length(beta_G2));

        for itest = 1:6
            %%
            s_train = [];
            s_test = [];
        
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
    
                test_idx = (itest-1) * 1080+(1:1080);
                train_idx = setdiff(1:10800, test_idx);
                s_train = [s_train s(:, train_idx)];
                s_test = [s_test s(:, test_idx)];
            end
            
            for ibeta = 1:length(beta_G2)
                beta_G = beta_G2(ibeta);
                edgeStruct.beta_G = beta_G;

                [jijFinal, hirFinal] = ...
                    learn_maxent_pairwise_pseudocount_regularized_l2 ...
                    (s_train, lambda, nStates, nodeMap, ...
                    edgeMap, edgeStruct, plotFlag);

                jijFinal3(:, itest, ibeta) = jijFinal;
                hirFinal3(:, itest, ibeta) = hirFinal;
            end
        end

    case '1day'

        jijFinal3 = zeros(nNodes*(nNodes-1)/2, df-di+1, length(beta_G2));
        hirFinal3 = zeros(nNodes*(nStates-1), df-di+1, length(beta_G2));

        for itest = 1:(df-di+1)
            s2 = [];
  
            for dayId = di:df     
                fn = ['etho_' nstrain '_12hr_dt2s_d' num2str(dayId) '.mat'];
                load(fn,'s');
                s = s';
                s = int64(s);
                s = remove_tunnel(s);
                
                s = s(:,1801:12600);                
                s2 = [s2 s];
            end
            
            tf = size(s2, 2);
            test_idx = (itest-1) * 10800+(1:10800);
            train_idx = setdiff(1:tf, test_idx);
            s_train = s2(:, train_idx);
            s_test = s2(:, test_idx);

            %%

            for ibeta = 1:length(beta_G2)
                beta_G = beta_G2(ibeta);
                edgeStruct.beta_G = beta_G;

                [jijFinal, hirFinal] = ...
                    learn_maxent_pairwise_pseudocount_regularized_l2 ...
                    (s_train, lambda, nStates, nodeMap, ...
                    edgeMap, edgeStruct, plotFlag);

                jijFinal3(:, itest, ibeta) = jijFinal;
                hirFinal3(:, itest, ibeta) = hirFinal;
            end
        end

    otherwise
        error('ttscheme can only be 1h or 1day.')
end

%% also do one learning on the combined training and test set
edgeStruct.beta_G = 0;
% plotFlag = 1;
[jijFinal, hirFinal] = learn_maxent_pairwise_pseudocount_regularized_l2 ...
                    ([s_train s_test], lambda, nStates, nodeMap, ...
                    edgeMap, edgeStruct, plotFlag);
jijFinal_alldata = jijFinal;
hirFinal_alldata = hirFinal;

%%
% fp = ['sme_pairwise_l2_' nstrain '_di' num2str(di) ...
%     '_df' num2str(df) '_' ttscheme '.mat'];


fp = ['sme_pairwise_l2_' nstrain_save '_di' num2str(di) ...
    '_df' num2str(df) '_' ttscheme '_beta0.mat'];

save(fp,'jijFinal3','hirFinal3', 'jijFinal_alldata','hirFinal_alldata');


end


