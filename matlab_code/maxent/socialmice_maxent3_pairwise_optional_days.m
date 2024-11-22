function socialmice_maxent3_pairwise_optional_days...
    (di, df, ttscheme_number)
% Learn pairwise maximum entropy models with training / test split, and
% with various strength of L2 regularization on the pairwise interactions,
% for aggregate data from day di to day df.
%%
addpath(genpath('../public_code'))
%%
nstrain = 'male_before_timp_180511'; nNodes = 15;
%%

log10_beta_G2 = -1:0.5:1;
beta_G2 = 10.^log10_beta_G2;

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
            s_train = [];
            s_test = [];
        
            for dayId = di:df
                    
                fn = ['etho_' nstrain '_12hr_dt2s_d' num2str(dayId) '.mat'];
                
                load(fn,'s');
                s = s';
                s = int64(s);
                s = remove_tunnel(s);
                
                s = s(:,1801:12600);
    
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

%%
fp = ['sme_pairwise_l2_' nstrain '_di' num2str(di) ...
    '_df' num2str(df) '_' ttscheme '.mat'];

save(fp,'jijFinal3','hirFinal3');

end


