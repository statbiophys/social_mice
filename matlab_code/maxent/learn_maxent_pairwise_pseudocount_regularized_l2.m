function [jijFinal, hirFinal, cijExp, cijExp_bsstd, ...
    mirExp, mirExp_bsstd, ...
    cijMod, mirMod, ...
    meanj2, stdj2] = ...
    learn_maxent_pairwise_pseudocount_regularized_l2 ...
    (dsSignal, lambda, nStates, nodeMap, ...
    edgeMap, edgeStruct, plotFlag)
% implement L2 regularization on the pairwise interactions for pairwise
% maximum entropy model

%%
maxIter = 2e4; 
burnIn = 2e3;
nHist = 1; 

edgeStruct.maxIter = maxIter;
[nNodes, nSamples] = size(dsSignal);
edgeStruct.useMex = 1;
%%
suffStat = UGM_MRF_computeSuffStat...
    (dsSignal',nodeMap,edgeMap,edgeStruct);


nEdges = edgeStruct.nEdges;
cijExp = (suffStat((nStates-1)*nNodes+(1:nEdges)) + lambda/nStates)...
    /(lambda+nSamples);
mirExp = (suffStat(1:(nStates-1)*nNodes) + lambda/nStates)...
/(lambda+nSamples);



mirExpMat = zeros(nNodes,nStates);
mirExpMat(:,1:nStates-1) = reshape(mirExp,nNodes,nStates-1);
mirExpMat(:,end) = 1-sum(mirExpMat(:,1:nStates-1),2);



%% run bs half to extract cthres and mthres
nbs = 50;

mirExp_bsmaster = nan(nbs, nNodes*(nStates-1)) ; %zeros(0);
cijExp_bsmaster = nan(nbs, nEdges); %zeros(0);



for ibs = 1:nbs
    bsSignal = bootstrapData(dsSignal',200)';
    nBsSamples = size(bsSignal,2);
    suffStat = UGM_MRF_computeSuffStat...
        (bsSignal',nodeMap,edgeMap,edgeStruct);

    mirExp_bsmaster(ibs,:) = ...
        (suffStat(1:(nStates-1)*nNodes) + lambda/nStates)/(lambda+nSamples);
    cijExp_bsmaster(ibs,:) = ...
        (suffStat((nStates-1)*nNodes+(1:nEdges)) + lambda/nStates)...
        /(lambda+nSamples);
end


mirExp_bsstd = nanstd(mirExp_bsmaster)/sqrt(nSamples/nBsSamples);
cijExp_bsstd = nanstd(cijExp_bsmaster)/sqrt(nSamples/nBsSamples);


mthres = sqrt(mean(mirExp_bsstd.^2));
cthres = sqrt(mean(cijExp_bsstd.^2));

%%

hirInitial = log(mirExpMat./repmat(mirExpMat(:,nStates),1,nStates));
hirInitial = reshape(hirInitial(:,1:nStates-1),nNodes*(nStates-1),1);

jijInitial = zeros(nEdges,1);


[nodePot,edgePot] = ...
    UGM_MRF_makePotentials...
    ([hirInitial;jijInitial],...
    nodeMap,edgeMap,edgeStruct);

subGibbs = UGM_Sample_Gibbs...
    (nodePot,edgePot,edgeStruct,burnIn);

%%
suffStat = UGM_MRF_computeSuffStat...
    (subGibbs',nodeMap,edgeMap,edgeStruct);

mirMod = suffStat(1:(nStates-1)*nNodes)/maxIter;
cijMod = suffStat((nStates-1)*nNodes+(1:nEdges)) / maxIter;



%% update J_ij and h_ir: optimal choice of \Delta J_ij and \Delta h_ir
%%
jijOld = jijInitial;
hirOld = hirInitial;

cijModMaster = zeros(0);
mirModMaster = zeros(0);

count=0;


merrormaster=[];
cerrormaster=[];

merrormaster = [merrormaster sqrt(mean(mirExp.^2))];
cerrormaster = [cerrormaster sqrt(mean(cijExp.^2))];

meanj2 = [];
stdj2 = [];

mcFlag = 1;

% beta_h = mirExp_bsstd;
% beta_J = cijExp_bsstd;
beta_G = edgeStruct.beta_G; %10; %1; %cijkExp_bsstd;


%% now, update J and h together 
% for triplet, do we want to use pairwise information? or not? 
% do we want to only use beta for the triplets? 
mcFlag = 1;
count = 0;

%
step_size = 0.1; 

while mcFlag
    %% J and h update using gradient descent

    %%
    
    hirNew = hirOld - (mirMod - mirExp) * step_size;
    jijNew = jijOld - (cijMod - cijExp + beta_G * jijOld) * step_size;

    
    %%         
  
    [nodePot,edgePot] = UGM_MRF_makePotentials...
        ([hirNew; jijNew], ...
        nodeMap,edgeMap,edgeStruct);
    subGibbs = UGM_Sample_Gibbs...
        (nodePot,edgePot,edgeStruct,burnIn);
    

    suffStat = UGM_MRF_computeSuffStat...
        (subGibbs',nodeMap,edgeMap,edgeStruct);
    
    mirMod = suffStat(1:(nStates-1)*nNodes)/maxIter;
    cijMod = suffStat((nStates-1)*nNodes+(1:nEdges)) / maxIter;
        
    merror = sqrt(mean((mirMod-mirExp).^2));
    cerror = sqrt(mean((cijMod-cijExp).^2));

   
 
    %%
    if plotFlag
        figure(1)
        subplot(2,2,1)
        plot(hirOld, hirNew, 'o');
        hold on
        plot([-2 2],[-2 2],'r')
        hold off
        axis square
        axis([-1 1 -1 1]*1.2*(max(abs([hirOld;hirNew]))) )
        xlabel('h, old');
        ylabel('h, new')
        
        subplot(2,2,2)
        plot(jijOld, jijNew, 'o');
        hold on
        plot([-1 1],[-1 1],'r')
        hold off
        axis square
        axis([-1 1 -1 1]*1.2*(max(abs([jijOld;jijNew]))) )
        xlabel('J, old');
        ylabel('J, new')
        
        subplot(2,2,3)
        plot(mirExp,mirMod,'o')
        hold on
        plot([0 1],[0 1],'r')
        hold off
        axis square
        axis([0.1 0.8 0.1 0.8]);
        xlabel('h_{ir}, data')
        ylabel('h_{ir}, model')
        
        subplot(2,2,4)
        plot(cijExp,cijMod,'o')
        hold on
        plot([0 0.5],[0 0.5],'r')
        hold off
        axis square
        xlabel('C_{ij}, data')
        ylabel('C_{ij}, model')

        
        
    end

     %%
    count = count + 1;       
    
    hirOld = hirNew;
    jijOld = jijNew;
  
    merrormaster = [merrormaster sqrt(mean((mirMod-mirExp).^2))];
    cerrormaster = [cerrormaster sqrt(mean((cijMod-cijExp).^2))];

    meanj2 = [meanj2 mean(jijOld)];
    stdj2 = [stdj2 std(jijOld)];

    %% plot
    if plotFlag
        figure(5)
        subplot(1,2,1)         
        plot(merrormaster)
        hold on
        plot(cerrormaster)
        hold off

        subplot(1,2,2)
        plot(meanj2)
        hold on
        plot(stdj2)
        hold off        
    end

    %%

    if count > 1000%500
        mcFlag = 0;
    end

    %%
    if count > 50
        if beta_G ~= 0
            if merror<mthres && ...
                    abs(meanj2(end) - meanj2(end-50)) + ...
                abs(stdj2(end) - stdj2(end-50)) < 0.0025
                mcFlag = 0;
            end
        else
            if cerror<cthres && merror<mthres && ...
                    abs(meanj2(end) - meanj2(end-50)) + ...
                abs(stdj2(end) - stdj2(end-50)) < 0.0025
                mcFlag = 0;
            end
        end
    end
end

disp(['Update h and J  together took ' num2str(count) ' steps']);

hirFinal = hirOld;
jijFinal = jijOld;

%% visualize
if plotFlag
    figure(13)
    errorbar(cijMod, cijExp, cijExp_bsstd,'.')
    hold on
    errorbar(mirMod, mirExp, mirExp_bsstd,'.')
    plot([0 1],[0 1],'color',[0 0 0]+0.65)
    hold off
    axis square
    axis([0 1 0 1])

    %%
    figure(14)
    subplot(1,2,1)
    errorbar(mirMod, mirExp, mirExp_bsstd,'.')
    hold on
    plot([0 1],[0 1],'color',[0 0 0]+0.65)
    hold off
    axis square
    axis([0 1 0 1])

    subplot(1,2,2)
    errorbar(cijMod, cijExp, cijExp_bsstd,'.')
    hold on
    plot([0 1],[0 1],'color',[0 0 0]+0.65)
    hold off
    axis square
    axis([0 1 0 1])

    
end
end

