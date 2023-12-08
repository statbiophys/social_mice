function [jijFinal, hirFinal, cijExp, cijExp_bsstd, ...
    mirExp, mirExp_bsstd, cijMod, mirMod, ...
    meanj2, stdj2] = ...
    learn_maxent_gd_pseudocount(dsSignal, lambda, nStates, nodeMap, edgeMap, ...
    edgeStruct, plotFlag, jijInitial, hirInitial)
% Learn the maximum entropy model using gradient descend

maxIter = 2e5;
burnIn = 2e3;
nHist = 5;

edgeStruct.maxIter = maxIter;
[nNodes, nSamples] = size(dsSignal);
% nEdges = edgeStruct.nEdges;
nEdges = max(edgeMap(:))-max(nodeMap(:));

%%
suffStat = UGM_MRF_computeSuffStat(dsSignal',nodeMap,edgeMap,edgeStruct);
cijExp = (suffStat((nStates-1)*nNodes+1:end) + lambda/nStates)/(lambda+nSamples) ;
mirExp = (suffStat(1:(nStates-1)*nNodes) + lambda/nStates)/(lambda+nSamples) ;

mirExpMat = zeros(nNodes,nStates);
mirExpMat(:,1:nStates-1) = reshape(mirExp,nNodes,nStates-1);
mirExpMat(:,end) = 1-sum(mirExpMat(:,1:nStates-1),2);

%% run bs half to extract cthres and mthres
nbs = 50;
cijExp_bsmaster = zeros(0);
mirExp_bsmaster = zeros(0);
ccijExp_bsmaster = zeros(0);

for ibs = 1:nbs
    bsSignal = bootstrapData(dsSignal',200)';
    nBsSamples = size(bsSignal,2);
    suffStat = UGM_MRF_computeSuffStat(bsSignal',nodeMap,edgeMap,edgeStruct);

    mirExp_bsmaster(ibs,:) = ...
        (suffStat(1:(nStates-1)*nNodes) + lambda/nStates ) ...
        / (lambda + nBsSamples) ;
    cijExp_bsmaster(ibs,:) = ...
        (suffStat((nStates-1)*nNodes+1:end) + lambda/nStates ) ...
        / (lambda + nBsSamples) ;
end

cijExp_bsstd = nanstd(cijExp_bsmaster)/sqrt(nSamples/nBsSamples);
mirExp_bsstd = nanstd(mirExp_bsmaster)/sqrt(nSamples/nBsSamples);
cthres = sqrt(mean(cijExp_bsstd.^2));
mthres = sqrt(mean(mirExp_bsstd.^2));


mirExp_bsstd(mirExp_bsstd<0.005) = 0.005; % thresholding data variability
%%

if nargin < 8
    hirInitial = log(mirExpMat./repmat(mirExpMat(:,nStates),1,nStates));
    hirInitial = reshape(hirInitial(:,1:nStates-1),nNodes*(nStates-1),1);
    
    cijExpMat_temp = recoverUpt(cijExp, nNodes) + eye(nNodes);
    jijInitialMat_temp = inv(cijExpMat_temp);
    jijInitial = -makeUpt(jijInitialMat_temp);
end


[nodePot,edgePot] = UGM_MRF_makePotentials([hirInitial;jijInitial],...
    nodeMap,edgeMap,edgeStruct);
subGibbs = UGM_Sample_Gibbs(nodePot, edgePot, edgeStruct, burnIn);
suffStat = UGM_MRF_computeSuffStat(subGibbs',nodeMap,edgeMap,edgeStruct);
cijMod = suffStat((nStates-1)*nNodes+1:end)/maxIter;

mirMod = suffStat(1:(nStates-1)*nNodes)/maxIter;
%%
jijOld = jijInitial;
hirOld = hirInitial;

%%
if plotFlag
    figure(3)
    errorbar(cijMod, cijExp, cijExp_bsstd,'.')
    hold on
    errorbar(mirMod, mirExp, mirExp_bsstd,'.')
    plot([0 1],[0 1],'color',[0 0 0]+0.65)
    hold off
    axis square
    axis([0 1 0 1])
end

%% Gradient descend
step_size = 0.25; 

mcFlag = 1;
count = 0;


cerrormaster = [];
merrormaster = [];

meanj2 = [];
stdj2 = [];

%%
while mcFlag
    if count > 5000 % put in some stop mechanism
        jijOld = nan(size(jijOld));
        hirOld = nan(size(hirOld));
        break
    end
    count = count + 1
    jijNew = jijOld - (cijMod - cijExp) * step_size;
    hirNew = hirOld - (mirMod - mirExp) * step_size;

    [nodePot,edgePot] = UGM_MRF_makePotentials([hirNew;jijNew],...
        nodeMap,edgeMap,edgeStruct);
    subGibbs = UGM_Sample_Gibbs(nodePot, edgePot, edgeStruct, burnIn);
    suffStat = UGM_MRF_computeSuffStat(subGibbs',nodeMap,edgeMap,edgeStruct);
    cijMod = suffStat((nStates-1)*nNodes+1:end)/maxIter;
    mirMod = suffStat(1:(nStates-1)*nNodes)/maxIter;


    jijOld = jijNew;
    hirOld = hirNew; 

    disp(['J = ' num2str(mean(jijOld)) ...
          '; c_exp = ' num2str(mean(cijExp)) ...
          ', c_mod = ' num2str(mean(cijMod))]) 
      
      
    cerror = sqrt(mean((cijMod-cijExp).^2));
    merror = sqrt(mean((mirMod-mirExp).^2));
    
    cerrormaster = [cerrormaster cerror];
    merrormaster = [merrormaster merror];

    meanj2 = [meanj2 mean(jijOld)];
    stdj2 = [stdj2 std(jijOld)];
    


    if count > 50
        if all(abs(cijMod-cijExp)<cijExp_bsstd') && ...
                merror<mthres && ...
                abs(meanj2(end) - meanj2(end-50)) + ...
            abs(stdj2(end) - stdj2(end-50)) < 0.0025
            mcFlag = 0;
        end
    end

   
    if plotFlag
        figure(3)
        subplot(1,3,1)
        errorbar(mirMod, mirExp, mirExp_bsstd,'.')
        hold on
        errorbar(cijMod, cijExp, cijExp_bsstd,'.')
        plot([0 1],[0 1],'color',[0 0 0]+0.65)
        hold off
        axis square
        axis([0 1 0 1])
        
        subplot(1,3,2)
        plot(merrormaster)
        hold on
        plot(cerrormaster)
        hold off

        subplot(1,3,3)
        plot(meanj2)
        hold on
        plot(stdj2)
        hold off

        figure(5)
        imagesc(subGibbs)
    end
    
    pause(0.5)
end
  



jijFinal = jijOld;
hirFinal = hirOld;


end














