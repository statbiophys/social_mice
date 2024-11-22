function [jijFinal, hirFinal, gijkFinal, cijExp, cijExp_bsstd, ...
    mirExp, mirExp_bsstd, cijkExp, cijkExp_bsstd, ...
    cijMod, mirMod, cijkMod, ...
    meanj2, stdj2] = ...
    learn_maxent_triplet_pseudocount_regularized_l2 ...
    (dsSignal, lambda, nStates, nodeMap, ...
    edgeMap, tripletMap, tripletStruct, plotFlag)
% Implement L2 regularization on the triplet interactions for the triplet
% maximum entropy model

%%
maxIter = 2e4; 
burnIn = 2e3;
nHist = 1; 

tripletStruct.maxIter = maxIter;
[nNodes, nSamples] = size(dsSignal);
tripletStruct.useMex = 1;
%%
suffStat = xcUGM_MRF_computeSuffStat_triplet...
    (dsSignal',nodeMap,edgeMap,tripletMap,tripletStruct);


nEdges = tripletStruct.nEdges;
cijExp = (suffStat((nStates-1)*nNodes+(1:nEdges)) + lambda/nStates)...
    /(lambda+nSamples);
mirExp = (suffStat(1:(nStates-1)*nNodes) + lambda/nStates)...
/(lambda+nSamples);
cijkExp = (suffStat((nStates-1)*nNodes+nEdges+1:end) + lambda/nStates)...
    /(lambda+nSamples);


mirExpMat = zeros(nNodes,nStates);
mirExpMat(:,1:nStates-1) = reshape(mirExp,nNodes,nStates-1);
mirExpMat(:,end) = 1-sum(mirExpMat(:,1:nStates-1),2);



%% run bs half to extract cthres and mthres
nbs = 50;

mirExp_bsmaster = nan(nbs, nNodes*(nStates-1)) ;
cijExp_bsmaster = nan(nbs, nEdges);
nTriplets = tripletStruct.nTriplets;
cijkExp_bsmaster = nan(nbs, nTriplets) ; 


for ibs = 1:nbs
    bsSignal = bootstrapData(dsSignal',200)';
    nBsSamples = size(bsSignal,2);
    suffStat = xcUGM_MRF_computeSuffStat_triplet...
        (bsSignal',nodeMap,edgeMap,tripletMap,tripletStruct);

    mirExp_bsmaster(ibs,:) = ...
        (suffStat(1:(nStates-1)*nNodes) + lambda/nStates)/(lambda+nSamples);
    cijExp_bsmaster(ibs,:) = ...
        (suffStat((nStates-1)*nNodes+(1:nEdges)) + lambda/nStates)...
        /(lambda+nSamples);
    cijkExp_bsmaster(ibs,:) = ...
        (suffStat((nStates-1)*nNodes+nEdges+1:end) + lambda/nStates)...
        /(lambda+nSamples);

end


mirExp_bsstd = nanstd(mirExp_bsmaster)/sqrt(nSamples/nBsSamples);
cijExp_bsstd = nanstd(cijExp_bsmaster)/sqrt(nSamples/nBsSamples);
cijkExp_bsstd = nanstd(cijkExp_bsmaster)/sqrt(nSamples/nBsSamples);


mthres = sqrt(mean(mirExp_bsstd.^2));
cthres = sqrt(mean(cijExp_bsstd.^2));
ctriplet_thres = sqrt(mean(cijkExp_bsstd.^2));

%%

hirInitial = log(mirExpMat./repmat(mirExpMat(:,nStates),1,nStates));
hirInitial = reshape(hirInitial(:,1:nStates-1),nNodes*(nStates-1),1);

jijInitial = zeros(size(cijExp));

gijkInitial = zeros(size(cijkExp));

[nodePot,edgePot,tripletPot] = ...
    xcUGM_MRF_makePotentials_triplet...
    ([hirInitial;jijInitial;gijkInitial],...
    nodeMap,edgeMap,tripletMap,tripletStruct);

subGibbs = xcUGM_Sample_Gibbs_triplet...
    (nodePot,edgePot,tripletPot,tripletStruct,burnIn);

%%
suffStat = xcUGM_MRF_computeSuffStat_triplet...
    (subGibbs',nodeMap,edgeMap,tripletMap,tripletStruct);

mirMod = suffStat(1:(nStates-1)*nNodes)/maxIter;
cijMod = suffStat((nStates-1)*nNodes+(1:nEdges)) / maxIter;
cijkMod = suffStat((nStates-1)*nNodes+nEdges+1:end)/maxIter;



%% update J_ij and h_ir: optimal choice of \Delta J_ij and \Delta h_ir
%%
jijOld = jijInitial;
hirOld = hirInitial;
gijkOld = gijkInitial;

cijModMaster=zeros(0);
mirModMaster=zeros(0);
cijkModMaster=zeros(0);

count=0;


merrormaster=[];
cerrormaster=[];
cijkerrormaster=[];

merrormaster = [merrormaster sqrt(mean(mirExp.^2))];
cerrormaster = [cerrormaster sqrt(mean(cijExp.^2))];
cijkerrormaster = [cijkerrormaster sqrt(mean(cijkExp.^2))];

meanj2 = [];
stdj2 = [];

mcFlag = 1;

beta_G = tripletStruct.beta_G; %10; %1; %cijkExp_bsstd;


%% now, update J and h together 
mcFlag = 1;
count = 0;

%
step_size = 0.25;
step_size_G = 0.1; 

while mcFlag
    %% J and h update using gradient descent

    %%
    jijNew = jijOld - (cijMod - cijExp) * step_size;
    hirNew = hirOld - (mirMod - mirExp) * step_size;
    gijkNew = gijkOld - (cijkMod - cijkExp + beta_G * gijkOld) * step_size_G; 

    
    %%         
  
    [nodePot,edgePot,tripletPot] = xcUGM_MRF_makePotentials_triplet...
        ([hirNew; jijNew; gijkNew], ...
        nodeMap,edgeMap,tripletMap,tripletStruct);
    subGibbs = xcUGM_Sample_Gibbs_triplet...
        (nodePot,edgePot,tripletPot,tripletStruct,burnIn);
    

    suffStat = xcUGM_MRF_computeSuffStat_triplet...
        (subGibbs',nodeMap,edgeMap,tripletMap,tripletStruct);
    
    mirMod = suffStat(1:(nStates-1)*nNodes)/maxIter;
    cijMod = suffStat((nStates-1)*nNodes+(1:nEdges)) / maxIter;
    cijkMod = suffStat((nStates-1)*nNodes+nEdges+1:end)/maxIter;
        
    merror = sqrt(mean((mirMod-mirExp).^2));
    cerror = sqrt(mean((cijMod-cijExp).^2));
    cijkerror = sqrt(mean((cijkMod-cijkExp).^2));

   
 
    %%
    if plotFlag
        figure(1)
        subplot(2,3,1)
        plot(hirOld, hirNew, 'o');
        hold on
        plot([-2 2],[-2 2],'r')
        hold off
        axis square
        axis([-1 1 -1 1]*1.2*(max(abs([hirOld;hirNew]))) )
        xlabel('h, old');
        ylabel('h, new')
        
        subplot(2,3,2)
        plot(jijOld, jijNew, 'o');
        hold on
        plot([-1 1],[-1 1],'r')
        hold off
        axis square
        axis([-1 1 -1 1]*1.2*(max(abs([jijOld;jijNew]))) )
        xlabel('J, old');
        ylabel('J, new')

        subplot(2,3,3)
        plot(gijkOld, gijkNew, 'o');
        hold on
        plot([-1 1],[-1 1],'r')
        hold off
        axis square
        axis([-1 1 -1 1]*1.2*(max(abs([gijkOld;gijkNew]))) )
        xlabel('G, old');
        ylabel('G, new')
        
        subplot(2,3,4)
        plot(mirExp,mirMod,'o')
        hold on
        plot([0 1],[0 1],'r')
        hold off
        axis square
        axis([0.1 0.8 0.1 0.8]);
        xlabel('h_{ir}, data')
        ylabel('h_{ir}, model')
        
        subplot(2,3,5)
        plot(cijExp,cijMod,'o')
        hold on
        plot([0 0.5],[0 0.5],'r')
        hold off
        axis square
        xlabel('C_{ij}, data')
        ylabel('C_{ij}, model')

        subplot(2,3,6)
        plot(cijkExp,cijkMod,'o')
        hold on
        plot([0 0.3],[0 0.3],'r')
        hold off
        axis square
        xlabel('C_{ijk}, data')
        ylabel('C_{ijk}, model')
        
    end

     %%
    count = count + 1;       
    
    hirOld = hirNew;
    jijOld = jijNew;
    gijkOld = gijkNew;
  
    merrormaster = [merrormaster sqrt(mean((mirMod-mirExp).^2))];
    cerrormaster = [cerrormaster sqrt(mean((cijMod-cijExp).^2))];
    cijkerrormaster = [cijkerrormaster sqrt(mean((cijkMod-cijkExp).^2))];

    meanj2 = [meanj2 mean(jijOld)];
    stdj2 = [stdj2 std(jijOld)];

    %% plot
    if plotFlag
        figure(5)
        subplot(1,2,1)         
        plot(merrormaster)
        hold on
        plot(cerrormaster)
        plot(cijkerrormaster)
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
        if cerror<cthres && ...
                merror<mthres && ...
                abs(meanj2(end) - meanj2(end-50)) + ...
            abs(stdj2(end) - stdj2(end-50)) < 0.0025
            mcFlag = 0;
        end
    end
end

disp(['Update h and J and G together took ' num2str(count) ' steps']);

hirFinal = hirOld;
jijFinal = jijOld;
gijkFinal = gijkOld;

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
    subplot(1,3,1)
    errorbar(mirMod, mirExp, mirExp_bsstd,'.')
    hold on
    plot([0 1],[0 1],'color',[0 0 0]+0.65)
    hold off
    axis square
    axis([0 1 0 1])

    subplot(1,3,2)
    errorbar(cijMod, cijExp, cijExp_bsstd,'.')
    hold on
    plot([0 1],[0 1],'color',[0 0 0]+0.65)
    hold off
    axis square
    axis([0 1 0 1])

    subplot(1,3,3)
    errorbar(cijkMod, cijkExp, cijkExp_bsstd,'.')
    hold on
    plot([0 1],[0 1],'color',[0 0 0]+0.65)
    hold off
    axis square
    axis([0 1 0 1])
    
end



end



function [mycijk] = compute_cijk(mys)

%%
nNodes = length(mys);
%%
mycijk = nan(nNodes*(nNodes-1)*(nNodes-2)/6,1);
count = 0;
for iNode = 1:nNodes-2
    for jNode = iNode+1:nNodes-1
        for kNode = jNode+1:nNodes
            count = count + 1;
            mycijk(count) = (mys(iNode) == mys(jNode))*((mys(iNode) == mys(kNode)));
        end
    end
end


end

