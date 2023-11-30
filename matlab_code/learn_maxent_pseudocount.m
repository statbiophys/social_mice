function [jijFinal, hirFinal, cijExp, cijExp_bsstd, ...
    mirExp, mirExp_bsstd, cijMod, mirMod, ...
    meanj2, stdj2] = ...
    learn_maxent_pseudocount(dsSignal, lambda, nStates, nodeMap, ...
    edgeMap, edgeStruct, plotFlag)
% Sections of code from inversePottsUGM_h.m, make this code easier to
% access;
% If we want to change model (the underlying structure of the maxent
% model), we just change edgeMap and nodeMap. That's it.
% xc 2020-1-13

maxIter = 2e5; %2e4;
burnIn = 2e3;
nHist = 5;

edgeStruct.maxIter = maxIter;
[nNodes, nSamples] = size(dsSignal);
% nEdges = edgeStruct.nEdges;
nEdges = max(edgeMap(:))-max(nodeMap(:));


%%
suffStat = UGM_MRF_computeSuffStat(dsSignal',nodeMap,edgeMap,edgeStruct);
% cijExp=suffStat((nStates-1)*nNodes+1:end)/nSamples;
% 
% mirExp=suffStat(1:(nStates-1)*nNodes)/nSamples;
% lambda = 1;
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

%     mirExp_bsmaster(ibs,:)=suffStat(1:(nStates-1)*nNodes)/nBsSamples;

%     if numel(cijExp) > 1  
%         cijExp_bsmaster(ibs,:)=suffStat((nStates-1)*nNodes+1:end)/nBsSamples;
%         ccijExp_bsmaster(ibs,:)= cm_to_cc2(cijExp_bsmaster(ibs,:),...
%             mirExp_bsmaster(ibs,:),15,5);
%     else
%         cijExp_bsmaster(ibs) = ...
%             suffStat((nStates-1)*nNodes+1:end)/nBsSamples...
%             /(nNodes*(nNodes-1)/2);
%     end
end


cijExp_bsstd = nanstd(cijExp_bsmaster)/sqrt(nSamples/nBsSamples);
mirExp_bsstd = nanstd(mirExp_bsmaster)/sqrt(nSamples/nBsSamples);
%ccijExp_bsstd = nanstd(ccijExp_bsmaster)/sqrt(nSamples/nBsSamples);
cthres = sqrt(mean(cijExp_bsstd.^2));
mthres = sqrt(mean(mirExp_bsstd.^2));

%%

hirInitial = log(mirExpMat./repmat(mirExpMat(:,nStates),1,nStates));
hirInitial = reshape(hirInitial(:,1:nStates-1),nNodes*(nStates-1),1);

% jijInitial = zeros(nEdges,1);
jijInitial = zeros(size(cijExp));
% jijInitial = ones(size(cijExp))*0.1;

[nodePot,edgePot] = UGM_MRF_makePotentials([hirInitial;jijInitial],...
    nodeMap,edgeMap,edgeStruct);
subGibbs = UGM_Sample_Gibbs(nodePot, edgePot, edgeStruct, burnIn);
suffStat = UGM_MRF_computeSuffStat(subGibbs',nodeMap,edgeMap,edgeStruct);
cijMod = suffStat((nStates-1)*nNodes+1:end)/maxIter;
% if numel(cijMod) == 1
%     cijMod = cijMod/(nNodes*(nNodes-1)/2);
% end
mirMod = suffStat(1:(nStates-1)*nNodes)/maxIter;

%% update J_ij and h_ir: optimal choice of \Delta J_ij and \Delta h_ir
%%
jijOld = jijInitial;
hirOld = hirInitial;
jijHist0 = jijInitial;
hirHist0 = hirInitial;

cijModMaster=zeros(0);
mirModMaster=zeros(0);

count=0;

cerrormaster=[];
merrormaster=[];
cerrormaster = [cerrormaster sqrt(mean(cijExp.^2))];
merrormaster = [merrormaster sqrt(mean(mirExp.^2))];

meanj2 = [];
stdj2 = [];

mcFlag = 1;

%% New method: first update all h's, then update the J's. This prevent artifitial frustration

mthresFlag = 1;
while mthresFlag
    for iHist = 1:nHist
        %% h update
        disp(['h update. count = ' num2str(count)]);

        deltaLMaster_h = zeros(nNodes*(nStates-1),1);
        deltaHMaster = zeros(nNodes*(nStates-1),1);

        for iSite = 1:nNodes*(nStates-1)
            mm = mirMod(iSite);
            me = mirExp(iSite);
            deltaH = log(me*(1.-mm)/mm/(1.-me));
            deltaHMaster(iSite) = deltaH;
            deltaLMaster_h(iSite) = ...
                -deltaH * me + log(1 + (exp(deltaH) - 1) * mm);
        end

        [deltaLMin, index] = min(deltaLMaster_h);
        hirNew = hirOld;
        hirNew(index) = hirOld(index) + deltaHMaster(index);
        jijNew = jijOld;
%         disp(['hir_old = ' num2str(hirOld(index)) ', hir_new = ' num2str(hirNew(index))]);
%         disp(['mir_model_old = ' num2str(mirMod(index)) ', mir_data_old = ' num2str(mirExp(index))]);

        %% histogram sampling to update c and m, when both h and j changes

        
        
        hist_estimate = xc_UGM_histogram_method(subGibbs',edgeStruct,edgeMap,...
            nodeMap,[hirHist0;jijHist0],[hirNew;jijNew]);
        cijModNew = hist_estimate(nNodes*(nStates-1)+1:end);
%         if numel(cijModNew) == 1
%             cijModNew = cijModNew/(nNodes*(nNodes-1)/2);
%         end
        mirModNew = hist_estimate(1:nNodes*(nStates-1));
        
        count = count + 1;
        cijModMaster(:,count)=cijModNew;
        mirModMaster(:,count)=mirModNew;

%         disp(['cijMod = ' num2str(cijModNew)])
        cijMod=cijModNew;
        mirMod=mirModNew;

        hirOld=hirNew;
        
%         cerrormaster = [cerrormaster sqrt(mean((cijMod-cijExp).^2))];
%         merrormaster = [merrormaster sqrt(mean((mirMod-mirExp).^2))];
    

%         pause;
    end
    
    jijHist0 = jijOld;
    hirHist0 = hirOld;
    
    [nodePot,edgePot] = UGM_MRF_makePotentials([hirHist0;jijHist0],...
        nodeMap,edgeMap,edgeStruct);
    subGibbs = UGM_Sample_Gibbs(nodePot, edgePot, edgeStruct, burnIn);

    suffStat = UGM_MRF_computeSuffStat(subGibbs',nodeMap,edgeMap,edgeStruct);
    mirMod = suffStat(1:(nStates-1)*nNodes)/maxIter;
    merror = sqrt(mean((mirMod-mirExp).^2));
    
    disp(['mthres = ' num2str(mthres) 'merror = ' num2str(merror)]);
    if merror<mthres
        mthresFlag=0;
    end  
end

disp(['Update h took ' num2str(count) ' steps']);

%% now, update J and h together
mcFlag = 1;
count = 0;
while mcFlag
%     count = count + 1
%     if count > 1000
%         mcFlag = 0;
%     end
    for iHist = 1:nHist
        %% j update
        
        disp(['J update. count = ' num2str(count) ', iHist = ' num2str(iHist)]);
        

        deltaLMaster_j = zeros(nEdges,1);
        deltaJMaster = zeros(nEdges,1);

        for iEdge = 1:nEdges
            cm = cijMod(iEdge);
            ce = cijExp(iEdge);
            deltaJ = log(ce*(1.-cm)/cm/(1.-ce));
            deltaJMaster(iEdge) = deltaJ;
            deltaLMaster_j(iEdge) = ...
                - deltaJ * ce + log(1 + (exp(deltaJ) - 1) * cm);
        end

        [deltaLMin, index] = min(deltaLMaster_j);
        jijNew = jijOld;
        jijNew(index) = jijOld(index) + deltaJMaster(index);
%         disp(['jij_old = ' num2str(jijOld(index)) ', jij_new = ' num2str(jijNew(index))]);
%         disp(['cij_model_old = ' num2str(cijMod(index)) ', cij_data_old = ' num2str(cijExp(index))]);
        

        %% h update

        deltaLMaster_h = zeros(nNodes*(nStates-1),1);
        deltaHMaster = zeros(nNodes*(nStates-1),1);

        for iSite = 1:nNodes*(nStates-1)
            mm = mirMod(iSite);
            me = mirExp(iSite);
            deltaH = log(me*(1.-mm)/mm/(1.-me));
            deltaHMaster(iSite) = deltaH;
            deltaLMaster_h(iSite) = ...
                -deltaH * me + log(1 + (exp(deltaH) - 1) * mm);
        end

        [deltaLMin, index] = min(deltaLMaster_h);
        hirNew = hirOld;
        hirNew(index) = hirOld(index) + deltaHMaster(index);
%         disp(['hir_old = ' num2str(hirOld(index)) ', hir_new = ' num2str(hirNew(index))]);
%         disp(['mir_model_old = ' num2str(mirMod(index)) ', mir_data_old = ' num2str(mirExp(index))]);

        %% histogram sampling to update c and m, when both h and j changes

        % test with new method to see if we get the same answer. 2018-10-30

        hist_estimate = xc_UGM_histogram_method(subGibbs',edgeStruct,edgeMap,...
            nodeMap,[hirHist0;jijHist0],[hirNew;jijNew]);
        cijModNew = hist_estimate(nNodes*(nStates-1)+1:end);
%         if numel(cijModNew) == 1
%             cijModNew = cijModNew/(nNodes*(nNodes-1)/2);
%         end
        mirModNew = hist_estimate(1:nNodes*(nStates-1));
        
        count = count + 1;
%         cijModMaster(:,count)=cijModNew;
%         mirModMaster(:,count)=mirModNew;
       

        if plotFlag
            figure(1)
            subplot(2,2,1)
            plot(deltaLMaster_j);
            hold on
            plot(deltaJMaster);
            hold off
            set(gca,'fontsize',14);
            xlabel('pairs');
            axis square
            
            subplot(2,2,2)
            plot(deltaLMaster_h);
            hold on
            plot(deltaHMaster);
            hold off
            set(gca,'fontsize',14);
            xlabel('local field');
            axis square
            
            subplot(2,2,3)
            plot(cijExp,cijMod,'o')
%             hold on
            axis square
%             axis([0 0.25 0 0.25])
%             axis([0.2 0.7 0.2 0.7]);
            set(gca,'fontsize',14);
            
            subplot(2,2,4)
            plot(mirExp,mirMod,'o')
%             hold on
            axis square
            axis([0.1 0.8 0.1 0.8]);
            set(gca,'fontsize',14);
            
        end
        
        cijMod=cijModNew;
        mirMod=mirModNew;

        jijOld=jijNew;
        hirOld=hirNew;
        
        cerrormaster = [cerrormaster sqrt(mean((cijMod-cijExp).^2))];
        merrormaster = [merrormaster sqrt(mean((mirMod-mirExp).^2))];
    
        meanj2 = [meanj2 mean(jijOld)];
        stdj2 = [stdj2 std(jijOld)];
        
        if plotFlag
%             cerrormaster = [cerrormaster sqrt(mean((cijMod-cijExp).^2))];
%             merrormaster = [merrormaster sqrt(mean((mirMod-mirExp).^2))];
%         
%             meanj2 = [meanj2 mean(jijOld)];
%             stdj2 = [stdj2 std(jijOld)];

            figure(5)
            subplot(1,2,1)
            plot(cerrormaster)
            hold on
            plot(merrormaster)
            hold off

            subplot(1,2,2)
            plot(meanj2)
            hold on
            plot(stdj2)
            hold off
            pause(0.1);
            
        end
%         pause;
    end

%     if plotFlag
%         figure(1)
%         subplot(2,2,1)
%         plot(deltaLMaster_j);
%         hold on
%         plot(deltaJMaster);
%         hold off
%         set(gca,'fontsize',14);
%         xlabel('pairs');
%         axis square
%         
%         subplot(2,2,2)
%         plot(deltaLMaster_h);
%         hold on
%         plot(deltaHMaster);
%         hold off
%         set(gca,'fontsize',14);
%         xlabel('local field');
%         axis square
%         
%         subplot(2,2,3)
%         plot(cijExp,cijMod,'o')
% %             hold on
%         axis square
% %             axis([0 0.25 0 0.25])
% %             axis([0.2 0.7 0.2 0.7]);
%         set(gca,'fontsize',14);
%         
%         subplot(2,2,4)
%         plot(mirExp,mirMod,'o')
% %             hold on
%         axis square
%         axis([0.1 0.8 0.1 0.8]);
%         set(gca,'fontsize',14);
%         
%     end
% 
%     if plotFlag
%         cerrormaster = [cerrormaster sqrt(mean((cijMod-cijExp).^2))];
%         merrormaster = [merrormaster sqrt(mean((mirMod-mirExp).^2))];
%     
%         meanj2 = [meanj2 mean(jijOld)];
%         stdj2 = [stdj2 std(jijOld)];
% 
%         figure(5)
%         subplot(1,2,1)
%         plot(cerrormaster)
%         hold on
%         plot(merrormaster)
%         hold off
% 
%         subplot(1,2,2)
%         plot(meanj2)
%         hold on
%         plot(stdj2)
%         hold off
%         pause(0.5);
%             
%     end

    
    jijHist0 = jijOld;
    hirHist0 = hirOld;
    

    [nodePot,edgePot] = UGM_MRF_makePotentials([hirHist0;jijHist0],...
        nodeMap,edgeMap,edgeStruct);
    subGibbs = UGM_Sample_Gibbs(nodePot, edgePot, edgeStruct, burnIn);

    suffStat = UGM_MRF_computeSuffStat(subGibbs',nodeMap,edgeMap,edgeStruct);
    cijMod = suffStat((nStates-1)*nNodes+1:end)/maxIter;
%     if numel(cijMod) == 1
%         cijMod = cijMod/(nNodes*(nNodes-1)/2);
%     end
    mirMod = suffStat(1:(nStates-1)*nNodes)/maxIter;
    
    cerror = sqrt(mean((cijMod-cijExp).^2));
    merror = sqrt(mean((mirMod-mirExp).^2));
    
%     if (cerror<cthres)&&(merror<mthres)
%     if all(abs(cijMod-cijExp)<cthres) && all(abs(mirMod-mirExp)<mthres)
%     if all(abs(cijMod-cijExp)<cijExp_bsstd') && ...
%                 all(abs(mirMod-mirExp)<mirExp_bsstd')
%         mcFlag=0;
%     end
    
    if count > 100
        if all(abs(cijMod-cijExp)<cijExp_bsstd') && ...
                merror<mthres && ...
                abs(meanj2(end) - meanj2(end-100)) + ...
            abs(stdj2(end) - stdj2(end-100)) < 0.01 %0.005
            mcFlag = 0;
        end
    end

    
end

disp(['Update h and J together took ' num2str(count) ' steps']);

jijFinal = jijHist0;
hirFinal = hirHist0;

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
end


end