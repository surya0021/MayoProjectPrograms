function displayFigure6(analysisType)
%analysis Type : 1= LDA ; 2= Regularized LDA

if ~exist('analysisType','var');      analysisType = 2;   end
    
folderSourceString = 'E:\Mayo\Data\populationSizeLDA';
commonStr = '_TargetOnset-0.5_0_tapers23_alpha8_12_gamma42_78';
analysisTypeStr{1}='LDA';   analysisTypeStr{2}='RLDA';
if analysisType==1
    dataTypeStr{1} = 'FR';      dataTypeStr{2} = 'Gamma';       dataTypeStr{3} = 'Alpha';
else
    dataTypeStr{1} = 'FR';      dataTypeStr{2} = 'Gamma';       dataTypeStr{3} = 'Alpha';       dataTypeStr{4}='FRGammaAlpha';
end
colorName =gray(5);
textPos = zeros(1,46);
for i=1:length(dataTypeStr)
    fileNameStr = (fullfile(folderSourceString,[dataTypeStr{i} analysisTypeStr{analysisType} '_Folds5' commonStr '.mat']));
    load(fileNameStr)
    projection = mProj(1,:); %#ok<NODEF>
    nSession = cell2mat(cellfun(@length,projection,'un',0));
    meanProj = cell2mat(cellfun(@mean,projection,'un',0));
    semProj = cell2mat(cellfun(@(x) std(x)/sqrt(length(x)),projection,'un',0));
    errbar(1:46,meanProj,semProj,colorName(i,:),1);     hold 'on';
    p(i) = plot(1:46,meanProj,'Color',colorName(i,:),'Marker','o','LineWidth',1); %#ok<AGROW>
    textPos = max(textPos,meanProj+semProj);
    clear mProj projection
        
end
for j=1:46
    text(j-0.2,textPos(j)+0.05,num2str(nSession(j)));
end
legend(gca,p,'Firing Rate','Gamma','Alpha','FR+Gamma+Alpha','Location','northwest')

xlabel('Population Size')
ylabel('Z-Score of projection')
end
    
    
function errbar(pos,x,y,color,lineWidth)
for i=1:length(pos)
    line([pos(i) pos(i)],[x(i)-y(i) x(i)+y(i)],'Color',color,'LineWidth',lineWidth)
end
end
