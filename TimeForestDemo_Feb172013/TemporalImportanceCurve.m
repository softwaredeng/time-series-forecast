function [] = TemporalImportanceCurve(allForest)

temp1 =what;
dir0=temp1.path;
dir1=[dir0 '\' 'dataset\'];
filename = dir(dir1);fileL=size(filename);fileL=fileL(1);

allData = cell(fileL-2,4);
for i= 1:fileL
    str = filename(i).name;
    sameIf1 = strcmp(str,'.');
    sameIf2 = strcmp(str,'..');
    if(sameIf1+sameIf2>0)%not a valid name
        continue;
    end    
     
   trainDir= [dir1 str '\' str '_TRAIN'];
   TRAIN= load(trainDir); % Onl
   testDir= [dir1 str '\' str '_TEST'];
   TEST= load(testDir); % Onl
   
   cls = TRAIN(:,1); TRAIN(:,1)=[]; cls0 = cls; clsTrain = cls;
   clsTest = TEST(:,1); TEST(:,1)=[];
   nrowTrian = size(TRAIN,1);
   timeL = size(TRAIN,2);
   uniqCls = unique(cls);
   for ii = 1: size(uniqCls,1)
       cls(find(cls0==uniqCls(ii,1))) = ii;
   end        
   uniqCls = unique(cls);
   
   if(size(uniqCls,1)~=2)
       continue;
   end
   
   
    forestCell = allForest(i-2,2); %parameterCell = cell(ntree,1);       
    forestThis = forestCell{:};
    ntree = size(forestThis,1); 

    fea =[];
    for itree=1:ntree
    tree=forestThis{itree};
    tempR=[];
    tempR=TTreeFea(tree,tempR);%[tree.stat tree.entropy tree.split tree.p1 tree.p2 tree.ndata tree.depth tree.minChild];
    fea = [fea;tempR];
    end

    fea = [fea fea(:,5)-fea(:,4)]; %add window size
    
    impV = zeros(3,timeL);
    for iii=1:size(fea,1)
        impV(fea(iii,1)+1,fea(iii,4):fea(iii,5))=impV(fea(iii,1)+1,fea(iii,4):fea(iii,5))+fea(iii,2);
    end
    
    figure,hold on;
    colV = {'red','blue','black'};
    for iii=1:size(impV,1)
    b=plot(impV(iii,:),colV{iii},'LineWidth',2);
    end    
    legend('Mean','Std. Dev.','Slope');     
    xlim([1 size(TRAIN,2)-1]);
    set(gca,'fontsize',14);
    set(get(b,'Parent'),'FontSize',14)    
    title('Temporal Importance Curves')

    c=figure;
    hold on;
    inx = find(cls==min(cls));     
        for jj=1:min(100,size(inx,1))
            temp=inx(jj);
           a=plot(TRAIN(temp,2:end),'-m','MarkerSize',6,'LineWidth',2);
        end        
     inx = find(cls==max(cls));     
        for jj=1:min(100,size(inx,1))
            temp=inx(jj);
           b=plot(TRAIN(temp,2:end),'--g', 'MarkerSize',6,'LineWidth',1);
        end  
    legend([a,b],'class 1','class 2');  
    xlim([1 size(TRAIN,2)-1])
    title(str)
    set(gca,'fontsize',14);
    set(get(b,'Parent'),'FontSize',14)         
end
    