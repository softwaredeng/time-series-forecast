% --- the results may be different from each run due to randomness
function [allErr,allForest]=TForest()

temp1 =what;dir0=temp1.path;dir1=[dir0 '\' 'dataset\'];
filename = dir(dir1);fileL=size(filename);fileL=fileL(1);
alphaV = 1; % 0: entropy gain; 1: Entrance gain

nRows = (fileL-2);
allForest=cell(fileL-2,2);
allErr=cell(nRows,4);
allErr(1,1)={'Data'};allErr(1,2)={'Error'};
allErr(1,3)={'Training Time'};allErr(1,4)={'Testing Time'};

%fileL
isOpen = matlabpool('size') > 0;
if isOpen~=1 
    %matlabpool open 4
end

counter=1;
for alpha=alphaV
    
%for i= 3:fileL
for i= 3:fileL
    forest = [];
    R=[];
    str = filename(i).name;
    sameIf1 = strcmp(str,'.');
    sameIf2 = strcmp(str,'..');
    if(sameIf1+sameIf2>0)%not a valid name
        continue;
    end
 
    
   testDir= [dir1 str '\' str '_TEST'];
   trainDir= [dir1 str '\' str '_TRAIN'];
   TEST = load(testDir); % Only these two lines need to be changed to test a different dataset. %
   TRAIN= load(trainDir); % Only these two lines need to be changed to test a different dataset. %
 
       
    X=TRAIN;X(:,1)=[];
    TL = size(TRAIN,2);
    cls=TRAIN(:,1);
    
    minWin = 1;  maxWin = TL;
    
    ntree=500;
    pre_l=1000000; %no limit on depth
    num = 20;nVar=num;nMean=num;nSlope=num;%here nVar is substitute the second derivative, i.e. smoothness
    
    sampleModeWSZ = 1; %window size sample: 0: log2;  1: sqrt; 2: all
    sampleModePos = 1; %position sample:              1: sqrt; 2: all
    
    tic;
    forest = cell(ntree,1);      
    %parfor: multicore computing. replace to "for" for single core
    parfor itree=1:ntree
        inx = randsample(size(X,1),ceil(size(X,1)*2/2),1);%1: with replacement; 0: without replacement  
        depth=0;
        [tree1] = TMakeTree(depth,pre_l,X(inx,:),cls(inx,:),nVar,nMean,nSlope,sampleModeWSZ,sampleModePos,alpha,minWin,maxWin);
        forest{itree}=tree1;
    end
    trainT = toc;
    
    tic;
    TEST1 = TEST;
    clsTest=TEST1(:,1);
    TEST1(:,1)=[];
    
    %---batch processing test cases
    clsBatch=TReadForestBatchTest(forest,ntree,TEST1);
    err = 1-sum((clsBatch-clsTest(:,1))==0)/size(TEST1,1);
    testT = toc;
    
    if(size(alphaV,2)==1)allForest(counter,1)={str};allForest(counter,2)={forest};end
    allErr(counter+1,1)={str};allErr(counter+1,2)={err};
    allErr(counter+1,3)={trainT};allErr(counter+1,4)={testT};
    % err
    allErr;
    counter=counter+1;
end
end


%matlabpool close 

isOpen = matlabpool('size') > 0;
if isOpen==1 
    %matlabpool close 
end
