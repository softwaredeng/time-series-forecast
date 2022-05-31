
function [clsBatch,clsV]=TReadForestBatchTest(forest,ntree,X)
%length=length+1;
clsV =zeros(size(X,1),ntree);
index = 1:size(X,1);
X=[X index'];
for itree=1:ntree %if there is no parall computing, here then change parfor to for
 tree=forest{itree};
 cls = zeros(size(X,1),1)-100;
 clsV(:,itree)=TReadTreeBatchTest(tree,X,cls);
 %clcV(itree);
end

clsBatch = mode(clsV,2);
%clc = mode(clcV);