
function [cls]=TReadTreeBatchTest(tree,X,cls)
%length=length+1;

% if the node is a terminal, then return
if(tree.terminal==1)  
    index = X(:,end);
    cls(index,:)=tree.class; 
return  
end


tempX = X(:,tree.p1:tree.p2);
if(tree.stat==0)%mean
    testStat = mean(tempX');    
end
 
if(tree.stat==2)%slope
    %temp = polyfit([tree.p1:tree.p2],X(tree.p1:tree.p2),1);
    %testStat = temp(1);
    p=[tree.p1:tree.p2];
    testStat=( p*tempX'-sum(p)*(mean(tempX')) )/(p*p'-sum(p)*mean(p));
end

if(tree.stat==1)%std
    %temp = polyfit([tree.p1:tree.p2],X(tree.p1:tree.p2),1);
    %testStat = temp(1);
    testStat = sqrt(var(tempX')); 
end

Xleft=X(testStat<=tree.split,:);
Xright=X(testStat>tree.split,:);

if(size(Xleft,1)>0)
[cls]=TReadTreeBatchTest(tree.childl,Xleft,cls);
end
if(size(Xright,1)>0)
[cls]=TReadTreeBatchTest(tree.childr,Xright,cls);
end


%return;
end
