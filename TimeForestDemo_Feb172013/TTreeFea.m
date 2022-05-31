
function [R]=TTreeFea(tree,R)
%length=length+1;

% if the node is a terminal, then return
if(tree.terminal==1)  
    found=1; 
    clc=tree.class; 
return ; 
end
temp = [tree.stat tree.entropy tree.split tree.p1 tree.p2 tree.ndata tree.depth tree.minChild];
R=[R;temp];
R=TTreeFea(tree.childl,R);
R=TTreeFea(tree.childr,R);
%length=length+1;
%tree.split=bestV;
%tree.entropy=bestIG;
%tree.length=length;
%tree.p1 = bestP1;
%tree.p2 = bestP2;
%tree.stat = bestTpe;

end
