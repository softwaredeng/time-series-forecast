
function [tree] = TMakeTree(depth,pre_l,X,cls,nVar,nMean,nSlope,sampleModeWSZ,sampleModePos,alpha,minWin,maxWin)
depth=depth+1;
tree.terminal=0;%not a terminal node

% if the node's are pure, this will also include the leaf nodes only
% consists of 1 data
 if(size(unique(cls),1)==1)
     tree.terminal=1; % ternimal node
     tree.depth=depth;
     %tree.data=data0;
     tree.ndata=size(X,1);
     tree.oriEntropy = IGEntropy(cls);
     tree.class=mode(cls);%the most frequent class; if all equal frequency, then select the first class
     return ;
 end
 
 if(depth>=pre_l)
    tree.terminal=1; % ternimal node
    tree.depth=depth;
    tree.ndata=size(X,1);
    tree.oriEntropy = IGEntropy(cls);
    tree.class=mode(cls);%the most frequent class; if all equal frequency, then select the first class
    %tree.wholesize=wholesize;
    return ;
 end
 
bestIG = 0; bestType=-1; %0:mean; 1:slope 2:var; 
bestV = -1000; bestP1 = -1; best2=-1; 
bestWinsz = -1;

%  X0=X;
%  X=X(1:50,1:20);
%  cls = cls(1:50);
% temp1 =what;
% dir0=temp1.path;
% fid = fopen([dir0 '\x.txt'], 'w');
% for i = 1:size(X,1)
%     for j = 1:size(X,2)
%         fprintf(fid, '%f,', X(i,j));
%     end
% end
% fclose(fid);
% fid = fopen([dir0 '\y.txt'], 'w');
% for i = 1:size(cls,1)
%         fprintf(fid, '%f,', cls(i));
% end
% fclose(fid);
% size(X)
% (mean(X(:,1:18)'))'

outDouble=zeros(1,2)-1;outInt=zeros(1,4)-1;
X1 = X'; vX = X1(:);% each time seires is in the same block
nrow=size(X,1);ncol = size(X,2);
%tic
randSeed = round(rand*10000000); 
%sampleMode = 3; %1: sprt win size, all posotions; 2: sqrt win size, sqrt position  3: all win size, all position
%sampleModeWSZ = 2; %window size sample: 0: log2;  1: sqrt; 2: all
%sampleModePos = 2; %position sample:              1: sqrt; 2: all

[outDouble,outInt]=slideWin(vX',cls',nrow,ncol,nMean,nSlope,nVar,randSeed,sampleModeWSZ,sampleModePos,alpha,minWin,maxWin);

if(outDouble(end)>0)
    debug = 1;
end

if (outDouble(2)==-10000)
    a=1
    display 'no entropy gain'
    tree.terminal=1; % ternimal node
    tree.depth=depth;
    tree.ndata=size(X,1);
    tree.oriEntropy = IGEntropy(cls);
    tree.class=mode(cls);%the most frequent class; if all equal frequency, then select the first class
    return ;
end

bestIG=outDouble(1);bestV=outDouble(2);
bestGainOrig =outDouble(3);%the original information gain
bestType=outInt(4);
bestP1=outInt(2);bestP2=outInt(3);
bestWin = outInt(1);


%[bestIG1,bestType1,bestV1,bestP11,bestP21]=IGSlideWin(X,cls,nVar,nMean,nSlope);
%toc
%     tempX = X(:,1:18);    
%     slopeThis=zeros(1,size(tempX,1));
%     p=[1:18];
%     polyfit(p,tempX(1,:),1);
%     slopeThis=( p*tempX'-sum(p)*(mean(tempX'))  )/(p*p'-sum(p)*mean(p));
 if(bestGainOrig<=0)%if there is no entropy gain
    tree.terminal=1; % ternimal node
    tree.depth=depth;
    tree.ndata=size(X,1);
    tree.oriEntropy = IGEntropy(cls);
    tree.class=mode(cls);%the most frequent class; if all equal frequency, then select the first class
    %tree.wholesize=wholesize;
    return ;
 end
 
tempX = X(:,bestP1:bestP2);
if(bestType==0)%mean
  if(size(tempX,2)==1)    
      display 'mean based on 1 point'
      stat = tempX';
  else
      stat = mean(tempX'); 
  end     
end

if(bestType==2)%slope
    %display 'slope used'
    if(size(tempX,2)==1)
        display 'slope based on 1 point'
        stat = tempX';
    else        
    slopeThis=zeros(1,size(tempX,1));
    p=[bestP1:bestP2];
    slopeThis=( p*tempX'-sum(p)*(mean(tempX')) )/(p*p'-sum(p)*mean(p));    
    %toc;
  
    %tic;
    %slopeThis1=[];
    %for j = 1: size(tempX,1)
    %temp = polyfit([bestP1:bestP2],tempX(j,:),1);
    %slopeThis1 = [slopeThis1 temp(1)];
    %end
    %toc;
    %[slopeThis;slopeThis1]    
    %debug = atan(slopeThis);
   stat = atan(slopeThis);
   end
end
if(bestType==1)%var
  if(size(tempX,2)==1)
        display 'var based on 1 point'
        stat = tempX';
  else
  stat = sqrt(var(tempX'));  
  %stat = (mean(abs(diff(diff(tempX')))));
  end
end



Xleft=X(stat<bestV+10^(-10),:);
Xright=X(stat>bestV,:);

if (size(Xleft,1)==0)
    nullSplit=1
end
if (size(Xright,1)==0)
    nullSplit=1
end

%clsleft=cls(stat<=bestV);
clsleft=cls(stat<bestV+10^(-10));
clsright=cls(stat>bestV);

if(size(clsleft,1)*size(clsright,1)==0)
    debug=1
end


tree.split=bestV;
tree.entropy=bestIG;
tree.depth=depth;
%tree.feature=fea;
tree.ndata=size(X,1);
tree.p1 = bestP1;
tree.p2 = bestP2;
tree.stat = bestType;
tree.split = bestV;
tree.oriEntropy = IGEntropy(cls);
tree.minChild = min(size(clsleft,1),size(clsright,1));
%tree.left_ndata=size(dataleft,1);tree.right_ndata=size(dataright,1);
%tree.wholesize=wholesize;


[tree.childl]=TMakeTree(depth,pre_l,Xleft,clsleft,nVar,nMean,nSlope,sampleModeWSZ,sampleModePos,alpha,minWin,maxWin);
[tree.childr]=TMakeTree(depth,pre_l,Xright,clsright,nVar,nMean,nSlope,sampleModeWSZ,sampleModePos,alpha,minWin,maxWin);

