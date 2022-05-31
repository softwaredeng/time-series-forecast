clear;


M = 1000;
V = [];
for t = 1:1000
    V = [V t*(M-t+1)];
end


b=plot(V,'LineWidth',2);
set(gca,'fontsize',14);
set(get(b,'Parent'),'FontSize',14)  