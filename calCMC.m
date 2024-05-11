
clear variables;clc;
files = dir('Track_*.mat');
%{
N = length(files);
k=0;
for i = 1:N
    k=k+1;
    load(files(i).name,'pos','fps');
    bacpos = pos;
    bact = ceil(bacpos(:,3)/fps);
    t = unique(bact);
    CMC = zeros(size(t));
    cellnum = zeros(size(t));
    parfor j=1:length(t)
         CMC(j) = mean(-(bacpos(bact==t(j),1)-705)/705);
         cellnum(j) = sum(bact==t(j));
    end
    save(files(i).name,'t','CMC','-append');
    eval(['CMC',num2str(k) ,'= CMC;']);
end
clearvars -except files;
%}
CMCall = [];
yall = [];
xall = [];
tall = [];
for i = 1:length(files)
    load(files(i).name,'CMC','t','x','y');
    CMCall = [CMCall,CMC];
    tall = [tall,t];
    
    yall = [yall,y'];
    xall = [xall,x'];
end
mCMC = mean(CMCall,2);
semCMC = std(CMCall,1,2)/sqrt(size(CMCall,2));
mt = mean(tall,2);

subplot(211)
shadedErrorBar(mt,mCMC,semCMC,'b');
set(gca,'fontsize',16);
ylim([0.05 0.55]);
xlabel('Time (s)');
ylabel('CMC');

my = mean(yall,2);
semy = std(yall,1,2)/sqrt(size(yall,2));
mx = mean(xall,2);
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
[fitresult, gof] = fit( mx, my, 'exp1', opts );

subplot(212)
hold on;
errorbar(mx,my,semy,'bo','markerfacecolor','b','markersize',6.0);
xi = 0:1:(1410*5.5/20);
plot(xi,fitresult(xi),'r-','linewidth',2.0);
hold off;

set(gca,'fontsize',16);
xlabel('x (\mum)');
ylabel('Probability Density (\mum)^{-1}');


