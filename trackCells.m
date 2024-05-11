
% 批量追踪细菌
%{
clear variables;clc;
files = dir('HCB*.mat');
for ifile=1:length(files)
    tic
    disp(['Tracking ...',num2str(ifile),'/',num2str(length(files))]);
    
    load(files(ifile).name,'bacpos','fps');
    %%%%%-----bactierium position-----%%%%%
    
    tracks = trackcell2D(bacpos,fps,5.5,50,20);  
    
    save(['Track_',files(ifile).name],'bacpos','fps','tracks');
    toc
clearvars -except files ifile;
end
clear variables;clc;
%}

%追踪后轨迹处理，剔除不动的细菌（v<5um/s）
file = dir('Track_*.mat');
for ifile = 1:length(file)
    load(file(ifile).name,'tracks','fps');
    
% 计算平均速度.单个点的数据定义速度为-1.
clc;
for i = 1:length(tracks)
    if length(tracks(i).x)>1
        x = tracks(i).x*5.5/20;
        y = tracks(i).y*5.5/20;
        v(i) = mean(sqrt(diff(y).^2+diff(x).^2)*fps);
    else
        v(i)=-1;
    end
end
% 去除单点与速度为0的轨迹点
tracks(v<=0) = [];
v(v<=0) = [];
clearvars -except tracks v ifile file rr fps;

% 画出轨迹图
%{
clc;
t = tracks(v>=5);
hold on;
for i = 1:length(t)
    plot(t(i).x,t(i).y,'-.');
end
hold off;
clearvars -except tracks l v stdp;
%}

t = tracks(v>5);
pos = [cell2mat({t(:).x}'),cell2mat({t(:).y}'),cell2mat({t(:).frame}')];
xp = pos(pos(:,3)>=299*fps,1)*5.5/20;
[c,x] = hist(xp,30);
y = c./sum(c)/(mean(diff(x)));

opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
[fitresult, gof] = fit( x', y', 'exp1', opts );

figure;
hold on;
bar(x,y);
xi = 0:1:(1410*5.5/20);
plot(xi,fitresult(xi),'r-','linewidth',2.0);
hold off;
title(num2str(fitresult.b,'%.4f'));

rr(ifile) = fitresult.b;
save(file(ifile).name,'pos','fitresult','x','y','-append');
clearvars -except ifile file rr;
end