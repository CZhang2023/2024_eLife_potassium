% serial vvtt to cwbiast
clear variables;
files = dir('*.mat');
N = length(files);
% figure;
for i=1:N
    load(files(i).name,'tt','vv');
    if sum(vv>0)<sum(vv<0)
        vv = -vv;
    end
    
    disp(['数据处理中...',num2str(i),'/',num2str(N)]);
    fps = 1./mean(diff(tt));
    
    tw = 20;                        %time window in seconds.
    ts = 1;                         %time slide in seconds.
    fw = ceil(fps*tw);              %time window in frames.
    fs = ceil(ts*fps);              %time slide in frames.
    
    cwbias = zeros(floor((length(vv)-fw)/fs)+1,1);
    timex = cwbias;
    sv = cwbias;
    
    flag = zeros(fw,1);
    for j=1:length(cwbias)
        vvp = vv(1+(j-1)*fs:(j-1)*fs+fw);
        [y,x] = hist(vvp(vvp>0),0:0.1:100);
        sv(j) = x(find(y==max(y),1,'first'));
        
        high = sv(j);
        low = -high;
        flag=binswitch(vv(1+(j-1)*fs:(j-1)*fs+fw),low*2/3,high*2/3); 
        [ccw,cw] = interval(flag,tt(1+(j-1)*fs:(j-1)*fs+fw));
        ccw(ccw==0) = [];
        cw(cw==0) = [];
        cwbias(j) = sum(cw)/(sum(cw)+sum(ccw));
        timex(j) = mean(tt(1+(j-1)*fs:(j-1)*fs+fw));
    end
    %{
    clf;
    yyaxis left;
    plot(tt,vv);
    hold on;plot(timex,sv,'r-','linewidth',1.5);hold off;
    yyaxis right;
    plot(timex,cwbias,'g-','linewidth',1.5);
    
    saveas(gcf,[files(i).name(1:end-3),'jpg']);
    %}
    
    save(files(i).name,'timex','cwbias','tw','ts','sv','-append');
    clearvars -except i N files;
end