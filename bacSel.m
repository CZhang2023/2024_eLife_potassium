clear variables;clc;

%
% Step 1: 预选范围和阈值.
files = dir('*.avi');
centp = zeros(length(files),1);
thp = zeros(length(files),1);
for ifile = 1:length(files)
    disp(['Step1/3 视频阈值区域选取...',num2str(ifile),'/',num2str(length(files))]);

    obj = VideoReader(files(ifile).name);
    I = read(obj,1);
    I(:,:,2:3) = [];
    [~, rect] = imcrop(I);
    centp(ifile) = rect(1)+fix(rect(3)/2); % the center of channel in pixel.
    close;
    
    J = I(:,(centp(ifile)-705):(centp(ifile)+705));
    se=strel('disk',13);
    J1 = imsubtract(imadd(J,imtophat(J,se)), imbothat(J,se)); % 消除不均匀背景
    I1 = imadjust(J1);
    
    judge=0;
    while judge==0
        imshow(I1);
        th = input('threshold selection:\n');
        bw = imbinarize(I1,th);
        se = strel('disk',2);
        bw = imopen(bw,se);
        bw = imclose(bw,se);
        imshow(bw);
        judge = input('是否存储当前阈值(1/0)？\n');
        close all;
    end
    thp(ifile) = th;
clearvars -except files ifile thp centp;
end
save('thresholdandrange.mat','files','thp','centp');


% Step 2: 批量处理图像.
clear variables;clc;
load('thresholdandrange.mat');
for ifile = 1:length(files)
    tic
    disp(['Step2/3 视频处理...',num2str(ifile),'/',num2str(length(files))]);
    
    obj = VideoReader(files(ifile).name); %#ok<*TNMLP>
    N = obj.Duration*obj.FrameRate;
    frames = zeros([1000, 1411, 9000],'logical');
   
    fps = obj.FrameRate;
    centofchannel = fix(centp(ifile));
    th = thp(ifile);
    for i = 1:9000
        I = read(obj,i); %#ok<*VIDREAD>
        I(:,:,2:3) = [];
        J = I(:,(centofchannel-705):(centofchannel+705));
        se=strel('disk',13);
        J1 = imsubtract(imadd(J,imtophat(J,se)), imbothat(J,se)); % 消除不均匀背景
        I1 = imadjust(J1);
        bw = ~imbinarize(I1,th);
        se = strel('disk',2);
        bw = imopen(bw,se);
        bw = imclose(bw,se);
        frames(:,:,i) = bw;
%         clearvars -except files ifile obj i N thp centp frames fps centofchannel th ;
    end
    save([files(ifile).name(1:end-3),'mat'],'frames','fps','th','centofchannel');
    toc
clearvars -except files ifile thp centp;
end
%}

% Step 3: 批量定位细胞.
clear variables;clc;
files = dir('HCB*.mat');
for ifile=1:length(files)
    tic
    disp(['Step3/3 细胞定位...',num2str(ifile),'/',num2str(length(files))]);
    
    load(files(ifile).name,'frames');
    %%%%%-----bactierium position-----%%%%%
    frame(1:size(frames,3)) = struct('Pos',[]);
    for i=1:size(frames,3)
        Labtemp = bwlabeln(frames(:,:,i),8);
        Strtemp = regionprops(Labtemp, 'Area');
        BW = ismember(Labtemp, find([Strtemp.Area]<=100&[Strtemp.Area]>=5));
        Labtemp = bwlabeln(BW,8);
        Strtemp = regionprops(Labtemp, 'Centroid');
        frame(i).Pos(:,1:2) = cell2mat((cell({Strtemp(:).Centroid}))');
        frame(i).Pos(:,3) = i;
        clear *temp;
    end
    bacpos = cell2mat((cell({frame(:).Pos}))');
    save(files(ifile).name,'bacpos','-append');
    toc
clearvars -except files ifile;
end
clear variables;clc;