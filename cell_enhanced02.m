clear;close all;clc
tic
runtime = [];
baseFileName = 'data/A172/';
imagename = 'A172_Phase_A7_1_00d00h00m_1.tif';
maskname =strcat(imagename(1:end-4),'_instance_masks.png');
%  maskname =strcat(imagename(1:end-4),'_cp_masks.png');
fullFileName = fullfile(baseFileName,imagename);
fullFileName2 = fullfile(baseFileName,maskname);
input = imread(fullFileName);
mask = imread(fullFileName2);
% input = imread('20221116 4hpt_BF-7_ch00.tif');
if ndims(input) > 2
    input = rgb2gray(input);    
end
% figure,histogram(input)
fsize = 15;
figure,imshow(input,[])
title(strcat('input',imagename(1:end-4)),'FontSize',fsize)
% print(figure(1),'resultpic3_C7_1/C7_1_input.png','-dpng')
input01 = single(input);
%% 该双树复小波变换将图像在4个尺度上分解，每个尺度上产生为6个高通子带和1个低通子带
lev = 4;
[Yl,Yh,Yscale] = dualtree2(input,'Level',lev,'LevelOneFilter','legall','FilterLength',6);

Yh02 = Yh;
% Level：分解的尺度，最大值为log2(min(H,W))
% LevelOneFilter：用于一级分析的双正交滤波器，选用legall对高频第一级分析得到的噪声最少；
% FileterLength：用于2级及更高级别的正交Hilbert Q-shift分析滤波器的长度，值越小，噪声越低；
% 最后一个级别的低频分量已经包含了大部分能量信息；
% 高频分量包含了图像的大部分细节信息。
% Yl是最后一个级别（尺度）的低频分量，而Yh包含了每个级别（尺度）的高频分量
% Yl = Yscale{1,1};
% figure % 低通效果
% imagesc(Yl)
% colormap gray
% % set(gca,'XTick',[]),set(gca,'YTick',[])
% figure % 高通效果
% montage(imag(Yh{1,1}))
%% scaling（低通）系数增强方法
Yscale02 = Yscale;
for m = 1:lev
    Yl = Yscale{m,1};
    strel1 = strel('disk',2); % 创新点在结构元素上
    Yl_b = imbothat(Yl,strel1); % 对低频部分进行底帽变换，因为原图中是需要增强细胞内部的暗部分：Iclose-I
    Yl_t = imtophat(Yl,strel1); % 增强图像中亮的部分：I-Iopen
    Yl02 = Yl + Yl_t - Yl_b ; % 原图 + 顶帽变换 - 底帽变换：result = I + (I - Io) - (Ic - I) = 3I - Io - Ic
    Yscale02{m,1} = Yl02;
end
% figure% 低频增强后的效果
% imagesc(Yl02)
% colormap gray
%% wavelet系数增强方法：创新点在阈值选取
% Yhgains = zeros(lev,6); 
% Yhgains02 = Yhgains; Yhgains02(:,[1 3 5]) = 1;
% Yhgains03 = rand(lev,6);
Yhgains04 = ones(lev,6);
% Yhgains05 = ones(1,6);
% 定义结构元素序列
B0 = strel('disk',3); % 影响后续分割结果
B1 = strel('disk',3); % 影响后续分割结果
B2 = strel('disk',2);
B3 = strel('disk',1);
B = {B0,B1,B2,B3};
[~,~,zidai] = size(Yh{1,1}); % 每个尺度高通子带的个数
sigma_f = 0.5; % 高斯滤波器的标准差：影响后续分割结果中的突刺
fs = 1; % 高斯滤波器的大小：2*ceil（2*sigma_f）+1
fd = 'spatial'; % 滤波域：可选择为'spatial'或者'frequency'
for m = 1:lev
    for n = 1:zidai
        switch n
            case 1 % HL
                rYh = real(Yh{m,1}(:,:,n)); % 每个尺度每个子带的实部
                rYh02 = denoise(rYh); % 每个尺度每个子带的实部进行降噪
                rYh02 = imgaussfilt(rYh02,sigma_f,"FilterSize",fs,'FilterDomain',fd); % 高频信号处理
                iYh = imag(Yh{m,1}(:,:,n)); % 每个尺度每个自带的虚部
                iYh02 = Lenhanced(iYh,B{1,m}); % 低频信号增强处理
                iYh02 = imgaussfilt(iYh02,sigma_f,"FilterSize",fs,'FilterDomain',fd);
                Yh02{m,1}(:,:,n) = rYh02 + 1i*iYh02;
            case 2 % HH
                rYh = real(Yh{m,1}(:,:,n));
                rYh02 = denoise(rYh);
                rYh02 = imgaussfilt(rYh02,sigma_f,"FilterSize",fs,'FilterDomain',fd);
                iYh = imag(Yh{m,1}(:,:,n));
                iYh02 = denoise(iYh);
                iYh02 = imgaussfilt(iYh02,sigma_f,"FilterSize",fs,'FilterDomain',fd);
                Yh02{m,1}(:,:,n) = rYh02 + 1i*iYh02;
            case 3 % LH
                rYh = real(Yh{m,1}(:,:,n));
                rYh02 = Lenhanced(rYh,B{1,m});
                rYh02 = imgaussfilt(rYh02,sigma_f,"FilterSize",fs,'FilterDomain',fd);
                iYh = imag(Yh{m,1}(:,:,n));
                iYh02 = denoise(iYh);
                iYh02 = imgaussfilt(iYh02,sigma_f,"FilterSize",fs,'FilterDomain',fd);
                Yh02{m,1}(:,:,n) = rYh02 + 1i*iYh02;
            case 4 % LH
                rYh = real(Yh{m,1}(:,:,n));
                rYh02 = Lenhanced(rYh,B{1,m});
                rYh02 = imgaussfilt(rYh02,sigma_f,"FilterSize",fs,'FilterDomain',fd);
                iYh = imag(Yh{m,1}(:,:,n));
                iYh02 = denoise(iYh);
                iYh02 = imgaussfilt(iYh02,sigma_f,"FilterSize",fs,'FilterDomain',fd);
                Yh02{m,1}(:,:,n) = rYh02 + 1i*iYh02;
            case 5 % HH
                rYh = real(Yh{m,1}(:,:,n));
                rYh02 = denoise(rYh);
                rYh02 = imgaussfilt(rYh02,sigma_f,"FilterSize",fs,'FilterDomain',fd);
                iYh = imag(Yh{m,1}(:,:,n));
                iYh02 = denoise(iYh);
                iYh02 = imgaussfilt(iYh02,sigma_f,"FilterSize",fs,'FilterDomain',fd);
                Yh02{m,1}(:,:,n) = rYh02 + 1i*iYh02;
            case 6 % HL
                rYh = real(Yh{m,1}(:,:,n));
                rYh02 = denoise(rYh);
                rYh02 = imgaussfilt(rYh02,sigma_f,"FilterSize",fs,'FilterDomain',fd);
                iYh = imag(Yh{m,1}(:,:,n));
                iYh02 = Lenhanced(iYh,B{1,m});
                iYh02 = imgaussfilt(iYh02,sigma_f,"FilterSize",fs,'FilterDomain',fd);
                Yh02{m,1}(:,:,n) = rYh02 + 1i*iYh02;
        end
    end
end
output = idualtree2(Yl02,Yh02,'DetailGain',Yhgains04,'LowpassGain',1);
% output = 2*output - input01;
figure,imshow(output,[])
title('output','FontSize',fsize)

% imwrite(output,'A172_Phase_C7_1_00d00h00m_1_output.tif')

output = round(output); % 将其转换为整数型
output = uint8(output); % 将其转换为8位整型


% 图像增强定量比较：均方误差immse，峰值信噪比psnr，结构相似性ssim
% err01 = immse(input,output);
% fprintf('\n 均方误差为%0.4f\n',err01)
% [peaknr,snr] = psnr(input,output);
% fprintf('\n 峰值信噪比为%0.4f\n',peaknr)
% fprintf('\n 信噪比为%0.4f\n',snr)
% [ssimval,ssimmap] = ssim(input,output);
% fprintf('\n 全局结构相似性为%0.4f\n',ssimval)

% figure,imshow(ssimmap,[])
% title("Local SSIM Map with Global SSIM Value: "+num2str(ssimval))
% runtime(1) = toc 

%% 使用非锐化掩膜锐化图像
% output01 = imsharpen(output,'Radius',2,'Amount',1);
% figure,imshow(output01,[])
%% 对比度调整
o_scale = stretchlim(output); % 对比度拉伸的调整范围
output01 = imadjust(output,o_scale); % 调整对比度
figure,imshow(output01,[]),title('对比度调整后')
% imwrite(output01,'A172_Phase_C7_1_00d00h00m_1_output01.tif')

% [peaknr,snr] = psnr(uint8(input),output01);
% fprintf('\n 峰值信噪比为%0.4f\n',peaknr)
% fprintf('\n 信噪比为%0.4f\n',snr)
% figure,subplot(131),imshow(input),title('input')
% hold on
% subplot(132),imshow(output),title('output')
% hold on
% subplot(133),imshow(output01),title('imadjust')
%% 双边滤波
% output01 = imbilatfilt(output);

%% 梯度阈值分割（EGT）
% [output03,output03_mask] = segmentImage2(output);
% figure,imshow(input),hold on
% visboundaries(output03)
% 判定去除小面积条件
area_min=150;
min_object_size=area_min; % 低于最小值的物体将会被去除
min_hole_size=322; % 大于这个值的孔将会被填充
treshold_finetune=0; % 贪婪系数，值越大，包含的区域越多
[output03,threshold] = my_EGT_Segmentation(output, min_object_size,min_hole_size,treshold_finetune);
% imwrite(output03,'A172_Phase_C7_1_00d00h00m_1_output003.tif')
%% 定义一个圆形滤波器，使得分割结果光滑
% dr = 5;
% for k = 1:dr
%     d = fspecial('disk',k);
%     temp = imfilter(output03,d,'replicate','corr');
%     figure,imshow(temp,[])
% end
d = fspecial('disk',5); % 定义圆形滤波器
output03 = imfilter(output03,d,'replicate','corr');

% imwrite(output03,'A172_Phase_C7_1_00d00h00m_1_output003_2.tif')
% output03_1 = labeloverlay(input,bwlabel(output03),'Transparency',0.5);
% imwrite(output03_1,'A172_Phase_C7_1_00d00h00m_1_output03_1.tif')
%% 非细胞区域去除
output03 = imfill(output03,'holes'); % 填充孔洞
output03_p = regionprops(output03, {'Area', 'MajorAxisLength', 'MinorAxisLength'});
output03_p_a = cat(1,output03_p.Area); % 面积序列
output03_p_a_m = median(output03_p_a); % 面积的中位数
output03_p_a_u = mean(output03_p_a,'all'); % 面积的均值
output03_l = bwlabel(output03); % 标签图像
a_number = max(output03_l,[],'all'); % 区域面积个数
output033 = zeros(size(output03)); % 结果掩膜
for m = 1:a_number
    output03_m = zeros(size(output03)); % 标签掩膜
    [a_r,a_c] = find(output03_l==m);
    output03_m(sub2ind(size(input),a_r,a_c)) = 1;
    output03_m_p = regionprops(output03_m, {'Area', 'MajorAxisLength', 'MinorAxisLength'});
    % 只有长轴/短轴的比值大于3且面积小于全部区域的面积均值时，该区域才会置为0
    while output03_m_p.Area >= output03_p_a_u || round(output03_m_p.MajorAxisLength) / round(output03_m_p.MinorAxisLength) <= 3
        output033(sub2ind(size(output03),a_r,a_c)) = 1;
        break
    end
end
%

%% 标记分水岭
Yl03 = Yscale{3,1};
% imwrite(Yl03,'A172_Phase_C7_1_00d00h00m_1_Yl03.png')
L = imsegkmeans(single(Yl03),2,'NumAttempts',2);
% 防止聚类的时候，将细胞目标作为背景
N = find(L==2);
[Yl03_x,Yl03_y] = size(Yl03);
if numel(N) > Yl03_x*Yl03_y/2
    L = imcomplement(L);
end
bw = L == max(L,[],'all'); % 取最大值
bw2 = labeloverlay(uint16(Yl03),bwlabel(bw),'Transparency',0.2);

%
bw_cen = regionprops(bw,'Centroid');
bw_centroids = cat(1,bw_cen.Centroid);
marked_centroids = round(bw_centroids.*4);
mask01 = zeros(size(input));
mask01(sub2ind(size(input),marked_centroids(:,2),marked_centroids(:,1))) = 1;
mask01 = imdilate(mask01,strel('disk',5,8)); % 获取种子点
D1 = -bwdist(~mask01);
l = watershed(D1); % 分水岭变换


% output04 = output03; % 未经过处理
% output04(l==0) = 0;
% output04 = bwareaopen(output04,70);
% figure,imshow(input)
% hold on
% visboundaries(output04)

%% 原图直接处理
% [output02,output02_masked] = segmentImage(input);
% output04 = logical(output02); % 经过面积去除的操作
% output04(l==0) = 0;
% output04 = bwareaopen(output04,70);
% figure,imshow(input),hold on
% visboundaries(output04)
% title('原图直接处理')


output04 = logical(output033); % 经过面积去除的操作
%% 加一个判定条件，如果种子点的个数远远超过两倍的经过EGT分割后的区域个数，则不采用该种子点
[~,N_egt] = bwlabel(output033); % EGT分割后的区域个数
[~,N_seed] = bwlabel(mask01); % 种子点个数
if N_seed < 2*N_egt
    output04(l==0) = 0;
end



%% 
output04 = bwareaopen(output04,150);
% 分割效果是提升了，但是IoU降低了将近0.10
output04_2 = imerode(output04,strel('disk',1,0));

d1 = fspecial('disk',3);  % 值太大会导致重新粘连
output04_3 = imfilter(output04_2,d1,'replicate','corr'); % 分割后使用盘形滤波器平滑后重叠或者粘连细胞问题根本解决不了

output04_3 = bwareaopen(output04_3,150);
% imwrite(output04_3,'A172_Phase_C7_1_00d00h00m_1_output04.png')
% output04_33 = labeloverlay(input,bwlabel(output04_3),"Transparency",0.5);
% figure,imshow(output04_33,[])
% imwrite(output04_33,'resultvis/C7_1.png')

figure,imshow(input),hold on
visboundaries(output04_3)
title('Proposed')
jaccard1 = jaccard(output04,logical(mask));
% fprintf('GT和所提出算法的jaccard系数是%.4f\n',jaccard1)
[score1,precision1,recall1] = bfscore(output04,logical(mask));
% fprintf('GT和所提算法的精确率Precision为%.2f\n召回率Recall为%.2f\nF1-score为%.2f\n', ...
%     precision1,recall1,score1)

runtime = toc;
% print(figure(4),'resultpic3_A7_4/A7_4_proposed.png','-dpng'); % 将图像保存为向量图格式PDF

% 
% figure,subplot(221),imshow(input),title('original')
% ,hold on
% subplot(222),imshow(output01),title('output'),hold on
% subplot(223),imshow(output04_3),title('EGT-output')
% subplot(224),imshow(input),hold on,
% visboundaries(output04_3),title('mask')
% output04 = output02;
% output04(l==0) = 0;
% output04 = bwareaopen(output04,70);
% figure,imshow(input)
% hold on
% visboundaries(output04)
% %% 超像素分割
% [L01,N01] = superpixels(output,2000);
% BW = boundarymask(L01);
% figure,imshow(imoverlay(input,BW,'cyan'),'InitialMagnification',100)
% figure,imshow(L01,[])
%% EGT

area_min=250;
min_object_size=area_min;
min_hole_size=322;
treshold_finetune=2.2;
res_EGT = EGT_Segmentation(input, min_object_size,min_hole_size,treshold_finetune);
figure;imshow(input,[]);hold on;
visboundaries(res_EGT)
title('EGT')
jaccard2 = jaccard(res_EGT,logical(mask));
% fprintf('EGT的jaccard为%.4f\n',jaccard2)
[score2,precision2,recall2] = bfscore(res_EGT,logical(mask));
% fprintf('GT和EGT算法的精确率Precision为%.2f\n召回率Recall为%.2f\nF1-score为%.2f\n', ...
%     precision2,recall2,score2)
% print(figure(5),'resultpic3_A7_4/A7_4_EGT.png','-dpng'); % 将图像保存为向量图格式PDF
% res_EGT(l==0) = 0;
% figure;imshow(input,[]);hold on;
% visboundaries(res_EGT)
% title('EGT2')


rmin=20;%estimated celll radius range and area
rmax=56;
area_min=250;
min_and_max_of_img1=[0 2.7198]; % minimal and maximal intensity of first image (for normalization)
%% simple treshold
% t=0.03;
% I_norm=mat2gray(input,min_and_max_of_img1);
% res_simple_treshold=I_norm>t;
% figure;imshow(input,[]);hold on;
% visboundaries(res_simple_treshold)
% title('simple treshold')

%% otsu treshold 
% I_norm=mat2gray(input);
% t = graythresh(I_norm);
% res_otsu_treshold=I_norm>t;
% figure;imshow(input,[]);hold on;
% visboundaries(res_otsu_treshold)
% title('otsu treshold')

%% poisson treshold
I_norm=imadjust(input);%works better after elimination of extreme values
t=poisson_tresh(I_norm);
res_poisson_treshold=I_norm>t;
figure;imshow(input,[]);hold on;
visboundaries(res_poisson_treshold)
title('poisson treshold')
jaccard3 = jaccard(res_poisson_treshold,logical(mask));
% fprintf('poisson的jaccard为%.4f\n',jaccard3)
[score3,precision3,recall3] = bfscore(res_poisson_treshold,logical(mask));
% fprintf('GT和poisson算法的精确率Precision为%.2f\n召回率Recall为%.2f\nF1-score为%.2f\n', ...
%     precision3,recall3,score3)
% print(figure(6),'resultpic3_A7_4/A7_4_poisson.png','-dpng'); % 将图像保存为向量图格式PDF
%% sJuneau
% window_size=3;
% t=0.01;
% min_object_size=area_min;
% I_norm=mat2gray(input,min_and_max_of_img1);
% [res_juneau] = juneau(I_norm,window_size,t, min_object_size);
% figure;imshow(input,[]);hold on;
% visboundaries(res_juneau)
% title('sJuneau')


%% level-set Chan-Vese
smoothness=0.2;
aditional_force=0.17;
I_norm=mat2gray(input,min_and_max_of_img1);
initialization=imdilate(res_EGT,strel('disk',5));
res_chanvese=activecontour(I_norm,initialization,500,'Chan-Vese','ContractionBias',-aditional_force,'SmoothFactor',smoothness);
figure;imshow(input,[]);hold on;
visboundaries(res_chanvese)
title('sLS-Chan-Vese')
jaccard4 = jaccard(res_chanvese,logical(mask));
% fprintf('chanvese的jaccard为%.4f\n',jaccard4)
[score4,precision4,recall4] = bfscore(res_chanvese,logical(mask));
% fprintf('GT和chanvese算法的精确率Precision为%.2f\n召回率Recall为%.2f\nF1-score为%.2f\n', ...
%     precision4,recall4,score4)
% print(figure(7),'resultpic3_A7_4/A7_4_chanvese.png','-dpng'); % 将图像保存为向量图格式PDF
%% sJuneau
% window_size=5;
% t=0.02;
% min_object_size=area_min;
% I_norm=mat2gray(input,min_and_max_of_img1);
% [res_juneau] = juneau(I_norm,window_size,t, min_object_size);
% figure;imshow(input,[]);hold on;
% visboundaries(res_juneau)
% title('sJuneau')

%% ilastik
% ilastik_path = 'E:\Documents\Research_Topic\Image_segmentation\Code\双树复小波\data\A172_ilastik\'; % ilastik图像文件路径
% ilastik_pic_name = imagename(1:end-4); % 提取所处理图像的前缀名字
% ilastik_pic_name2 = strcat(ilastik_pic_name,'_ilastik.png'); % 获取该图像的完整名称
% ilastik_pic = imread(strcat(ilastik_path,ilastik_pic_name2));% 读取该路径下的图片
% figure,imshow(input,[]),hold on
% visboundaries(ilastik_pic)
% title('ilastik')
% 
% jaccard5 = jaccard(logical(mask),ilastik_pic);
% fprintf('ilastik的jaccard为%.4f\n',jaccard5)
% [score5,precision5,recall5] = bfscore(ilastik_pic,logical(mask));
% fprintf('GT和ilastik算法的精确率Precision为%.2f\n召回率Recall为%.2f\nF1-score为%.2f\n', ...
%     precision5,recall5,score5)

%% Jaccard
jaccard_path = 'data\A172_jaccard\'; % ilastik图像文件路径
jaccard_pic_name = imagename(1:end-4); % 提取所处理图像的前缀名字
jaccard_pic_name2 = strcat(jaccard_pic_name,'_instance_masks.png'); % 获取该图像的完整名称
jaccard_pic = imread(strcat(jaccard_path,jaccard_pic_name2));% 读取该路径下的图片
figure,imshow(input,[]),hold on
visboundaries(jaccard_pic)
title('jaccard')

jaccard5 = jaccard(jaccard_pic,logical(mask));
% fprintf('ilastik的jaccard为%.4f\n',jaccard5)
[score5,precision5,recall5] = bfscore(jaccard_pic,logical(mask));
% fprintf('GT和jaccard算法的精确率Precision为%.2f\n召回率Recall为%.2f\nF1-score为%.2f\n', ...
%     precision5,recall5,score5)
% print(figure(8),'resultpic3_A7_4/A7_4_mask.png','-dpng'); % 将图像保存为向量图格式PDF


%% 将mask转换为标签mask
if max(mask,[],'all') == 1
    mask = bwlabel(mask);
end
mask2 = false(size(input));
figure,imshow(input),title('mask'),hold on
result = input;
% 防止粘连现象显示不出
for k = 1:max(mask,[],'all')
    mask1 = mask;
    mask1(mask1~=k) = 0;
    mask1 = bwareaopen(mask1,100);
    mask2 = mask2 + mask1;
    visboundaries(mask1) % 利用visbounaries来显示结果
    % 保存结果
%     mask1_e = edge(mask1,'canny');
%     result = labeloverlay(result,mask1_e,"Colormap",[1 0 0],'transparency',0);
end
% print(figure(9),'resultpic3_A7_4/A7_4_mask.png','-dpng'); % 将图像保存为向量图格式PDF
%     imwrite(result,'A172.bmp') 
%     print(figure(13), 'result.eps', '-depsc'); % 将图像保存为向量图格式EPS、PDF
%     print(figure(13),'result.pdf','-dpdf'); % 将图像保存为向量图格式PDF
% print(figure(4),'resultpic/A7-propoese.png','-dpng');
% print(figure(5),'resultpic/A7-EGT.png','-dpng');
% print(figure(6),'resultpic/A7-poisson.png','-dpng');
% print(figure(7),'resultpic/A7-chan.png','-dpng');
% print(figure(8),'resultpic/A7-mask.png','-dpng');
% jaccard(logical(mask2),output03)

%% 计数指标
% Ngt = max(mask,[],'all');
% [~,Nj] = bwlabel(jaccard_pic);
% [~,Negt] = bwlabel(res_EGT);
% [~,Ncv] = bwlabel(res_chanvese);
% [~,Np] = bwlabel(output04_3);
% fprintf('GT为%d\n jaccard为%d\n EGT为%d\n CV为%d\n Proposed为%d\n',Ngt,Nj,Negt,Ncv,Np);
%% 分割评价指标：扯犊子，用mask乘以算法结果，那还用个屁的算法结果啊，直接套用mask来计算不是更好吗？？？
% result_noNaN = [];
% result_NaN20 = [];
% result_noNaN0 = [];
% %%% proposed
% pr = immultiply(mask,output04_3);
% rr1 = jaccard(double(pr),double(mask));
% rr2 = rr1;
% rr1 = rr1(~isnan(rr1)); % 去掉NaN值
% result_noNaN(1) = mean(rr1);
% rr2(isnan(rr2)) = 0; % 替换NaN值为0
% result_NaN20(1) = mean(rr2);
% rr3 = rr1;
% rr3(rr3==0) = []; % 去掉0和NaN值
% result_noNaN0(1) = mean(rr3);
% 
% %%% EGT
% er = immultiply(mask,res_EGT);
% err1 = jaccard(double(er),double(mask));
% err2 = err1;
% err1 = err1(~isnan(err1)); % 去掉NaN值
% result_noNaN(2) = mean(err1);
% err2(isnan(err2)) = 0; % 替换NaN值为0
% result_NaN20(2) = mean(err2);
% err3 = err1;
% err3(err3==0) = []; % 去掉0和NaN值
% result_noNaN0(2) = mean(err3);
% 
% %%% CV
% cr = immultiply(mask,res_chanvese);
% crr1 = jaccard(double(cr),double(mask));
% crr2 = crr1;
% crr1 = crr1(~isnan(crr1)); % 去掉NaN值
% crr3 = crr1;
% result_noNaN(3) = mean(crr1);
% crr2(isnan(crr2)) = 0; % 替换NaN值为0
% result_NaN20(3) = mean(crr2);
% crr3(crr3==0) = []; % 去掉0和NaN值
% result_noNaN0(3) = mean(crr3);
% 
% %%% Jaccard
% 
% jr = immultiply(mask,jaccard_pic);
% jrr1 = jaccard(double(jr),double(mask));
% jrr2 = jrr1;
% jrr1 = jrr1(~isnan(jrr1)); % 去掉NaN值
% jrr3 = jrr1;
% result_noNaN(4) = mean(jrr1);
% jrr2(isnan(jrr2)) = 0; % 替换NaN值为0
% result_NaN20(4) = mean(jrr2);
% jrr3(jrr3==0) = []; % 去掉0和NaN值
% result_noNaN0(4) = mean(jrr3);


%% 保存mask和算法中的每个细胞
%%% 找到mask中每一个标签的位置并保存
% for cell_v = 1:max(mask,[],'all')
%     biaoqian = double(zeros(size(mask)));
%     [biaoqian_r,biaoqian_c,biaoqian_v] = find(mask == cell_v);
%     biaoqian_s = find(mask == cell_v);  % 对应值所在位置的索引
%     biaoqian(biaoqian_s)= biaoqian_v;
%     biaoqian_folder = imagename(1:end-4); % 文件夹名称
%     if isfolder(biaoqian_folder) == 0
%         mkdir(biaoqian_folder);
% %     else
% %         disp('文件夹已存在');
%     end
%     imwrite(biaoqian,strcat(biaoqian_folder,'\',imagename(1:end-4),'-',num2str(cell_v),'.png'))
%     % figure,imshow(biaoqian,[])
% end
% 
% %%% 保存DTCWT-EGT算法的每个细胞结果
% dwegt_r = bwlabel(output04_3);
% for dwegt_v = 1:max(dwegt_r,[],'all')
%     dwegt_bq = double(zeros(size(dwegt_r)));
%     [dwegt_bq_r,dwegt_bq_c,dwegt_bq_v] = find(dwegt_r == dwegt_v);
%     dwegt_s = find(dwegt_r == dwegt_v);
%     dwegt_bq(dwegt_s) = dwegt_bq_v;
%     dwegt_folder = strcat(imagename(1:end-4),'_dwegt');
%     if isfolder(dwegt_folder) == 0
%         mkdir(dwegt_folder);
% %     else
% %         disp('文件夹已存在');
%     end
%     imwrite(dwegt_bq,strcat(dwegt_folder,'\',imagename(1:end-4),'-dwegt-',num2str(dwegt_v),'.png'))
% 
% end
% 
% %%% 保存EGT算法的每个细胞结果
% egt_r = bwlabel(res_EGT);
% for egt_v = 1:max(egt_r,[],'all')
%     egt_bq = double(zeros(size(egt_r)));
%     [egt_bq_r,egt_bq_c,egt_bq_v] = find(egt_r == egt_v);
%     egt_s = find(egt_r == egt_v);
%     egt_bq(egt_s) = egt_bq_v;
%     egt_folder = strcat(imagename(1:end-4),'_egt');
%     if isfolder(egt_folder) == 0
%         mkdir(egt_folder);
% %     else
% %         disp('文件夹已存在');
% 
%     end
%     imwrite(egt_bq,strcat(egt_folder,'\',imagename(1:end-4),'-egt-',num2str(egt_v),'.png'))
% 
% end
% 
% %%% 保存Jaccard算法的每个细胞结果
% jaccard_r = bwlabel(jaccard_pic);
% for jaccard_v = 1:max(jaccard_r,[],'all')
%     jaccard_bq = double(zeros(size(jaccard_r)));
%     [jaccard_bq_r,jaccard_bq_c,jaccard_bq_v] = find(jaccard_r == jaccard_v);
%     jaccard_s = find(jaccard_r == jaccard_v);
%     jaccard_bq(jaccard_s) = jaccard_bq_v;
%     jaccard_folder = strcat(imagename(1:end-4),'_jaccard');
%     if isfolder(jaccard_folder) == 0
%         mkdir(jaccard_folder);
% %     else 
% %         disp('文件夹已存在');
%     end
%     imwrite(jaccard_bq,strcat(jaccard_folder,'\',imagename(1:end-4),'-jaccard-',num2str(jaccard_v),'.png'))
% 
% end
% 
% %%% 保存CV算法的每个细胞结果
% cv_r = bwlabel(res_chanvese);
% for cv_v = 1:max(cv_r,[],'all')
%     cv_bq = double(zeros(size(cv_r)));
%     [cv_bq_r,cv_bq_c,cv_bq_v] = find(cv_r == cv_v);
%     cv_s = find(cv_r == cv_v);
%     cv_bq(cv_s) = cv_bq_v;
%     cv_folder = strcat(imagename(1:end-4),'_cv');
%     if isfolder(cv_folder) == 0
%         mkdir(cv_folder);
% %     else
% %         disp('文件夹已存在');
%     end
%     imwrite(cv_bq,strcat(cv_folder,'\',imagename(1:end-4),'-cv-',num2str(cv_v),'.png'))
% 
% end

