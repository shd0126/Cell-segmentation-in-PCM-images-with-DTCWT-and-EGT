function outputArg1 = Lenhanced(Yl,strel1)
%UNTITLED 原图+顶帽变换-底帽变换
%   Yl输入图像，strel1结构元素
Yl_b = imbothat(Yl,strel1); % 对低频部分进行底帽变换，因为原图中是需要增强细胞内部的暗部分：Iclose-I
Yl_t = imtophat(Yl,strel1); % 增强图像中亮的部分：I-Iopen
outputArg1 = Yl + Yl_t - Yl_b ; % 原图 + 顶帽变换 - 底帽变换：result = I + (I - Io) - (Ic - I) = 3I - Io - Ic
end