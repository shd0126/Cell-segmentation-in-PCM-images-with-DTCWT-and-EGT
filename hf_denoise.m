function output = denoise(input)
% 适应的阈值去噪函数
sigma_n2 = median(abs(input),'all')/0.6745; % 求噪声的方差
sigma_f2 = var(input,0,'all'); % 每个尺度每个子带虚部f的方差
sigma_s = sqrt(max(sigma_f2 - sigma_n2,0)); % 每个尺度每个子带虚部不含噪声的标准差
output = sign(input).*max(abs(input)-sqrt(2)*sigma_n2/sigma_s,0);  % 求最后的输出结果
end

