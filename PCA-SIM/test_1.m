%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Main file for Experiment for SIM.
% Version 2.0 -
% Related Reference:
% last modified on 06/17/2022
% by Jiaming Qian, Yu Cao and Chao Zuo (zuochao@njust.edu.cn,jiaming_qian@njust.edu.cn)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear;
clc;
close all

%% Set internal parameter
pixel_size=6.26;                                                             % Pixel size
mag=100;                                                                    % Magnification
NA=1.45;                                                                     % NA,别动了这个值就刚好
lamda=0.488;                                                                % Fluorescent wavelength      
sub_optimization=0;


%% Read raw image

%for train
% output_folder_WF = 'D:\毕业设计\doc from article 1\Supplementary code and data\Test\trainA\';  % 指定保存文件夹
% output_folder_WI = 'D:\毕业设计\doc from article 1\Supplementary code and data\Test\trainB\';  % 指定保存文件夹


%for test
output_folder_WF = 'D:\毕业设计\doc from article 1\Supplementary code and data\Test\testA\';  % 指定保存文件夹
output_folder_WI = 'D:\毕业设计\doc from article 1\Supplementary code and data\Test\testB\';  % 指定保存文件夹



for i = 46:55

    output_filename_WF = sprintf('cell%02d.tif',i);       % 输出文件名
    output_path_WF = fullfile(output_folder_WF, output_filename_WF);
    
    output_filename_WI = sprintf('cell%02d.tif',i);       % 输出文件名
    output_path_WI = fullfile(output_folder_WI, output_filename_WI);
    
    input_folder = 'D:\毕业设计\dataset1\Microtubules\';
    input_folder = fullfile(input_folder,sprintf('Cell_%03d',i));
    input_file = 'RawSIMData_level_06.mrc';
    input_path = fullfile(input_folder,input_file);
    [header,Iraw] = read_mrc(input_path);
    
    % 重塑为三维数组
    Iraw = reshape(Iraw, [header(1), header(2), header(3)]);
    
    % 处理复数数据（若mode=4）
    if header(4) == 4
        Iraw = complex(Iraw(1:2:end), Iraw(2:2:end));
        Iraw = reshape(Iraw, [header(1), header(2), header(3)]);
    end
    
    % 调整维度顺序
    Iraw = permute(Iraw, [2, 1, 3]);
    
    % % 显示九张切片
    % figure;
    % num_slices = double(header(3)); % Z方向总层数
    % selected_slices = round(linspace(1, num_slices, 9)); % 均匀选择9个索引
    % 
    % for i = 1:9
    %     subplot(3, 3, i);
    %     slice_idx = selected_slices(i);
    %     imshow(Iraw(:, :, slice_idx), []);
    %     title(sprintf('Z = %d', slice_idx));
    %     axis off; % 隐藏坐标轴
    % 
    % end
    % 
    % % 添加总标题
    % sgtitle('均匀分布的Z层切片');
    
    Iraw = double(Iraw);
    param=parameter_set(Iraw,pixel_size,NA,lamda,mag);                          % parameters setting
    IIraw=edgefilter(Iraw);                                                     % edgefilter
    IIraw=RLdeconv(Iraw,param.psf,1);                                           % RLdeconv
    WFimage=edgefilter(WF_double(IIraw));                                       % double size Rawimage
    % figure;imshow(WFimage,[]);title('WF');
    max_value_WF = max(WFimage(:));
    disp(['最大像素值: ', num2str(max_value_WF)]);
    % 旋转 180 度（顺时针）
    WFimage_out = imrotate(WFimage, -180);
    
    % 水平翻转
    WFimage_out = flip(WFimage_out, 2);
    WFimage_out = uint16((WFimage_out / max_value_WF) * 65535);
    imwrite(WFimage_out, output_path_WF);
    %% PCA
    sub_optimization=1;
    Filter_size=11;%default= 11
    Mask_size=3;
    par=Parameter_estimation_PCA(IIraw,param,Filter_size,Mask_size);            % parameter estimation PCA
    
    %% COR
    % sub_optimization=1;
    % par=Parameter_estimation_COR(IIraw,param);                                % parameter estimation COR
    
    %% POP
    % par=Parameter_estimation_POP(IIraw,param);
    
    %% ACR
    % par=Parameter_estimation_ACR(IIraw,param);
    
    %% IRT
    % par=Parameter_estimation_IRT(IIraw,param);
    
    %% Wiener reconstruction
    Wiener=edgefilter(Wiener_recon(IIraw,param,par,sub_optimization));        % Wiener_reconstruction
    % figure;imshow(Wiener,[]);title('Wiener');
    max_value_WI = max(Wiener(:));
    disp(['最大像素值: ', num2str(max_value_WI)]);
    % 旋转 180 度（顺时针）
    Wiener_out = imrotate(Wiener, -180);
    
    % 水平翻转
    Wiener_out = flip(Wiener_out, 2);
    Wiener_out = uint16((Wiener_out / max_value_WI) * 65535);
    imwrite(Wiener_out, output_path_WI);
end
%% HiFi reconstruction
% HiFi=edgefilter(Hifi_recon(IIraw,param,par,sub_optimization));
% figure;imshow(HiFi,[]),colormap hot;title('HiFi')                           % HiFi_reconstruction
% figure;imshow(log(abs(fftshift(fft2(HiFi)))),[])

%% TV reconstruction
% TV=edgefilter(TV_recon(IIraw,param,par,sub_optimization));
% figure;imshow(TV,[]),colormap hot;title('TV')                             % TV_reconstruction
