import cv2
from skimage import util
import torch
import random
import numpy as np
from tifffile import tifffile

from options.test_options import TestOptions
from models import create_model

device = torch.device("cuda" if torch.cuda.is_available() else "cpu")

def gaussian_weighted_crop(image, model):
    """
    采用高斯加权融合，减少拼接痕迹。
    """
    crop_size=256
    stride=44

    _,h,w = image.shape

    # **2. 创建高斯加权矩阵**
    x = np.linspace(-1, 1, crop_size)#生成一个从-1到1均匀分布的数组，用于后续创建高斯加权矩阵
    y = np.linspace(-1, 1, crop_size)
    xx, yy = np.meshgrid(x, y)#将这两个一维数组扩展为二维网格坐标，xx 表示 x 方向的坐标矩阵，yy 表示 y 方向的坐标矩阵。
    gaussian_kernel = np.exp(- (xx**2 + yy**2)/0.5)  # 2D 高斯
    gaussian_kernel = torch.tensor(gaussian_kernel, dtype=torch.float32)
    gaussian_kernel /= gaussian_kernel.max()  # 归一化到 0~1
    gaussian_kernel = gaussian_kernel.unsqueeze(0).to(device)  # (1, 256, 256)

    # **3. 进行滑动窗口预测**
    output_map = torch.zeros(1, h, w).to(device)
    count_map = torch.zeros(1, h, w).to(device)

    for i in range(0,18):
        for j in range(0,18):
            crop = image[:,i*stride:i*stride+crop_size,j*stride:j*stride+crop_size]
            crop = crop.unsqueeze(0)
            with torch.no_grad():
                pred = model.netG(crop)
                pred = pred.squeeze(0)

            weighted_pred = pred * gaussian_kernel
            output_map[:,i*stride:i*stride+crop_size,j*stride:j*stride+crop_size] += weighted_pred
            count_map[:,i*stride:i*stride+crop_size,j*stride:j*stride+crop_size] += gaussian_kernel

    # **4. 计算加权平均**
    output_map /= count_map
    return output_map



# **加载模型**
opt = TestOptions().parse()  # get test options
opt.num_threads = 0   # test code only supports num_threads = 0
opt.batch_size = 1    # test code only supports batch_size = 1
opt.serial_batches = True  # disable data shuffling; comment this line if results on randomly chosen images are needed.
opt.no_flip = True    # no flip; comment this line if results on flipped images are needed.
opt.display_id = -1   # no visdom display; the test code saves the results to a HTML file.
model = create_model(opt)      # create a model given opt.model and other options
model.setup(opt)               # regular setup: load and print networks; create schedulers
if opt.eval:
    model.eval()

# **加载图像**
import os
def find_dic(directory):
    return [f for f in os.listdir(directory) if os.path.isfile(os.path.join(directory, f))]

directory_path = './datasets/Microtubules-1/testA/'
# directory_path = './datasets/ER-1/testA/'
file_names = find_dic(directory_path)

for i in range(len(file_names)):
    image_path = directory_path + file_names[i]
    print('正在处理',image_path)
    image = cv2.imread(image_path, cv2.IMREAD_UNCHANGED)
    image = torch.tensor(image,dtype = torch.float32)[None]/65535.0

    #获取结果
    result = gaussian_weighted_crop(image, model)

    image = image.cpu().float().numpy()  # 转换成numpy
    image = ((image + 1) / 2.0 * 65535.0).astype(np.uint16)  # 转换回0-65535

    #保存原始图像
    save_path_true = f'./Myresult/cell{i + 46}_true.tif'
    # 去除多余的维度
    image_2d = np.squeeze(image, axis=0)
    # 保存为16位TIFF
    tifffile.imwrite(save_path_true, image_2d, dtype=np.uint16)
    print('已保存', save_path_true)

    result = result.cpu().float().numpy()#转换成numpy
    result =((result+1)/2.0 * 65535.0).astype(np.uint16)#转换回0-65535

    #保存预测图像
    save_path_fake = f'./Myresult/cell{i+46}_fake.tif'
    # 去除多余的维度
    result_2d = np.squeeze(result, axis=0)
    # 保存为16位TIFF
    tifffile.imwrite(save_path_fake, result_2d, dtype=np.uint16)
    print('已保存', save_path_fake)
