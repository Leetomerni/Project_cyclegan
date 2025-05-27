import os

import cv2
import numpy as np
import torch
from matplotlib import pyplot as plt

from data.base_dataset import BaseDataset, get_transform
from data.image_folder import make_dataset
from PIL import Image
import random
from skimage import util
def add_noise(image):

    # 归一化到 [0, 1] 区间
    image_norm = image / 65535.0
    gaussian_var = random.uniform(0.001, 0.01)

    # 添加高斯噪声
    gaussian_noisy = util.random_noise(image_norm, mode='gaussian', var=gaussian_var)

    # 还原为16位图像
    noisy_image = np.clip(gaussian_noisy * 65535, 0, 65535).astype(np.uint16)
    return noisy_image
class UnalignedDataset(BaseDataset):
    """
    This dataset class can load unaligned/unpaired datasets.

    It requires two directories to host training images from domain A '/path/to/data/trainA'
    and from domain B '/path/to/data/trainB' respectively.
    You can train the model with the dataset flag '--dataroot /path/to/data'.
    Similarly, you need to prepare two directories:
    '/path/to/data/testA' and '/path/to/data/testB' during test time.
    """

    def __init__(self, opt):
        """Initialize this dataset class.

        Parameters:
            opt (Option class) -- stores all the experiment flags; needs to be a subclass of BaseOptions
        """
        BaseDataset.__init__(self, opt)
        self.dir_A = os.path.join(opt.dataroot, opt.phase + 'A')  # create a path '/path/to/data/trainA'
        self.dir_B = os.path.join(opt.dataroot, opt.phase + 'B')  # create a path '/path/to/data/trainB'
        #dir_A = ./datasets/F-actin/trainA, dir_B = ./datasets/F-actin/trainB
        self.A_paths = sorted(make_dataset(self.dir_A, opt.max_dataset_size))   # load images from '/path/to/data/trainA'
        self.B_paths = sorted(make_dataset(self.dir_B, opt.max_dataset_size))    # load images from '/path/to/data/trainB'
        self.A_size = len(self.A_paths)  # get the size of dataset A
        self.B_size = len(self.B_paths)  # get the size of dataset B
        #print('dataset A: %d, dataset B: %d' % (self.A_size, self.B_size))
        btoA = self.opt.direction == 'BtoA'#direction 默认是 AtoB
        #print('direction:',btoA)--false
        input_nc = self.opt.output_nc if btoA else self.opt.input_nc       # get the number of channels of input image
        output_nc = self.opt.input_nc if btoA else self.opt.output_nc      # get the number of channels of output image
        #print('通道',input_nc,output_nc)
        self.transform_A = get_transform(self.opt, grayscale=(input_nc == 1))
        self.transform_B = get_transform(self.opt, grayscale=(output_nc == 1))
        #数据增强次数
        self.n_aug = self.opt.n_aug
        #是否添加噪声
        self.noise = self.opt.noise

    def __getitem__(self, index):
        """Return a data point and its metadata information.

        Parameters:
            index (int)      -- a random integer for data indexing

        Returns a dictionary that contains A, B, A_paths and B_paths
            A (tensor)       -- an image in the input domain
            B (tensor)       -- its corresponding image in the target domain
            A_paths (str)    -- image paths
            B_paths (str)    -- image paths
        """
        original_idx = index//self.n_aug
        A_path = self.A_paths[original_idx % self.A_size]  # make sure index is within then range
        if self.opt.serial_batches:   # make sure index is within then range
            index_B = original_idx % self.B_size
        else:   # randomize the index for domain B to avoid fixed pairs.
            index_B = random.randint(0, self.B_size - 1)
        B_path = self.B_paths[index_B]
        #从文件中读取图像
        A_img = cv2.imread(A_path, cv2.IMREAD_UNCHANGED)
        B_img = cv2.imread(B_path, cv2.IMREAD_UNCHANGED)
        #对低分辨率图像添加噪声
        if self.noise == True:
            A_img = add_noise(A_img)
        #归一化
        A_img = torch.tensor(A_img,dtype=torch.float32)[None]/65535.0
        B_img = torch.tensor(B_img,dtype=torch.float32)[None]/65535.0
        # 应用transforms
        A = self.transform_A(A_img)
        B = self.transform_B(B_img)

        return {'A': A, 'B': B, 'A_paths': A_path, 'B_paths': B_path}

    def __len__(self):
        """Return the total number of images in the dataset.

        As we have two datasets with potentially different number of images,
        we take a maximum of
        """
        return max(self.A_size, self.B_size)*self.n_aug
