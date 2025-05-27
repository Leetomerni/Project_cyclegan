# Project_cyclegan

这里包含了四个项目文件夹和一个环境配置文件，PCA-SIM是Matlab项目，其余三个是Python项目，环境配置是conda环境。

### 配置conda环境

见environment.yml，这个环境是三个Python项目公用的。

### 数据集

使用BioSR数据集，读者需自行前往https://figshare.com/articles/dataset/BioSR/13264793下载。

### 关于PCA-SIM

这个文件夹用于处理低分辨率图像，将原始数据的多帧低分辨率图像处理为一张低分辨率和一张高分辨率图像。

具体使用方法见PCA-SIM/Guide，本项目只需要使用test部分的代码，即PCA-SIM/test-1。按照指导配置好Matlab工具包之后，只需要运行test-1文件即可，需要自己替换数据集的路径。

原始论文链接在https://elight.springeropen.com/articles/10.1186/s43593-022-00035-x，论文也提供了更全面的开源代码包，本项目只上传了用得到的部分。

### 关于CycleGAN

这个文件夹是原始CycleGAN的代码，源自Github项目https://github.com/junyanz/pytorch-CycleGAN-and-pix2pix，本项目对源代码进行了适应性修改。具体使用方法如下：

##### 训练

打开train.py，在参数中配置：

```
--dataroot
./datasets/Microtubules
--name
Microtubules
--model
cycle_gan
--n_epochs_decay
50
--n_epochs
50
--n_aug
30
--noise
False
```

dataroot是训练数据集的路径；name是训练产生的文件夹的存储名字；model是表示CycleGAN模型；n_epochs_decay表示学习率不变的批次数量；n_epochs表示学习率线性下降的批次数量；noise表示是否添加噪声进行训练。

训练产生的所有文件在CycleGAN/checkpoints/“数据集名字”

##### 测试

打开Mytest.py，在参数中配置：

```
--dataroot
./datasets/Microtubules/testA--
--name
Microtubules_pretrained
--num_test
10
--model
test
--no_dropout
```

dataroot是测试数据集的路径，只需要低分辨率图像；name是训练产生的pth文件的名字，要先从CycleGAN/checkpoints/“数据集名字”下取出”lasted_netG_A“文件，把这个文件复制到CycleGAN/checkpoints/“数据集名字_pretrained”下，并且改名为“lasted_netG”；num_test表示用到数据集中多少个图像；model是表示测试模型。

测试产生的文件会存储在自己设置的路径里，需要读者在Mytest.py文件中自行修改。

### 关于MyCycleGAN

这是作者提出的基于物理过程优化的CycleGAN模型，具体使用方法和CycleGAN完全相同。

### 关于Django_pro

这是基于Django的网页系统，该项目依托作者给出的conda环境，如果是读者自己配制的虚拟环境，需要在虚拟环境里面安装Django包。

运行方法为：

进入Django_pro/SRpro，打开cmd，并且执行以下命令：

```
conda activate "虚拟环境名字"
python manage.py runserver
```

之后点击cmd给出的链接，即可访问。

### 致谢

感谢项目的贡献者开源