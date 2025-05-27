from .base_options import BaseOptions


class TestOptions(BaseOptions):
    """This class includes test options.

    It also includes shared options defined in BaseOptions.
    """


    def parse(self):
        """Parse our options and set up gpu device without command line arguments."""
        # 模拟 argparse.Namespace 对象
        opt = type('Namespace', (), {})()
        # 手动设置测试所需参数
        opt.aspect_ratio = 1.0
        opt.batch_size = 1
        opt.checkpoints_dir = './checkpoints'
        opt.crop_size = 256
        opt.dataroot = './datasets/Microtubules-1/testA'
        opt.dataset_mode = 'single'
        opt.direction = 'AtoB'
        opt.display_winsize = 256
        opt.epoch = 'latest'
        opt.eval = False
        opt.gpu_ids = '0'
        opt.init_gain = 0.02
        opt.init_type = 'normal'
        opt.input_nc = 1
        opt.output_nc = 1
        opt.isTrain = False
        opt.load_iter = 0
        opt.load_size = 256
        opt.max_dataset_size = float("inf")
        opt.model = 'test'
        opt.n_layers_D = 3
        opt.name='Microtubules-1_pretrained'
        opt.ngf = 64
        opt.ndf = 64
        opt.netD = 'basic'
        opt.netG = 'resnet_9blocks'
        opt.norm = 'instance'
        opt.no_dropout = True
        opt.no_flip = False
        opt.num_test = 10
        opt.num_threads = 8
        opt.phase = 'test'
        opt.preprocess = 'crop'
        opt.results_dir = './results/'
        opt.serial_batches= False
        opt.verbose = False
        opt.suffix = ''
        opt.use_wandb = False
        opt.wandb_project_name = 'CycleGAN-and-pix2pix'

        # GPU 设置

        str_ids = opt.gpu_ids.split(',')
        opt.gpu_ids = []
        for str_id in str_ids:
            id = int(str_id)
            if id >= 0:
                opt.gpu_ids.append(id)

        if len(opt.gpu_ids) > 0:
            import torch
            torch.cuda.set_device(opt.gpu_ids[0])

        self.opt = opt
        return self.opt