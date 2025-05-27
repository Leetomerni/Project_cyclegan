from django.shortcuts import render
from django.http import FileResponse, Http404
from .forms import ImageUploadForm
import os
import sys
from django.conf import settings
from PIL import Image
# Create your views here.
sys.path.append(r"D:\Django_pro")#改成自己的项目路径，不要有中文
from forDjangotest import process_image  # 导入你写的函数
def index(request):
    if request.method == 'POST':
        form = ImageUploadForm(request.POST, request.FILES)
        selected_option = request.POST.get('option')
        if selected_option == 'option1':
            construction_model = 1
        else :
            construction_model = 2

        if form.is_valid():
            image = form.cleaned_data['image']
            input_path = os.path.join(settings.MEDIA_ROOT, image.name)
            output_path = os.path.join(settings.MEDIA_ROOT, 'processed_' + image.name)

            with open(input_path, 'wb+') as f:
                for chunk in image.chunks():
                    f.write(chunk)

            process_image(input_path, output_path, construction_model)

            image_name_base = os.path.splitext(image.name)[0]
            image_name_base = image_name_base + '.png'

            image_in = Image.open(input_path)
            image_in_png = os.path.join(settings.MEDIA_ROOT ,image_name_base)
            image_in.save(image_in_png)
            image_out = Image.open(output_path)
            image_out_png = os.path.join(settings.MEDIA_ROOT ,'processed_' + image_name_base)
            image_out.save(image_out_png)

            original_url = os.path.join(settings.MEDIA_URL, image_name_base)
            processed_url = os.path.join(settings.MEDIA_URL, 'processed_' + image_name_base)
            return render(request, 'index.html', {
                'form': form,
                'original_image': original_url,
                'processed_image':processed_url,
                'processed_filename': os.path.basename(output_path)  # 用于下载
            })
    else:
        form = ImageUploadForm()
    return render(request, 'index.html', {'form': form})
def download_result(request, filename):
    file_path = os.path.join(settings.MEDIA_ROOT, filename)
    if os.path.exists(file_path):
        return FileResponse(open(file_path, 'rb'), as_attachment=True, filename=filename)
    else:
        raise Http404("文件未找到")