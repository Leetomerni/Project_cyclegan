from django.urls import path
from . import views

urlpatterns = [
    path('', views.index, name='index'),
    path('download/<str:filename>/', views.download_result, name='download_result'),
]
