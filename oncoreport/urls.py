from django.urls import path
from . import views

urlpatterns = [
    path('', views.home, name ='home'),
    path('upload_input_report', views.upload_input_report, name='upload_input_report'),
    path('import_onco_report', views.import_onco_report, name='import_onco_report'),
    path('search_sample_tumortype', views.search_sample_tumortype, name='search_sample_tumortype'),
    path('index/<runid>', views.index, name='index'),
    path('index/<runid>/<sample>/<tumor_type>', views.index_runid_sample, name='index_runid_sample'),
    path('edit/<int:id>', views.edit),
    path('update/<int:id>', views.update),
    path('preview_report', views.preview_report, name='preview_report'),
]