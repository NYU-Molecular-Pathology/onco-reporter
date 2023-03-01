from django.contrib import admin
from .models import Oncoreport, Tumortypes, PMKBdb, Allaberrations

# Register your models here.
admin.site.register(Tumortypes),
admin.site.register(Oncoreport),
admin.site.register(PMKBdb),
admin.site.register(Allaberrations),
