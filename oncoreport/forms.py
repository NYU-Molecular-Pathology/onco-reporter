from django import forms
from .models import Allaberrations

class AllaberrationsForm(forms.ModelForm):
    class Meta:
        model = Allaberrations
        fields = ['gene','variants','tier','variant_type','vaf','depth','transcript','locus','exon','copy_number','variant_pmkb','pmkb_tumortype']

        widgets = {
            'gene': forms.TextInput(attrs={ 'class': 'form-control'}),
            'variants': forms.TextInput(attrs={ 'class': 'form-control'}),
            'tier': forms.TextInput(attrs={ 'class': 'form-control'}),
            'variant_type': forms.TextInput(attrs={ 'class': 'form-control'}),
            'vaf': forms.TextInput(attrs={ 'class': 'form-control'}),
            'depth': forms.TextInput(attrs={ 'class': 'form-control'}),
            'transcript': forms.TextInput(attrs={ 'class': 'form-control'}),
            'locus': forms.TextInput(attrs={ 'class': 'form-control'}),
            'exon': forms.TextInput(attrs={ 'class': 'form-control'}),
            'copy_number': forms.TextInput(attrs={ 'class': 'form-control'}),
            'variant_pmkb': forms.TextInput(attrs={ 'class': 'form-control'}),
            'pmkb_tumortype': forms.TextInput(attrs={ 'class': 'form-control'}),
        } 