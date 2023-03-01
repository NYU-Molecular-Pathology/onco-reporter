from django.db import models
from django.conf import settings

# Create your models here.
class Tumortypes(models.Model):
    site_of_origin = models.CharField(blank=False, max_length=255)
    diagnosis = models.CharField(blank=False, max_length=255)
    tumor_type = models.CharField(blank=False, default="NA", max_length=255)

    def __str__(self):
        return self.site_of_origin


class PMKBdb(models.Model):
    tumor_type = models.TextField(blank=False)
    gene = models.CharField(blank=False, max_length=255) 
    variant = models.CharField(blank=False, max_length=255)
    tier = models.CharField(blank=False,max_length=25)
    interpretations = models.TextField(blank=False)
    citations = models.TextField(blank=False)
    
    def __str__(self):
        return self.tumor_type


class Oncoreport(models.Model):
    runid = models.CharField(blank=False,default='',max_length=255)
    sample = models.CharField(blank=False, max_length=255)
    locus = models.CharField(blank=False, max_length=255)
    genes = models.CharField(blank=False, max_length=255)
    variant_type = models.CharField(blank=False, default='NA',max_length=255)
    exon = models.IntegerField(default=0)
    transcript = models.CharField(default='NA',max_length=255)
    coding = models.CharField(blank=False, max_length=255)
    variant_effect = models.CharField(default='NA',max_length=255)
    frequency = models.DecimalField(max_digits=6,decimal_places=2,default=0)
    amino_acid_change = models.CharField(default='NA',max_length=255)
    aa = models.CharField(default='NA',max_length=255)
    copy_number = models.DecimalField(max_digits=6,decimal_places=2,default=0)
    coverage = models.IntegerField(default=0)
    class Meta:
        db_table="oncoreport"

class Allaberrations(models.Model):
    runid = models.CharField(blank=False,default='',max_length=255)
    sample = models.CharField(blank=False, max_length=255)
    gene = models.CharField(blank=False, max_length=255)
    variant_type = models.CharField(blank=False, default='NA',max_length=255)
    transcript = models.CharField(default='NA',max_length=255)
    variants = models.CharField(default='NA',max_length=255)
    vaf = models.DecimalField(blank=False,max_digits=6,decimal_places=2,default=0)
    aa = models.CharField(default='NA',max_length=255)
    exon = models.IntegerField(blank=False,default=0)
    copy_number = models.DecimalField(max_digits=6,decimal_places=2,default=0)
    depth = models.IntegerField(blank=False,default=0)
    tier = models.CharField(default='1/2',max_length=25)
    locus = models.CharField(blank=False,default='NA',max_length=255) 
    variant_pmkb = models.TextField(default='NA')
    pmkb_tumortype = models.TextField(default='NA')
    interpretations = models.TextField(default='NA')
    citations = models.TextField(default='NA')
    status = models.TextField(default='NA')



