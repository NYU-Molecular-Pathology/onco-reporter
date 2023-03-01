import csv
from oncoreport.models import PMKBdb,Tumortypes

## Import backend PMKB data into PMKBdb table ##
with open('oncomine_backend_v4.csv') as csvfile1:
    reader = csv.DictReader(csvfile1)
    for row in reader:
        p = PMKBdb(tumor_type=row['tumor_type'], gene=row['gene'], variant=row['variant'], tier=row['tier'],interpretations=row['interpretations'], citations=row['citations'])
        p.save()
print("Successfully Imported PMKB data!!")

## Import backend Tumortype data into Tumortypes table ##
with open('Site_Diagnosis_TumorType.csv') as csvfile2:
    reader = csv.DictReader(csvfile2)
    for row in reader:
        t = Tumortypes(site_of_origin=row['Site_of_Origin'], diagnosis=row['Diagnosis'], tumor_type=row['Tumor_type_forComment'])
        t.save()
print("Successfully Imported Tumortypes data!!")