from django.shortcuts import render, redirect
from django.urls import reverse
import pandas as pd
import numpy as np
import json
from . import process_input
from .forms import AllaberrationsForm
from .models import Tumortypes, Oncoreport, PMKBdb, Allaberrations

# Create your views here.
def home(request):
    return render(request, 'oncoreport/home.html', {} )

def upload_input_report(request):
    return render(request, 'oncoreport/upload_input_report.html', {})

def import_onco_report(request):
    csv_file = request.FILES["csv_file"]
    csv_file_name = request.FILES['csv_file'].name
    name = csv_file_name.split(".csv")
    filename = name[0]
    input_df = pd.read_csv(csv_file)
    input_df['Runid'] = filename

    if Oncoreport.objects.filter(runid=filename).exists():
        tumortypes_db = pd.DataFrame(list(Tumortypes.objects.all().values()))
        key_values = zip(tumortypes_db.site_of_origin,tumortypes_db.diagnosis)
        site_diagnosis_dict = dict()
        ## key = site_diagnosis[0](site of origin), value = site_diagnosis[1](list of diagnosis tumor types)
        for site_diagnosis in key_values:
            if site_diagnosis[0] in site_diagnosis_dict:
                site_diagnosis_dict[site_diagnosis[0]].append(site_diagnosis[1])
            else:
                site_diagnosis_dict[site_diagnosis[0]] = [site_diagnosis[1]]

        all_data = pd.DataFrame(list(Oncoreport.objects.filter(runid=filename).values()))
        samples = all_data['sample'].unique()

        return render(request, 'oncoreport/import_onco_report.html',{'samples':samples, 'site_diagnosis_dict':site_diagnosis_dict})
    else:
        ## Strip of all samples except TMs ##
        input_df_TMonly = input_df.loc[(input_df['Sample'].str.startswith("TM",na=False))]
        input_df_TMonly_setindex = input_df_TMonly.reset_index(drop=True)
        input_df_TMonly_setindex['Genes'] = input_df_TMonly_setindex['Genes'].fillna("NEGATIVE")
        oncomine_report = input_df_TMonly_setindex

        ## Load to database for input oncomine report ##
        row_iter_all = oncomine_report.iterrows()

        objs_all = [
            Oncoreport(
                    runid = row['Runid'],
                    sample = row['Sample'],
                    locus = row['Locus'],
                    genes = row['Genes'],
                    variant_type = row['Type'],
                    ## convert NA in exon,freq,copynumber columns of report to 0 ##
                    exon = row['Exon'],
                    transcript = row['Transcript'],
                    coding = row['Coding'],
                    variant_effect = row['Variant.Effect'],
                    frequency = row['Frequency'],
                    amino_acid_change = row['Amino.Acid.Change'],
                    aa = row['AA'],
                    copy_number = row['Copy.Number'],
                    coverage = row['Coverage']
                )
                for index, row in row_iter_all
                ]

        Oncoreport.objects.bulk_create(objs_all)

        tumortypes_db = pd.DataFrame(list(Tumortypes.objects.all().values()))
        key_values = zip(tumortypes_db.site_of_origin,tumortypes_db.diagnosis)
        site_diagnosis_dict = dict()
        ## key = site_diagnosis[0](site of origin), value = site_diagnosis[1](list of diagnosis tumor types)
        for site_diagnosis in key_values:
            if site_diagnosis[0] in site_diagnosis_dict:
                site_diagnosis_dict[site_diagnosis[0]].append(site_diagnosis[1])
            else:
                site_diagnosis_dict[site_diagnosis[0]] = [site_diagnosis[1]]

        samples = oncomine_report['Sample'].unique()

        return render(request, 'oncoreport/import_onco_report.html',{'samples':samples, 'site_diagnosis_dict':site_diagnosis_dict})

def search_sample_tumortype(request):
    sample = request.POST['Sample']
    site_of_origin = request.POST['siteDropdown']
    diagnosis = request.POST['diagnosisDropdown']

    oncomine_report = pd.DataFrame(list(Oncoreport.objects.all().values()))
    PMKB = pd.DataFrame(list(PMKBdb.objects.all().values()))
    tumortype_db = pd.DataFrame(list(Tumortypes.objects.filter(site_of_origin=site_of_origin,diagnosis=diagnosis).values()))
        # to get single value instead of whole object #
    tumor_type = tumortype_db['tumor_type'].values[0]

    if Allaberrations.objects.filter(sample=sample).exists():
        abberations = pd.DataFrame(list(Allaberrations.objects.filter(sample=sample).values()))
        runid_abberations = abberations.runid.unique()
        runid = runid_abberations[0]
        json_records = abberations.to_json(orient ='records') 
        data_df = []
        data_df = json.loads(json_records)
        return render(request, 'oncoreport/search_sample_tumortype.html', {'data_df':data_df,'tumor_sample':sample,'tumor_type':tumor_type,'runid':runid})
    else:
        input_all_abberations = process_input.main(oncomine_report,PMKB,sample,tumor_type)
        row_iter_all = input_all_abberations.iterrows()

        aberrations_all = [
            Allaberrations(
                    runid = row['runid'],
                    sample = row['sample'],
                    gene = row['Gene'],
                    variant_type = row['variant_type'],
                    transcript = row['transcript'],
                    variants = row['Variants'],
                    vaf = row['frequency'],
                    aa = row['aa'],
                    exon = row['exon'],
                    copy_number = row['copy_number'],
                    depth = row['coverage'],
                    tier = row['tier'],
                    locus = row['locus'],
                    variant_pmkb = row['variant'],
                    pmkb_tumortype = row['tumor_type'],
                    interpretations = row['interpretations'],
                    citations = row['citations'],
                )
                for index, row in row_iter_all
                ]
    
        Allaberrations.objects.bulk_create(aberrations_all)
    
        abberations = pd.DataFrame(list(Allaberrations.objects.filter(sample=sample).values()))
        runid_abberations = abberations.runid.unique()
        runid = runid_abberations[0]
        json_records = abberations.to_json(orient ='records') 
        data_df = []
        data_df = json.loads(json_records)
        return render(request, 'oncoreport/search_sample_tumortype.html', {'data_df':data_df,'tumor_sample':sample,'tumor_type':tumor_type,'runid':runid})

def index(request, runid):
    tumortypes_db = pd.DataFrame(list(Tumortypes.objects.all().values()))
    key_values = zip(tumortypes_db.site_of_origin,tumortypes_db.diagnosis)
    site_diagnosis_dict = dict()
    ## key = site_diagnosis[0](site of origin), value = site_diagnosis[1](list of diagnosis tumor types)
    for site_diagnosis in key_values:
        if site_diagnosis[0] in site_diagnosis_dict:
            site_diagnosis_dict[site_diagnosis[0]].append(site_diagnosis[1])
        else:
            site_diagnosis_dict[site_diagnosis[0]] = [site_diagnosis[1]]

        all_data = pd.DataFrame(list(Oncoreport.objects.filter(runid=runid).values()))
        samples = all_data['sample'].unique()
    return render(request, 'oncoreport/import_onco_report.html',{'samples':samples, 'site_diagnosis_dict':site_diagnosis_dict})

def index_runid_sample(request, runid, sample, tumor_type):
    allaberrations = Allaberrations.objects.filter(runid=runid,sample=sample)
    selected_sample = sample
    tumor_type = tumor_type
    runid = runid

    return render(request,"oncoreport/show_persample.html",{'allaberrations':allaberrations,'tumor_sample':selected_sample,'tumor_type':tumor_type, 'runid':runid})

def edit(request, id):
    allaberrations = Allaberrations.objects.get(id=id)
    return render(request, 'oncoreport/edit.html', {'allaberrations': allaberrations})

def update(request, id):
    allaberrations_id = Allaberrations.objects.get(id=id)
    runid = allaberrations_id.runid
    sample = allaberrations_id.sample
    tumor_type = allaberrations_id.pmkb_tumortype

    allaberrations = Allaberrations.objects.filter(id=allaberrations_id.id,runid=runid).first()
    form = AllaberrationsForm(request.POST, instance = allaberrations)
    if form.is_valid():
        form.save()
        return redirect('index_runid_sample',runid=runid,sample=sample,tumor_type=tumor_type)
    return render(request, 'oncoreport/edit.html', {'allaberrations': allaberrations} ) 

def preview_report(request):
    selected_variants = request.POST.getlist('selected_values')
    sample = request.POST['tumorsample']
    tumor_sample = sample.split("-B")[0]
    tumor_type = request.POST['tumortype']

    if not selected_variants:
        ## Get the run id ##
        df = pd.DataFrame(list(Allaberrations.objects.filter(sample=sample).values()))
        ## Get the runid to redirect to index page with runid ##
        variants_runid = df.runid.unique()
        runid = variants_runid[0]
        return render(request, 'oncoreport/negative_report.html', {'tumor_sample': tumor_sample, 'runid': runid, 'tumor_type': tumor_type})
    else:    
        ## Get selected checkbox variants ##
        res = [sub.replace('/', '') for sub in selected_variants]
        ## query the abberations if id is present is res list ##
        selected_variants_df = pd.DataFrame(list(Allaberrations.objects.filter(id__in=res).values()))
        selected_variants_df['sample'] = tumor_sample
        ## Get the runid to redirect to index page with runid ##
        variants_runid = selected_variants_df.runid.unique()
        runid = variants_runid[0]

        ## Need to categorize df based on DNA results,RNA results and Tiers ##
        variants_table = selected_variants_df[['gene','variants','tier','variant_type','vaf','depth','transcript','locus','exon','copy_number','interpretations','citations']]

        ## DNA results SNV and INDELS ##
        DNA_results = variants_table[(variants_table['variant_type'] == 'SNV') | (variants_table['variant_type'] == 'INDEL') | (variants_table['variant_type'] == 'CNV') | (variants_table['variant_type'] == 'MNV')]
        DNA_results_snv_indels_only = DNA_results[(DNA_results['variant_type'] == 'SNV') | (DNA_results['variant_type'] == 'INDEL') | (DNA_results['variant_type'] == 'MNV')]
        DNA_results_selected_cols = DNA_results_snv_indels_only[['gene','variants','tier','variant_type','vaf','depth','transcript','locus','exon','copy_number']]
        DNA_results_selected_cols_final = DNA_results_selected_cols.copy()
        ## Need only coding and protein change in table for Variant ##
        DNA_results_selected_cols_final['variants'] = DNA_results_selected_cols_final['variants'].str.split(", ").str[1:3].str.join(', ')
        DNA_results_selected_cols_df = DNA_results_selected_cols_final.sort_values('tier')
        DNA_results_final = DNA_results_selected_cols_df.to_dict('records')

        ## DNA results copy number changes ##
        DNA_results_cnvs = variants_table[(variants_table['variant_type'] == 'CNV')]
        DNA_results_selected_CNV = DNA_results_cnvs[['gene','variants','tier','variant_type','vaf','depth','transcript','locus','exon','copy_number','interpretations','citations']]
        DNA_results_selected_CNV['CNV'] = DNA_results_selected_CNV['copy_number'].apply(lambda x: 'Amplification' if x > 8 else 'Gain')
        ## Need only gene and CNV type in table for Variant ##
        DNA_results_selected_CNV['variants'] = DNA_results_selected_CNV['gene'] + ", " + DNA_results_selected_CNV['CNV']
        DNA_results_CNV_final_df = DNA_results_selected_CNV[['gene','CNV']]
        DNA_results_CNV_final = DNA_results_CNV_final_df.to_dict('records')

        ## combine snv, indel and cnv ##
        DNA_all = pd.concat([DNA_results_snv_indels_only,DNA_results_selected_CNV])
        ## DNA Tiers ##
        # change variants for cnv to gene, gain or amp #
        DNA_Tier1or2 = DNA_all[(DNA_all['tier'] == '1/2')]
        DNA_Tier1or2_final= DNA_Tier1or2.to_dict('records')
        DNA_Tier3 = DNA_all[(DNA_all['tier'] == '3')]
        DNA_Tier3_final = DNA_Tier3.to_dict('records')

        ## RNA results ##
        RNA_results = variants_table[(variants_table['variant_type'] == 'FUSION') | (variants_table['variant_type'] == 'RNAExonVariant')]
        RNA_results_final = RNA_results.to_dict('records')

        all_values = {'tumor_sample': tumor_sample,'DNA_results_final': DNA_results_final, 'DNA_results_CNV_final': DNA_results_CNV_final,'DNA_Tier1or2_final':DNA_Tier1or2_final,'DNA_Tier3_final': DNA_Tier3_final,'RNA_results_final': RNA_results_final, 'runid': runid, 'tumor_type': tumor_type}
        return render(request, 'oncoreport/insert_selected_variants.html', all_values)
