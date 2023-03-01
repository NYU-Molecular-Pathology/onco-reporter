import pandas as pd
import re

def process_files(onco_report_df_DNA,PMKB,tumor_type,type):
    onco_values = pd.DataFrame()
    pmkb_values = pd.DataFrame()

    reshaped_df_orig = onco_report_df_DNA
    reshaped_df = reshaped_df_orig.copy()
    
    reshaped_df['Gene'] = reshaped_df['genes'] ## to get fusion gene original name

    ## For RNA samples example: BCR(2)-ABL1(4) to BCR(2) and ABL1(4) as separate rows ##
    reshaped_df_final = reshaped_df.set_index(reshaped_df.columns.drop('genes',1).tolist()).genes.str.split(" - ", expand=True).stack().reset_index().rename(columns={0:'genes'}).loc[:,reshaped_df.columns]

    ## Add Variants column to get oncomine report values ##
    if type == "DNA":
        reshaped_df['Variants'] = reshaped_df[['genes','coding','aa']].apply(lambda x: ', '.join(x), axis=1)
        onco_values = onco_values.append(reshaped_df)
    else:
        reshaped_df_RNA = reshaped_df_final
        reshaped_df_RNA['Variants'] = "NA" # No need to combine variants (no AA change)
        reshaped_df_RNA['genes'] = reshaped_df_RNA['genes'].str.replace(r" ?\([^)]+\)","",regex=True) ## Remove exon numbers and brackets from gene names
        onco_values = onco_values.append(reshaped_df_RNA)
    
    PMKB_tumortype_final = PMKB.loc[(PMKB['tumor_type'] == tumor_type)]
    PMKB_tumortype = PMKB_tumortype_final.copy()
    
    AA_function = lambda exon, gene, varianttype, AA: "Exon"+" " +str(int(exon)) if gene == "EGFR" and varianttype == "INDEL" else ( "Amplification" if varianttype == "CNV" else AA.replace('p.',''))

    for gene,AA,varianttype,exon in zip(reshaped_df_final.genes,reshaped_df_final.aa,reshaped_df_final.variant_type,reshaped_df_final.exon):
        AA_final = str(AA_function(exon,gene,varianttype,AA))
        if type == "DNA":
            PMKB_result = (PMKB_tumortype.loc[(PMKB_tumortype['gene'] == gene) & (PMKB_tumortype['tumor_type'] == tumor_type) & (PMKB_tumortype['variant'].str.startswith(AA_final))]) if gene in PMKB_tumortype.values and tumor_type in PMKB_tumortype.values and any(PMKB_tumortype.variant.str.startswith(AA_final)) else ((PMKB_tumortype.loc[(PMKB_tumortype['gene'] == gene) & (PMKB_tumortype['tumor_type'] == tumor_type) & (PMKB_tumortype['variant'] == "Default")] if gene in PMKB_tumortype.values and tumor_type in PMKB_tumortype.values and any(~PMKB_tumortype.variant.str.startswith(AA_final)) else((PMKB.loc[(PMKB['gene'] == gene) & (PMKB['tumor_type'] == "Other") & (PMKB['variant'].str.startswith(AA_final))] if gene not in PMKB_tumortype.values and any(PMKB.variant.str.startswith(AA_final)) else(PMKB.loc[(PMKB['gene'] == gene) & (PMKB['tumor_type'] == "Other")])))))
            pmkb_values = pmkb_values.append(PMKB_result)
        else:
            remove_exons_gene = re.sub(r" ?\([^)]+\)", "", gene)
            gene_final = remove_exons_gene
            PMKB_result = (PMKB_tumortype.loc[(PMKB_tumortype['gene'] == gene_final) & (PMKB_tumortype['variant'].str.contains("Fusion"))]) if gene_final in PMKB_tumortype.values and any(PMKB_tumortype.variant.str.contains("Fusion")) else (PMKB.loc[(PMKB['gene'] == gene) & (PMKB['tumor_type'] == "Other")])
            pmkb_values = pmkb_values.append(PMKB_result)

    pmkb_values_nodups = pmkb_values.drop_duplicates(subset=['gene','variant','interpretations','citations'])
    pmkb_values_final = pmkb_values_nodups.copy()
    pmkb_values_final['Type_pmkb'] = pmkb_values_final['variant'].apply(lambda variant_type: 'CNV' if variant_type == "Amplification" else ( "INDEL" if "deletion" in variant_type or "insertion" in variant_type else "SNV"))

    merged_table_dna = pd.merge(onco_values,pmkb_values_final,how='left',left_on=['genes','variant_type'], right_on=['gene','Type_pmkb'])
    merged_table_rna = pd.merge(onco_values,pmkb_values_final,how='left',left_on='genes', right_on='gene')
    merged_table_rna_final = merged_table_rna.drop_duplicates(subset=["locus"],keep="last")
    return merged_table_dna,merged_table_rna_final


def all_results(onco_report,PMKB,selected_sample,tumor_type):
    onco_report_persample = onco_report[(onco_report['sample'] == selected_sample)]
    ## If Sample has SNV,INDEL then combine DNA+RNA 
    onco_report_persample_type = onco_report_persample.variant_type.unique()
    type_status = any(x in onco_report_persample_type for x in ['SNV', 'INDEL', 'CNV', 'MNV'])
    ## Checking if Genes contians NEGATIVE ##
    negative_found = onco_report_persample[onco_report_persample['genes'].str.contains('NEGATIVE')]
    
    if type_status == True: 
        ## Get both DNA and RNA results in other cases which have both info ##
        ### Get DNA results table ###
        onco_report_df_DNA = onco_report_persample[(onco_report_persample['variant_type'] != 'RNAExonVariant') & (onco_report_persample['variant_type'] != 'FUSION')]

        DNA_results_processed = process_files(onco_report_df_DNA,PMKB,tumor_type,"DNA")[0]
        DNA_results_final = DNA_results_processed[['runid','sample','Gene','variant_type','transcript','Variants','frequency','exon','aa','copy_number','coverage','locus','tumor_type','tier','interpretations','citations','variant']]
        DNA_results = DNA_results_final.drop_duplicates()

        onco_report_df_RNA = onco_report_persample.loc[(onco_report_persample['variant_type'] == 'RNAExonVariant') | (onco_report_persample['variant_type'] == 'FUSION')]
        if onco_report_df_RNA.empty:
            DNA_RNA = DNA_results
            abberations_all = DNA_RNA.drop_duplicates()
            return abberations_all
        else:
            RNA_results_processed = process_files(onco_report_df_RNA,PMKB,tumor_type,"RNA")[1]
            RNA_results_final = RNA_results_processed[['runid','sample','Gene','variant_type','transcript','Variants','frequency','exon','aa','copy_number','coverage','locus','tumor_type','tier','interpretations','citations','variant']]
            RNA_results = RNA_results_final.drop_duplicates()
        
            DNA_RNA = DNA_results.append(RNA_results)
            abberations_all = DNA_RNA.drop_duplicates()
            return abberations_all

    elif  negative_found['genes'].count() > 0:
        runid = onco_report_persample.runid.unique()
        onco_negative_result = pd.DataFrame(columns=['runid','sample','Gene','variant_type','transcript','Variants','frequency','exon','aa','copy_number','coverage','locus','tumor_type','tier','interpretations','citations','variant'])
        onco_negative_result = onco_negative_result.append({'runid': runid[0], 'sample': selected_sample, 'Gene': "negative",'variant_type' : "negative",'transcript': "negative",'Variants': "negative",'frequency': 0.0,'exon':0,'aa': "negative",'copy_number' : 0.0,'coverage': 0,'locus' : "negative",'tumor_type': tumor_type,'tier': 0,'interpretations': "negative",'citations': "negative",'variant': "negative"}, ignore_index=True)
        return onco_negative_result

    else:
        onco_report_df_RNA = onco_report_persample[(onco_report_persample['variant_type'] == 'RNAExonVariant') | (onco_report_persample['variant_type'] == 'FUSION')]
        RNA_results_processed = process_files(onco_report_df_RNA,PMKB,tumor_type,"RNA")[1]
        RNA_results_final = RNA_results_processed[['runid','sample','Gene','variant_type','transcript','Variants','frequency','exon','aa','copy_number','coverage','locus','tumor_type','tier','interpretations','citations','variant']]
        RNA_results = RNA_results_final.drop_duplicates()
        return RNA_results


def main(onco_report,PMKB,sample,tumor_type):
    abberation_per_sample = all_results(onco_report,PMKB,sample,tumor_type)
    abberation_per_sample_nona = abberation_per_sample.fillna({'vaf': 0.0, 'exon':0, 'tier':'1/2'})
    #abberation_per_sample_nona.to_csv("abberation_per_sample_nona.csv")
    input_all_abberations = abberation_per_sample_nona[['runid','sample','Gene','variant_type','transcript','Variants','frequency','aa','exon','copy_number','coverage','tier','locus','variant','tumor_type','interpretations','citations']]
    input_all_abberations_tiering = input_all_abberations.copy()
    input_all_abberations_final = input_all_abberations_tiering.reset_index(drop=True)
    return input_all_abberations_final
