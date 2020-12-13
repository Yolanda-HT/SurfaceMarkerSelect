#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 21 13:03:45 2020

@author: yolandatiao
"""

import csv
import pandas as pd
import os
import numpy as np

def greekToEnglish(inStr):
    gTe_dict = {"α": "a", "β":"b", "γ":"g", "δ":"d", "ζ": "z"}
    for g, e in gTe_dict.items():
        inStr = inStr.replace(g, e)
    return(inStr)


#####---------- Unique gene name to alternative gene name reference
if False:
    wk_dir = "/Volumes/Yolanda/Exp391_Acute-Chronic_SC/z_Resources"
    os.chdir(wk_dir)
    
    ref_file = "ENSEMBL_99_Mosue_Genename_GRCm38_p6.txt"
    out_name = "ENSEMBL_99_Mosue_Genename_GRCm38_p6_uniqGenames.csv"
    ref_tb = pd.read_table(ref_file)
    ref_tb['Gene name'] = [str(i).replace(", ", ",") for i in list(ref_tb["Gene name"])]
    ref_tb['Gene name'] = [str(i).replace(",", "|") for i in list(ref_tb["Gene name"])]
    ref_tb['Gene name'] = [str(i).replace("|", "") if str(i).endswith("|") else i for i in list(ref_tb["Gene name"])]
    ref_tb['Gene Synonym'] = [str(i).replace(", ", ",") for i in list(ref_tb["Gene Synonym"])]
    ref_tb['Gene Synonym'] = [str(i).replace(",", "|") for i in list(ref_tb["Gene Synonym"])]
    ref_tb['Gene Synonym'] = [str(i).replace("|", "") if str(i).endswith("|") else i for i in list(ref_tb["Gene Synonym"])]
    drop_cols = ['Gene stable ID', 'Gene stable ID version', 
                 'Transcript stable ID', 'Transcript stable ID version']
    for col in drop_cols:
        del ref_tb[col]
    ref_tb = ref_tb.drop_duplicates()
    gene_name_uniq = list(set(list(ref_tb['Gene name'])))
    gene_sn_list = []
    for gene in gene_name_uniq:
        gene_tb = ref_tb[ref_tb["Gene name"] == gene]
        gene_sn = list(gene_tb["Gene Synonym"])
        gene_sn_str = "|".join(gene_sn)
        gene_sn = gene_sn_str.split("|")
        gene_sn = list(set(gene_sn))
        gene_sn_str = "|".join(gene_sn)
        gene_sn_list.append(gene_sn_str)
    ref_tb_uniq = pd.DataFrame({"gene_name": gene_name_uniq, "gene_synonym":gene_sn_list})    
    ref_tb_uniq.to_csv(out_name, index=False)
    
    ref_file = "ENSEMBL_99_Human_Genename_GRCh38_p13.txt"
    out_name = "ENSEMBL_99_Human_Genename_GRCh38_p13_uniqGenames.csv"
    ref_tb = pd.read_table(ref_file)
    ref_tb['Gene name'] = [str(i).replace(", ", ",") for i in list(ref_tb["Gene name"])]
    ref_tb['Gene name'] = [str(i).replace(",", "|") for i in list(ref_tb["Gene name"])]
    ref_tb['Gene name'] = [str(i).replace("|", "") if str(i).endswith("|") else i for i in list(ref_tb["Gene name"])]
    ref_tb['Gene Synonym'] = [str(i).replace(", ", ",") for i in list(ref_tb["Gene Synonym"])]
    ref_tb['Gene Synonym'] = [str(i).replace(",", "|") for i in list(ref_tb["Gene Synonym"])]
    ref_tb['Gene Synonym'] = [str(i).replace("|", "") if str(i).endswith("|") else i for i in list(ref_tb["Gene Synonym"])]
    drop_cols = ['Gene stable ID', 'Gene stable ID version', 
                 'Transcript stable ID', 'Transcript stable ID version']
    for col in drop_cols:
        del ref_tb[col]
    ref_tb = ref_tb.drop_duplicates()
    gene_name_uniq = list(set(list(ref_tb['Gene name'])))
    gene_sn_list = []
    for gene in gene_name_uniq:
        gene_tb = ref_tb[ref_tb["Gene name"] == gene]
        gene_sn = list(gene_tb["Gene Synonym"])
        gene_sn_str = "|".join(gene_sn)
        gene_sn = gene_sn_str.split("|")
        gene_sn = list(set(gene_sn))
        gene_sn_str = "|".join(gene_sn)
        gene_sn_list.append(gene_sn_str)
    ref_tb_uniq = pd.DataFrame({"gene_name": gene_name_uniq, "gene_synonym":gene_sn_list})    
    ref_tb_uniq.to_csv(out_name, index=False)

#####---------- CD format unify
if False:
    wk_dir = "/Volumes/Yolanda/Exp391_Acute-Chronic_SC/z_Resources"
    os.chdir(wk_dir)
    
    cd_file = "BD_mouse_CDs.csv"
    out_file = "BD_mouse_CDs_fmt.csv"
    with open(out_file, "w") as fout:
        with open(cd_file, "r") as fin:
            rfin = csv.reader(fin, delimiter=",")
            wfout = csv.writer(fout, delimiter=",")
            for row in rfin:
                row_alt = row[1]
                row_alt = row_alt.replace("; ", "|")
                row_alt = row_alt.replace(";", "|")
                row_alt = greekToEnglish(row_alt)
                row_alt = "".join(row_alt.split())
                wfout.writerow([row[0], row_alt])
    
    cd_file = "BD_human_CDs.csv"
    out_file = "BD_human_CDs_fmt.csv"
    with open(out_file, "w") as fout:
        with open(cd_file, "r") as fin:
            rfin = csv.reader(fin, delimiter=",")
            wfout = csv.writer(fout, delimiter=",")
            for row in rfin:
                row_alt = row[1]
                row_alt = row_alt.replace("; ", "|")
                row_alt = row_alt.replace(";", "|")
                row_alt = greekToEnglish(row_alt)
                row_alt = "".join(row_alt.split())
                wfout.writerow([row[0], row_alt])

#####---------- CD to gene name
def cd_to_gn(cd_file, ref_file):
    #cd_file = "LEGENDSCREEN_MOUSE_fmt.csv"
    #ref_file = "ENSEMBL_99_Mosue_Genename_GRCm38_p6_uniqGenames.csv"
    
    cd_tb = pd.read_csv(cd_file)
    ref_tb = pd.read_csv(ref_file)
    
    cd_tb_allGeneNames = list(cd_tb["CD"])
    cd_tb_altnames = [str(i) for i in list(cd_tb["AlternativeName"])]
    cd_tb_altnames = "|".join(cd_tb_altnames).split("|")
    cd_tb_allGeneNames = cd_tb_allGeneNames + cd_tb_altnames
    cd_tb_allGeneNames = [i for i in cd_tb_allGeneNames if i !='']
    cd_tb_allGeneNames = [i for i in cd_tb_allGeneNames if i !='nan']
    cd_tb_allGeneNames = list(set(cd_tb_allGeneNames))
    cd_tb_allGeneNames_lower = [i.lower() for i in cd_tb_allGeneNames]
    
    
    # Keep genes that with name / synoname in CD table
    ref_use = []
    for i in range(0, len(ref_tb)):
        i_names = [ref_tb["gene_name"][i].lower()]
        i_syns = str(ref_tb["gene_synonym"][i]).lower().split("|")
        i_names += i_syns
        i_use = [True if x in cd_tb_allGeneNames_lower else False for x in i_names]
        if any(i_use):
            ref_use.append("y")
        else:
            ref_use.append("n")
    ref_tb["use"] = ref_use
    ref_tb = ref_tb[ref_tb["use"] == "y"]
    ref_tb_gns = list(ref_tb["gene_name"])
    ref_tb_sns = list(ref_tb['gene_synonym'])
    ref_tb_gns_lower = [str(i).lower() for i in ref_tb_gns]
    ref_tb_sns_lower = [str(i).lower() for i in ref_tb_sns]
    
    # Matching CD table to gene reference table
    gene_names = []
    gene_sns = []
    for i in range(0, len(cd_tb)):
    #for i in range(0, 5):
        i_cd = cd_tb["CD"][i].lower()
        i_alt = str(cd_tb["AlternativeName"][i]).lower().split("|")
        
        i_ref_idx = None
        # First order match: CD match with gene_name
        if i_cd in ref_tb_gns_lower:
            i_ref_idx = ref_tb_gns_lower.index(i_cd)
        # Second order match: CD match with gene_synonym
        else:
            for x in range(0, len(ref_tb_sns_lower)):
                x_sn = ref_tb_sns_lower[x].split("|")
                if i_cd in x_sn:
                    i_ref_idx = x
                    break
        if i_ref_idx == None:
            for j in i_alt:
                if (j != "nan"):
                    # Third order match: CD alternative name with gene_name
                    if j in ref_tb_gns_lower:
                        i_ref_idx = ref_tb_gns_lower.index(j)
                        break
                    # Last match: CD alternative name with gene_synonym
                    else:
                        for x in range(0, len(ref_tb_sns_lower)):
                            x_sn = ref_tb_sns_lower[x].split("|")
                            if j in x_sn:
                                i_ref_idx = x
                                break
                    if i_ref_idx != None:
                        break
        if i_ref_idx != None:
            gene_names.append(ref_tb_gns[i_ref_idx])
            gene_sns.append(ref_tb_sns[i_ref_idx])
        else:
            gene_names.append("")
            gene_sns.append("")
    cd_tb["gene_name"] = gene_names
    cd_tb["gene_synonym"] = gene_sns
    cd_tb.to_csv(cd_file.replace(".csv", "_GN.csv"), index=False)


if False:
    # Mouse
    cd_mm = "/Volumes/Yolanda/Exp391_Acute-Chronic_SC/z_Resources/BD_mouse_CDs_fmt.csv"
    ref_mm = "/Volumes/Yolanda/Exp391_Acute-Chronic_SC/z_Resources/ENSEMBL_99_Mosue_Genename_GRCm38_p6_uniqGenames.csv"
    cd_to_gn(cd_mm, ref_mm)
    
    # Human
    cd_h = "/Volumes/Yolanda/Exp391_Acute-Chronic_SC/z_Resources/BD_human_CDs_fmt.csv"
    ref_h = "/Volumes/Yolanda/Exp391_Acute-Chronic_SC/z_Resources/ENSEMBL_99_Human_Genename_GRCh38_p13_uniqGenames.csv"
    cd_to_gn(cd_h, ref_h)

#####---------- Biolegend name to genename
wk_dir = "/Volumes/Yolanda/Exp391_Acute-Chronic_SC/z_Resources"
os.chdir(wk_dir)

###----- Mouse
if False:
    biolegend_file = "LEGENDSCREEN_MOUSE.csv"
    list_file = biolegend_file.replace(".csv", "_list.csv")
    correct_file = biolegend_file.replace(".csv", "_list_corrected.csv")
    use_to_match_file = biolegend_file.replace(".csv", "_fmt.csv")
    
    biolegend_name_list = []
    with open(biolegend_file, "r") as fin:
        rfin = csv.reader(fin, delimiter=",")
        for row in rfin:
            row_use = ["".join(x.split()) for x in row]
            biolegend_name_list += row_use
    
    biolegend_name_list = [i for i in biolegend_name_list if "Ctrl" not in i]
    biolegend_name_list = [i for i in biolegend_name_list if "control" not in i]
    biolegend_name_list = [i for i in biolegend_name_list if "Blank" not in i]
    biolegend_name_list = list(set(biolegend_name_list))
    biolegend_name_list = [greekToEnglish(i) for i in biolegend_name_list]
    
    biolegend_name_tb = pd.DataFrame({"marker_name": biolegend_name_list})
    biolegend_name_tb.to_csv(list_file, index=False)
    
    ### Manually remove "/" --> corrected file
    biolegend_name_tb = pd.read_csv(correct_file)
    marker_name = list(biolegend_name_tb["marker_name"])
    cd_name = []
    alt_name = []
    for i in marker_name:
        if "(" in i:
            cd_i = i.replace(")", "").split("(")[0]
            alt_i = i.replace(")", "").split("(")[1]
            cd_name.append(cd_i)
            alt_name.append(alt_i)
        elif "-" in i:
            cd_name.append(i)
            alt_name.append(i.replace("-", ""))
        else:
            cd_name.append(i)
            alt_name.append("")
    biolegend_name_tb_fmt = pd.DataFrame({"CD": cd_name, "AlternativeName": alt_name})
    biolegend_name_tb_fmt.to_csv(use_to_match_file, index=False)


###----- Human
if False:
    biolegend_file = "LEGENDSCREEN_HUMAN.csv"
    list_file = biolegend_file.replace(".csv", "_list.csv")
    correct_file = biolegend_file.replace(".csv", "_list_corrected.csv")
    use_to_match_file = biolegend_file.replace(".csv", "_fmt.csv")
    
    biolegend_name_list = []
    with open(biolegend_file, "r") as fin:
        rfin = csv.reader(fin, delimiter=",")
        for row in rfin:
            row_use = ["".join(x.split()) for x in row]
            biolegend_name_list += row_use
    
    biolegend_name_list = [i for i in biolegend_name_list if "Ctrl" not in i]
    biolegend_name_list = [i for i in biolegend_name_list if "control" not in i]
    biolegend_name_list = [i for i in biolegend_name_list if "Blank" not in i]
    biolegend_name_list = list(set(biolegend_name_list))
    biolegend_name_list = [greekToEnglish(i) for i in biolegend_name_list]
    
    biolegend_name_tb = pd.DataFrame({"marker_name": biolegend_name_list})
    biolegend_name_tb.to_csv(list_file, index=False)
    
    ### Manually remove "/" --> corrected file
    biolegend_name_tb = pd.read_csv(correct_file)
    marker_name = list(biolegend_name_tb["marker_name"])
    cd_name = []
    alt_name = []
    for i in marker_name:
        if "(" in i:
            cd_i = i.replace(")", "").split("(")[0]
            alt_i = i.replace(")", "").split("(")[1]
            cd_name.append(cd_i)
            alt_name.append(alt_i)
        elif "-" in i:
            cd_name.append(i)
            alt_name.append(i.replace("-", ""))
        else:
            cd_name.append(i)
            alt_name.append("")
    biolegend_name_tb_fmt = pd.DataFrame({"CD": cd_name, "AlternativeName": alt_name})
    biolegend_name_tb_fmt.to_csv(use_to_match_file, index=False)

if False:
    # Mouse
    biolegend_mm = "LEGENDSCREEN_MOUSE_fmt.csv"
    ref_mm = "/Volumes/Yolanda/Exp391_Acute-Chronic_SC/z_Resources/ENSEMBL_99_Mosue_Genename_GRCm38_p6_uniqGenames.csv"
    cd_to_gn(biolegend_mm, ref_mm)
    
    # Human
    biolegend_h = "LEGENDSCREEN_HUMAN_fmt.csv"
    ref_h = "/Volumes/Yolanda/Exp391_Acute-Chronic_SC/z_Resources/ENSEMBL_99_Human_Genename_GRCh38_p13_uniqGenames.csv"
    cd_to_gn(biolegend_h, ref_h)

###----- Fill in gene names that were not matched with BD matahced gene names
if False:
    biolegend_gn = "LEGENDSCREEN_MOUSE_fmt_GN.csv"
    bd_gn = "BD_mouse_CDs_fmt_GN.csv"
    
    biolegend_gn_add = biolegend_gn.replace(".csv", "_add.csv")
    biolegend_gn_tb = pd.read_csv(biolegend_gn)
    bd_gn_tb = pd.read_csv(bd_gn)
    
    biolegend_gn_tb_matched = biolegend_gn_tb[ biolegend_gn_tb["gene_name"].notna()]
    biolegend_gn_tb_notmatched = biolegend_gn_tb[ biolegend_gn_tb["gene_name"].isna()]
    
    del biolegend_gn_tb_notmatched["gene_name"]
    del biolegend_gn_tb_notmatched["gene_synonym"]
    del biolegend_gn_tb_notmatched["AlternativeName"]
    
    bd_gn_tb.columns
    biolegend_gn_tb_notmatched.columns
    
    biolegend_gn_tb_notmatched = biolegend_gn_tb_notmatched.set_index("CD").join(
            bd_gn_tb.set_index("CD"), how="left")
    biolegend_gn_tb_notmatched = biolegend_gn_tb_notmatched.reset_index()
    
    biolegend_gn_tb_notmatched["CD"]
    biolegend_gn_tb_notmatched["gene_name"]
    biolegend_gn_tb = biolegend_gn_tb_matched.append(biolegend_gn_tb_notmatched)
    
    biolegend_gn_tb.to_csv(biolegend_gn_add, index=False)


###----- Create final list of usable markers
wkdir = "/Volumes/Yolanda/Exp391_Acute-Chronic_SC/z_Resources"
os.chdir(wkdir)

bd_file = "BD_mouse_CDs_fmt_GN_filled.csv"
biolegend_file = "LEGENDSCREEN_MOUSE_fmt_GN_add_filled.csv"

bd_tb = pd.read_csv(bd_file)
biolegend_tb = pd.read_csv(biolegend_file)

use_markers = list(bd_tb["gene_name"]) + list(biolegend_tb["gene_name"])
use_markers = [i for i in use_markers if str(i) != "nan"]
use_markers = list(set(use_markers))

out_file = "MM_MARKERS.csv"
with open(out_file, "w") as fout:
    wfout = csv.writer(fout, delimiter=",")
    wfout.writerow(["gene_name"])
    for i in use_markers:
        wfout.writerow([i])






















