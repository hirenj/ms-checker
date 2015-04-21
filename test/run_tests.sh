#!/bin/bash

/Applications/nwjs.app/Contents/MacOS/nwjs .. --debug=true --write-peptides=true "/Scratch/MSF-test/HepG2 TCL T2 (SV)/130114_HepG2_(041012)_SC(L)_T2KO(M)_95%_IEF_1_12_100%_4h_comb_1x2x3x4xGalNAc_Sprot_UniProt_Tyr.msf" > t2_tcl_checked.tsv

/Applications/nwjs.app/Contents/MacOS/nwjs .. --debug=true --write-peptides=true "/Scratch/MSF-test/HepG2 sec T2 (SBL)/130116_HepG2_sec_(051012)_SC(L)_T2KO(M)_95%_IEF_1_12_100%_4h_comb_1x2x3x4xGalNAc_Sprot_UniProt_Tyr.msf" > t2_sec_checked.tsv

extract_excel_sheet "/Scratch/MSF-test/HepG2 TCL T2 (SV)/130114_HepG2_(041012)_SC(L)_T2KO(M)_95%_IEF_1_12_100%_4h_final.xlsx" 2 | awk -F$'\t' 'BEGIN { print "uniprot" FS "sequence" FS "quant" FS "original_quant" FS "glyco" FS "intensity" } NR < 2 { next } { print $3 FS $15 FS $11 FS $10 FS $17 FS $37 }' > t2_tcl_manual.tsv
extract_excel_sheet "/Scratch/MSF-test/HepG2 sec T2 (SBL)/130116_HepG2_sec_(051012)_SC(L)_T2KO(M)_95%_IEF_1_12_100%_4h_final_SBL_mod (2).xlsx" 2 | awk -F$'\t' 'BEGIN { print "uniprot" FS "sequence" FS "quant" FS "original_quant" FS "glyco" FS "intensity" } NR < 2 { next }  { print $3 FS $15 FS $11 FS $10 FS $17 FS $37 }' > t2_sec_manual.tsv


Rscript analysis.R t2_tcl_manual.tsv t2_tcl_checked.tsv t2_tcl_results.xls
Rscript analysis.R t2_sec_manual.tsv t2_sec_checked.tsv t2_sec_results.xls


# select * from SpectrumHeaders LEFT join MassPeaks USING(MassPeakID) LEFT JOIN Peptides USING(SpectrumID) LEFT JOIN PrecursorIonQuanResultsSearchSpectra ON Peptides.SpectrumID = PrecursorIonQuanResultsSearchSpectra.SearchSpectrumID LEFT JOIN EventAnnotations USING(QuanResultID) LEFT JOIN Events using(EventID) where Peptides.Sequence = 'QFTSSTSYNRGDSTFESKSYK' and Peptides.ConfidenceLevel = 3 and Peptides.SearchEngineRank = 1 and EventAnnotations.EventID is not null