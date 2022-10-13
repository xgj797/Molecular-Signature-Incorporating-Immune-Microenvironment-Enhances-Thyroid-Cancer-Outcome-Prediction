
cd '/scratch/weissvl/shengq2/20210624_human_wntpathway_rnaseq/deseq2_proteincoding_genetable_WebGestalt/result'

set -o pipefail


 
if [[ ! -s /scratch/weissvl/shengq2/20210624_human_wntpathway_rnaseq/deseq2_proteincoding_genetable_WebGestalt/result/malignant_vs_benign/Project_malignant_vs_benign_pathway_KEGG/enrichment_results_malignant_vs_benign_pathway_KEGG.txt || ! -d /scratch/weissvl/shengq2/20210624_human_wntpathway_rnaseq/deseq2_proteincoding_genetable_WebGestalt/result/malignant_vs_benign/Project_malignant_vs_benign_pathway_KEGG/enrichment_results_malignant_vs_benign_pathway_KEGG.txt ]]; then
  cd /scratch/weissvl/shengq2/20210624_human_wntpathway_rnaseq/deseq2_proteincoding_genetable_WebGestalt/result/malignant_vs_benign 
  if [[ -f /scratch/weissvl/shengq2/20210624_human_wntpathway_rnaseq/deseq2_proteincoding_genetable/result/malignant_vs_benign_min5_fdr0.05_DESeq2_sig_genename.txt ]]; then
    if [[ -s /scratch/weissvl/shengq2/20210624_human_wntpathway_rnaseq/deseq2_proteincoding_genetable/result/malignant_vs_benign_min5_fdr0.05_DESeq2_sig_genename.txt ]]; then
      R --vanilla -f /scratch/weissvl/shengq2/20210624_human_wntpathway_rnaseq/deseq2_proteincoding_genetable_WebGestalt/result/weiss_human_wnt.WebGestaltR.r --args hsapiens malignant_vs_benign /scratch/weissvl/shengq2/20210624_human_wntpathway_rnaseq/deseq2_proteincoding_genetable/result/malignant_vs_benign_min5_fdr0.05_DESeq2_sig_genename.txt . genesymbol genome
      rm */*/*.zip
    else
      echo "Empty gene file" > malignant_vs_benign.empty
    fi 
  else
    echo "Gene file not exist: /scratch/weissvl/shengq2/20210624_human_wntpathway_rnaseq/deseq2_proteincoding_genetable/result/malignant_vs_benign_min5_fdr0.05_DESeq2_sig_genename.txt" .
  fi
fi

 
if [[ ! -s /scratch/weissvl/shengq2/20210624_human_wntpathway_rnaseq/deseq2_proteincoding_genetable_WebGestalt/result/metastatic_vs_primary/Project_metastatic_vs_primary_pathway_KEGG/enrichment_results_metastatic_vs_primary_pathway_KEGG.txt || ! -d /scratch/weissvl/shengq2/20210624_human_wntpathway_rnaseq/deseq2_proteincoding_genetable_WebGestalt/result/metastatic_vs_primary/Project_metastatic_vs_primary_pathway_KEGG/enrichment_results_metastatic_vs_primary_pathway_KEGG.txt ]]; then
  cd /scratch/weissvl/shengq2/20210624_human_wntpathway_rnaseq/deseq2_proteincoding_genetable_WebGestalt/result/metastatic_vs_primary 
  if [[ -f /scratch/weissvl/shengq2/20210624_human_wntpathway_rnaseq/deseq2_proteincoding_genetable/result/metastatic_vs_primary_min5_fdr0.05_DESeq2_sig_genename.txt ]]; then
    if [[ -s /scratch/weissvl/shengq2/20210624_human_wntpathway_rnaseq/deseq2_proteincoding_genetable/result/metastatic_vs_primary_min5_fdr0.05_DESeq2_sig_genename.txt ]]; then
      R --vanilla -f /scratch/weissvl/shengq2/20210624_human_wntpathway_rnaseq/deseq2_proteincoding_genetable_WebGestalt/result/weiss_human_wnt.WebGestaltR.r --args hsapiens metastatic_vs_primary /scratch/weissvl/shengq2/20210624_human_wntpathway_rnaseq/deseq2_proteincoding_genetable/result/metastatic_vs_primary_min5_fdr0.05_DESeq2_sig_genename.txt . genesymbol genome
      rm */*/*.zip
    else
      echo "Empty gene file" > metastatic_vs_primary.empty
    fi 
  else
    echo "Gene file not exist: /scratch/weissvl/shengq2/20210624_human_wntpathway_rnaseq/deseq2_proteincoding_genetable/result/metastatic_vs_primary_min5_fdr0.05_DESeq2_sig_genename.txt" .
  fi
fi

 
if [[ ! -s /scratch/weissvl/shengq2/20210624_human_wntpathway_rnaseq/deseq2_proteincoding_genetable_WebGestalt/result/poor_outcome_yes_vs_no/Project_poor_outcome_yes_vs_no_pathway_KEGG/enrichment_results_poor_outcome_yes_vs_no_pathway_KEGG.txt || ! -d /scratch/weissvl/shengq2/20210624_human_wntpathway_rnaseq/deseq2_proteincoding_genetable_WebGestalt/result/poor_outcome_yes_vs_no/Project_poor_outcome_yes_vs_no_pathway_KEGG/enrichment_results_poor_outcome_yes_vs_no_pathway_KEGG.txt ]]; then
  cd /scratch/weissvl/shengq2/20210624_human_wntpathway_rnaseq/deseq2_proteincoding_genetable_WebGestalt/result/poor_outcome_yes_vs_no 
  if [[ -f /scratch/weissvl/shengq2/20210624_human_wntpathway_rnaseq/deseq2_proteincoding_genetable/result/poor_outcome_yes_vs_no_min5_fdr0.05_DESeq2_sig_genename.txt ]]; then
    if [[ -s /scratch/weissvl/shengq2/20210624_human_wntpathway_rnaseq/deseq2_proteincoding_genetable/result/poor_outcome_yes_vs_no_min5_fdr0.05_DESeq2_sig_genename.txt ]]; then
      R --vanilla -f /scratch/weissvl/shengq2/20210624_human_wntpathway_rnaseq/deseq2_proteincoding_genetable_WebGestalt/result/weiss_human_wnt.WebGestaltR.r --args hsapiens poor_outcome_yes_vs_no /scratch/weissvl/shengq2/20210624_human_wntpathway_rnaseq/deseq2_proteincoding_genetable/result/poor_outcome_yes_vs_no_min5_fdr0.05_DESeq2_sig_genename.txt . genesymbol genome
      rm */*/*.zip
    else
      echo "Empty gene file" > poor_outcome_yes_vs_no.empty
    fi 
  else
    echo "Gene file not exist: /scratch/weissvl/shengq2/20210624_human_wntpathway_rnaseq/deseq2_proteincoding_genetable/result/poor_outcome_yes_vs_no_min5_fdr0.05_DESeq2_sig_genename.txt" .
  fi
fi

