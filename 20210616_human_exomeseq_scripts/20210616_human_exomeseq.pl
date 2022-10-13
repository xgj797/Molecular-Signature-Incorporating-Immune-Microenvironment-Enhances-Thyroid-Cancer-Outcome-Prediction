#!/usr/bin/perl
use strict;
use warnings;

use CQS::Global;
use CQS::FileUtils;
use CQS::SystemUtils;
use CQS::ConfigUtils;
use CQS::ClassFactory;
use CQS::PerformExomeSeq;
use Data::Dumper;

my $intervals = "/data/cqs/references/exomeseq/IDT/Exome-IDT_V1.hg38.slop50.bed.interval_list";
my $bait_intervals = "/data/cqs/references/exomeseq/IDT/Exome-IDT_V1.hg38.slop50.bed.interval_list";
my $cover_bed = "/data/cqs/references/exomeseq/IDT/Exome-IDT_V1.hg38.slop50.bed";

my $def = {
  task_name => "human_exomeseq",

  target_dir => "/scratch/weissvl/shengq2/20210616_human_exomeseq",
  #target_dir => "/workspace/shengq2/vivian_weiss/20210616_human_exomeseq",
  
  email      => "quanhu.sheng.1\@vumc.org",
  
  __ignore_samples => [qw(P4601_CP_01 
P4601_CP_02
P4601_CP_03
P4601_CP_04
P4601_CP_05
P4601_CP_06
P4601_CP_07
P4601_CP_08
P4601_CP_09
P4601_CP_10
P4601_CP_11
P4601_CP_12
P4601_CP_13
P4601_CP_14
P4601_CP_15
P4601_CP_16
P4601_CP_17
P4601_CP_19
P4601_CP_20
P4601_CP_21
P4601_CP_22
P4601_CP_23
P4601_CP_24
P4601_CP_25
P4601_CP_26
P4601_CP_27
P4601_CP_28
P4601_CP_29
P4601_CP_30
P4601_CP_31
P4601_CP_32
P4601_CP_33
P4601_CP_34
P4601_CP_35
P4601_CP_36
P4601_CP_37
P4601_CP_38
P4601_CP_39
P4601_CP_40
P4601_CP_41
P4601_CP_42
P4601_CP_43
P4601_CP_44
P4601_CP_45
P4601_CP_46
P4601_CP_47
P4601_CP_48
P4601_CP_49
P4601_CP_50
P1896_AC_027
P1809_AC_031
P1809_AC_032
P1809_AC_033
P1809_AC_034
P1809_AC_035
P1809_AC_036
P1809_AC_037
P1809_AC_038
P1809_AC_039
P1809_AC_040
P1809_AC_041
P1809_AC_042
P1809_AC_043
P1809_AC_044
P1809_AC_045
P1809_AC_046
P1809_AC_047
P1809_AC_048
P1809_AC_049
P1809_AC_050
P1809_AC_122
P1809_AC_123
P1809_AC_124
P1809_AC_125
P1809_AC_126
P1809_AC_127
P1809_AC_128
P1809_AC_129
P1809_AC_130
P1809_AC_131
P1809_AC_132
P1809_AC_133
P1809_AC_134
P1809_AC_135
)],
  files      => {
    #file_def.py -i /data/h_vivian_weiss/DNAseq/1809/ -n "(1809....\\d+)" -f "*.gz" -a -p
    'P1809_AC_003' => [ '/data/h_vivian_weiss/DNAseq/1809/1809-AC-3-GGATATATCC-TGTTCACCAT_S3_R1_001.fastq.gz', '/data/h_vivian_weiss/DNAseq/1809/1809-AC-3-GGATATATCC-TGTTCACCAT_S3_R2_001.fastq.gz' ],
    'P1809_AC_004' => [ '/data/h_vivian_weiss/DNAseq/1809/1809-AC-4-AAGCGCGCTT-CACCTGTTGC_S4_R1_001.fastq.gz', '/data/h_vivian_weiss/DNAseq/1809/1809-AC-4-AAGCGCGCTT-CACCTGTTGC_S4_R2_001.fastq.gz' ],
    'P1809_AC_005' => [ '/data/h_vivian_weiss/DNAseq/1809/1809-AC-5-CAACGAGAGC-GTAGGTGGTG_S5_R1_001.fastq.gz', '/data/h_vivian_weiss/DNAseq/1809/1809-AC-5-CAACGAGAGC-GTAGGTGGTG_S5_R2_001.fastq.gz' ],
    'P1809_AC_006' => [ '/data/h_vivian_weiss/DNAseq/1809/1809-AC-6-TGGTAGAGAT-ACGAACAACA_S6_R1_001.fastq.gz', '/data/h_vivian_weiss/DNAseq/1809/1809-AC-6-TGGTAGAGAT-ACGAACAACA_S6_R2_001.fastq.gz' ],
    'P1809_AC_008' => [ '/data/h_vivian_weiss/DNAseq/1809/1809-AC-8-TTCGGTGTGA-ACGCCTTGTT_S7_R1_001.fastq.gz', '/data/h_vivian_weiss/DNAseq/1809/1809-AC-8-TTCGGTGTGA-ACGCCTTGTT_S7_R2_001.fastq.gz' ],
    'P1809_AC_010' => [ '/data/h_vivian_weiss/DNAseq/1809/1809-AC-10-CCTAACACAG-GTATTCCACC_S8_R1_001.fastq.gz', '/data/h_vivian_weiss/DNAseq/1809/1809-AC-10-CCTAACACAG-GTATTCCACC_S8_R2_001.fastq.gz' ],
    'P1809_AC_011' => [ '/data/h_vivian_weiss/DNAseq/1809/1809-AC-11-GGTAGAATTA-TGTAAGGTGG_S9_R1_001.fastq.gz', '/data/h_vivian_weiss/DNAseq/1809/1809-AC-11-GGTAGAATTA-TGTAAGGTGG_S9_R2_001.fastq.gz' ],
    'P1809_AC_012' => [ '/data/h_vivian_weiss/DNAseq/1809/1809-AC-12-AACGAGGCCG-CACGGAACAA_S10_R1_001.fastq.gz', '/data/h_vivian_weiss/DNAseq/1809/1809-AC-12-AACGAGGCCG-CACGGAACAA_S10_R2_001.fastq.gz' ],
    'P1809_AC_013' => [ '/data/h_vivian_weiss/DNAseq/1809/1809-AC-13-TATCCGAGGC-CTCAACGCTT_S11_R1_001.fastq.gz', '/data/h_vivian_weiss/DNAseq/1809/1809-AC-13-TATCCGAGGC-CTCAACGCTT_S11_R2_001.fastq.gz' ],
    'P1809_AC_014' => [ '/data/h_vivian_weiss/DNAseq/1809/1809-AC-14-CGCTTAGAAT-TCTGGTATCC_S12_R1_001.fastq.gz', '/data/h_vivian_weiss/DNAseq/1809/1809-AC-14-CGCTTAGAAT-TCTGGTATCC_S12_R2_001.fastq.gz' ],
    'P1809_AC_015' => [ '/data/h_vivian_weiss/DNAseq/1809/1809-AC-15-ACGGTCCAAC-GAGTCATAGG_S13_R1_001.fastq.gz', '/data/h_vivian_weiss/DNAseq/1809/1809-AC-15-ACGGTCCAAC-GAGTCATAGG_S13_R2_001.fastq.gz' ],
    'P1809_AC_016' => [ '/data/h_vivian_weiss/DNAseq/1809/1809-AC-16-GTAACTTGGT-AGACTGCGAA_S14_R1_001.fastq.gz', '/data/h_vivian_weiss/DNAseq/1809/1809-AC-16-GTAACTTGGT-AGACTGCGAA_S14_R2_001.fastq.gz' ],
    'P1809_AC_017' => [ '/data/h_vivian_weiss/DNAseq/1809/1809-AC-17-TTATACGCGA-AGCTGGAATG_S15_R1_001.fastq.gz', '/data/h_vivian_weiss/DNAseq/1809/1809-AC-17-TTATACGCGA-AGCTGGAATG_S15_R2_001.fastq.gz' ],
    'P1809_AC_018' => [ '/data/h_vivian_weiss/DNAseq/1809/1809-AC-18-CCGCGTATAG-GATCAAGGCA_S16_R1_001.fastq.gz', '/data/h_vivian_weiss/DNAseq/1809/1809-AC-18-CCGCGTATAG-GATCAAGGCA_S16_R2_001.fastq.gz' ],
    'P1809_AC_019' => [ '/data/h_vivian_weiss/DNAseq/1809/1809-AC-19-GACGTCTGCA-TCAGTCTCGT_S17_R1_001.fastq.gz', '/data/h_vivian_weiss/DNAseq/1809/1809-AC-19-GACGTCTGCA-TCAGTCTCGT_S17_R2_001.fastq.gz' ],
    'P1809_AC_020' => [ '/data/h_vivian_weiss/DNAseq/1809/1809-AC-20-AGTACTCATG-CTGACTCTAC_S18_R1_001.fastq.gz', '/data/h_vivian_weiss/DNAseq/1809/1809-AC-20-AGTACTCATG-CTGACTCTAC_S18_R2_001.fastq.gz' ],
    'P1809_AC_021' => [ '/data/h_vivian_weiss/DNAseq/1809/1809-AC-21-CATTAACTGA-CGAAGATTCT_S19_R1_001.fastq.gz', '/data/h_vivian_weiss/DNAseq/1809/1809-AC-21-CATTAACTGA-CGAAGATTCT_S19_R2_001.fastq.gz' ],
    'P1809_AC_022' => [ '/data/h_vivian_weiss/DNAseq/1809/1809-AC-22-TGCCGGTCAG-TAGGAGCCTC_S20_R1_001.fastq.gz', '/data/h_vivian_weiss/DNAseq/1809/1809-AC-22-TGCCGGTCAG-TAGGAGCCTC_S20_R2_001.fastq.gz' ],
    'P1809_AC_023' => [ '/data/h_vivian_weiss/DNAseq/1809/1809-AC-23-TCAGATTAAC-CCTCTACATG_S21_R1_001.fastq.gz', '/data/h_vivian_weiss/DNAseq/1809/1809-AC-23-TCAGATTAAC-CCTCTACATG_S21_R2_001.fastq.gz' ],
    'P1809_AC_024' => [ '/data/h_vivian_weiss/DNAseq/1809/1809-AC-24-CTGAGCCGGT-TTCTCGTGCA_S22_R1_001.fastq.gz', '/data/h_vivian_weiss/DNAseq/1809/1809-AC-24-CTGAGCCGGT-TTCTCGTGCA_S22_R2_001.fastq.gz' ],
    'P1809_AC_025' => [ '/data/h_vivian_weiss/DNAseq/1809/1809-AC-25-GTGACGGAGC-GCGTTGGTAT_S23_R1_001.fastq.gz', '/data/h_vivian_weiss/DNAseq/1809/1809-AC-25-GTGACGGAGC-GCGTTGGTAT_S23_R2_001.fastq.gz' ],
    'P1809_AC_026' => [ '/data/h_vivian_weiss/DNAseq/1809/1809-AC-26-ACAGTAAGAT-ATACCAACGC_S24_R1_001.fastq.gz', '/data/h_vivian_weiss/DNAseq/1809/1809-AC-26-ACAGTAAGAT-ATACCAACGC_S24_R2_001.fastq.gz' ],
    'P1809_AC_027' => [ '/data/h_vivian_weiss/DNAseq/1809/1809-AC-27-CACTCAATTC-AATAGCTGAG_S25_R1_001.fastq.gz', '/data/h_vivian_weiss/DNAseq/1809/1809-AC-27-CACTCAATTC-AATAGCTGAG_S25_R2_001.fastq.gz' ],
    'P1809_AC_028' => [ '/data/h_vivian_weiss/DNAseq/1809/1809-AC-28-TGTCTGGCCT-GGCGATCAGA_S26_R1_001.fastq.gz', '/data/h_vivian_weiss/DNAseq/1809/1809-AC-28-TGTCTGGCCT-GGCGATCAGA_S26_R2_001.fastq.gz' ],
    'P1809_AC_029' => [ '/data/h_vivian_weiss/DNAseq/1809/1809-AC-29-GAACGCAATA-CAATATAGGT_S27_R1_001.fastq.gz', '/data/h_vivian_weiss/DNAseq/1809/1809-AC-29-GAACGCAATA-CAATATAGGT_S27_R2_001.fastq.gz' ],
    'P1809_AC_030' => [ '/data/h_vivian_weiss/DNAseq/1809/1809-AC-30-AGGTATGGCG-TGGCGCGAAC_S28_R1_001.fastq.gz', '/data/h_vivian_weiss/DNAseq/1809/1809-AC-30-AGGTATGGCG-TGGCGCGAAC_S28_R2_001.fastq.gz' ],
    'P1809_AC_031' => [ '/data/h_vivian_weiss/DNAseq/1809/1809-AC-31-CGTACAGGAA-ATGGTTGACT_S29_R1_001.fastq.gz', '/data/h_vivian_weiss/DNAseq/1809/1809-AC-31-CGTACAGGAA-ATGGTTGACT_S29_R2_001.fastq.gz' ],
    'P1809_AC_032' => [ '/data/h_vivian_weiss/DNAseq/1809/1809-AC-32-TACGTGAAGG-GCAACCAGTC_S30_R1_001.fastq.gz', '/data/h_vivian_weiss/DNAseq/1809/1809-AC-32-TACGTGAAGG-GCAACCAGTC_S30_R2_001.fastq.gz' ],
    'P1809_AC_033' => [ '/data/h_vivian_weiss/DNAseq/1809/1809-AC-33-CCGGTTCCTA-TCCATTGCCG_S31_R1_001.fastq.gz', '/data/h_vivian_weiss/DNAseq/1809/1809-AC-33-CCGGTTCCTA-TCCATTGCCG_S31_R2_001.fastq.gz' ],
    'P1809_AC_034' => [ '/data/h_vivian_weiss/DNAseq/1809/1809-AC-34-TTAACCTTCG-CTTGCCATTA_S32_R1_001.fastq.gz', '/data/h_vivian_weiss/DNAseq/1809/1809-AC-34-TTAACCTTCG-CTTGCCATTA_S32_R2_001.fastq.gz' ],
    'P1809_AC_035' => [ '/data/h_vivian_weiss/DNAseq/1809/1809-AC-35-ATCAGTACCA-AGCGAATTAG_S33_R1_001.fastq.gz', '/data/h_vivian_weiss/DNAseq/1809/1809-AC-35-ATCAGTACCA-AGCGAATTAG_S33_R2_001.fastq.gz' ],
    'P1809_AC_036' => [ '/data/h_vivian_weiss/DNAseq/1809/1809-AC-36-GCTGACGTTG-GATAGGCCGA_S34_R1_001.fastq.gz', '/data/h_vivian_weiss/DNAseq/1809/1809-AC-36-GCTGACGTTG-GATAGGCCGA_S34_R2_001.fastq.gz' ],
    'P1809_AC_037' => [ '/data/h_vivian_weiss/DNAseq/1809/1809-AC-37-ATACTTGTTC-GTGCAGACAG_S35_R1_001.fastq.gz', '/data/h_vivian_weiss/DNAseq/1809/1809-AC-37-ATACTTGTTC-GTGCAGACAG_S35_R2_001.fastq.gz' ],
    'P1809_AC_038' => [ '/data/h_vivian_weiss/DNAseq/1809/1809-AC-38-GCGTCCACCT-ACATGAGTGA_S36_R1_001.fastq.gz', '/data/h_vivian_weiss/DNAseq/1809/1809-AC-38-GCGTCCACCT-ACATGAGTGA_S36_R2_001.fastq.gz' ],
    'P1809_AC_039' => [ '/data/h_vivian_weiss/DNAseq/1809/1809-AC-39-TGCGCATAGC-TAACCGTAAT_S37_R1_001.fastq.gz', '/data/h_vivian_weiss/DNAseq/1809/1809-AC-39-TGCGCATAGC-TAACCGTAAT_S37_R2_001.fastq.gz' ],
    'P1809_AC_040' => [ '/data/h_vivian_weiss/DNAseq/1809/1809-AC-40-CATATGCGAT-CGGTTACGGC_S38_R1_001.fastq.gz', '/data/h_vivian_weiss/DNAseq/1809/1809-AC-40-CATATGCGAT-CGGTTACGGC_S38_R2_001.fastq.gz' ],
    'P1809_AC_041' => [ '/data/h_vivian_weiss/DNAseq/1809/1809-AC-41-TGATGGCTAC-TACCGCCTCG_S39_R1_001.fastq.gz', '/data/h_vivian_weiss/DNAseq/1809/1809-AC-41-TGATGGCTAC-TACCGCCTCG_S39_R2_001.fastq.gz' ],
    'P1809_AC_042' => [ '/data/h_vivian_weiss/DNAseq/1809/1809-AC-42-CAGCAATCGT-CGTTATTCTA_S40_R1_001.fastq.gz', '/data/h_vivian_weiss/DNAseq/1809/1809-AC-42-CAGCAATCGT-CGTTATTCTA_S40_R2_001.fastq.gz' ],
    'P1809_AC_043' => [ '/data/h_vivian_weiss/DNAseq/1809/1809-AC-43-CACAACTTAA-TTCATAAGGT_S41_R1_001.fastq.gz', '/data/h_vivian_weiss/DNAseq/1809/1809-AC-43-CACAACTTAA-TTCATAAGGT_S41_R2_001.fastq.gz' ],
    'P1809_AC_044' => [ '/data/h_vivian_weiss/DNAseq/1809/1809-AC-44-TGTGGTCCGG-CCTGCGGAAC_S42_R1_001.fastq.gz', '/data/h_vivian_weiss/DNAseq/1809/1809-AC-44-TGTGGTCCGG-CCTGCGGAAC_S42_R2_001.fastq.gz' ],
    'P1809_AC_045' => [ '/data/h_vivian_weiss/DNAseq/1809/1809-AC-45-ATACATCACA-AACACTGTTA_S43_R1_001.fastq.gz', '/data/h_vivian_weiss/DNAseq/1809/1809-AC-45-ATACATCACA-AACACTGTTA_S43_R2_001.fastq.gz' ],
    'P1809_AC_046' => [ '/data/h_vivian_weiss/DNAseq/1809/1809-AC-46-GCGTGCTGTG-GGTGTCACCG_S44_R1_001.fastq.gz', '/data/h_vivian_weiss/DNAseq/1809/1809-AC-46-GCGTGCTGTG-GGTGTCACCG_S44_R2_001.fastq.gz' ],
    'P1809_AC_047' => [ '/data/h_vivian_weiss/DNAseq/1809/1809-AC-47-ACCTAAGACC-ACATATCCAG_S45_R1_001.fastq.gz', '/data/h_vivian_weiss/DNAseq/1809/1809-AC-47-ACCTAAGACC-ACATATCCAG_S45_R2_001.fastq.gz' ],
    'P1809_AC_048' => [ '/data/h_vivian_weiss/DNAseq/1809/1809-AC-48-GTTCGGAGTT-GTGCGCTTGA_S46_R1_001.fastq.gz', '/data/h_vivian_weiss/DNAseq/1809/1809-AC-48-GTTCGGAGTT-GTGCGCTTGA_S46_R2_001.fastq.gz' ],
    'P1809_AC_049' => [ '/data/h_vivian_weiss/DNAseq/1809/1809-AC-49-CGAATCTATA-TGGAGTACTT_S47_R1_001.fastq.gz', '/data/h_vivian_weiss/DNAseq/1809/1809-AC-49-CGAATCTATA-TGGAGTACTT_S47_R2_001.fastq.gz' ],
    'P1809_AC_050' => [ '/data/h_vivian_weiss/DNAseq/1809/1809-AC-50-TAGGCTCGCG-CAAGACGTCC_S48_R1_001.fastq.gz', '/data/h_vivian_weiss/DNAseq/1809/1809-AC-50-TAGGCTCGCG-CAAGACGTCC_S48_R2_001.fastq.gz' ],
    'P1809_AC_051' => [ '/data/h_vivian_weiss/DNAseq/1809/1809-AC-51-GGCCTCCAAG-ATTGCGCGGT_S49_R1_001.fastq.gz', '/data/h_vivian_weiss/DNAseq/1809/1809-AC-51-GGCCTCCAAG-ATTGCGCGGT_S49_R2_001.fastq.gz' ],
    'P1809_AC_052' => [ '/data/h_vivian_weiss/DNAseq/1809/1809-AC-52-AATTCTTGGA-GCCATATAAC_S50_R1_001.fastq.gz', '/data/h_vivian_weiss/DNAseq/1809/1809-AC-52-AATTCTTGGA-GCCATATAAC_S50_R2_001.fastq.gz' ],
    'P1809_AC_053' => [ '/data/h_vivian_weiss/DNAseq/1809/1809-AC-53-GCCAATCCTC-CAGTGGCACT_S51_R1_001.fastq.gz', '/data/h_vivian_weiss/DNAseq/1809/1809-AC-53-GCCAATCCTC-CAGTGGCACT_S51_R2_001.fastq.gz' ],
    'P1809_AC_054' => [ '/data/h_vivian_weiss/DNAseq/1809/1809-AC-54-ATTGGCTTCT-TGACAATGTC_S52_R1_001.fastq.gz', '/data/h_vivian_weiss/DNAseq/1809/1809-AC-54-ATTGGCTTCT-TGACAATGTC_S52_R2_001.fastq.gz' ],
    'P1809_AC_055' => [ '/data/h_vivian_weiss/DNAseq/1809/1809-AC-55-TTAGTGAGTC-CGACCTAACG_S53_R1_001.fastq.gz', '/data/h_vivian_weiss/DNAseq/1809/1809-AC-55-TTAGTGAGTC-CGACCTAACG_S53_R2_001.fastq.gz' ],
    'P1809_AC_056' => [ '/data/h_vivian_weiss/DNAseq/1809/1809-AC-56-CCGACAGACT-TAGTTCGGTA_S54_R1_001.fastq.gz', '/data/h_vivian_weiss/DNAseq/1809/1809-AC-56-CCGACAGACT-TAGTTCGGTA_S54_R2_001.fastq.gz' ],
    'P1809_AC_057' => [ '/data/h_vivian_weiss/DNAseq/1809/1809-AC-57-TAGAGAATAC-GCCGCACTCT_S55_R1_001.fastq.gz', '/data/h_vivian_weiss/DNAseq/1809/1809-AC-57-TAGAGAATAC-GCCGCACTCT_S55_R2_001.fastq.gz' ],
    'P1809_AC_058' => [ '/data/h_vivian_weiss/DNAseq/1809/1809-AC-58-CGAGAGGCGT-ATTATGTCTC_S56_R1_001.fastq.gz', '/data/h_vivian_weiss/DNAseq/1809/1809-AC-58-CGAGAGGCGT-ATTATGTCTC_S56_R2_001.fastq.gz' ],
    'P1809_AC_059' => [ '/data/h_vivian_weiss/DNAseq/1809/1809-AC-59-TCGGCGGTTA-AGAACCGAGT_S57_R1_001.fastq.gz', '/data/h_vivian_weiss/DNAseq/1809/1809-AC-59-TCGGCGGTTA-AGAACCGAGT_S57_R2_001.fastq.gz' ],
    'P1809_AC_060' => [ '/data/h_vivian_weiss/DNAseq/1809/1809-AC-60-CTAATAACCG-GAGGTTAGAC_S58_R1_001.fastq.gz', '/data/h_vivian_weiss/DNAseq/1809/1809-AC-60-CTAATAACCG-GAGGTTAGAC_S58_R2_001.fastq.gz' ],
    'P1809_AC_061' => [ '/data/h_vivian_weiss/DNAseq/1809/1809-AC-61-AGGCGTTCGC-GTTAATTACG_S59_R1_001.fastq.gz', '/data/h_vivian_weiss/DNAseq/1809/1809-AC-61-AGGCGTTCGC-GTTAATTACG_S59_R2_001.fastq.gz' ],
    'P1809_AC_062' => [ '/data/h_vivian_weiss/DNAseq/1809/1809-AC-62-GAATACCTAT-ACCGGCCGTA_S60_R1_001.fastq.gz', '/data/h_vivian_weiss/DNAseq/1809/1809-AC-62-GAATACCTAT-ACCGGCCGTA_S60_R2_001.fastq.gz' ],
    'P1809_AC_063' => [ '/data/h_vivian_weiss/DNAseq/1809/1809-AC-63-CCACGCTGAA-CCTCGTGCGT_S61_R1_001.fastq.gz', '/data/h_vivian_weiss/DNAseq/1809/1809-AC-63-CCACGCTGAA-CCTCGTGCGT_S61_R2_001.fastq.gz' ],
    'P1809_AC_064' => [ '/data/h_vivian_weiss/DNAseq/1809/1809-AC-64-TTGTATCAGG-TTCTACATAC_S62_R1_001.fastq.gz', '/data/h_vivian_weiss/DNAseq/1809/1809-AC-64-TTGTATCAGG-TTCTACATAC_S62_R2_001.fastq.gz' ],
    'P1809_AC_065' => [ '/data/h_vivian_weiss/DNAseq/1809/1809-AC-65-GGACCAGTGG-ACTTCAAGCG_S63_R1_001.fastq.gz', '/data/h_vivian_weiss/DNAseq/1809/1809-AC-65-GGACCAGTGG-ACTTCAAGCG_S63_R2_001.fastq.gz' ],
    'P1809_AC_066' => [ '/data/h_vivian_weiss/DNAseq/1809/1809-AC-66-AAGTTGACAA-GTCCTGGATA_S64_R1_001.fastq.gz', '/data/h_vivian_weiss/DNAseq/1809/1809-AC-66-AAGTTGACAA-GTCCTGGATA_S64_R2_001.fastq.gz' ],
    'P1809_AC_067' => [ '/data/h_vivian_weiss/DNAseq/1809/1809-AC-67-GCTTGAACGC-AATCACCAGC_S65_R1_001.fastq.gz', '/data/h_vivian_weiss/DNAseq/1809/1809-AC-67-GCTTGAACGC-AATCACCAGC_S65_R2_001.fastq.gz' ],
    'P1809_AC_068' => [ '/data/h_vivian_weiss/DNAseq/1809/1809-AC-68-ATCCAGGTAT-GGCTGTTGAT_S66_R1_001.fastq.gz', '/data/h_vivian_weiss/DNAseq/1809/1809-AC-68-ATCCAGGTAT-GGCTGTTGAT_S66_R2_001.fastq.gz' ],
    'P1809_AC_069' => [ '/data/h_vivian_weiss/DNAseq/1809/1809-AC-69-GTCTCGCCAC-CAATCGGCTG_S67_R1_001.fastq.gz', '/data/h_vivian_weiss/DNAseq/1809/1809-AC-69-GTCTCGCCAC-CAATCGGCTG_S67_R2_001.fastq.gz' ],
    'P1809_AC_070' => [ '/data/h_vivian_weiss/DNAseq/1809/1809-AC-70-ACTCTATTGT-TGGCTAATCA_S68_R1_001.fastq.gz', '/data/h_vivian_weiss/DNAseq/1809/1809-AC-70-ACTCTATTGT-TGGCTAATCA_S68_R2_001.fastq.gz' ],
    'P1809_AC_071' => [ '/data/h_vivian_weiss/DNAseq/1809/1809-AC-71-CTTATTGGCC-TTAAGACAAG_S69_R1_001.fastq.gz', '/data/h_vivian_weiss/DNAseq/1809/1809-AC-71-CTTATTGGCC-TTAAGACAAG_S69_R2_001.fastq.gz' ],
    'P1809_AC_072' => [ '/data/h_vivian_weiss/DNAseq/1809/1809-AC-72-TCCGCCAATT-CCGGAGTGGA_S70_R1_001.fastq.gz', '/data/h_vivian_weiss/DNAseq/1809/1809-AC-72-TCCGCCAATT-CCGGAGTGGA_S70_R2_001.fastq.gz' ],
    'P1809_AC_073' => [ '/data/h_vivian_weiss/DNAseq/1809/1809-AC-73-TCAATGGAGA-AACTTATCCT_S71_R1_001.fastq.gz', '/data/h_vivian_weiss/DNAseq/1809/1809-AC-73-TCAATGGAGA-AACTTATCCT_S71_R2_001.fastq.gz' ],
    'P1809_AC_074' => [ '/data/h_vivian_weiss/DNAseq/1809/1809-AC-74-CTGGCAAGAG-GGTCCGCTTC_S72_R1_001.fastq.gz', '/data/h_vivian_weiss/DNAseq/1809/1809-AC-74-CTGGCAAGAG-GGTCCGCTTC_S72_R2_001.fastq.gz' ],
    'P1809_AC_075' => [ '/data/h_vivian_weiss/DNAseq/1809/1809-AC-75-AGAGAACCTA-CTCGCTTCGG_S73_R1_001.fastq.gz', '/data/h_vivian_weiss/DNAseq/1809/1809-AC-75-AGAGAACCTA-CTCGCTTCGG_S73_R2_001.fastq.gz' ],
    'P1809_AC_076' => [ '/data/h_vivian_weiss/DNAseq/1809/1809-AC-76-GAGAGGTTCG-TCTATCCTAA_S74_R1_001.fastq.gz', '/data/h_vivian_weiss/DNAseq/1809/1809-AC-76-GAGAGGTTCG-TCTATCCTAA_S74_R2_001.fastq.gz' ],
    'P1809_AC_077' => [ '/data/h_vivian_weiss/DNAseq/1809/1809-AC-77-ACTTGTTATC-GAACTTCCTT_S75_R1_001.fastq.gz', '/data/h_vivian_weiss/DNAseq/1809/1809-AC-77-ACTTGTTATC-GAACTTCCTT_S75_R2_001.fastq.gz' ],
    'P1809_AC_078' => [ '/data/h_vivian_weiss/DNAseq/1809/1809-AC-78-GTCCACCGCT-AGGTCCTTCC_S76_R1_001.fastq.gz', '/data/h_vivian_weiss/DNAseq/1809/1809-AC-78-GTCCACCGCT-AGGTCCTTCC_S76_R2_001.fastq.gz' ],
    'P1809_AC_079' => [ '/data/h_vivian_weiss/DNAseq/1809/1809-AC-79-TAGCATAACC-TATGGAGATT_S77_R1_001.fastq.gz', '/data/h_vivian_weiss/DNAseq/1809/1809-AC-79-TAGCATAACC-TATGGAGATT_S77_R2_001.fastq.gz' ],
    'P1809_AC_080' => [ '/data/h_vivian_weiss/DNAseq/1809/1809-AC-80-CGATGCGGTT-CGCAAGAGCC_S78_R1_001.fastq.gz', '/data/h_vivian_weiss/DNAseq/1809/1809-AC-80-CGATGCGGTT-CGCAAGAGCC_S78_R2_001.fastq.gz' ],
    'P1809_AC_081' => [ '/data/h_vivian_weiss/DNAseq/1809/1809-AC-81-CGCTGTCTCA-ACCACGACAT_S79_R1_001.fastq.gz', '/data/h_vivian_weiss/DNAseq/1809/1809-AC-81-CGCTGTCTCA-ACCACGACAT_S79_R2_001.fastq.gz' ],
    'P1809_AC_082' => [ '/data/h_vivian_weiss/DNAseq/1809/1809-AC-82-TATCACTCTG-GTTGTAGTGC_S80_R1_001.fastq.gz', '/data/h_vivian_weiss/DNAseq/1809/1809-AC-82-TATCACTCTG-GTTGTAGTGC_S80_R2_001.fastq.gz' ],
    'P1809_AC_083' => [ '/data/h_vivian_weiss/DNAseq/1809/1809-AC-83-GTTGGATGAA-TCATAGATTG_S81_R1_001.fastq.gz', '/data/h_vivian_weiss/DNAseq/1809/1809-AC-83-GTTGGATGAA-TCATAGATTG_S81_R2_001.fastq.gz' ],
    'P1809_AC_084' => [ '/data/h_vivian_weiss/DNAseq/1809/1809-AC-84-ACCAAGCAGG-CTGCGAGCCA_S82_R1_001.fastq.gz', '/data/h_vivian_weiss/DNAseq/1809/1809-AC-84-ACCAAGCAGG-CTGCGAGCCA_S82_R2_001.fastq.gz' ],
    'P1809_AC_085' => [ '/data/h_vivian_weiss/DNAseq/1809/1809-AC-85-GCACCACCAA-CAGCACGGAG_S83_R1_001.fastq.gz', '/data/h_vivian_weiss/DNAseq/1809/1809-AC-85-GCACCACCAA-CAGCACGGAG_S83_R2_001.fastq.gz' ],
    'P1809_AC_086' => [ '/data/h_vivian_weiss/DNAseq/1809/1809-AC-86-ATGTTGTTGG-TGATGTAAGA_S84_R1_001.fastq.gz', '/data/h_vivian_weiss/DNAseq/1809/1809-AC-86-ATGTTGTTGG-TGATGTAAGA_S84_R2_001.fastq.gz' ],
    'P1809_AC_087' => [ '/data/h_vivian_weiss/DNAseq/1809/1809-AC-87-AGCCTCAGGC-CCGTTCAAGG_S85_R1_001.fastq.gz', '/data/h_vivian_weiss/DNAseq/1809/1809-AC-87-AGCCTCAGGC-CCGTTCAAGG_S85_R2_001.fastq.gz' ],
    'P1809_AC_088' => [ '/data/h_vivian_weiss/DNAseq/1809/1809-AC-88-GATTCTGAAT-TTACCTGGAA_S86_R1_001.fastq.gz', '/data/h_vivian_weiss/DNAseq/1809/1809-AC-88-GATTCTGAAT-TTACCTGGAA_S86_R2_001.fastq.gz' ],
    'P1809_AC_089' => [ '/data/h_vivian_weiss/DNAseq/1809/1809-AC-89-CAATTCTCAC-AAGCATCTTG_S87_R1_001.fastq.gz', '/data/h_vivian_weiss/DNAseq/1809/1809-AC-89-CAATTCTCAC-AAGCATCTTG_S87_R2_001.fastq.gz' ],
    'P1809_AC_090' => [ '/data/h_vivian_weiss/DNAseq/1809/1809-AC-90-TGGCCTCTGT-GGATGCTCCA_S88_R1_001.fastq.gz', '/data/h_vivian_weiss/DNAseq/1809/1809-AC-90-TGGCCTCTGT-GGATGCTCCA_S88_R2_001.fastq.gz' ],
    'P1809_AC_091' => [ '/data/h_vivian_weiss/DNAseq/1809/1809-AC-91-GGCATAGGTG-AACAGACGGC_S89_R1_001.fastq.gz', '/data/h_vivian_weiss/DNAseq/1809/1809-AC-91-GGCATAGGTG-AACAGACGGC_S89_R2_001.fastq.gz' ],
    'P1809_AC_092' => [ '/data/h_vivian_weiss/DNAseq/1809/1809-AC-92-AATGCGAACA-GGTGAGTAAT_S90_R1_001.fastq.gz', '/data/h_vivian_weiss/DNAseq/1809/1809-AC-92-AATGCGAACA-GGTGAGTAAT_S90_R2_001.fastq.gz' ],
    'P1809_AC_093' => [ '/data/h_vivian_weiss/DNAseq/1809/1809-AC-93-GCTAATAGGA-CTAGTGCTCT_S91_R1_001.fastq.gz', '/data/h_vivian_weiss/DNAseq/1809/1809-AC-93-GCTAATAGGA-CTAGTGCTCT_S91_R2_001.fastq.gz' ],
    'P1809_AC_094' => [ '/data/h_vivian_weiss/DNAseq/1809/1809-AC-94-ATCGGCGAAG-TCGACATCTC_S92_R1_001.fastq.gz', '/data/h_vivian_weiss/DNAseq/1809/1809-AC-94-ATCGGCGAAG-TCGACATCTC_S92_R2_001.fastq.gz' ],
    'P1809_AC_095' => [ '/data/h_vivian_weiss/DNAseq/1809/1809-AC-95-TGGACCGCCA-CGTCTCATAT_S93_R1_001.fastq.gz', '/data/h_vivian_weiss/DNAseq/1809/1809-AC-95-TGGACCGCCA-CGTCTCATAT_S93_R2_001.fastq.gz' ],
    'P1809_AC_096' => [ '/data/h_vivian_weiss/DNAseq/1809/1809-AC-96-CAAGTTATTG-TACTCTGCGC_S94_R1_001.fastq.gz', '/data/h_vivian_weiss/DNAseq/1809/1809-AC-96-CAAGTTATTG-TACTCTGCGC_S94_R2_001.fastq.gz' ],
    'P1809_AC_097' => [ '/data/h_vivian_weiss/DNAseq/1809/1809-AC-97-CTCGAATATA-TCAGAAGGCG_S95_R1_001.fastq.gz', '/data/h_vivian_weiss/DNAseq/1809/1809-AC-97-CTCGAATATA-TCAGAAGGCG_S95_R2_001.fastq.gz' ],
    'P1809_AC_098' => [ '/data/h_vivian_weiss/DNAseq/1809/1809-AC-98-TCTAGGCGCG-CTGAGGAATA_S96_R1_001.fastq.gz', '/data/h_vivian_weiss/DNAseq/1809/1809-AC-98-TCTAGGCGCG-CTGAGGAATA_S96_R2_001.fastq.gz' ],
    'P1809_AC_099' => [ '/data/h_vivian_weiss/DNAseq/1809/1809-AC-99-GCACGGTACC-GACACCATGT_S97_R1_001.fastq.gz', '/data/h_vivian_weiss/DNAseq/1809/1809-AC-99-GCACGGTACC-GACACCATGT_S97_R2_001.fastq.gz' ],
    'P1809_AC_100' => [ '/data/h_vivian_weiss/DNAseq/1809/1809-AC-100-ATGTAACGTT-AGTGTTGCAC_S98_R1_001.fastq.gz', '/data/h_vivian_weiss/DNAseq/1809/1809-AC-100-ATGTAACGTT-AGTGTTGCAC_S98_R2_001.fastq.gz' ],
    'P1809_AC_101' => [ '/data/h_vivian_weiss/DNAseq/1809/1809-AC-101-AGTGTTGCAC-ATGTAACGTT_S99_R1_001.fastq.gz', '/data/h_vivian_weiss/DNAseq/1809/1809-AC-101-AGTGTTGCAC-ATGTAACGTT_S99_R2_001.fastq.gz' ],
    'P1809_AC_102' => [ '/data/h_vivian_weiss/DNAseq/1809/1809-AC-102-GACACCATGT-GCACGGTACC_S100_R1_001.fastq.gz', '/data/h_vivian_weiss/DNAseq/1809/1809-AC-102-GACACCATGT-GCACGGTACC_S100_R2_001.fastq.gz' ],
    'P1809_AC_103' => [ '/data/h_vivian_weiss/DNAseq/1809/1809-AC-103-CTGAGGAATA-TCTAGGCGCG_S101_R1_001.fastq.gz', '/data/h_vivian_weiss/DNAseq/1809/1809-AC-103-CTGAGGAATA-TCTAGGCGCG_S101_R2_001.fastq.gz' ],
    'P1809_AC_104' => [ '/data/h_vivian_weiss/DNAseq/1809/1809-AC-104-TCAGAAGGCG-CTCGAATATA_S102_R1_001.fastq.gz', '/data/h_vivian_weiss/DNAseq/1809/1809-AC-104-TCAGAAGGCG-CTCGAATATA_S102_R2_001.fastq.gz' ],
    'P1809_AC_105' => [ '/data/h_vivian_weiss/DNAseq/1809/1809-AC-105-TACTCTGCGC-CAAGTTATTG_S103_R1_001.fastq.gz', '/data/h_vivian_weiss/DNAseq/1809/1809-AC-105-TACTCTGCGC-CAAGTTATTG_S103_R2_001.fastq.gz' ],
    'P1809_AC_106' => [ '/data/h_vivian_weiss/DNAseq/1809/1809-AC-106-CGTCTCATAT-TGGACCGCCA_S104_R1_001.fastq.gz', '/data/h_vivian_weiss/DNAseq/1809/1809-AC-106-CGTCTCATAT-TGGACCGCCA_S104_R2_001.fastq.gz' ],
    'P1809_AC_107' => [ '/data/h_vivian_weiss/DNAseq/1809/1809-AC-107-TCGACATCTC-ATCGGCGAAG_S105_R1_001.fastq.gz', '/data/h_vivian_weiss/DNAseq/1809/1809-AC-107-TCGACATCTC-ATCGGCGAAG_S105_R2_001.fastq.gz' ],
    'P1809_AC_108' => [ '/data/h_vivian_weiss/DNAseq/1809/1809-AC-108-CTAGTGCTCT-GCTAATAGGA_S106_R1_001.fastq.gz', '/data/h_vivian_weiss/DNAseq/1809/1809-AC-108-CTAGTGCTCT-GCTAATAGGA_S106_R2_001.fastq.gz' ],
    'P1809_AC_109' => [ '/data/h_vivian_weiss/DNAseq/1809/1809-AC-109-GGTGAGTAAT-AATGCGAACA_S107_R1_001.fastq.gz', '/data/h_vivian_weiss/DNAseq/1809/1809-AC-109-GGTGAGTAAT-AATGCGAACA_S107_R2_001.fastq.gz' ],
    'P1809_AC_110' => [ '/data/h_vivian_weiss/DNAseq/1809/1809-AC-110-AACAGACGGC-GGCATAGGTG_S108_R1_001.fastq.gz', '/data/h_vivian_weiss/DNAseq/1809/1809-AC-110-AACAGACGGC-GGCATAGGTG_S108_R2_001.fastq.gz' ],
    'P1809_AC_111' => [ '/data/h_vivian_weiss/DNAseq/1809/1809-AC-111-GGATGCTCCA-TGGCCTCTGT_S109_R1_001.fastq.gz', '/data/h_vivian_weiss/DNAseq/1809/1809-AC-111-GGATGCTCCA-TGGCCTCTGT_S109_R2_001.fastq.gz' ],
    'P1809_AC_112' => [ '/data/h_vivian_weiss/DNAseq/1809/1809-AC-112-AAGCATCTTG-CAATTCTCAC_S110_R1_001.fastq.gz', '/data/h_vivian_weiss/DNAseq/1809/1809-AC-112-AAGCATCTTG-CAATTCTCAC_S110_R2_001.fastq.gz' ],
    'P1809_AC_113' => [ '/data/h_vivian_weiss/DNAseq/1809/1809-AC-113-TTACCTGGAA-GATTCTGAAT_S111_R1_001.fastq.gz', '/data/h_vivian_weiss/DNAseq/1809/1809-AC-113-TTACCTGGAA-GATTCTGAAT_S111_R2_001.fastq.gz' ],
    'P1809_AC_114' => [ '/data/h_vivian_weiss/DNAseq/1809/1809-AC-114-CCGTTCAAGG-AGCCTCAGGC_S112_R1_001.fastq.gz', '/data/h_vivian_weiss/DNAseq/1809/1809-AC-114-CCGTTCAAGG-AGCCTCAGGC_S112_R2_001.fastq.gz' ],
    'P1809_AC_115' => [ '/data/h_vivian_weiss/DNAseq/1809/1809-AC-115-TGATGTAAGA-ATGTTGTTGG_S113_R1_001.fastq.gz', '/data/h_vivian_weiss/DNAseq/1809/1809-AC-115-TGATGTAAGA-ATGTTGTTGG_S113_R2_001.fastq.gz' ],
    'P1809_AC_116' => [ '/data/h_vivian_weiss/DNAseq/1809/1809-AC-116-CAGCACGGAG-GCACCACCAA_S114_R1_001.fastq.gz', '/data/h_vivian_weiss/DNAseq/1809/1809-AC-116-CAGCACGGAG-GCACCACCAA_S114_R2_001.fastq.gz' ],
    'P1809_AC_117' => [ '/data/h_vivian_weiss/DNAseq/1809/1809-AC-117-CTGCGAGCCA-ACCAAGCAGG_S115_R1_001.fastq.gz', '/data/h_vivian_weiss/DNAseq/1809/1809-AC-117-CTGCGAGCCA-ACCAAGCAGG_S115_R2_001.fastq.gz' ],
    'P1809_AC_118' => [ '/data/h_vivian_weiss/DNAseq/1809/1809-AC-118-TCATAGATTG-GTTGGATGAA_S116_R1_001.fastq.gz', '/data/h_vivian_weiss/DNAseq/1809/1809-AC-118-TCATAGATTG-GTTGGATGAA_S116_R2_001.fastq.gz' ],
    'P1809_AC_119' => [ '/data/h_vivian_weiss/DNAseq/1809/1809-AC-119-GTTGTAGTGC-TATCACTCTG_S117_R1_001.fastq.gz', '/data/h_vivian_weiss/DNAseq/1809/1809-AC-119-GTTGTAGTGC-TATCACTCTG_S117_R2_001.fastq.gz' ],
    'P1809_AC_120' => [ '/data/h_vivian_weiss/DNAseq/1809/1809-AC-120-ACCACGACAT-CGCTGTCTCA_S118_R1_001.fastq.gz', '/data/h_vivian_weiss/DNAseq/1809/1809-AC-120-ACCACGACAT-CGCTGTCTCA_S118_R2_001.fastq.gz' ],
    'P1809_AC_121' => [ '/data/h_vivian_weiss/DNAseq/1809/1809-AC-121-CGCAAGAGCC-CGATGCGGTT_S119_R1_001.fastq.gz', '/data/h_vivian_weiss/DNAseq/1809/1809-AC-121-CGCAAGAGCC-CGATGCGGTT_S119_R2_001.fastq.gz' ],
    'P1809_AC_122' => [ '/data/h_vivian_weiss/DNAseq/1809/1809-AC-122-TATGGAGATT-TAGCATAACC_S120_R1_001.fastq.gz', '/data/h_vivian_weiss/DNAseq/1809/1809-AC-122-TATGGAGATT-TAGCATAACC_S120_R2_001.fastq.gz' ],
    'P1809_AC_123' => [ '/data/h_vivian_weiss/DNAseq/1809/1809-AC-123-AGGTCCTTCC-GTCCACCGCT_S121_R1_001.fastq.gz', '/data/h_vivian_weiss/DNAseq/1809/1809-AC-123-AGGTCCTTCC-GTCCACCGCT_S121_R2_001.fastq.gz' ],
    'P1809_AC_124' => [ '/data/h_vivian_weiss/DNAseq/1809/1809-AC-124-GAACTTCCTT-ACTTGTTATC_S122_R1_001.fastq.gz', '/data/h_vivian_weiss/DNAseq/1809/1809-AC-124-GAACTTCCTT-ACTTGTTATC_S122_R2_001.fastq.gz' ],
    'P1809_AC_125' => [ '/data/h_vivian_weiss/DNAseq/1809/1809-AC-125-TCTATCCTAA-GAGAGGTTCG_S123_R1_001.fastq.gz', '/data/h_vivian_weiss/DNAseq/1809/1809-AC-125-TCTATCCTAA-GAGAGGTTCG_S123_R2_001.fastq.gz' ],
    'P1809_AC_126' => [ '/data/h_vivian_weiss/DNAseq/1809/1809-AC-126-CTCGCTTCGG-AGAGAACCTA_S124_R1_001.fastq.gz', '/data/h_vivian_weiss/DNAseq/1809/1809-AC-126-CTCGCTTCGG-AGAGAACCTA_S124_R2_001.fastq.gz' ],
    'P1809_AC_127' => [ '/data/h_vivian_weiss/DNAseq/1809/1809-AC-127-GGTCCGCTTC-CTGGCAAGAG_S125_R1_001.fastq.gz', '/data/h_vivian_weiss/DNAseq/1809/1809-AC-127-GGTCCGCTTC-CTGGCAAGAG_S125_R2_001.fastq.gz' ],
    'P1809_AC_128' => [ '/data/h_vivian_weiss/DNAseq/1809/1809-AC-128-AACTTATCCT-TCAATGGAGA_S126_R1_001.fastq.gz', '/data/h_vivian_weiss/DNAseq/1809/1809-AC-128-AACTTATCCT-TCAATGGAGA_S126_R2_001.fastq.gz' ],
    'P1809_AC_129' => [ '/data/h_vivian_weiss/DNAseq/1809/1809-AC-129-CCGGAGTGGA-TCCGCCAATT_S127_R1_001.fastq.gz', '/data/h_vivian_weiss/DNAseq/1809/1809-AC-129-CCGGAGTGGA-TCCGCCAATT_S127_R2_001.fastq.gz' ],
    'P1809_AC_130' => [ '/data/h_vivian_weiss/DNAseq/1809/1809-AC-130-TTAAGACAAG-CTTATTGGCC_S128_R1_001.fastq.gz', '/data/h_vivian_weiss/DNAseq/1809/1809-AC-130-TTAAGACAAG-CTTATTGGCC_S128_R2_001.fastq.gz' ],
    'P1809_AC_131' => [ '/data/h_vivian_weiss/DNAseq/1809/1809-AC-131-TGGCTAATCA-ACTCTATTGT_S129_R1_001.fastq.gz', '/data/h_vivian_weiss/DNAseq/1809/1809-AC-131-TGGCTAATCA-ACTCTATTGT_S129_R2_001.fastq.gz' ],
    'P1809_AC_132' => [ '/data/h_vivian_weiss/DNAseq/1809/1809-AC-132-CAATCGGCTG-GTCTCGCCAC_S130_R1_001.fastq.gz', '/data/h_vivian_weiss/DNAseq/1809/1809-AC-132-CAATCGGCTG-GTCTCGCCAC_S130_R2_001.fastq.gz' ],
    'P1809_AC_133' => [ '/data/h_vivian_weiss/DNAseq/1809/1809-AC-133-GGCTGTTGAT-ATCCAGGTAT_S131_R1_001.fastq.gz', '/data/h_vivian_weiss/DNAseq/1809/1809-AC-133-GGCTGTTGAT-ATCCAGGTAT_S131_R2_001.fastq.gz' ],
    'P1809_AC_134' => [ '/data/h_vivian_weiss/DNAseq/1809/1809-AC-134-AATCACCAGC-GCTTGAACGC_S132_R1_001.fastq.gz', '/data/h_vivian_weiss/DNAseq/1809/1809-AC-134-AATCACCAGC-GCTTGAACGC_S132_R2_001.fastq.gz' ],
    'P1809_AC_135' => [ '/data/h_vivian_weiss/DNAseq/1809/1809-AC-135-GTCCTGGATA-AAGTTGACAA_S133_R1_001.fastq.gz', '/data/h_vivian_weiss/DNAseq/1809/1809-AC-135-GTCCTGGATA-AAGTTGACAA_S133_R2_001.fastq.gz' ],
    #file_def.py -i /data/h_vivian_weiss/DNAseq/1896/ -n "(........\\d+)" -f "*.gz" -a -p
    'P1896_AC_01' => [ '/data/h_vivian_weiss/DNAseq/1896/1896-AC-1-GGATATATCC-TGTTCACCAT_S44_R1_001.fastq.gz', '/data/h_vivian_weiss/DNAseq/1896/1896-AC-1-GGATATATCC-TGTTCACCAT_S44_R2_001.fastq.gz' ],
    'P1896_AC_02' => [ '/data/h_vivian_weiss/DNAseq/1896/1896-AC-2-TTATACGCGA-AGCTGGAATG_S45_R1_001.fastq.gz', '/data/h_vivian_weiss/DNAseq/1896/1896-AC-2-TTATACGCGA-AGCTGGAATG_S45_R2_001.fastq.gz' ],
    'P1896_AC_03' => [ '/data/h_vivian_weiss/DNAseq/1896/1896-AC-3-GAACGCAATA-CAATATAGGT_S46_R1_001.fastq.gz', '/data/h_vivian_weiss/DNAseq/1896/1896-AC-3-GAACGCAATA-CAATATAGGT_S46_R2_001.fastq.gz' ],
    'P1896_AC_04' => [ '/data/h_vivian_weiss/DNAseq/1896/1896-AC-4-TGATGGCTAC-TACCGCCTCG_S47_R1_001.fastq.gz', '/data/h_vivian_weiss/DNAseq/1896/1896-AC-4-TGATGGCTAC-TACCGCCTCG_S47_R2_001.fastq.gz' ],
    'P1896_AC_05' => [ '/data/h_vivian_weiss/DNAseq/1896/1896-AC-5-GCCAATCCTC-CAGTGGCACT_S48_R1_001.fastq.gz', '/data/h_vivian_weiss/DNAseq/1896/1896-AC-5-GCCAATCCTC-CAGTGGCACT_S48_R2_001.fastq.gz' ],
    'P1896_AC_06' => [ '/data/h_vivian_weiss/DNAseq/1896/1896-AC-6-GGACCAGTGG-ACTTCAAGCG_S49_R1_001.fastq.gz', '/data/h_vivian_weiss/DNAseq/1896/1896-AC-6-GGACCAGTGG-ACTTCAAGCG_S49_R2_001.fastq.gz' ],
    'P1896_AC_07' => [ '/data/h_vivian_weiss/DNAseq/1896/1896-AC-7-ACTTGTTATC-GAACTTCCTT_S50_R1_001.fastq.gz', '/data/h_vivian_weiss/DNAseq/1896/1896-AC-7-ACTTGTTATC-GAACTTCCTT_S50_R2_001.fastq.gz' ],
    'P1896_AC_08' => [ '/data/h_vivian_weiss/DNAseq/1896/1896-AC-8-CAATTCTCAC-AAGCATCTTG_S51_R1_001.fastq.gz', '/data/h_vivian_weiss/DNAseq/1896/1896-AC-8-CAATTCTCAC-AAGCATCTTG_S51_R2_001.fastq.gz' ],
    'P1896_AC_09' => [ '/data/h_vivian_weiss/DNAseq/1896/1896-AC-9-AAGCGCGCTT-CACCTGTTGC_S52_R1_001.fastq.gz', '/data/h_vivian_weiss/DNAseq/1896/1896-AC-9-AAGCGCGCTT-CACCTGTTGC_S52_R2_001.fastq.gz' ],
    'P1896_AC_10' => [ '/data/h_vivian_weiss/DNAseq/1896/1896-AC-10-CCGCGTATAG-GATCAAGGCA_S53_R1_001.fastq.gz', '/data/h_vivian_weiss/DNAseq/1896/1896-AC-10-CCGCGTATAG-GATCAAGGCA_S53_R2_001.fastq.gz' ],
    'P1896_AC_11' => [ '/data/h_vivian_weiss/DNAseq/1896/1896-AC-11-AGGTATGGCG-TGGCGCGAAC_S54_R1_001.fastq.gz', '/data/h_vivian_weiss/DNAseq/1896/1896-AC-11-AGGTATGGCG-TGGCGCGAAC_S54_R2_001.fastq.gz' ],
    'P1896_AC_12' => [ '/data/h_vivian_weiss/DNAseq/1896/1896-AC-12-CAGCAATCGT-CGTTATTCTA_S55_R1_001.fastq.gz', '/data/h_vivian_weiss/DNAseq/1896/1896-AC-12-CAGCAATCGT-CGTTATTCTA_S55_R2_001.fastq.gz' ],
    'P1896_AC_13' => [ '/data/h_vivian_weiss/DNAseq/1896/1896-AC-13-ATTGGCTTCT-TGACAATGTC_S56_R1_001.fastq.gz', '/data/h_vivian_weiss/DNAseq/1896/1896-AC-13-ATTGGCTTCT-TGACAATGTC_S56_R2_001.fastq.gz' ],
    'P1896_AC_14' => [ '/data/h_vivian_weiss/DNAseq/1896/1896-AC-14-AAGTTGACAA-GTCCTGGATA_S57_R1_001.fastq.gz', '/data/h_vivian_weiss/DNAseq/1896/1896-AC-14-AAGTTGACAA-GTCCTGGATA_S57_R2_001.fastq.gz' ],
    'P1896_AC_15' => [ '/data/h_vivian_weiss/DNAseq/1896/1896-AC-15-GTCCACCGCT-AGGTCCTTCC_S58_R1_001.fastq.gz', '/data/h_vivian_weiss/DNAseq/1896/1896-AC-15-GTCCACCGCT-AGGTCCTTCC_S58_R2_001.fastq.gz' ],
    'P1896_AC_16' => [ '/data/h_vivian_weiss/DNAseq/1896/1896-AC-16-TGGCCTCTGT-GGATGCTCCA_S59_R1_001.fastq.gz', '/data/h_vivian_weiss/DNAseq/1896/1896-AC-16-TGGCCTCTGT-GGATGCTCCA_S59_R2_001.fastq.gz' ],
    'P1896_AC_17' => [ '/data/h_vivian_weiss/DNAseq/1896/1896-AC-17-CAACGAGAGC-GTAGGTGGTG_S60_R1_001.fastq.gz', '/data/h_vivian_weiss/DNAseq/1896/1896-AC-17-CAACGAGAGC-GTAGGTGGTG_S60_R2_001.fastq.gz' ],
    'P1896_AC_18' => [ '/data/h_vivian_weiss/DNAseq/1896/1896-AC-18-GACGTCTGCA-TCAGTCTCGT_S61_R1_001.fastq.gz', '/data/h_vivian_weiss/DNAseq/1896/1896-AC-18-GACGTCTGCA-TCAGTCTCGT_S61_R2_001.fastq.gz' ],
    'P1896_AC_19' => [ '/data/h_vivian_weiss/DNAseq/1896/1896-AC-19-CGTACAGGAA-ATGGTTGACT_S62_R1_001.fastq.gz', '/data/h_vivian_weiss/DNAseq/1896/1896-AC-19-CGTACAGGAA-ATGGTTGACT_S62_R2_001.fastq.gz' ],
    'P1896_AC_20' => [ '/data/h_vivian_weiss/DNAseq/1896/1896-AC-20-CACAACTTAA-TTCATAAGGT_S63_R1_001.fastq.gz', '/data/h_vivian_weiss/DNAseq/1896/1896-AC-20-CACAACTTAA-TTCATAAGGT_S63_R2_001.fastq.gz' ],
    'P1896_AC_21' => [ '/data/h_vivian_weiss/DNAseq/1896/1896-AC-21-TTAGTGAGTC-CGACCTAACG_S64_R1_001.fastq.gz', '/data/h_vivian_weiss/DNAseq/1896/1896-AC-21-TTAGTGAGTC-CGACCTAACG_S64_R2_001.fastq.gz' ],
    'P1896_AC_22' => [ '/data/h_vivian_weiss/DNAseq/1896/1896-AC-22-GCTTGAACGC-AATCACCAGC_S65_R1_001.fastq.gz', '/data/h_vivian_weiss/DNAseq/1896/1896-AC-22-GCTTGAACGC-AATCACCAGC_S65_R2_001.fastq.gz' ],
    'P1896_AC_23' => [ '/data/h_vivian_weiss/DNAseq/1896/1896-AC-23-TAGCATAACC-TATGGAGATT_S66_R1_001.fastq.gz', '/data/h_vivian_weiss/DNAseq/1896/1896-AC-23-TAGCATAACC-TATGGAGATT_S66_R2_001.fastq.gz' ],
    'P1896_AC_24' => [ '/data/h_vivian_weiss/DNAseq/1896/1896-AC-24-GGCATAGGTG-AACAGACGGC_S67_R1_001.fastq.gz', '/data/h_vivian_weiss/DNAseq/1896/1896-AC-24-GGCATAGGTG-AACAGACGGC_S67_R2_001.fastq.gz' ],
    'P1896_AC_25' => [ '/data/h_vivian_weiss/DNAseq/1896/1896-AC-25-TGGTAGAGAT-ACGAACAACA_S68_R1_001.fastq.gz', '/data/h_vivian_weiss/DNAseq/1896/1896-AC-25-TGGTAGAGAT-ACGAACAACA_S68_R2_001.fastq.gz' ],
    'P1896_AC_26' => [ '/data/h_vivian_weiss/DNAseq/1896/1896-AC-26-AGTACTCATG-CTGACTCTAC_S69_R1_001.fastq.gz', '/data/h_vivian_weiss/DNAseq/1896/1896-AC-26-AGTACTCATG-CTGACTCTAC_S69_R2_001.fastq.gz' ],
    'P1896_AC_27' => [ '/data/h_vivian_weiss/DNAseq/1896/1896-AC-27-TACGTGAAGG-GCAACCAGTC_S70_R1_001.fastq.gz', '/data/h_vivian_weiss/DNAseq/1896/1896-AC-27-TACGTGAAGG-GCAACCAGTC_S70_R2_001.fastq.gz' ],
    'P1896_AC_28' => [ '/data/h_vivian_weiss/DNAseq/1896/1896-AC-28-TGTGGTCCGG-CCTGCGGAAC_S71_R1_001.fastq.gz', '/data/h_vivian_weiss/DNAseq/1896/1896-AC-28-TGTGGTCCGG-CCTGCGGAAC_S71_R2_001.fastq.gz' ],
    'P1896_AC_29' => [ '/data/h_vivian_weiss/DNAseq/1896/1896-AC-29-CCGACAGACT-TAGTTCGGTA_S72_R1_001.fastq.gz', '/data/h_vivian_weiss/DNAseq/1896/1896-AC-29-CCGACAGACT-TAGTTCGGTA_S72_R2_001.fastq.gz' ],
    'P1896_AC_30' => [ '/data/h_vivian_weiss/DNAseq/1896/1896-AC-30-ATCCAGGTAT-GGCTGTTGAT_S73_R1_001.fastq.gz', '/data/h_vivian_weiss/DNAseq/1896/1896-AC-30-ATCCAGGTAT-GGCTGTTGAT_S73_R2_001.fastq.gz' ],
    'P1896_AC_31' => [ '/data/h_vivian_weiss/DNAseq/1896/1896-AC-31-CGATGCGGTT-CGCAAGAGCC_S74_R1_001.fastq.gz', '/data/h_vivian_weiss/DNAseq/1896/1896-AC-31-CGATGCGGTT-CGCAAGAGCC_S74_R2_001.fastq.gz' ],
    'P1896_AC_32' => [ '/data/h_vivian_weiss/DNAseq/1896/1896-AC-32-AATGCGAACA-GGTGAGTAAT_S75_R1_001.fastq.gz', '/data/h_vivian_weiss/DNAseq/1896/1896-AC-32-AATGCGAACA-GGTGAGTAAT_S75_R2_001.fastq.gz' ],
    'P1896_AC_33' => [ '/data/h_vivian_weiss/DNAseq/1896/1896-AC-33-TTCGGTGTGA-ACGCCTTGTT_S76_R1_001.fastq.gz', '/data/h_vivian_weiss/DNAseq/1896/1896-AC-33-TTCGGTGTGA-ACGCCTTGTT_S76_R2_001.fastq.gz' ],
    'P1896_AC_34' => [ '/data/h_vivian_weiss/DNAseq/1896/1896-AC-34-CATTAACTGA-CGAAGATTCT_S77_R1_001.fastq.gz', '/data/h_vivian_weiss/DNAseq/1896/1896-AC-34-CATTAACTGA-CGAAGATTCT_S77_R2_001.fastq.gz' ],
    'P1896_AC_35' => [ '/data/h_vivian_weiss/DNAseq/1896/1896-AC-35-CCGGTTCCTA-TCCATTGCCG_S78_R1_001.fastq.gz', '/data/h_vivian_weiss/DNAseq/1896/1896-AC-35-CCGGTTCCTA-TCCATTGCCG_S78_R2_001.fastq.gz' ],
    'P1896_AC_36' => [ '/data/h_vivian_weiss/DNAseq/1896/1896-AC-36-ATACATCACA-AACACTGTTA_S79_R1_001.fastq.gz', '/data/h_vivian_weiss/DNAseq/1896/1896-AC-36-ATACATCACA-AACACTGTTA_S79_R2_001.fastq.gz' ],
    'P1896_AC_37' => [ '/data/h_vivian_weiss/DNAseq/1896/1896-AC-37-TAGAGAATAC-GCCGCACTCT_S80_R1_001.fastq.gz', '/data/h_vivian_weiss/DNAseq/1896/1896-AC-37-TAGAGAATAC-GCCGCACTCT_S80_R2_001.fastq.gz' ],
    'P1896_AC_38' => [ '/data/h_vivian_weiss/DNAseq/1896/1896-AC-38-GTCTCGCCAC-CAATCGGCTG_S81_R1_001.fastq.gz', '/data/h_vivian_weiss/DNAseq/1896/1896-AC-38-GTCTCGCCAC-CAATCGGCTG_S81_R2_001.fastq.gz' ],
    'P1896_AC_39' => [ '/data/h_vivian_weiss/DNAseq/1896/1896-AC-39-CGCTGTCTCA-ACCACGACAT_S82_R1_001.fastq.gz', '/data/h_vivian_weiss/DNAseq/1896/1896-AC-39-CGCTGTCTCA-ACCACGACAT_S82_R2_001.fastq.gz' ],
    #file_def.py -i /data/h_vivian_weiss/DNAseq/1988/ -n "(........\\d+)" -f "*.gz" -a -p
    'P1988_AC_01' => [ '/data/h_vivian_weiss/DNAseq/1988/1988-AC-1-GCTAATAGGA-CTAGTGCTCT_S83_R1_001.fastq.gz', '/data/h_vivian_weiss/DNAseq/1988/1988-AC-1-GCTAATAGGA-CTAGTGCTCT_S83_R2_001.fastq.gz' ],
    'P1988_AC_02' => [ '/data/h_vivian_weiss/DNAseq/1988/1988-AC-2-CCTAACACAG-GTATTCCACC_S84_R1_001.fastq.gz', '/data/h_vivian_weiss/DNAseq/1988/1988-AC-2-CCTAACACAG-GTATTCCACC_S84_R2_001.fastq.gz' ],
    'P1988_AC_03' => [ '/data/h_vivian_weiss/DNAseq/1988/1988-AC-3-TGCCGGTCAG-TAGGAGCCTC_S85_R1_001.fastq.gz', '/data/h_vivian_weiss/DNAseq/1988/1988-AC-3-TGCCGGTCAG-TAGGAGCCTC_S85_R2_001.fastq.gz' ],
    'P1988_AC_04' => [ '/data/h_vivian_weiss/DNAseq/1988/1988-AC-4-TTAACCTTCG-CTTGCCATTA_S86_R1_001.fastq.gz', '/data/h_vivian_weiss/DNAseq/1988/1988-AC-4-TTAACCTTCG-CTTGCCATTA_S86_R2_001.fastq.gz' ],
    'P1988_AC_05' => [ '/data/h_vivian_weiss/DNAseq/1988/1988-AC-5-GCGTGCTGTG-GGTGTCACCG_S87_R1_001.fastq.gz', '/data/h_vivian_weiss/DNAseq/1988/1988-AC-5-GCGTGCTGTG-GGTGTCACCG_S87_R2_001.fastq.gz' ],
    'P1988_AC_06' => [ '/data/h_vivian_weiss/DNAseq/1988/1988-AC-6-CGAGAGGCGT-ATTATGTCTC_S88_R1_001.fastq.gz', '/data/h_vivian_weiss/DNAseq/1988/1988-AC-6-CGAGAGGCGT-ATTATGTCTC_S88_R2_001.fastq.gz' ],
    'P1988_AC_07' => [ '/data/h_vivian_weiss/DNAseq/1988/1988-AC-7-ACTCTATTGT-TGGCTAATCA_S89_R1_001.fastq.gz', '/data/h_vivian_weiss/DNAseq/1988/1988-AC-7-ACTCTATTGT-TGGCTAATCA_S89_R2_001.fastq.gz' ],
    'P1988_AC_08' => [ '/data/h_vivian_weiss/DNAseq/1988/1988-AC-8-TATCACTCTG-GTTGTAGTGC_S90_R1_001.fastq.gz', '/data/h_vivian_weiss/DNAseq/1988/1988-AC-8-TATCACTCTG-GTTGTAGTGC_S90_R2_001.fastq.gz' ],
    'P1988_AC_09' => [ '/data/h_vivian_weiss/DNAseq/1988/1988-AC-9-ATCGGCGAAG-TCGACATCTC_S91_R1_001.fastq.gz', '/data/h_vivian_weiss/DNAseq/1988/1988-AC-9-ATCGGCGAAG-TCGACATCTC_S91_R2_001.fastq.gz' ],
    'P1988_AC_10' => [ '/data/h_vivian_weiss/DNAseq/1988/1988-AC-10-GGTAGAATTA-TGTAAGGTGG_S92_R1_001.fastq.gz', '/data/h_vivian_weiss/DNAseq/1988/1988-AC-10-GGTAGAATTA-TGTAAGGTGG_S92_R2_001.fastq.gz' ],
    'P1988_AC_11' => [ '/data/h_vivian_weiss/DNAseq/1988/1988-AC-11-TCAGATTAAC-CCTCTACATG_S93_R1_001.fastq.gz', '/data/h_vivian_weiss/DNAseq/1988/1988-AC-11-TCAGATTAAC-CCTCTACATG_S93_R2_001.fastq.gz' ],
    'P1988_AC_12' => [ '/data/h_vivian_weiss/DNAseq/1988/1988-AC-12-ATCAGTACCA-AGCGAATTAG_S94_R1_001.fastq.gz', '/data/h_vivian_weiss/DNAseq/1988/1988-AC-12-ATCAGTACCA-AGCGAATTAG_S94_R2_001.fastq.gz' ],
    'P1988_AC_13' => [ '/data/h_vivian_weiss/DNAseq/1988/1988-AC-13-ACCTAAGACC-ACATATCCAG_S95_R1_001.fastq.gz', '/data/h_vivian_weiss/DNAseq/1988/1988-AC-13-ACCTAAGACC-ACATATCCAG_S95_R2_001.fastq.gz' ],
    'P1988_AC_14' => [ '/data/h_vivian_weiss/DNAseq/1988/1988-AC-14-TCGGCGGTTA-AGAACCGAGT_S96_R1_001.fastq.gz', '/data/h_vivian_weiss/DNAseq/1988/1988-AC-14-TCGGCGGTTA-AGAACCGAGT_S96_R2_001.fastq.gz' ],
    'P1988_AC_15' => [ '/data/h_vivian_weiss/DNAseq/1988/1988-AC-15-CTTATTGGCC-TTAAGACAAG_S97_R1_001.fastq.gz', '/data/h_vivian_weiss/DNAseq/1988/1988-AC-15-CTTATTGGCC-TTAAGACAAG_S97_R2_001.fastq.gz' ],
    'P1988_AC_16' => [ '/data/h_vivian_weiss/DNAseq/1988/1988-AC-16-GTTGGATGAA-TCATAGATTG_S98_R1_001.fastq.gz', '/data/h_vivian_weiss/DNAseq/1988/1988-AC-16-GTTGGATGAA-TCATAGATTG_S98_R2_001.fastq.gz' ],
    'P1988_AC_17' => [ '/data/h_vivian_weiss/DNAseq/1988/1988-AC-17-TGGACCGCCA-CGTCTCATAT_S99_R1_001.fastq.gz', '/data/h_vivian_weiss/DNAseq/1988/1988-AC-17-TGGACCGCCA-CGTCTCATAT_S99_R2_001.fastq.gz' ],
    'P1988_AC_18' => [ '/data/h_vivian_weiss/DNAseq/1988/1988-AC-18-AACGAGGCCG-CACGGAACAA_S100_R1_001.fastq.gz', '/data/h_vivian_weiss/DNAseq/1988/1988-AC-18-AACGAGGCCG-CACGGAACAA_S100_R2_001.fastq.gz' ],
    'P1988_AC_19' => [ '/data/h_vivian_weiss/DNAseq/1988/1988-AC-19-CTGAGCCGGT-TTCTCGTGCA_S101_R1_001.fastq.gz', '/data/h_vivian_weiss/DNAseq/1988/1988-AC-19-CTGAGCCGGT-TTCTCGTGCA_S101_R2_001.fastq.gz' ],
    'P1988_AC_20' => [ '/data/h_vivian_weiss/DNAseq/1988/1988-AC-20-GCTGACGTTG-GATAGGCCGA_S102_R1_001.fastq.gz', '/data/h_vivian_weiss/DNAseq/1988/1988-AC-20-GCTGACGTTG-GATAGGCCGA_S102_R2_001.fastq.gz' ],
    'P1988_AC_21' => [ '/data/h_vivian_weiss/DNAseq/1988/1988-AC-21-GTTCGGAGTT-GTGCGCTTGA_S103_R1_001.fastq.gz', '/data/h_vivian_weiss/DNAseq/1988/1988-AC-21-GTTCGGAGTT-GTGCGCTTGA_S103_R2_001.fastq.gz' ],
    'P1988_AC_22' => [ '/data/h_vivian_weiss/DNAseq/1988/1988-AC-22-CTAATAACCG-GAGGTTAGAC_S104_R1_001.fastq.gz', '/data/h_vivian_weiss/DNAseq/1988/1988-AC-22-CTAATAACCG-GAGGTTAGAC_S104_R2_001.fastq.gz' ],
    'P1988_AC_23' => [ '/data/h_vivian_weiss/DNAseq/1988/1988-AC-23-TCCGCCAATT-CCGGAGTGGA_S105_R1_001.fastq.gz', '/data/h_vivian_weiss/DNAseq/1988/1988-AC-23-TCCGCCAATT-CCGGAGTGGA_S105_R2_001.fastq.gz' ],
    'P1988_AC_24' => [ '/data/h_vivian_weiss/DNAseq/1988/1988-AC-24-ACCAAGCAGG-CTGCGAGCCA_S106_R1_001.fastq.gz', '/data/h_vivian_weiss/DNAseq/1988/1988-AC-24-ACCAAGCAGG-CTGCGAGCCA_S106_R2_001.fastq.gz' ],
    'P1988_AC_25' => [ '/data/h_vivian_weiss/DNAseq/1988/1988-AC-25-CAAGTTATTG-TACTCTGCGC_S107_R1_001.fastq.gz', '/data/h_vivian_weiss/DNAseq/1988/1988-AC-25-CAAGTTATTG-TACTCTGCGC_S107_R2_001.fastq.gz' ],
    'P1988_AC_26' => [ '/data/h_vivian_weiss/DNAseq/1988/1988-AC-26-TATCCGAGGC-CTCAACGCTT_S108_R1_001.fastq.gz', '/data/h_vivian_weiss/DNAseq/1988/1988-AC-26-TATCCGAGGC-CTCAACGCTT_S108_R2_001.fastq.gz' ],
    #file_def.py -i /data/h_vivian_weiss/3544/ -n "(........\\d+)" -f "*.gz" -a -p
    'P3544_SB_01' => [ '/data/h_vivian_weiss/3544/3544-SB-1-AGTGTTGC-ATGTAACG_S11_R1_001.fastq.gz', '/data/h_vivian_weiss/3544/3544-SB-1-AGTGTTGC-ATGTAACG_S11_R2_001.fastq.gz' ],
    'P3544_SB_02' => [ '/data/h_vivian_weiss/3544/3544-SB-2-GACACCAT-GCACGGTA_S12_R1_001.fastq.gz', '/data/h_vivian_weiss/3544/3544-SB-2-GACACCAT-GCACGGTA_S12_R2_001.fastq.gz' ],
    'P3544_SB_03' => [ '/data/h_vivian_weiss/3544/3544-SB-3-CTGAGGAA-TCTAGGCG_S16_R1_001.fastq.gz', '/data/h_vivian_weiss/3544/3544-SB-3-CTGAGGAA-TCTAGGCG_S16_R2_001.fastq.gz' ],
    'P3544_SB_04' => [ '/data/h_vivian_weiss/3544/3544-SB-4-TCAGAAGG-CTCGAATA_S25_R1_001.fastq.gz', '/data/h_vivian_weiss/3544/3544-SB-4-TCAGAAGG-CTCGAATA_S25_R2_001.fastq.gz' ],
    'P3544_SB_05' => [ '/data/h_vivian_weiss/3544/3544-SB-5-TACTCTGC-CAAGTTAT_S24_R1_001.fastq.gz', '/data/h_vivian_weiss/3544/3544-SB-5-TACTCTGC-CAAGTTAT_S24_R2_001.fastq.gz' ],
    'P3544_SB_06' => [ '/data/h_vivian_weiss/3544/3544-SB-6-CGTCTCAT-TGGACCGC_S29_R1_001.fastq.gz', '/data/h_vivian_weiss/3544/3544-SB-6-CGTCTCAT-TGGACCGC_S29_R2_001.fastq.gz' ],
    'P3544_SB_07' => [ '/data/h_vivian_weiss/3544/3544-SB-7-TCGACATC-ATCGGCGA_S21_R1_001.fastq.gz', '/data/h_vivian_weiss/3544/3544-SB-7-TCGACATC-ATCGGCGA_S21_R2_001.fastq.gz' ],
    'P3544_SB_08' => [ '/data/h_vivian_weiss/3544/3544-SB-8-CTAGTGCT-GCTAATAG_S15_R1_001.fastq.gz', '/data/h_vivian_weiss/3544/3544-SB-8-CTAGTGCT-GCTAATAG_S15_R2_001.fastq.gz' ],
    'P3544_SB_09' => [ '/data/h_vivian_weiss/3544/3544-SB-9-GGTGAGTA-AATGCGAA_S32_R1_001.fastq.gz', '/data/h_vivian_weiss/3544/3544-SB-9-GGTGAGTA-AATGCGAA_S32_R2_001.fastq.gz' ],
    'P3544_SB_10' => [ '/data/h_vivian_weiss/3544/3544-SB-10-AACAGACG-GGCATAGG_S22_R1_001.fastq.gz', '/data/h_vivian_weiss/3544/3544-SB-10-AACAGACG-GGCATAGG_S22_R2_001.fastq.gz' ],
    'P3544_SB_11' => [ '/data/h_vivian_weiss/3544/3544-SB-11-GGATGCTC-TGGCCTCT_S18_R1_001.fastq.gz', '/data/h_vivian_weiss/3544/3544-SB-11-GGATGCTC-TGGCCTCT_S18_R2_001.fastq.gz' ],
    'P3544_SB_12' => [ '/data/h_vivian_weiss/3544/3544-SB-12-AAGCATCT-CAATTCTC_S14_R1_001.fastq.gz', '/data/h_vivian_weiss/3544/3544-SB-12-AAGCATCT-CAATTCTC_S14_R2_001.fastq.gz' ],
    'P3544_SB_13' => [ '/data/h_vivian_weiss/3544/3544-SB-13-TTACCTGG-GATTCTGA_S26_R1_001.fastq.gz', '/data/h_vivian_weiss/3544/3544-SB-13-TTACCTGG-GATTCTGA_S26_R2_001.fastq.gz' ],
    'P3544_SB_14' => [ '/data/h_vivian_weiss/3544/3544-SB-14-CCGTTCAA-AGCCTCAG_S40_R1_001.fastq.gz', '/data/h_vivian_weiss/3544/3544-SB-14-CCGTTCAA-AGCCTCAG_S40_R2_001.fastq.gz' ],
    'P3544_SB_15' => [ '/data/h_vivian_weiss/3544/3544-SB-15-TGATGTAA-ATGTTGTT_S13_R1_001.fastq.gz', '/data/h_vivian_weiss/3544/3544-SB-15-TGATGTAA-ATGTTGTT_S13_R2_001.fastq.gz' ],
    'P3544_SB_16' => [ '/data/h_vivian_weiss/3544/3544-SB-16-CAGCACGG-GCACCACC_S31_R1_001.fastq.gz', '/data/h_vivian_weiss/3544/3544-SB-16-CAGCACGG-GCACCACC_S31_R2_001.fastq.gz' ],
    'P3544_SB_17' => [ '/data/h_vivian_weiss/3544/3544-SB-17-CTGCGAGC-ACCAAGCA_S35_R1_001.fastq.gz', '/data/h_vivian_weiss/3544/3544-SB-17-CTGCGAGC-ACCAAGCA_S35_R2_001.fastq.gz' ],
    'P3544_SB_18' => [ '/data/h_vivian_weiss/3544/3544-SB-18-TCATAGAT-GTTGGATG_S45_R1_001.fastq.gz', '/data/h_vivian_weiss/3544/3544-SB-18-TCATAGAT-GTTGGATG_S45_R2_001.fastq.gz' ],
    'P3544_SB_19' => [ '/data/h_vivian_weiss/3544/3544-SB-19-GTTGTAGT-TATCACTC_S41_R1_001.fastq.gz', '/data/h_vivian_weiss/3544/3544-SB-19-GTTGTAGT-TATCACTC_S41_R2_001.fastq.gz' ],
    'P3544_SB_20' => [ '/data/h_vivian_weiss/3544/3544-SB-20-ACCACGAC-CGCTGTCT_S43_R1_001.fastq.gz', '/data/h_vivian_weiss/3544/3544-SB-20-ACCACGAC-CGCTGTCT_S43_R2_001.fastq.gz' ],
    'P3544_SB_21' => [ '/data/h_vivian_weiss/3544/3544-SB-21-CGCAAGAG-CGATGCGG_S42_R1_001.fastq.gz', '/data/h_vivian_weiss/3544/3544-SB-21-CGCAAGAG-CGATGCGG_S42_R2_001.fastq.gz' ],
    'P3544_SB_22' => [ '/data/h_vivian_weiss/3544/3544-SB-22-TATGGAGA-TAGCATAA_S30_R1_001.fastq.gz', '/data/h_vivian_weiss/3544/3544-SB-22-TATGGAGA-TAGCATAA_S30_R2_001.fastq.gz' ],
    'P3544_SB_23' => [ '/data/h_vivian_weiss/3544/3544-SB-23-AGGTCCTT-GTCCACCG_S44_R1_001.fastq.gz', '/data/h_vivian_weiss/3544/3544-SB-23-AGGTCCTT-GTCCACCG_S44_R2_001.fastq.gz' ],
    'P3544_SB_24' => [ '/data/h_vivian_weiss/3544/3544-SB-24-GAACTTCC-ACTTGTTA_S37_R1_001.fastq.gz', '/data/h_vivian_weiss/3544/3544-SB-24-GAACTTCC-ACTTGTTA_S37_R2_001.fastq.gz' ],
    'P3544_SB_25' => [ '/data/h_vivian_weiss/3544/3544-SB-25-TCTATCCT-GAGAGGTT_S20_R1_001.fastq.gz', '/data/h_vivian_weiss/3544/3544-SB-25-TCTATCCT-GAGAGGTT_S20_R2_001.fastq.gz' ],
    'P3544_SB_26' => [ '/data/h_vivian_weiss/3544/3544-SB-26-CTCGCTTC-AGAGAACC_S34_R1_001.fastq.gz', '/data/h_vivian_weiss/3544/3544-SB-26-CTCGCTTC-AGAGAACC_S34_R2_001.fastq.gz' ],
    'P3544_SB_27' => [ '/data/h_vivian_weiss/3544/3544-SB-27-GGTCCGCT-CTGGCAAG_S23_R1_001.fastq.gz', '/data/h_vivian_weiss/3544/3544-SB-27-GGTCCGCT-CTGGCAAG_S23_R2_001.fastq.gz' ],
    'P3544_SB_28' => [ '/data/h_vivian_weiss/3544/3544-SB-28-AACTTATC-TCAATGGA_S19_R1_001.fastq.gz', '/data/h_vivian_weiss/3544/3544-SB-28-AACTTATC-TCAATGGA_S19_R2_001.fastq.gz' ],
    'P3544_SB_29' => [ '/data/h_vivian_weiss/3544/3544-SB-29-CCGGAGTG-TCCGCCAA_S27_R1_001.fastq.gz', '/data/h_vivian_weiss/3544/3544-SB-29-CCGGAGTG-TCCGCCAA_S27_R2_001.fastq.gz' ],
    'P3544_SB_30' => [ '/data/h_vivian_weiss/3544/3544-SB-30-TTAAGACA-CTTATTGG_S28_R1_001.fastq.gz', '/data/h_vivian_weiss/3544/3544-SB-30-TTAAGACA-CTTATTGG_S28_R2_001.fastq.gz' ],
    'P3544_SB_31' => [ '/data/h_vivian_weiss/3544/3544-SB-31-TGGCTAAT-ACTCTATT_S39_R1_001.fastq.gz', '/data/h_vivian_weiss/3544/3544-SB-31-TGGCTAAT-ACTCTATT_S39_R2_001.fastq.gz' ],
    'P3544_SB_32' => [ '/data/h_vivian_weiss/3544/3544-SB-32-CAATCGGC-GTCTCGCC_S17_R1_001.fastq.gz', '/data/h_vivian_weiss/3544/3544-SB-32-CAATCGGC-GTCTCGCC_S17_R2_001.fastq.gz' ],
    'P3544_SB_33' => [ '/data/h_vivian_weiss/3544/3544-SB-33-GGCTGTTG-ATCCAGGT_S36_R1_001.fastq.gz', '/data/h_vivian_weiss/3544/3544-SB-33-GGCTGTTG-ATCCAGGT_S36_R2_001.fastq.gz' ],
    'P3544_SB_34' => [ '/data/h_vivian_weiss/3544/3544-SB-34-AATCACCA-GCTTGAAC_S33_R1_001.fastq.gz', '/data/h_vivian_weiss/3544/3544-SB-34-AATCACCA-GCTTGAAC_S33_R2_001.fastq.gz' ],
    'P3544_SB_35' => [ '/data/h_vivian_weiss/3544/3544-SB-35-GTCCTGGA-AAGTTGAC_S38_R1_001.fastq.gz', '/data/h_vivian_weiss/3544/3544-SB-35-GTCCTGGA-AAGTTGAC_S38_R2_001.fastq.gz' ],
    #file_def.py -i /data/h_vivian_weiss/4271/ -n "(........\\d+)" -f "*.gz" -a -p
    'P4271_CP_01' => [ '/data/h_vivian_weiss/4271/4271-CP-1-AGTGTTGC-ATGTAACG_S99_R1_001.fastq.gz', '/data/h_vivian_weiss/4271/4271-CP-1-AGTGTTGC-ATGTAACG_S99_R2_001.fastq.gz' ],
    'P4271_CP_02' => [ '/data/h_vivian_weiss/4271/4271-CP-2-GACACCAT-GCACGGTA_S100_R1_001.fastq.gz', '/data/h_vivian_weiss/4271/4271-CP-2-GACACCAT-GCACGGTA_S100_R2_001.fastq.gz' ],
    'P4271_CP_03' => [ '/data/h_vivian_weiss/4271/4271-CP-3-CTGAGGAA-TCTAGGCG_S101_R1_001.fastq.gz', '/data/h_vivian_weiss/4271/4271-CP-3-CTGAGGAA-TCTAGGCG_S101_R2_001.fastq.gz' ],
    'P4271_CP_04' => [ '/data/h_vivian_weiss/4271/4271-CP-4-TCAGAAGG-CTCGAATA_S109_R1_001.fastq.gz', '/data/h_vivian_weiss/4271/4271-CP-4-TCAGAAGG-CTCGAATA_S109_R2_001.fastq.gz' ],
    'P4271_CP_05' => [ '/data/h_vivian_weiss/4271/4271-CP-5-TACTCTGC-CAAGTTAT_S110_R1_001.fastq.gz', '/data/h_vivian_weiss/4271/4271-CP-5-TACTCTGC-CAAGTTAT_S110_R2_001.fastq.gz' ],
    'P4271_CP_06' => [ '/data/h_vivian_weiss/4271/4271-CP-6-CGTCTCAT-TGGACCGC_S111_R1_001.fastq.gz', '/data/h_vivian_weiss/4271/4271-CP-6-CGTCTCAT-TGGACCGC_S111_R2_001.fastq.gz' ],
    'P4271_CP_07' => [ '/data/h_vivian_weiss/4271/4271-CP-7-TCGACATC-ATCGGCGA_S89_R1_001.fastq.gz', '/data/h_vivian_weiss/4271/4271-CP-7-TCGACATC-ATCGGCGA_S89_R2_001.fastq.gz' ],
    'P4271_CP_08' => [ '/data/h_vivian_weiss/4271/4271-CP-8-CTAGTGCT-GCTAATAG_S112_R1_001.fastq.gz', '/data/h_vivian_weiss/4271/4271-CP-8-CTAGTGCT-GCTAATAG_S112_R2_001.fastq.gz' ],
    'P4271_CP_09' => [ '/data/h_vivian_weiss/4271/4271-CP-9-GGTGAGTA-AATGCGAA_S90_R1_001.fastq.gz', '/data/h_vivian_weiss/4271/4271-CP-9-GGTGAGTA-AATGCGAA_S90_R2_001.fastq.gz' ],
    'P4271_CP_10' => [ '/data/h_vivian_weiss/4271/4271-CP-10-AACAGACG-GGCATAGG_S102_R1_001.fastq.gz', '/data/h_vivian_weiss/4271/4271-CP-10-AACAGACG-GGCATAGG_S102_R2_001.fastq.gz' ],
    'P4271_CP_11' => [ '/data/h_vivian_weiss/4271/4271-CP-11-GGATGCTC-TGGCCTCT_S113_R1_001.fastq.gz', '/data/h_vivian_weiss/4271/4271-CP-11-GGATGCTC-TGGCCTCT_S113_R2_001.fastq.gz' ],
    'P4271_CP_12' => [ '/data/h_vivian_weiss/4271/4271-CP-12-AAGCATCT-CAATTCTC_S103_R1_001.fastq.gz', '/data/h_vivian_weiss/4271/4271-CP-12-AAGCATCT-CAATTCTC_S103_R2_001.fastq.gz' ],
    'P4271_CP_13' => [ '/data/h_vivian_weiss/4271/4271-CP-13-TTACCTGG-GATTCTGA_S114_R1_001.fastq.gz', '/data/h_vivian_weiss/4271/4271-CP-13-TTACCTGG-GATTCTGA_S114_R2_001.fastq.gz' ],
    'P4271_CP_14' => [ '/data/h_vivian_weiss/4271/4271-CP-14-CCGTTCAA-AGCCTCAG_S104_R1_001.fastq.gz', '/data/h_vivian_weiss/4271/4271-CP-14-CCGTTCAA-AGCCTCAG_S104_R2_001.fastq.gz' ],
    'P4271_CP_15' => [ '/data/h_vivian_weiss/4271/4271-CP-15-TGATGTAA-ATGTTGTT_S91_R1_001.fastq.gz', '/data/h_vivian_weiss/4271/4271-CP-15-TGATGTAA-ATGTTGTT_S91_R2_001.fastq.gz' ],
    'P4271_CP_16' => [ '/data/h_vivian_weiss/4271/4271-CP-16-CAGCACGG-GCACCACC_S81_R1_001.fastq.gz', '/data/h_vivian_weiss/4271/4271-CP-16-CAGCACGG-GCACCACC_S81_R2_001.fastq.gz' ],
    'P4271_CP_17' => [ '/data/h_vivian_weiss/4271/4271-CP-17-CTGCGAGC-ACCAAGCA_S82_R1_001.fastq.gz', '/data/h_vivian_weiss/4271/4271-CP-17-CTGCGAGC-ACCAAGCA_S82_R2_001.fastq.gz' ],
    'P4271_CP_18' => [ '/data/h_vivian_weiss/4271/4271-CP-18-TCATAGAT-GTTGGATG_S92_R1_001.fastq.gz', '/data/h_vivian_weiss/4271/4271-CP-18-TCATAGAT-GTTGGATG_S92_R2_001.fastq.gz' ],
    'P4271_CP_19' => [ '/data/h_vivian_weiss/4271/4271-CP-19-GTTGTAGT-TATCACTC_S93_R1_001.fastq.gz', '/data/h_vivian_weiss/4271/4271-CP-19-GTTGTAGT-TATCACTC_S93_R2_001.fastq.gz' ],
    'P4271_CP_20' => [ '/data/h_vivian_weiss/4271/4271-CP-20-ACCACGAC-CGCTGTCT_S105_R1_001.fastq.gz', '/data/h_vivian_weiss/4271/4271-CP-20-ACCACGAC-CGCTGTCT_S105_R2_001.fastq.gz' ],
    'P4271_CP_21' => [ '/data/h_vivian_weiss/4271/4271-CP-21-CGCAAGAG-CGATGCGG_S94_R1_001.fastq.gz', '/data/h_vivian_weiss/4271/4271-CP-21-CGCAAGAG-CGATGCGG_S94_R2_001.fastq.gz' ],
    'P4271_CP_22' => [ '/data/h_vivian_weiss/4271/4271-CP-22-TATGGAGA-TAGCATAA_S106_R1_001.fastq.gz', '/data/h_vivian_weiss/4271/4271-CP-22-TATGGAGA-TAGCATAA_S106_R2_001.fastq.gz' ],
    'P4271_CP_23' => [ '/data/h_vivian_weiss/4271/4271-CP-23-AGGTCCTT-GTCCACCG_S115_R1_001.fastq.gz', '/data/h_vivian_weiss/4271/4271-CP-23-AGGTCCTT-GTCCACCG_S115_R2_001.fastq.gz' ],
    'P4271_CP_24' => [ '/data/h_vivian_weiss/4271/4271-CP-24-GAACTTCC-ACTTGTTA_S95_R1_001.fastq.gz', '/data/h_vivian_weiss/4271/4271-CP-24-GAACTTCC-ACTTGTTA_S95_R2_001.fastq.gz' ],
    'P4271_CP_25' => [ '/data/h_vivian_weiss/4271/4271-CP-25-TCTATCCT-GAGAGGTT_S116_R1_001.fastq.gz', '/data/h_vivian_weiss/4271/4271-CP-25-TCTATCCT-GAGAGGTT_S116_R2_001.fastq.gz' ],
    'P4271_CP_26' => [ '/data/h_vivian_weiss/4271/4271-CP-26-CTCGCTTC-AGAGAACC_S117_R1_001.fastq.gz', '/data/h_vivian_weiss/4271/4271-CP-26-CTCGCTTC-AGAGAACC_S117_R2_001.fastq.gz' ],
    'P4271_CP_27' => [ '/data/h_vivian_weiss/4271/4271-CP-27-GGTCCGCT-CTGGCAAG_S107_R1_001.fastq.gz', '/data/h_vivian_weiss/4271/4271-CP-27-GGTCCGCT-CTGGCAAG_S107_R2_001.fastq.gz' ],
    'P4271_CP_28' => [ '/data/h_vivian_weiss/4271/4271-CP-28-AACTTATC-TCAATGGA_S83_R1_001.fastq.gz', '/data/h_vivian_weiss/4271/4271-CP-28-AACTTATC-TCAATGGA_S83_R2_001.fastq.gz' ],
    'P4271_CP_29' => [ '/data/h_vivian_weiss/4271/4271-CP-29-CCGGAGTG-TCCGCCAA_S118_R1_001.fastq.gz', '/data/h_vivian_weiss/4271/4271-CP-29-CCGGAGTG-TCCGCCAA_S118_R2_001.fastq.gz' ],
    'P4271_CP_30' => [ '/data/h_vivian_weiss/4271/4271-CP-30-TTAAGACA-CTTATTGG_S84_R1_001.fastq.gz', '/data/h_vivian_weiss/4271/4271-CP-30-TTAAGACA-CTTATTGG_S84_R2_001.fastq.gz' ],
    'P4271_CP_32' => [ '/data/h_vivian_weiss/4271/4271-CP-32-CAATCGGC-GTCTCGCC_S85_R1_001.fastq.gz', '/data/h_vivian_weiss/4271/4271-CP-32-CAATCGGC-GTCTCGCC_S85_R2_001.fastq.gz' ],
    'P4271_CP_33' => [ '/data/h_vivian_weiss/4271/4271-CP-33-GGCTGTTG-ATCCAGGT_S86_R1_001.fastq.gz', '/data/h_vivian_weiss/4271/4271-CP-33-GGCTGTTG-ATCCAGGT_S86_R2_001.fastq.gz' ],
    'P4271_CP_34' => [ '/data/h_vivian_weiss/4271/4271-CP-34-AATCACCA-GCTTGAAC_S96_R1_001.fastq.gz', '/data/h_vivian_weiss/4271/4271-CP-34-AATCACCA-GCTTGAAC_S96_R2_001.fastq.gz' ],
    'P4271_CP_35' => [ '/data/h_vivian_weiss/4271/4271-CP-35-GTCCTGGA-AAGTTGAC_S108_R1_001.fastq.gz', '/data/h_vivian_weiss/4271/4271-CP-35-GTCCTGGA-AAGTTGAC_S108_R2_001.fastq.gz' ],
    'P4271_CP_36' => [ '/data/h_vivian_weiss/4271/4271-CP-36-ACTTCAAG-GGACCAGT_S97_R1_001.fastq.gz', '/data/h_vivian_weiss/4271/4271-CP-36-ACTTCAAG-GGACCAGT_S97_R2_001.fastq.gz' ],
    'P4271_CP_37' => [ '/data/h_vivian_weiss/4271/4271-CP-37-TTCTACAT-TTGTATCA_S87_R1_001.fastq.gz', '/data/h_vivian_weiss/4271/4271-CP-37-TTCTACAT-TTGTATCA_S87_R2_001.fastq.gz' ],
    'P4271_CP_38' => [ '/data/h_vivian_weiss/4271/4271-CP-38-CCTCGTGC-CCACGCTG_S88_R1_001.fastq.gz', '/data/h_vivian_weiss/4271/4271-CP-38-CCTCGTGC-CCACGCTG_S88_R2_001.fastq.gz' ],
    'P4271_CP_39' => [ '/data/h_vivian_weiss/4271/4271-CP-39-ACCGGCCG-GAATACCT_S98_R1_001.fastq.gz', '/data/h_vivian_weiss/4271/4271-CP-39-ACCGGCCG-GAATACCT_S98_R2_001.fastq.gz' ],
    #file_def.py -i /data/h_vivian_weiss/4601/ -n "(........\\d+)" -f "*.gz" -a -p
    'P4601_CP_01' => [ '/data/h_vivian_weiss/4601/4601-CP-1_S1_L005_R1_001.fastq.gz', '/data/h_vivian_weiss/4601/4601-CP-1_S1_L005_R2_001.fastq.gz' ],
    'P4601_CP_02' => [ '/data/h_vivian_weiss/4601/4601-CP-2_S6_L005_R1_001.fastq.gz', '/data/h_vivian_weiss/4601/4601-CP-2_S6_L005_R2_001.fastq.gz' ],
    'P4601_CP_03' => [ '/data/h_vivian_weiss/4601/4601-CP-3_S30_L005_R1_001.fastq.gz', '/data/h_vivian_weiss/4601/4601-CP-3_S30_L005_R2_001.fastq.gz' ],
    'P4601_CP_04' => [ '/data/h_vivian_weiss/4601/4601-CP-4_S7_L005_R1_001.fastq.gz', '/data/h_vivian_weiss/4601/4601-CP-4_S7_L005_R2_001.fastq.gz' ],
    'P4601_CP_05' => [ '/data/h_vivian_weiss/4601/4601-CP-5_S8_L005_R1_001.fastq.gz', '/data/h_vivian_weiss/4601/4601-CP-5_S8_L005_R2_001.fastq.gz' ],
    'P4601_CP_06' => [ '/data/h_vivian_weiss/4601/4601-CP-6_S9_L005_R1_001.fastq.gz', '/data/h_vivian_weiss/4601/4601-CP-6_S9_L005_R2_001.fastq.gz' ],
    'P4601_CP_07' => [ '/data/h_vivian_weiss/4601/4601-CP-7_S42_L005_R1_001.fastq.gz', '/data/h_vivian_weiss/4601/4601-CP-7_S42_L005_R2_001.fastq.gz' ],
    'P4601_CP_08' => [ '/data/h_vivian_weiss/4601/4601-CP-8_S18_L005_R1_001.fastq.gz', '/data/h_vivian_weiss/4601/4601-CP-8_S18_L005_R2_001.fastq.gz' ],
    'P4601_CP_09' => [ '/data/h_vivian_weiss/4601/4601-CP-9_S10_L005_R1_001.fastq.gz', '/data/h_vivian_weiss/4601/4601-CP-9_S10_L005_R2_001.fastq.gz' ],
    'P4601_CP_10' => [ '/data/h_vivian_weiss/4601/4601-CP-10_S11_L005_R1_001.fastq.gz', '/data/h_vivian_weiss/4601/4601-CP-10_S11_L005_R2_001.fastq.gz' ],
    'P4601_CP_11' => [ '/data/h_vivian_weiss/4601/4601-CP-11_S12_L005_R1_001.fastq.gz', '/data/h_vivian_weiss/4601/4601-CP-11_S12_L005_R2_001.fastq.gz' ],
    'P4601_CP_12' => [ '/data/h_vivian_weiss/4601/4601-CP-12_S2_L005_R1_001.fastq.gz', '/data/h_vivian_weiss/4601/4601-CP-12_S2_L005_R2_001.fastq.gz' ],
    'P4601_CP_13' => [ '/data/h_vivian_weiss/4601/4601-CP-13_S19_L005_R1_001.fastq.gz', '/data/h_vivian_weiss/4601/4601-CP-13_S19_L005_R2_001.fastq.gz' ],
    'P4601_CP_14' => [ '/data/h_vivian_weiss/4601/4601-CP-14_S13_L005_R1_001.fastq.gz', '/data/h_vivian_weiss/4601/4601-CP-14_S13_L005_R2_001.fastq.gz' ],
    'P4601_CP_15' => [ '/data/h_vivian_weiss/4601/4601-CP-15_S3_L005_R1_001.fastq.gz', '/data/h_vivian_weiss/4601/4601-CP-15_S3_L005_R2_001.fastq.gz' ],
    'P4601_CP_16' => [ '/data/h_vivian_weiss/4601/4601-CP-16_S14_L005_R1_001.fastq.gz', '/data/h_vivian_weiss/4601/4601-CP-16_S14_L005_R2_001.fastq.gz' ],
    'P4601_CP_17' => [ '/data/h_vivian_weiss/4601/4601-CP-17_S20_L005_R1_001.fastq.gz', '/data/h_vivian_weiss/4601/4601-CP-17_S20_L005_R2_001.fastq.gz' ],
    #'P4601_CP_18' => [ '/data/h_vivian_weiss/4601/4601-CP-18_S31_L005_R1_001.fastq.gz', '/data/h_vivian_weiss/4601/4601-CP-18_S31_L005_R2_001.fastq.gz' ],
    'P4601_CP_19' => [ '/data/h_vivian_weiss/4601/4601-CP-19_S32_L005_R1_001.fastq.gz', '/data/h_vivian_weiss/4601/4601-CP-19_S32_L005_R2_001.fastq.gz' ],
    'P4601_CP_20' => [ '/data/h_vivian_weiss/4601/4601-CP-20_S43_L005_R1_001.fastq.gz', '/data/h_vivian_weiss/4601/4601-CP-20_S43_L005_R2_001.fastq.gz' ],
    'P4601_CP_21' => [ '/data/h_vivian_weiss/4601/4601-CP-21_S17_L005_R1_001.fastq.gz', '/data/h_vivian_weiss/4601/4601-CP-21_S17_L005_R2_001.fastq.gz' ],
    'P4601_CP_22' => [ '/data/h_vivian_weiss/4601/4601-CP-22_S33_L005_R1_001.fastq.gz', '/data/h_vivian_weiss/4601/4601-CP-22_S33_L005_R2_001.fastq.gz' ],
    'P4601_CP_23' => [ '/data/h_vivian_weiss/4601/4601-CP-23_S34_L005_R1_001.fastq.gz', '/data/h_vivian_weiss/4601/4601-CP-23_S34_L005_R2_001.fastq.gz' ],
    'P4601_CP_24' => [ '/data/h_vivian_weiss/4601/4601-CP-24_S15_L005_R1_001.fastq.gz', '/data/h_vivian_weiss/4601/4601-CP-24_S15_L005_R2_001.fastq.gz' ],
    'P4601_CP_25' => [ '/data/h_vivian_weiss/4601/4601-CP-25_S21_L005_R1_001.fastq.gz', '/data/h_vivian_weiss/4601/4601-CP-25_S21_L005_R2_001.fastq.gz' ],
    'P4601_CP_26' => [ '/data/h_vivian_weiss/4601/4601-CP-26_S35_L005_R1_001.fastq.gz', '/data/h_vivian_weiss/4601/4601-CP-26_S35_L005_R2_001.fastq.gz' ],
    'P4601_CP_27' => [ '/data/h_vivian_weiss/4601/4601-CP-27_S22_L005_R1_001.fastq.gz', '/data/h_vivian_weiss/4601/4601-CP-27_S22_L005_R2_001.fastq.gz' ],
    'P4601_CP_28' => [ '/data/h_vivian_weiss/4601/4601-CP-28_S54_L005_R1_001.fastq.gz', '/data/h_vivian_weiss/4601/4601-CP-28_S54_L005_R2_001.fastq.gz' ],
    'P4601_CP_29' => [ '/data/h_vivian_weiss/4601/4601-CP-29_S36_L005_R1_001.fastq.gz', '/data/h_vivian_weiss/4601/4601-CP-29_S36_L005_R2_001.fastq.gz' ],
    'P4601_CP_30' => [ '/data/h_vivian_weiss/4601/4601-CP-30_S44_L005_R1_001.fastq.gz', '/data/h_vivian_weiss/4601/4601-CP-30_S44_L005_R2_001.fastq.gz' ],
    'P4601_CP_31' => [ '/data/h_vivian_weiss/4601/4601-CP-31_S23_L005_R1_001.fastq.gz', '/data/h_vivian_weiss/4601/4601-CP-31_S23_L005_R2_001.fastq.gz' ],
    'P4601_CP_32' => [ '/data/h_vivian_weiss/4601/4601-CP-32_S45_L005_R1_001.fastq.gz', '/data/h_vivian_weiss/4601/4601-CP-32_S45_L005_R2_001.fastq.gz' ],
    'P4601_CP_33' => [ '/data/h_vivian_weiss/4601/4601-CP-33_S46_L005_R1_001.fastq.gz', '/data/h_vivian_weiss/4601/4601-CP-33_S46_L005_R2_001.fastq.gz' ],
    'P4601_CP_34' => [ '/data/h_vivian_weiss/4601/4601-CP-34_S24_L005_R1_001.fastq.gz', '/data/h_vivian_weiss/4601/4601-CP-34_S24_L005_R2_001.fastq.gz' ],
    'P4601_CP_35' => [ '/data/h_vivian_weiss/4601/4601-CP-35_S37_L005_R1_001.fastq.gz', '/data/h_vivian_weiss/4601/4601-CP-35_S37_L005_R2_001.fastq.gz' ],
    'P4601_CP_36' => [ '/data/h_vivian_weiss/4601/4601-CP-36_S25_L005_R1_001.fastq.gz', '/data/h_vivian_weiss/4601/4601-CP-36_S25_L005_R2_001.fastq.gz' ],
    'P4601_CP_37' => [ '/data/h_vivian_weiss/4601/4601-CP-37_S47_L005_R1_001.fastq.gz', '/data/h_vivian_weiss/4601/4601-CP-37_S47_L005_R2_001.fastq.gz' ],
    'P4601_CP_38' => [ '/data/h_vivian_weiss/4601/4601-CP-38_S26_L005_R1_001.fastq.gz', '/data/h_vivian_weiss/4601/4601-CP-38_S26_L005_R2_001.fastq.gz' ],
    'P4601_CP_39' => [ '/data/h_vivian_weiss/4601/4601-CP-39_S48_L005_R1_001.fastq.gz', '/data/h_vivian_weiss/4601/4601-CP-39_S48_L005_R2_001.fastq.gz' ],
    'P4601_CP_40' => [ '/data/h_vivian_weiss/4601/4601-CP-40_S27_L005_R1_001.fastq.gz', '/data/h_vivian_weiss/4601/4601-CP-40_S27_L005_R2_001.fastq.gz' ],
    'P4601_CP_41' => [ '/data/h_vivian_weiss/4601/4601-CP-41_S49_L005_R1_001.fastq.gz', '/data/h_vivian_weiss/4601/4601-CP-41_S49_L005_R2_001.fastq.gz' ],
    'P4601_CP_42' => [ '/data/h_vivian_weiss/4601/4601-CP-42_S55_L005_R1_001.fastq.gz', '/data/h_vivian_weiss/4601/4601-CP-42_S55_L005_R2_001.fastq.gz' ],
    'P4601_CP_43' => [ '/data/h_vivian_weiss/4601/4601-CP-43_S56_L005_R1_001.fastq.gz', '/data/h_vivian_weiss/4601/4601-CP-43_S56_L005_R2_001.fastq.gz' ],
    'P4601_CP_44' => [ '/data/h_vivian_weiss/4601/4601-CP-44_S50_L005_R1_001.fastq.gz', '/data/h_vivian_weiss/4601/4601-CP-44_S50_L005_R2_001.fastq.gz' ],
    'P4601_CP_45' => [ '/data/h_vivian_weiss/4601/4601-CP-45_S38_L005_R1_001.fastq.gz', '/data/h_vivian_weiss/4601/4601-CP-45_S38_L005_R2_001.fastq.gz' ],
    'P4601_CP_46' => [ '/data/h_vivian_weiss/4601/4601-CP-46_S39_L005_R1_001.fastq.gz', '/data/h_vivian_weiss/4601/4601-CP-46_S39_L005_R2_001.fastq.gz' ],
    'P4601_CP_47' => [ '/data/h_vivian_weiss/4601/4601-CP-47_S57_L005_R1_001.fastq.gz', '/data/h_vivian_weiss/4601/4601-CP-47_S57_L005_R2_001.fastq.gz' ],
    'P4601_CP_48' => [ '/data/h_vivian_weiss/4601/4601-CP-48_S58_L005_R1_001.fastq.gz', '/data/h_vivian_weiss/4601/4601-CP-48_S58_L005_R2_001.fastq.gz' ],
    'P4601_CP_49' => [ '/data/h_vivian_weiss/4601/4601-CP-49_S51_L005_R1_001.fastq.gz', '/data/h_vivian_weiss/4601/4601-CP-49_S51_L005_R2_001.fastq.gz' ],
    'P4601_CP_50' => [ '/data/h_vivian_weiss/4601/4601-CP-50_S64_L005_R1_001.fastq.gz', '/data/h_vivian_weiss/4601/4601-CP-50_S64_L005_R2_001.fastq.gz' ],
    'P4601_CP_51' => [ '/data/h_vivian_weiss/4601/4601-CP-51_S4_L005_R1_001.fastq.gz', '/data/h_vivian_weiss/4601/4601-CP-51_S4_L005_R2_001.fastq.gz' ],
    #'P4601_CP_52' => [ '/data/h_vivian_weiss/4601/4601-CP-52_S29_L005_R1_001.fastq.gz', '/data/h_vivian_weiss/4601/4601-CP-52_S29_L005_R2_001.fastq.gz' ],
    'P4601_CP_53' => [ '/data/h_vivian_weiss/4601/4601-CP-53_S59_L005_R1_001.fastq.gz', '/data/h_vivian_weiss/4601/4601-CP-53_S59_L005_R2_001.fastq.gz' ],
    'P4601_CP_54' => [ '/data/h_vivian_weiss/4601/4601-CP-54_S28_L005_R1_001.fastq.gz', '/data/h_vivian_weiss/4601/4601-CP-54_S28_L005_R2_001.fastq.gz' ],
    'P4601_CP_55' => [ '/data/h_vivian_weiss/4601/4601-CP-55_S52_L005_R1_001.fastq.gz', '/data/h_vivian_weiss/4601/4601-CP-55_S52_L005_R2_001.fastq.gz' ],
    'P4601_CP_56' => [ '/data/h_vivian_weiss/4601/4601-CP-56_S40_L005_R1_001.fastq.gz', '/data/h_vivian_weiss/4601/4601-CP-56_S40_L005_R2_001.fastq.gz' ],
    'P4601_CP_57' => [ '/data/h_vivian_weiss/4601/4601-CP-57_S65_L005_R1_001.fastq.gz', '/data/h_vivian_weiss/4601/4601-CP-57_S65_L005_R2_001.fastq.gz' ],
    'P4601_CP_58' => [ '/data/h_vivian_weiss/4601/4601-CP-58_S53_L005_R1_001.fastq.gz', '/data/h_vivian_weiss/4601/4601-CP-58_S53_L005_R2_001.fastq.gz' ],
    'P4601_CP_59' => [ '/data/h_vivian_weiss/4601/4601-CP-59_S16_L005_R1_001.fastq.gz', '/data/h_vivian_weiss/4601/4601-CP-59_S16_L005_R2_001.fastq.gz' ],
    'P4601_CP_60' => [ '/data/h_vivian_weiss/4601/4601-CP-60_S60_L005_R1_001.fastq.gz', '/data/h_vivian_weiss/4601/4601-CP-60_S60_L005_R2_001.fastq.gz' ],
    'P4601_CP_61' => [ '/data/h_vivian_weiss/4601/4601-CP-61_S41_L005_R1_001.fastq.gz', '/data/h_vivian_weiss/4601/4601-CP-61_S41_L005_R2_001.fastq.gz' ],
    'P4601_CP_62' => [ '/data/h_vivian_weiss/4601/4601-CP-62_S61_L005_R1_001.fastq.gz', '/data/h_vivian_weiss/4601/4601-CP-62_S61_L005_R2_001.fastq.gz' ],
    'P4601_CP_63' => [ '/data/h_vivian_weiss/4601/4601-CP-63_S66_L005_R1_001.fastq.gz', '/data/h_vivian_weiss/4601/4601-CP-63_S66_L005_R2_001.fastq.gz' ],
    'P4601_CP_64' => [ '/data/h_vivian_weiss/4601/4601-CP-64_S67_L005_R1_001.fastq.gz', '/data/h_vivian_weiss/4601/4601-CP-64_S67_L005_R2_001.fastq.gz' ],
    'P4601_CP_65' => [ '/data/h_vivian_weiss/4601/4601-CP-65_S62_L005_R1_001.fastq.gz', '/data/h_vivian_weiss/4601/4601-CP-65_S62_L005_R2_001.fastq.gz' ],
    'P4601_CP_66' => [ '/data/h_vivian_weiss/4601/4601-CP-66_S68_L005_R1_001.fastq.gz', '/data/h_vivian_weiss/4601/4601-CP-66_S68_L005_R2_001.fastq.gz' ],
    'P4601_CP_67' => [ '/data/h_vivian_weiss/4601/4601-CP-67_S69_L005_R1_001.fastq.gz', '/data/h_vivian_weiss/4601/4601-CP-67_S69_L005_R2_001.fastq.gz' ],
    'P4601_CP_68' => [ '/data/h_vivian_weiss/4601/4601-CP-68_S70_L005_R1_001.fastq.gz', '/data/h_vivian_weiss/4601/4601-CP-68_S70_L005_R2_001.fastq.gz' ],
    #'P4601_CP_69' => [ '/data/h_vivian_weiss/4601/4601-CP-69_S5_L005_R1_001.fastq.gz', '/data/h_vivian_weiss/4601/4601-CP-69_S5_L005_R2_001.fastq.gz' ],
    'P4601_CP_70' => [ '/data/h_vivian_weiss/4601/4601-CP-70_S71_L005_R1_001.fastq.gz', '/data/h_vivian_weiss/4601/4601-CP-70_S71_L005_R2_001.fastq.gz' ],
    'P4601_CP_71' => [ '/data/h_vivian_weiss/4601/4601-CP-71_S72_L005_R1_001.fastq.gz', '/data/h_vivian_weiss/4601/4601-CP-71_S72_L005_R2_001.fastq.gz' ],
    'P4601_CP_72' => [ '/data/h_vivian_weiss/4601/4601-CP-72_S73_L005_R1_001.fastq.gz', '/data/h_vivian_weiss/4601/4601-CP-72_S73_L005_R2_001.fastq.gz' ],
    'P4601_CP_73' => [ '/data/h_vivian_weiss/4601/4601-CP-73_S63_L005_R1_001.fastq.gz', '/data/h_vivian_weiss/4601/4601-CP-73_S63_L005_R2_001.fastq.gz' ],
    #file_def.py -i /data/h_vivian_weiss/6121/ -n "(........\\d+)" -f "*.gz" -a -p
    'P6121_CP_21' => [ '/data/h_vivian_weiss/6121/6121-CP-21_S1_L005_R1_001.fastq.gz', '/data/h_vivian_weiss/6121/6121-CP-21_S1_L005_R2_001.fastq.gz' ],
    'P6121_CP_22' => [ '/data/h_vivian_weiss/6121/6121-CP-22_S1_L005_R1_001.fastq.gz', '/data/h_vivian_weiss/6121/6121-CP-22_S1_L005_R2_001.fastq.gz' ],
    'P6121_CP_23' => [ '/data/h_vivian_weiss/6121/6121-CP-23_S1_L005_R1_001.fastq.gz', '/data/h_vivian_weiss/6121/6121-CP-23_S1_L005_R2_001.fastq.gz' ],
    'P6121_CP_24' => [ '/data/h_vivian_weiss/6121/6121-CP-24_S1_L005_R1_001.fastq.gz', '/data/h_vivian_weiss/6121/6121-CP-24_S1_L005_R2_001.fastq.gz' ],
    'P6121_CP_25' => [ '/data/h_vivian_weiss/6121/6121-CP-25_S1_L005_R1_001.fastq.gz', '/data/h_vivian_weiss/6121/6121-CP-25_S1_L005_R2_001.fastq.gz' ],
    'P6121_CP_26' => [ '/data/h_vivian_weiss/6121/6121-CP-26_S1_L005_R1_001.fastq.gz', '/data/h_vivian_weiss/6121/6121-CP-26_S1_L005_R2_001.fastq.gz' ],
    'P6121_CP_27' => [ '/data/h_vivian_weiss/6121/6121-CP-27_S1_L005_R1_001.fastq.gz', '/data/h_vivian_weiss/6121/6121-CP-27_S1_L005_R2_001.fastq.gz' ],
    'P6121_CP_28' => [ '/data/h_vivian_weiss/6121/6121-CP-28_S1_L005_R1_001.fastq.gz', '/data/h_vivian_weiss/6121/6121-CP-28_S1_L005_R2_001.fastq.gz' ],
    'P6121_CP_29' => [ '/data/h_vivian_weiss/6121/6121-CP-29_S1_L005_R1_001.fastq.gz', '/data/h_vivian_weiss/6121/6121-CP-29_S1_L005_R2_001.fastq.gz' ],
    'P6121_CP_30' => [ '/data/h_vivian_weiss/6121/6121-CP-30_S1_L005_R1_001.fastq.gz', '/data/h_vivian_weiss/6121/6121-CP-30_S1_L005_R2_001.fastq.gz' ],
    'P6121_CP_31' => [ '/data/h_vivian_weiss/6121/6121-CP-31_S1_L005_R1_001.fastq.gz', '/data/h_vivian_weiss/6121/6121-CP-31_S1_L005_R2_001.fastq.gz' ],
    'P6121_CP_32' => [ '/data/h_vivian_weiss/6121/6121-CP-32_S1_L005_R1_001.fastq.gz', '/data/h_vivian_weiss/6121/6121-CP-32_S1_L005_R2_001.fastq.gz' ],
    'P6121_CP_33' => [ '/data/h_vivian_weiss/6121/6121-CP-33_S1_L005_R1_001.fastq.gz', '/data/h_vivian_weiss/6121/6121-CP-33_S1_L005_R2_001.fastq.gz' ],
    'P6121_CP_34' => [ '/data/h_vivian_weiss/6121/6121-CP-34_S1_L005_R1_001.fastq.gz', '/data/h_vivian_weiss/6121/6121-CP-34_S1_L005_R2_001.fastq.gz' ],
    'P6121_CP_35' => [ '/data/h_vivian_weiss/6121/6121-CP-35_S1_L005_R1_001.fastq.gz', '/data/h_vivian_weiss/6121/6121-CP-35_S1_L005_R2_001.fastq.gz' ],
    'P6121_CP_36' => [ '/data/h_vivian_weiss/6121/6121-CP-36_S1_L005_R1_001.fastq.gz', '/data/h_vivian_weiss/6121/6121-CP-36_S1_L005_R2_001.fastq.gz' ],
    'P6121_CP_37' => [ '/data/h_vivian_weiss/6121/6121-CP-37_S1_L005_R1_001.fastq.gz', '/data/h_vivian_weiss/6121/6121-CP-37_S1_L005_R2_001.fastq.gz' ],
    'P6121_CP_38' => [ '/data/h_vivian_weiss/6121/6121-CP-38_S1_L005_R1_001.fastq.gz', '/data/h_vivian_weiss/6121/6121-CP-38_S1_L005_R2_001.fastq.gz' ],
    'P6121_CP_39' => [ '/data/h_vivian_weiss/6121/6121-CP-39_S1_L005_R1_001.fastq.gz', '/data/h_vivian_weiss/6121/6121-CP-39_S1_L005_R2_001.fastq.gz' ],
    'P6121_CP_40' => [ '/data/h_vivian_weiss/6121/6121-CP-40_S1_L005_R1_001.fastq.gz', '/data/h_vivian_weiss/6121/6121-CP-40_S1_L005_R2_001.fastq.gz' ],
  },

  ####################################################################################################
  perform_gatk4_pairedfastq2bam => 0,

  perform_gatk4_refine => 1,

  perform_muTect2indel => 0,

  #muTect1
  perform_muTect => 0,
  muTect_option => '--min_qscore 20 --filter_reads_with_N_cigar',
  mutect_walltime => "22",
  mutect_memory => "40gb",

  #muTect2
  'perform_muTect2' => 1,
  'perform_mutect2_by_wdl' => 0,
  'Mutect2.run_orientation_bias_mixture_model_filter' => "true",
  'mutect2_walltime' => "22",
  'mutect2_memory' => "40gb",

  'perform_mutect2_pon' => 0,
  'pon_samples' => [qw(
P4601_CP_01 
P4601_CP_02
P4601_CP_03
P4601_CP_04
P4601_CP_05
P4601_CP_06
P4601_CP_07
P4601_CP_08
P4601_CP_09
P4601_CP_10
P4601_CP_11
P4601_CP_12
P4601_CP_13
P4601_CP_14
P4601_CP_15
P4601_CP_16
P4601_CP_17
P4601_CP_19
P4601_CP_20
P4601_CP_21
P4601_CP_22
P4601_CP_23
P4601_CP_24
P4601_CP_25
P4601_CP_26
P4601_CP_27
P4601_CP_28
P4601_CP_29
P4601_CP_30
P4601_CP_31
P4601_CP_32
P4601_CP_33
P4601_CP_34
P4601_CP_35
P4601_CP_36
P4601_CP_37
P4601_CP_38
P4601_CP_39
P4601_CP_40
P4601_CP_41
P4601_CP_42
P4601_CP_43
P4601_CP_44
P4601_CP_45
P4601_CP_46
P4601_CP_47
P4601_CP_48
P4601_CP_49
P4601_CP_50
P1896_AC_027
P1809_AC_031
P1809_AC_032
P1809_AC_033
P1809_AC_034
P1809_AC_035
P1809_AC_036
P1809_AC_037
P1809_AC_038
P1809_AC_039
P1809_AC_040
P1809_AC_041
P1809_AC_042
P1809_AC_043
P1809_AC_044
P1809_AC_045
P1809_AC_046
P1809_AC_047
P1809_AC_048
P1809_AC_049
P1809_AC_050
P1809_AC_122
P1809_AC_123
P1809_AC_124
P1809_AC_125
P1809_AC_126
P1809_AC_127
P1809_AC_128
P1809_AC_129
P1809_AC_130
P1809_AC_131
P1809_AC_132
P1809_AC_133
P1809_AC_134
P1809_AC_135
)],

  perform_gatk4_callvariants  => 1,
  perform_check_fastq_duplicate => 0,
  perform_gatk_callvariants => 0,
  perform_cnv_gatk4_cohort => 1,
  callvariants_vqsr_mode => 1,
  perform_IBS => 0,
  perform_multiqc => 0, 
  perform_target_coverage => 0,
  bait_intervals_file => $bait_intervals,
  target_intervals_file => $intervals,
  merge_fastq => 0,
  is_paired_end    => 1,
  perform_cutadapt => 1,
  cutadapt_option  => "-q 20 -a AGATCGGAAGAGCACACGTC -A AGATCGGAAGAGCGTCGTGT",
  min_read_length  => 30,
  covered_bed => $cover_bed,
  perform_cutadapt_validate => 0,

  max_thread => 8,
  bwa_option => "-K 100000000 -v 3",
  bwa_walltime => "22",
  bwa_memory => "40g",
  bwa_use_tmp_folder => 1,
  use_sambamba => 0,

  use_tmp_folder => 0,

  aligner_scatter_count => 0,
  perform_bam_validation => 0,
  bam_validation_walltime => 8,

  perform_cnv_gatk4_somatic => 0,

  perform_CNV_Radar => 0,
  CNVRadar_docker_command => "singularity exec -e /data/cqs/softwares/singularity/cnvradar.v1.2.1.sif",

  perform_cra_gatk4_somatic => 0,
  #default padding is 250. Since our intervals file already has 50 bases padding, we need only 200 in cra.
  cra_padding => 200,
  ref_fasta_main_dict => "/data/cqs/references/broad/hg38/v0/Homo_sapiens_assembly38.main.dict",

  annotation_genes=> "BRAF RAS NTRK2 CDKN2A CDKN2B NF1 KMT2D RB1 MMR ARID2 ATM",
};

if ($def->{target_dir} =~ /workspace/){
  $singularity_prefix_str = "singularity exec -B /workspace,/scratch,/data,/home,/tmp";
}

my $config = performExomeSeq_gatk_hg38( $def, 1 );
#performTask($config, "CNV_Radar_01_roi");
#performTask($config, "gatk4_cra_02_CreateReadCountPanelOfNormals");
#performTask($config, "gatk4_cra_03_DenoiseReadCounts");
#performTask($config, "gatk4_cra_04_PlotDenoisedCopyRatios");
# performTask($config, "gatk4_cra_05_CollectAllelicCounts");
# performTask($config, "gatk4_cra_06_ModelSegments");
# performTask($config, "gatk4_cra_07_CallCopyRatioSegments");
# performTask($config, "gatk4_cra_08_PlotModeledSegments");
#performTask($config, "gatk4_cra_09_MergeModelSegments");

#print(Dumper(keys %$config));


if(0){
  for my $key (sort keys %$config){
    if ($key =~ /_CNV_/) {
      print($key . "\n");
      performTask($config, $key);
    }
  }
  performTask($config, "sequencetask");
}

if(0){
  for my $key (sort keys %$config){
    if ($key =~ /_geneannotation/) {
      print($key . "\n");
      #performTask($config, $key);
    }
  }
  performTask($config, "sequencetask");
}

if(0){
  for my $key (sort keys %$config){
    if ($key =~ /bwa_g4_refine_gatk4_SNV_0/) {
      performTask($config, $key);
    }
  }
  performTask($config, "sequencetask");
}

#performTask($config, "sequencetask");


#performTask($config, 'bwa_g4_refine_gatk4_SNV_05_filter_geneannotation');
#performTask($config, 'bwa_g4_refine_CrosscheckFingerprints');
#performTask($config, "fastqc_raw");
#performTask($config, "cutadapt");
#performTask($config, "fastqc_post_trim");

#$config->{bwa}{target_dir} = $config->{bwa}{target_dir} . ".v2";
#performTask($config, "bwa");
#performTask($config, "bwa_bam_validation");

#$config->{cutadapt_validate}{target_dir} = $config->{cutadapt_validate}{target_dir} . ".v3";
#performTask($config, "cutadapt_validate");
#performTask($config, "bwa_g4_refine");
#performTask($config, "bwa_g4_refine_bam_validation");
#performTask($config, "PON_muTect2_01_call");
#performTask($config, "PON_muTect2_01_call_wdl");
#performTask($config, "PON_muTect2_02_table");
#performTask($config, "muTect2_01_call");
#performTask($config, "muTect2_01_call_wdl");
#performTask($config, "muTect2_04_annovar");

1;
