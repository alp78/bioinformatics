import os
import imp
import time
import datetime as dt
from os import listdir
from os.path import isfile, join
from datetime import datetime
from pathlib import Path
from multiprocessing import cpu_count, Manager
from concurrent.futures import ProcessPoolExecutor
import trimmRE
import collapser
import samfilter_pysam
import bwamem
import exactmatch_fpALU
import intersection_fpALU
import misseq
import intersection_replib
import trash_deleting
import metacluster
from osmgmt import concat_logs_and_stats

MAPPER = 'bowtie' #bwa or bowtie
NCORE = 1

#PREPROCESSING
PRIMER = 'AGCCACCGCGC'
SHIFT = 5
MIST = 1
AD1 = 'CGTGCTGCGG'
AD2 = 'AGGCGTCGTGCG'
BLEN = 10
RE = 'CCGGCC'
RESTRICT_SITE = {'TCGA': 'CGA', 'CTAG': 'TAG'}
IS_SHORT_FLANK = False
CHAOS = True #False
TRIMM_N = 0
TRIMM_POLY_N = False
POLY_N_R1 = () # poly_n_win_size_r1, poly_n_th_r1, poly_n_shift_r1
POLY_N_R2 = () # poly_n_win_size_r2, poly_n_th_r2, poly_n_shift_r2
SKIP_SHORT_READS = 50
MID_MIST_SHORT_READS = (3, 2) # (r1 ,r2) max mismatches in AD2 or PRIMER which will remove from reads
END_MIST_SHORT_READS = (6, 6) # (r1 ,r2) min length of AD2 or PRIMER at the reads ends
PLACE_OF_SEARCH_TAIL = (None, None) # (r1 ,r2)
MIN_SEQ_LEN_AFTER_TRIMM = (25, 25) # (r1 ,r2)

#SAMFILTER
FLAGS = [99, 83, 147, 163]# + [355, 339, 403, 419]

R_SITE_DISTANCE = 1000

#FILTERS
FIX_WINDOW = 20
MIN_READ = 1
MAX_DIST = 1000
MIN_DIST = 0
INSWINDOW_FIX = 40
SHIFT_RESTRICT_SITE = 3 
SHIFT_MISS_PRIMER = 18
SHIFT_MISS_RE = 8

#THRESHOLDS
RE_HAMMING = 2
FLANK_ERRORS = 5
REPEAT = 0.8
MISS_PRIMER = 4
MISS_RE = 2
TEMPLATE_SWITCH_MD = 30

#METACLUSTERING
WINDOW = 100
MAIN_FLANK_LEN = 11

INPUTDIR = '/mnt/nfs/home/244205/000000-My_Documents/VM-home/alu_snakemake/raw_fastq'
OUTPUTDIR = '/mnt/nfs/home/244205/000000-My_Documents/VM-home/alu_snakemake/output'

out_preproc = os.path.abspath(OUTPUTDIR) + '/preprocessing/'
out_map = os.path.abspath(OUTPUTDIR) + '/mapping/'
out_table = os.path.abspath(OUTPUTDIR) + '/table/'
out_collapse = os.path.abspath(OUTPUTDIR) + '/collapse_table/'
out_ematch = os.path.abspath(OUTPUTDIR) + '/ematch_table/'
out_notfp = os.path.abspath(OUTPUTDIR) + '/notfpALU_table/'
out_fix = os.path.abspath(OUTPUTDIR) + '/fix_ins_table/'
out_filter = os.path.abspath(OUTPUTDIR) + '/filter_rsite/'
out_primer = os.path.abspath(OUTPUTDIR) + '/filter_primer/'
out_re = os.path.abspath(OUTPUTDIR) + '/filter_re/'
out_replib = os.path.abspath(OUTPUTDIR) + '/filter_replib/'
out_premeta = os.path.abspath(OUTPUTDIR) + '/pre_metatable/'
out_meta = os.path.abspath(OUTPUTDIR) + '/metatable/'

# unique names of samples prefixes
SAMPLES = list(set([filename.split('.')[0].split('_')[0] for filename in listdir(INPUTDIR) if filename[-2:]=='gz']))

full_start = time.time()

rule all:
    input:
        meta_out = expand('%smetatable_{switch}read.txt' % out_meta, switch = ['human', 'pc'])

rule preprocess_fastq:
    input:
        input_r1 = expand('%s/{sample}_{lane}.fastq.gz' % INPUTDIR, sample = SAMPLES, lane ='R1'),
        input_r2 = expand('%s/{sample}_{lane}.fastq.gz' % INPUTDIR, sample = SAMPLES, lane ='R2')
    output:
        fastq = expand('%s{sample}_{lane}_{switch}.fastq' % out_preproc, sample = SAMPLES, lane = ['R1', 'R2'], switch = ['good', 'bad']),
        meta = expand('%s{sample}_meta.txt' % out_preproc, sample = SAMPLES),
        logs = '%slogfile_trim.log' % out_preproc,
        stats = '%sstatistics_trim.csv' % out_preproc     
    run:
        imp.reload(trimmRE)
        combo_list = list(zip(sorted(list(input.input_r1)), sorted(list(input.input_r2))))
        manager = Manager()
        folders_q = manager.Queue()
        number_of_cpus = cpu_count()
        for d in combo_list:
            folders_q.put(d)
        start = time.time()
        print('')
        print('Trimming and sorting reads...')
        with ProcessPoolExecutor(max_workers=number_of_cpus) as step1:
            for in_dir in range(folders_q.qsize()):
                step1.submit(trimmRE.main, folders_q, INPUTDIR, out_preproc, PRIMER, AD1, AD2, BLEN, SHIFT, MIST, RESTRICT_SITE, RE,
                            IS_SHORT_FLANK, CHAOS, NCORE, TRIMM_N, TRIMM_POLY_N, POLY_N_R1, POLY_N_R2, SKIP_SHORT_READS,
                            MID_MIST_SHORT_READS, END_MIST_SHORT_READS, PLACE_OF_SEARCH_TAIL, MIN_SEQ_LEN_AFTER_TRIMM)
        concat_logs_and_stats(out_preproc)
        end = time.time()
        print('')
        print('Done.')
        print(f'{str(dt.timedelta(seconds=int(end - start)))}')
        print('')


rule map_fastq:
    input:
        input_r1 = expand('%s{sample}_{lane}_good.fastq' % out_preproc, sample = SAMPLES, lane ='R1'),
        input_r2 = expand('%s{sample}_{lane}_good.fastq' % out_preproc, sample = SAMPLES, lane ='R2')
    output:
        sams = expand('%s{sample}.sam' % out_map, sample = SAMPLES),
        logs = '%slogfile_map.log' % out_map
    run:
        imp.reload(bwamem)
        combo_list = list(zip(sorted(list(input.input_r1)), sorted(list(input.input_r2))))
        if MAPPER == 'bowtie':
            mapper_execline = 'bowtie2 -p 4'
            refway = '/mnt/nfs/shared/999993-Bioda/projects/kaja_transposons/reference/bowtie2_index/hg38'
        else:
            if MAPPER != 'bwa':
                mapper_execline = './bwa-0.7.13/bwa mem -t 4'
                refway = '/mnt/nfs/shared/999993-Bioda/projects/kaja_transposons/reference/bowtie2_index/hg38'
        bwamem.main(outputdir = out_map,
                    refway = refway,
                    bwaline = mapper_execline,
                    conform_files = combo_list)


rule filter_sam:
    input:
        input_sam = expand('%s{sample}.sam' % out_map, sample = SAMPLES)
    output:
        text_out = expand('%s{sample}.txt' % out_table, sample = SAMPLES)
    run:
        imp.reload(samfilter_pysam)
        sorted_sams = sorted(list(input.input_sam))
        samfilter_pysam.main(inputdir = out_map,
                            outputdir = out_table,
                            flags = FLAGS,
                            n_core = NCORE,
                            onlyfiles = sorted_sams)      


rule collapse_reads:
    input:
        in_tables = expand('%s{sample}.txt' % out_table, sample = SAMPLES)
    output:
        collapse_reads = expand('%s{sample}_{switch}read.txt' % out_collapse, sample = SAMPLES, switch = ['human', 'pc'])
    run:
        imp.reload(collapser)
        filtered_sams = sorted(list(input.in_tables))
        collapser.main(outputdir = out_collapse,
                    target_re = RE,
                    n_core = 1,
                    onlyfiles = filtered_sams)  


rule exact_match_and_remove_poly_fix:
    input:
        in_collapse = expand('%s{sample}_{switch}read.txt' % out_collapse, sample = SAMPLES, switch = 'human')
    output:
        ematch_out = expand('%s{sample}_ematch.txt' % out_ematch, sample = SAMPLES),
        notfp_out = expand('%s{sample}.txt' % out_notfp, sample = SAMPLES),
        fix_out = expand('%s{sample}.txt' % out_fix, sample = SAMPLES)
    run:
        imp.reload(exactmatch_fpALU)
        sorted_collapse = sorted(list(input.in_collapse))
        exactmatch_fpALU.main(inputdir = out_collapse,
                            outputdir = out_ematch,
                            replibrary = '/mnt/nfs/shared/999993-Bioda/projects/kaja_transposons/moscow/repdata/Alu_replibrary_hg38.txt',
                            refway = '/mnt/nfs/shared/999993-Bioda/projects/kaja_transposons/reference/hg38.fa',
                            restrict_site = RESTRICT_SITE,
                            max_dist = MAX_DIST,
                            min_dist = 100,
                            min_read = MIN_READ,
                            inswindow = FIX_WINDOW,
                            direction = 5,
                            n_core = NCORE,
                            onlyfiles = sorted_collapse)

        imp.reload(intersection_fpALU)
        sorted_collapse = sorted(list(input.in_collapse))
        intersection_fpALU.main(inputdir = out_collapse,
                                outputdir = out_notfp,
                                outputdir_fix = out_fix,
                                fix_ins = ["AluYa5"],
                                replib_inputdir = out_ematch,
                                inswindow = INSWINDOW_FIX,
                                n_core = NCORE,
                                onlyfiles = sorted_collapse)


rule find_restrict_site:
    input:
        notfp_in = expand('%s{sample}.txt' % out_notfp, sample = SAMPLES)
    output:
        filter_out = expand('%s{sample}.txt' % out_filter, sample = SAMPLES)
    run:
        imp.reload(misseq)
        sorted_notfp = sorted(list(input.notfp_in))
        print('')
        print('Finding restrict site..')
        misseq.main(inputdir = out_notfp,
                    outputdir = out_filter,
                    refway = '/mnt/nfs/shared/999993-Bioda/projects/kaja_transposons/reference/hg38.fa',
                    mseq = RESTRICT_SITE,
                    mname = "R_SITE",
                    shift = SHIFT_RESTRICT_SITE,
                    n_core = NCORE,
                    onlyfiles = sorted_notfp)
        

rule find_flank_primer:
    input:
        filter_in = expand('%s{sample}.txt' % out_filter, sample = SAMPLES)
    output:
        primer_out = expand('%s{sample}.txt' % out_primer, sample = SAMPLES)
    run:
        imp.reload(misseq)
        sorted_filter = sorted(list(input.filter_in))
        print('')
        print('Finding primer in flank..')
        misseq.main(inputdir = out_filter,
                    outputdir = out_primer,
                    refway = '/mnt/nfs/shared/999993-Bioda/projects/kaja_transposons/reference/hg38.fa',
                    mseq = PRIMER,
                    mname = 'MISS_P_HAMMING',
                    shift = SHIFT_MISS_PRIMER,
                    n_core = NCORE,
                    onlyfiles = sorted_filter)


rule find_flank_re:
    input:
        primer_in = expand('%s{sample}.txt' % out_primer, sample = SAMPLES)
    output:
        re_out = expand('%s{sample}.txt' % out_re, sample = SAMPLES)
    run:
        imp.reload(misseq)
        sorted_primer = sorted(list(input.primer_in))
        print('')
        print('Finding part of RE in flank..')
        misseq.main(inputdir = out_primer,
                    outputdir = out_re,
                    refway = '/mnt/nfs/shared/999993-Bioda/projects/kaja_transposons/reference/hg38.fa',
                    mseq = None,
                    mname = 'MISS_RE_HAMMING',
                    shift = SHIFT_MISS_RE,
                    n_core = NCORE,
                    onlyfiles = sorted_primer)


rule intersect_repeats:
    input:
        re_in = expand('%s{sample}.txt' % out_re, sample = SAMPLES)
    output:
        replib_out = expand('%s{sample}.txt' % out_replib, sample = SAMPLES)
    run:
        imp.reload(intersection_replib)
        sorted_re = sorted(list(input.re_in))
        print('')
        print('Intersecting reads with repeats...')
        intersection_replib.main(inputdir = out_re,
                                outputdir = out_replib,
                                repeatway = '/mnt/nfs/shared/999993-Bioda/projects/kaja_transposons/moscow/repdata/repeats_hg38.tabular',
                                n_core = 1,
                                onlyfiles = sorted_re)


rule clean_dataset:
    input:
        replib_in = expand('%s{sample}.txt' % out_replib, sample = SAMPLES)
    output:
        premeta_out = expand('%s{sample}.txt' % out_premeta, sample = SAMPLES)
    run:
        imp.reload(trash_deleting)
        sorted_replib = sorted(list(input.replib_in))
        trash_deleting.main(inputdir = out_replib,
                    outputdir = out_premeta,
                    re_hamming = RE_HAMMING,
                    flank_errors = FLANK_ERRORS,
                    rsite = 'R_SITE',
                    repeat = REPEAT,
                    m_primer = MISS_PRIMER,
                    primer_name = 'P',
                    m_re = MISS_RE,
                    re_name = 'RE',
                    onlyfiles = sorted_replib)


rule meta_cluster:
    input:
        premeta_in = expand('%s{sample}.txt' % out_premeta, sample = SAMPLES),
        collapse_in = expand('%s{sample}_{switch}read.txt' % out_collapse, sample = SAMPLES, switch='pc')
    output:
        meta_out = expand('%smetatable_{switch}read.txt' % out_meta, switch = ['human', 'pc'])
    run:
        imp.reload(metacluster)
        sorted_premeta = sorted(list(input.premeta_in))
        sorted_pc = sorted(list(input.collapse_in))
        metacluster.main(inputdir = out_premeta,
                        outputdir = out_meta,
                        pcdir = out_collapse,
                        target_re = RE,
                        window = WINDOW,
                        blen = BLEN,
                        onlyfiles = sorted_premeta,
                        pcfiles = sorted_pc)
        full_end = time.time()
        print('')
        print('Full process terminated.')
        print(f'{str(dt.timedelta(seconds=int(full_end - full_start)))}')
        print('')

# snakemake --dag | dot -Tsvg > dag.svg