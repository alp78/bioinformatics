"""
Reference files:
FASTA: https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz
GTF: ftp://ftp.ensembl.org/pub/release-99/gtf/homo_sapiens/Homo_sapiens.GRCh38.99.chr.gtf.gz
"""

import pandas as pd
from Bio import SeqIO
import mygene
from gtfparse import read_gtf
import pybedtools
import logging
from os import path, remove
from sys import exit
from termcolor import colored

def make_utr_tables(gtf, _types, logger):
    logger.info(colored('Making UTRs tables...', 'blue'))
    feature = ''
    gtf_df = read_gtf(gtf)
    for _type in _types:
        if _type == 5:
            feature = 'five_prime_utr'
        if _type == 3:
            feature = 'three_prime_utr'
        utr_df = gtf_df[gtf_df['feature']==feature]
        utr_unique = utr_df.drop_duplicates(subset=['gene_id', 'source'], keep='first')
        utr_unique.reset_index(drop=True, inplace=True)
        utr_unique = utr_unique.loc[:, ['source', 'gene_id', 'feature', 'seqname', 'strand', 'start', 'end']]
        utr_unique.to_csv(f'utr{_type}.txt', sep='\t', index = False)

def convert_id(_id, logger):
    logger.info(colored('Converting ID...', 'blue'))
    mg = mygene.MyGeneInfo()
    entrez_id = mg.query(_id, species='human')['hits'][0]['entrezgene']
    gene_symbol = mg.query(entrez_id, species='human')['hits'][0]['symbol']
    ensembl_id = mg.getgene(entrez_id, 'ensembl')['ensembl']['gene']
    return ensembl_id, gene_symbol

def get_utr(ensembl_id, gene_symbol, _type, logger):
    logger.info(colored(f'Getting {_type}-UTR...', 'blue'))
    utr_unique = pd.read_csv(f'utr{_type}.txt', sep='\t', skiprows=(0), header=(0), low_memory=False)
    choices = utr_unique[utr_unique['gene_id']==ensembl_id]['source'].to_list()
    if len(choices) > 0:
        temp_df = utr_unique[utr_unique['gene_id']==ensembl_id]
        if 'ensembl' in choices:
            location = temp_df[temp_df['source']=='ensembl'].values.tolist()[0][-6:]
        elif 'ensembl_havana' in choices:
            location = temp_df[temp_df['source']=='ensembl_havana'].values.tolist()[0][-6:]
        else:
            location = temp_df[temp_df['source']=='havana'].values.tolist()[0][-6:]
        location[2] = 'chr' + location[2]
        location[0] = gene_symbol
    return location

def make_utrs_fasta(ref_fasta_full_path, location_list, logger):
    logger.info(colored('Cooking FASTA file...', 'blue'))
    recs_list = []
    for location in location_list:
        strand = location[3]
        loc = pybedtools.BedTool(f"""{location[2]} {location[4]} {location[5]}""", from_string=True)
        fasta = pybedtools.example_filename(ref_fasta_full_path)
        seq = loc.sequence(fi=fasta)
        out_fasta = seq.save_seqs('temp_UTR.fasta')
        for rec in SeqIO.parse('temp_UTR.fasta', 'fasta'):
            seq = rec.seq.upper()
            if strand == '-':
                rev_seq = seq.reverse_complement()
                rec.seq = rev_seq
            rec.id = f'{location[0]}|{location[1]}|{location[2]}|{location[3]}|{location[4]}|{location[5]}'
            rec.description = ''
            recs_list.append(rec)
    recs_gen = (rec for rec in recs_list)
    SeqIO.write(recs_gen, 'UTR.fasta', "fasta")
    remove('temp_UTR.fasta')

# PROGRAM START
def get_utrs(gtf, ref_fasta_full_path, _id):
    logger = logging.getLogger()
    logger.setLevel(logging.INFO)

    location_list = []
    ensembl_id = ''
    gene_symbol = ''
    utr3 = []
    utr5 = []

    if not (path.isfile('utr3.txt') and path.isfile('utr5.txt')):
        try:
            types = [5,3]
            make_utr_tables(gtf, types, logger)
        except:
            logger.info(colored('Error in parsing GTF file', 'red'))
            del logger
            exit(0)

    try:
        ensembl_id, gene_symbol = convert_id(_id, logger)
    except:
        logger.info(colored('Error in converting ID', 'red'))
        del logger
        exit(0)

    if ensembl_id:
        try:
            utr3 = get_utr(ensembl_id, gene_symbol, 3, logger)
        except:
            logger.info(colored('Error in retrieving 3-UTR', 'red'))
            del logger
            exit(0)
        try:
            utr5 = get_utr(ensembl_id, gene_symbol, 5, logger)
        except:
            logger.info(colored('Error in retrieving 5-UTR', 'red'))
            del logger
            exit(0)
    else:
        logger.info(colored('Error in converting ID', 'red'))
        del logger
        exit(0)

    location_list.append(utr5)
    location_list.append(utr3)

    if (len(location_list) == 2 and len(utr3) == 6 and len(utr5) == 6):
        try:
            make_utrs_fasta(ref_fasta_full_path, location_list, logger)
        except:
            logger.info(colored('Error in making fasta file', 'red'))
            del logger
            exit(0)
    else:
        logger.info(colored('Error in location list', 'red'))
        del logger
        exit(0)

    logger.info(colored('All done!', 'green'))
    for location in location_list:
        print(location)
    del logger
