import sys
import numpy as np
import deeptools.countReadsPerBin as crpb
from collections import namedtuple
from hexcaller.src.utils import *

"""call_he.py: Given a coverage file from either chromosome pairs or from a sample versus reference
comparison, detect bias in CNV calls."""


def load_he(filename):
    """Load HE gene pairs file
    args:
        filename: filename to load

    returns:
        dict: dictionary of HE pairs
    """
    he_pairs = {}

    with open(filename, 'r') as infile:
        for line in infile:
            dat = line.strip().split('\t')
            if len(dat) == 8:
                if dat[0] in covs.keys():
                    covs[dat[0]].append({
                        'chr': dat[1],
                        'position': int(dat[2]),
                        'scoreA': float(dat[4]),
                        'scoreB': float(dat[6]),
                        'info': line.strip()
                    })
                else:
                    covs[dat[0]] = [{
                        'chr': dat[1],
                        'position': int(dat[2]),
                        'scoreA': float(dat[4]),
                        'scoreB': float(dat[6]),
                        'info': line.strip()
                    }]
    return covs



def is_hit(he_array, i, seed_size):
    """Determine if an he_array contains a hit in (i,i+seed_size) index pos.
    args:
        he_array (list): array of he_event strings
        i: start index
        seed_size: range to scan

    returns:
        boolean: return if a hit is found
    """

    if he_array[i] == "norm/norm":
        return False
    for j in range(i, i+seed_size):
        #print(f"j is {j}")
        if he_array[i] != he_array[j]:
            return False
    return True


def call_he_type(zscore_A, zscore_B):
    """Given a pair of z-scores return a CNV event type
    args:
        zscore_A (float): z-score observed in the first sample
        zscoreB (float): z-score observed in the reference

    returns:
        string: return CNV event type
    """

    cutoff = 1
    if zscore_A < -1*cutoff and zscore_B < -1*cutoff:
        return "del/del"
    if zscore_A < -1*cutoff and abs(zscore_B) <= cutoff:
        return "del/norm"
    if zscore_A < -1*cutoff and zscore_B > cutoff:
        return "del/dup"
    if abs(zscore_A) <= cutoff and zscore_B < -1*cutoff:
        return "norm/del"
    if abs(zscore_A) <= cutoff and abs(zscore_B) <= cutoff:
        return "norm/norm"
    if abs(zscore_A) <= cutoff and zscore_B > cutoff:
        return "norm/dup"
    if zscore_A > cutoff and zscore_B < -1*cutoff:
        return "dup/del"
    if zscore_A > cutoff and abs(zscore_B) <= cutoff:
        return "dup/norm"
    if zscore_A > cutoff and zscore_B > cutoff:
        return "dup/dup"
    return "undefined"


def normalize(covs):
    """Given raw coverage input normalize data and detect CNV events.
    args:
        covs: raw coverage loaded from input file

    """

    for sample in covs.keys():
        z_scores = normalize_sample(sample, covs[sample])
        call_he_events(sample, z_scores, covs[sample])

    #one-sided filtering of z-score to remove extreme high cov.
    #filtered_entries = (z_scores < 3).all(axis=1)
    #disable for now.
    #new_df = df[filtered_entries]
    #print(f"{data}")
    #print(f"{z_scores}")


def normalize_sample(sample, sample_cov):
    """Given raw coverage of a sample return normalized values
    args:
        sample_cov: raw coverage of a sample

    returns:
        array: return an array of normalized values
    """

    gene_A_covs = []
    gene_B_covs = []
    log_message(f"normalizing sample: {sample}")

    for gene in sample_cov:
        gene_A_covs.append(gene['scoreA'])
        gene_B_covs.append(gene['scoreB'])
    observed_covs = [gene_A_covs, gene_B_covs]
    data = np.array(observed_covs)
    mean = np.mean(data)
    stdev = np.std(data)
    log_message(f"sample= {sample}, mean= {mean:.2f}, stdev= {stdev:.2f}")
    z_scores = (data - mean) / stdev

    return z_scores


def call_he_events(sample, z_scores, sample_cov):
    gene_count = len(z_scores[0])
    he_array = []
    for i in range(gene_count):
        zscore_A = z_scores[0][i]
        zscore_B = z_scores[1][i]
        he_type = call_he_type(zscore_A, zscore_B)
        he_array.append(he_type)
    seeds = find_seeds(he_array)
    he_events = extend_seeds(seeds, sample_cov)
    summarize_cnvs(sample, he_events, z_scores, sample_cov)


def summarize_cnvs(sample, he_events, z_scores, sample_cov):
    summary = []
    cur_event = {}
    cur_event_type = 'undef'
    for i, gene in enumerate(sample_cov):
        #generate event summary
        if len(he_events[i]):
            if cur_event_type == he_events[i]:
                # update end position, gene_counts
                cur_event['end'] = gene['position']
                cur_event['count'] += 1
            else:
                cur_event_type = he_events[i]
                cur_event = {
                    'sample': sample,
                    'chr': gene['chr'],
                    'event': cur_event_type,
                    'start': gene['position'],
                    'end': gene['position'],
                    'count': 1
                }
        else:
            if len(cur_event):
                summary.append(cur_event)
            cur_event = []
            cur_event_type = 'undef'
    print_cnv_summary(sample, summary)
    print_event(sample, he_events, z_scores, sample_cov)


def print_cnv_summary(sample, summary):
    summary_file = open(sample+".hexcaller.summary.txt","w")
    for cnv in summary:
        if cnv['end'] - cnv['start'] > 20000:
            summary_file.write(f"{cnv}")


def find_seeds(he_array):
    seed_size = 5
    gene_count = len(he_array) - seed_size
    he_seeds = [''] * len(he_array)
    for i, he_type in enumerate(he_array):
        #print(f"{i} {he_type}")
        if i <= gene_count:
            if is_hit(he_array, i, seed_size):
                #print(f"{i} is hit")
                for j in range(i,i+seed_size):
                    he_seeds[j] = he_type
#    print(f"{he_events}")
    return he_seeds

def print_event(sample, he_events, z_scores, sample_cov):
    event_file = open(sample+".hexcaller.txt","w")
    header_line = ",".join(['sample','group','pos','geneA','covA','geneB','covB','ZA','ZB','HeEvent'])
    event_file.write(header_line+"\n")

    for i, gene in enumerate(sample_cov):
        event_file.write(f"{gene['info']},{z_scores[0][i]},{z_scores[1][i]},{he_events[i]}\n")


def extend_seeds(he_array, sample_cov):
    #skip
    return he_array

def run_bedtools(he_pair, bam_file):
    subprocess.run(["bedtools", "coverage", "-mean", "-a "+he_pair,
                    "-b "+bam_file], stdout=f)

def load_covs_from_bam(he_gene_pair, input_bam_file):
    """Load cov from bam
    """
    log_message(f"open {input_bam_file}")
    cr = crpb.CountReadsPerBin([input_bam_file],binLength=1, stepSize=1, numberOfProcessors=10)

    headers = ['chrA','g1_start', 'g1_end', 'g1_name', 'chrB', 'g2_start', 'g2_end', 'g2_name']
    covs = {}
    num_processed = 0
    sample_name = input_bam_file.rpartition(".bam")[0]
    log_message(f"assign sample_name: {sample_name}")
    covs[sample_name] = []
    with open(he_gene_pair, 'r') as infile:
        for line in infile:
            dat = line.strip().split('\t')
            if len(dat) == 8:
                he_pair = dict(zip(headers,dat))
                g1_cov = cr.count_reads_in_region(he_pair['chrA'],int(he_pair['g1_start']),int(he_pair['g1_end']))
                scoreA = np.mean(g1_cov[0])
                g2_cov = cr.count_reads_in_region(he_pair['chrB'],int(he_pair['g2_start']),int(he_pair['g2_end']))
                scoreB = np.mean(g2_cov[0])

                covs[sample_name].append({
                        'chr': dat[0],
                        'position': int(dat[1]),
                        'scoreA': scoreA,
                        'scoreB': scoreB,
                        'info': ','.join([sample_name, he_pair['chrA'], he_pair['g1_start'],
                                          he_pair['g1_name'],str(scoreA), he_pair['g2_name'], str(scoreB)])
                    })
                num_processed += 1
                if num_processed % 1000 == 0:
                    log_message(f"calculated mean cov for {num_processed} pairs")
            else:
                log_error("invalid HE file format: expected 8-col tab-delimited file.")
    return covs


def call_cnv_from_cov(he_gene_pair, input_bam_file):
    covs = load_covs_from_bam(he_gene_pair, input_bam_file)
    #covs = load_covs(input_cov_file)
    norm_cov = normalize(covs)

