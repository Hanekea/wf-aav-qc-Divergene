#!/usr/bin/env python
"""Create reference with the various input sequences."""

import json
from pathlib import Path
import subprocess

import numpy as np
import pandas as pd

from .util import wf_parser  # noqa: ABS101


def argparser():
    """Create argument parser."""
    parser = wf_parser("contam")

    parser.add_argument(
        '--bam_info',
        help="The output of seqkit bam",
        type=Path)
    parser.add_argument(
        '--sample_id',
        help="sample ID")
    parser.add_argument(
        '--transgene_fasta',
        help="The transgene plasmid fasta sequence")
    parser.add_argument(
        '--ref_ids',
        help="Json file with reference:contig mappings")
    parser.add_argument(
        '--n_reads',
        help="Total number of input reads", type=int)
    parser.add_argument(
        '--contam_class_counts', help="Contamination counts output", type=Path)

    return parser


def main(args):
    """Run main entry point."""
    # Extract the reference sequence names from the FATA files,
    # seqkit adds a newline to the output, so remove that with `rstrip`
    transgene_plasmid_name = subprocess.check_output(
        ['seqkit', 'seq', '--name', '--only-id', args.transgene_fasta],
        encoding='UTF-8').rstrip()

    with open(args.ref_ids, 'r') as json_fh:
        ref_ids = json.load(json_fh)

    # Read the per-alignment read summaries
    df_bam = pd.read_csv(
        args.bam_info,
        sep='\t',
        usecols=['Read', 'Ref', 'ReadLen'],
        dtype={
            'Read': str,
            'Ref': str,
            'ReadLen': np.uint32
            }
    )
    # Assign reference category to alignments
    df_bam['contam_class'] = None
    for ref, contigs in ref_ids.items():
        df_bam.loc[df_bam.Ref.isin(contigs), 'contam_class'] = ref
    df_bam.loc[df_bam.Ref == transgene_plasmid_name, 'contam_class'] = 'Transgene'

    # Calculate mapped/unmapped as a percentage of reads not alignments.
    # Note `seqkit bam` does not return info for unammoed reads, so we need to get the
    # total number of input reads separately.
    n_mapped_reads = df_bam.Read.nunique()
    n_input_reads = args.n_reads
    n_unmapped_reads = n_input_reads - n_mapped_reads
    unmapped_pct = 100 / n_input_reads * n_unmapped_reads
    mapped_pct = 100 - unmapped_pct

    df_contam_class = pd.DataFrame(df_bam['contam_class'].value_counts())
    df_contam_class = df_contam_class.rename(columns={'count': 'Number of alignments'})
    df_contam_class['Percentage of alignments'] = (
            100 / len(df_bam) * df_contam_class['Number of alignments']).round(2).T
    df_contam_class.index.name = 'Reference'
    df_contam_class.loc['Unmapped'] = [n_unmapped_reads, unmapped_pct]
    df_contam_class.loc['Mapped'] = [n_mapped_reads, mapped_pct]
    df_contam_class['sample_id'] = args.sample_id

    df_lengths = df_bam[['ReadLen', 'contam_class']]
    df_lengths['sample_id'] = args.sample_id

    df_contam_class.to_csv(args.contam_class_counts, sep='\t')

    '''START OF NEW MODIFICATIONS'''
    #map reads to contam class it aligns with 
    read_to_classes = df_bam.groupby('Read')['contam_class'].apply(set)
    vector = 'Transgene'

    map_tovector = 0
    map_tononvector = 0
    map_toboth = 0 

    for classes in read_to_classes:
        has_vector = vector in classes
        has_nonvector = any(c != map_tovector for c in classes if pd.notna(c))

        if has_vector and not has_nonvector:
            map_tovector += 1
        elif not has_vector and has_nonvector:
            map_tononvector += 1
        elif has_vector and has_nonvector:
            map_toboth += 1
        else:
            continue

    total_reads_classified = map_tovector + map_tononvector + map_toboth

    #creates dataframe for output file
    df_vector_class = pd.DataFrame({
        'Reference': ['Map to Vector', 'Map to Non-Vector', 'Map to Both'],
        'Number of Alignments': [map_tovector, map_tononvector, map_toboth],
        'Alignment Percentages': [
            100 * map_tovector / total_reads_classified if total_reads_classified else 0,
            100 * map_tononvector / total_reads_classified if total_reads_classified else 0,
            100 * map_toboth / total_reads_classified if total_reads_classified else 0,
        ],
        'sample_id': args.sample_id
    })

    #write to TSV file
    vector_out_path = args.contam_class_counts.parent / f"{args.sample_id}_vector_vs_nonvector.tsv"
    df_vector_class.to_csv(args.vector_out_path, sep='\t', index=False)

    print(f"Vector/non-vector classification saved to {vector_out_path}")
