from Bio import AlignIO
import os
import click
import logging
from rich.logging import RichHandler
import matplotlib.pyplot as plt
import numpy as np

@click.group(context_settings={'help_option_names': ['-h', '--help']})
@click.option('--verbosity', '-v', type=click.Choice(['info', 'debug']), default='info', help="Verbosity level, default = info.")
def cli(verbosity):
    """
    This script contains multiple functions for the CoBi course
    """
    logging.basicConfig(
        format='%(message)s',
        handlers=[RichHandler()],
        datefmt='%H:%M:%S',
        level=verbosity.upper())
    pass

def basicalninfo(aln):
    Align_obj = AlignIO.read(aln, "fasta")
    total_length = Align_obj.get_alignment_length()
    gap_ncol = 0;iden_ncol = 0
    for col in range(total_length):
        col_seq = Align_obj[:,col]
        if '-' in col_seq: gap_ncol+=1
        elif len(set(col_seq))==1: iden_ncol+=1
    return total_length,gap_ncol,iden_ncol

@cli.command(context_settings={'help_option_names': ['-h', '--help']})
@click.argument('aln', type=click.Path(exists=True))
def getalninfo(aln):
    """
    get basic msa information
    """
    total_length,gap_ncol,iden_ncol = basicalninfo(aln)
    print("Gap column: {}/{:.2f}".format(gap_ncol,gap_ncol/total_length))
    print("Identical column: {}/{:.2f}".format(iden_ncol,iden_ncol/total_length))

@cli.command(context_settings={'help_option_names': ['-h', '--help']})
@click.argument('alns', nargs=-1, type=click.Path(exists=True))
@click.option('--output', '-o', default='MSA.summary.svg', show_default=True, help='output filename')
def summarizealns(alns,output):
    gap_ncols_p = [];iden_ncols_p = [];total_lengths = []
    for aln in alns:
        total_length,gap_ncol,iden_ncol = basicalninfo(aln)
        gap_ncols_p += [gap_ncol/total_length]
        iden_ncols_p += [iden_ncol/total_length]
        total_lengths += [total_length]
    fig, ax = plt.subplots(1,1,figsize=(4, 4))
    ax.scatter(total_lengths,gap_ncols_p,s=20,marker='o',color='r',label='Gap percentage')
    ax.hlines(np.mean(gap_ncols_p),np.min(total_lengths),np.max(total_lengths),color='r', linestyle='--', lw=2,label='Mean gap percentage')
    ax.scatter(total_lengths,iden_ncols_p,s=20,marker='o',color='b',label='Identity percentage')
    ax.hlines(np.mean(iden_ncols_p),np.min(total_lengths),np.max(total_lengths),color='b', linestyle='--', lw=2,label='Mean identity percentage')
    ax.set_xlabel("MSA length",fontsize=10)
    ax.set_ylabel("Gap/Identity percentage",fontsize=10)
    ax.grid(ls=":")
    ax.legend(loc=0,fontsize=10,frameon=False)
    fig.tight_layout()
    fig.savefig(output,format ='svg', bbox_inches='tight')
    plt.close()


if __name__ == '__main__':
	cli()
