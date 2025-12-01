from Bio import AlignIO
import os
import click
import logging
from rich.logging import RichHandler


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
    logging.info("Proper Initiation")
    pass

@cli.command(context_settings={'help_option_names': ['-h', '--help']})
@click.argument('aln', type=click.Path(exists=True))
def getalninfo(aln):
    """
    get basic msa information
    """
    Align_obj = AlignIO.read(aln, "fasta")
    total_length = Align_obj.get_alignment_length()
    gap_ncol = 0
    iden_ncol = 0
    for col in range(total_length):
        col_seq = Align_obj[:,col]
        if '-' in col_seq: gap_ncol+=1
        elif len(set(col_seq))==1:
            iden_ncol+=1
    logging.info("Gap column: {}/{:.2f}".format(gap_ncol,gap_ncol/total_length))
    logging.info("Identical column: {}/{:.2f}".format(iden_ncol,iden_ncol/total_length))

if __name__ == '__main__':
	cli()
