from Bio import AlignIO,Phylo
import os
import click
import logging
from rich.logging import RichHandler
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

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
@click.option('--output', '-o', default='msa.summary.svg', show_default=True, help='output filename')
def summarizealns(alns,output):
    """
    summarize msas
    """
    gap_ncols_p = [];iden_ncols_p = [];total_lengths = []
    for aln in alns:
        total_length,gap_ncol,iden_ncol = basicalninfo(aln)
        gap_ncols_p += [gap_ncol/total_length*100]
        iden_ncols_p += [iden_ncol/total_length*100]
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

@cli.command(context_settings={'help_option_names': ['-h', '--help']})
@click.argument('treefile', type=click.Path(exists=True))
@click.option('--output', '-o', default=None, show_default=True, help='output filename')
@click.option('--outgroup', '-og', default='Zygnema_circumcarinatum_SAG_698-1b', show_default=True, help='outgroup name')
def reroot(treefile,output,outgroup):
    """
    reroot gene tree on algal outgroup
    """
    tree = Phylo.read(treefile,'newick')
    tipnames = [tip.name for tip in tree.get_terminals()]
    outgroupname = [tipname for tipname in tipnames if tipname.startswith(outgroup)][0]
    tree.root_with_outgroup(outgroupname)
    Phylo.draw_ascii(tree)
    outputname = treefile+'.reroot' if output is None else output
    Phylo.write(tree,outputname,format='newick')

@cli.command(context_settings={'help_option_names': ['-h', '--help']})
@click.argument('alns', nargs=-1, type=click.Path(exists=True))
@click.option('--gsmap', '-m', default=None, show_default=True,help='gsmap')
@click.option('--output', '-o', default=None, show_default=True, help='output filename')
def Concat(alns,gsmap,output):
    """
    Concatenate msas
    """
    seqs = {}
    fname = gsmap+'.Concat' if output is None else output
    if gsmap is not None:
        Gsmap = {}
        with open(gsmap,'r') as f:
            for line in f.readlines():
                gene, sp = line.strip().split(' ')
                Gsmap[gene] = sp
    splist = list(set(Gsmap.values()))
    for aln in alns:
        aln_object = AlignIO.read(aln,'fasta')
        for j in range(len(aln_object)):
            spn = aln_object[j].id
            if gsmap is not None: spn = Gsmap[spn]
            sequence = aln_object[j].seq
            if seqs.get(spn) is None:
                seqs[spn] = sequence
            else:
                seqs[spn] = seqs[spn] + sequence
    with open (fname,"w") as f:
        for spn,seq in seqs.items():
            f.write('>{}\n{}\n'.format(spn,seq))

def kde_mode(kde_x, kde_y):
    maxy_iloc = np.argmax(kde_y)
    mode = kde_x[maxy_iloc]
    return mode, max(kde_y)

def get_totalH(Hs):
    CHF = 0
    for i in Hs: CHF = CHF + i
    return CHF

def addvvline(ax,xvalue,color,lstyle,labell):
    if labell == '': ax.axvline(xvalue,color=color, ls=lstyle, lw=1)
    else: ax.axvline(xvalue,color=color, ls=lstyle, lw=1, label='{}: {:.1f}'.format(labell,xvalue))
    return ax

def drawdistribution(X):
    lower,upper = min(X)*0.9,max(X)*1.1
    fig, ax = plt.subplots()
    Hs, Bins, patches = ax.hist(X, bins = np.arange(lower,upper,1), color='gray', alpha=1, rwidth=0.8)
    kde_x = np.linspace(lower,upper,num=500)
    kde_y = stats.gaussian_kde(X,bw_method='scott').pdf(kde_x)
    CHF = get_totalH(Hs)
    scaling = CHF*1
    ax.plot(kde_x, kde_y*scaling, color='k',alpha=0.8, ls = '--',lw = 1,label='KDE curve')
    ax.legend(loc=1,fontsize=10,frameon=False)
    ax.set_xlabel("MLE site rate", fontsize = 10)
    ax.set_ylabel("Number of sites", fontsize = 10)
    inset_ax = fig.add_axes([2/5, 2/5, 8/16, 8/16])
    X_below5 = X[X<=5]
    lower,upper = min(X_below5)*0.9,max(X_below5)*1.1
    Hs, Bins, patches = inset_ax.hist(X_below5, bins = np.arange(lower,upper,0.1), color='gray', alpha=1, rwidth=0.8)
    inset_ax.set_xlabel("MLE site rate", fontsize = 10)
    inset_ax.set_ylabel("Number of sites", fontsize = 10)
    return fig,ax

@cli.command(context_settings={'help_option_names': ['-h', '--help']})
@click.argument('ratefile', type=click.Path(exists=True))
@click.option('--output', '-o', default=None, show_default=True, help='output filename')
def ratedis(ratefile,output):
    """
    draw site rate distribution
    """
    df = pd.read_csv(ratefile,skiprows=6,header=0,index_col=None,sep='\t')
    rates = df['Rate'].to_numpy()
    fig,ax = drawdistribution(rates)
    fig.tight_layout()
    fig.savefig(output,format ='svg', bbox_inches='tight')
    plt.close()

@cli.command(context_settings={'help_option_names': ['-h', '--help']})
@click.argument('treefile', type=click.Path(exists=True))
@click.option('--gsmap', '-m', default=None, show_default=True,help='gsmap')
@click.option('--output', '-o', default=None, show_default=True, help='output filename')
def renamegenetree(treefile,gsmap,output):
    """
    Rename gene names to species names
    """
    Gsmap = {}
    with open(gsmap,'r') as f:
        for line in f.readlines():
            gene, sp = line.strip().split(' ')
            Gsmap[gene] = sp
    tree = Phylo.read(treefile,'newick')
    for tip in tree.get_terminals():
        tip.name = Gsmap[tip.name]
    Phylo.write(tree,output,format='newick')

def plotlinearregression(X,Y,output):
    data_xy = [(x,y) for x,y in zip(X,Y)]
    data_xy = sorted(data_xy, key=lambda x:x[0])
    sorted_X,sorted_Y = np.array([x for x,y in data_xy]),np.array([y for x,y in data_xy])
    fig, ax = plt.subplots(1,1,figsize=(4, 4))
    ax.scatter(sorted_X,sorted_Y,color='k',label='Original data points',alpha=0.5,zorder=1)
    ax.set_xlabel("sCF",fontsize=10)
    ax.set_ylabel("gCF",fontsize=10)
    slope, intercept, r_value, p_value, std_err = stats.linregress(sorted_X,sorted_Y)
    y_fit = intercept + slope * sorted_X
    ax.plot(sorted_X, y_fit, color='k', label='Fitted linear function',zorder=2)
    se = np.sqrt(np.sum((sorted_Y - y_fit)**2) / (len(sorted_Y) - 2))
    t_value = stats.t.ppf(1 - 0.025, df=len(sorted_X) - 2)
    ci = t_value * se * np.sqrt(1/len(sorted_X) + (sorted_X - np.mean(sorted_X))**2 / np.sum((sorted_X - np.mean(sorted_X))**2))
    ax.fill_between(sorted_X, y_fit-ci, y_fit+ci, color='gray', alpha=0.5, label='95% confidence band',zorder=3)
    ax.plot([],[],color='white',label='R-squared: {:.4f}'.format(r_value**2))
    ax.plot([],[],color='white',label='Slope: {:.4f}'.format(slope))
    ax.plot([],[],color='white',label='P-value: {:.4f}'.format(p_value))
    order_list = ['Fitted linear function','95% confidence band','Original data points','R-squared: {:.4f}'.format(r_value**2),'Slope: {:.4f}'.format(slope),'P-value: {:.4f}'.format(p_value)]
    handles, labels = ax.get_legend_handles_labels()
    labels, handles = zip(*sorted(zip(labels, handles), key=lambda x: order_list.index(x[0])))
    ax.legend(handles,labels,fontsize=10,loc=0,frameon=False)
    fig.tight_layout()
    fig.savefig(output,format ='svg', bbox_inches='tight')
    plt.close()

@cli.command(context_settings={'help_option_names': ['-h', '--help']})
@click.argument('scffile', type=click.Path(exists=True))
@click.argument('gcffile', type=click.Path(exists=True))
@click.option('--output', '-o', default=None, show_default=True, help='output filename')
def scfgcfconstrast(scffile,gcffile,output):
    """
    compare scf with gcf
    """
    df_gcf = pd.read_csv(gcffile,skiprows=17,header=0,index_col=None,sep='\t')
    df_scf = pd.read_csv(scffile,skiprows=14,header=0,index_col=None,sep='\t')
    gcf = df_gcf['gCF'].to_list()
    scf = df_scf['sCF'].to_list()
    plotlinearregression(scf,gcf,output)

if __name__ == '__main__':
	cli()
