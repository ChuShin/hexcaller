import click
from hexcaller.src.call_cnv import call_cnv_from_cov

@click.command()
@click.argument('he_gene_pair', type=click.Path(exists=True))
@click.argument('cov_file', type=click.Path(exists=True))
def call_cnv(he_gene_pair,cov_file):
    """call CNV based on read coverage."""
    call_cnv_from_cov(he_gene_pair,cov_file)
