"""Command-line interface for ProbeDesign."""

import click
from pathlib import Path

from .core import design_probes
from .output import write_output_files, format_probes_table


@click.group()
@click.version_option()
def main():
    """ProbeDesign - Design oligonucleotide probes for FISH experiments."""
    pass


@main.command()
@click.argument('input_file', type=click.Path(exists=True))
@click.option(
    '-n', '--probes', 'n_probes',
    default=48,
    type=int,
    help='Number of probes to design (default: 48)'
)
@click.option(
    '-l', '--oligo-length',
    default=20,
    type=int,
    help='Length of each oligonucleotide (default: 20)'
)
@click.option(
    '-s', '--spacer-length',
    default=2,
    type=int,
    help='Minimum gap between probes (default: 2)'
)
@click.option(
    '-g', '--target-gibbs',
    default=-23.0,
    type=float,
    help='Target Gibbs free energy in kcal/mol (default: -23)'
)
@click.option(
    '--allowable-gibbs',
    default='-26,-20',
    type=str,
    help='Allowable Gibbs FE range as min,max (default: -26,-20)'
)
@click.option(
    '-o', '--output',
    default=None,
    type=str,
    help='Output file prefix (default: derived from input filename)'
)
@click.option(
    '--quiet', '-q',
    is_flag=True,
    help='Suppress output to stdout'
)
@click.option(
    '--species',
    default='human',
    type=click.Choice(['human', 'mouse', 'elegans', 'drosophila', 'rat']),
    help='Species for masking databases (default: human)'
)
@click.option(
    '--pseudogene-mask/--no-pseudogene-mask',
    default=False,
    help='Mask regions that align to pseudogenes (default: off)'
)
@click.option(
    '--genome-mask/--no-genome-mask',
    default=False,
    help='Mask repetitive genomic regions (default: off)'
)
@click.option(
    '--index-dir',
    default=None,
    type=click.Path(exists=True, file_okay=False),
    help='Directory containing bowtie indexes (default: auto-detect)'
)
@click.option(
    '--repeatmask-file',
    default=None,
    type=click.Path(exists=True),
    help='FASTA file with N\'s marking repeat regions (for manual repeat masking)'
)
@click.option(
    '--repeatmask/--no-repeatmask',
    default=False,
    help='Run RepeatMasker to automatically mask repeat regions (default: off)'
)
def design(
    input_file: str,
    n_probes: int,
    oligo_length: int,
    spacer_length: int,
    target_gibbs: float,
    allowable_gibbs: str,
    output: str,
    quiet: bool,
    species: str,
    pseudogene_mask: bool,
    genome_mask: bool,
    index_dir: str,
    repeatmask_file: str,
    repeatmask: bool,
):
    """Design probes for a target sequence.

    INPUT_FILE is a FASTA file containing the target sequence.

    Examples:

        probedesign design input.fa

        probedesign design input.fa -n 6 -l 20 -o MyProbes

        probedesign design input.fa --target-gibbs -23 --allowable-gibbs -26,-20

        probedesign design input.fa --pseudogene-mask --species human
    """
    # Parse allowable Gibbs range
    try:
        min_gibbs, max_gibbs = map(float, allowable_gibbs.split(','))
    except ValueError:
        raise click.BadParameter(
            f"Invalid format for --allowable-gibbs: {allowable_gibbs}. "
            "Expected format: min,max (e.g., -26,-20)"
        )

    # Determine output prefix
    if output is None:
        output = Path(input_file).stem

    # Handle --repeatmask flag: run RepeatMasker and use output as repeatmask_file
    actual_repeatmask_file = repeatmask_file
    if repeatmask:
        if repeatmask_file:
            raise click.BadParameter(
                "Cannot use both --repeatmask and --repeatmask-file. "
                "Choose one or the other."
            )
        if not quiet:
            click.echo(f"Running RepeatMasker on {input_file}...")
        try:
            from .masking import run_repeatmasker
            masked_file = run_repeatmasker(Path(input_file), species=species)
            actual_repeatmask_file = str(masked_file)
            if not quiet:
                click.echo(f"RepeatMasker output: {masked_file}")
        except FileNotFoundError as e:
            raise click.ClickException(str(e))
        except RuntimeError as e:
            raise click.ClickException(str(e))

    # Run probe design
    if not quiet:
        click.echo(f"Designing probes for {input_file}...")

    result = design_probes(
        input_file=input_file,
        n_probes=n_probes,
        oligo_length=oligo_length,
        spacer_length=spacer_length,
        target_gibbs=target_gibbs,
        allowable_gibbs=(min_gibbs, max_gibbs),
        output_name=output,
        species=species,
        pseudogene_mask=pseudogene_mask,
        genome_mask=genome_mask,
        index_dir=index_dir,
        repeatmask_file=actual_repeatmask_file,
    )

    if not result.probes:
        click.echo("ERROR: Could not find any valid probes.", err=True)
        raise SystemExit(1)

    # Write output files
    write_output_files(result, output, mask_seqs=result.mask_strings)

    # Print summary
    if not quiet:
        click.echo(f"\nFound {len(result.probes)} probes (score: {result.score:.4f})")
        click.echo(f"\nOutput files:")
        click.echo(f"  {output}_oligos.txt")
        click.echo(f"  {output}_seq.txt")
        click.echo(f"\nProbes:")
        click.echo(format_probes_table(result))


@main.command()
@click.argument('input_file', type=click.Path(exists=True))
@click.option(
    '-l', '--oligo-length',
    default=20,
    type=int,
    help='Length of oligonucleotide to analyze (default: 20)'
)
def analyze(input_file: str, oligo_length: int):
    """Analyze a sequence for probe design feasibility.

    Shows the distribution of Gibbs free energies across the sequence.
    """
    from .fasta import read_fasta, sequences_to_single_string
    from .sequence import clean_sequence, has_invalid_chars
    from .thermodynamics import gibbs_rna_dna

    headers, seqs = read_fasta(input_file)
    seq = clean_sequence(sequences_to_single_string(seqs))

    click.echo(f"Sequence length: {len(seq)}")
    click.echo(f"Analyzing with oligo length: {oligo_length}")

    gibbs_values = []
    invalid_count = 0

    for i in range(len(seq) - oligo_length + 1):
        oligo = seq[i:i + oligo_length]
        if has_invalid_chars(oligo):
            invalid_count += 1
        else:
            try:
                gibbs_values.append(gibbs_rna_dna(oligo))
            except KeyError:
                invalid_count += 1

    if gibbs_values:
        click.echo(f"\nGibbs FE distribution:")
        click.echo(f"  Min: {min(gibbs_values):.2f} kcal/mol")
        click.echo(f"  Max: {max(gibbs_values):.2f} kcal/mol")
        click.echo(f"  Mean: {sum(gibbs_values)/len(gibbs_values):.2f} kcal/mol")

    click.echo(f"\nValid positions: {len(gibbs_values)}")
    click.echo(f"Invalid positions: {invalid_count}")


if __name__ == '__main__':
    main()
