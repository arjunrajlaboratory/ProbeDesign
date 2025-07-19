"""Command-line interface for oligo probe design."""

import argparse
import sys
import os
from .core.probe_designer import design_probes


def main():
    """Main command-line interface function."""
    parser = argparse.ArgumentParser(
        description="Design oligo probes for single molecule RNA FISH",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Basic usage - design 48 probes for human sequence
  python -m oligo_probe_design input.fasta
  
  # Design 32 mouse probes with custom parameters
  python -m oligo_probe_design input.fasta --n-oligos 32 --species mouse --target-tm -25
  
  # Disable masking options
  python -m oligo_probe_design input.fasta --no-pseudogene-mask --no-genome-mask
  
  # Save output to specific files
  python -m oligo_probe_design input.fasta --output my_probes
        """
    )
    
    # Required arguments
    parser.add_argument('input_file', 
                       help='Input FASTA file containing target sequence')
    
    # Basic design parameters
    parser.add_argument('--n-oligos', type=int, default=48,
                       help='Number of oligos to design (default: 48)')
    parser.add_argument('--oligo-length', type=int, default=20,
                       help='Length of each oligo (default: 20)')
    parser.add_argument('--spacer-length', type=int, default=2,
                       help='Minimum spacer between oligos (default: 2)')
    
    # Thermodynamic parameters
    parser.add_argument('--target-tm', type=float, default=-23.0,
                       help='Target free energy in kcal/mol (default: -23.0)')
    parser.add_argument('--target-range', nargs=2, type=float, 
                       default=[-26.0, -20.0], metavar=('MIN', 'MAX'),
                       help='Allowable free energy range in kcal/mol (default: -26.0 -20.0)')
    
    # Species and masking
    parser.add_argument('--species', default='human',
                       choices=['human', 'mouse', 'drosophila', 'elegans', 'rat', 'cow'],
                       help='Species for masking databases (default: human)')
    
    # Masking options
    parser.add_argument('--no-pseudogene-mask', action='store_true',
                       help='Disable pseudogene masking')
    parser.add_argument('--no-genome-mask', action='store_true', 
                       help='Disable genome masking')
    parser.add_argument('--gc-run-mask', action='store_true',
                       help='Enable GC run masking (off by default)')
    parser.add_argument('--gc-mask', action='store_true',
                       help='Enable GC composition masking (off by default)')
    parser.add_argument('--no-repeat-mask', action='store_true',
                       help='Disable repeat masking (on by default)')
    
    # Output options
    parser.add_argument('--output', '-o',
                       help='Output file prefix (default: input filename without extension)')
    parser.add_argument('--quiet', '-q', action='store_true',
                       help='Suppress progress messages')
    
    args = parser.parse_args()
    
    # Validate input file
    if not os.path.exists(args.input_file):
        print(f"Error: Input file '{args.input_file}' not found", file=sys.stderr)
        return 1
    
    # Set default output filename if not provided
    if not args.output:
        base_name = os.path.splitext(os.path.basename(args.input_file))[0]
        args.output = base_name
    
    # Validate parameters
    if args.n_oligos <= 0:
        print("Error: Number of oligos must be positive", file=sys.stderr)
        return 1
    
    if args.oligo_length < 10 or args.oligo_length > 50:
        print("Error: Oligo length must be between 10 and 50", file=sys.stderr)
        return 1
    
    if args.spacer_length < 0:
        print("Error: Spacer length must be non-negative", file=sys.stderr)
        return 1
    
    if args.target_range[0] >= args.target_range[1]:
        print("Error: Target range minimum must be less than maximum", file=sys.stderr)
        return 1
    
    # Display parameters
    if not args.quiet:
        print("Oligo Probe Design")
        print("=" * 40)
        print(f"Input file: {args.input_file}")
        print(f"Number of oligos: {args.n_oligos}")
        print(f"Oligo length: {args.oligo_length}")
        print(f"Spacer length: {args.spacer_length}")
        print(f"Target Tm: {args.target_tm} kcal/mol")
        print(f"Target range: {args.target_range[0]} to {args.target_range[1]} kcal/mol")
        print(f"Species: {args.species}")
        print()
        
        print("Masking options:")
        print(f"  Pseudogene masking: {'disabled' if args.no_pseudogene_mask else 'enabled'}")
        print(f"  Genome masking: {'disabled' if args.no_genome_mask else 'enabled'}")
        print(f"  GC run masking: {'enabled' if args.gc_run_mask else 'disabled'}")
        print(f"  GC composition masking: {'enabled' if args.gc_mask else 'disabled'}")
        print(f"  Repeat masking: {'disabled' if args.no_repeat_mask else 'enabled'}")
        print()
    
    try:
        # Run probe design
        results = design_probes(
            input_file=args.input_file,
            n_oligos=args.n_oligos,
            oligo_length=args.oligo_length,
            spacer_length=args.spacer_length,
            target_tm=args.target_tm,
            target_range=args.target_range,
            species=args.species,
            pseudogene_mask=not args.no_pseudogene_mask,
            genome_mask=not args.no_genome_mask,
            gc_run_mask=args.gc_run_mask,
            gc_mask=args.gc_mask,
            repeat_mask=not args.no_repeat_mask,
            output_file=args.output,
            verbose=not args.quiet
        )
        
        if not results['oligos']:
            print("No suitable oligos were found. Try adjusting parameters.", file=sys.stderr)
            return 1
        
        if not args.quiet:
            print()
            print("Design completed successfully!")
            print(f"Designed {len(results['oligos'])} oligos")
            print(f"Output saved with prefix: {args.output}")
            
    except KeyboardInterrupt:
        print("\nOperation cancelled by user", file=sys.stderr)
        return 1
    except Exception as e:
        print(f"Error during probe design: {e}", file=sys.stderr)
        if not args.quiet:
            import traceback
            traceback.print_exc()
        return 1
    
    return 0


if __name__ == '__main__':
    sys.exit(main())