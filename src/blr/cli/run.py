"""
Run the BLR pipeline

This is a small wrapper around Snakemake that sets some default parameters
"""
import sys
import logging
from importlib_resources import path as resource_path
from snakemake import snakemake

logger = logging.getLogger(__name__)


def add_arguments(parser):
    arg = parser.add_argument
    arg('--dryrun', '-n', default=False, action='store_true',
        help='Do not execute anything')
    # TODO
    # Change the default for --cores (use available_cpu_count() perhaps)
    arg('--cores', '--jobs', '-j', metavar='N', type=int, default=4,
        help='Run on at most N CPU cores in parallel. '
        'Default: %(default)s')
    arg('--keepgoing', '-k', default=False, action='store_true',
        help='If one job fails, finish the others.')
    arg("--dag", default=False, action='store_true',
        help="Print the dag in the graphviz dot language (requires graphviz to be installed). Default: %(default)s. "
             "To get output to pdf file, pipe output into dot as follows: blr run --dag | dot -Tpdf > dag.pdf")
    arg('targets', nargs='*', default=[],
        help='File(s) to create. If omitted, the full pipeline is run.')


def main(args):
    # snakemake sets up its own logging, and this cannot be easily changed
    # (setting keep_logger=True crashes), so remove our own log handler
    # for now
    logger.root.handlers = []
    with resource_path('blr', 'Snakefile') as snakefile_path:
        success = snakemake(
            snakefile_path,
            snakemakepath='snakemake',
            dryrun=args.dryrun,
            cores=args.cores,
            keepgoing=args.keepgoing,
            printshellcmds=True,
            printdag=args.dag,
            targets=args.targets if args.targets else None,
        )
    sys.exit(0 if success else 1)
