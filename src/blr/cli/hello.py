"""
Hello
"""
import logging

logger = logging.getLogger(__name__)


def add_arguments(parser):
    arg = parser.add_argument
    arg('--value', default="world")


def main(args):
    logger.info('Hello, %s', args.value)

