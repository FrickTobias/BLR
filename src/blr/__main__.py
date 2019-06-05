import sys
import logging


logger = logging.getLogger(__name__)


def main() -> int:
    logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
    logger.warning("A warning!")
    logger.info("Hello")

    return 0


if __name__ == "__main__":
    sys.exit(main())
