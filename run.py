from eri import call
import hydra
from hydra.experimental import compose, initialize
from omegaconf import DictConfig, OmegaConf
import os
import logging
import sys

class Logger:
    def __init__(self, filename="Default.log"):
        self.terminal = sys.stdout
        self.log = open(filename, "a")
 
    def write(self, message):
        self.terminal.write(message)
        self.log.write(message)
 
    def flush(self):
        pass

def setup_logging():
    """Enable pretty logging and sets the level to DEBUG."""
    logging.addLevelName(logging.DEBUG, "D")
    logging.addLevelName(logging.INFO, "I")
    logging.addLevelName(logging.WARNING, "W")
    logging.addLevelName(logging.ERROR, "E")
    logging.addLevelName(logging.CRITICAL, "C")

    formatter = logging.Formatter(
        fmt=("%(levelname)s%(asctime)s" " [%(module)s:%(lineno)d] %(message)s"),
        datefmt="%m%d %H:%M:%S",
    )

    # create console handler and set level to info
    console_handler = logging.StreamHandler()
    console_handler.setLevel(logging.INFO)
    console_handler.setFormatter(formatter)

    # create logger and set level to debug
    logger = logging.getLogger()
    logger.handlers = []
    logger.setLevel(logging.DEBUG)
    logger.propagate = False
    logger.addHandler(console_handler)

    return logger


def init(cfg):
    logger = setup_logging()
    print_log(cfg, logger)
    return logger

def exit(cfg):
    pass

def print_log(cfg, logger):
    sys.stdout = Logger('stdout.txt')
    logger.info(OmegaConf.to_yaml(cfg))

@hydra.main(config_path="conf", config_name="config")
def run(cfg:DictConfig):
    logger = init(cfg)
    call(cfg, logger)
    exit(cfg)

def run_debug():
    # global initialization
    with initialize(config_path="conf", job_name="test_app"):
        cfg = compose(config_name="config")
        logger = init(cfg)
        call(cfg, logger)
        exit(cfg)

if __name__ == "__main__":
    # run_debug()
    run()
   
