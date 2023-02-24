import logging
import os
from logging import Logger


def get_module_logger(mod_name: str) -> Logger:
    log_level = os.getenv("LOG_LEVEL", "WARNING")
    logger = logging.getLogger(mod_name)
    handler = logging.StreamHandler()
    formatter = logging.Formatter("%(asctime)s %(name)-12s %(levelname)-8s %(message)s")
    handler.setFormatter(formatter)
    logger.addHandler(handler)
    logger.info(f"log_level provided from env vars is {log_level}")
    if log_level == "DEBUG":
        logger.setLevel(logging.DEBUG)
    elif log_level == "INFO":
        logger.setLevel(logging.INFO)
    elif log_level == "WARNING":
        logger.setLevel(logging.WARNING)
    else:
        logger.setLevel(logging.WARNING)
    return logger
