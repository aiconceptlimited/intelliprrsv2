"""Centralized logging for publication pipeline."""

import logging
import sys
from pathlib import Path
from datetime import datetime


def setup_logging(
    log_file: Path = None,
    level: str = 'INFO',
    console: bool = True
) -> logging.Logger:
    """Set up logging with console and file output."""
    
    logger = logging.getLogger('intelliprrsv2_paper')
    logger.setLevel(getattr(logging, level.upper()))
    
    # Remove existing handlers
    logger.handlers.clear()
    
    # Console handler
    if console:
        console_handler = logging.StreamHandler(sys.stdout)
        console_handler.setLevel(getattr(logging, level.upper()))
        console_format = '%(asctime)s | %(levelname)-8s | %(message)s'
        console_handler.setFormatter(logging.Formatter(console_format))
        logger.addHandler(console_handler)
    
    # File handler
    if log_file:
        log_file.parent.mkdir(parents=True, exist_ok=True)
        file_handler = logging.FileHandler(log_file)
        file_handler.setLevel(logging.DEBUG)
        file_format = '%(asctime)s | %(levelname)-8s | %(name)s | %(message)s'
        file_handler.setFormatter(logging.Formatter(file_format))
        logger.addHandler(file_handler)
    
    return logger


def get_logger(name: str = None) -> logging.Logger:
    """Get the configured logger."""
    return logging.getLogger('intelliprrsv2_paper')
