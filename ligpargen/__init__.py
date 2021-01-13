import logging
logger = logging.getLogger(__name__)

console_handler = logging.StreamHandler()
console_handler.setLevel(logging.ERROR)

logger.addHandler(console_handler)
