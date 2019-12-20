import logging

def setup_custom_logger(name, pack_output):
    formatter = logging.Formatter(fmt='%(asctime)s - %(levelname)s - %(module)s - %(message)s')

    handler = logging.StreamHandler()
    handler.setFormatter(formatter)

    logging.basicConfig(filename=pack_output+'qc_data_collector.log', level=logging.INFO)
    logger = logging.getLogger(name)
    logger.addHandler(handler)

    return logger