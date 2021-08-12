from .eri_cal import generate_h2o_data

def call(cfg, logger):
    if cfg.task == 'collect_data':
        generate_h2o_data(cfg, logger)