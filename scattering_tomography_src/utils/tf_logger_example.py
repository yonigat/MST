from utils.tf_logger import Logger
import numpy as np

logger = Logger('/home/adamgeva/tmp')
for i in range(1000):
    logger.log_histogram('test_hist', np.random.rand(50) * (i + 1), i)
