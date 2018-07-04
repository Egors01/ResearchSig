import time
import os
from tools.set_env import Environment

class Logger:
    def __init__(self, source_name='default run log', filename='default_log'):
        self.log_path = SetUpEnv().log_path
        if not os.path.exists(self.log_path):
            os.makedirs(self.log_path)
        self.log_path = os.path.join(self.log_path, filename+'.log')
        self.time_start = time.time()
        self.source_name = source_name
        self.time_tick = time.time()
        self.space = ' ' * len(self.source_name)
        f_log = open(self.log_path, 'w')
        f_log.write(self.source_name + ' log: Started at ' + str(time.asctime(time.localtime())) + '\n')
        f_log.close()
        pass

    def tick(self):
        self.tick = time.time()

    def print_message(self, message, print_time=False):
        f_log = open(self.log_path, 'a')
        f_log.write(self.source_name + ' log: ' + message + ' \n')
        if print_time:
            diff = time.time() - self.time_tick
            min = round(diff // 60, 2)
            sec = round(diff % 60, 2)
            f_log.write(self.space + ' passed time: ' + str(min) + ' min ' + str(sec) + ' s\n')
        f_log.close()

    def print_end(self):
        f_log = open(self.log_path, 'a')
        f_log.write(self.source_name + ' log: process ended at ' + str(time.asctime(time.localtime())) + '\n')
        diff = time.time() - self.time_start
        min = round(diff // 60, 2)
        sec = round(diff % 60,2)

        f_log.write(self.space + ' time passed from start ' + str(min) + ' min ' + str(sec) + ' s\n')
        f_log.close()



