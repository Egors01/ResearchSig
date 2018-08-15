import time
import os

from constants import LOG_PATH
from constants import RUN_ID
from tools.set_env import Environment


class Logger:
    def __init__(self, source_name='main_thread', filename='run_log', msg='', reset=False):
        self.log_path = LOG_PATH
        if not os.path.exists(self.log_path):
            os.makedirs(self.log_path)
        filename+='_'+str(RUN_ID)
        source_name+=' run_id '+str(RUN_ID)
        self.log_path = os.path.join(self.log_path, filename + '.log')
        self.time_start = time.time()
        self.logger_source_name = source_name
        self.space = ' ' * len(self.logger_source_name)
        self.timers = Timer()
        if not os.path.isfile(self.log_path) or reset:
            f_log = open(self.log_path, 'w')
            f_log.write(self.logger_source_name + ' log: Started at ' + str(time.asctime(time.localtime())) + '\n')

        elif msg != '':
            f_log = open(self.log_path, 'a')
            f_log.write(self.logger_source_name + ' log: ' + msg + '\n')
            f_log.close()
        else:
            pass

    def tick(self, timer_name=0):
        self.timers.make_new_tick(timer_name)

    def message(self, msg, source_name='', print_time=False, timer_name=0):
        f_log = open(self.log_path, 'a')
        if source_name == '':
            source_name = self.logger_source_name
        f_log.write(source_name + ' log: ' + msg + ' \n')
        if print_time:
            diff = time.time() - self.timers.get_timer_by_name(timer_name)
            min = round(diff // 60, 2)
            sec = round(diff % 60, 2)
            f_log.write(source_name + ' log: passed time: ' + str(min) + ' min ' + str(sec) + ' s\n')
        f_log.close()

    def print_end(self):
        f_log = open(self.log_path, 'a')
        f_log.write(self.logger_source_name + ' log: all parts ended at ' + str(time.asctime(time.localtime())) + '\n')
        diff = time.time() - self.time_start
        min = round(diff // 60, 2)
        sec = round(diff % 60, 2)

        f_log.write(self.logger_source_name+ ' log: time passed from start ' + str(min) + ' min ' + str(sec) + ' s\n')
        f_log.close()


class Timer():
    def __init__(self):
        self.ticks = {0: time.time()}

    def make_new_tick(self, tick_name):
        self.ticks[tick_name] = time.time()

    def get_timer_by_name(self, tick_name):
        return self.ticks[tick_name]
