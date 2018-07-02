from tools.logger import Logger
from tools.set_env import SetUpEnv
if __name__ == '__main__':
    env = SetUpEnv()
    var_log = Logger('variant_call','vclog')
    for i in range(1, 20000000):
        a = 8 * i
    var_log.tick()
    var_log.print_message('in the middle ',print_time=True)
    for i in range(1, 20000000):
        a = 8 * i
    var_log.print_message('at the end ', print_time=True)
    var_log.print_end()