from core_scripts.variantcaller import variant_call_pairs


PAIRS = [["nomLeu3", "ponAbe2", "hg38"]]
chr_names = ['example', 'example2', 'example3']
variant_call_pairs(PAIRS,chr_names)

    # env = Environment()
    # var_log = Logger('variant_call','vclog')
    # for i in range(1, 20000000):
    #     a = 8 * i
    # var_log.tick()
    # var_log.message('in the middle ', print_time=True)
    # for i in range(1, 20000000):
    #     a = 8 * i
    # var_log.message('at the end ', print_time=True)
    # var_log.print_end()

