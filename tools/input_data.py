from constants import SPECIES_NAMES

class InputParameters(object):
    def __init__(self):
        self.run_ids_and_pairs = {}
    def get_pairs_by_run_id(self, run_id):
        run_id = str(run_id)
        return self.run_ids_and_pairs.get(run_id)

    def add_pair(self, first, second, ref, run_id):
        run_id=str(run_id)
        for sp_name in [first, second, ref]:
            if sp_name in SPECIES_NAMES:
                pass
            else:
                Exception('Invalid species name in pair {}\n run id {}'.format(
                    sp_name, run_id))
        self.run_ids_and_pairs.setdefault(run_id, []).append(
            Pair(first=first, second=second, ref=ref))

class Pair():
    def __init__(self, first, second, ref):
        self.species_1 = first
        self.species_2 = second
        self.reference = ref
    def string(self):
        return self.species_1+' '+self.species_2+' ref '+self.reference
