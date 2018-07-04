from constants import OS_TYPE,PRIMARY_PATH
import os


class Environment(object):
    def __init__(self,):
        if OS_TYPE=='LINUX':
            self.primary_path = PRIMARY_PATH
            self.variant_path = os.path.join(self.primary_path,'data','variants')
            self.raw_maf_path = os.path.join(self.primary_path,'data','maf_data')
            self.fasta_chrs_path = os.path.join(self.primary_path,'data', 'chr_raw_fasta_data')
            self.log_path = os.path.join(self.primary_path,'logs')
            self.resourses_path =os.path.join(self.primary_path,'resourses')
            attributes= [attrib for attrib in dir(self) if '__' not in attrib]
            for path in self.__dict__.values():
                if not os.path.exists(path):
                    os.makedirs(path)
            self.ref_startpospath = os.path.join(self.primary_path,'resourses','start_positions.csv')
            self.species_names_path = os.path.join(self.primary_path,'resourses','species_names.csv')
    def get_primary_path(self):
        return self.primary_path

    def get_ref_startpospath(self):
        return self.ref_startpospath

Environment()
def import_create_print():
    df={'A':[1,2,3],'B':[3,4,5]}
    print(df)
