from constants import OS_TYPE
import os
import inspect


class SetUpEnv(object):
    def __init__(self,):
        if OS_TYPE=='LINUX':
            self.primary_path = os.getcwd()
            self.variant_path = os.path.join(self.primary_path,'variants')
            self.chrs_path = os.path.join(self.primary_path,'chr_raw_fasta_data')
            self.info_path = os.path.join(self.primary_path,'supplementary_data')
            self.log_path = os.path.join(self.primary_path,'logs')
            attributes= [attrib for attrib in dir(self) if '__' not in attrib]
            for path in self.__dict__.values():
                if not os.path.exists(path):
                    os.makedirs(path)

    def get_primary_path(self):
        return self.primary_path



def import_create_print():
    df={'A':[1,2,3],'B':[3,4,5]}
    print(df)
