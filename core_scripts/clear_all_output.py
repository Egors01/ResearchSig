from constants import *
import os, shutil

varinats_folders_list = [VCF_VARIANT_PATH,
                         FILT_VCF_VARIANT_PATH,
                         LOG_PATH,
                         MATRICES_PATH]
for folder in varinats_folders_list:
    for the_file in os.listdir(folder):
        file_path = os.path.join(folder, the_file)
        try:
            if os.path.isfile(file_path):
                os.unlink(file_path)
        except Exception as e:
            print(e)