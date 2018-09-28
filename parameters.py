from tools.input_data import InputParameters




input_data = InputParameters()
input_data.add_pair(first="hg38", second="panTro4", ref="inner.18", run_id=999)


input_data.add_pair(first="panPan1", second=   "hg38", ref="ponAbe2", run_id=0)
input_data.add_pair(first="gorGor3", second=   "hg38", ref="ponAbe2", run_id=0)
input_data.add_pair(first="ponAbe2", second=   "hg38", ref="nomLeu3", run_id=0)
input_data.add_pair(first="panTro4", second=   "hg38", ref="ponAbe2", run_id=0)
input_data.add_pair(first="nomLeu3", second=   "hg38", ref="ponAbe2", run_id=0)
input_data.add_pair(first="panTro4", second="panPan1", ref="ponAbe2", run_id=0)
input_data.add_pair(first="panTro4", second="gorGor3", ref="ponAbe2", run_id=0)
input_data.add_pair(first="gorGor3", second="panPan1", ref="ponAbe2", run_id=0)
input_data.add_pair(first="gorGor3", second="ponAbe2", ref="nomLeu3", run_id=0)
input_data.add_pair(first="nomLeu3", second="ponAbe2", ref="rheMac3", run_id=0)
input_data.add_pair(first="nomLeu3", second="ponAbe2", ref="papAnu2", run_id=0)
input_data.add_pair(first="nomLeu3", second="ponAbe2", ref="nasLar1", run_id=0)



#input_data.add_pair(first="nomLeu3", second="hg38", ref="rheMac3", run_id=5)
#input_data.add_pair(first="nomLeu3", second="panTro4", ref="rheMac3", run_id=5)


input_data.add_pair(first="nasLar1", second="rhiRox1", ref="hg38", run_id=6)
input_data.add_pair(first="nasLar1", second="rhiRox1", ref="rheMac3", run_id=6)
input_data.add_pair(first="nasLar1", second="rhiRox1", ref="chlSab2", run_id=6)


input_data.add_pair(first="macFas5", second="rheMac3", ref="papAnu2", run_id=8)
input_data.add_pair(first="rheMac3", second="papAnu2", ref="chlSab2", run_id=8)
input_data.add_pair(first="papAnu2", second="chlSab2", ref="rhiRox1", run_id=8)

input_data.add_pair(first="papAnu2", second="chlSab2", ref="nasLar1", run_id=8)
input_data.add_pair(first="chlSab2", second="rhiRox1", ref="nomLeu3", run_id=8)
input_data.add_pair(first="chlSab2", second="nasLar1", ref="nomLeu3", run_id=8)


input_data.add_pair(first="chlSab2", second="rhiRox1", ref="hg38", run_id=8)
input_data.add_pair(first="chlSab2", second="nasLar1", ref="hg38", run_id=8)
input_data.add_pair(first="chlSab2", second="rhiRox1", ref="nomLeu3", run_id=8)


input_data.add_pair(first="calJac3", second="saiBol1", ref="nomLeu3", run_id=8)
input_data.add_pair(first="calJac3", second="saiBol1", ref="nasLar1", run_id=8)
input_data.add_pair(first="calJac3", second="saiBol1", ref="hg38", run_id=8)

#by 18 sept nomLeu outgroup
input_data.add_pair(first="hg38", second="ponAbe2", ref="nomLeu3", run_id=1990)
input_data.add_pair(first="panTro4", second="ponAbe2", ref="nomLeu3", run_id=1990)
input_data.add_pair(first="panPan1", second="ponAbe2", ref="nomLeu3", run_id=1990)
input_data.add_pair(first="gorGor3", second="ponAbe2", ref="nomLeu3", run_id=1990)

input_data.add_pair(first="hg38", second="gorGor3", ref="nomLeu3", run_id=1991)
input_data.add_pair(first="panTro4", second="gorGor3", ref="nomLeu3", run_id=1991)
input_data.add_pair(first="panPan1", second="gorGor3", ref="nomLeu3", run_id=1991)

input_data.add_pair(first="hg38", second="panPan1", ref="nomLeu3", run_id=1992)
input_data.add_pair(first="panTro4", second="panPan1", ref="nomLeu3", run_id=1992)

input_data.add_pair(first="hg38", second="panTro4", ref="nomLeu3", run_id=1992)




input_data.add_pair(first="panTro4", second="hg38", ref="nomLeu3", run_id=1993)
#input_data.add_pair(first="ponAbe2", second="hg38", ref="nomLeu3", run_id=1990)
input_data.add_pair(first="panPan1", second="hg38", ref="nomLeu3", run_id=1993)
input_data.add_pair(first="gorGor3", second="hg38", ref="nomLeu3", run_id=1993)


input_data.add_pair(first="hg38", second="panTro4", ref="inner.18", run_id=28090)
input_data.add_pair(first="panPan1", second="panTro4", ref="inner.17", run_id=28090)
input_data.add_pair(first="panPan1", second="hg38", ref="inner.18", run_id=28092)
input_data.add_pair(first="hg38", second="gorGor3", ref="inner.19", run_id=28091)
input_data.add_pair(first="hg38", second="ponAbe2", ref="inner.20", run_id=28091)
# PAIRS = [["ponAbe2", "panPan1", "panTro4"],
#          ["ponAbe2", "panPan1", "hg38"],
#          ["ponAbe2", "gorGor3", "panPan1"],
#          ["nomLeu3", "ponAbe2", "gorGor3"],
#          ["rheMac3", "nomLeu3", "ponAbe2"],
#
#          ["ponAbe2", "gorGor3", "hg38"],
#
#          ["nomLeu3", "ponAbe2", "hg38"],
#          ["rheMac3", "nomLeu3", "hg38"]]
