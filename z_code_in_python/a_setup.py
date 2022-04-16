import platform
import os
import pandas as pd
import numpy as np

# set wd by os
win_path = "D:\\MyZone\\Dropbox\\pool_80_for_paper\\16s_dairyCattle\\code_in_python"
mac_path = "/Users/lzhao13/Dropbox/pool_80_for_paper/16s_dairyCattle/"
if platform.system() == "Windows":
    os.chdir(win_path)
else:
    os.chdir(mac_path)

# costant
N = 90
MIN_DEPTH = 12637

# variable
col_index_list_name = ["HC","PWD","BH","SCU","FH","MLLH","FC","MLC","PLLC","FOC","CUC","Hos"]
col_index_list = [
    list(range(1, 10)), list(range(10, 16)),    # HC, PWD
    list(range(16, 24)), list(range(24, 32)),    # BH, SCU
    list(range(32, 40)), list(range(40, 48)),    # FH, MLLH
    list(range(48, 56)), list(range(56, 64)),    # FC, MLC
    list(range(64, 72)), list(range(72, 75)),    # PLLC, FOC
    list(range(75, 83)), list(range(83, 91))]    # CUC, Hos

mu = ["01-HC", "02-PWH", "03-BH", "04-SCU", "05-FH", "06-MLLH", 
      "07-FC", "08-MLC", "09-PLLC", "10-FOC", "11-CUC", "12-Hos"]
sample_size = [len(x) for x in col_index_list]
# sample_size = [9,6,8,8,8,8,8,8,8,3,8,8]
stage_arr = np.repeat(mu, sample_size)
stage = stage_arr.tolist()

production_tracks = ("HC.001", "HC.002", "HC.003", "HC.004", "HC.005", "HC.006", "HC.007", "HC.008", "HC.009",
                     "PWH.001", "PWH.002", "PWH.003", "PWH.004", "PWH.005", "PWH.006",
                     "BH.001", "BH.002", "BH.003", "BH.004", "BH.005", "BH.006", "BH.007", "BH.008",
                     "SCU.001", "SCU.002", "SCU.003", "SCU.004", "SCU.005", "SCU.006", "SCU.007", "SCU.008",
                     "FH.001", "FH.002", "FH.003", "FH.004", "FH.005", "FH.006", "FH.007", "FH.008",
                     "MLLH.001","MLLH.002", "MLLH.003", "MLLH.004", "MLLH.005", "MLLH.006", "MLLH.007", "MLLH.008",
                     "FC.001", "FC.002","FC.003", "FC.004", "FC.005", "FC.006", "FC.007", "FC.008",
                     "MLC.001", "MLC.002", "MLC.003", "MLC.004", "MLC.005", "MLC.006", "MLC.007", "MLC.008",
                     "PLLC.001","PLLC.002","PLLC.003", "PLLC.004","PLLC.005","PLLC.006","PLLC.007","PLLC.008",
                     "FOC.001", "FOC.002", "FOC.003",
                     "CUC.001", "CUC.002", "CUC.003", "CUC.004", "CUC.005","CUC.006", "CUC.007", "CUC.008",
                     "Hos.001", "Hos.002",  "Hos.003",  "Hos.004",  "Hos.005",  "Hos.006", "Hos.007", "Hos.008",
                     "Kingdom", "Phylum","Class","Order","Family","Genus")    # correct order
# palette
my_colros = ("steelblue2","blue3","seagreen3","green4","dimgray","coral4",
             "pink2","gold","darkorange2","deeppink","purple","red1")  
