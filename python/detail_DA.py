import numpy as np
import sys

import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties

def get_time( fn ):

   skip = True
#   skip = False

   TIME = { "NCYCLE":0 }
   dcl_ocnt = 0
   dcl_acnt = 0

   f = open(fn, encoding='UTF-8')
   lines = open(fn, encoding='UTF-8').readlines()
   for l in range(len(lines)):
      data = lines[l].split(" ,")

      if data[0] == "FINALIZE" or data[0] == "INITIALIZE" or data[0] == "INIT_LETKF" or \
         data[0] == "SET_GRID" or data[0] == "DEALLOCATE" or data[0] == "INITIALIZE_OTHERS":
         continue 

      if data[0] == "DCL_OBS" and dcl_ocnt == 0:
         dcl_ocnt = 1
         continue
      if data[0] == "DCL_ANAL" and dcl_acnt == 0:
         dcl_acnt = 1
         continue

      if not skip:
#         if data[0] == "FINALIZE": 
#            break

         # count # of cycles
         if data[0] == "SCALE":
            TIME["NCYCLE"] += 1

         if not data[0] in TIME.keys():
            TIME[data[0]] = float(data[1])
         else:
            TIME[data[0]] += float(data[1])


      # skip 1st cycle
      if data[0] == "SEND":
         skip = False

   return(TIME)

def draw_graph( TIME ):

   fig, ((ax)) = plt.subplots(1, 1, figsize=(4.5, 4.5))
   fig.subplots_adjust(left=0.2, bottom=0.1, right=0.5, top=0.9, wspace=0.2, hspace=0.3)

#   print(TIME.keys())
#   print(len(TIME.keys()))
#   sys.exit()

   left = 0.5
#   ax.bar(left, TIME["SCALE"]/TIME["NCYCLE"])

   k_list = list(TIME.keys())
   print(k_list)

   xlabels = ["Nov2019"]

#   cc_l = ["k"]*len(k_list)
#   cc_l[1] = "deepskyblue"
#   cc_l[2] = "limegreen"
#   cc_l[3] = "r"
#   cc_l[4] = "y"
#   cc_l[5] = "gray"
   bottom = 0.0
   for k in range(1, len(k_list)):

      if k_list[k] == "FINALIZE":
         continue

#   ax.bar(left, TIME[k_list[1]]/TIME["NCYCLE"])
      val = TIME[k_list[k]]/TIME["NCYCLE"]
      if k_list[k] == "READ_OBS":
         val -= TIME["DCL_OBS"]/TIME["NCYCLE"]
      if k_list[k] == "WRITE_ANAL":
         val -= TIME["DCL_ANAL"]/TIME["NCYCLE"]

      ax.bar(left, val, bottom=bottom, label=k_list[k]) #, color=cc_l[k])

      # bottom for the next
      bottom += val

   handles, labels = ax.get_legend_handles_labels() 
   ax.legend(handles[::-1], labels[::-1],
             bbox_to_anchor=(1.05, 1), loc='upper left')

#   ax.set_xticks(left, xlabels)
   labels = [item.get_text() for item in ax.get_xticklabels()]
   labels[1] = 'Nov2019'
   ax.set_xticklabels(labels)

   ax.set_xlim(0,1)
   ax.set_ylim(0,15)
   ax.set_ylabel("Time (s)")

   plt.show()
   sys.exit()

   for key in TIME.keys():
  
      if key == "NCYCLE":
         continue
      print(key)


####
fn = "/data15/honda/SCALE-LETKF/AIP_SAFE/HPC20191128/D4_250M_mercator_6528_allsngl_ldt2.0/exp/2781300_cycle_20190610080000/DA_TIMER.txt"
#fn = "/data15/honda/SCALE-LETKF/AIP_SAFE/HPC20190829/D4_250M_NP4096_600S_SMTH_HPC0829_6528/exp/2656030_cycle_20190610080000/DA_TIMER.txt"


TIME = get_time( fn )
print (TIME)
draw_graph( TIME )
