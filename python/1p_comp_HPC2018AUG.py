import numpy as np
import sys
import os

top = "/data_honda01/honda/SCALE-LETKF/scale-5.3.3/OUTPUT/AIP_REALTIME"
top = "/data15/honda/SCALE-LETKF/AIP_SAFE/HPC20190829"
exp6528 = "D4_250M_NP4096_600S_SMTH_HPC0829_6528"
exp4032 = "D4_250M_NP4096_600S_SMTH_HPC0829_4032"
exp8192 = "D4_250M_NP4096_600S_SMTH_HPC0829_8192"

fn_DA6528 = os.path.join(top,exp6528,"exp/2656030_cycle_20190610080000/DA.txt")
fn_FCST6528 = os.path.join(top,exp6528,"exp/2656030_cycle_20190610080000/FCST.txt")

fn_DA8192 = os.path.join(top,exp8192,"exp/2656036_cycle_20190610080000/DA.txt")

fn_DA4032 = os.path.join(top,exp4032,"exp/2655604_cycle_20190610080000/DA.txt")
fn_FCST4032 = os.path.join(top,exp4032,"exp/2655604_cycle_20190610080000/FCST.txt")

quick = True
#quick = False

def read_txt(fn):

   try:
      f = open(fn)
   except:
      print("Failed to open")
      print(fn)
      sys.exit()
   
   lines = f.read().split('\n')

   data = []
   for line in lines:
      if line == "":
         continue
      data.append(float(line))

   return(np.array(data))

def main():

   data_DA4032  = read_txt(fn_DA4032)
   data_FCST4032  = read_txt(fn_FCST4032)

   data_DA6528 = read_txt(fn_DA6528)
   data_FCST6528 = read_txt(fn_FCST6528)

   data_DA8192 = read_txt(fn_DA8192)

   data_DA = data_DA6528
   data_FCST = data_FCST6528

   import matplotlib.pyplot as plt

   fig, (ax1,ax2) = plt.subplots(2, 1, figsize=(7.0, 12.0))
   fig.subplots_adjust(left=0.15, bottom=0.07, right=0.98, top=0.92, wspace=0.2, hspace=0.0)

   ymin_DA = 0.0
   ymax_DA = 60.0

   xvals_DA = np.arange(data_DA.size) + 1
   xvals_FCST = np.arange(data_FCST.size) + 1

   xmin_DA = 0.0 #np.min(xvals_DA)
   xmax_DA = 3.0 #np.max(xvals_DA)

#   ax1.plot(xvals_DA, data_DA, marker='o',markersize=6,color='r',alpha=0.8)

   ax1.tick_params(labelsize=16)
   ax1.tick_params(axis='x',which='both', bottom=False, top=False, labelbottom=False) 
   ax2.tick_params(labelsize=16)


   data_DA = [np.mean(data_DA4032), np.mean(data_DA6528), np.mean(data_DA8192)]
   data_FCST = [np.mean(data_FCST4032), np.mean(data_FCST6528), np.nan]
   names = ['4032 nodes\n(20 DA cycles +\n11 forecasts)', '6528 nodes\n(64 cycles +\n50 forecasts)', '8192 nodes\n(17 DA cycles)']
   xvals = [0.5,1.5,2.5]

   ax1.bar(xvals, data_DA, color='k',alpha=0.5,align='center')
   ax2.bar(xvals, data_FCST, color='b',alpha=0.5,align='center')

   ax2.set_xticks(xvals)
   ax2.set_xticklabels(names, minor=False, fontsize=14)



   ax1.set_ylim(ymin_DA,ymax_DA)

   ymin_FCST = 0.0
   ymax_FCST = 900.0

   ax2.set_ylim(ymax_FCST,ymin_FCST)

   ax1.set_xlim(xmin_DA,xmax_DA)
   ax2.set_xlim(xmin_DA,xmax_DA)

   ax1.set_ylabel("Data assimilation cycle (s)", fontsize=18)
   ax2.set_ylabel("30-min forecast (s)", fontsize=18)


   tit = "Computation time"# (6528 nodes) \n64 DA cycles + 50 forecasts"
   fig.suptitle(tit, fontsize=20)

   odir = "png/HPC201908"
   ofig = "1p_scalability.png"

   if not quick:

      os.makedirs(odir, exist_ok=True)
      plt.savefig(os.path.join(odir,ofig), bbox_inches="tight", pad_inches = 0.1)
      plt.clf()
      plt.close('all')
   else:
      plt.show()
      plt.clf()
      plt.close('all')

main()

