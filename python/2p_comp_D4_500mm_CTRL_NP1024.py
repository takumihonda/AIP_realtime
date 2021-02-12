import numpy as np
import sys
import os


#top = "/data_ballantine02/miyoshi-t/honda/SCALE-LETKF/AIP_D4_VERIFY/20201117/D4_500m_CTRL_NP1024/exp/3660515_cycle_20190824150000"
top = "/data_ballantine02/miyoshi-t/honda/SCALE-LETKF/AIP_D4_VERIFY/20201117/D4_500m_CTRL_NP1024/exp/3662073_cycle_20190824150000"
fn_DA = os.path.join(top, "DA.txt")
fn_FCST = os.path.join(top, "FCST.txt")

quick = False
quick = True

def read_txt(fn):

   try:
      f = open(fn)
   except:
      print("Failed to open")
      sys.exit()
   
   lines = f.read().split('\n')

   data = []
   for line in lines:
      if line == "":
         continue
      data.append(float(line))

   return(np.array(data))

def main(fn_DA, fn_FCST):

   data_DA = read_txt(fn_DA)

   data_FCST = read_txt(fn_FCST)

   import matplotlib.pyplot as plt

   fig, (ax1,ax2) = plt.subplots(2, 1, figsize=(14.0, 13.0))
   fig.subplots_adjust(left=0.08, bottom=0.05, right=0.98, top=0.92, wspace=0.2, hspace=0.0)

   ymin_DA = 0.0
   ymax_DA = 60.0

   xvals_DA = np.arange(data_DA.size) + 1
   xvals_FCST = np.arange(data_FCST.size) + 1

   xmin_DA = 0.0 #np.min(xvals_DA)
   xmax_DA = np.max(xvals_DA)

#   ax1.plot(xvals_DA, data_DA, marker='o',markersize=6,color='r',alpha=0.8)
   ax1.bar(xvals_DA, data_DA, color='k',alpha=0.5,
           align='center')
   ax1.set_ylim(ymin_DA,ymax_DA)

   ymin_FCST = 0.0
   ymax_FCST = 900.0
#   ymax_FCST = 600.0

#   ax2.plot(xvals_FCST, data_FCST, marker='o',markersize=6)
   ax2.bar(xvals_FCST, data_FCST, color='b',alpha=0.5, 
           align="center")#, marker='o',markersize=6)
   ax2.set_ylim(ymax_FCST,ymin_FCST)

   xmax_DA = 64
   xmax_DA = 120
   ax1.set_xlim(xmin_DA,xmax_DA)
   ax2.set_xlim(xmin_DA,xmax_DA)

   ax2.set_xlabel("# of cycle/forecasts", fontsize=18)

   ax1.set_ylabel("Data assimilation cycle (s)", fontsize=18)
   ax2.set_ylabel("30-min forecast (s)", fontsize=18)


   ax1.hlines(xmin=xmin_DA, xmax=xmax_DA, y=30.0, colors='k',
              alpha=0.5, linestyle='dashed')

   print("DA:",np.mean(data_DA[1:]))
   print("FCST:",np.mean(data_FCST))

#   ax1.hlines(xmin=xmin_DA, xmax=xmax_DA, y=np.mean(data_DA), colors='k',
#              alpha=0.5, linestyle='dashed')
  

   ax1.tick_params(labelsize=16)
   ax1.tick_params(axis='x',which='both', bottom=False, top=False,labelbottom=False) 
   ax2.tick_params(labelsize=16)

   tit = "Computation time (6528 nodes) \n64 DA cycles + 50 forecasts"
   fig.suptitle(tit, fontsize=20)

   odir = "png/HPC201911"
   ofig = "6528.png"

   if not quick:

      os.makedirs(odir, exist_ok=True)
      plt.savefig(os.path.join(odir,ofig), bbox_inches="tight", pad_inches = 0.1)
      plt.clf()
      plt.close('all')
   else:
      plt.show()
      plt.clf()
      plt.close('all')

main(fn_DA, fn_FCST)

