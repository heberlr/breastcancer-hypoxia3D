from pyMCDS import pyMCDS
import numpy as np
import matplotlib.pyplot as plt
import scipy.ndimage

initial_index = 1;
last_index = 1;
Lcell_size = 5;
Dcell_size = 1;
Ecell_size = 10;
live_count = np.zeros( last_index+1 );
dead_count = np.zeros( last_index+1 );
endo_count = np.zeros( last_index+1 );
tip_count = np.zeros( last_index+1 );
times = np.zeros( last_index+1 );

for n in range( initial_index,last_index+1 ):
  filename='output'+"%08i"%n+'.xml'
  filenameOut='output/output'+"%08i"%n+'.png'
  mcds=pyMCDS(filename,'output')
  times[n]= mcds.get_time()
  
  cx = mcds.data['discrete_cells']['position_x'];
  cy = mcds.data['discrete_cells']['position_y'];
  cycle = mcds.data['discrete_cells']['cycle_model']
  cycle = cycle.astype( int )
  current_phase = mcds.data['discrete_cells']['current_phase']
  current_phase = current_phase.astype(int)
  elapsed_time_in_phase = mcds.data['discrete_cells']['elapsed_time_in_phase']
  cell_type = mcds.data['discrete_cells']['cell_type']
  cell_type = cell_type.astype(int)
  live = np.argwhere( (cycle < 100) & (cell_type==0) ).flatten()
  dead = np.argwhere( (cycle >= 100) & (cell_type==0) ).flatten()
  endo = np.argwhere( cell_type==1 ).flatten()
  tip = np.argwhere( cell_type==2 ).flatten()
  live_count[n] = len(live)
  dead_count[n] = len(dead)
  endo_count[n] = len(endo)
  tip_count[n] = len(tip)
  
  proteinsX = mcds.data['discrete_cells']['proteins_x']
  proteinsY = mcds.data['discrete_cells']['proteins_y']
  rgb_var = np.vstack((0.5*proteinsX[live],0.5*proteinsY[live],0*proteinsX[live])).T
  
  RadiusSize = 2000.0

  figure, axes = plt.subplots(nrows=1, ncols=3,figsize=(15,5))
  figure.suptitle( '#LiveCell:'+str("%04i"%(live_count[n]))+'  #DeadCell:'+str("%04i"%(dead_count[n]))+ '--  #StalkCell:'+str("%04i"%(endo_count[n]))+ '  #TipCell:'+str("%04i"%(tip_count[n]))+'  Time:' +str("%8.2f"%(n)) + ' hours', size=14)
  plt.subplot(131)
  o2 = mcds.get_concentrations( 'oxygen' );
  X1,Y1 = mcds.get_2D_mesh();
  if (o2.max() > 0):
    v1 = np.linspace(0, o2.max(), 10, endpoint=True)
    plt.contourf(X1,Y1,o2[:,:,0],v1,cmap='binary');
    x = plt.colorbar(ticks=v1)
  else:
    plt.contourf(X1,Y1,o2[:,:,0],cmap='binary');
    x = plt.colorbar()
  plt.xlim(X1.min(), X1.max())
  plt.ylim(Y1.min(), Y1.max())

  
  plt.subplot(132)
  plt.scatter( cx[live],cy[live],c=rgb_var,s=Lcell_size);
  plt.scatter( cx[dead],cy[dead],c='moccasin',s=Dcell_size );
  plt.scatter( cx[endo],cy[endo],c='r',s=Ecell_size );
  plt.scatter( cx[tip],cy[tip],c='b',s=Ecell_size );
  plt.xlim(X1.min(), X1.max())
  plt.ylim(Y1.min(), Y1.max())
  plt.subplots_adjust(left=0.06,right=0.96,bottom=0.08,top=0.82,wspace=0.26)
  
  plt.subplot(133)
  plt.plot(X1[25,:],o2[25,:,2]) 
  plt.show()
  #figure.savefig(filenameOut)
  #plt.close()
