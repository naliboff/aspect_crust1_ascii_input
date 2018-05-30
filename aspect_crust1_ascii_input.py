# Author:  John Naliboff
#
# Overview
#   A python script to process CRUST1.0 data and write ascii files
#   that can be used in ASPECT. At present, the code can only
#   write files for compositional initial conditions. 
#

# Input file information
#   Running the code requires specifying various parameters in a separate python
#   file that is specified when the code is executed. These parameters include
#   information about the geographical extent for CRUST1.0 data extraction,
#   model height/base/reference surface, vertical sampling resolution, coordinate
#   system (cartesian or spherical) and model name.
#        Description                   Units    name
#     1. Minimum longitude            (degrees; lon1)
#     2. Maximum longitude            (degrees; lon2)
#     3. Minimum latitude             (degrees; lat1)
#     4. Maximum latitude             (degrees; lat2)
#     5. Maximum colatitude           (degres; 90. - lat1)
#     6. Minimum colatitude           (degrees; 90. - lat2)
#     7. Model height                 (meters; top)
#     8. Model reference surface      (meters; ref)
#     9. Model base                   (meters; bottom)
#    10. Radius at base of model      (meters; radb)
#    11. Vertical sampling resolution (meters; res)
#    12. Coordinate system            (none; crs)
#    13. Model name                   (none; name)

# General notes:
#   1. No optimized for efficient numpy array operations
#   2. Layers for 'sticky air' and the asthenosphere are currently added and 
#      defined based on the input parameters for model height and depth.

# Instructions (run time)
#   The code should be executed from a terminal with the following general syntax:
#     python aspect_crust1_ascii_input.py full/path/to/input/file/ input_file_name
#
#   For example, the command to run 'tests/test.py' from the current directory is
#     python aspect_crust1_ascii_input.py "$PWD"/examples/ example

# Load modules
import numpy as np; import sys; import pickle; import os;

#------------------------------------------------------------------------------

# Main function
def main():


  # Get input file variables
  input_file_directory = str(sys.argv[1])
  input_file_name = str(sys.argv[2])

  sys.path.insert(0, input_file_directory)

  # Get user input variables
  exec('from ' + input_file_name + ' import *')

  # Caculate geographic grid coordinate values from input file parameters
  rad_pts_asp,rad_grd_asp,lon_pts_asp,col_pts_asp,geo_pts_asp = grid_values(top,bot,res,lon1,lon2,col1,col2)

  # Generate crust1 longitude-colatitude coordinate array
  lon_col_cr1 = crust1_coord()

  # Generate compositional field ascii data file
  compositional_data(crs,top,bot,ref,lon_col_cr1,lon1,lon2,col1,col2,name,rad_pts_asp, \
                       lon_pts_asp,col_pts_asp,rad_grd_asp,input_file_directory)
    
  # Remove .pyc file
  os.system('rm ' + input_file_directory + '/*.pyc')


#----------------------------------------------------------------------------

def grid_values(top,bot,res,lon1,lon2,col1,col2):

  # Calculate number of radial points
  rad_pts_asp = int((top-bot)/res + 1.0);

  # Create zeros array for radial positions
  rad_grd_asp = np.zeros([rad_pts_asp,1]);

  # Fill zeros array with radial positions
  for i in range(rad_pts_asp):
    rad_grd_asp[i] = bot + i*res

  # Number of user defined longitude and colatitude points
  lon_pts_asp = int(lon2 - lon1 + 1.0);
  col_pts_asp = int(col1 - col2 + 1.0);

  # Total number of user defined longitude and colatitude (geographic) points
  geo_pts_asp = lon_pts_asp * col_pts_asp

  # Return values
  return rad_pts_asp, rad_grd_asp, lon_pts_asp, col_pts_asp, geo_pts_asp

#----------------------------------------------------------------------------

def compositional_data(crs,top,bot,ref,lon_col_cr1,lon1,lon2,col1,col2,name,rad_pts_asp, \
                       lon_pts_asp,col_pts_asp,rad_grd_asp, input_file_directory):

  # Load crustal thickness bounds (depth (km, convert to m) at top of 9 layers:
  #   water, ice, sed1, sed2, sed3, crust1, crust2, crust3, mantle  
  dep_lay_cr1 = np.loadtxt('data/crust1/crust1.bnds')*1.e3;

  # Load layer density (g/cm3 convert to kg/m3) for 9 layers:
  #   water, ice, sed1, sed2, sed3, crust1, crust2, crust3, mantle  
  den_lay_cr1 = np.loadtxt('data/crust1/crust1.rho')*1.e3;

  # Make depth and density (10 kg/m3) layer for air
  dep_lay_air = np.zeros([dep_lay_cr1.shape[0],1]) + (top-ref);
  den_lay_air = np.zeros([dep_lay_cr1.shape[0],1]) + 10.;

  # Add air depth and density layer arrays (1 row) to depth and density layer arrays (air is now top row)
  dep_lay_cr1 = np.concatenate((dep_lay_air,dep_lay_cr1),axis=1); 
  den_lay_cr1 = np.concatenate((den_lay_air,den_lay_cr1),axis=1);

  # Make depth and density layers for "asthenosphere" (base of lithosphere)
  dep_lay_ast = np.zeros([dep_lay_cr1.shape[0],1]) + (bot-ref);
  den_lay_ast = np.zeros([den_lay_cr1.shape[0],1]) + 3300.;

  # Add asthenosphere depth layer array (1 row) to depth layer array
  dep_lay_cr1 = np.concatenate((dep_lay_cr1,dep_lay_ast),axis=1);
  den_lay_cr1 = np.concatenate((den_lay_cr1,den_lay_ast),axis=1);

  # Convert layer depths to radius
  rad_lay_cr1 = dep_lay_cr1 + ref;

  # Create variable that combines coordinate ("cor"), layer radii ("rad") and densities ("den")
  cor_rad_den_cr1 = np.concatenate((lon_col_cr1,rad_lay_cr1,den_lay_cr1),axis=1);

  # Find which CRUST1 points lie within defined geographical bounds
  inds = np.where( (cor_rad_den_cr1[:,0]>=lon1)
                 & (cor_rad_den_cr1[:,0]<=lon2)
                 & (cor_rad_den_cr1[:,1]<=col1)
                 & (cor_rad_den_cr1[:,1]>=col2) )

  # Create new layer radius and density array based on defined geographic bounds
  tmp = cor_rad_den_cr1[inds,:]; cor_rad_den_asp = tmp[0,:,:]; 

  # Write ASPECT composition ascii input file
  if crs == 'car':
    write_car_com_output(lon1,col2,name,rad_pts_asp,lon_pts_asp,col_pts_asp, \
                                rad_grd_asp,cor_rad_den_asp,input_file_directory)
  

  elif crs == 'sph':
    write_sph_com_output(radb,name,rad_pts_asp,lon_pts_asp,col_pts_asp, \
                               rad_grd_asp,cor_rad_den_asp, input_file_directory)

#----------------------------------------------------------------------------

def write_sph_com_output(radb,name,rad_pts_asp,lon_pts_asp,col_pts_asp, \
                               rad_grd_asp,cor_rad_den_asp,input_file_directory):

  # Open file to write results to
  outfile=open(input_file_directory+'/crust1_'+name+'.txt','w')

  # Write header line
  print >> outfile, '# POINTS: %-i %i %i'% (rad_pts_asp,lon_pts_asp,col_pts_asp)

  # Convert coordinates from degrees to radians
  cor_rad_den_asp[:,0:2] = np.radians(cor_rad_den_asp[:,0:2])

  # Sort coordinate arrays so that colatitude is in ascending order
  inds = np.lexsort((cor_rad_asp[:,0],cor_rad_den_asp[:,1])); 
  cor_rad_den_asp = cor_rad_den_asp[inds,:];
  
  # Loop through lon/colat points
  for i in range(lon_pts_asp*col_pts_asp):
  
    # Loop through radial points
    for j in range(rad_pts_asp):
      
      # Create variable containing current row of cor_rad_asp variable
      cur = np.copy(cor_rad_den_asp[i,:]);
      
      # Calculate thickness of each layer at current point
      thk = cur[2:12] - cur[3:13];
      
      # If the thickness of a layer is 0, set its radius to 1.e9 meters
      for k in range(thk.shape[0]):
        if thk[k]==0.: cur[k+2]=1.e9
      
      # Find index of closest composition to current depth
      inds = np.argmin(np.absolute(cur[2:12]-rad_grd_asp[j]));
      
      # Define variables for output longitiude, colatitude and density
      lon = cur[0]; col = cur[1]; den = cur[13+inds];
      
      # Write values to a file
      print >> outfile, '%-12.2f' % (rad_grd_asp[j]+ radb),
      print >> outfile, '%-9.6f'  % (lon),
      print >> outfile, '%-9.6f'  % (col),
      print >> outfile, '%3.1f'   % (inds),
      print >> outfile, '%10.1f'  % (den)

  # Close output file
  outfile.close()

#------------------------------------------------------------------------------

def write_car_com_output(lon1,col2,name,rad_pts_asp,lon_pts_asp,col_pts_asp, \
                                rad_grd_asp,cor_rad_den_asp,input_file_directory):

  # Open data file
  outfile=open(input_file_directory+'/crust1_'+name+'.txt','w')

  # Write header line for different cases
  # 2D longitudinal profile
  if lon_pts_asp>1 and col_pts_asp==1:
    print >> outfile, '# POINTS: %-i %i'% (lon_pts_asp,rad_pts_asp)
  # 2D latitudinal profile
  elif lon_pts_asp==1 and col_pts_asp>1:
    print >> outfile, '# POINTS: %-i %i'% (col_pts_asp,rad_pts_asp)
  # 3D
  elif lon_pts_asp>1 and col_pts_asp>1:
    print >> outfile, '# POINTS: %-i %i %i'% (lon_pts_asp,col_pts_asp,rad_pts_asp)

  # Loop through radial points
  for j in range(rad_pts_asp):

    # Loop through lon/colat points
    for i in range(lon_pts_asp*col_pts_asp):

      # Create variable containing current row of cor_rad_asp variable
      cur = np.copy(cor_rad_den_asp[i,:]); 

      # Calculate thickness of each layer at current point
      thk = cur[2:12] - cur[3:13];

      # If the thickness of a layer is 0, set its radius to 1.e9 meters
      for k in range(thk.shape[0]):
        if thk[k]==0.: cur[k+2]=1.e9
      
      # Find index of closest composition to current depth
      inds = np.argmin(np.absolute(cur[2:12]-rad_grd_asp[j]));
      
      # Define variables for output longitiude, colatitude and density
      lon = cur[0]; col = cur[1]; den = cur[13+inds];

      # Write values to a file
      if lon_pts_asp>1:
        print >> outfile, '%-12.2f' % ((lon-lon1)*111.e3),
      if col_pts_asp>1:
        print >> outfile, '%-12.2f' % ((col-col2)*111.e3),
      print >> outfile, '%-12.2f' % (rad_grd_asp[j]),
      print >> outfile, '%3.1f' % (inds),
      print >> outfile, '%10.1f' % (den)

  # Close output file
  outfile.close()

#------------------------------------------------------------------------------

def crust1_coord():
  
  # Coordinate parameters
  lon_min_cr1 = -179.5    # Minimum longitude
  lon_max_cr1 =  179.5    # Maximum longitude
  lon_pts_cr1 =  360      # Total longitude points (1 degree spacing)
  lat_min_cr1 =  -89.5    # Minimum latitude
  lat_max_cr1 =   89.5    # Maximum latitude
  lat_pts_cr1 =   180     # Total latitude points (1 degree spacing)

  # Longitude coordinate array 
  lon_val_cr1 = np.linspace(lon_min_cr1,lon_max_cr1,lon_pts_cr1)

  # Latitude coordinate array
  lat_val_cr1 = np.linspace(lat_max_cr1,lat_min_cr1,lat_pts_cr1)

  # Zeros array for CRUST1 format longitude-latitude array
  lon_lat_cr1 = np.zeros(shape=(lon_pts_cr1*lat_pts_cr1,2))

  # Fill CRUST1 format zeros array with longitude-latitude values
  count = 0
  for i in range(lat_pts_cr1):
    for j in range(lon_pts_cr1):
      lon_lat_cr1[count,0] = lon_val_cr1[j]
      lon_lat_cr1[count,1] = lat_val_cr1[i]
      count = count + 1

  # Create array for longitude and colatitude values in CRUST1 array format
  lon_col_cr1 = np.copy(lon_lat_cr1[:,:])

  # Modify longitude to range from 0-360 degrees 
  for i in range(lon_pts_cr1*lat_pts_cr1):
    if lon_col_cr1[i,0] < 0.:
      lon_col_cr1[i,0] = lon_col_cr1[i,0] + 360.;

  # Convert latitude to colatitude
  lon_col_cr1[:,1] = 90. - lon_col_cr1[:,1]

  # Return array
  return lon_col_cr1

#------------------------------------------------------------------------------

# Call main function
main()

#
