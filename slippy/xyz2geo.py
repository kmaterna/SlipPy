import numpy as np 
import collections

Map_boundaries_tuple=collections.namedtuple("Map_boundaries_tuple",["lon_0",
    "lat_0","llcrnrlat","llcrnrlon","urcrnrlat","urcrnrlon","proj","resolution"]);


def geo2xyz(lons, lats, lon0, lat0):
  # returns units of meters
  if len(np.atleast_1d(lons))==1:
  	x_array = (lons-lon0)*111000.0*(np.cos(np.deg2rad(lat0)));
  	y_array = (lats-lat0)*111000.0;
  else:
    x_array = np.zeros(np.shape(lons));
    y_array = np.zeros(np.shape(lats));

    for i in range(len(lons)):
      x_array[i]=(lons[i]-lon0)*111000.0*(np.cos(np.deg2rad(lat0)));
      y_array[i]=(lats[i]-lat0)*111000.0;

  return x_array,y_array;


def xyz2geo(x, y, lon0, lat0):
  # takes units of meters
  if len(np.atleast_1d(x))==1:
    lon_array = lon0 + (x/(111000.0*np.cos(np.deg2rad(lat0))));
    lat_array = lat0 + (y/111000.0);
  else:
    lon_array = np.zeros(np.shape(x));
    lat_array = np.zeros(np.shape(y));

    for i in range(len(x)):
      lon_array[i] = lon0 + (x[i]/(111000.0*np.cos(np.deg2rad(lat0))));
      lat_array[i] = lat0 + (y[i]/111000.0);

  return lon_array, lat_array;


def geodetic_to_cartesian(pos_geo, collection):
  ''' 
  Parameters
  ----------
    pos_geo : (...,D) array
      array of geodetic positions. The first and second component of 
      the last axis are longitude and latitude.  The last axis can 
      have additional components (e.g. height) and they will be 
      returned unaltered.  

  Returns
  -------
    pos_cart : (...,D) array
  '''    
  pos_cart = np.array(pos_geo,copy=True)
  pos_geo = np.asarray(pos_geo)
  pos_x,pos_y = geo2xyz(pos_geo[...,0],pos_geo[...,1], collection.lon_0, collection.lat_0)  
  pos_x = np.asarray(pos_x)
  pos_y = np.asarray(pos_y)
  pos_cart[...,0] = pos_x
  pos_cart[...,1] = pos_y
  return pos_cart



def cartesian_to_geodetic(pos_cart,collection):
  ''' 
  Parameters
  ----------
    pos_cart : (...,D) array
      array of cartesian positions 
    
  Returns
  -------
    pos_geo : (...,D) array  
  '''
  pos_geo = np.array(pos_cart,copy=True)
  pos_cart = np.asarray(pos_cart)
  pos_lon,pos_lat = xyz2geo(pos_cart[...,0],pos_cart[...,1],collection.lon_0, collection.lat_0)  
  pos_lon = np.asarray(pos_lon)
  pos_lat = np.asarray(pos_lat)
  pos_geo[...,0] = pos_lon
  pos_geo[...,1] = pos_lat
  return pos_geo



def create_default_basemap(lon_lst,lat_lst, proj='M5i', resolution='i'):
  ''' 
  creates a named tuple that bounds lat_lst and lon_lst
  '''
  print("Creating collection of map boundaries, for plotting with GMT")
  if (len(lon_lst) == 0) | (len(lat_lst) == 0):
    return Map_boundaries_tuple(lon_0=-90.0, lat_0=41.0, llcrnrlat = 26.0, llcrnrlon=-128.0, urcrnrlat=48.0, urcrnrlon=-53.0, proj='M5i');
    
  lon_buff = (max(lon_lst) - min(lon_lst))/20.0
  lat_buff = (max(lat_lst) - min(lat_lst))/20.0
  if lon_buff < 0.15:
    lon_buff = 0.15

  if lat_buff < 0.15:
    lat_buff = 0.15

  llcrnrlon = min(lon_lst) - lon_buff
  llcrnrlat = min(lat_lst) - lat_buff
  urcrnrlon = max(lon_lst) + lon_buff
  urcrnrlat = max(lat_lst) + lat_buff
  lon_0 = (llcrnrlon + urcrnrlon)/2.0
  lat_0 = (llcrnrlat + urcrnrlat)/2.0
  return Map_boundaries_tuple(lon_0=lon_0, lat_0=lat_0, llcrnrlat = llcrnrlat, 
    llcrnrlon=llcrnrlon, urcrnrlat=urcrnrlat, urcrnrlon=urcrnrlon, proj=proj, resolution=resolution); 

