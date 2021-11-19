import numpy as np

TrackName="DT173"

mask = np.loadtxt("mask_"+TrackName+".xyz")
inc = np.loadtxt("incidenceAngle_"+TrackName+".xyz")
inc[mask[:,2]==0,2]=0
masked_inc = inc[inc[:,2]!=0]

azm = np.loadtxt("azimuthAngle_"+TrackName+".xyz")
azm[mask[:,2]==0,2]=0
masked_azm = azm[azm[:,2]!=0]


#losx=np.cos(masked_inc[:,2])*-1*np.cos(masked_azm[:,2])
#losy=np.cos(masked_inc[:,2])*-1*np.sin(masked_azm[:,2])
#losz=-1*np.sin(masked_inc[:,2])

losx=np.sin(masked_inc[:,2]*np.pi/180)*np.cos(masked_azm[:,2]*np.pi/180)
losy=np.sin(masked_inc[:,2]*np.pi/180)*np.sin(masked_azm[:,2]*np.pi/180)
losz=np.cos(masked_inc[:,2]*np.pi/180)


lon = masked_azm[:,0]
lat = masked_azm[:,1]

lonlat_losx=np.vstack([lon, lat, losx]).transpose()
lonlat_losy=np.vstack([lon, lat, losy]).transpose()
lonlat_losz=np.vstack([lon, lat, losz]).transpose()

np.savetxt("masked_losx_"+TrackName+".xyz", lonlat_losx, delimiter = ' ',fmt='%g')
np.savetxt("masked_losy_"+TrackName+".xyz", lonlat_losy, delimiter = ' ',fmt='%g')
np.savetxt("masked_losz_"+TrackName+".xyz", lonlat_losz, delimiter = ' ',fmt='%g')

