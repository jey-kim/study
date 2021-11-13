def distance_calculation(lat1,lon1,lat2,lon2):
    '''

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % <distnace_calculation> is a function that calculates           %
    % the angular distance and azimuth between two coordinates       %
    % on the surface of a spherical body.                            %
    % ** See Haversine Formula for details **                        %
    % The syntax of this code is this :                              %
    % ang_dist, az = distance_calculation(lat1,lon1,lat2,lon2)       %
    % The inputs must be in [degree] and the outputs are in [degree] %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    Jeonghyeop Kim (jeonghyeop.kim@gmail.com)
    11/13/2021
    '''

    import numpy as np

    lat1 = lat1*np.pi/180
    lat2 = lat2*np.pi/180

    lon1 = lon1*np.pi/180;
    lon2 = lon2*np.pi/180;

    angular_distance = (180/np.pi)*np.arccos(np.sin(lat1)*np.sin(lat2) \
                          + np.cos(lat1)*np.cos(lat2)*np.cos(lon1-lon2))

    X = np.sin(lon2-lon1)*np.cos(lat2)
    Y = np.cos(lat1)*np.sin(lat2)-np.sin(lat1)*np.cos(lat2)*np.cos(lon2-lon1)
    azimuth = np.arctan2(X,Y)*180/np.pi


    if azimuth < 0:
        azimuth = 360 + azimuth
    

    
    return angular_distance, azimuth   

