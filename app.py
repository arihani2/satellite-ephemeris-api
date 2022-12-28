from flask import Flask, redirect
from skyfield.api import EarthSatellite
from skyfield.api import load, wgs84
import numpy as np


app = Flask(__name__)

#Redirects user to the Center for the Protection of Dark and Quiet Sky homepage
@app.route('/')
def root():
    return redirect('https://cps.iau.org/')

@app.route('/tle/<tle>/<latitude>/<longitude>')
def read_tle_string(tle, latitude, longitude):
    '''
    Returns the Right Ascension and Declination relative to the observer's coordinates
    for the given satellite's Two Line Element Data Set at epoch

    Parameters
    ---------
    tle: 'str'
        The Two Line Element set
    latitude: 'float'
        The observers latitude coordinate (positive value represents north, negative value represents south)
    longitude: 'float'
        The observers longitude coordinate (positive value represents east, negatie value represents west)

    Returns
    -------
    Right Ascension: 'float'
        The right ascension of the satellite in degrees relative to the observer
    Declination: 'float'
        The declination of the satellite in degrees relative to the observer
    '''
    #Get rid of ASCII representation for space
    tle = tle.replace("%20", ' ')
    #Retrieve the two lines
    u,w = tle[:69], tle[70:]

    #Cast the latitude and longitude to floats (request parses as a string)
    lat = float(latitude)
    lon = float(longitude)

    #This is the skyfield implementation
    ts = load.timescale()
    satellite = EarthSatellite(u,w,ts = ts)

    #Get current position and find topocentric ra and dec
    currPos = wgs84.latlon(lat, lon)
    # Both times t are correct, one returns relative at runtime, other returns relative at satellite epoch
    # t = ts.now()
    t = ts.ut1_jd(satellite.model.jdsatepoch)
    difference = satellite - currPos
    topocentric = difference.at(t)
    ra, dec, distance = topocentric.radec()

    return {
        'Right Ascension': ra.dstr(warn = False),
        'Declination': dec.dstr()
    }

@app.route('/tle_file/<file_path>')
def read_tle_from_file(file_path:str):
    '''
    Returns dictionary of relative Right Ascension and Declination for all
    epoch times in Two Line Element data set text file

    ***This function is not complete. Will need a file upload of some sort to properly
    parse the data. This is a future project.***

    Parameters
    ---------
    file_path: 'file_path'
        The Two Line Element text file. This is the relative path

    Returns
    -------
    Dictionary of values
    '''
    to_return = []
    f = open(file_path)
    words = f.read()
    single_lines = words.split('\n')
    lat_lon = single_lines[0]
    for i in range(1, len(single_lines)-1, 2):
        u = single_lines[i]
        w = single_lines[i+1]
        TLE = f'{u} {w}, {lat_lon}'
        to_return.append(read_tle_string(TLE))

    return to_return


@app.route('/pos/<pos_string>')
def get_ephemeris(pos_string: str):
    '''
    Returns the geocentric Right Ascension and Declination of the orbiting 
    mass given a position vector

    Parameters
    ---------
    pos_string: 'str'
        The position and velocity vector. Values are in order of (x,y,z, dx/dt, dy/dt, dz/dt).
        Values are comma separated.

    Returns
    -------
    Right Ascension: 'float'
        The geocentric right ascension of the satellite in degrees
    Declination: 'float'
        The geocentric declination of the satellite in degrees
    '''
    pos_vector = pos_string.split(',')
    position = np.array([float(element) for element in pos_vector])
    ra, dec = icrf2radec(position)

    return {
        'Right Ascension': ra,
        'Declination': dec
    }


def icrf2radec(pos, deg=True):
    """
    Convert ICRF xyz to Right Ascension and Declination.
    Geometric states on unit sphere, no light travel time/aberration correction.
    
    Parameters:
    -----------
    pos ... real, dim=[n, 3], 3D vector of unit length (ICRF)
    deg ... True: angles in degrees, False: angles in radians
    Returns:
    --------
    ra ... Right Ascension [deg]
    dec ... Declination [deg]
    """
    norm=np.linalg.norm
    array=np.array
    arctan2=np.arctan2
    arcsin=np.arcsin
    rad2deg=np.rad2deg
    modulo=np.mod
    pix2=2.*np.pi
    
    if(pos.ndim>1):
        r=norm(pos,axis=1)
        xu=pos[:,0]/r
        yu=pos[:,1]/r
        zu=pos[:,2]/r
    else:
        r=norm(pos)
        xu=pos[0]/r
        yu=pos[1]/r
        zu=pos[2]/r
    
    phi=arctan2(yu,xu)
    delta=arcsin(zu)
    
    if(deg):
        ra = modulo(rad2deg(phi)+360,360)
        dec = rad2deg(delta)
    else:
        ra = modulo(phi+pix2,pix2)
        dec = delta
    
    return ra, dec