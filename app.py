from flask import Flask, redirect
from astropy.coordinates import SkyCoord
from skyfield.api import EarthSatellite
from skyfield.api import load, wgs84
import numpy as np


app = Flask(__name__)

#Redirects user to the Center for the Protection of Dark and Quiet Sky homepage
@app.route('/')
def root():
    return redirect('https://cps.iau.org/')

@app.route('/tle/<tle>/<latitude>/<longitude>/<julian_date>')
def read_tle_string(tle, latitude, longitude, julian_date):
    '''
    Returns the Right Ascension and Declination relative to the observer's coordinates
    for the given satellite's Two Line Element Data Set at inputted Julian Date.

    **Please note, for the most accurate results, an inputted Julian Date close to the TLE epoch is necessary.

    Parameters
    ---------
    tle: 'str'
        The Two Line Element set
    latitude: 'float'
        The observers latitude coordinate (positive value represents north, negative value represents south)
    longitude: 'float'
        The observers longitude coordinate (positive value represents east, negatie value represents west)
    julian_date: 'float
        UT1 Universal Time Julian Date. An input of 0 will use the TLE epoch.

    Returns
    -------
    Right Ascension: 'float'
        The right ascension of the satellite relative to observer coordinates in ICRS reference frame in degrees. Range of response is [0,360)
    Declination: 'float'
        The declination of the satellite relative to observer coordinates in ICRS reference frame in degrees. Range of response is [-90,90]
    '''
    #Get rid of ASCII representation for space
    tle = tle.replace("%20", ' ')
    #Retrieve the two lines
    u,w = tle[:69], tle[70:]

    #Cast the latitude, longitude, and jd to floats (request parses as a string)
    lat = float(latitude)
    lon = float(longitude)
    jd = float(julian_date)

    #This is the skyfield implementation
    ts = load.timescale()
    satellite = EarthSatellite(u,w,ts = ts)

    #Get current position and find topocentric ra and dec
    currPos = wgs84.latlon(lat, lon)
    # t = ts.now()
    # Set time to satellite epoch if input jd is 0, otherwise time is inputted jd
    if jd == 0: t = ts.ut1_jd(satellite.model.jdsatepoch)
    else: t = ts.ut1_jd(jd)

    difference = satellite - currPos
    topocentric = difference.at(t)
    ra, dec, distance = topocentric.radec()

    deg = SkyCoord(ra.hstr(), dec.dstr())

    return {
        'Right Ascension': deg.ra.degree,
        'Declination': deg.dec.degree
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


@app.route('/pos/<x>/<y>/<z>')
def get_ephemeris(x, y, z):
    '''
    Returns the geocentric Right Ascension and Declination of the orbiting 
    mass given the geocentric position vector

    Parameters
    ---------
    x: 'float'
        The x position of the orbiting mass in km
    y: 'float'
        The y position of the orbiting mass in km
    z: 'float'
        The z position of the orbiting mass in km
        

    Returns
    -------
    Right Ascension: 'float'
        The geocentric right ascension of the satellite in degrees
    Declination: 'float'
        The geocentric declination of the satellite in degrees
    '''
    position = np.array([float(x), float(y), float(z)])
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