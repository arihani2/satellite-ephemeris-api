import unittest
import requests
from astropy.coordinates import SkyCoord
from datetime import datetime
from random import randint
import julian


class TestAPI(unittest.TestCase):
    def setUp(self) -> None:
        self.RADEC_TOLERANCE = 0.01
        self.AZALT_TOLERANCE = 0.1
        self.jd = julian.to_jd(datetime.now(), fmt='jd')
        self.lat = 51.4773207 #Greenwich Latitude
        self.long = 0.0 #Greenwich Longitude
        self.BASE_URL = 'http://apexgroup.web.illinois.edu/ephemeris/'
        self.JPL_URL = 'https://ssd.jpl.nasa.gov/api/horizons.api/'
        self.NORAD_IDS = [55269, 55270, 55271, 55272, 55273, 55274, 55275, 55276, 55277, 55278,\
                          55279, 55280, 55281, 55282, 55283, 55284, 55285, 55286, 55287, 55288]

    def testTLE(self):
        for id in self.NORAD_IDS:
            #Get TLE from Celestrak
            CELESTRAK = f'https://celestrak.org/NORAD/elements/gp.php?CATNR={id}&FORMAT=2LE'
            tle = requests.get(CELESTRAK).text

            #Get TLE information in proper format for each API call
            jpl_tle = tle.replace('\r\n', '%0A').strip()
            tle = tle.replace('\r\n', ' ').strip()

            #Get RA, DEC, Altitude and Azimuth from JPL Horizons API
            jpl_response = requests.get(f'{self.JPL_URL}?format=json&COMMAND=\'TLE\'&TLE={jpl_tle}\
            &MAKE_EPHEM=\'YES\'&TLIST=\'{self.jd}\'&EPHEM_TYPE=\'OBSERVER\'&CENTER=\'000@399\'&ANG_FORMAT=\'DEG\'').json()
            lines = jpl_response['result'].splitlines()
            index = lines.index('$$SOE')
            jpl_response_data = lines[index+1].split()
            correct_ra = float(jpl_response_data[2])
            correct_dec = float(jpl_response_data[3])
            correct_az = float(jpl_response_data[8])
            correct_alt = float(jpl_response_data[9])

            #Get RA, DEC, Altitude, and Azimuth from Apex API
            response = requests.get(f'{self.BASE_URL}tle/{tle}/{self.lat}/{self.long}/{self.jd}').json()
            ra = response['Right Ascension']
            dec = response['Declination']
            alt = response['Altitude']
            az = response['Azimuth']

            #Compare RA, DEC, Altitude, and Azimuth
            self.assertAlmostEqual(ra, correct_ra, delta=self.RADEC_TOLERANCE)
            self.assertAlmostEqual(dec, correct_dec, delta=self.RADEC_TOLERANCE)
            self.assertAlmostEqual(alt, correct_alt, delta=self.AZALT_TOLERANCE)
            self.assertAlmostEqual(az, correct_az, delta=self.AZALT_TOLERANCE)
        print(f'\n{len(self.NORAD_IDS)} TLE tests passed\n')

    def testPosition(self):
        for id in self.NORAD_IDS:
            #Get TLE formatted for each API call
            CELESTRAK = f'https://celestrak.org/NORAD/elements/gp.php?CATNR={id}&FORMAT=2LE'
            tle = requests.get(CELESTRAK).text
            jpl_tle = tle.replace('\r\n', '%0A').strip()

            #Retrieve vector response from JPL Horizons API
            #USE SKYFIELD TO GET CORRECT RA AND DEC possibly 
            jpl_vec = requests.get(f'{self.JPL_URL}?format=json&EPHEM_TYPE=\'VECTORS\'&COMMAND=\'TLE\'&TLE={jpl_tle}\
            &MAKE_EPHEM=\'YES\'&TLIST=\'{self.jd}\'&TLIST_TYPE=\'JD\'&VEC_TABLE=1&TIME_TYPE=\'UT\'&REF_PLANE=\'FRAME\'').json()
            index_vec = jpl_vec['result'].splitlines().index('$$SOE')
            str_resp = (jpl_vec['result'].splitlines()[index_vec+2].split('E'))
            x_str = str_resp[0]
            y_str = str_resp[1]
            z_str = str_resp[2]
            x_pow = float(str_resp[1][1:3])
            y_pow = float(str_resp[2][1:3])
            z_pow = float(str_resp[3][1:3])
            x = float(x_str[4:]) * 10 ** x_pow
            y = float(y_str[7:]) * 10 ** y_pow
            z = float(z_str[7:]) * 10 ** z_pow

            #Retrieve coordinate response from JPL Horizons API
            jpl_coords = requests.get(f'{self.JPL_URL}?format=json&COMMAND=\'TLE\'&TLE={jpl_tle}\
            &MAKE_EPHEM=\'YES\'&TLIST=\'{self.jd}\'&EPHEM_TYPE=\'OBSERVER\'&CENTER=\'500\'&ANG_FORMAT=\'DEG\'').json()
            lines = jpl_coords['result'].splitlines()
            index = lines.index('$$SOE')
            jpl_response_data = lines[index+1].split()
            correct_ra = float(jpl_response_data[2])
            correct_dec = float(jpl_response_data[3])

            #Get RA and DEC from Apex API
            response = requests.get(f'{self.BASE_URL}pos/{x},{y},{z}').json() 
            ra = response['Right Ascension']
            dec = response['Declination']

            #Compare RA and DEC
            self.assertAlmostEqual(ra, correct_ra, delta=self.RADEC_TOLERANCE)
            self.assertAlmostEqual(dec, correct_dec, delta=self.RADEC_TOLERANCE)
        print(f'\n{len(self.NORAD_IDS)} Position tests passed\n')
    
    def testElevation(self):
        for id in self.NORAD_IDS:
            #Get TLE from Celestrak
            CELESTRAK = f'https://celestrak.org/NORAD/elements/gp.php?CATNR={id}&FORMAT=2LE'
            tle = requests.get(CELESTRAK).text
            elevation = randint(0, 50000)
            elevation_km = elevation / 1000

            #Get TLE information in proper format for each API call
            jpl_tle = tle.replace('\r\n', '%0A').strip()
            tle = tle.replace('\r\n', ' ').strip()

            #Get RA and DEC from JPL Horizons API
            jpl_response = requests.get(f'{self.JPL_URL}?format=json&COMMAND=\'TLE\'&TLE={jpl_tle}\
            &MAKE_EPHEM=\'YES\'&TLIST=\'{self.jd}\'&EPHEM_TYPE=\'OBSERVER\'&CENTER=\'coord\'\
            &SITE_COORD=\'{self.long},{self.lat},{elevation_km}\'&COORD_TYPE=\'GEODETIC\'&ANG_FORMAT=\'DEG\'').json()
            lines = jpl_response['result'].splitlines()
            index = lines.index('$$SOE')
            jpl_response_data = lines[index+1].split()
            correct_ra = float(jpl_response_data[2])
            correct_dec = float(jpl_response_data[3])

            #Get RA and DEC from Apex API
            response = requests.get(f'{self.BASE_URL}tle/{tle}/{self.lat}/{self.long}/{self.jd}/&elevation={elevation}').json()
            ra = response['Right Ascension']
            dec = response['Declination']

            #Compare RA and DEC
            self.assertAlmostEqual(ra, correct_ra, delta=self.RADEC_TOLERANCE)
            self.assertAlmostEqual(dec, correct_dec, delta=self.RADEC_TOLERANCE)
        print(f'\n{len(self.NORAD_IDS)} Elevation tests passed\n')



if __name__ == '__main__':
    unittest.main()