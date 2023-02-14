import unittest
import requests
from astropy.coordinates import SkyCoord
from datetime import datetime
import julian


class TestAPI(unittest.TestCase):
    def setUp(self) -> None:
        self.TOLERANCE = 0.01
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

            #Get RA and DEC from JPL Horizons API
            jpl_response = requests.get(f'{self.JPL_URL}?format=json&COMMAND=\'TLE\'&TLE={jpl_tle}\
            &MAKE_EPHEM=\'YES\'&TLIST=\'{self.jd}\'&EPHEM_TYPE=\'OBSERVER\'&CENTER=\'000@399\'').json()
            lines = jpl_response['result'].splitlines()
            index = lines.index('$$SOE')
            jpl_response_data = lines[index+1].split()
            ra_str = f'{jpl_response_data[2]}h{jpl_response_data[3]}m{jpl_response_data[4]}s'
            dec_str = f'{jpl_response_data[5]}d{jpl_response_data[6]}m{jpl_response_data[7]}s'
            coords = SkyCoord(ra_str, dec_str)
            correct_ra = coords.ra.deg
            correct_dec = coords.dec.deg

            #Get RA and DEC from Apex API
            response = requests.get(f'{self.BASE_URL}tle/{tle}/{self.lat}/{self.long}/{self.jd}').json()
            ra = response['Right Ascension']
            dec = response['Declination']

            #Compare RA and DEC
            self.assertAlmostEqual(ra, correct_ra, delta=self.TOLERANCE)
            self.assertAlmostEqual(dec, correct_dec, delta=self.TOLERANCE)
        print(f'{len(self.NORAD_IDS)} TLE tests passed')

    def testPosition(self):
        for id in self.NORAD_IDS:
            # print(id)
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
            &MAKE_EPHEM=\'YES\'&TLIST=\'{self.jd}\'&EPHEM_TYPE=\'OBSERVER\'&CENTER=\'500\'').json()
            lines = jpl_coords['result'].splitlines()
            index = lines.index('$$SOE')
            jpl_response_data = lines[index+1].split()
            ra_str = f'{jpl_response_data[2]}h{jpl_response_data[3]}m{jpl_response_data[4]}s'
            dec_str = f'{jpl_response_data[5]}d{jpl_response_data[6]}m{jpl_response_data[7]}s'
            coords = SkyCoord(ra_str, dec_str)
            correct_ra = coords.ra.deg
            correct_dec = coords.dec.deg

            #Get RA and DEC from Apex API
            response = requests.get(f'{self.BASE_URL}pos/{x}/{y}/{z}').json() #breaks here
            ra = response['Right Ascension']
            dec = response['Declination']

            #Compare RA and DEC
            self.assertAlmostEqual(ra, correct_ra, delta=self.TOLERANCE)
            self.assertAlmostEqual(dec, correct_dec, delta=self.TOLERANCE)
        print(f'{len(self.NORAD_IDS)} Position tests passed')



if __name__ == '__main__':
    unittest.main()