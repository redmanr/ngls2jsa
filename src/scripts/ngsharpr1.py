#!/usr/bin/env python2.7

import argparse
from collections import OrderedDict
from ConfigParser import SafeConfigParser
import csv
from datetime import datetime
import logging
import pyfits
import os.path
import re
import shutil
import subprocess
from subprocess import CalledProcessError
import sys

from tools4caom2.logger import logger
from tools4caom2.tapclient import tapclient
from tools4caom2.utdate_string import utdate_string

from tools4caom2.__version__ import version as tools4caom2version
from jcmt2caom2.__version__ import version as jcmt2caom2version


"""
Custom code to prepare files from the JLS Nearby Galaxy Survey (NGS)
for ingestion into the JSA.
"""

planeURI_cache = {}
 
def rewrite_fits(insdf, outfits, project_name, dprcinst, workdir, tap, log):
    """
    Rewrite a single sdf file into FITS, setting custom headers as needed.
    
    The NGS did not originally preserve the membership headers OBSCNT, OBSn
    when converting from sdf to fits, but the headers are still present in 
    the sdf files.  This code therefore does the conversion again, using the
    same conversion that is done by the jsawrapdr pipeline script.  It calls
    the Starlink programs ndfcopy, fitsmod, and ndf2fits, so Starlink 
    must be installed and the kappa and convert packages configured.
    
    Arguments:
    insdf: input file, which will be in Starlink NDF format
    outfits: output file, which will be in FITS format
    project_name: used to label things, esp. in observationID
    dprcinst: data processing recipe instance identifier
    workdir: absolute path to the working directory (cwd by default)
    tap: a tools4caom2.tapclient object
    log: a tools4caom2.logger object
    """
    global planeURI_cache
    
    log.console('PROGRESS: ' + insdf)
    
    # Making a copy of the sdf file updates the provenance structure, and avoids 
    # accidentally corrupting the original file.  Beware of complications 
    # if the files are very large.
    myindir, mysdfile = os.path.split(insdf)
    sdfcopy = os.path.join(workdir, 'copy_' + mysdfile)
    
    mydir, myfile = os.path.split(outfits)
    
    # Find some useful Starlink commands
    if 'KAPPA_DIR' not in os.environ:
        log.console(' run kappa command before proceeding', logging.ERROR)
    if 'CONVERT_DIR' not in os.environ:
        log.console(' run convert command before proceeding', logging.ERROR)

    ndfcopy = os.path.abspath(
                    os.path.expandvars('$KAPPA_DIR/ndfcopy'))
    fitsmod = os.path.abspath(
                    os.path.expandvars('$KAPPA_DIR/fitsmod'))
    ndf2fits = os.path.abspath(
                    os.path.expandvars('$CONVERT_DIR/ndf2fits'))
    
    # ndfcopy will update the PROVENANCE structure to avoid needless repetition
    ndfcopy_cmd = [ndfcopy, insdf, sdfcopy]
    log.file(' '.join(ndfcopy_cmd))
    try:
        output = subprocess.check_output(ndfcopy_cmd,
                                         stderr=subprocess.STDOUT)
    except subprocess.CalledProcessError as e:
        log.console('ndfcopy command failed: \n' + output,
                    logging.ERROR)
    
    # We need to add the PRODUCT header so that ndf2fits will write the 
    # membership and provenance headers. 
    # Extract the product from the insdf file name.
    # The insdf filename looks like one of
    # NGC0210_nearest_v2-0.sdf
    # NGC0210_nearest_totint20_v2-0.sdf
    # NGC0210_nearest_totint20_noise_v2-0.sdf
    # i.e. 2, 3 or 4 tokens followed by a version number.  The product is
    # "reduced" if there are two tokens, the third token if there are three, 
    # and the dash-separated concatenation of the third and fourth if there are 
    # four tokens. 
    name_token = mysdfile.split('_')
    if len(name_token) == 3:
        product = 'reduced'
    elif len(name_token) == 4:
        product = name_token[2]
    elif len(name_token) == 5:
        product = name_token[2] + '-' + name_token[3]
    else:
        log.console('name_token = ' + repr(name_token) + 
                    ' does not have 3-5 tokens',
                    logging.ERROR)

    # The set of science_products guides how we will sort files into planes.
    # NGS had originally intended two planes, for totint20 and reduced+others.
    # However totint20 dors not have energy metadata, so we will start with just
    # one plane.
    science_product = 'reduced'
    # if product == 'totint20':
    #     science_product = product

    # fitswrite will add the product header that is needed for the provenance 
    # to be written.  
    fitsmod_cmd = [fitsmod,
                   'edit=write',
                   'mode=interface',
                   'position=!',
                   sdfcopy,
                   'product',
                   'value=' + product,
                   'comment="product"']
    log.file(' '.join(fitsmod_cmd))
    try:
        output = subprocess.check_output(fitsmod_cmd,
                                         stderr=subprocess.STDOUT)
    except subprocess.CalledProcessError as e:
        log.console('fitsmod command failed: \n' + output,
                    logging.ERROR)
    
    # Now convert the sdfcopy into outfits
    comp = 'd'
    if product == 'reduced':
        comp = 'dv'
    # Convert to a CADC-compliant FITS file
    ndf2fits_cmd = [ndf2fits,
                    sdfcopy,
                    outfits,
                    'provenance=cadc',
                    'proexts',
                    'profits',
                    'prohis',
                    'duplex',
                    'checksum',
                    'encoding="fits-wcs(cd)"',
                    'comp=' + comp]
    log.file(' '.join(ndf2fits_cmd))
    try:
        output = subprocess.check_output(ndf2fits_cmd,
                                         stderr=subprocess.STDOUT)
    except subprocess.CalledProcessError as e:
        log.console('ndf2fits command failed: \n' + output,
                    logging.ERROR)
    
    hdulist = pyfits.open(outfits, mode='update')
    head = hdulist[0].header
    
    # Gather all the modified headers into headerdict, then update them in head.
    # We could edit head directly, but this approach allows the code to be 
    # modified for other sources of header metadata more easily.
    headerdict = {}
    
    # NGS is one of the JCMT Legacy Surveys
    headerdict['INSTREAM'] = 'JCMTLS'
    
    instrument = head['INSTRUME']

    # Observations will be distinguished by source and instrument configuration
    if instrument == 'HARP':
        restfreqstrs = {'co3-2': '345796MHz'}
        restfreqnums = {'co3-2': 345795989900.0}
        transition = 'unknown'
        if 'MOLECULE' in head and 'TRANSITI' in head:
            transition = re.sub(r'[^-0-9A-Za-z]', 
                                '', 
                                head['MOLECULE'] + head['TRANSITI']).lower()
        # ASN_ID fills Observation.observationID
        headerdict['ASN_ID'] = '-'.join([project_name,
                                         re.sub(r'\s', '', 
                                                head['OBJECT']).lower(),
                                         transition])
        # Set the resfreq string used in the PRODID header.
        restfreqstr = 'unknownHz'
        restfreq = None
        if transition in restfreqstrs:
            restfreqstr = restfreqstrs[transition]
            restfreq = restfreqnums[transition]
        else:
            log.console('transition = ' + transition + ' is not in ' +
                        repr(restfreqstrs.keys()),
                        logging.ERROR)
        
        bwmode = head['BWMODE']
        headerdict['PRODID'] = '-'.join([science_product,
                                         restfreqstr,
                                         bwmode])
    elif instrument == 'SCUBA-2':
        headerdict['ASN_ID'] = '-'.join([project_name,
                                         re.sub(r'\s', '', 
                                                head['OBJECT']).lower(),
                                         'continuum'])
        filter = str(head['FILTER']) + 'um'
        headerdict['PRODID'] = '-'.join([science_product,
                                         filter]) 
    
    headerdict['ASN_TYPE'] = 'custom'
    # headerdict['MBRCNT'] = 0 # number of membership URIs
    # headerdict['MBR1'] = <membership URI 1>
    headerdict['OBS-TYPE'] = 'science'
    
    # If defined, head['PROJECT'] is the observing project, which will be 
    # by the publication project for JLS products but can be recovered through
    # the membership.
    headerdict['PROJECT'] = project_name
    headerdict['SURVEY'] = 'NGS'

    # The PI of the whole project, not the PI of the project for which raw data
    # was collected.
    headerdict['PI'] = 'Christine Wilson'
    
    # Ambiguous, it may be that no document has this title.  Check with NGS.
    headerdict['TITLE'] = 'Nearby Galaxies Legacy Survey'
    
    # We use the instrument name to set the observationID
    # headerdict['INSTRUME'] already set correctly
    # headerdict['INBEAM'] already set correctly
    # headerdict['BACKEND'] already set correctly
    # headerdict['SW_MODE'] already set correctly
    # headerdict['SCAN_PAT'] already set correctly
    # headerdict['OBS_SB'] already set correctly
    # headerdict['SB_MODE'] already set correctly
    
    # Oddly, not set by default
    headerdict['TELESCOP'] = 'JCMT'
    headerdict['OBSGEO-X'] = -5464588.652191697
    headerdict['OBSGEO-Y'] = -2493003.0215722183
    headerdict['OBSGEO-Z'] = 2150655.6609171447
    
    # headerdict['OBJECT'] already set correctly
    headerdict['TARGTYPE'] = 'OBJECT' # or 'FIELD'
    # headerdict['ZSOURCE'] = <redshift in BARYCENT frame>
    # headerdict['TARGKEYW'] = <target keyword string>
    headerdict['MOVING'] = False
    # headerdict['OBSRA'] = <target RA in ICRS>
    # headerdict['OBSDEC'] = <target Dec in ICRS>
    # headerdict['RADESYS'] = <RA/Dec system>
    # headerdict['EQUINOX'] = <equinox of coordinates>
    # headerdict['FILTER'] = <characteristic wavelength>
    
    # if instrument == 'HARP':
    #     headerdict['RESTFREQ'] = restfreq
    
    # headerdict['BWMODE'] = already set correctly
    # headerdict['SUBSYSNR'] = ignored
    # headerdict['RECIPE'] = already set correctly
    # headerdict['PROCVERS'] = <data processing software version>
    # headerdict['ENGVERS'] = <data processing engine version>
    
    # DataProductType is a crude classification of the shape of the data
    if product in ['reduced', '20kms']:
        headerdict['DATAPROD'] = 'cube'
    else:
        headerdict['DATAPROD'] = 'image'

    # ProductType is a crude classification of the nature of the data 
    # in each extension of a FITS file
    if product == science_product:
        headerdict['PRODTYPE'] = '0=science,auxiliary'
    else:
        headerdict['PRODTYPE'] = 'auxiliary'

    # CalibrationLevel is a crude classification of the degree of processing 
    headerdict['CALLEVEL'] = 'calibrated'

    # Ask who gets the credit
    headerdict['PRODUCER'] = 'NGS'
    
    # Fill the data processing recipe instance identifier, filled in this project
    # with the project name, galaxy class and galaxy ID.
    headerdict['DPRCINST'] = dprcinst
    
    # A suitable, representative processing datetime
    headerdict['DPDATE'] = '2010-07-07T00:00:00'
    
    # The membership headers OBSCNT, OBS1, ... are set correctly, but the 
    # provenance headers PRVCNT, PRV1, ... are mostly set to temporary
    # files that do not exist.  Delete them from the header and insert
    # input headers INPCNT, INP1, ... derived from a TAP query on the assumption
    # that all of the input planes have productID like 'raw-%'.  Complain if
    # any inputs do not match this pattern.
    if 'PRVCNT' in head and int(head['PRVCNT']) > 0:
        for n in range(int(head['PRVCNT'])):
            prvn = 'PRV' + str(n + 1)
            del head[prvn]
        del head['PRVCNT']
    
    if 'OBSCNT' in head and int(head['OBSCNT']) > 0:
        inpcnt = 0
        for n in range(int(head['OBSCNT'])):
            obsn = 'OBS' + str(n + 1)
            # Construct planeURI's for the input planes using a TAP query
            # to find all the raw planes in the raw observations.
            # The OBSn headers actually record the obsid_subsysnr values from
            # the ACSIS, SCUBA2 and FILES tables, but the latter will not be 
            # accessible at sites other than the JAC.  Formally there is no ICD
            # that allows us to convert obsid_subsysnr into obsid except to
            # look up the value in the FILES table, but as a practical 
            # alternative, we can split the obsid_subsysnr into parts
            #     instrument, obsnum, dateobs, subsysnr = \
            #         head[obsn].split('_')
            # where the dateobs uniquely identifies the observation.
            # We can then use a TAP query to find the actual observationID and 
            # productID for each raw plane.
            if obsn in head:
                raw_instr = None
                raw_obsnum = None
                raw_dateobs = None
                raw_subsysnr = None
                obsnval = head[obsn]
                
                # Have we already found this planeURI?
                if obsnval in planeURI_cache:
                    inpcnt += 1
                    inpnn = 'INP' + str(inpcnt)
                    headerdict[inpnn] = planeURI_cache[obsnval]
                else:
                    try:
                        raw_instr, raw_obsnum, raw_dateobs, raw_subsysnr = \
                            head[obsn].split('_')
                    except:
                        pass
                    if not (raw_instr is None or raw_dateobs is None):
                        obsn_pat = "'" + raw_instr + "%" + raw_dateobs + "'"
                        
                        tapcmd = '\n'.join([
                            "SELECT DISTINCT",
                            "       Observation.observationID,",
                            "       Plane.productID",
                            "FROM caom2.Observation AS Observation",
                            "    INNER JOIN caom2.Plane AS Plane",
                            "        ON Observation.obsID = Plane.obsID",
                            "WHERE Observation.collection = 'JCMT'",
                            "    AND Observation.observationID LIKE " + obsn_pat,
                            "    AND Plane.productID LIKE 'raw%'"])
                        results = tap.query(tapcmd)
                        if results:
                            for raw_obsid, raw_prodid in results:
                                inpcnt += 1
                                planeURI = '/'.join(['caom:JCMT',
                                                     raw_obsid,
                                                     raw_prodid])
                                inpnn = 'INP' + str(inpcnt)
                                headerdict[inpnn] = planeURI
                                planeURI_cache[obsnval] = planeURI
                headerdict['INPCNT'] = inpcnt
    
    # Are there any new keywords in the headerdict
    newkeys = False
    for key in headerdict:
        if key not in head:
            newkeys = True
    
    # If so, add a comment to label the section containing new keys
    # Existing keywords will be updated in-place
    if newkeys:
        endcard = len(head)
        head.update('', '', comment='JSA Headers', after=endcard)
    
    # update FITS headers with those supplied in headerdict
    for key in sorted(headerdict.keys(), reverse=True):
        if key in head:
            head.update(key, headerdict[key])
        else:
            head.update(key, headerdict[key], after=endcard)

    hdulist.flush()
    hdulist.close()
    os.remove(sdfcopy)
    # os.remove(fitscopy)

def fix_name(outdir, prefix, filename):
    """
    Compose a new name from the prefix and basename of the file.
    """
    dirpath, basename = os.path.split(filename)
    file_id = os.path.splitext(basename)[0].lower()
    return os.path.join(outdir, dirpath, prefix + '_' + file_id + '.fits')

def readfilelist(rootdir, subdirpath, filelist, log):
    """
    Construct a list of file names rooted at indir by reading names from indir
    and calling readfilelist recursively for each directory.
    """  
    dirlist = []
    if subdirpath:
        readdir = os.path.join(rootdir, subdirpath)
    else:
        readdir = rootdir
    
    for f in os.listdir(readdir):
        log.file('examine: ' + f)
        filename = os.path.join(rootdir, subdirpath, f)
        
        if os.path.isfile(filename):
            # Append the path to the file to filelist
            filelist.append(os.path.join(subdirpath, f))
        
        elif os.path.isdir(filename):
            # Append the subdirectory to the list awaitin recursion
            dirlist.append(os.path.join(subdirpath, f))
    for d in dirlist:
        # recurse into the subdirectory d
        readfilelist(rootdir, d, filelist, log)

def is_ingestible(filename):
    """
    Return True if the extension is an sdf file, False otherwise
    """
    return os.path.splitext(filename)[1] == '.sdf'

def run():
    """
    The run() method for ngs_prepare_files.  In addition to its function in
    preparing NGS files, it serves as an example of how other projects
    can prepare their own files.
    
    The code has been simplified in comparison to prepare_files.py by the
    elimination of command line and csv file methods for setting file headers.
    """
    progname = os.path.basename(os.path.splitext(sys.path[1])[0])

    ap = argparse.ArgumentParser(progname)

    ap.add_argument('--proxy',
                    default='~/.ssl/cadcproxy.pem',
                    help='path to CADC proxy')

    ap.add_argument('--indir',
                    required=True,
                    help='existing release directory')
    ap.add_argument('--outdir',
                    required=True,
                    help='new release directory to which files will be written')
    # default prefix is the same as the NGS project name
    ap.add_argument('--prefix',
                    required=True,
                    help='prefix for ingestible file names')

    ap.add_argument('--workdir',
                    default='.',
                    help='directory to hold working files (default=cwd)')
    
    ap.add_argument('--log',
                    default='ngs_prepare_files_' + utdate_string() + '.log',
                    help='(optional) name of log file')

    # verbosity
    ap.add_argument('--debug', '-d',
                    action='store_true',
                    help='run in debug mode')

    a = ap.parse_args()

    loglevel = logging.INFO
    if a.debug:
        loglevel = logging.DEBUG

    with logger(a.log, loglevel).record() as log:
        # Report all command line arguments
        log.file(progname)
        for attr in dir(a):
            if attr != 'id' and attr[0] != '_':
                log.file('%-15s= %s' % 
                                 (attr, str(getattr(a, attr))))
        log.console('log = ' + a.log)

        proxy = os.path.abspath(
                    os.path.expandvars(
                        os.path.expanduser(a.proxy)))
        
        tap = tapclient(log, proxy)
        
        workdir = os.path.abspath(
                    os.path.expandvars(
                        os.path.expanduser(a.workdir)))

        # Convert a.indir into an abspath and verify that it is a directory
        if not a.indir:
            log.console('specify --indir as the path to the input directory',
                        logging.ERROR)
        a.indir = os.path.abspath(
                    os.path.expandvars(
                        os.path.expanduser(a.indir)))
        if not os.path.isdir(a.indir):
            log.console('indir = ' + a.indir + 
                        ' is not a directory',
                        logging.ERROR)

        # Convert a.outdir into an abspath and verify that it is a directory
        if not a.outdir:
            log.console('specify both --indir and --outdir, since it is '
                        'forbidden to overwrite the original files',
                        logging.ERROR)
        a.outdir = os.path.abspath(
                   os.path.expandvars(
                       os.path.expanduser(a.outdir)))
        if not os.path.isdir(a.outdir):
            log.console('output directory ' + a.outdir + 
                        ' is not a directory',
                        logging.ERROR)
        
        # filelist contains a list of file paths relative to a.indir.
        filelist = []
        readfilelist(a.indir, '', filelist, log)
        
        for infile in filelist:
            # Be sure the directory path exists before creating the FITS file
            dirpath = os.path.join(a.outdir, 
                                    os.path.dirname(infile))
            if not os.path.isdir(dirpath):
                os.makedirs(dirpath)
            
            # Existing FITDS files are defective, so skip them
            if os.path.splitext(infile)[1] == '.fits':
                continue
    
            inpath = os.path.join(a.indir, infile)
            if is_ingestible(inpath):
                # Data files are always in a dirctory called Data in the NGS.
                # The galaxy class and object name are the preceding two
                # directories.
                dirparts = inpath.split('/')
                dprcinst = ''
                i = -1
                for part in dirparts:
                    i += 1
                    if part == 'Data':
                        break
                if i > 1:
                    dprcinst = '-'.join([a.prefix, 
                                         dirparts[i-2], 
                                         dirparts[i-1]])
                if not dprcinst:
                    log.console('could not form dprcinst from ' + 
                                repr(dirparts),
                                logging.ERROR) 
                                
                # Add the prefix to fits files generated from sdf files,
                # but not to other files that will simply be copied.
                newfile = fix_name(a.outdir, a.prefix, infile)
                
                rewrite_fits(inpath, 
                             newfile, 
                             a.prefix,
                             dprcinst,
                             workdir,
                             tap,
                             log)
            else:
                newfile = os.path.join(a.outdir, infile)
                shutil.copyfile(inpath, newfile)

if __name__ == '__main__':
    run()