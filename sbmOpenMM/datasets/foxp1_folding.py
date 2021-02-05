import os
import shutil
from .downloader import get_file_from_url

def download_dataset(output_dir, dcd_only=False, data_only=False, overwrite=False):
    """
    Download the dataset in the specified location.

    Parameters
    ==========
    output_dir : str
        Path to the output directory where to store the files.
    dcd_only : bool
        Get only DCD trajectories files?
    data_only : bool
        Get only energy DATA files?
    overwrite : bool
        Delete and replace all dataset files?
    """

    dcd = True
    data = True
    if dcd_only:
        if data_only:
            raise ValueError('Only one option must be given (dcd_only of data_only?)')
        data = False
        dcd = True

    if data_only:
        if dcd_only:
            raise ValueError('Only one option must be given (dcd_only of data_only?)')
        data = True
        dcd = False

    # Create output directory if it does not exists
    if os.path.exists(output_dir):
        if overwrite:
            shutil.rmtree(output_dir)
            os.mkdir(output_dir)
    else:
        os.mkdir(output_dir)

    if data:
        data_paths = pathsToDATAfiles()
        for p in data_paths:
            get_file_from_url(data_paths[p], output_dir+'/'+p)

    if dcd:
        dcd_paths = pathsToDCDfiles()
        for p in dcd_paths:
            get_file_from_url(dcd_paths[p], output_dir+'/'+p)

def pathsToDCDfiles():
    """
    Method to contain paths to DCD files in the dataset
    """
    doi = 'https://dataverse.csuc.cat/api/access/datafile/:persistentId?persistentId=doi:10.34810/data31/'
    # Get DCD files permanent links
    dcd_suffix = '_trajectory.dcd'
    dcd_paths = { '01':'32','02':'33','03':'37','04':'36','05':'35','06':'39',
                  '07':'38','08':'40','09':'41','10':'42','11':'43','12':'44',
                  '13':'45','14':'60','15':'46'}
    return {p+dcd_suffix:doi+dcd_paths[p] for p in dcd_paths}

def pathsToDATAfiles():
    """
    Method to contain paths to DATA files in the dataset
    """
    doi = 'https://dataverse.csuc.cat/api/access/datafile/:persistentId?persistentId=doi:10.34810/data31/'
    # Get DATA energy files permanent links
    data_suffix = '_energies.data'
    data_paths = { '01':'24','02':'16','03':'6','04':'18','05':'3','06':'11',
                   '07':'28','08':'20','09':'4','10':'7','11':'26','12':'1',
                   '13':'8','14':'14','15':'23'}
    return {p+data_suffix:doi+data_paths[p] for p in data_paths}
