import requests
from tqdm import tqdm

def get_file_from_url(url, filename):

    r = requests.get(url, stream=True)
    
    file_size = int(r.headers.get('content-length', 0))
    initial_pos = 0

    with open(filename, 'wb') as f:
        with tqdm(total=file_size, unit='B',
                  unit_scale=True, unit_divisor=1024,
                  desc=filename, initial=initial_pos,
                  ascii=True, miniters=1) as pbar:
            for chunk in r.iter_content(32 * 1024):
                f.write(chunk)
                pbar.update(len(chunk))

#     urllib.request.urlretrieve(url, filename, reporthook=progress)
#     print('\n')
#
# def progress(count, blockSize, totalSize):
#   percent = int(count*blockSize*100/totalSize)
#   sys.stdout.write("\r" + "Progress" + "...%d%%" % percent)
#   sys.stdout.flush()
