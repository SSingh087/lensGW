### calling method

#### install pycbc
```
import sys
!{sys.executable} -m pip install pycbc ligo-common emcee==2.2.1 --no-cache-dir
```
<hr>

#### install lenstronomy
  - upload version 1.7.0 on drive (unzip using `!unzip`) / upload unzipped version
  - mount drive
```
    from google.colab import drive
    drive.mount('/content/drive')
  ```
  - use `cd` comand to locate the file Eg:`cd drive/MyDrive/lenstronomy/lenstronomy-1.7.0`
  - `!python setup.py install`
  - return to home `cd ..`
 <HR>

#### install lensGW
  - fork lensGW 
  - clone lensGW `!git clone https://username:password@github.com/SSingh087/test.git`
  - `cd test`
  - `!python setup.py install`
  - return to home `cd ..`
  
<hr>

- install pycbc-lensGW
  - clone pycbc-lensGW `!git clone https://github.com/SSingh087/lensGW-PyCBC-plugin.git`
  - `!python setup.py install`
  - restart runtime

<hr>

### Run this sample for testing

```
from pycbc import waveform
from lgw import *
waveform.add_custom_waveform('lensed', lensed_gw_td, 'time', force=True)
hp_tilde_lensed, hc_tilde_lensed = waveform.get_td_waveform(
                approximant="lensed", y0=0.1, y1=0.7937005, zS=2.0, zL=0.5, mL=1e8,
                ml=1e5, l0=0.5, l1=0, lens_model_list=['POINT_MASS'],
                optim='True', mass1=30, mass2=30, delta_t=1.0/16384, f_lower=20)
```
