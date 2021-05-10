### calling method

- install pycbc
- install lensGW
- install pycbc-lensGW

```
from pycbc import waveform
from lgw import *
waveform.add_custom_waveform('lensed', lensed_gw_td, 'time', force=True)
hp_tilde_lensed, hc_tilde_lensed = waveform.get_td_waveform(
                approximant="lensed", y0=0.1, y1=0.7937005, zS=2.0, zL=0.5, mL=1e8,
                ml=1e5, l0=0.5, l1=0, lens_model_list=['POINT_MASS'],
                optim='True', mass1=30, mass2=30, delta_t=1.0/16384, f_lower=20)
```
