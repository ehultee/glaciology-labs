## Helper functions to clean up GEOL 362 notebooks
## Based on Clubes de Ciencia helper funcs
## Original: 5 July 2019  EHU
## New: 24 Aug 2021
import xarray as xr
import pandas as pd
import numpy as np
from oggm import utils

def ice_to_freshwater(icevol, rho_ice=900, rho_water=1000):
    """Cleanly convert volume of glacial ice (km3) to equivalent volume fresh water (liter).
    Arguments:
        icevol = volume of ice to convert, in km3
        rho_ice = density of glacial ice (default 900 kg/m3)
        rho_water = density of freshwater (default 1000 kg/m3)
        """
    km3_to_ltr = 1E12
    water_vol_km3 = icevol * rho_ice / rho_water
    return water_vol_km3 * km3_to_ltr


def read_run_results(gdir, filesuffix=None):
    """Reads the output diagnostics of a simulation and puts the data in a pandas dataframe.
    
    Parameters
    ----------
    gdir : the glacier directory
    filesuffix : the file identifier 
    
    Returns
    -------
    a pandas Dataframe with monthly temp and precip
    """
    
    with xr.open_dataset(gdir.get_filepath('model_diagnostics', filesuffix=filesuffix)) as ds:
        ds = ds.load()
      
    # Lemgth needs filtering
    ts = ds.length_m.to_series()
    ts = ts.rolling(12*3).min()
    ts.iloc[0:12*3] = ts.iloc[12*3]
    
    # Volume change
    delta_vol = np.append(ds.volume_m3.data[1:] - ds.volume_m3.data[0:-1], [0])
    
    if ds.calendar_month[0] == 10 and gdir.cenlat < 0:
        # this is to cover up a bug in OGGM
        _, m = utils.hydrodate_to_calendardate(ds.hydro_year.data, ds.hydro_month.data, start_month=4)
        ds.calendar_month[:] = m
    
    odf = pd.DataFrame()
    odf['length_m'] = ts
    odf['volume_m3'] = ds.volume_m3
    odf['delta_water_m3'] = delta_vol * 0.9
    odf['month'] = ds.calendar_month
    
    return odf


def read_climate_statistics(gdir):
    """Reads the annual cycle of climate for [1985-2015] at the glacier terminus elevation.
    
    Parameters
    ----------
    gdir : the glacier directory
    
    Returns
    -------
    a pandas Dataframe with monthly average temp and precip
    """
    
    with xr.open_dataset(gdir.get_filepath('climate_monthly')) as ds:
        ds = ds.load()
        
    ds = ds.sel(time=slice('1985', '2015'))
        
    dsm = ds.groupby('time.month').mean(dim='time')
    odf = pd.DataFrame()
    odf['temp_celcius'] = dsm.temp.to_series()
    odf['prcp_mm_mth'] = dsm.prcp.to_series()
    
    # We correct for altitude difference
    d = utils.glacier_statistics(gdir)
    odf['temp_celcius'] += (ds.ref_hgt - d['flowline_min_elev']) * 0.0065
    
    return odf

def plot_xz_bed(x, bed, ax=None, ylim=None):
    """This function implements a glacier bed, prepared axes and a legend in
    altitude vs. distance along a glacier plot.  Based on function of the same 
    name in OGGM-Edu, but adds explicit axes argument.
    Parameters
    ----------
    x : ndarray
        distance along glacier (all steps in km)
    bed : ndarray
        bed rock
        
    Parameters (Optional)
    ----------
    ax : matplotlib axes instance on which to plot
        If None, calls plt.gca()
    ylim : tuple, y-limits of plot
        If None, calls ax.get_ylim()
    """
    if ax is None:
        ax = plt.gca()
    if ylim is None:
        ylim = ax.get_ylim()
       
    ax.plot(x, bed, color='k', label='Bedrock', linestyle=':', linewidth=1.5)
    ax.set_xlabel('Distance along glacier [km]')
    ax.set_ylabel('Altitude [m]')
    ax.set_ylim(ylim)
    ax.legend(loc='best', frameon=False)
