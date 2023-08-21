# Paris50C
Exploring CMIP6 to find episodes of tmax exceeding 50°C in the Paris area
Exploring episodes of extreme heat in the CMIP6 archive

by Pascal Yiou (LSCE), August 2023.

This suite of files contains code (bash and R) to extract extreme events (when a variable exceeds a threshold over a given region) in the CMIP6 archive.
The code is used in a paper submitted to Climate Services, by Yiou et al. (2023). It is based on a report (in French) of the GREC Francilien, on the possibility of exceeding 50°C in the Paris area: https://cdn.paris.fr/paris/2023/06/27/scenarios_paris50c_note-grec_2023-guwL.pdf

The codes have been developed by Pascal Yiou (LSCE, IPSL, U Paris Saclay). The codes are distributed "as is" under a CeCILL license:
http://www.cecill.info/
They can be downloaded and used for free for academic purposes.
For commercial uses, please contact Pascal Yiou (pascal.yiou@lsce.ipsl.fr)

The archive also comes with ERA5 and EOBS time series that were extracted with
the Climate Explorer: https://climexp.knmi.nl/

List of code files, in use order:
1. extract_var_CMIP6_IdF_EN.sh: bash code to extract a climate variable, e.g. tasmax, tasmin, tas from the CMIP6 archive, for all scenarios, all models, all runs. The code extracts a predefined region (in longitude and latitude). Many parameters of extraction can be changed (variable, region). The region to be extracted needs to be coherent with the extracted data from ERA5 and EOBS. This code is devised to run in BATCH mode, on the IPSL computing server (spiritx). It uses the file architecture of CMIP6.
2. GWD_var_CMIP6_EN.sh: bash code to extract yearly Global Surface Temprature (GST) from all simulations in the CMIP6 archive. This code is devised to run in BATCH mode, on the IPSL computing server (spiritx).  It uses the file architecture of CMIP6.
3. T_kstest_CMIP6_EN.R: R code to determine which CMIP6 model yield a similar probability distribution of summer tmax, to ERA5 and EOBS, for the Paris area (Ile de France).
4. Textr-CMIP6_EN.R: R code to identify models that have runs with episodes of temperature exceeding T0 (T0=48°C) in the Paris area.
5. Textr-CMIP6_diags_EN.R: R code to produce simple diagnostics of extremes (T>T0=48°C) in CMIP6 simulations.

List of data files:
1. era5_t2m_year_GWD.nc: global mean surface temperature from ERA5 (extracted from the Climate Explorer)
2. iensembles_025_tx_1.75-3.25E_48.25-49.25N_n_max.nc: tmax over the Ile de France region from EOBS (extracted from the Climate Explorer)
3. iera5_tmax_daily_eu_1.75-3.25E_48.25-49.25N_n_max.nc: tmax over the Ile de France region from ERA5 (extracted from the Climate Explorer)
4. GRID_WEUR: an ASCII file that describes a ERA5 grid over Western Europe.
