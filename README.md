## README

## rSTEMMUS_SCOPE

Codes to run the STEMMUS_SCOPE model in MATLAB from R


### To install rSTEMMUS_SCOPE, use:

``` r
devtools::install_github("AlbyDR/rSTEMMUS.SCOPE")
library(rSTEMMUS.SCOPE)
```

It may need to install the rhdf5 package before using BiocManager as below.

``` r
install.packages("BiocManager")
BiocManager::install("rhdf5")
```


MATLAB R2015b or superior is required to run SCOPE and the SCOPE code need to be downloaded and unzipped in a directory of your choice (suggestion "D:/model/STEMMUS_SCOPE/"). The SCOPE code is available at <https://github.com/EcoExtreML/STEMMUS_SCOPE>. Only src folder is needed.

<br/>

**The required data are divided into:**

```         
  1.1 Meteorological and vegetation properties time series inputs
  
  1.2 Site specific characteristics

  1.3 Soil Initial Conditions
  
  1.4 Soil Properties

  1.5 Constants and model settings
  
```

<br/>

#### 1.1 Time Series inputs (vectors .dat)
------------------------------------------------------------------------

| symbol | variable                     |   unit    | observation         |
|:-------|:-----------------------------|:---------:|:--------------------|
| rain\_ | Precipitation                | [cm s-1]  | 1 cm/s = 36000 mm/h |
| Ta\_   | Air Temperature              |   [°C]    |                     |
| RH\_   | Relative Humidity            |    [%]    | 0 \<= RH \<= 100    |
| p\_    | Atmospheric Pressure         |   [hPa]   |                     |
| u\_    | Wind speed                   |  [m s-1]  | u \>= 0.05          |
| CO2\_  | Carbon Dioxide Concentration | [mg m-3]  |                     |
| Rin\_  | Incoming Shortwave Radiation |  [W m-2]  |                     |
| Rli\_  | Incoming Longwave Radiation  |  [W m-2]  |                     |
| LAI\_  | Leaf Area Index              | [m2 m-2]  | LAI \>= 0.01        |
| ea\_   | Air Vapor Pressure           |   [hPa]   | eq-01               |
| VPD\_  | Vapor Pressure Deficit       |   [hPa]   | eq-02               |
| tts\_  | Zenith Solar Angle           |    [-]    |                     |
| t\_    | Timestamp (doy_float)        |    [-]    |                     |
| year\_ | Year                         | [integer] | decimal Julian day  |
| Cab\_  | Chlorophyll ab               | [ug cm-2] | Calendar year       |
| hc\_   | Canopy Height                |    [m]    | hc \>= 0.01         |

### $`ea = 6.107*10^{7.5 X Ta \choose 237.3 + Ta}* {RH \choose 100}`$ **(eq-01)** 

### $`VPD = 6.107*10^{7.5 X Ta \choose 237.3 + Ta}* 1 - {RH\choose 100}`$ **(eq-02)**

<br/>

#### 1.2 Site specific characteristics
------------------------------------------------------------------------

| symbol | variable | unit | observation |
|:---------------|:---------------------|:--------------:|:------------------|
| sitename | Name on the site | string | 2 chr - 3 chr (DE-C01) |
| Dur_tot | Number of timestamps | double |  |
| DELT | Timestep size in seconds | double | 60x30min or 60x60min hourly |
| latitude | Latitude (x) | degree |  |
| longitude | Longitude (y) | degree |  |
| elevation | Altitude (DEM) | [m] |  |
| reference_height | Measurement Height (z) | [m] |  |
| IGBP_veg_long | Long name IGBP vegetation class | string |  |
| canopy_height | Canopy Height (hc) | [m] |  |

<br/>

#### 1.3 Soil Initial Conditions (soil_init.mat)
------------------------------------------------------------------------

| symbol | variable | unit | observation |
|:-----------|:-----------|:----------:|:------------------------------------|
| SWC | Initial soil water content | m3 m-3 | “volumetric soil water” layer 1 to 4 (swvl1, swvl2, swvl3, swvl4) |
| Ts | Initial soil temperature | °C | "Skin temperature” (skt) and “Soil temperature” level 1 to 4 (stl1, stl2, stl3, stl4) |

note: data from CDS ERA5 Land

<br/>

#### 1.4 Soil Properties (soil_parameters.mat)
------------------------------------------------------------------------

| symbol | variable | unit | observation |
|:--------------|:-----------------|:-------------:|:----------------------|
| SaturatedK | Saturated hydraulic conductivity | [cm s-1] | 1x6 ^[1]^ [Ks]/(24\*3600) cm d-1 |
| ks0 | First element of SaturatedK vector | [cm s-1] | 1x1 ^[1]^ |
| porosity | Porosity | [m3 m-3] | 1x6 ^[1]^ [thetas] |
| theta_s0 | First element of porosity vector | [m3 m-3] | 1x1 ^[1]^ |
| SaturatedMC | Saturated SWC | [m3 m-3] | 1x6^[1]^ [thetas] |
| ResidualMC | Residual SWC | [m3 m-3] | 1x6 ^[1]^ [thetar] |
| Coefficient_Alpha | Coefficient Alpha | [cm-1] | 1x6 ^[1]^ [alpha] |
| Coefficient_n | Coefficient n | [-] | 1x6 ^[1]^ [n] |
| fieldMC | Field Capacity | [m3 m-3] | 6x1 ^[2]^ |
| FOS | Sandy Fraction | fraction (/100) | 6x1x1^3^ [SAND1/2] |
| FOC | Clay | fraction (/100) | 6x1x1^3^ [CLAY1/2] |
| MSOC | Organic Fraction (Carbon) | fraction (/10000) | 6x1x1^3^ [OC1/2] |
| fmax |  | [-] | 1x1 surfdata |
| Coef_Lamda | Lambda per depth | [-] | 6x1x1^4^ Lambda folder |

<sup>[1]</sup> Derived from the **PTF_SoilGrids_Schaap** datasets (*n, alpha, Ks, thetas, thetar) from the valid depths: 0, 5, 15, 30, 60, 100 and 200 cm (sl1* to sl7), excluding *sl3 (15cm).*

<sup>[2]</sup> Calculated field capacity form field moisture content\
field_moisture_content = theta_r + (theta_s - theta_r) / (1 + (alpha \* phi_fc) \^ coef_n) \^ (1 - 1/coef_n)\
phi_fc = 341.9 - soil water potential at field) capacity (cm)

<sup>3</sup> Both layers (1, 2) are combined and the values from the depth 1,3,5,6,7,8 are used per variable

<sup>4</sup> Only the Lambda file layers l1, l3, l5, l6, l7, l8 (depth_indices) were combined and used

<br/>

#### 1.5 Constants and model settings

Use the function of the family "info", "check" and "change" to get more information about which constants parameters and model setting from STEMMUS and SCOPE can be changed to calibrate the model for the site characteristics.
