# README

# rSTEMMUS_SCOPE
Codes to run the STEMMUS_SCOPE model in MATLAB from R

### To install rSCOPE, use:
devtools::install_github("AlbyDR/rSTEMMUS.SCOPE")

library(rSTEMMUS.SCOPE)

It may need to install the rhdf5 package before using BiocManager

install.packages("BiocManager")

BiocManager::install("rhdf5")

### Time Series inputs (vectors - .dat)

| symbol | variable                     |    unit     | observation         |
|:-------|:-----------------------------|:-----------:|:--------------------|
| rain\_ | Precipitation                | \[cm s-1\]  | 1 cm/s = 36000 mm/h |
| Ta\_   | Air Temperature              |   \[°C\]    |                     |
| RH\_   | Relative Humidity            |    \[%\]    | 0 \<= RH \<= 100    |
| p\_    | Atmospheric Pressure         |   \[hPa\]   |                     |
| u\_    | Wind speed                   |  \[m s-1\]  | u \>= 0.05          |
| CO2\_  | Carbon Dioxide Concentration | \[mg m-3\]  |                     |
| Rin\_  | Incoming Shortwave Radiation |  \[W m-2\]  |                     |
| Rli\_  | Incoming Longwave Radiation  |  \[W m-2\]  |                     |
| LAI\_  | Leaf Area Index              | \[m2 m-2\]  | LAI \>= 0.01        |
| ea\_   | Air Vapor Pressure           |   \[hPa\]   | eq-01               |
| VPD\_  | Vapor Pressure Deficit       |   \[hPa\]   | eq-02               |
| tts\_  | Zenith Solar Angle           |    \[-\]    |                     |
| t\_    | Timestamp (doy_float)        |    \[-\]    |                     |
| year\_ | Year                         | \[integer\] | decimal Julian day  |
| Cab\_  | Chlorophyll ab               | \[ug cm-2\] | Calendar year       |
| hc\_   | Canopy Height                |    \[m\]    | hc \>= 0.01         |

#### $$
ea = 6.107*10^{7.5 X Ta \choose 237.3 + Ta} * {RH \choose 100}
$$

#### $$
VPD = 6.107*10^{7.5 X Ta \choose 237.3 + Ta} * 1 - {RH\choose 100}
$$

``` r
1 + 1
```

    [1] 2

## Input required per location

**Location and Vegetation characteristics**

| symbol | variable | unit | observation |
|:---|:---|:--:|:---|
| sitename | Name on the site | string | 2 chr - 3 chr (DE-C01) |
| Dur_tot | Number of timestamps | double |  |
| DELT | Timestep size in seconds | double | 60x30min or 60x60min hourly |
| latitude | Latitude (x) | degree |  |
| longitude | Longitude (y) | degree |  |
| elevation | Altitude (DEM) | \[m\] |  |
| reference_height | Measurement Height (z) | \[m\] |  |
| IGBP_veg_long | Long name IGBP vegetation class | string |  |
| canopy_height | Canopy Height (hc) | \[m\] |  |

``` r
1 + 2
```

    [1] 3

**Soil Initial Conditions** (soil_init.mat)

| symbol | variable                   |    unit    | observation |
|:-------|:---------------------------|:----------:|:------------|
| SWC    | Initial soil water content | \[m3 m-3\] | 1\)         |
| Ts     | Initial soil temperature   |   \[°C\]   | 2\)         |

note: data from CDS ERA5 Land:  
1) “volumetric soil water” layer 1 to 4 (swvl1,swvl2,swvl3,swvl4), and  
2) “Skin temperature” (skt) and “Soil temperature” level 1 to 4
(stl1,stl2,stl3,stl4).

``` r
 1 + 3
```

    [1] 4

**Soil Properties** (soil_parameters.mat)

<table style="width:98%;">
<colgroup>
<col style="width: 16%" />
<col style="width: 30%" />
<col style="width: 16%" />
<col style="width: 34%" />
</colgroup>
<thead>
<tr class="header">
<th style="text-align: left;">symbol</th>
<th style="text-align: left;">variable</th>
<th style="text-align: center;">unit</th>
<th style="text-align: left;">observation</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td style="text-align: left;">SaturatedK</td>
<td style="text-align: left;">Saturated hydraulic conductivity</td>
<td style="text-align: center;">[cm s-1]</td>
<td style="text-align: left;"><p>1x6 <sup><span
class="smallcaps">1</span></sup> [Ks]/(24*3600)</p>
<p>cm d-1</p></td>
</tr>
<tr class="even">
<td style="text-align: left;">ks0</td>
<td style="text-align: left;">First element of SaturatedK vector</td>
<td style="text-align: center;">[cm s-1]</td>
<td style="text-align: left;">1x1 <sup><span
class="smallcaps">1</span></sup></td>
</tr>
<tr class="odd">
<td style="text-align: left;">porosity</td>
<td style="text-align: left;">Porosity</td>
<td style="text-align: center;">[m3 m-3]</td>
<td style="text-align: left;">1x6 <sup><span
class="smallcaps">1</span></sup> [thetas]</td>
</tr>
<tr class="even">
<td style="text-align: left;">theta_s0</td>
<td style="text-align: left;">First element of porosity vector</td>
<td style="text-align: center;">[m3 m-3]</td>
<td style="text-align: left;">1x1 <sup><span
class="smallcaps">1</span></sup></td>
</tr>
<tr class="odd">
<td style="text-align: left;">SaturatedMC</td>
<td style="text-align: left;">Saturated SWC</td>
<td style="text-align: center;">[m3 m-3]</td>
<td style="text-align: left;">1x6<sup><span
class="smallcaps">1</span></sup> [thetas]</td>
</tr>
<tr class="even">
<td style="text-align: left;">ResidualMC</td>
<td style="text-align: left;">Residual SWC</td>
<td style="text-align: center;">[m3 m-3]</td>
<td style="text-align: left;">1x6 <sup><span
class="smallcaps">1</span></sup> [thetar]</td>
</tr>
<tr class="odd">
<td style="text-align: left;">Coefficient_Alpha</td>
<td style="text-align: left;">Coefficient Alpha</td>
<td style="text-align: center;">?</td>
<td style="text-align: left;">1x6 <sup><span
class="smallcaps">1</span></sup> [alpha]</td>
</tr>
<tr class="even">
<td style="text-align: left;">Coefficient_n</td>
<td style="text-align: left;">Coefficient n</td>
<td style="text-align: center;">[-]</td>
<td style="text-align: left;">1x6 <sup><span
class="smallcaps">1</span></sup> [n]</td>
</tr>
<tr class="odd">
<td style="text-align: left;">fieldMC</td>
<td style="text-align: left;">Field Capacity</td>
<td style="text-align: center;">[m3 m-3]</td>
<td style="text-align: left;">6x1 <sup><span
class="smallcaps">2</span></sup></td>
</tr>
<tr class="even">
<td style="text-align: left;">FOS</td>
<td style="text-align: left;">Sandy Fraction</td>
<td style="text-align: center;">fraction (/100)</td>
<td style="text-align: left;">6x1x1<sup>3</sup> [SAND1/2]</td>
</tr>
<tr class="odd">
<td style="text-align: left;">FOC</td>
<td style="text-align: left;">Clay</td>
<td style="text-align: center;">fraction (/100)</td>
<td style="text-align: left;">6x1x1<sup>3</sup> [CLAY1/2]</td>
</tr>
<tr class="even">
<td style="text-align: left;">MSOC</td>
<td style="text-align: left;">Organic Fraction (Carbon)</td>
<td style="text-align: center;">fraction (/10000)</td>
<td style="text-align: left;">6x1x1<sup>3</sup> [OC1/2]</td>
</tr>
<tr class="odd">
<td style="text-align: left;">fmax</td>
<td style="text-align: left;"></td>
<td style="text-align: center;">?</td>
<td style="text-align: left;">1x1 surfdata</td>
</tr>
<tr class="even">
<td style="text-align: left;">Coef_Lamda</td>
<td style="text-align: left;">Lambda per depth</td>
<td style="text-align: center;">?</td>
<td style="text-align: left;">6x1x1<sup>4</sup> Lambda folder</td>
</tr>
</tbody>
</table>

<sup><span class="smallcaps">1</span></sup> Derived from the
**PTF_SoilGrids_Schaap** datasets (*n, alpha, Ks, thetas, thetar) from
the valid depths: 0, 5, 15, 30, 60, 100 and 200 cm (sl1* to sl7),
excluding *sl3 (15cm).*

<sup><span class="smallcaps">2</span></sup> Calculated field capacity
form field moisture content  
field_moisture_content = theta_r + (theta_s - theta_r) / (1 + (alpha \*
phi_fc) ^ coef_n) ^ (1 - 1/coef_n)  
phi_fc = 341.9 - soil water potential at field) capacity (cm)

<sup>3</sup> Both layers (1, 2) are combined and the values from the
depth 1,3,5,6,7,8 are used per variable

<sup>4</sup> Only the Lambda file layers l1, l3, l5, l6, l7, l8
(depth_indices) were combined and used

``` r
1 + 4
```

    [1] 5
