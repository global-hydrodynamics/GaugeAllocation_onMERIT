Allocate GRDC gauges on 1min MERIT Hydro.

---------------------------------------
 Dai Yamazaki (U-Tokyo)
 2022 April 5th
 Distributed under CC-BY 4.0 license.
---------------------------------------

[1] Code:
src/allocate_GRDC.f90

Calculate optimum gauge location on MERIT Hydro 1min river map, for each gauge (ID, lat, lon, area)

[2] Input Datas

<2a: 1min river map>
- Please use 1min river map of CaMa-Flood. A sample map is included in the package (1min_v400)
- 1min map should be linked as 1min_river (automatically done in s01-GRDC.sh script)

<2b: list of gauges>
GRDC_input.txt
- List of gauges with data (ID, Lat, Lon, Area [km2])

[3] Output:
GRDC_alloc.txt
- List of gauges with allocated data (ID, Lat, Lon, Area [km2])

GaugeID  basinID  lat_ori  lon_ori  area_ori  lat_MERIT  lon_MERIT  area_MERIT  diff  error  ix  iy  bsn_lat  bsn_lon

-- "basin ID" is river basin ID of CaMa-Flood map (basin.bin)
-- "lat_ori  lon_ori  area_ori" are valuse in input data
-- "lat_MERIT  lon_MERIT  area_MERIT" valus of allocated pixel of 1min river map
-- "diff error" are area difference and relative error
-- "ix iy" are ix,iy coordinate of allocated pixel (fortran id 1:nx, 1:ny)
-- "bsn_lat  bsn_lon" are center lat/lon of the allocated river basin

[4] How to use:
The code automatically finds the optimam pixel to allocate gauge.

However, gauge is sometimes not appropriatelly allocated due to
 -- errors in the gauge metadata (lat, lon, area)
 -- errors in the MERIT Hydro river network.

 Pleaes check the output data carefully (the developepr did this process on the Excel file included)

If large error is found: please modify the lat,lon,area data of the input data (or exclude the gauge from the list)




****** 

In addition to river gauge allocation, similar code for allocating virtual stations of Satellite Altimetry is included.

******