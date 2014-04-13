--This query produces the sample dataset for the POC
--It returns the lat, lon and gravity in the grid bounded
--by -35.100S, -35.900S, 148.6E & 149.4E at 0.02 deg resolution
--execution time on the dataset is 450ms on server as built.
SELECT lat, lon, (grav/1000)
FROM fine_data 
WHERE lat > -35900 
AND lat < -35100 
AND lon > 148600 
AND lon < 149400
AND lon % 20 = 1
AND lat % 20 =-1;

--The input parameters are shown below. Lat and Lon bounds
--are self explanatory. Returning the data at a lower resolution
--is achieved through the use of modulo division. The dataset is 
--at 0.002 deg resolution so setting the resolution parameter 
--will result in [RESOLUTION]/(1000*0.002) degree spacing between
--points. eg a value of 200 will give you 0.1 degree grid spacing.
--Query execution time is approximately inversely proportional to 
--the square of the resolution parameter for a given area. 
SELECT lat, lon, (grav/1000)
FROM fine_data 
WHERE lat > [LAT_MIN] 
AND lat < [LAT_MAX] 
AND lon > [LON_MIN] 
AND lon < [ON_MAX]
AND lon % [RESOLUTION] =  1
AND lat % [RESOLUTION] = -1;

