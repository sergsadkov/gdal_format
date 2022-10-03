# gdal_format
Python TUI for performing operations with raster files using GDAL functions with minimum code writing. Might make Python scripts and QGIS plugins a little shorter.

It is inspired with my 3 years of working as a GIS Developer in Python when I faced a lot of routine operations similar to each other. I needed to write the same functions (with the same typos) over and over again. So I made several functions to make my life easier. Here I tried to put 'em all together hoping it would be useful for anybody else.

Python files:

  data_properties.py -- DataProperties class containing all features for a new raster file
  
  check_properties.py -- CheckProperties class containing mismatches for an existing raster file
  
  save_data.py -- functions to create new files:
      
      saveRaster -- save an existing raster with the new properties
      
      saveRGB -- save an RGB image, with or without alpha-channel
  
  temp_paths.py -- creates paths for temp files in a specified directory deleting old data from there
