library(gdalUtils)
library(raster)

data_folder <- 'D:/Documents and Settings/azvoleff/OneDrive - Conservation International Foundation/Data'

load_as_vrt <- function(folder, pattern, band=FALSE, raster=TRUE) {
    vrt_file <- tempfile(fileext='.vrt')
    files <- list.files(folder, pattern=pattern)
    if (length(files) == 0) {
        stop('No files found')
    }
    if (band) {
        gdalbuildvrt(paste0(folder, '/', files), vrt_file, b=band)
        r <- raster(vrt_file)
    } else {
        gdalbuildvrt(paste0(folder, '/', files), vrt_file)
        r <- raster::stack(vrt_file)
    }
    if (raster) {
        return(r)
    } else {
        return(vrt_file)
    }
}

population <- load_as_vrt(file.path(data_folder, 'Degradation_Paper', 'GEE_Rasters'), 'stack_pop_2000_2015_cnt[-.0-9]*tif')
names(population) <- c('pop_2000', 'pop_2005', 'pop_2010', 'pop_2015')
NAvalue(population) <- -32768
population_growth <- overlay(population$pop_2000, population$pop_2015,
                             fun=function(pop_2000, pop_2015) {
                                r = ((pop_2015/pop_2000)^(1/15)-1) * 100
                                # Set growth to zero when 2000 and 2015 
                                # populations are equal
                                r[pop_2000 == pop_2015] <- 0
                                # Set growth rate to maximum of other 
                                # pixels for areas that grew from zero 
                                # population in 2000.
                                isinfs <- is.infinite(r)
                                r[isinfs] <- 0
                                r[isinfs] <- max(r, na.rm=TRUE)
                                r
                             })
writeRaster(population_growth, filename='population_growth.tif', 
            overwrite=TRUE, options="COMPRESS=LZW", datatype="INT2S")

fc00_09 <- load_as_vrt(file.path(data_folder, 'Degradation_Paper', 'GEE_Rasters'), 'fc00_09_ha[-.0-9]*tif')
NAvalue(fc00_09) <- -32768
fc_2000 <- stack(fc00_09[[1]])
fc10_18 <- load_as_vrt(file.path(data_folder, 'Degradation_Paper', 'GEE_Rasters'), 'fc10_18_ha[-.0-9]*tif')
NAvalue(fc10_18) <- -32768
fc_2015 <- stack(fc10_18[[9]])

fc_change <- overlay(fc_2000, fc_2015,
    fun=function(fc_2000, fc_2015) {
        r = ((fc_2015 - fc_2000) / fc_2000) * 100
        # Set growth to zero when 2000 and 2015 covers are equal
        r[fc_2000 == fc_2015] <- 0
        # Set growth rate to maximum of other pixels for areas that grew from 
        # zero forest cover in 2000.
        isinfs <- is.infinite(r)
        r[isinfs] <- 0
        r[isinfs] <- max(r, na.rm=TRUE)
        r
    })
writeRaster(fc_change, filename='fc_change.tif', 
            overwrite=TRUE, options="COMPRESS=LZW", datatype="INT2S")

biomass <- load_as_vrt(file.path(data_folder, 'Degradation_Paper', 'GEE_Rasters'), 'biomass_above_below_tons[-.0-9]*tif')
NAvalue(biomass) <- -32768
names(biomass) <- c('agb', 'bgb')
overlay(biomass, fun=function(agb, bgb) {
        agb + bgb
    },
    filename='biomass.tif', 
    overwrite=TRUE,
    options="COMPRESS=LZW",
    datatype="INT2S")
