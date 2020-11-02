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
fc10_19 <- load_as_vrt(file.path(data_folder, 'Degradation_Paper', 'GEE_Rasters'), 'fc10_19_ha[-.0-9]*tif')
NAvalue(fc10_19) <- -32768
fc <- stack(fc00_09, fc10_19)
b <- brick(fc, values=FALSE)  
names(b) <- c(paste0('fcc_0', seq(0, 9)),
              paste0('fcc_', seq(10, 19)))
b <- writeStart(b, filename="fc_change.tif", overwrite=TRUE,
                options="COMPRESS=LZW", datatype="INT2S")
bs <- raster::blockSize(b)
n_blocks <- bs$n
block_num <- 500
pb <- pbCreate(bs$n)
for (block_num in 1:n_blocks) {
    v <- getValuesBlock(fc, row=bs$row[block_num], nrows=bs$nrows[block_num])
    # fc_2000 is the first column of v, with other columns being subsequent 
    # years. So divide each column by that first column and then convert to 
    # percent change in forest cover
    v <- ((v  - v[, 1])/ v[, 1]) * 100
    # Set growth rate to maximum of other pixels for areas that grew from zero 
    # forest cover in 2000.
    isinfs <- is.infinite(v)
    v[isinfs] <- 0
    v[isinfs] <- max(v, na.rm=TRUE)
    b <- writeValues(b, v, bs$row[block_num])
    pbStep(pb, block_num)
}
b <- writeStop(b)
pbClose(pb)

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
