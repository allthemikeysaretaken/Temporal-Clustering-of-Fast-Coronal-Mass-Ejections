from active_region_timelapse import *

## SPECIFY DIRECTORIES EXPLICITLY
dir_home = '/Users/.../'
dir_figs = '{}Figures/'.format(dir_home)
dir_data = '{}Data/'.format(dir_home)
dir_fits = '{}MichelsonDopplerImager/'.format(dir_data)
dir_active_regions = '{}ActiveRegions/'.format(dir_figs)
##

## GENERATE TIME-LAPSE OF ACTIVE REGIONS
fps = 5
img_extension = '.png'
dpi = 800
mov_extension = '.mkv'
for cmap in (None, 'hot'):
    MDI = MichelsonDopplerImager(rdirectory=dir_fits, sdirectory=dir_active_regions)
    MDI.convert_from_fits(cmap, img_extension, dpi)
    MDI.save_timelapse(fps, cmap, mov_extension, codec='mpeg4', search_extension=img_extension)
##
