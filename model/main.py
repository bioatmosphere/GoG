import os
import sys
import subprocess

import numpy as np


# Now import the compiled Fortran modules
os.chdir("vegetation/src/build/")

import constants

print(constants.codename)



# import Parameters
# import Input
# import Output
# import Species
# import Site
# import Sitelist
# import GenusGroups
# import Model



# # Prepare input files
# initialize_inputFiles(filelist)

# # Prepare site and species data
# call initialize_sitelist(sites,species_present)

# # Prepare output files
# call initialize_outputFiles(species_present)

# # write runtime vars to screen
# call drawBanner(size(sites),species_present)

# # Start timing
# call cpu_time(start_time)

# for sndx in range(numsites):

#     # make sure the site exists
#     if ( sites(sndx)%site_wmo == rnvalid ) then
#         write(*,*) '             No site or climate file for site ',         &
#                                 sites(sndx)%site_id
#         write(*,*) '             Skipping site ',sites(sndx)%site_name
#         write(*,*) 
#         cycle
#     endif

#     # load climate and site specific vars,then adjust for altitude if requested
#     call set_site_rng_seed(fixed_seed)
#     call set_site_climate(same_climate,fixed_seed)

#     # skip this site if no species are present
#     if ( size(sites(sndx)%species) .eq. 0 )  then
#         write(*,*) '              No species present in site ',              &
#                                 sites(sndx)%site_id
#         write(*,*) '              Skipping site ',sites(sndx)%site_name
#         write(*,*) '             '
#         cycle
#     endif
#     call showProgress(sites(sndx))

#     # run the model
#     for year in range(numyears):
    
#         call BioGeoClimate(sites(sndx),year)
#         # We will write soil/CN/clim data for last year's trees but after BioGeo
#         call write_soil_data(sites(sndx),year)
#         call write_site_data(sites(sndx),year)
#         # Trees !
#         call Canopy(sites(sndx))
#         call Growth(sites(sndx))
#         call Mortality(sites(sndx))
#         call Renewal(sites(sndx))

#         if spinup == True:
#             if ( year < spinup_yrs ) then
#                 print_interval=10*year_print_interval
#             else
#                 print_interval=year_print_interval
#             endif
#         else
#             print_interval=year_print_interval
        

#         if ((mod(year,print_interval).eq.0) .or. year==numyears) then
#             call write_genus_data(sites(sndx),species_present,year)
#             call write_species_data(sites(sndx),species_present,year)
#             if (tree_level_data):
#                 call write_tree_data(sites(sndx),year)
        

#     call cpu_time(total_time)
#     write(*,*) 'Cumulative time :',total_time-start_time
#     write(*,'(a80)') &
#     '=============================================================================='