FC = gfortran

FFLAGS = -g -fdefault-real-8 -Wall -Wextra -Warray-temporaries -Wconversion -fimplicit-none -fbacktrace \
-ffree-line-length-0 -fcheck=all -ffpe-trap=zero,overflow,underflow -finit-real=nan

#FFLAGS = 

LDFLAGS = 

SRCS = mod_setup.f90 mod_const.f90 mod_surf1.f90 mod_soil.f90 master.f90 main.f90 solar_angle.f90 thermo_functions.f90 \
downward_radiation.f90  diffusion_coeff_louis.f90 diffusion_coeff_obrien.f90 init_surf1.f90 inv_mat_tri.f90 \
vertical_diffusion.f90 pbl_height.f90 init_surf1.f90 soil_def.f90 soil_prop.f90 rs_soil.f90 rs_veg.f90 \
drag_coeff_z0h.f90 fluxes.f90 energy_budget.f90 water_budget.f90  vdfppcfls.f90  \
lw_emissivities.f90 lw_radiation.f90 coriolis.f90 lw_radiation2.f90 lw_radiation3.f90 soilt_vertical_diffusion.f90 \
soilw_vertical_diffusion.f90 soil_prop_diff.f90

OBJS = ${SRCS:.f90=.o}

all: main

main: $(OBJS)
	$(FC) $(LDFLAGS) -o $@ $^

%.o: %.f90
	$(FC) $(FFLAGS) -c $<

clean:
	rm -f $(OBJS) *.mod



