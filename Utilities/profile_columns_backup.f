   zone       ! numbers start with 1 at the surface
   mass       ! m/Msun. mass coordinate of outer boundary of cell.
   logR       ! log10(radius/Rsun) at outer boundary of zone
   logT       ! log10(temperature) at center of zone
   logRho     ! log10(density) at center of zone
   logP       ! log10(pressure) at center of zone
   !x_mass_fraction_H
   !y_mass_fraction_He
   !z_mass_fraction_metals

   ! Stuff for reynolds criterion + christensen relations
   mlt_mixing_length ! mixing length for mlt (cm)
   conv_vel ! convection velocity (cm/sec)
   dlnRho_dlnT_const_Pgas
   cp ! specific heat at constant total pressure
   grav ! gravitational acceleration (cm sec^2)

   ! Grads
   grada ! dlnT_dlnP at constant S
   gradr ! dlnT/dlnP required for purely radiative transport
