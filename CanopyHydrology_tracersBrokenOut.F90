! This version just breaks out tracer updates. I don't love that this still has loops over
! tracers mixed with the science code, but this way everything done to the bulk is visible
! inline.

subroutine CanopyHydrology

  ! Top-level associate statement would go here. All associates here would refer to bulk
  ! quantities. So, when you see a bare variable in the code (not accessed via %), it
  ! refers to bulk water.

  ! Note about filters: I'm pretty sure that I'm missing some settings that need to be
  ! done outside of the soil filter.

  ! ------------------------------------------------------------------------
  ! Compute patch-level precipitation inputs for bulk and all tracers
  ! ------------------------------------------------------------------------

  do i = water_inst%bulk_and_tracers_beg, water_inst%bulk_and_tracers_end
     associate( &
          waterflux_inst => water_inst%bulk_and_tracers(i)%waterflux_inst, &
          wateratm2lnd_inst => water_inst%bulk_and_tracers(i)%wateratm2lnd_inst &
          )

     do fp = 1, num_soilp
        p = filter_soilp(fp)
        c = patch%column(p)

        waterflux_inst%qflx_liq_above_canopy_patch(p) = &
             wateratm2lnd_inst%forc_rain_col(c) + &
             waterflux_inst%qflx_irrig_sprinkler_patch(p)

        wateratm2lnd_inst%forc_snow_patch(p) = wateratm2lnd_inst%forc_snow_col(c)
     end do
     end associate
  end do

  ! ------------------------------------------------------------------------
  ! Compute canopy interception and throughfall for bulk
  ! ------------------------------------------------------------------------

  ! Note: Previously this used a nolakep filter and set fluxes for all special landunits
  ! other than lakes. I'm pretty sure the fluxes ended up as 0 everywhere other than
  ! veg/crop. Here we use a soilp filter; this should give the same answers as long as we
  ! initialize the relevant fluxes to 0 in initCold.
  do fp = 1, num_soilp
     p = filter_soilp(fp)
     check_point_for_interception_and_excess(p) = &
          (frac_veg_nosno(p) == 1 .and. (forc_snow(p) + qflx_liq_above_canopy(p)) > 0._r8)
     if (check_point_for_interception_and_excess(p)) then
        ! Coefficient of interception
        if (use_clm5_fpi) then 
           fpi = interception_fraction * tanh(elai(p) + esai(p))
        else
           fpi = 0.25_r8*(1._r8 - exp(-0.5_r8*(elai(p) + esai(p))))
        end if

        fpisnow = (1._r8 - exp(-0.5_r8*(elai(p) + esai(p))))  ! max interception of 1

        ! Direct throughfall
        qflx_through_snow(p) = forc_snow(p) * (1._r8-fpisnow)
        qflx_through_rain(p) = qflx_liq_above_canopy(p) * (1._r8-fpi)

        ! Canopy interception
        qflx_intercepted_snow(p) = forc_snow(p) * fpisnow
        qflx_intercepted_rain(p) = qflx_liq_above_canopy(p) * fpi

     else
        ! Note: setting qflx_through_snow = forc_snow and qflx_through_rain =
        ! qflx_liq_above_canopy could change answers from the earlier logic if either of
        ! these could ever be negative.
        qflx_through_snow(p) = forc_snow(p)
        qflx_through_rain(p) = qflx_liq_above_canopy(p)
        qflx_intercepted_snow(p) = 0._r8
        qflx_intercepted_rain(p) = 0._r8
     end if
  end do

  call TracerCanopyInterceptionAndThroughfall(bounds, num_soilp, filter_soilp, &
       water_inst)

  ! ------------------------------------------------------------------------
  ! Update snocan and liqcan based on interception, for bulk and all tracers
  ! ------------------------------------------------------------------------

  do i = water_inst%bulk_and_tracers_beg, water_inst%bulk_and_tracers_end
     associate( &
          waterstate_inst => water_inst%bulk_and_tracers(i)%waterstate_inst, &
          waterflux_inst => water_inst%bulk_and_tracers(i)%waterflux_inst &
          )

     do fp = 1, num_soilp
        p = filter_soilp(fp)

        waterstate_inst%snocan_patch(p) = max(0._r8, &
             waterstate_inst%snocan_patch(p) + &
             dtime * waterflux_inst%qflx_intercepted_snow_patch(p))

        waterstate_inst%liqcan_patch(p) = max(0._r8, &
             waterstate_inst%liqcan_patch(p) + &
             dtime * waterflux_inst%qflx_intercepted_rain_patch(p))
     end do
     end associate
  end do

  ! ------------------------------------------------------------------------
  ! Compute runoff from canopy due to exceeding maximum storage, for bulk
  ! ------------------------------------------------------------------------

  do fp = 1, num_soilp
     p = filter_soilp(fp)
     qflx_snocanfall(p) = 0._r8
     qflx_liqcanfall(p) = 0._r8

     if (check_point_for_interception_and_excess(p)) then
        liqcanmx = dewmx(p) * (elai(p) + esai(p))
        qflx_liqcanfall(p) = max((liqcan(p) - liqcanmx)/dtime, 0._r8)
        snocanmx = 60._r8*dewmx(p) * (elai(p) + esai(p))  ! 6*(LAI+SAI)
        qflx_snocanfall(p) = max((snocan(p) - snocanmx)/dtime, 0._r8)
     end if
  end do

  call TracerCanopyExcess(bounds, num_soilp, filter_soilp, water_inst)

  ! ------------------------------------------------------------------------
  ! Update snocan and liqcan based on canfall, for bulk and all tracers
  ! ------------------------------------------------------------------------

  do i = water_inst%bulk_and_tracers_beg, water_inst%bulk_and_tracers_end
     associate( &
          waterstate_inst => water_inst%bulk_and_tracers(i)%waterstate_inst, &
          waterflux_inst => water_inst%bulk_and_tracers(i)%waterflux_inst &
          )

     do fp = 1, num_soilp
        p = filter_soilp(fp)

        waterstate_inst%liqcan_patch(p) = waterstate_inst%liqcan_patch(p) - &
             dtime * waterflux_inst%qflx_liqcanfall_patch(p)

        waterstate_inst%snocan_patch(p) = waterstate_inst%snocan_patch(p) - &
             dtime * waterflux_inst%qflx_snocanfall_patch(p)
     end do
     end associate
  end do

  ! ------------------------------------------------------------------------
  ! Compute snow unloading for bulk
  ! ------------------------------------------------------------------------

  do fp = 1, num_soilp
     p = filter_soilp(fp)
     c = patch%column(p)
     qflx_snow_unload(p) = 0._r8
     if (frac_veg_nosno(p) == 1 .and. snocan(p) > 0._r8) then
        qflx_snow_temp_unload = max(0._r8,snocan(p)*(forc_t(c)-270.15_r8)/1.87e5_r8)
        qflx_snow_wind_unload = 0.5_r8*snocan(p)*forc_wind(g)/1.56e5_r8
        qflx_snow_unload(p) = min(qflx_snow_temp_unload + qflx_snow_wind_unload, snocan(p))
     end if
  end do

  call TracerSnowUnloading(bounds, num_soilp, filter_soilp, water_inst)

  ! ------------------------------------------------------------------------
  ! Update snocan based on snow unloading, for bulk and all tracers
  ! ------------------------------------------------------------------------

  do i = water_inst%bulk_and_tracers_beg, water_inst%bulk_and_tracers_end
     associate( &
          waterstate_inst => water_inst%bulk_and_tracers(i)%waterstate_inst, &
          waterflux_inst => water_inst%bulk_and_tracers(i)%waterflux_inst &
          )

     do fp = 1, num_soilp
        p = filter_soilp(fp)

        waterstate_inst%snocan_patch(p) = &
             waterstate_inst%snocan_patch(p) - &
             dtime * waterflux_inst%qflx_snow_unload_patch(p)
     end do
     end associate
  end do

  ! ------------------------------------------------------------------------
  ! Compute summed fluxes onto ground, for bulk and all tracers
  ! ------------------------------------------------------------------------

  do i = water_inst%bulk_and_tracers_beg, water_inst%bulk_and_tracers_end
     associate( &
          waterflux_inst => water_inst%bulk_and_tracers(i)%waterflux_inst &
          )

     do fp = 1, num_soilp
        p = filter_soilp(fp)

        waterflux_inst%qflx_snow_grnd_patch(p) = &
             waterflux_inst%qflx_through_snow_patch(p) + &
             waterflux_inst%qflx_snocanfall_patch(p) + &
             waterflux_inst%qflx_snow_unload_patch(p)

        waterflux_inst%qflx_rain_grnd_patch(p) = &
             waterflux_inst%qflx_through_rain_patch(p) + &
             waterflux_inst%qflx_liqcanfall_patch(p) + &
             waterflux_inst%qflx_irrig_drip_patch(p)
     end do
     end associate
  end do

end subroutine CanopyHydrology

subroutine TracerCanopyInterceptionAndThroughfall(bounds, num_soilp, filter_soilp, &
     water_inst)
  ! Calculate canopy interception and throughfall for all tracers
  type(bounds_type), intent(in) :: bounds
  integer, intent(in) :: num_soilp
  integer, intent(in) :: filter_soilp(:)
  type(water_type), intent(inout) :: water_inst

  do i = water_inst%tracers_beg, water_inst%tracers_end
     associate( &
          bulk   => water_inst%bulk_and_tracers(i_bulk), &
          tracer => water_inst%bulk_and_tracers(i) &
          )

     call CalcTracerFromBulk( &
          lb            = begp, &
          num_pts       = num_soilp, &
          filter_pts    = filter_soilp, &
          bulk_source   = bulk  %wateratm2lnd_inst%forc_snow_patch(begp:endp), &
          bulk_val      = bulk  %waterflux_inst%qflx_through_snow_patch(begp:endp), &
          tracer_source = tracer%wateratm2lnd_inst%forc_snow_patch(begp:endp), &
          tracer_val    = tracer%waterflux_inst%qflx_through_snow_patch(begp:endp))

     call CalcTracerFromBulk( &
          lb            = begp, &
          num_pts       = num_soilp, &
          filter_pts    = filter_soilp, &
          bulk_source   = bulk  %wateratm2lnd_inst%forc_snow_patch(begp:endp), &
          bulk_val      = bulk  %waterflux_inst%qflx_intercepted_snow_patch(begp:endp), &
          tracer_source = tracer%wateratm2lnd_inst%forc_snow_patch(begp:endp), &
          tracer_val    = tracer%waterflux_inst%qflx_intercepted_snow_patch(begp:endp))

     call CalcTracerFromBulk( &
          lb            = begp, &
          num_pts       = num_soilp, &
          filter_pts    = filter_soilp, &
          bulk_source   = bulk  %waterflux_inst%qflx_liq_above_canopy_patch(begp:endp), &
          bulk_val      = bulk  %waterflux_inst%qflx_through_rain_patch(begp:endp), &
          tracer_source = tracer%waterflux_inst%qflx_liq_above_canopy_patch(begp:endp), &
          tracer_val    = tracer%waterflux_inst%qflx_through_rain_patch(begp:endp))

     call CalcTracerFromBulk( &
          lb            = begp, &
          num_pts       = num_soilp, &
          filter_pts    = filter_soilp, &
          bulk_source   = bulk  %waterflux_inst%qflx_liq_above_canopy_patch(begp:endp), &
          bulk_val      = bulk  %waterflux_inst%qflx_intercepted_rain_patch(begp:endp), &
          tracer_source = tracer%waterflux_inst%qflx_liq_above_canopy_patch(begp:endp), &
          tracer_val    = tracer%waterflux_inst%qflx_intercepted_rain_patch(begp:endp))

     end associate
  end do
end subroutine TracerCanopyInterceptionAndThroughfall

subroutine TracerCanopyExcess(bounds, num_soilp, filter_soilp, water_inst)
  ! Calculate runoff from canopy due to exceeding maximum storage, for all tracers
  type(bounds_type), intent(in) :: bounds
  integer, intent(in) :: num_soilp
  integer, intent(in) :: filter_soilp(:)
  type(water_type), intent(inout) :: water_inst

  do i = water_inst%tracers_beg, water_inst%tracers_end
     associate( &
          bulk   => water_inst%bulk_and_tracers(i_bulk), &
          tracer => water_inst%bulk_and_tracers(i) &
          )

     call CalcTracerFromBulk( &
          lb            = begp, &
          num_pts       = num_soilp, &
          filter_pts    = filter_soilp, &
          bulk_source   = bulk  %waterstate_inst%liqcan_patch(begp:endp), &
          bulk_val      = bulk  %waterflux_inst%qflx_liqcanfall_patch(begp:endp), &
          tracer_source = tracer%waterstate_inst%liqcan_patch(begp:endp), &
          tracer_val    = tracer%waterflux_inst%qflx_liqcanfall_patch(begp:endp))

     call CalcTracerFromBulk( &
          lb            = begp, &
          num_pts       = num_soilp, &
          filter_pts    = filter_soilp, &
          bulk_source   = bulk  %waterstate_inst%snocan_patch(begp:endp), &
          bulk_val      = bulk  %waterflux_inst%qflx_snocanfall_patch(begp:endp), &
          tracer_source = tracer%waterstate_inst%snocan_patch(begp:endp), &
          tracer_val    = tracer%waterflux_inst%qflx_snocanfall_patch(begp:endp))
     end associate
  end do
end subroutine TracerCanopyExcess

subroutine TracerSnowUnloading(bounds, num_soilp, filter_soilp, water_inst)
  ! Compute snow unloading for all tracers
  type(bounds_type), intent(in) :: bounds
  integer, intent(in) :: num_soilp
  integer, intent(in) :: filter_soilp(:)
  type(water_type), intent(inout) :: water_inst

  do i = water_inst%tracers_beg, water_inst%tracers_end
     associate( &
          bulk   => water_inst%bulk_and_tracers(i_bulk), &
          tracer => water_inst%bulk_and_tracers(i) &
          )

     call CalcTracerFromBulk( &
          lb            = begp, &
          num_pts       = num_soilp, &
          filter_pts    = filter_soilp, &
          bulk_source   = bulk  %waterstate_inst%snocan_patch(begp:endp), &
          bulk_val      = bulk  %waterflux_inst%qflx_snow_unload_patch(begp:endp), &
          tracer_source = tracer%waterstate_inst%snocan_patch(begp:endp), &
          tracer_val    = tracer%waterflux_inst%qflx_snow_unload_patch(begp:endp))
     end associate
  end do
end subroutine TracerSnowUnloading
