subroutine CanopyHydrology

  ! Note about filters: I'm pretty sure that I'm missing some settings that need to be
  ! done outside of the soil filter.

  do i = water_inst%bulk_and_tracers_beg, water_inst%bulk_and_tracers_end
     ! A possible improvement for this and other calls that operate on bulk or
     ! bulk+tracers would be passing individual arrays rather than whole instances. This
     ! way, you could trace the data flow through the routine.
     call PrecipInputs(num_soilp, filter_soilp, &
          water_inst%bulk_and_tracers(i)%waterflux_inst, &
          water_inst%bulk_and_tracers(i)%wateratm2lnd_inst)
  end do

  call BulkCanopyInterceptionAndThroughfall(bounds, num_soilp, filter_soilp, &
       ! and some other inputs...
       water_inst%waterfluxbulk_inst, &
       check_point_for_interception_and_excess = check_point_for_interception_and_excess(begp:endp))

  do i = water_inst%tracers_beg, water_inst%tracers_end
     call TracerCanopyInterceptionAndThroughfall(bounds, num_soilp, filter_soilp, &
          water_inst%wateratm2lndbulk_inst, water_inst%waterfluxbulk_inst, &
          water_inst%bulk_and_tracers(i)%wateratm2lnd_inst, &
          water_inst%bulk_and_tracers(i)%waterflux_inst)
  end do

  do i = water_inst%bulk_and_tracers_beg, water_inst%bulk_and_tracers_end
     call AddInterceptionToCanopy(num_soilp, filter_soilp, &
          water_inst%bulk_and_tracers(i)%waterflux_inst, &
          water_inst%bulk_and_tracers(i)%waterstate_inst)
  end do

  call BulkCanopyExcess(bounds, num_soilp, filter_soilp, &
       waterstatebulk_inst, &  ! and some other inputs...
       waterfluxbulk_inst, &
       check_point_for_interception_and_excess = check_point_for_interception_and_excess(begp:endp))

  do i = water_inst%tracers_beg, water_inst%tracers_end
     call TracerCanopyExcess(bounds, num_soilp, filter_soilp, &
          water_inst%waterstatebulk_inst, water_inst%waterfluxbulk_inst, &
          water_inst%bulk_and_tracers(i)%waterstate_inst, &
          water_inst%bulk_and_tracers(i)%waterflux_inst)
  end do

  do i = water_inst%bulk_and_tracers_beg, water_inst%bulk_and_tracers_end
     call RemoveCanfallFromCanopy(num_soilp, filter_soilp, &
          water_inst%bulk_and_tracers(i)%waterflux_inst, &
          water_inst%bulk_and_tracers(i)%waterstate_inst)
  end do

  call BulkSnowUnloading(num_soilp, filter_soilp, &
       waterstatebulk_inst, &  ! and some other inputs...
       waterfluxbulk_inst)

  do i = water_inst%tracers_beg, water_inst%tracers_end
     call TracerSnowUnloading(bounds, num_soilp, filter_soilp, &
          water_inst%waterstatebulk_inst, water_inst%waterfluxbulk_inst, &
          water_inst%bulk_and_tracers(i)%waterstate_inst, &
          water_inst%bulk_and_tracers(i)%waterflux_inst)
  end do

  do i = water_inst%bulk_and_tracers_beg, water_inst%bulk_and_tracers_end
     call RemoveSnowUnloading(num_soilp, filter_soilp, &
          water_inst%bulk_and_tracers(i)%waterflux_inst, &
          water_inst%bulk_and_tracers(i)%waterstate_inst)
  end do

  do i = water_inst%bulk_and_tracers_beg, water_inst%bulk_and_tracers_end
     call FluxesOntoGround(num_soilp, filter_soilp, &
          water_inst%bulk_and_tracers(i)%waterflux_inst)
  end do

end subroutine CanopyHydrology

subroutine PrecipInputs(num_soilp, filter_soilp, waterflux_inst, wateratm2lnd_inst)
  ! Compute patch-level precipitation inputs for bulk water or one tracer
  integer, intent(in) :: num_soilp
  integer, intent(in) :: filter_soilp(:)
  class(waterflux_type), intent(inout) :: waterflux_inst
  class(wateratm2lnd_type), intent(inout) :: wateratm2lnd_inst

  ! Associate statements go here

  do fp = 1, num_soilp
     p = filter_soilp(fp)
     c = patch%column(p)

     qflx_liq_above_canopy_patch(p) = forc_rain_col(c) + qflx_irrig_sprinkler_patch(p)

     forc_snow_patch(p) = forc_snow_col(c)
  end do
end subroutine PrecipInputs

subroutine BulkCanopyInterceptionAndThroughfall(bounds, num_soilp, filter_soilp, &
     ! and some other inputs...
     waterfluxbulk_inst, check_point_for_interception_and_excess)
  ! Compute canopy interception and throughfall for bulk water
  type(bounds_type), intent(in) :: bounds
  integer, intent(in) :: num_soilp
  integer, intent(in) :: filter_soilp(:)
  ! And some other inputs...
  type(waterfluxbulk_type), intent(inout) :: waterfluxbulk_inst
  logical, intent(inout) :: check_point_for_interception_and_excess( bounds%begp: )

  ! Associate statements go here

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
end subroutine BulkCanopyInterceptionAndThroughfall

subroutine TracerCanopyInterceptionAndThroughfall(bounds, num_soilp, filter_soilp, &
     wateratm2lndbulk_inst, waterfluxbulk_inst, &
     wateratm2lnd_tracer_inst, waterflux_tracer_inst)
  ! Calculate canopy interception and throughfall for one tracer
  type(bounds_type), intent(in) :: bounds
  integer, intent(in) :: num_soilp
  integer, intent(in) :: filter_soilp(:)
  type(wateratm2lndbulk_type), intent(in) :: wateratm2lndbulk_inst
  type(waterfluxbulk_type), intent(in) :: waterfluxbulk_inst
  type(wateratm2lnd_type), intent(in) :: wateratm2lnd_tracer_inst
  type(waterflux_type), intent(inout) :: waterflux_tracer_inst

  call CalcTracerFromBulk( &
       lb            = begp, &
       num_pts       = num_soilp, &
       filter_pts    = filter_soilp, &
       bulk_source   = wateratm2lndbulk_inst%forc_snow(begp:endp), &
       bulk_val      = waterfluxbulk_inst%qflx_through_snow(begp:endp), &
       tracer_source = wateratm2lnd_tracer_inst%forc_snow_patch(begp:endp), &
       tracer_val    = waterflux_tracer_inst%qflx_through_snow_patch(begp:endp))

  call CalcTracerFromBulk( &
       lb            = begp, &
       num_pts       = num_soilp, &
       filter_pts    = filter_soilp, &
       bulk_source   = wateratm2lndbulk_inst%forc_snow(begp:endp), &
       bulk_val      = waterfluxbulk_inst%qflx_intercepted_snow(begp:endp), &
       tracer_source = wateratm2lnd_tracer_inst%forc_snow_patch(begp:endp), &
       tracer_val    = waterflux_tracer_inst%qflx_intercepted_snow_patch(begp:endp))

  call CalcTracerFromBulk( &
       lb            = begp, &
       num_pts       = num_soilp, &
       filter_pts    = filter_soilp, &
       bulk_source   = waterfluxbulk_inst%qflx_liq_above_canopy(begp:endp), &
       bulk_val      = waterfluxbulk_inst%qflx_through_rain(begp:endp), &
       tracer_source = waterflux_tracer_inst%qflx_liq_above_canopy_patch(begp:endp), &
       tracer_val    = waterflux_tracer_inst%qflx_through_rain_patch(begp:endp))

  call CalcTracerFromBulk( &
       lb            = begp, &
       num_pts       = num_soilp, &
       filter_pts    = filter_soilp, &
       bulk_source   = waterfluxbulk_inst%qflx_liq_above_canopy(begp:endp), &
       bulk_val      = waterfluxbulk_inst%qflx_intercepted_rain(begp:endp), &
       tracer_source = waterflux_tracer_inst%qflx_liq_above_canopy_patch(begp:endp), &
       tracer_val    = waterflux_tracer_inst%qflx_intercepted_rain_patch(begp:endp))
end subroutine TracerCanopyInterceptionAndThroughfall

subroutine AddInterceptionToCanopy(num_soilp, filter_soilp, waterflux_inst, waterstate_inst)
  ! Update snocan and liqcan based on interception, for bulk or one tracer
  integer, intent(in) :: num_soilp
  integer, intent(in) :: filter_soilp(:)
  class(waterflux_type), intent(in) :: waterflux_inst
  class(waterstate_type), intent(inout) :: waterstate_inst

  ! Associates go here

  do fp = 1, num_soilp
     p = filter_soilp(fp)

     snocan_patch(p) = max(0._r8, snocan_patch(p) + dtime * qflx_intercepted_snow_patch(p))
     liqcan_patch(p) = max(0._r8, liqcan_patch(p) + dtime * qflx_intercepted_rain_patch(p))
  end do
end subroutine AddInterceptionToCanopy

subroutine BulkCanopyExcess(bounds, num_soilp, filter_soilp, &
     waterstatebulk_inst, &  ! and some other inputs...
     waterfluxbulk_inst, &
     check_point_for_interception_and_excess)
  ! Compute runoff from canopy due to exceeding maximum storage, for bulk
  type(bounds_type), intent(in) :: bounds
  integer, intent(in) :: num_soilp
  integer, intent(in) :: filter_soilp(:)
  type(waterstatebulk_type), intent(in) :: waterstatebulk_inst
  ! And some other inputs...
  type(waterfluxbulk_type), intent(inout) :: waterfluxbulk_inst
  logical, intent(in) :: check_point_for_interception_and_excess( bounds%begp: )

  ! Associates go here

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

end subroutine BulkCanopyExcess

subroutine TracerCanopyExcess(bounds, num_soilp, filter_soilp, &
     waterstatebulk_inst, waterfluxbulk_inst, &
     waterstate_tracer_inst, waterflux_tracer_inst)
  ! Calculate runoff from canopy due to exceeding maximum storage, for one tracer
  type(bounds_type), intent(in) :: bounds
  integer, intent(in) :: num_soilp
  integer, intent(in) :: filter_soilp(:)
  type(waterstatebulk_type), intent(in) :: waterstatebulk_inst
  type(waterfluxbulk_type), intent(in) :: waterfluxbulk_inst
  type(waterstate_type), intent(in) :: waterstate_tracer_inst
  type(waterflux_type), intent(inout) :: waterflux_tracer_inst

  call CalcTracerFromBulk( &
       lb            = begp, &
       num_pts       = num_soilp, &
       filter_pts    = filter_soilp, &
       bulk_source   = waterstatebulk_inst%liqcan(begp:endp), &
       bulk_val      = waterfluxbulk_inst%qflx_liqcanfall(begp:endp), &
       tracer_source = waterstate_tracer_inst%liqcan_patch(begp:endp), &
       tracer_val    = waterflux_tracer_inst%qflx_liqcanfall_patch(begp:endp))

  call CalcTracerFromBulk( &
       lb            = begp, &
       num_pts       = num_soilp, &
       filter_pts    = filter_soilp, &
       bulk_source   = waterstatebulk_inst%snocan(begp:endp), &
       bulk_val      = waterfluxbulk_inst%qflx_snocanfall(begp:endp), &
       tracer_source = waterstate_tracer_inst%snocan_patch(begp:endp), &
       tracer_val    = waterflux_tracer_inst%qflx_snocanfall_patch(begp:endp))
end subroutine TracerCanopyExcess

subroutine RemoveCanfallFromCanopy(num_soilp, filter_soilp, waterflux_inst, waterstate_inst)
  ! Update snocan and liqcan based on canfall, for bulk or one tracer
  integer, intent(in) :: num_soilp
  integer, intent(in) :: filter_soilp(:)
  class(waterflux_type), intent(in) :: waterflux_inst
  class(waterstate_type), intent(inout) :: waterstate_inst

  ! Associates go here

  do fp = 1, num_soilp
     p = filter_soilp(fp)

     liqcan_patch(p) = liqcan_patch(p) - dtime * qflx_liqcanfall_patch(p)
     snocan_patch(p) = snocan_patch(p) - dtime * qflx_snocanfall_patch(p)
  end do
end subroutine RemoveCanfallFromCanopy

subroutine BulkSnowUnloading(num_soilp, filter_soilp, &
     waterstatebulk_inst, &  ! and some other inputs...
     waterfluxbulk_inst)
  ! Compute snow unloading for bulk
  integer, intent(in) :: num_soilp
  integer, intent(in) :: filter_soilp(:)
  type(waterstatebulk_type), intent(in) :: waterstatebulk_inst
  ! And some other inputs...
  type(waterfluxbulk_type), intent(inout) :: waterfluxbulk_inst

  ! Associates go here

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
end subroutine BulkSnowUnloading

subroutine TracerSnowUnloading(bounds, num_soilp, filter_soilp, &
     waterstatebulk_inst, waterfluxbulk_inst, &
     waterstate_tracer_inst, waterflux_tracer_inst)
  ! Compute snow unloading for one tracer
  type(bounds_type), intent(in) :: bounds
  integer, intent(in) :: num_soilp
  integer, intent(in) :: filter_soilp(:)
  type(waterstatebulk_type), intent(in) :: waterstatebulk_inst
  type(waterfluxbulk_type), intent(in) :: waterfluxbulk_inst
  type(waterstate_type), intent(in) :: waterstate_tracer_inst
  type(waterflux_type), intent(inout) :: waterflux_tracer_inst

  call CalcTracerFromBulk( &
       lb            = begp, &
       num_pts       = num_soilp, &
       filter_pts    = filter_soilp, &
       bulk_source   = waterstatebulk_inst%snocan(begp:endp), &
       bulk_val      = waterfluxbulk_inst%qflx_snow_unload(begp:endp), &
       tracer_source = waterstate_tracer_inst%snocan_patch(begp:endp), &
       tracer_val    = waterflux_tracer_inst%qflx_snow_unload_patch(begp:endp))
end subroutine TracerSnowUnloading

subroutine RemoveSnowUnloading(num_soilp, filter_soilp, waterflux_inst, waterstate_inst)
  ! Update snocan based on snow unloading, for bulk or one tracer
  integer, intent(in) :: num_soilp
  integer, intent(in) :: filter_soilp(:)
  class(waterflux_type), intent(in) :: waterflux_inst
  class(waterstate_type), intent(inout) :: waterstate_inst

  ! Associates go here

  do fp = 1, num_soilp
     p = filter_soilp(fp)

     snocan_patch(p) = snocan_patch(p) - dtime * qflx_snow_unload_patch(p)
  end do
end subroutine RemoveSnowUnloading

subroutine FluxesOntoGround(num_soilp, filter_soilp, waterflux_inst)
  ! Compute summed fluxes onto ground, for bulk or one tracer
  integer, intent(in) :: num_soilp
  integer, intent(in) :: filter_soilp(:)
  class(waterflux_type), intent(inout) :: waterflux_inst

  ! Associates go here

  do fp = 1, num_soilp
     p = filter_soilp(fp)

     qflx_snow_grnd_patch(p) = &
          qflx_through_snow_patch(p) + &
          qflx_snocanfall_patch(p) + &
          qflx_snow_unload_patch(p)

     qflx_rain_grnd_patch(p) = &
          qflx_through_rain_patch(p) + &
          qflx_liqcanfall_patch(p) + &
          qflx_irrig_drip_patch(p)
  end do
end subroutine FluxesOntoGround
