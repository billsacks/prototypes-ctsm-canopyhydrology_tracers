! This implementation passes individual arrays explicitly for ALL routines. This adds long
! argument lists for the tracer-only routines, but keeps the tracer-only routines
! consistent with the bulk-only or bulk-and-tracer routines. (I have only implemented this
! for the first few routines.) In addition to the consistency, I like this because I think
! this will make it more likely for the sources (for state updates & tracer updates) to
! stay in sync: I think it will be easier for someone to notice that they need to change
! the argument list for the tracer flux update routine.
!
! Note that, for the argument names to the tracer flux updating routine, I'm using
! prefixes of bulk_ and trac_. Using prefixes makes it more obvious if you accidentally
! use a tracer variable where bulk should go, or vice versa. Using 'trac_' makes the
! tracer variables the same number of characters as the bulk, letting you see that things
! are consistent at a glance.

subroutine CanopyHydrology(bounds, num_soilp, filter_soilp, patch, water_inst)
  logical  :: check_point_for_interception_and_excess(bounds%begp:bounds%endp)
  real(r8) :: qflx_liq_above_canopy_patch(bounds%begp:bounds%endp) ! liquid water input above canopy (rain plus irrigation) [mm/s]
  real(r8) :: tracer_qflx_liq_above_canopy_patch(bounds%begp:bounds%endp) ! For one tracer: liquid water input above canopy (rain plus irrigation) [mm/s]
  real(r8) :: forc_snow_patch(bounds%begp:bounds%endp)
  real(r8) :: tracer_forc_snow_patch(bounds%begp:bounds%endp)

  associate( &
       b_wateratm2lnd_inst => water_inst%wateratm2lndbulk_inst, &
       b_waterflux_inst    => water_inst%waterfluxbulk_inst, &
       begp => bounds%begp, &
       endp => bounds%endp, &
       begc => bounds%begc, &
       endc => bounds%endc  &
       )

  ! Note about filters: I'm pretty sure that I'm missing some settings that need to be
  ! done outside of the soil filter.

  ! Compute patch-level precipitation inputs for bulk water
  call SumFlux_TopOfCanopyInputs(bounds, num_soilp, filter_soilp, &
       ! Inputs
       forc_rain             = b_wateratm2lnd_inst%forc_rain_col(begc:endc), &
       qflx_irrig_sprinkler  = b_waterflux_inst%qflx_irrig_sprinkler_patch(begp:endp), &
       forc_snow_col         = b_wateratm2lnd_inst%forc_snow_col(begc:endc), &
       ! Outputs
       qflx_liq_above_canopy = qflx_liq_above_canopy_patch(begp:endp), &
       forc_snow_patch       = forc_snow_patch(begp:endp))

  ! Compute canopy interception and throughfall for bulk water
  call BulkFlux_CanopyInterceptionAndThroughfall(bounds, num_soilp, filter_soilp, &
       ! Inputs
       frac_veg_nosno        = canopystate_inst%frac_veg_nosno_patch(begp:endp), &
       elai                  = canopystate_inst%elai_patch(begp:endp), &
       esai                  = canopystate_inst%esai_patch(begp:endp), &
       forc_snow             = forc_snow_patch(begp:endp), &
       qflx_liq_above_canopy = qflx_liq_above_canopy_patch(begp:endp), &
       ! Outputs
       qflx_through_snow     = b_waterflux_inst%qflx_through_snow_patch(begp:endp), &
       qflx_through_rain     = b_waterflux_inst%qflx_through_rain_patch(begp:endp), &
       qflx_intercepted_snow = b_waterflux_inst%qflx_intercepted_snow_patch(begp:endp), &
       qflx_intercepted_rain = b_waterflux_inst%qflx_intercepted_rain_patch(begp:endp), &
       check_point_for_interception_and_excess = check_point_for_interception_and_excess(begp:endp))

  ! Calculate canopy interception and throughfall for each tracer
  do i = water_inst%tracers_beg, water_inst%tracers_end
     associate(w => water_inst%bulk_and_tracers(i))
     call SumFlux_TopOfCanopyInputs(bounds, num_soilp, filter_soilp, &
          ! Inputs
          forc_rain             = w%wateratm2lnd_inst%forc_rain_col(begc:endc), &
          qflx_irrig_sprinkler  = w%waterflux_inst%qflx_irrig_sprinkler_patch(begp:endp), &
          forc_snow_col         = w%wateratm2lnd_inst%forc_snow_col(begc:endc), &
          ! Outputs
          qflx_liq_above_canopy = tracer_qflx_liq_above_canopy_patch(begp:endp), &
          forc_snow_patch       = tracer_forc_snow_patch(begp:endp))

     call TracerFlux_CanopyInterceptionAndThroughfall(bounds, num_soilp, filter_soilp, &
          ! Inputs
          bulk_forc_snow             = forc_snow_patch(begp:endp), &
          bulk_qflx_liq_above_canopy = qflx_liq_above_canopy_patch(begp:endp), &
          bulk_qflx_through_snow     = b_waterflux_inst%qflx_through_snow_patch(begp:endp), &
          bulk_qflx_intercepted_snow = b_waterflux_inst%qflx_intercepted_snow_patch(begp:endp), &
          bulk_qflx_through_rain     = b_waterflux_inst%qflx_through_rain_patch(begp:endp), &
          bulk_qflx_intercepted_rain = b_waterflux_inst%qflx_intercepted_rain_patch(begp:endp), &
          trac_forc_snow             = tracer_forc_snow_patch(begp:endp), &
          trac_qflx_liq_above_canopy = tracer_qflx_liq_above_canopy_patch(begp:endp), &
          ! Outputs
          trac_qflx_through_snow     = w%waterflux_inst%qflx_through_snow_patch(begp:endp), &
          trac_qflx_intercepted_snow = w%waterflux_inst%qflx_intercepted_snow_patch(begp:endp), &
          trac_qflx_through_rain     = w%waterflux_inst%qflx_through_rain_patch(begp:endp), &
          trac_qflx_intercepted_rain = w%waterflux_inst%qflx_intercepted_rain_patch(begp:endp))
     end associate
  end do

  ! Update snocan and liqcan based on interception, for bulk water and each tracer
  do i = water_inst%bulk_and_tracers_beg, water_inst%bulk_and_tracers_end
     associate(w => water_inst%bulk_and_tracers(i))
     call UpdateState_AddInterceptionToCanopy(bounds, num_soilp, filter_soilp, &
          ! Inputs
          qflx_intercepted_snow = w%waterflux_inst%qflx_intercepted_snow_patch(begp:endp), &
          qflx_intercepted_rain = w%waterflux_inst%qflx_intercepted_rain_patch(begp:endp), &
          ! Outputs
          snocan                = w%waterstate_inst%snocan_patch(begp:endp), &
          liqcan                = w%waterstate_inst%liqcan_patch(begp:endp))
     end associate
  end do

  ! Compute runoff from canopy due to exceeding maximum storage, for bulk
  call BulkFlux_CanopyExcess(bounds, num_soilp, filter_soilp, &
       waterstatebulk_inst, &  ! and some other inputs...
       waterfluxbulk_inst, &
       check_point_for_interception_and_excess = check_point_for_interception_and_excess(begp:endp))

  ! Calculate runoff from canopy due to exceeding maximum storage, for each tracer
  do i = water_inst%tracers_beg, water_inst%tracers_end
     call TracerFlux_CanopyExcess(bounds, num_soilp, filter_soilp, &
          water_inst%waterstatebulk_inst, water_inst%waterfluxbulk_inst, &
          water_inst%bulk_and_tracers(i)%waterstate_inst, &
          water_inst%bulk_and_tracers(i)%waterflux_inst)
  end do

  ! Update snocan and liqcan based on canfall, for bulk water and each tracer
  do i = water_inst%bulk_and_tracers_beg, water_inst%bulk_and_tracers_end
     call UpdateState_RemoveCanfallFromCanopy(num_soilp, filter_soilp, &
          water_inst%bulk_and_tracers(i)%waterflux_inst, &
          water_inst%bulk_and_tracers(i)%waterstate_inst)
  end do

  ! Compute snow unloading for bulk
  call BulkFlux_SnowUnloading(num_soilp, filter_soilp, &
       waterstatebulk_inst, &  ! and some other inputs...
       waterfluxbulk_inst)

  ! Compute snow unloading for each tracer
  do i = water_inst%tracers_beg, water_inst%tracers_end
     call TracerFlux_SnowUnloading(bounds, num_soilp, filter_soilp, &
          water_inst%waterstatebulk_inst, water_inst%waterfluxbulk_inst, &
          water_inst%bulk_and_tracers(i)%waterstate_inst, &
          water_inst%bulk_and_tracers(i)%waterflux_inst)
  end do

  ! Update snocan based on snow unloading, for bulk water and each tracer
  do i = water_inst%bulk_and_tracers_beg, water_inst%bulk_and_tracers_end
     call UpdateState_RemoveSnowUnloading(num_soilp, filter_soilp, &
          water_inst%bulk_and_tracers(i)%waterflux_inst, &
          water_inst%bulk_and_tracers(i)%waterstate_inst)
  end do

  ! Compute summed fluxes onto ground, for bulk water and each tracer
  do i = water_inst%bulk_and_tracers_beg, water_inst%bulk_and_tracers_end
     call SumFlux_FluxesOntoGround(num_soilp, filter_soilp, &
          water_inst%bulk_and_tracers(i)%waterflux_inst)
  end do

  end associate

end subroutine CanopyHydrology

subroutine SumFlux_TopOfCanopyInputs(bounds, num_soilp, filter_soilp, &
     forc_rain, qflx_irrig_sprinkler, forc_snow_col, &
     qflx_liq_above_canopy, forc_snow_patch)
  ! Compute patch-level precipitation inputs for bulk water or one tracer
  type(bounds_type), intent(in) :: bounds
  integer, intent(in) :: num_soilp
  integer, intent(in) :: filter_soilp(:)
  real(r8), intent(in) :: forc_rain( bounds%begc: )
  real(r8), intent(in) :: qflx_irrig_sprinkler( bounds%begp: )
  real(r8), intent(in) :: forc_snow_col( bounds%begc: )
  real(r8), intent(inout) :: qflx_liq_above_canopy( bounds%begp: )
  real(r8), intent(inout) :: forc_snow_patch( bounds%begp: )

  SHR_ASSERT_FL((ubound(forc_rain, 1) == bounds%endc), sourcefile, __LINE__)
  SHR_ASSERT_FL((ubound(qflx_irrig_sprinkler, 1) == bounds%endp), sourcefile, __LINE__)
  SHR_ASSERT_FL((ubound(forc_snow_col, 1) == bounds%endc), sourcefile, __LINE__)
  SHR_ASSERT_FL((ubound(qflx_liq_above_canopy, 1) == bounds%endp), sourcefile, __LINE__)
  SHR_ASSERT_FL((ubound(forc_snow_patch, 1) == bounds%endp), sourcefile, __LINE__)

  do fp = 1, num_soilp
     p = filter_soilp(fp)
     c = patch%column(p)

     qflx_liq_above_canopy_patch(p) = forc_rain_col(c) + qflx_irrig_sprinkler_patch(p)
     forc_snow_patch(p) = forc_snow_col(c)
  end do
end subroutine SumFlux_TopOfCanopyInputs

subroutine BulkFlux_CanopyInterceptionAndThroughfall(bounds, num_soilp, filter_soilp, &
     frac_veg_nosno, elai, esai, forc_snow, qflx_liq_above_canopy, &
     qflx_through_snow, qflx_through_rain, &
     qflx_intercepted_snow, qflx_intercepted_rain, &
     check_point_for_interception_and_excess)
  ! Compute canopy interception and throughfall for bulk water
  type(bounds_type), intent(in) :: bounds
  integer, intent(in) :: num_soilp
  integer, intent(in) :: filter_soilp(:)
  real(r8), intent(in) :: frac_veg_nosno( bounds%begp: )
  real(r8), intent(in) :: elai( bounds%begp: )
  real(r8), intent(in) :: esai( bounds%begp: )
  real(r8), intent(in) :: forc_snow( bounds%begp: )
  real(r8), intent(in) :: qflx_liq_above_canopy( bounds%begp: )

  real(r8), intent(inout) :: qflx_through_snow( bounds%begp: )
  real(r8), intent(inout) :: qflx_through_rain( bounds%begp: )
  real(r8), intent(inout) :: qflx_intercepted_snow( bounds%begp: )
  real(r8), intent(inout) :: qflx_intercepted_rain( bounds%begp: )
  logical , intent(inout) :: check_point_for_interception_and_excess( bounds%begp: )

  SHR_ASSERT_FL((ubound(frac_veg_nosno, 1) == bounds%endp), sourcefile, __LINE__)
  SHR_ASSERT_FL((ubound(elai, 1) == bounds%endp), sourcefile, __LINE__)
  SHR_ASSERT_FL((ubound(esai, 1) == bounds%endp), sourcefile, __LINE__)
  SHR_ASSERT_FL((ubound(forc_snow, 1) == bounds%endp), sourcefile, __LINE__)
  SHR_ASSERT_FL((ubound(qflx_liq_above_canopy, 1) == bounds%endp), sourcefile, __LINE__)
  SHR_ASSERT_FL((ubound(qflx_through_snow, 1) == bounds%endp), sourcefile, __LINE__)
  SHR_ASSERT_FL((ubound(qflx_through_rain, 1) == bounds%endp), sourcefile, __LINE__)
  SHR_ASSERT_FL((ubound(qflx_intercepted_snow, 1) == bounds%endp), sourcefile, __LINE__)
  SHR_ASSERT_FL((ubound(qflx_intercepted_rain, 1) == bounds%endp), sourcefile, __LINE__)
  SHR_ASSERT_FL((ubound(check_point_for_interception_and_excess, 1) == bounds%endp), sourcefile, __LINE__)

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
end subroutine BulkFlux_CanopyInterceptionAndThroughfall

subroutine TracerFlux_CanopyInterceptionAndThroughfall(bounds, num_soilp, filter_soilp, &
     bulk_forc_snow, bulk_qflx_liq_above_canopy, &
     bulk_qflx_through_snow, bulk_qflx_intercepted_snow, &
     bulk_qflx_through_rain, bulk_qflx_intercepted_rain, &
     trac_forc_snow, trac_qflx_liq_above_canopy, &
     trac_qflx_through_snow, trac_qflx_intercepted_snow, &
     trac_qflx_through_rain, trac_qflx_intercepted_rain)
  ! Calculate canopy interception and throughfall for one tracer
  type(bounds_type), intent(in) :: bounds
  integer, intent(in) :: num_soilp
  integer, intent(in) :: filter_soilp(:)
  real(r8), intent(in) :: bulk_forc_snow( bounds%begp: )
  real(r8), intent(in) :: bulk_qflx_liq_above_canopy( bounds%begp: )
  real(r8), intent(in) :: bulk_qflx_through_snow( bounds%begp: )
  real(r8), intent(in) :: bulk_qflx_intercepted_snow( bounds%begp: )
  real(r8), intent(in) :: bulk_qflx_through_rain( bounds%begp: )
  real(r8), intent(in) :: bulk_qflx_intercepted_rain( bounds%begp: )
  real(r8), intent(in) :: trac_forc_snow( bounds%begp: )
  real(r8), intent(in) :: trac_qflx_liq_above_canopy( bounds%begp: )
  real(r8), intent(inout) :: trac_qflx_through_snow( bounds%begp: )
  real(r8), intent(inout) :: trac_qflx_intercepted_snow( bounds%begp: )
  real(r8), intent(inout) :: trac_qflx_through_rain( bounds%begp: )
  real(r8), intent(inout) :: trac_qflx_intercepted_rain( bounds%begp: )

  SHR_ASSERT_FL((ubound(bulk_forc_snow, 1) == bounds%endp), sourcefile, __LINE__)
  SHR_ASSERT_FL((ubound(bulk_qflx_liq_above_canopy, 1) == bounds%endp), sourcefile, __LINE__)
  SHR_ASSERT_FL((ubound(bulk_qflx_through_snow, 1) == bounds%endp), sourcefile, __LINE__)
  SHR_ASSERT_FL((ubound(bulk_qflx_intercepted_snow, 1) == bounds%endp), sourcefile, __LINE__)
  SHR_ASSERT_FL((ubound(bulk_qflx_through_rain, 1) == bounds%endp), sourcefile, __LINE__)
  SHR_ASSERT_FL((ubound(bulk_qflx_intercepted_rain, 1) == bounds%endp), sourcefile, __LINE__)
  SHR_ASSERT_FL((ubound(trac_forc_snow, 1) == bounds%endp), sourcefile, __LINE__)
  SHR_ASSERT_FL((ubound(trac_qflx_liq_above_canopy, 1) == bounds%endp), sourcefile, __LINE__)
  SHR_ASSERT_FL((ubound(trac_qflx_through_snow, 1) == bounds%endp), sourcefile, __LINE__)
  SHR_ASSERT_FL((ubound(trac_qflx_intercepted_snow, 1) == bounds%endp), sourcefile, __LINE__)
  SHR_ASSERT_FL((ubound(trac_qflx_through_rain, 1) == bounds%endp), sourcefile, __LINE__)
  SHR_ASSERT_FL((ubound(trac_qflx_intercepted_rain, 1) == bounds%endp), sourcefile, __LINE__)

  call CalcTracerFromBulk( &
       lb            = begp, &
       num_pts       = num_soilp, &
       filter_pts    = filter_soilp, &
       bulk_source   = bulk_forc_snow(begp:endp), &
       bulk_val      = bulk_qflx_through_snow(begp:endp), &
       tracer_source = trac_forc_snow(begp:endp), &
       tracer_val    = trac_qflx_through_snow(begp:endp))

  call CalcTracerFromBulk( &
       lb            = begp, &
       num_pts       = num_soilp, &
       filter_pts    = filter_soilp, &
       bulk_source   = bulk_forc_snow(begp:endp), &
       bulk_val      = bulk_qflx_intercepted_snow(begp:endp), &
       tracer_source = trac_forc_snow(begp:endp), &
       tracer_val    = trac_qflx_intercepted_snow(begp:endp))

  call CalcTracerFromBulk( &
       lb            = begp, &
       num_pts       = num_soilp, &
       filter_pts    = filter_soilp, &
       bulk_source   = bulk_qflx_liq_above_canopy(begp:endp), &
       bulk_val      = bulk_qflx_through_rain(begp:endp), &
       tracer_source = trac_qflx_liq_above_canopy(begp:endp), &
       tracer_val    = trac_qflx_through_rain(begp:endp))

  call CalcTracerFromBulk( &
       lb            = begp, &
       num_pts       = num_soilp, &
       filter_pts    = filter_soilp, &
       bulk_source   = bulk_qflx_liq_above_canopy(begp:endp), &
       bulk_val      = bulk_qflx_intercepted_rain(begp:endp), &
       tracer_source = trac_qflx_liq_above_canopy(begp:endp), &
       tracer_val    = trac_qflx_intercepted_rain(begp:endp))
end subroutine TracerFlux_CanopyInterceptionAndThroughfall

subroutine UpdateState_AddInterceptionToCanopy(bounds, num_soilp, filter_soilp, &
     qflx_intercepted_snow, qflx_intercepted_rain, snocan, liqcan)
  ! Update snocan and liqcan based on interception, for bulk or one tracer
  integer, intent(in) :: num_soilp
  integer, intent(in) :: filter_soilp(:)
  real(r8), intent(in) :: qflx_intercepted_snow( bounds%begp: )
  real(r8), intent(in) :: qflx_intercepted_rain( bounds%begp: )
  real(r8), intent(in) :: snocan( bounds%begp: )
  real(r8), intent(in) :: liqcan( bounds%begp: )

  SHR_ASSERT_FL((ubound(qflx_intercepted_snow, 1) == bounds%endp), sourcefile, __LINE__)
  SHR_ASSERT_FL((ubound(qflx_intercepted_rain, 1) == bounds%endp), sourcefile, __LINE__)
  SHR_ASSERT_FL((ubound(snocan, 1) == bounds%endp), sourcefile, __LINE__)
  SHR_ASSERT_FL((ubound(liqcan, 1) == bounds%endp), sourcefile, __LINE__)

  do fp = 1, num_soilp
     p = filter_soilp(fp)

     snocan(p) = max(0._r8, snocan(p) + dtime * qflx_intercepted_snow(p))
     liqcan(p) = max(0._r8, liqcan(p) + dtime * qflx_intercepted_rain(p))
  end do
end subroutine UpdateState_AddInterceptionToCanopy

subroutine BulkFlux_CanopyExcess(bounds, num_soilp, filter_soilp, &
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

end subroutine BulkFlux_CanopyExcess

subroutine TracerFlux_CanopyExcess(bounds, num_soilp, filter_soilp, &
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
end subroutine TracerFlux_CanopyExcess

subroutine UpdateState_RemoveCanfallFromCanopy(num_soilp, filter_soilp, waterflux_inst, waterstate_inst)
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
end subroutine UpdateState_RemoveCanfallFromCanopy

subroutine BulkFlux_SnowUnloading(num_soilp, filter_soilp, &
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
end subroutine BulkFlux_SnowUnloading

subroutine TracerFlux_SnowUnloading(bounds, num_soilp, filter_soilp, &
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
end subroutine TracerFlux_SnowUnloading

subroutine UpdateState_RemoveSnowUnloading(num_soilp, filter_soilp, waterflux_inst, waterstate_inst)
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
end subroutine UpdateState_RemoveSnowUnloading

subroutine SumFlux_FluxesOntoGround(num_soilp, filter_soilp, waterflux_inst)
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
end subroutine SumFlux_FluxesOntoGround

