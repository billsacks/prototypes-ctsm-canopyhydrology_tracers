subroutine CanopyHydrology

  do f = 1, num_nolakep
     p = filter_nolakep(f)
     g = patch%gridcell(p)
     l = patch%landunit(p)
     c = patch%column(p)

     ! Canopy interception and precipitation onto ground surface
     ! Add precipitation to leaf water

     if (lun%itype(l)==istsoil .or. lun%itype(l)==istwet .or. lun%urbpoi(l) .or. &
          lun%itype(l)==istcrop) then

        qflx_snocanfall(p) = 0._r8      ! rate of just snow canopy fall
        qflx_liqcanfall(p) = 0._r8
        qflx_snowindunload(p) = 0._r8
        qflx_snotempunload(p) = 0._r8
        snounload(p)=0._r8
        qflx_through_snow(p) = 0._r8 ! rain precipitation direct through canopy
        qflx_through_rain(p) = 0._r8 ! snow precipitation direct through canopy
        qflx_prec_intr(p) = 0._r8    ! total intercepted precipitation


        if (col%itype(c) /= icol_sunwall .and. col%itype(c) /= icol_shadewall) then
           if (frac_veg_nosno(p) == 1 .and. (forc_rain(c) + forc_snow(c) + qflx_irrig_sprinkler(p)) > 0._r8) then

              ! total liquid water inputs above canopy
              qflx_liq_above_canopy(p) = forc_rain(c)+ qflx_irrig_sprinkler(p)

              ! Coefficient of interception
              if(use_clm5_fpi) then 
                 fpi = interception_fraction * tanh(elai(p) + esai(p))
              else
                 fpi = 0.25_r8*(1._r8 - exp(-0.5_r8*(elai(p) + esai(p))))
              endif

              snocanmx = 60._r8*dewmx(p) * (elai(p) + esai(p))  ! 6*(LAI+SAI)
              liqcanmx = dewmx(p) * (elai(p) + esai(p))

              fpisnow = (1._r8 - exp(-0.5_r8*(elai(p) + esai(p))))  ! max interception of 1
              ! Direct throughfall
              qflx_through_snow(p) = forc_snow(c) * (1._r8-fpisnow)
              qflx_through_rain(p) = qflx_liq_above_canopy(p) * (1._r8-fpi)

              ! Intercepted precipitation [mm/s]
              qflx_prec_intr(p) = forc_snow(c)*fpisnow + qflx_liq_above_canopy(p)*fpi
              ! storage of intercepted snowfall, rain, and dew
              snocan(p) = max(0._r8, snocan(p) + dtime*forc_snow(c)*fpisnow)
              liqcan(p) = max(0._r8, liqcan(p) + dtime*qflx_liq_above_canopy(p)*fpi)

              ! Initialize rate of canopy runoff and snow falling off canopy
              qflx_snocanfall(p) = 0._r8
              qflx_liqcanfall(p) = 0._r8
              qflx_snowindunload(p) = 0._r8
              qflx_snotempunload(p) = 0._r8
              snounload(p)=0._r8

              if (forc_t(c) > tfrz) then ! Above freezing (Use t_veg?)
                 xliqrun = (liqcan(p) - liqcanmx)/dtime
                 if (xliqrun > 0._r8) then
                    qflx_liqcanfall(p) = xliqrun
                    liqcan(p) = liqcanmx
                 end if
              else ! Below freezing
                 xsnorun = (snocan(p) - snocanmx)/dtime
                 if (xsnorun > 0._r8) then ! exceeds snow capacity
                    qflx_snocanfall(p) = xsnorun
                    snocan(p) = snocanmx
                 end if
              end if
           end if
        end if

     else if (lun%itype(l)==istice_mec) then

        qflx_through_snow(p) = 0._r8
        qflx_through_rain(p) = 0._r8
        qflx_prec_intr(p)    = 0._r8
        snocan(p)            = 0._r8
        liqcan(p)            = 0._r8
        qflx_snocanfall(p)    = 0._r8
        qflx_liqcanfall(p)    = 0._r8
        qflx_snowindunload(p) = 0._r8 
        qflx_snotempunload(p) = 0._r8 
        snounload(p)=0._r8

     end if

     ! Precipitation onto ground (kg/(m2 s))

     if (col%itype(c) /= icol_sunwall .and. col%itype(c) /= icol_shadewall) then
        if (frac_veg_nosno(p) == 0) then
           qflx_prec_grnd_snow(p) = forc_snow(c)
           qflx_prec_grnd_rain(p) = forc_rain(c) + qflx_irrig_sprinkler(p)
        else
           qflx_snowindunload(p)=0._r8 
           qflx_snotempunload(p)=0._r8 
           snounload(p)=0._r8
           if (snocan(p) > 0._r8) then
              qflx_snotempunload(p) = max(0._r8,snocan(p)*(forc_t(c)-270.15_r8)/1.87e5_r8) 
              qflx_snowindunload(p) = 0.5_r8*snocan(p)*forc_wind(g)/1.56e5_r8 
              snounload(p) = (qflx_snowindunload(p)+qflx_snotempunload(p))*dtime ! total canopy unloading in timestep
              if ( snounload(p) > snocan(p) ) then ! Limit unloading to snow in canopy
                 snounload(p) = snocan(p)
              end if
              snocan(p) = snocan(p) - snounload(p)
           endif
           qflx_prec_grnd_snow(p) = qflx_through_snow(p) + qflx_snocanfall(p)  + snounload(p)/dtime
           qflx_prec_grnd_rain(p) = qflx_through_rain(p) + qflx_liqcanfall(p) 
        end if
        ! Urban sunwall and shadewall have no intercepted precipitation
     else
        qflx_prec_grnd_snow(p) = 0.
        qflx_prec_grnd_rain(p) = 0.
     end if

     ! Add irrigation water directly onto ground (bypassing canopy interception)
     ! Note that it's still possible that (some of) this irrigation water will runoff (as runoff is computed later)
     qflx_prec_grnd_rain(p) = qflx_prec_grnd_rain(p) + qflx_irrig_drip(p)

     qflx_prec_grnd(p) = qflx_prec_grnd_snow(p) + qflx_prec_grnd_rain(p)

     qflx_snow_grnd_patch(p) = qflx_prec_grnd_snow(p)           ! ice onto ground (mm/s)
     qflx_rain_grnd(p)       = qflx_prec_grnd_rain(p)           ! liquid water onto ground (mm/s)

  end do ! (end patch loop)
end subroutine CanopyHydrology
