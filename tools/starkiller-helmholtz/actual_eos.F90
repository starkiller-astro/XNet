module actual_eos_module

    use xnet_types, only: dp
    use eos_type_module

    implicit none
    private

    character (len=64), public :: eos_name = "helmholtz"

    ! Runtime parameters
    logical, allocatable :: do_coulomb
    logical, allocatable :: input_is_constant

    !..for the tables, in general
    integer, parameter, private :: imax = 541, jmax = 201
    integer, allocatable :: itmax, jtmax
    real(dp), allocatable :: d(:), t(:)

    real(dp), allocatable :: tlo, thi, tstp, tstpi
    real(dp), allocatable :: dlo, dhi, dstp, dstpi

    real(dp), allocatable :: ttol, dtol

    !..for the helmholtz free energy tables
    real(dp), allocatable :: f(:,:), fd(:,:),                &
                             ft(:,:), fdd(:,:), ftt(:,:),    &
                             fdt(:,:), fddt(:,:), fdtt(:,:), &
                             fddtt(:,:)

    !..for the pressure derivative with density ables
    real(dp), allocatable :: dpdf(:,:), dpdfd(:,:),          &
                             dpdft(:,:), dpdfdt(:,:)

    !..for chemical potential tables
    real(dp), allocatable :: ef(:,:), efd(:,:),              &
                             eft(:,:), efdt(:,:)

    !..for the number density tables
    real(dp), allocatable :: xf(:,:), xfd(:,:),              &
                             xft(:,:), xfdt(:,:)

    !..for storing the differences
    real(dp), allocatable :: dt_sav(:), dt2_sav(:),          &
                             dti_sav(:), dt2i_sav(:),        &
                             dd_sav(:), dd2_sav(:),          &
                             ddi_sav(:), dd2i_sav(:)

    integer, parameter :: max_newton = 100

    ! 2006 CODATA physical constants

    ! Math constants
    real(dp), parameter :: pi       = 3.1415926535897932384_dp

    ! Physical constants
    real(dp), parameter :: h       = 6.6260689633e-27_dp
    real(dp), parameter :: hbar    = 0.5_dp * h/pi
    real(dp), parameter :: qe      = 4.8032042712e-10_dp
    real(dp), parameter :: avo_eos = 6.0221417930e23_dp
    real(dp), parameter :: clight  = 2.99792458e10_dp
    real(dp), parameter :: kerg    = 1.380650424e-16_dp
    real(dp), parameter :: ev2erg_eos  = 1.60217648740e-12_dp
    real(dp), parameter :: kev     = kerg/ev2erg_eos
    real(dp), parameter :: amu     = 1.66053878283e-24_dp
    real(dp), parameter :: me_eos  = 9.1093821545e-28_dp
    real(dp), parameter :: rbohr   = hbar*hbar/(me_eos * qe * qe)
    real(dp), parameter :: fine    = qe*qe/(hbar*clight)

    real(dp), parameter :: ssol    = 5.67051e-5_dp
    real(dp), parameter :: asol    = 4.0_dp * ssol / clight
    real(dp), parameter :: weinlam = h*clight/(kerg * 4.965114232_dp)
    real(dp), parameter :: weinfre = 2.821439372_dp*kerg/h

    ! Astronomical constants
    real(dp), parameter :: ly      = 9.460528e17_dp
    real(dp), parameter :: pc      = 3.261633_dp * ly

    ! Some other useful combinations of the constants
    real(dp), parameter :: sioncon = (2.0_dp * pi * amu * kerg)/(h*h)
    real(dp), parameter :: forth   = 4.0_dp/3.0_dp
    real(dp), parameter :: forpi   = 4.0_dp * pi
    real(dp), parameter :: forthpi = forth * pi
    real(dp), parameter :: kergavo = kerg * avo_eos
    real(dp), parameter :: ikavo   = 1.0_dp/kergavo
    real(dp), parameter :: asoli3  = asol/3.0_dp
    real(dp), parameter :: light2  = clight * clight

    ! Constants used for the Coulomb corrections
    real(dp), parameter :: a1    = -0.898004_dp
    real(dp), parameter :: b1    =  0.96786_dp
    real(dp), parameter :: c1    =  0.220703_dp
    real(dp), parameter :: d1    = -0.86097_dp
    real(dp), parameter :: e1    =  2.5269_dp
    real(dp), parameter :: a2    =  0.29561_dp
    real(dp), parameter :: b2    =  1.9885_dp
    real(dp), parameter :: c2    =  0.288675_dp
    real(dp), parameter :: onethird = 1.0_dp/3.0_dp
    real(dp), parameter :: esqu = qe * qe

    !$acc declare &
    !$acc create(tlo, thi, dlo, dhi) &
    !$acc create(tstp, tstpi, dstp, dstpi) &
    !$acc create(ttol, dtol) &
    !$acc create(itmax, jtmax, d, t) &
    !$acc create(f, fd, ft, fdd, ftt, fdt, fddt, fdtt, fddtt) &
    !$acc create(dpdf, dpdfd, dpdft, dpdfdt) &
    !$acc create(ef, efd, eft, efdt, xf, xfd, xft, xfdt)  &
    !$acc create(dt_sav, dt2_sav, dti_sav, dt2i_sav) &
    !$acc create(dd_sav, dd2_sav, ddi_sav, dd2i_sav) &
    !$acc create(do_coulomb, input_is_constant)

    public :: actual_eos, actual_eos_init, actual_eos_finalize, eos_supports_input_type
    public :: xnet_actual_eos, actual_eos_eta, actual_eos_cv

contains


    function eos_supports_input_type(input) result(supported)

        implicit none

        integer, intent(in) :: input
        logical :: supported

        if (input == eos_input_rt .or. &
            input == eos_input_rp .or. &
            input == eos_input_rh .or. &
            input == eos_input_re .or. &
            input == eos_input_tp .or. &
            input == eos_input_th .or. &
            input == eos_input_ps .or. &
            input == eos_input_ph) then

            supported = .true.

         else

            supported = .false.

         endif

       end function eos_supports_input_type



    !  Frank Timmes Helmholtz based Equation of State
    !  http://cococubed.asu.edu/

    !..given a temperature temp [K], density den [g/cm**3], and a composition
    !..characterized by abar and zbar, this routine returns most of the other
    !..thermodynamic quantities. of prime interest is the pressure [erg/cm**3],
    !..specific thermal energy [erg/gr], the entropy [erg/g/K], along with
    !..their derivatives with respect to temperature, density, abar, and zbar.
    !..other quantites such the normalized chemical potential eta (plus its
    !..derivatives), number density of electrons and positron pair (along
    !..with their derivatives), adiabatic indices, specific heats, and
    !..relativistically correct sound speed are also returned.
    !..
    !..this routine assumes planckian photons, an ideal gas of ions,
    !..and an electron-positron gas with an arbitrary degree of relativity
    !..and degeneracy. interpolation in a table of the helmholtz free energy
    !..is used to return the electron-positron thermodynamic quantities.
    !..all other derivatives are analytic.
    !..
    !..references: cox & giuli chapter 24 ; timmes & swesty apj 1999

    subroutine actual_eos(input, state)

        !$acc routine seq

        implicit none

        !..input arguments
        integer,      intent(in   ) :: input
        type (eos_t), intent(inout) :: state

        !..rows to store EOS data
        real(dp) :: temp_row, &
                    den_row, &
                    abar_row, &
                    zbar_row, &
                    ye_row, &
                    etot_row, &
                    ptot_row, &
                    cv_row, &
                    cp_row,  &
                    xne_row, &
                    xnp_row, &
                    etaele_row, &
                    detadt_row, &
                    pele_row, &
                    ppos_row, &
                    dpd_row,  &
                    dpt_row, &
                    dpa_row, &
                    dpz_row,  &
                    ded_row, &
                    det_row, &
                    dea_row,  &
                    dez_row,  &
                    stot_row, &
                    dsd_row, &
                    dst_row, &
                    htot_row, &
                    dhd_row, &
                    dht_row, &
                    dpe_row, &
                    dpdr_e_row, &
                    gam1_row, &
                    cs_row

        !..declare local variables

        logical :: single_iter, double_iter, converged
        integer :: var, dvar, var1, var2, iter
        real(dp) :: v_want
        real(dp) :: v1_want, v2_want
        real(dp) :: xnew, xtol, dvdx, smallx, error, v
        real(dp) :: v1, v2, dv1dt, dv1dr, dv2dt,dv2dr, delr, error1, error2, told, rold, tnew, rnew, v1i, v2i

        real(dp) :: x,y,z,zz,zzi,deni,tempi,xni,dxnidd,dxnida, &
                    dpepdt,dpepdd,deepdt,deepdd,dsepdd,dsepdt, &
                    dpraddd,dpraddt,deraddd,deraddt,dpiondd,dpiondt, &
                    deiondd,deiondt,dsraddd,dsraddt,dsiondd,dsiondt, &
                    dse,dpe,dsp,kt,ktinv,prad,erad,srad,pion,eion, &
                    sion,xnem,pele,eele,sele,pres,ener,entr,dpresdd, &
                    dpresdt,denerdd,denerdt,dentrdd,dentrdt,cv,cp, &
                    gam1,gam2,gam3,chit,chid,nabad,sound,etaele, &
                    detadt,detadd,xnefer,dxnedt,dxnedd,s, &
                    temp,den,abar,zbar,ytot1,ye,din


        !..for the interpolations
        integer :: iat,jat
        real(dp) :: free,df_d,df_t,df_tt,df_dt
        real(dp) :: xt,xd,mxt,mxd,fi(36), &
                    si0t,si1t,si2t,si0mt,si1mt,si2mt, &
                    si0d,si1d,si2d,si0md,si1md,si2md, &
                    dsi0t,dsi1t,dsi2t,dsi0mt,dsi1mt,dsi2mt, &
                    dsi0d,dsi1d,dsi2d,dsi0md,dsi1md,dsi2md, &
                    ddsi0t,ddsi1t,ddsi2t,ddsi0mt,ddsi1mt,ddsi2mt

        !..for the coulomb corrections
        real(dp) :: dsdd,dsda,lami,inv_lami,lamida,lamidd,     &
                    plasg,plasgdd,plasgdt,plasgda,plasgdz,     &
                    ecoul,decouldd,decouldt,decoulda,decouldz, &
                    pcoul,dpcouldd,dpcouldt,dpcoulda,dpcouldz, &
                    scoul,dscouldd,dscouldt,dscoulda,dscouldz

        real(dp) :: p_temp, e_temp
        real(dp) :: smallt, smalld

        call eos_get_small_temp(smallt)
        call eos_get_small_dens(smalld)

        temp_row = state % T
        den_row  = state % rho
        abar_row = state % abar
        zbar_row = state % zbar
        ye_row   = state % y_e

        ! Initial setup for iterations

        single_iter = .false.
        double_iter = .false.

        if (input .eq. eos_input_rt) then

          ! Nothing to do here.

        elseif (input .eq. eos_input_rh) then

          single_iter = .true.
          v_want = state % h
          var  = ienth
          dvar = itemp

        elseif (input .eq. eos_input_tp) then

          single_iter = .true.
          v_want = state % p
          var  = ipres
          dvar = idens

        elseif (input .eq. eos_input_rp) then

          single_iter = .true.
          v_want = state % p
          var  = ipres
          dvar = itemp

        elseif (input .eq. eos_input_re) then

          single_iter = .true.
          v_want = state % e
          var  = iener
          dvar = itemp

        elseif (input .eq. eos_input_ps) then

          double_iter = .true.
          v1_want = state % p
          v2_want = state % s
          var1 = ipres
          var2 = ientr

        elseif (input .eq. eos_input_ph) then

          double_iter = .true.
          v1_want = state % p
          v2_want = state % h
          var1 = ipres
          var2 = ienth

        elseif (input .eq. eos_input_th) then

          single_iter = .true.
          v_want = state % h
          var  = ienth
          dvar = idens

        endif

        ptot_row = 0.0_dp
        dpt_row = 0.0_dp
        dpd_row = 0.0_dp
        dpa_row = 0.0_dp
        dpz_row = 0.0_dp
        dpe_row = 0.0_dp
        dpdr_e_row = 0.0_dp

        etot_row = 0.0_dp
        det_row = 0.0_dp
        ded_row = 0.0_dp
        dea_row = 0.0_dp
        dez_row = 0.0_dp

        stot_row = 0.0_dp
        dst_row = 0.0_dp
        dsd_row = 0.0_dp

        htot_row = 0.0_dp
        dhd_row = 0.0_dp
        dht_row = 0.0_dp

        pele_row = 0.0_dp
        ppos_row = 0.0_dp

        xne_row = 0.0_dp
        xnp_row = 0.0_dp

        etaele_row = 0.0_dp
        detadt_row = 0.0_dp

        cv_row = 0.0_dp
        cp_row = 0.0_dp
        cs_row = 0.0_dp
        gam1_row = 0.0_dp

        converged = .false.

        if (input .eq. eos_input_rt) converged = .true.

        do iter = 1, max_newton

           temp  = temp_row
           den   =  den_row
           abar  = abar_row
           zbar  = zbar_row

           ytot1 = 1.0_dp / abar
           ye    = ye_row
           din   = ye * den

           !..initialize
           deni    = 1.0_dp/den
           tempi   = 1.0_dp/temp
           kt      = kerg * temp
           ktinv   = 1.0_dp/kt

           !..radiation section:
           prad    = asoli3 * temp * temp * temp * temp
           dpraddd = 0.0_dp
           dpraddt = 4.0_dp * prad*tempi

           erad    = 3.0_dp * prad*deni
           deraddd = -erad*deni
           deraddt = 3.0_dp * dpraddt*deni

           srad    = (prad*deni + erad)*tempi
           dsraddd = (dpraddd*deni - prad*deni*deni + deraddd)*tempi
           dsraddt = (dpraddt*deni + deraddt - srad)*tempi

           !..ion section:
           xni     = avo_eos * ytot1 * den
           dxnidd  = avo_eos * ytot1
           dxnida  = -xni * ytot1

           pion    = xni * kt
           dpiondd = dxnidd * kt
           dpiondt = xni * kerg

           eion    = 1.5_dp * pion*deni
           deiondd = (1.5_dp * dpiondd - eion)*deni
           deiondt = 1.5_dp * dpiondt*deni

           x       = abar*abar*sqrt(abar) * deni/avo_eos
           s       = sioncon * temp
           z       = x * s * sqrt(s)
           y       = log(z)
           sion    = (pion*deni + eion)*tempi + kergavo * ytot1 * y
           dsiondd = (dpiondd*deni - pion*deni*deni + deiondd)*tempi &
                - kergavo * deni * ytot1
           dsiondt = (dpiondt*deni + deiondt)*tempi -  &
                (pion*deni + eion) * tempi*tempi  &
                + 1.5_dp * kergavo * tempi*ytot1
           x       = avo_eos*kerg/abar

           !..electron-positron section:
           !..assume complete ionization
           xnem    = xni * zbar

           !..enter the table with ye*den
           din = ye*den

           !..hash locate this temperature and density
           jat = int((log10(temp) - tlo)*tstpi) + 1
           jat = max(1,min(jat,jtmax-1))
           iat = int((log10(din) - dlo)*dstpi) + 1
           iat = max(1,min(iat,itmax-1))

           !..access the table locations only once
           fi(1)  = f(iat,jat)
           fi(2)  = f(iat+1,jat)
           fi(3)  = f(iat,jat+1)
           fi(4)  = f(iat+1,jat+1)
           fi(5)  = ft(iat,jat)
           fi(6)  = ft(iat+1,jat)
           fi(7)  = ft(iat,jat+1)
           fi(8)  = ft(iat+1,jat+1)
           fi(9)  = ftt(iat,jat)
           fi(10) = ftt(iat+1,jat)
           fi(11) = ftt(iat,jat+1)
           fi(12) = ftt(iat+1,jat+1)
           fi(13) = fd(iat,jat)
           fi(14) = fd(iat+1,jat)
           fi(15) = fd(iat,jat+1)
           fi(16) = fd(iat+1,jat+1)
           fi(17) = fdd(iat,jat)
           fi(18) = fdd(iat+1,jat)
           fi(19) = fdd(iat,jat+1)
           fi(20) = fdd(iat+1,jat+1)
           fi(21) = fdt(iat,jat)
           fi(22) = fdt(iat+1,jat)
           fi(23) = fdt(iat,jat+1)
           fi(24) = fdt(iat+1,jat+1)
           fi(25) = fddt(iat,jat)
           fi(26) = fddt(iat+1,jat)
           fi(27) = fddt(iat,jat+1)
           fi(28) = fddt(iat+1,jat+1)
           fi(29) = fdtt(iat,jat)
           fi(30) = fdtt(iat+1,jat)
           fi(31) = fdtt(iat,jat+1)
           fi(32) = fdtt(iat+1,jat+1)
           fi(33) = fddtt(iat,jat)
           fi(34) = fddtt(iat+1,jat)
           fi(35) = fddtt(iat,jat+1)
           fi(36) = fddtt(iat+1,jat+1)

           !..various differences
           xt  = max( (temp - t(jat))*dti_sav(jat), 0.0_dp)
           xd  = max( (din - d(iat))*ddi_sav(iat), 0.0_dp)
           mxt = 1.0_dp - xt
           mxd = 1.0_dp - xd

           !..the six density and six temperature basis functions
           si0t =   psi0(xt)
           si1t =   psi1(xt)*dt_sav(jat)
           si2t =   psi2(xt)*dt2_sav(jat)

           si0mt =  psi0(mxt)
           si1mt = -psi1(mxt)*dt_sav(jat)
           si2mt =  psi2(mxt)*dt2_sav(jat)

           si0d =   psi0(xd)
           si1d =   psi1(xd)*dd_sav(iat)
           si2d =   psi2(xd)*dd2_sav(iat)

           si0md =  psi0(mxd)
           si1md = -psi1(mxd)*dd_sav(iat)
           si2md =  psi2(mxd)*dd2_sav(iat)

           !..derivatives of the weight functions
           dsi0t =   dpsi0(xt)*dti_sav(jat)
           dsi1t =   dpsi1(xt)
           dsi2t =   dpsi2(xt)*dt_sav(jat)

           dsi0mt = -dpsi0(mxt)*dti_sav(jat)
           dsi1mt =  dpsi1(mxt)
           dsi2mt = -dpsi2(mxt)*dt_sav(jat)

           dsi0d =   dpsi0(xd)*ddi_sav(iat)
           dsi1d =   dpsi1(xd)
           dsi2d =   dpsi2(xd)*dd_sav(iat)

           dsi0md = -dpsi0(mxd)*ddi_sav(iat)
           dsi1md =  dpsi1(mxd)
           dsi2md = -dpsi2(mxd)*dd_sav(iat)

           !..second derivatives of the weight functions
           ddsi0t =   ddpsi0(xt)*dt2i_sav(jat)
           ddsi1t =   ddpsi1(xt)*dti_sav(jat)
           ddsi2t =   ddpsi2(xt)

           ddsi0mt =  ddpsi0(mxt)*dt2i_sav(jat)
           ddsi1mt = -ddpsi1(mxt)*dti_sav(jat)
           ddsi2mt =  ddpsi2(mxt)

           !     ddsi0d =   ddpsi0(xd)*dd2i_sav(iat)
           !     ddsi1d =   ddpsi1(xd)*ddi_sav(iat)
           !     ddsi2d =   ddpsi2(xd)

           !     ddsi0md =  ddpsi0(mxd)*dd2i_sav(iat)
           !     ddsi1md = -ddpsi1(mxd)*ddi_sav(iat)
           !     ddsi2md =  ddpsi2(mxd)


           !..the free energy
           free  = h5( fi, &
                si0t,   si1t,   si2t,   si0mt,   si1mt,   si2mt, &
                si0d,   si1d,   si2d,   si0md,   si1md,   si2md)

           !..derivative with respect to density
           df_d  = h5( fi, &
                si0t,   si1t,   si2t,   si0mt,   si1mt,   si2mt, &
                dsi0d,  dsi1d,  dsi2d,  dsi0md,  dsi1md,  dsi2md)

           !..derivative with respect to temperature
           df_t = h5( fi, &
                dsi0t,  dsi1t,  dsi2t,  dsi0mt,  dsi1mt,  dsi2mt, &
                si0d,   si1d,   si2d,   si0md,   si1md,   si2md)

           !..derivative with respect to density**2
           !     df_dd = h5( &
           !               si0t,   si1t,   si2t,   si0mt,   si1mt,   si2mt, &
           !               ddsi0d, ddsi1d, ddsi2d, ddsi0md, ddsi1md, ddsi2md)

           !..derivative with respect to temperature**2
           df_tt = h5( fi, &
                ddsi0t, ddsi1t, ddsi2t, ddsi0mt, ddsi1mt, ddsi2mt, &
                si0d,   si1d,   si2d,   si0md,   si1md,   si2md)

           !..derivative with respect to temperature and density
           df_dt = h5( fi, &
                dsi0t,  dsi1t,  dsi2t,  dsi0mt,  dsi1mt,  dsi2mt, &
                dsi0d,  dsi1d,  dsi2d,  dsi0md,  dsi1md,  dsi2md)

           !..now get the pressure derivative with density, chemical potential, and
           !..electron positron number densities
           !..get the interpolation weight functions
           si0t   =  xpsi0(xt)
           si1t   =  xpsi1(xt)*dt_sav(jat)

           si0mt  =  xpsi0(mxt)
           si1mt  =  -xpsi1(mxt)*dt_sav(jat)

           si0d   =  xpsi0(xd)
           si1d   =  xpsi1(xd)*dd_sav(iat)

           si0md  =  xpsi0(mxd)
           si1md  =  -xpsi1(mxd)*dd_sav(iat)

           !..derivatives of weight functions
           dsi0t  = xdpsi0(xt)*dti_sav(jat)
           dsi1t  = xdpsi1(xt)

           dsi0mt = -xdpsi0(mxt)*dti_sav(jat)
           dsi1mt = xdpsi1(mxt)

           dsi0d  = xdpsi0(xd)*ddi_sav(iat)
           dsi1d  = xdpsi1(xd)

           dsi0md = -xdpsi0(mxd)*ddi_sav(iat)
           dsi1md = xdpsi1(mxd)

           !..look in the pressure derivative only once
           fi(1)  = dpdf(iat,jat)
           fi(2)  = dpdf(iat+1,jat)
           fi(3)  = dpdf(iat,jat+1)
           fi(4)  = dpdf(iat+1,jat+1)
           fi(5)  = dpdft(iat,jat)
           fi(6)  = dpdft(iat+1,jat)
           fi(7)  = dpdft(iat,jat+1)
           fi(8)  = dpdft(iat+1,jat+1)
           fi(9)  = dpdfd(iat,jat)
           fi(10) = dpdfd(iat+1,jat)
           fi(11) = dpdfd(iat,jat+1)
           fi(12) = dpdfd(iat+1,jat+1)
           fi(13) = dpdfdt(iat,jat)
           fi(14) = dpdfdt(iat+1,jat)
           fi(15) = dpdfdt(iat,jat+1)
           fi(16) = dpdfdt(iat+1,jat+1)

           !..pressure derivative with density
           dpepdd  = h3(   fi, &
                si0t,   si1t,   si0mt,   si1mt, &
                si0d,   si1d,   si0md,   si1md)
           dpepdd  = max(ye * dpepdd,0.0_dp)

           !..look in the electron chemical potential table only once
           fi(1)  = ef(iat,jat)
           fi(2)  = ef(iat+1,jat)
           fi(3)  = ef(iat,jat+1)
           fi(4)  = ef(iat+1,jat+1)
           fi(5)  = eft(iat,jat)
           fi(6)  = eft(iat+1,jat)
           fi(7)  = eft(iat,jat+1)
           fi(8)  = eft(iat+1,jat+1)
           fi(9)  = efd(iat,jat)
           fi(10) = efd(iat+1,jat)
           fi(11) = efd(iat,jat+1)
           fi(12) = efd(iat+1,jat+1)
           fi(13) = efdt(iat,jat)
           fi(14) = efdt(iat+1,jat)
           fi(15) = efdt(iat,jat+1)
           fi(16) = efdt(iat+1,jat+1)

           !..electron chemical potential etaele
           etaele  = h3( fi, &
                si0t,   si1t,   si0mt,   si1mt, &
                si0d,   si1d,   si0md,   si1md)

           !..derivative with respect to density
           x       = h3( fi, &
                si0t,   si1t,   si0mt,   si1mt, &
                dsi0d,  dsi1d,  dsi0md,  dsi1md)
           detadd  = ye * x

           !..derivative with respect to temperature
           detadt  = h3( fi, &
                dsi0t,  dsi1t,  dsi0mt,  dsi1mt, &
                si0d,   si1d,   si0md,   si1md)


           !..look in the number density table only once
           fi(1)  = xf(iat,jat)
           fi(2)  = xf(iat+1,jat)
           fi(3)  = xf(iat,jat+1)
           fi(4)  = xf(iat+1,jat+1)
           fi(5)  = xft(iat,jat)
           fi(6)  = xft(iat+1,jat)
           fi(7)  = xft(iat,jat+1)
           fi(8)  = xft(iat+1,jat+1)
           fi(9)  = xfd(iat,jat)
           fi(10) = xfd(iat+1,jat)
           fi(11) = xfd(iat,jat+1)
           fi(12) = xfd(iat+1,jat+1)
           fi(13) = xfdt(iat,jat)
           fi(14) = xfdt(iat+1,jat)
           fi(15) = xfdt(iat,jat+1)
           fi(16) = xfdt(iat+1,jat+1)

           !..electron + positron number densities
           xnefer   = h3( fi, &
                si0t,   si1t,   si0mt,   si1mt, &
                si0d,   si1d,   si0md,   si1md)

           !..derivative with respect to density
           x        = h3( fi, &
                si0t,   si1t,   si0mt,   si1mt, &
                dsi0d,  dsi1d,  dsi0md,  dsi1md)
           x = max(x,0.0_dp)
           dxnedd   = ye * x

           !..derivative with respect to temperature
           dxnedt   = h3( fi, &
                dsi0t,  dsi1t,  dsi0mt,  dsi1mt, &
                si0d,   si1d,   si0md,   si1md)


           !..the desired electron-positron thermodynamic quantities

           !..dpepdd at high temperatures and low densities is below the
           !..floating point limit of the subtraction of two large terms.
           !..since dpresdd doesn't enter the maxwell relations at all, use the
           !..bicubic interpolation done above instead of this one
           x       = din * din
           pele    = x * df_d
           dpepdt  = x * df_dt
           !     dpepdd  = ye * (x * df_dd + 2.0_dp * din * df_d)
           s       = dpepdd/ye - 2.0_dp * din * df_d

           x       = ye * ye
           sele    = -df_t * ye
           dsepdt  = -df_tt * ye
           dsepdd  = -df_dt * x

           eele    = ye*free + temp * sele
           deepdt  = temp * dsepdt
           deepdd  = x * df_d + temp * dsepdd

           !..coulomb section:
           !..initialize
           pcoul    = 0.0_dp
           dpcouldd = 0.0_dp
           dpcouldt = 0.0_dp
           dpcoulda = 0.0_dp
           dpcouldz = 0.0_dp
           ecoul    = 0.0_dp
           decouldd = 0.0_dp
           decouldt = 0.0_dp
           decoulda = 0.0_dp
           decouldz = 0.0_dp
           scoul    = 0.0_dp
           dscouldd = 0.0_dp
           dscouldt = 0.0_dp
           dscoulda = 0.0_dp
           dscouldz = 0.0_dp

           !..uniform background corrections only
           !..from yakovlev & shalybkov 1989
           !..lami is the average ion seperation
           !..plasg is the plasma coupling parameter
           z        = forth * pi
           s        = z * xni
           dsdd     = z * dxnidd
           dsda     = z * dxnida

           lami     = 1.0_dp/s**onethird
           inv_lami = 1.0_dp/lami
           z        = -onethird * lami
           lamidd   = z * dsdd/s
           lamida   = z * dsda/s

           plasg    = zbar*zbar*esqu*ktinv*inv_lami
           z        = -plasg * inv_lami
           plasgdd  = z * lamidd
           plasgda  = z * lamida
           plasgdt  = -plasg*ktinv * kerg
           plasgdz  = 2.0_dp * plasg/zbar

           !     TURN ON/OFF COULOMB
           if ( do_coulomb ) then
              !...yakovlev & shalybkov 1989 equations 82, 85, 86, 87
              if (plasg .ge. 1.0_dp) then
                 x        = plasg**(0.25_dp)
                 y        = avo_eos * ytot1 * kerg
                 ecoul    = y * temp * (a1*plasg + b1*x + c1/x + d1)
                 pcoul    = onethird * den * ecoul
                 scoul    = -y * (3.0_dp*b1*x - 5.0_dp*c1/x &
                      + d1 * (log(plasg) - 1.0_dp) - e1)

                 y        = avo_eos*ytot1*kt*(a1 + 0.25_dp/plasg*(b1*x - c1/x))
                 decouldd = y * plasgdd
                 decouldt = y * plasgdt + ecoul/temp
                 decoulda = y * plasgda - ecoul/abar
                 decouldz = y * plasgdz

                 y        = onethird * den
                 dpcouldd = onethird * ecoul + y*decouldd
                 dpcouldt = y * decouldt
                 dpcoulda = y * decoulda
                 dpcouldz = y * decouldz

                 y        = -avo_eos*kerg/(abar*plasg)* &
                      (0.75_dp*b1*x+1.25_dp*c1/x+d1)
                 dscouldd = y * plasgdd
                 dscouldt = y * plasgdt
                 dscoulda = y * plasgda - scoul/abar
                 dscouldz = y * plasgdz

                 !...yakovlev & shalybkov 1989 equations 102, 103, 104
              else if (plasg .lt. 1.0_dp) then
                 x        = plasg*sqrt(plasg)
                 y        = plasg**b2
                 z        = c2 * x - onethird * a2 * y
                 pcoul    = -pion * z
                 ecoul    = 3.0_dp * pcoul/den
                 scoul    = -avo_eos/abar*kerg*(c2*x -a2*(b2-1.0_dp)/b2*y)

                 s        = 1.5_dp*c2*x/plasg - onethird*a2*b2*y/plasg
                 dpcouldd = -dpiondd*z - pion*s*plasgdd
                 dpcouldt = -dpiondt*z - pion*s*plasgdt

                 s        = 3.0_dp/den
                 decouldd = s * dpcouldd - ecoul/den
                 decouldt = s * dpcouldt
                 decoulda = s * dpcoulda
                 decouldz = s * dpcouldz

                 s        = -avo_eos*kerg/(abar*plasg)* &
                      (1.5_dp*c2*x-a2*(b2-1.0_dp)*y)
                 dscouldd = s * plasgdd
                 dscouldt = s * plasgdt
                 dscoulda = s * plasgda - scoul/abar
                 dscouldz = s * plasgdz
              end if

              ! Disable Coulomb corrections if they cause
              ! the energy or pressure to go negative.

              p_temp = prad + pion + pele + pcoul
              e_temp = erad + eion + eele + ecoul

              if (p_temp .le. 0.0_dp .or. e_temp .le. 0.0_dp) then

                 pcoul    = 0.0_dp
                 dpcouldd = 0.0_dp
                 dpcouldt = 0.0_dp
                 dpcoulda = 0.0_dp
                 dpcouldz = 0.0_dp
                 ecoul    = 0.0_dp
                 decouldd = 0.0_dp
                 decouldt = 0.0_dp
                 decoulda = 0.0_dp
                 decouldz = 0.0_dp
                 scoul    = 0.0_dp
                 dscouldd = 0.0_dp
                 dscouldt = 0.0_dp
                 dscoulda = 0.0_dp
                 dscouldz = 0.0_dp

              end if
           end if

           !..sum all the components
           pres    = prad + pion + pele + pcoul
           ener    = erad + eion + eele + ecoul
           entr    = srad + sion + sele + scoul

           dpresdd = dpraddd + dpiondd + dpepdd + dpcouldd
           dpresdt = dpraddt + dpiondt + dpepdt + dpcouldt
           denerdd = deraddd + deiondd + deepdd + decouldd
           denerdt = deraddt + deiondt + deepdt + decouldt

           dentrdd = dsraddd + dsiondd + dsepdd + dscouldd
           dentrdt = dsraddt + dsiondt + dsepdt + dscouldt

           !..the temperature and density exponents (c&g 9.81 9.82)
           !..the specific heat at constant volume (c&g 9.92)
           !..the third adiabatic exponent (c&g 9.93)
           !..the first adiabatic exponent (c&g 9.97)
           !..the second adiabatic exponent (c&g 9.105)
           !..the specific heat at constant pressure (c&g 9.98)
           !..and relativistic formula for the sound speed (c&g 14.29)
           zz    = pres*deni
           zzi   = den/pres
           chit  = temp/pres * dpresdt
           chid  = dpresdd*zzi
           cv    = denerdt
           x     = zz * chit/(temp * cv)
           gam3  = x + 1.0_dp
           gam1  = chit*x + chid
           nabad = x/gam1
           gam2  = 1.0_dp/(1.0_dp - nabad)
           cp    = cv * gam1/chid
           z     = 1.0_dp + (ener + light2)*zzi
           sound = clight * sqrt(gam1/z)

           !..maxwell relations; each is zero if the consistency is perfect
           x   = den * den
           dse = temp*dentrdt/denerdt - 1.0_dp
           dpe = (denerdd*x + temp*dpresdt)/pres - 1.0_dp
           dsp = -dentrdd*x/dpresdt - 1.0_dp

           ptot_row = pres
           dpt_row = dpresdt
           dpd_row = dpresdd
           dpe_row = dpresdt / denerdt
           dpdr_e_row = dpresdd - dpresdt * denerdd / denerdt

           etot_row = ener
           det_row = denerdt
           ded_row = denerdd

           stot_row = entr
           dst_row = dentrdt
           dsd_row = dentrdd

           htot_row = ener + pres / den
           dhd_row = denerdd + dpresdd / den - pres / den**2
           dht_row = denerdt + dpresdt / den

           pele_row = pele
           ppos_row = 0.0_dp

           xne_row = xnefer
           xnp_row = 0.0_dp

           etaele_row = etaele
           detadt_row = detadt

           cv_row = cv
           cp_row = cp
           cs_row = sound
           gam1_row = gam1

           if (converged) then

              exit

           elseif (single_iter) then

              if (dvar .eq. itemp) then

                 x = temp_row
                 smallx = smallt
                 xtol = ttol

                 if (var .eq. ipres) then
                    v    = ptot_row
                    dvdx = dpt_row
                 elseif (var .eq. iener) then
                    v    = etot_row
                    dvdx = det_row
                 elseif (var .eq. ientr) then
                    v    = stot_row
                    dvdx = dst_row
                 elseif (var .eq. ienth) then
                    v    = htot_row
                    dvdx = dht_row
                 else
                    exit
                 endif

              else ! dvar == density

                 x = den_row
                 smallx = smalld
                 xtol = dtol

                 if (var .eq. ipres) then
                    v    = ptot_row
                    dvdx = dpd_row
                 elseif (var .eq. iener) then
                    v    = etot_row
                    dvdx = ded_row
                 elseif (var .eq. ientr) then
                    v    = stot_row
                    dvdx = dsd_row
                 elseif (var .eq. ienth) then
                    v    = htot_row
                    dvdx = dhd_row
                 else
                    exit
                 endif

              endif

              ! Now do the calculation for the next guess for T/rho

              xnew = x - (v - v_want) / dvdx

              ! Don't let the temperature/density change by more than a factor of two
              xnew = max(0.5 * x, min(xnew, 2.0 * x))

              ! Don't let us freeze/evacuate
              xnew = max(smallx, xnew)

              ! Store the new temperature/density

              if (dvar .eq. itemp) then
                 temp_row = xnew
              else
                 den_row  = xnew
              endif

              ! Compute the error from the last iteration

              error = abs( (xnew - x) / x )

              if (error .lt. xtol) converged = .true.

           elseif (double_iter) then

              ! Figure out which variables we're using

              told = temp_row
              rold = den_row

              if (var1 .eq. ipres) then
                 v1    = ptot_row
                 dv1dt = dpt_row
                 dv1dr = dpd_row
              elseif (var1 .eq. iener) then
                 v1    = etot_row
                 dv1dt = det_row
                 dv1dr = ded_row
              elseif (var1 .eq. ientr) then
                 v1    = stot_row
                 dv1dt = dst_row
                 dv1dr = dsd_row
              elseif (var1 .eq. ienth) then
                 v1    = htot_row
                 dv1dt = dht_row
                 dv1dr = dhd_row
              else
                 exit
              endif

              if (var2 .eq. ipres) then
                 v2    = ptot_row
                 dv2dt = dpt_row
                 dv2dr = dpd_row
              elseif (var2 .eq. iener) then
                 v2    = etot_row
                 dv2dt = det_row
                 dv2dr = ded_row
              elseif (var2 .eq. ientr) then
                 v2    = stot_row
                 dv2dt = dst_row
                 dv2dr = dsd_row
              elseif (var2 .eq. ienth) then
                 v2    = htot_row
                 dv2dt = dht_row
                 dv2dr = dhd_row
              else
                 exit
              endif

              ! Two functions, f and g, to iterate over
              v1i = v1_want - v1
              v2i = v2_want - v2

              !
              ! 0 = f + dfdr * delr + dfdt * delt
              ! 0 = g + dgdr * delr + dgdt * delt
              !

              ! note that dfi/dT = - df/dT
              delr = (-v1i*dv2dt + v2i*dv1dt) / (dv2dr*dv1dt - dv2dt*dv1dr)

              rnew = rold + delr

              tnew = told + (v1i - dv1dr*delr) / dv1dt

              ! Don't let the temperature or density change by more
              ! than a factor of two
              tnew = max(0.5_dp * told, min(tnew, 2.0_dp * told))
              rnew = max(0.5_dp * rold, min(rnew, 2.0_dp * rold))

              ! Don't let us freeze or evacuate
              tnew = max(smallt, tnew)
              rnew = max(smalld, rnew)

              ! Store the new temperature and density
              den_row  = rnew
              temp_row = tnew

              ! Compute the errors
              error1 = abs( (rnew - rold) / rold )
              error2 = abs( (tnew - told) / told )

              if (error1 .LT. dtol .and. error2 .LT. ttol) converged = .true.

           endif

        enddo

        state % T    = temp_row
        state % rho  = den_row

        state % p    = ptot_row
        state % dpdT = dpt_row
        state % dpdr = dpd_row


        state % dpde = dpe_row
        state % dpdr_e = dpdr_e_row

        state % e    = etot_row
        state % dedT = det_row
        state % dedr = ded_row

        state % s    = stot_row
        state % dsdT = dst_row
        state % dsdr = dsd_row

        state % h    = htot_row
        state % dhdR = dhd_row
        state % dhdT = dht_row

        state % pele = pele_row
        state % ppos = ppos_row

        state % xne = xne_row
        state % xnp = xnp_row

        state % eta = etaele_row
        state % detadt = detadt_row

        state % cv   = cv_row
        state % cp   = cp_row
        state % gam1 = gam1_row
        ! state % cs   = cs_row

        ! Take care of final housekeeping.

        ! Count the positron contribution in the electron quantities.

        state % xne  = state % xne  + state % xnp
        state % pele = state % pele + state % ppos

        ! Use the non-relativistic version of the sound speed, cs = sqrt(gam_1 * P / rho).
        ! This replaces the relativistic version that comes out of helmeos.

        state % cs = sqrt(state % gam1 * state % p / state % rho)

        if (input_is_constant) then

          if (input .eq. eos_input_rh) then

            state % h = v_want

          elseif (input .eq. eos_input_tp) then

            state % p = v_want

          elseif (input .eq. eos_input_rp) then

            state % p = v_want

          elseif (input .eq. eos_input_re) then

            state % e = v_want

          elseif (input .eq. eos_input_ps) then

            state % p = v1_want
            state % s = v2_want

          elseif (input .eq. eos_input_ph) then

            state % p = v1_want
            state % h = v2_want

          elseif (input .eq. eos_input_th) then

            state % h = v_want

          endif

        endif

    end subroutine actual_eos


    subroutine xnet_actual_eos(input, state)

        ! This is a pruned version of actual_eos that only calculates those 
        ! quantities needed by XNet: electron chemical potential, its derivative
        ! w.r.t. temperature, and specific heat.

        !$acc routine seq

        implicit none

        !..input arguments
        integer,      intent(in   ) :: input
        type (eos_t), intent(inout) :: state

        !..declare local variables
        real(dp) :: cv,etaele,detadt, &
                    temp,den,abar,zbar,ye

        temp  = state % T
        den   = state % rho
        abar  = state % abar
        zbar  = state % zbar

        ye    = state % y_e

        call actual_eos_cv(temp,den,abar,zbar,ye,cv)
        call actual_eos_eta(temp,den,ye,etaele,detadt)

        state % cv     = cv
        state % eta    = etaele
        state % detadt = detadt

    end subroutine xnet_actual_eos


    subroutine actual_eos_eta(temp,den,ye,etaele,detadt)

        ! This is a pruned version of actual_eos that only calculates those 
        ! quantities needed by XNet: electron chemical potential, its derivative
        ! w.r.t. temperature, and specific heat.

        !$acc routine seq

        implicit none

        !..input arguments
        real(dp),     intent(in ) :: temp,den,ye
        real(dp),     intent(out) :: etaele,detadt

        !..declare local variables
        real(dp) :: din

        !..for the interpolations
        integer  :: iat,jat
        real(dp) :: xt,xd,mxt,mxd,dfi(16), &
                    si0t,si1t,si2t,si0mt,si1mt,si2mt, &
                    si0d,si1d,si2d,si0md,si1md,si2md, &
                    dsi0t,dsi1t,dsi2t,dsi0mt,dsi1mt,dsi2mt

        din   = ye * den

        !..electron-positron section:
        !..hash locate this temperature and density
        jat = int((log10(temp) - tlo)*tstpi) + 1
        jat = max(1,min(jat,jtmax-1))
        iat = int((log10(din) - dlo)*dstpi) + 1
        iat = max(1,min(iat,itmax-1))

        !..various differences
        xt  = max( (temp - t(jat))*dti_sav(jat), 0.0_dp)
        xd  = max( (din - d(iat))*ddi_sav(iat), 0.0_dp)
        mxt = 1.0_dp - xt
        mxd = 1.0_dp - xd

        !..now get the pressure derivative with density, chemical potential, and
        !..electron positron number densities
        !..get the interpolation weight functions
        si0t   =  xpsi0(xt)
        si1t   =  xpsi1(xt)*dt_sav(jat)

        si0mt  =  xpsi0(mxt)
        si1mt  =  -xpsi1(mxt)*dt_sav(jat)

        si0d   =  xpsi0(xd)
        si1d   =  xpsi1(xd)*dd_sav(iat)

        si0md  =  xpsi0(mxd)
        si1md  =  -xpsi1(mxd)*dd_sav(iat)

        !..derivatives of weight functions
        dsi0t  = xdpsi0(xt)*dti_sav(jat)
        dsi1t  = xdpsi1(xt)

        dsi0mt = -xdpsi0(mxt)*dti_sav(jat)
        dsi1mt = xdpsi1(mxt)

        !..look in the electron chemical potential table only once
        dfi(1)  = ef(iat,jat)
        dfi(2)  = ef(iat+1,jat)
        dfi(3)  = ef(iat,jat+1)
        dfi(4)  = ef(iat+1,jat+1)
        dfi(5)  = eft(iat,jat)
        dfi(6)  = eft(iat+1,jat)
        dfi(7)  = eft(iat,jat+1)
        dfi(8)  = eft(iat+1,jat+1)
        dfi(9)  = efd(iat,jat)
        dfi(10) = efd(iat+1,jat)
        dfi(11) = efd(iat,jat+1)
        dfi(12) = efd(iat+1,jat+1)
        dfi(13) = efdt(iat,jat)
        dfi(14) = efdt(iat+1,jat)
        dfi(15) = efdt(iat,jat+1)
        dfi(16) = efdt(iat+1,jat+1)

        !..electron chemical potential etaele
        etaele  = h3( dfi, &
          si0t,   si1t,   si0mt,   si1mt, &
          si0d,   si1d,   si0md,   si1md)

        !..derivative with respect to temperature
        detadt  = h3( dfi, &
          dsi0t,  dsi1t,  dsi0mt,  dsi1mt, &
          si0d,   si1d,   si0md,   si1md)

    end subroutine actual_eos_eta


    subroutine actual_eos_cv(temp,den,abar,zbar,ye,cv)

        ! This is a pruned version of actual_eos that only calculates those 
        ! quantities needed by XNet: electron chemical potential, its derivative
        ! w.r.t. temperature, and specific heat.

        implicit none

        !$acc routine seq

        !..input arguments
        real(dp),     intent(in ) :: temp,den,abar,zbar,ye
        real(dp),     intent(out) :: cv

        !..declare local variables
        real(dp) :: x,y,z,deni,tempi,xni, &
                    deepdt,dsepdt, &
                    dpraddt,deraddt,dpiondt, &
                    deiondt, &
                    kt,ktinv,prad,erad,pion,eion, &
                    pele,eele,sele, &
                    s, &
                    ytot1,din


        !..for the interpolations
        integer  :: iat,jat
        real(dp) :: free,df_d,df_t,df_tt
        real(dp) :: xt,xd,mxt,mxd,fi(36),dfi(16), &
                    si0t,si1t,si2t,si0mt,si1mt,si2mt, &
                    si0d,si1d,si2d,si0md,si1md,si2md, &
                    dsi0t,dsi1t,dsi2t,dsi0mt,dsi1mt,dsi2mt, &
                    dsi0d,dsi1d,dsi2d,dsi0md,dsi1md,dsi2md, &
                    ddsi0t,ddsi1t,ddsi2t,ddsi0mt,ddsi1mt,ddsi2mt

        !..for the coulomb corrections
        real(dp) :: lami,inv_lami, &
                    plasg,plasgdt, &
                    ecoul,decouldt, &
                    pcoul,dpcouldt

        real(dp) :: p_temp, e_temp

        ytot1 = 1.0_dp/abar
        din   = ye * den

        !..initialize
        deni    = 1.0_dp/den
        tempi   = 1.0_dp/temp
        kt      = kerg * temp
        ktinv   = 1.0_dp/kt

        !..radiation section:
        prad    = asoli3 * temp * temp * temp * temp
        dpraddt = 4.0_dp * prad*tempi
        deraddt = 3.0_dp * dpraddt*deni

        !..ion section:
        xni     = avo_eos * ytot1 * den
        pion    = xni * kt
        dpiondt = xni * kerg
        deiondt = 1.5_dp * dpiondt*deni

        !..electron-positron section:
        !..hash locate this temperature and density
        jat = int((log10(temp) - tlo)*tstpi) + 1
        jat = max(1,min(jat,jtmax-1))
        iat = int((log10(din) - dlo)*dstpi) + 1
        iat = max(1,min(iat,itmax-1))

        !..various differences
        xt  = max( (temp - t(jat))*dti_sav(jat), 0.0_dp)
        xd  = max( (din - d(iat))*ddi_sav(iat), 0.0_dp)
        mxt = 1.0_dp - xt
        mxd = 1.0_dp - xd

        !..the six density and six temperature basis functions
        si0t =   psi0(xt)
        si1t =   psi1(xt)*dt_sav(jat)
        si2t =   psi2(xt)*dt2_sav(jat)

        si0mt =  psi0(mxt)
        si1mt = -psi1(mxt)*dt_sav(jat)
        si2mt =  psi2(mxt)*dt2_sav(jat)

        si0d =   psi0(xd)
        si1d =   psi1(xd)*dd_sav(iat)
        si2d =   psi2(xd)*dd2_sav(iat)

        si0md =  psi0(mxd)
        si1md = -psi1(mxd)*dd_sav(iat)
        si2md =  psi2(mxd)*dd2_sav(iat)

        !..derivatives of the weight functions
        !dsi0t =   dpsi0(xt)*dti_sav(jat)
        !dsi1t =   dpsi1(xt)
        !dsi2t =   dpsi2(xt)*dt_sav(jat)

        !dsi0mt = -dpsi0(mxt)*dti_sav(jat)
        !dsi1mt =  dpsi1(mxt)
        !dsi2mt = -dpsi2(mxt)*dt_sav(jat)

        dsi0d =   dpsi0(xd)*ddi_sav(iat)
        dsi1d =   dpsi1(xd)
        dsi2d =   dpsi2(xd)*dd_sav(iat)

        dsi0md = -dpsi0(mxd)*ddi_sav(iat)
        dsi1md =  dpsi1(mxd)
        dsi2md = -dpsi2(mxd)*dd_sav(iat)

        !..second derivatives of the weight functions
        ddsi0t =   ddpsi0(xt)*dt2i_sav(jat)
        ddsi1t =   ddpsi1(xt)*dti_sav(jat)
        ddsi2t =   ddpsi2(xt)

        ddsi0mt =  ddpsi0(mxt)*dt2i_sav(jat)
        ddsi1mt = -ddpsi1(mxt)*dti_sav(jat)
        ddsi2mt =  ddpsi2(mxt)

        !..access the table locations only once
        fi(1)  = f(iat,jat)
        fi(2)  = f(iat+1,jat)
        fi(3)  = f(iat,jat+1)
        fi(4)  = f(iat+1,jat+1)
        fi(5)  = ft(iat,jat)
        fi(6)  = ft(iat+1,jat)
        fi(7)  = ft(iat,jat+1)
        fi(8)  = ft(iat+1,jat+1)
        fi(9)  = ftt(iat,jat)
        fi(10) = ftt(iat+1,jat)
        fi(11) = ftt(iat,jat+1)
        fi(12) = ftt(iat+1,jat+1)
        fi(13) = fd(iat,jat)
        fi(14) = fd(iat+1,jat)
        fi(15) = fd(iat,jat+1)
        fi(16) = fd(iat+1,jat+1)
        fi(17) = fdd(iat,jat)
        fi(18) = fdd(iat+1,jat)
        fi(19) = fdd(iat,jat+1)
        fi(20) = fdd(iat+1,jat+1)
        fi(21) = fdt(iat,jat)
        fi(22) = fdt(iat+1,jat)
        fi(23) = fdt(iat,jat+1)
        fi(24) = fdt(iat+1,jat+1)
        fi(25) = fddt(iat,jat)
        fi(26) = fddt(iat+1,jat)
        fi(27) = fddt(iat,jat+1)
        fi(28) = fddt(iat+1,jat+1)
        fi(29) = fdtt(iat,jat)
        fi(30) = fdtt(iat+1,jat)
        fi(31) = fdtt(iat,jat+1)
        fi(32) = fdtt(iat+1,jat+1)
        fi(33) = fddtt(iat,jat)
        fi(34) = fddtt(iat+1,jat)
        fi(35) = fddtt(iat,jat+1)
        fi(36) = fddtt(iat+1,jat+1)

        !..the free energy
        !free  = h5( fi, &
        !   si0t,   si1t,   si2t,   si0mt,   si1mt,   si2mt, &
        !   si0d,   si1d,   si2d,   si0md,   si1md,   si2md)

        !..derivative with respect to density
        df_d  = h5( fi, &
           si0t,   si1t,   si2t,   si0mt,   si1mt,   si2mt, &
           dsi0d,  dsi1d,  dsi2d,  dsi0md,  dsi1md,  dsi2md)

        !..derivative with respect to temperature
        !df_t = h5( fi, &
        !   dsi0t,  dsi1t,  dsi2t,  dsi0mt,  dsi1mt,  dsi2mt, &
        !   si0d,   si1d,   si2d,   si0md,   si1md,   si2md)

        !..derivative with respect to temperature**2
        df_tt = h5( fi, &
           ddsi0t, ddsi1t, ddsi2t, ddsi0mt, ddsi1mt, ddsi2mt, &
           si0d,   si1d,   si2d,   si0md,   si1md,   si2md)

        !..the desired electron-positron thermodynamic quantities
        x       = din * din
        pele    = x * df_d

        !sele    = -df_t * ye
        dsepdt  = -df_tt * ye

        !eele    = ye*free + temp * sele
        deepdt  = temp * dsepdt

        !..coulomb section:
        !..uniform background corrections only
        !..from yakovlev & shalybkov 1989
        !..plasg is the plasma coupling parameter
        z        = forth * pi
        s        = z * xni

        lami     = 1.0_dp/s**onethird
        inv_lami = 1.0_dp/lami

        plasg    = zbar*zbar*esqu*ktinv*inv_lami
        plasgdt  = -plasg*ktinv * kerg

        !     TURN ON/OFF COULOMB
        !...yakovlev & shalybkov 1989 equations 82, 85, 86, 87
        if (plasg .ge. 1.0_dp) then
          x        = plasg**(0.25_dp)
          y        = avo_eos * ytot1 * kerg
          ecoul    = y * temp * (a1*plasg + b1*x + c1/x + d1)
          pcoul    = onethird * den * ecoul

          y        = avo_eos*ytot1*kt*(a1 + 0.25_dp/plasg*(b1*x - c1/x))
          decouldt = y * plasgdt + ecoul/temp

        !...yakovlev & shalybkov 1989 equations 102, 103, 104
        else if (plasg .lt. 1.0_dp) then
          x        = plasg*sqrt(plasg)
          y        = plasg**b2
          z        = c2 * x - onethird * a2 * y
          pcoul    = -pion * z
          ecoul    = 3.0_dp * pcoul/den

          s        = 1.5_dp*c2*x/plasg - onethird*a2*b2*y/plasg
          dpcouldt = -dpiondt*z - pion*s*plasgdt

          s        = 3.0_dp/den
          decouldt = s * dpcouldt
        end if

        ! Disable Coulomb corrections if they cause
        ! the energy or pressure to go negative.

        p_temp = prad + pion + pele + pcoul
        e_temp = 0.0_dp
        !e_temp = erad + eion + eele + ecoul

        if (p_temp .le. 0.0_dp ) then
          decouldt = 0.0_dp
        end if

        !..the specific heat at constant volume (c&g 9.92)
        cv = deraddt + deiondt + deepdt + decouldt

    end subroutine actual_eos_cv


    subroutine actual_eos_init

        use xnet_parallel, only: parallel_bcast, parallel_IOProcessor
        use xnet_util, only: xnet_terminate

        implicit none

        real(dp) :: dth, dt2, dti, dt2i
        real(dp) :: dd, dd2, ddi, dd2i
        real(dp) :: tsav, dsav
        integer :: i, j
        integer :: status

        ! Allocate managed module variables

        allocate(do_coulomb)
        allocate(input_is_constant)
        allocate(itmax)
        allocate(jtmax)
        allocate(d(imax))
        allocate(t(jmax))
        allocate(tlo)
        allocate(thi)
        allocate(tstp)
        allocate(tstpi)
        allocate(dlo)
        allocate(dhi)
        allocate(dstp)
        allocate(dstpi)
        allocate(ttol)
        allocate(dtol)
        allocate(f(imax,jmax))
        allocate(fd(imax,jmax))
        allocate(ft(imax,jmax))
        allocate(fdd(imax,jmax))
        allocate(ftt(imax,jmax))
        allocate(fdt(imax,jmax))
        allocate(fddt(imax,jmax))
        allocate(fdtt(imax,jmax))
        allocate(fddtt(imax,jmax))
        allocate(dpdf(imax,jmax))
        allocate(dpdfd(imax,jmax))
        allocate(dpdft(imax,jmax))
        allocate(dpdfdt(imax,jmax))
        allocate(ef(imax,jmax))
        allocate(efd(imax,jmax))
        allocate(eft(imax,jmax))
        allocate(efdt(imax,jmax))
        allocate(xf(imax,jmax))
        allocate(xfd(imax,jmax))
        allocate(xft(imax,jmax))
        allocate(xfdt(imax,jmax))
        allocate(dt_sav(jmax))
        allocate(dt2_sav(jmax))
        allocate(dti_sav(jmax))
        allocate(dt2i_sav(jmax))
        allocate(dd_sav(imax))
        allocate(dd2_sav(imax))
        allocate(ddi_sav(imax))
        allocate(dd2i_sav(imax))

        ! Read in the runtime parameters

        input_is_constant = .true.
        do_coulomb = .true.
        ttol = 1.0e-8_dp
        dtol = 1.0e-8_dp

        !..   read the helmholtz free energy table
        itmax = imax
        jtmax = jmax
        tlo   = 3.0_dp
        thi   = 13.0_dp
        tstp  = (thi - tlo)/float(jmax-1)
        tstpi = 1.0_dp/tstp
        dlo   = -12.0_dp
        dhi   = 15.0_dp
        dstp  = (dhi - dlo)/float(imax-1)
        dstpi = 1.0_dp/dstp

        do j=1,jmax
           tsav = tlo + (j-1)*tstp
           t(j) = 10.0_dp**(tsav)
           do i=1,imax
              dsav = dlo + (i-1)*dstp
              d(i) = 10.0_dp**(dsav)
           end do
        end do

        if (parallel_IOProcessor()) then

           !..   open the table
           open(unit=2,file='helm_table.dat',status='old',iostat=status,action='read')
           if (status > 0) then

              call xnet_terminate('actual_eos_init: Failed to open helm_table.dat')

           endif

           !...  read in the free energy table
           do j=1,jmax
              do i=1,imax
                 read(2,*) f(i,j),fd(i,j),ft(i,j),fdd(i,j),ftt(i,j),fdt(i,j), &
                      fddt(i,j),fdtt(i,j),fddtt(i,j)
              end do
           end do

           !..   read the pressure derivative with density table
           do j = 1, jmax
              do i = 1, imax
                 read(2,*) dpdf(i,j),dpdfd(i,j),dpdft(i,j),dpdfdt(i,j)
              end do
           end do

           !..   read the electron chemical potential table
           do j = 1, jmax
              do i = 1, imax
                 read(2,*) ef(i,j),efd(i,j),eft(i,j),efdt(i,j)
              end do
           end do

           !..   read the number density table
           do j = 1, jmax
              do i = 1, imax
                 read(2,*) xf(i,j),xfd(i,j),xft(i,j),xfdt(i,j)
              end do
           end do
        end if

        call parallel_bcast(f)
        call parallel_bcast(fd)
        call parallel_bcast(ft)
        call parallel_bcast(fdd)
        call parallel_bcast(ftt)
        call parallel_bcast(fdt)
        call parallel_bcast(fddt)
        call parallel_bcast(fdtt)
        call parallel_bcast(fddtt)
        call parallel_bcast(dpdf)
        call parallel_bcast(dpdfd)
        call parallel_bcast(dpdft)
        call parallel_bcast(dpdfdt)
        call parallel_bcast(ef)
        call parallel_bcast(efd)
        call parallel_bcast(eft)
        call parallel_bcast(efdt)
        call parallel_bcast(xf)
        call parallel_bcast(xfd)
        call parallel_bcast(xft)
        call parallel_bcast(xfdt)

        !..   construct the temperature and density deltas and their inverses
        do j = 1, jmax-1
           dth         = t(j+1) - t(j)
           dt2         = dth * dth
           dti         = 1.0_dp/dth
           dt2i        = 1.0_dp/dt2
           dt_sav(j)   = dth
           dt2_sav(j)  = dt2
           dti_sav(j)  = dti
           dt2i_sav(j) = dt2i
        end do
        do i = 1, imax-1
           dd          = d(i+1) - d(i)
           dd2         = dd * dd
           ddi         = 1.0_dp/dd
           dd2i        = 1.0_dp/dd2
           dd_sav(i)   = dd
           dd2_sav(i)  = dd2
           ddi_sav(i)  = ddi
           dd2i_sav(i) = dd2i
        end do

        if (parallel_IOProcessor()) then
           close(unit=2)
        endif

        ! Set up the minimum and maximum possible densities.

        mintemp = 10.d0**tlo
        maxtemp = 10.d0**thi
        mindens = 10.d0**dlo
        maxdens = 10.d0**dhi

        !$acc update device(mintemp, maxtemp, mindens, maxdens)

        !$acc update &
        !$acc device(tlo, thi, dlo, dhi) &
        !$acc device(tstp, tstpi, dstp, dstpi) &
        !$acc device(itmax, jtmax, d, t) &
        !$acc device(f, fd, ft, fdd, ftt, fdt, fddt, fdtt, fddtt) &
        !$acc device(dpdf, dpdfd, dpdft, dpdfdt) &
        !$acc device(ef, efd, eft, efdt, xf, xfd, xft, xfdt)  &
        !$acc device(dt_sav, dt2_sav, dti_sav, dt2i_sav) &
        !$acc device(dd_sav, dd2_sav, ddi_sav, dd2i_sav) &
        !$acc device(do_coulomb, input_is_constant)

    end subroutine actual_eos_init



    ! quintic hermite polynomial functions
    ! psi0 and its derivatives
    function psi0(z) result(psi0r)
      !$acc routine seq
      real(dp), intent(in) :: z
      real(dp) :: psi0r
      psi0r = z**3 * ( z * (-6.0_dp*z + 15.0_dp) -10.0_dp) + 1.0_dp
    end function psi0

    function dpsi0(z) result(dpsi0r)
      !$acc routine seq
      real(dp), intent(in) :: z
      real(dp) :: dpsi0r
      dpsi0r = z**2 * ( z * (-30.0_dp*z + 60.0_dp) - 30.0_dp)
    end function dpsi0

    function ddpsi0(z) result(ddpsi0r)
      !$acc routine seq
      real(dp), intent(in) :: z
      real(dp) :: ddpsi0r
      ddpsi0r = z* ( z*( -120.0_dp*z + 180.0_dp) -60.0_dp)
    end function ddpsi0

    ! psi1 and its derivatives
    function psi1(z) result(psi1r)
      !$acc routine seq
      real(dp), intent(in) :: z
      real(dp) :: psi1r
      psi1r = z* ( z**2 * ( z * (-3.0_dp*z + 8.0_dp) - 6.0_dp) + 1.0_dp)
    end function psi1

    function dpsi1(z) result(dpsi1r)
      !$acc routine seq
      real(dp), intent(in) :: z
      real(dp) :: dpsi1r
      dpsi1r = z*z * ( z * (-15.0_dp*z + 32.0_dp) - 18.0_dp) +1.0_dp
    end function dpsi1

    function ddpsi1(z) result(ddpsi1r)
      !$acc routine seq
      real(dp), intent(in) :: z
      real(dp) :: ddpsi1r
      ddpsi1r = z * (z * (-60.0_dp*z + 96.0_dp) -36.0_dp)
    end function ddpsi1

    ! psi2  and its derivatives
    function psi2(z) result(psi2r)
      !$acc routine seq
      real(dp), intent(in) :: z
      real(dp) :: psi2r
      psi2r = 0.5_dp*z*z*( z* ( z * (-z + 3.0_dp) - 3.0_dp) + 1.0_dp)
    end function psi2

    function dpsi2(z) result(dpsi2r)
      !$acc routine seq
      real(dp), intent(in) :: z
      real(dp) :: dpsi2r
      dpsi2r = 0.5_dp*z*( z*(z*(-5.0_dp*z + 12.0_dp) - 9.0_dp) + 2.0_dp)
    end function dpsi2

    function ddpsi2(z) result(ddpsi2r)
      !$acc routine seq
      real(dp), intent(in) :: z
      real(dp) :: ddpsi2r
      ddpsi2r = 0.5_dp*(z*( z * (-20.0_dp*z + 36.0_dp) - 18.0_dp) + 2.0_dp)
    end function ddpsi2


    ! biquintic hermite polynomial function
    function h5(fi,w0t,w1t,w2t,w0mt,w1mt,w2mt,w0d,w1d,w2d,w0md,w1md,w2md) result(h5r)
      !$acc routine seq
      real(dp), intent(in) :: fi(36)
      real(dp), intent(in) :: w0t,w1t,w2t,w0mt,w1mt,w2mt,w0d,w1d,w2d,w0md,w1md,w2md
      real(dp) :: h5r

      h5r =  fi(1)  *w0d*w0t   + fi(2)  *w0md*w0t &
           + fi(3)  *w0d*w0mt  + fi(4)  *w0md*w0mt &
           + fi(5)  *w0d*w1t   + fi(6)  *w0md*w1t &
           + fi(7)  *w0d*w1mt  + fi(8)  *w0md*w1mt &
           + fi(9)  *w0d*w2t   + fi(10) *w0md*w2t &
           + fi(11) *w0d*w2mt  + fi(12) *w0md*w2mt &
           + fi(13) *w1d*w0t   + fi(14) *w1md*w0t &
           + fi(15) *w1d*w0mt  + fi(16) *w1md*w0mt &
           + fi(17) *w2d*w0t   + fi(18) *w2md*w0t &
           + fi(19) *w2d*w0mt  + fi(20) *w2md*w0mt &
           + fi(21) *w1d*w1t   + fi(22) *w1md*w1t &
           + fi(23) *w1d*w1mt  + fi(24) *w1md*w1mt &
           + fi(25) *w2d*w1t   + fi(26) *w2md*w1t &
           + fi(27) *w2d*w1mt  + fi(28) *w2md*w1mt &
           + fi(29) *w1d*w2t   + fi(30) *w1md*w2t &
           + fi(31) *w1d*w2mt  + fi(32) *w1md*w2mt &
           + fi(33) *w2d*w2t   + fi(34) *w2md*w2t &
           + fi(35) *w2d*w2mt  + fi(36) *w2md*w2mt
    end function h5


    ! cubic hermite polynomial functions
    ! psi0 & derivatives
    function xpsi0(z) result(xpsi0r)
      !$acc routine seq
      real(dp), intent(in) :: z
      real(dp) :: xpsi0r
      xpsi0r = z * z * (2.0_dp*z - 3.0_dp) + 1.0
    end function xpsi0

    function xdpsi0(z) result(xdpsi0r)
      !$acc routine seq
      real(dp), intent(in) :: z
      real(dp) :: xdpsi0r
      xdpsi0r = z * (6.0_dp*z - 6.0_dp)
    end function xdpsi0


    ! psi1 & derivatives
    function xpsi1(z) result(xpsi1r)
      !$acc routine seq
      real(dp), intent(in) :: z
      real(dp) :: xpsi1r
      xpsi1r = z * ( z * (z - 2.0_dp) + 1.0_dp)
    end function xpsi1

    function xdpsi1(z) result(xdpsi1r)
      !$acc routine seq
      real(dp), intent(in) :: z
      real(dp) :: xdpsi1r
      xdpsi1r = z * (3.0_dp*z - 4.0_dp) + 1.0_dp
    end function xdpsi1

    ! bicubic hermite polynomial function
    function h3(dfi,w0t,w1t,w0mt,w1mt,w0d,w1d,w0md,w1md) result(h3r)
      !$acc routine seq
      real(dp), intent(in) :: dfi(16)
      real(dp), intent(in) :: w0t,w1t,w0mt,w1mt,w0d,w1d,w0md,w1md
      real(dp) :: h3r
      h3r =  dfi(1)  *w0d*w0t   +  dfi(2)  *w0md*w0t &
           + dfi(3)  *w0d*w0mt  +  dfi(4)  *w0md*w0mt &
           + dfi(5)  *w0d*w1t   +  dfi(6)  *w0md*w1t &
           + dfi(7)  *w0d*w1mt  +  dfi(8)  *w0md*w1mt &
           + dfi(9)  *w1d*w0t   +  dfi(10) *w1md*w0t &
           + dfi(11) *w1d*w0mt  +  dfi(12) *w1md*w0mt &
           + dfi(13) *w1d*w1t   +  dfi(14) *w1md*w1t &
           + dfi(15) *w1d*w1mt  +  dfi(16) *w1md*w1mt
    end function h3

    subroutine actual_eos_finalize

      implicit none

      ! Deallocate managed module variables

      deallocate(do_coulomb)
      deallocate(input_is_constant)
      deallocate(itmax)
      deallocate(jtmax)
      deallocate(d)
      deallocate(t)
      deallocate(tlo)
      deallocate(thi)
      deallocate(tstp)
      deallocate(tstpi)
      deallocate(dlo)
      deallocate(dhi)
      deallocate(dstp)
      deallocate(dstpi)
      deallocate(ttol)
      deallocate(dtol)
      deallocate(f)
      deallocate(fd)
      deallocate(ft)
      deallocate(fdd)
      deallocate(ftt)
      deallocate(fdt)
      deallocate(fddt)
      deallocate(fdtt)
      deallocate(fddtt)
      deallocate(dpdf)
      deallocate(dpdfd)
      deallocate(dpdft)
      deallocate(dpdfdt)
      deallocate(ef)
      deallocate(efd)
      deallocate(eft)
      deallocate(efdt)
      deallocate(xf)
      deallocate(xfd)
      deallocate(xft)
      deallocate(xfdt)
      deallocate(dt_sav)
      deallocate(dt2_sav)
      deallocate(dti_sav)
      deallocate(dt2i_sav)
      deallocate(dd_sav)
      deallocate(dd2_sav)
      deallocate(ddi_sav)
      deallocate(dd2i_sav)

    end subroutine actual_eos_finalize

end module actual_eos_module
