module rad
use rad_param
use rad_config
use rad_spectral
use rad_rte
use rad_fields, only: Srad, dimension_bc, smax
use rad_fields, only: k_g, k_s
implicit none
private

public :: rad_init ! inialization
public :: rad_final ! finalization, deallocation variables
public :: rad_calc ! calculate new radiative heat sources
public :: rad_write_src ! export rad sources
public :: S_Rpg, S_Rcg, S_Rps, S_Rcs ! rad heat sources for energy equation

type(configuration) :: config
contains

subroutine rad_init()
    use rad_fields, only: radiationOn, nq, ew, Tw
    implicit none
    include 'usrnlst.inc'
    radiationOn = RAD_ON
    if (.not.radiationOn) return

    Tw = RAD_T_W
    ew = RAD_EMIS_W
    ! create configuration data
    call nullconfig(config)
    config%SpectralModelName=RAD_SPECTRAL
    config%RTEModelName=RAD_RTE
    config%nQuad = nq
    config%skipSteps = RAD_SKIP
    config%nrr = RAD_NRR
    call to_upper(config)
    call rad_spectral_init(config)
    call rad_rte_init(config)
end subroutine rad_init

subroutine rad_final()
    use rad_fields, only: radiationOn
    if (.not.radiationOn) return
    call rad_spectral_final()
    call rad_rte_final()
end subroutine rad_final

subroutine rad_calc()
    use rad_fields, only: radiationOn
    implicit none
    if (.not.radiationOn) return
    ! spectral calculation
    call rad_spectral_calc()
    ! rte calculation
    call rad_rte_calc()
    ! source calculation
    call rad_spectral_srad()
end subroutine rad_calc

subroutine rad_write_src()
    use rad_fields, only: radiationOn, ReactionRates, nRR
    integer :: irr, m
    if (.not.radiationOn) return
    irr = config%nrr
    if (irr.eq.0) return
    if (irr.le.nRR) then
        ReactionRates(:,irr) = Srad(:,0)
    endif
    do m=1, smax
        irr = m+config%nrr
        if (irr .le. nRR) then
            ReactionRates(:,irr) = Srad(:,m)
        endif
    enddo
end subroutine

real(dp) function S_Rpg(ijk)
    implicit none
    integer ijk
    S_Rpg = 0.d0
end function

real(dp) function S_Rcg(ijk)
    use rad_fields, only: radiationOn
    implicit none
    integer ijk
    if (radiationOn) then
        S_Rcg = Srad(ijk,0)
    else
        S_Rcg = 0.d0
    endif
end function

real(dp) function S_Rps(ijk,m)
    implicit none
    integer ijk, m
    S_Rps = 0.d0
end function

real(dp) function S_Rcs(ijk, m)
    use rad_fields, only: radiationOn
    implicit none
    integer ijk, m
    if (radiationOn) then
        S_Rcs = Srad(ijk,m)
    else
        S_Rcs = 0.d0
    endif
end function

end module rad
