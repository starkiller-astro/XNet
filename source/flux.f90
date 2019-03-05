Module flux_data
!-----------------------------------------------------------------------------
! This module contains the data for calculating the net flux of each 
! matched reaction pair.
!-----------------------------------------------------------------------------
  Real(8),  Dimension(:), Allocatable    :: flx,flx_int 
  Real(8),  Dimension(:,:), Allocatable  :: dcflx
  Integer, Dimension(:), Allocatable    :: ifl_orig,ifl_term
  
! Threading Scope
!$OMP THREADPRIVATE(flx,flx_int)
  
End Module flux_data
  
Subroutine flux_init()
!-----------------------------------------------------------------------------
! This routine allocates the flux arrays and determines the double
! counting factors necessary for reactions with identical reactants.
!-----------------------------------------------------------------------------
  Use controls
  Use nuclear_data
  Use match_data
  Use flux_data
  Real(8), Dimension(5)  :: fact=1.0/(/1.0,2.0,6.0,24.0,4.0/)
  Integer                :: i,j,countf,countr
  
!$OMP PARALLEL DEFAULT(SHARED)
  Allocate(flx(mflx),flx_int(mflx))
!$OMP End PARALLEL
  Allocate(dcflx(2,mflx),ifl_orig(mflx),ifl_term(mflx))
  
! Check double counting of reactants for each flux
!$OMP PARALLEL DEFAULT(SHARED)
  If(idiag>=3) Write(lun_diag,'(a)') 'Flux Init'
!$OMP End PARALLEL
  Do i=1,mflx
    countf=1
    countr=1
    Do j=2,3
      If(nflx(j,i)/=0.and.nflx(j,i)==nflx(j-1,i)) countf=countf+1
    EndDo
    Do j=5,7
      If(nflx(j,i)/=0.and.nflx(j,i)==nflx(j-1,i)) countr=countr+1
    EndDo
  
! countr=3 can result from a+a+a+b or a+a+b+b, which has different factor
    If(countr==3.and.nflx(5,i)/=nflx(6,i)) countr=5
    dcflx(1,i)=fact(countf)
    dcflx(2,i)=fact(countr)
  
! Determine flux origin, ifl_orig, and terminus, ifl_term
    ifl_orig(i)=nflx(count(nflx(1:3,i)/=0),i)
    ifl_term(i)=nflx(count(nflx(4:7,i)/=0)+3,i)
!$OMP PARALLEL DEFAULT(SHARED)
    If(idiag>=3) Write(lun_diag,'(11i5,2f6.3)')&  
&     nflx(:,i),ifl_orig(i),ifl_term(i),countf,countr,dcflx(:,i)
!$OMP End PARALLEL
  EndDo
  Return
End Subroutine flux_init
  
Subroutine flux()
!-----------------------------------------------------------------------------
! This routine calculates the net flux of each matched reaction pair.
! Positive fluxes flow in the direction of increasing binding.
!-----------------------------------------------------------------------------
  Use controls
  Use conditions
  Use cross_sect_data
  Use nuclear_data
  Use abundances
  Use match_data
  Use flux_data
  Integer :: i,k,ifl,idcfl
  Real(8)  :: flt,flxmin
  flx=0.0
  Do i=1,nreac(1)
    ifl=abs(ifl1(i))
    idcfl=int(1.5-sign(0.5,Real(ifl1(i))))
    flt=csect1(i)*yt(n1i(1,i))
    flx(ifl)=flx(ifl)+dcflx(idcfl,ifl)*sign(flt,Real(ifl1(i),kind(flt)))
  EndDo
  Do i=1,nreac(2)
    ifl=abs(ifl2(i))
    idcfl=int(1.5-sign(0.5,Real(ifl2(i))))
    flt=csect2(i)*yt(n2i(1,i))*yt(n2i(2,i))
    flx(ifl)=flx(ifl)+dcflx(idcfl,ifl)*sign(flt,Real(ifl2(i),kind(flt)))
  EndDo
  Do i=1,nreac(3)
    ifl=abs(ifl3(i))
    idcfl=int(1.5-sign(0.5,Real(ifl3(i))))
    flt=csect3(i)*product(yt(n3i(1:3,i)))
    flx(ifl)=flx(ifl)+dcflx(idcfl,ifl)*sign(flt,Real(ifl3(i),kind(flt)))
  EndDo
  
! Since we cut off abundances less than ymin, we should similarly cut off 
! small fluxes
  flxmin=ymin/tdel
  Where(abs(flx)>flxmin)
      flx_int=flx_int+flx*tdel
  EndWhere
  If(idiag>=3) Then
    Write(lun_diag,'(a5,es10.3)') 'Flux',flxmin
    Write(lun_diag,'(i5,8a5,i5,2es10.3)') (k,nname(nflx(1:3,k)),' <-> ',&
&       nname(nflx(4:7,k)),iwflx(k),flx(k),flx_int(k),k=1,mflx)
  EndIf
  Return
End Subroutine flux
  
Subroutine flux_check()
!-----------------------------------------------------------------------------
! This routine compares the fluxes to the changes in abundance
!-----------------------------------------------------------------------------
  Use controls
  Use nuclear_data
  Use match_data
  Use flux_data
  Use abundances
  Use conditions
  Real(8), Dimension(0:ny) :: flx_test
  Integer i
  flx_test(0)=0.0
!    flx_test(1:ny)=-ydot*tdel ! Use only with convc>0
  flx_test(1:ny)=yo-y
  Do i=1,mflx
    flx_test(nflx(1:3,i))=flx_test(nflx(1:3,i))-flx(i)*tdel
    flx_test(nflx(4:7,i))=flx_test(nflx(4:7,i))+flx(i)*tdel
  EndDo
  If(idiag>=3) Then
    Write(lun_diag,"(a,i5,es11.3)") "Flux Check",mflx,tdel
    Write(lun_diag,"(a5,4es11.3)") &
&     (nname(i),y(i),yo(i)-y(i),ydot(i)*tdel,flx_test(i),i=1,ny)
  EndIf
  Return
End Subroutine flux_check
  
