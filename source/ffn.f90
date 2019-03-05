!*******************************************************************************
! FFN.f90 for XNet version 6 5/28/10
! This file contains reaction data structures for Fuller, Fowler, Neuman
! (FFN; 1982,1985) formated weak reactions and routines to read the data
! and calculate the reaction rates.
!*******************************************************************************

Module ffn_data
!-----------------------------------------------------------------------------
! This module contains the data to calculate FFN formatted weak reactions.
!-----------------------------------------------------------------------------
  Use cross_sect_data
  Real(8), Dimension(:,:), Allocatable :: ffnsum,ffnenu! dim(nffn,143)

End Module ffn_data

Subroutine read_ffn_data(nin,data_dir)
!-----------------------------------------------------------------------------
! This routine allocates and loads the data structures for FFN reaction rates
!-----------------------------------------------------------------------------
  Use ffn_data
  Integer, Intent(in) :: nin
  Integer             :: i,j
  Character (LEN=*) :: data_dir
  Character (LEN=1) :: null
  Open(3,FILE=trim(data_dir)//"/netweak",status='old')
  nffn=nin

  Allocate (ffnsum(nffn,143),ffnenu(nffn,143))
  Do i=1,nffn
    Read(3,"(a1)") null
    Read(3,"(9(f8.3))") (ffnsum(i,j),ffnenu(i,j),j=1,143)
  EndDo
End Subroutine read_ffn_data

Subroutine ffn_rate(t9,ene,rf,dlnrfdt9)
!-----------------------------------------------------------------------------
! This routine calculates the reaction rates for FFN weak rates
!-----------------------------------------------------------------------------
  Use controls
  Use ffn_data
  Real(8), Intent (in)  :: t9  ! Temperature in GK
  Real(8), Intent (in)  :: ene ! Electron Density (CGS)
  Real(8), Intent (out) :: rf(nffn) ! Temperature and density dependent FFN rates
  Real(8), Intent (out) :: dlnrfdt9(nffn) ! Temperature and density dependent log FFN rate derivatives
  Real(8)               :: r1(nffn),r2(nffn) ! temporary interpolation variables
  Real(8) :: tg(13) = (/ 0.01,0.1,0.2,0.4,0.7,1.0,1.5,2.,3.,5.,10.,30.,100./)
  Real(8) :: egl(11) = (/1.,2.,3.,4.,5.,6.,7.,8.,9.,10.,11./)
  Real(8) :: enl,dt,de,rdt,rde
  Integer i1(4)
  Integer i,le1,lt1
  Do i=1,13
    If(t9<=tg(i)) Exit
  EndDo
  lt1=i-1
  If(lt1<1) lt1=1
  enl=log10(ene)
  le1=int(enl)
  If(le1<1) le1=1
  dt=tg(lt1+1)-tg(lt1)
  rdt=(t9-tg(lt1))/dt
  de=egl(le1+1)-egl(le1)
  rde=(enl-egl(le1))/de
  i1(1)=13*(le1-1)+lt1
  i1(2)=i1(1)+1
  i1(3)=13*le1+lt1
  i1(4)=i1(3)+1
  r1(:)=ffnsum(:,i1(1))+(ffnsum(:,i1(2))-ffnsum(:,i1(1)))*rdt
  r2(:)=ffnsum(:,i1(3))+(ffnsum(:,i1(4))-ffnsum(:,i1(3)))*rdt
  rf(:)=r1(:)+(r2(:)-r1(:))*rde
  Where(rf(:).lt.-30.)
    rf(:)=0.
  ElseWhere
    rf(:)=10.**rf(:)
  EndWhere
  If(iheat>0) dlnrfdt9=log(10.)*(rde*(ffnsum(:,i1(4))-ffnsum(:,i1(3)))+(1-rde)*(ffnsum(:,i1(2))-ffnsum(:,i1(1))))/dt
  Return
End Subroutine ffn_rate
