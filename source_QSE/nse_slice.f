      module constants
!-----------------------------------------------------------------------------
!  These are fundamental constants in CGS, except energies in MeV, Temp in GK
!-----------------------------------------------------------------------------
      real*8, parameter :: pi  =3.1415926536, hbar=6.582122E-22
      real*8, parameter :: amu =1.036427E-18, bok =8.617385E-02
      real*8, parameter :: avn =6.022137E+23, e2  =1.439830E-13
      end module constants

      program nse_slice
      use nse_abundance
      use cnse_data
      use screen_nse
      use nuclear_data
      real(8) :: rho,ye,t9fnl
      character*80 :: data_dir,data_desc
      real(8) :: t9w(27)=(/100.,75.,50.,40.,35.,30.,28.,26.,24.,22.,20.,
     &         18.,16.,15.,14.,13.,12.,11.,10.,9.,8.,7.,6.,5.,4.,3.,2./)

!  For a given density and Ye, this program solves for NSE as a funtion of T9 
!  down to T9fnl.

      Write(6,'(a)') 'Rho?'
      Read(5,*) rho
      Write(6,'(a)') 'Ye?'
      Read(5,*) ye
      Write(6,'(a)') 'T9 stop?'
      Read(5,*) t9fnl
      Write(6,'(a,f6.4,a,es11.3,a,f7.3)') 
     &    'NSE for Ye=',ye,'and Rho=',rho,', stopping at T9=',t9fnl

!  Read nuclear dataset and allocate nse arrays
      data_dir='Data_NSE'
      call read_nuclear_data(data_dir,data_desc)
      zmax=maxval(zz)
      Allocate(ynse(ny),cnse(ny),he(zmax),intz(ny))
      Where(zz>1.0) 
        intz=int(zz)
      ElseWhere
        intz=1
      EndWhere

!  Open output files
      Open(12,file='nsl_out')
      Write(12,'(2es13.6)') rho,ye
      Open(14,file='nsl_ab')
      Open(15,file='nse_diag')

!  Descend to desired temperature starting from high temperature.
      call nse_descend(rho,ye,t9fnl,t9w)
 
      End                                                                       

