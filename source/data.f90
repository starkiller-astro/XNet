!-----------------------------------------------------------------------------  
! Data.f for Full_Net 4.10 2/26/07  
! This file contains the nuclear and reaction data structures and the 
! Subroutines which read in the data and Allocate the arrays.
!-----------------------------------------------------------------------------  

Module constants
!===============================================================================
! These are fundamental constants in CGS, except energies in MeV, Temp in GK
!===============================================================================
  Real(8), Parameter :: pi   =3.1415926536, hbar=6.582122E-22
  Real(8), Parameter :: amu  =1.036427E-18, bok =8.617385E-02
  Real(8), Parameter :: avn  =6.022137E+23, e2  =1.439830E-13
  Real(8), Parameter :: epmev=1.602177E-06, ele_en=.510999
  Real(8), parameter :: clt  =2.997925E+10
End Module constants

Module nuc_number
!-----------------------------------------------------------------------------
! This module contains ny, the number of nuclear species whose abundances
! are evolved by the network. The value of ny is read in by read_reaction_data
! or read_nuclear_data. nymx is only Used in the network setup.
!-----------------------------------------------------------------------------
  Integer             ::  ny
End Module nuc_number
  
Module nuclear_data
!-----------------------------------------------------------------------------
! This module contains the essential data for each included specie.
! aa, zz, & nn are the mass, proton and neutron numbers, and be is the 
! binding energy. 
! nname is the nuclear name (nname(0) is placeholder for non-nuclei). 
! Their array sizes are set in the routine read_nuclear_data. 
!-----------------------------------------------------------------------------
  Use nuc_number
  Real(8), Dimension(:), Allocatable   :: aa,zz,nn,be,mex
  Integer, Dimension(:), Allocatable   :: ia,iz,in
  Integer :: inmin,inmax,izmin,izmax
  Character(len=5), Dimension(:), Allocatable  :: nname
  
End Module nuclear_data
  
Module part_funct_data
!-----------------------------------------------------------------------------
! This module contains the nuclear partition function data.  The value of
! the partition function is interpolated over a grid of temperatures, t9i.
! g is the interpolation data (24,ny), gg is the current parition function,
! and angm is the J. gg(0) and angm(0) are placeholders for non-nuclei
! The array sizes are set in read_nuclear_data.
!-----------------------------------------------------------------------------
  Use nuc_number
  Real(8) ::  t9i(24)
  Real(8), Dimension(:), Allocatable   :: gg,angm
  Real(8), Dimension(:,:), Allocatable :: g
  Real(8), Dimension(:), Allocatable   :: dlngdt9
  
! Threading Scope
!$OMP THREADPRIVATE(gg,dlngdt9)
  
End Module part_funct_data
  
Module cross_sect_data
!-----------------------------------------------------------------------------
! This module contains the data needed to calculate the cross sections
! nreac(i) are the # of reactions with i reactants.  
! n{1,2,3}i list the nuclei affected by each reaction.  
! The rc variables are data for the temperature dependent part of the csect.
! The csect variables are the results, velocity integrated cross section * density dependence.
!-----------------------------------------------------------------------------
  Integer                              :: nreac(3)
  Integer, Dimension(:,:), Allocatable :: n1i,n2i,n3i
  Real(8), Dimension(:,:), Allocatable :: rc1,rc2,rc3    ! dim(7,nreac)
  Real(8), Dimension(:),   Allocatable :: csect1,csect2,csect3
  Real(8), Dimension(:),   Allocatable :: dcsect1dt9,dcsect2dt9,dcsect3dt9

! Reaction flags to indicate variations in how the rate is calculated
  Integer, Dimension(:),   Allocatable :: iwk1,iwk2,iwk3  ! Weak Reaction
  Integer, Dimension(:),   Allocatable :: ires1,ires2,ires3 ! Resonant React
  Integer, Dimension(:),   Allocatable :: irev1,irev2,irev3 ! Reverse React

! Additional rate or energy factors
  Real(8), Dimension(:),   Allocatable :: q1,q2,q3  ! reaction Q values
  Real(8), Dimension(:),   Allocatable :: rpf1,rpf2 ! ratios of Part. Functions.
  Real(8), Dimension(:),   Allocatable :: dlnrpf1dt9,dlnrpf2dt9 ! ratios of Part. Functions.
  Real(8), Dimension(:),   Allocatable :: h1,h2,h3  ! screening factors
  Real(8), Dimension(:),   Allocatable :: dh1dt9,dh2dt9,dh3dt9 ! d(screening factors)/dT

  Real(8), Dimension(:),   Allocatable :: z21,z22,z12,z12i,z31,z32,z33,z123 ! proton numbers of reactants
  Integer, Dimension(:),   Allocatable :: iz21,iz22,iz12,iz31,iz32,iz33,iz123
  
!-----------------------------------------------------------------------------
! Aside from the REACLIB formated data, this dataset includes pointers to sets
! of external reaction rates, in the form of indices incoded in rc{1,2,3}(1).
! FFN are electron and positron capture rates encoded in the data format of 
! Fuller, Fowler & Neuman (1982,1985).
! NNU are neutrino and antineutrino capture rates from Zinner & Langanke, 
! implemented by Carla Froehlich.
!-----------------------------------------------------------------------------
  Integer                              :: nffn ! The number of FFN reactions
  Integer                              :: nnnu ! The number of NNU reactions
  Integer, Dimension(:),   Allocatable :: iffn ! Pointers to FFN list
  Integer, Dimension(:),   Allocatable :: innu ! Pointers to NNU list

! Threading Scope
!$OMP THREADPRIVATE(csect1,csect2,csect3,rpf1,rpf2,h1,h2,h3)
!$OMP THREADPRIVATE(dcsect1dt9,dcsect2dt9,dcsect3dt9,dlnrpf1dt9,dlnrpf2dt9,dh1dt9,dh2dt9,dh3dt9)
  
End Module cross_sect_data
  
Module reac_rate_data
!-----------------------------------------------------------------------------
! This module contains the data necessary to calculate the reaction rates
! and to map to each species those reactions which affect it.
!-----------------------------------------------------------------------------
  Use nuc_number
  Integer                              ::  nan(3)
  Integer, Dimension(:,:), Allocatable :: la,le
  Integer, Dimension(:),   Allocatable :: mu1,mu2,mu3
  Integer, Dimension(:),   Allocatable :: n11,n21,n22,n31,n32,n33
  Real(8), Dimension(:),   Allocatable :: a1,a2,a3,b1,b2,b3 ! dim(nan(i))
  Real(8), Dimension(:),   Allocatable :: dr1dt9,dr2dt9,dr3dt9 ! dim(nan(i))
  
! Threading Scope
!$OMP THREADPRIVATE(b1,b2,b3,dr1dt9,dr2dt9,dr3dt9)
  
End Module reac_rate_data
  
Subroutine read_nuclear_data(data_dir,data_desc)
!-----------------------------------------------------------------------------  
! This routine reads, from the file netwinv, the nuclei included along 
! with the nuclear data which will be needed for later calculations.  This 
! data includes the atomic number, the number of protons and neutrons, and 
! the binding energy (calculated from the tabulated mass excess).  Also the 
! tabulations of the  partition functions, g, are read in for later 
! interpolation. Once the set of nuclear data is read in, it is assigned 
! to the proper nuclei.
!-----------------------------------------------------------------------------  
  Use controls
  Use nuclear_data
  Use part_funct_data
  Character (LEN=*),  INTENT(in)  :: data_dir
  Character (LEN=80), INTENT(out) :: data_desc
  Real(8) :: a,sp 
  Integer :: i,l,n,m,na,nb                                     
  Integer :: it9i(24)
  Character (LEN=5) nam                                 
    
! Read in the data description
  Open(13,file=trim(data_dir)//"/net_desc",status='old')
  Read(13,"(A)") data_desc
  Open(12,FILE=trim(data_dir)//"/netwinv",status='old')
  Read(12,"(i5)") ny
    
! Read in the partition function iteration grid, and fix Endpoints
  Read(12,"(24i3)") (it9i(i),i=1,24)                   
  Do i=1,24                                                    
    t9i(i)=it9i(i)*0.01     
  EndDo                                             
  t9i(24)=t9i(24)*10.                                                       
  Do n=1,ny
    Read(12,"(a5)") nam
  EndDo                     
    
! Set size of nuclear Parameter arrays and read in nuclear Parameters 
! and partition function interpolation table. nname(0), gg(0) and angm(0) 
! are placeholders for non-nuclei.
  Allocate (nname(0:ny))
  Allocate (aa(ny),zz(ny),nn(ny),be(ny),ia(ny),iz(ny),in(ny),mex(ny))
  Allocate (g(24,ny),angm(0:ny))
!$OMP PARALLEL
  Allocate (gg(0:ny))
  If(iheat>0) Allocate (dlngdt9(0:ny))
!$OMP End PARALLEL
  nname(0)=' === '
  angm(0)=0.0
  Do l=1,ny                                                  
    Read(12,"(a5,f12.3,2i4,f6.1,f15.8)") nam,aa(l),iz(l),in(l),sp,mex(l)
    Read(12,*) (g(m,l),m=1,24)     
    ia(l)=int(aa(l))                                                         
    zz(l)=real(iz(l))                                                
    nn(l)=real(in(l))                                                
    be(l)=8.07144*nn(l)+7.28899*zz(l)-mex(l)
    angm(l)=2.*sp+1.
    nname(l)=nam
  EndDo
  inmin=minval(in) ; inmax=maxval(in)
  izmin=iz(1)      ; izmax=iz(ny)                                                           
! Write(lun_diag,"(a5,4es12.4)") (nname(i),aa(i),zz(i),nn(i),be(i),i=1,ny)        
  Return                                                                    
  End                                                                       
  
Subroutine read_reaction_data(data_dir)    
!-----------------------------------------------------------------------------  
! This routine reads in the necessary reaction data
!-----------------------------------------------------------------------------  
  Use controls
  Use nuclear_data
  Use cross_sect_data
  Use reac_rate_data
  Character (LEN=*) :: data_dir
  Integer :: i,j,n,l
  Integer :: nr,idel,jdel,ijd,ndi,nds
  Character (LEN=1) :: null
  common /ceign/ nds,ndi,idel,jdel,ijd
    
! Read in nuclear set and numbers of reactions
  Open(4,file=trim(data_dir)//"/nets4",form='unformatted',status='old')
  Read(4) ny                                                               
  Read(4) (nname(i),i=1,ny)  
  Read(4) nffn,nnnu 
  Read(4) (nreac(i),i=1,3)   
    
! If there are FFN rates, read in the FFN data and set FFN array sizes
  If(nffn>0) Call read_ffn_data(nffn,data_dir)                          !FFN
    
! If there are NNU rates, read in the NNU data and set NNU array sizes
! If(nnnu>0) Call read_nnu_data(nnnu,data_dir)                          !NNU
    
! Read in the reaction cross section data
  Open(2,file=trim(data_dir)//"/nets3",form='unformatted',status='old')
    
! Allocate and read in reaction arrays for 1 reactant reactions
  nr=nreac(1)
  Allocate (rc1(7,nr),q1(nr))
  Allocate (n1i(4,nr),iwk1(nr),ires1(nr),irev1(nr),iffn(nr),innu(nr))
!$OMP PARALLEL      
  Allocate (csect1(nr),rpf1(nr),h1(nr))
  If(iheat>0) Allocate (dcsect1dt9(nr),dlnrpf1dt9(nr),dh1dt9(nr))
!$OMP End PARALLEL
  Do j=1,nr
! Read REACLIB data
    Read(2) n,(n1i(l,j),l=1,4),iwk1(j),ires1(j),irev1(j),(rc1(l,j),l=1,7),q1(j)  
    If(n/=j) Then
        Write(6,*) 'Error in nets3, 1',j,n
        Exit
    EndIf
! Calculate pointers to non-REACLIB data
    If(iwk1(j)==2.or.iwk1(j)==3) Then
      iffn(j)=int(rc1(1,j))
    ElseIf(iwk1(j)==7.or.iwk1(j)==8) Then
      innu(j)=int(rc1(1,j))
    Else
      iffn(j)=0 ; innu(j)=0
    Endif
  EndDo
    
! Allocate and read in reaction arrays for 2 reactant reactions
  nr=nreac(2)
  Allocate (rc2(7,nr),q2(nr))
  Allocate (n2i(5,nr),iwk2(nr),ires2(nr),irev2(nr))
!$OMP PARALLEL      
  Allocate (csect2(nr),rpf2(nr),h2(nr))
  If(iheat>0) Allocate (dcsect2dt9(nr),dlnrpf2dt9(nr),dh2dt9(nr))
!$OMP End PARALLEL
  Do j=1,nr
    Read(2)  n,(n2i(l,j),l=1,5),iwk2(j),ires2(j),irev2(j),(rc2(l,j),l=1,7),q2(j)                           
    If(n/=j) Then
        Write(6,*) 'Error in nets3, 2',j,n
        Exit
    EndIf
  EndDo
  If(iscrn>0) Then
    Allocate (z21(nr),z22(nr),iz21(nr),iz22(nr))
    z21=zz(n2i(1,:)) ; z22=zz(n2i(2,:))
    iz21=int(z21)    ; iz22=int(z22)
    Allocate (z12(nr),iz12(nr),z12i(nr))
    z12=z21*z22
    iz12=iz21+iz22
    z12i=((z21+z22)**1.86-z21**1.86-z22**1.86)
  EndIf
    
! Allocate and read in reaction arrays for 3 reactant reactions
  nr=nreac(3)
  Allocate (rc3(7,nr),q3(nr))
  Allocate (n3i(6,nr),iwk3(nr),ires3(nr),irev3(nr))
!$OMP PARALLEL      
  Allocate (csect3(nr),h3(nr))
  If(iheat>0) Allocate (dcsect3dt9(nr),dh3dt9(nr))
!$OMP End PARALLEL
  Do j=1,nr
    Read(2)  n,(n3i(l,j),l=1,6),iwk3(j),ires3(j),irev3(j),(rc3(l,j),l=1,7),q3(j)                                                
    If(n/=j) Then
        Write(6,*) 'Error in nets3, 3',j,n
        Exit
    EndIf
  EndDo
  If(iscrn>0) Then
    Allocate (z31(nr),z32(nr),z33(nr),iz31(nr),iz32(nr),iz33(nr))
    z31=zz(n3i(1,:)) ; z32=zz(n3i(2,:)) ; z33=zz(n3i(3,:))
    iz31=int(z31)    ; iz32=int(z32)    ; iz33=int(z33)
    Allocate (z123(nr),iz123(nr))
    z123=z31*z32+z31*z33+z32*z33
    iz123=iz31+iz32+iz33
  EndIf
    
! Allocate and read in the data linking nuclei to the reactions which 
! affect them.  Also read in the matrix sparseness descriptors.
  Allocate (la(3,ny),le(3,ny))
  Do i=1,ny                                                             
    Read(4)  n,la(1,i),le(1,i),la(2,i),le(2,i),la(3,i),le(3,i)
    If(n.ne.i) Then
        Write(6,*) 'Error in nets4',i,n
    EndIf
  EndDo
  Read(4) idel,jdel,ijd,ndi,nds
!    Write(lun_diag,*) 'ceign',idel,jdel,ijd,ndi,nds
    
! Create and fill extended reaction->nuclei arrays
  nan(1)=le(1,ny)
  Allocate (mu1(nan(1)),a1(nan(1)),n11(nan(1)))
!$OMP PARALLEL
  Allocate (b1(nan(1)))
  If(iheat>0) Allocate (dr1dt9(nan(1)))
!$OMP End PARALLEL
  Do j=1,nan(1)
    Read( 2)  a1(j),mu1(j) 
    n11(j)=n1i(1,mu1(j))
  EndDo                                                                  
  nan(2)=le(2,ny)
  Allocate (mu2(nan(2)),a2(nan(2)))
  Allocate (n21(nan(2)),n22(nan(2)))
!$OMP PARALLEL
  Allocate (b2(nan(2)))
  If(iheat>0) Allocate (dr2dt9(nan(2)))
!$OMP End PARALLEL
  Do j=1,nan(2)
    Read( 2)  a2(j),mu2(j)  
    n21(j)=n2i(1,mu2(j))
    n22(j)=n2i(2,mu2(j))
  EndDo                                                                     
  nan(3)=le(3,ny)
  Allocate (mu3(nan(3)),a3(nan(3)))
  Allocate (n31(nan(3)),n32(nan(3)),n33(nan(3)))
!$OMP PARALLEL
  Allocate (b3(nan(3)))
  If(iheat>0) Allocate (dr3dt9(nan(3)))
!$OMP End PARALLEL
  
  Do j=1,nan(3)
    Read( 2)  a3(j),mu3(j)
    n31(j)=n3i(1,mu3(j))
    n32(j)=n3i(2,mu3(j))
    n33(j)=n3i(3,mu3(j))
  EndDo                                                                  
  Return                                                                    
  End                                                                       
  
Subroutine index_from_name(name,index)
!-----------------------------------------------------------------------
! This Subroutine takes a nuclear name and finds the corresponding
! index for the current data set. index=0 indicates that the nuclear name
! is not found in the current set.
!-----------------------------------------------------------------------
  Use nuclear_data
  Character (len=5) :: name,sname
  Integer :: n,index,name_len
  
  name_len=len_trim(name)
  If(name_len==0) then
      index = 0
      return
  elseIf(name_len<5) then
      sname='     '
      sname((6-name_len):5)=name(1:name_len)
  else
      sname=name
  EndIf
  Do n=1,ny
    If(sname==nname(n)) Exit
  EndDo
  
  index=n
  
  Return
End Subroutine index_from_name
  
Subroutine name_ordered(name,num,max)
!-----------------------------------------------------------------------------
! This routine appEnds the number num (padded with "0"s up to the size of 
! the number max) to the Character variable name.
!-----------------------------------------------------------------------------
  Character (LEN=*)   :: name
  Character (LEN=9)   :: num_string,form
  Character (LEN=1)   :: lmax_string
  Integer             :: num,max,lnum,lmax
    
! Find Character length of max
  lmax=int(log10(float(max)))+1
  Write(unit=lmax_string,fmt='(I1)') lmax
    
! Construct Format Spec and write num as zero padded string
  form='(I'//lmax_string//'.'//lmax_string//')'
  Write(unit=num_string,fmt=form) num
    
! AppEnd num_ string to name
  name=trim(name)//trim(num_string)
  Return
End Subroutine name_ordered
  
