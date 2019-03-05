!-----------------------------------------------------------------------------  
!  Data.f for Full_Net 4.10 2/26/07  
!  This file contains the nuclear and reaction data structures and the 
!  subroutines which read in the data and allocate the arrays.
!-----------------------------------------------------------------------------  

      module nuc_number
!-----------------------------------------------------------------------------
!  This module contains ny, the number of nuclear species whose abundances
!  are evolved by the network. The value of ny is read in by read_reaction_data
!  or read_nuclear_data
!-----------------------------------------------------------------------------
      integer             ::  ny
      end module nuc_number

      module nuclear_data
!-----------------------------------------------------------------------------
!  This module contains the essential data for each included specie.
!  aa, zz, & nn are the mass, proton and neutron numbers, and be is the 
!  binding energy. 
!  nname is the nuclear name (nname(0) is placeholder for non-nuclei). 
!  Their array sizes are set in the routine read_nuclear_data. 
!-----------------------------------------------------------------------------
      use nuc_number
      real(8), dimension(:), allocatable   :: aa,zz,nn,be
      character(len=5), dimension(:), allocatable  :: nname
      end module nuclear_data

      module part_funct_data
!-----------------------------------------------------------------------------
!  This module contains the nuclear partition function data.  The value of
!  the partition function is interpolated over a grid of temperatures, t9i.
!  g is the interpolation data (24,ny), gg is the current parition function,
!  and angm is the J. gg(0) and angm(0) are placeholders for non-nuclei
!  The array sizes are set in read_nuclear_data.
!-----------------------------------------------------------------------------
      use nuc_number
      real(8) ::  t9i(24)
      real(8), dimension(:), allocatable   :: gg,angm
      real(8), dimension(:,:), allocatable :: g
      end module part_funct_data

      module cross_sect_data
!-----------------------------------------------------------------------------
!  This module contains the data needed to calculate the cross sections
!  nreac(i) are the # of reactions with i reactants.  
!  n?i list the nuclei affected by each reaction.  
!  The csect variables are the results.
!  The rc variables are data for the temperature dependant part of the csect.
!-----------------------------------------------------------------------------
      integer                              :: nreac(3)
      real(8),  dimension(:),   allocatable :: csect1,csect2,csect3
      real(8),  dimension(:,:), allocatable :: rc1,rc2,rc3   ! dim(7,nreac)
      integer, dimension(:,:), allocatable :: n1i,n2i,n3i
      integer, dimension(:),   allocatable :: iwk1,iwk2,iwk3    ! Weak Reaction
      integer, dimension(:),   allocatable :: ires1,ires2,ires3 ! Resonant React
      integer, dimension(:),   allocatable :: irev1,irev2,irev3 ! Reverse React
      integer, dimension(:),   allocatable :: iffn ! Maps reaction to FFN list
      real(8),  dimension(:),   allocatable :: rpf1,rpf2,h2,h3 ! ratios of PartF.
      real(8),  dimension(:),   allocatable :: q1,q2,q3 ! reaction Q values
      end module cross_sect_data

      module reac_rate_data
!-----------------------------------------------------------------------------
!  This module contains the data necessary to calculate the reaction rates
!  and to map to each species those reactions which affect it.
!-----------------------------------------------------------------------------
      use nuc_number
      integer   ::  nan(3)
      integer, dimension(:,:), allocatable :: la,le
      integer, dimension(:), allocatable :: mu1,mu2,mu3
      integer, dimension(:), allocatable :: n11,n21,n22,n31,n32,n33
      real(8),  dimension(:), allocatable :: a1,a2,a3,b1,b2,b3 ! dim(nan(i))
      end module reac_rate_data

      module ffn_data
!-----------------------------------------------------------------------------
!  This module contains the data to calculate weak reactions according to
!  Fuller, Fowler, Neuman (1982,1985).
!-----------------------------------------------------------------------------
      integer :: nffn ! The number of FFN reactions
      real(8), dimension(:),   allocatable :: rf,r1,r2       ! dim(nffn)
      real(8), dimension(:,:), allocatable :: ffnsum,ffnenu  ! dim(nffn,143)
      end module ffn_data

      subroutine read_nuclear_data(data_dir,data_desc)
!-----------------------------------------------------------------------------  
!  This routine reads, from the file netwinv, the nuclei included along 
!  with the nuclear data which will be needed for later calculations.  This 
!  data includes the atomic number, the number of protons and neutrons, and 
!  the binding energy (calculated from the tabulated mass excess).  Also the 
!  tabulations of the  partition functions, g, are read in for later 
!  interpolation. Once the set of nuclear data is read in, it is assigned 
!  to the proper nuclei.
!-----------------------------------------------------------------------------  
      use nuclear_data
      use part_funct_data
      character (LEN=*),  INTENT(in)  :: data_dir
      character (LEN=80), INTENT(out) :: data_desc
      real(8) ::  mex,a,sp 
      integer :: i,l,n,m,na,nb                                     
      integer :: it9i(24)
      character (LEN=5) nam                                 
  
!  Read in the data description
      Open(13,file=trim(data_dir)//"/net_desc",status='old')
      Read(13,"(A)") data_desc
      Open(12,FILE=trim(data_dir)//"/netwinv",status='old')
      Read(12,"(i5)") ny
  
!  Read in the partition function iteration grid, and fix endpoints
      Read(12,"(24i3)") (it9i(i),i=1,24)                   
      Do i=1,24                                                    
        t9i(i)=it9i(i)*0.01     
      Enddo                                             
      t9i(24)=t9i(24)*10.                                                       
      Do n=1,ny
        Read(12,"(a5)") nam
      Enddo                     
  
!  Set size of nuclear parameter arrays and read in nuclear parameters 
!  and partition function interpolation table. nname(0), gg(0) and angm(0) 
!  are placeholders for non-nuclei.
      Allocate (nname(0:ny))
      Allocate (aa(ny),zz(ny),nn(ny),be(ny))
      Allocate (g(24,ny),gg(0:ny),angm(0:ny))
      nname(0)=' === '
      angm(0)=0.0
      Do l=1,ny                                                  
        Read(12,"(a5,f12.3,2i4,f6.1,f10.3)") nam,a,na,nb,sp,mex
        Read(12,"(8f9.2)") (g(m,l),m=1,24)     
        aa(l)=a                                                         
        zz(l)=dble(float(na))                                                
        nn(l)=dble(float(nb))                                                
        be(l)=8.07144*nn(l)+7.28899*zz(l)-mex
        angm(l)=2.*sp+1.
        nname(l)=nam
      Enddo                                                           
!     Write(50,"(a5,4es12.4)") 
!    &      (nname(i),aa(i),zz(i),nn(i),be(i),i=1,ny)        
      Return                                                                    
      End                                                                       

      subroutine read_reaction_data(data_dir)    
!-----------------------------------------------------------------------------  
!  This routine reads in the necessary reaction data
!-----------------------------------------------------------------------------  
      use nuclear_data
      use ffn_data
      use cross_sect_data
      use reac_rate_data
      character (LEN=*) :: data_dir
      integer :: i,j,n,l
      integer :: nr,idel,jdel,ijd,ndi,nds
      character (LEN=1) :: null
      common /ceign/ nds,ndi,idel,jdel,ijd
  
!  Read in nuclear set and numbers of reactions
      Open(4,file=trim(data_dir)//"/nets4",form='unformatted',
     &       status='old')
      Read(4) ny                                                               
      Read(4) (nname(i),i=1,ny)  
      Read(4) nffn 
      Read(4) (nreac(i),i=1,3)   
  
!  If there are FFN rates, read in the FFN data and set FFN array sizes
      If(nffn>0) Then
          Open(3,FILE=trim(data_dir)//"/netweak",status='old')
          Allocate (rf(nffn),r1(nffn),r2(nffn))
          Allocate (ffnsum(nffn,143),ffnenu(nffn,143))
          Do i=1,nffn
            Read(3,"(a1)") null  
            Read(3,"(9(f8.3))") (ffnsum(i,j),ffnenu(i,j),j=1,143)  
          Enddo
      Endif
  
!  Read in the reaction cross section data
      Open(2,file=trim(data_dir)//"/nets3",form='unformatted',
     &       status='old')
  
!  Allocate and read in reaction arrays for 1 reactant reactions
      nr=nreac(1)
      Allocate (csect1(nr),rc1(7,nr),q1(nr),rpf1(nr))
      Allocate (n1i(3,nr),iwk1(nr),ires1(nr),irev1(nr),iffn(nr))
      Do j=1,nr
        Read(2) n,(n1i(l,j),l=1,3),iwk1(j),ires1(j),irev1(j),
     &          (rc1(l,j),l=1,7),q1(j)  
        If(n/=j) Then
            Write(6,*) 'Error in nets3, 1',j,n
            Exit
        Endif
      Enddo
  
!  Allocate and read in reaction arrays for 2 reactant reactions
      nr=nreac(2)
      Allocate (csect2(nr),rc2(7,nr),q2(nr),rpf2(nr),h2(nr))
      Allocate (n2i(4,nr),iwk2(nr),ires2(nr),irev2(nr))
      Do j=1,nr
        Read(2)  n,(n2i(l,j),l=1,4),iwk2(j),ires2(j),irev2(j),
     &         (rc2(l,j),l=1,7),q2(j)                           
        If(n/=j) Then
            Write(6,*) 'Error in nets3, 2',j,n
            Exit
        Endif
      Enddo
  
!  Allocate and read in reaction arrays for 3 reactant reactions
      nr=nreac(3)
      Allocate (csect3(nr),rc3(7,nr),q3(nr))
      Allocate (n3i(3,nr),iwk3(nr),ires3(nr),irev3(nr))
      Do j=1,nr
        Read(2)  n,(n3i(l,j),l=1,3),iwk3(j),ires3(j),irev3(j),
     &         (rc3(l,j),l=1,7),q3(j)                                                
        If(n/=j) Then
            Write(6,*) 'Error in nets3, 3',j,n
            Exit
        Endif
      Enddo
  
!  Allocate and read in the data linking nuclei to the reactions which 
!  affect them.  Also read in the matrix sparseness descriptors.
      Allocate (la(3,ny),le(3,ny))
      Do i=1,ny                                                             
        Read(4)  n,la(1,i),le(1,i),la(2,i),le(2,i),la(3,i),le(3,i)
        If(n.ne.i) Then
            Write(6,*) 'Error in nets4',i,n
        Endif
      Enddo
      Read(4) idel,jdel,ijd,ndi,nds
!     Write(50,*) 'ceign',idel,jdel,ijd,ndi,nds
  
!  Create and fill extended reaction->nuclei arrays
      nan(1)=le(1,ny)
      Allocate (mu1(nan(1)),a1(nan(1)),b1(nan(1)),n11(nan(1)))
      Do j=1,nan(1)
        Read( 2)  a1(j),mu1(j) 
        n11(j)=n1i(1,mu1(j))
      Enddo                                                                  
      nan(2)=le(2,ny)
      Allocate (mu2(nan(2)),a2(nan(2)),b2(nan(2)))
      Allocate (n21(nan(2)),n22(nan(2)))
      Do j=1,nan(2)
        Read( 2)  a2(j),mu2(j)  
        n21(j)=n2i(1,mu2(j))
        n22(j)=n2i(2,mu2(j))
      Enddo                                                                     
      nan(3)=le(3,ny)
      Allocate (mu3(nan(3)),a3(nan(3)),b3(nan(3)))
      Allocate (n31(nan(3)),n32(nan(3)),n33(nan(3)))
      Do j=1,nan(3)
        Read( 2)  a3(j),mu3(j)
        n31(j)=n3i(1,mu3(j))
        n32(j)=n3i(2,mu3(j))
        n33(j)=n3i(3,mu3(j))
      Enddo                                                                  
      Return                                                                    
      End                                                                       

      subroutine index_from_name(name,index)
!-----------------------------------------------------------------------
!  This subroutine takes a nuclear name and finds the corresponding
!  index for the current data set. index=0 indicates that the nuclear name
!  is not found in the current set.
!-----------------------------------------------------------------------
      use nuclear_data
      character (len=5) :: name,sname
      integer :: n,index,name_len
      name_len=len_trim(name)
      If(name_len<5) Then
          sname='     '
          sname((6-name_len):5)=name(1:name_len)
      Else
          sname=name
      Endif
      Do n=1,ny
        If(sname==nname(n)) Exit
      Enddo
      If(n>ny) Then
           index=0
      Else
           index=n
      Endif
      Return
      End

      Subroutine name_ordered(name,num,max)
!-----------------------------------------------------------------------------
!  This routine appends the number num (padded with "0"s up to the size of 
!  the number max) to the character variable name.
!-----------------------------------------------------------------------------
      character (LEN=*)   :: name
      character (LEN=9)   :: num_string,form
      character (LEN=1)   :: lmax_string
      integer             :: num,max,lnum,lmax
  
!  Find character length of max
      lmax=int(log10(float(max)))+1
      Write(unit=lmax_string,fmt='(I1)') lmax
  
!  Construct Format Spec and write num as zero padded string
      form='(I'//lmax_string//'.'//lmax_string//')'
      Write(unit=num_string,fmt=form) num
  
!  Append num_ string to name
      name=trim(name)//trim(num_string)
      Return
      End
