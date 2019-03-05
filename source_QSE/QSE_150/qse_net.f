!*******************************************************************************
!  Qse_Net 4.9 10/05/01
!  The subroutines in this file perform nucleosynthesis for a single Lagrangian
!  mass zone.  The handling of multiple zones, including multi-processing, is
!  carried out externally.
!
!  The data.f file contains the nuclear and reaction data structures and the
!  subroutines which read in the data and allocate the arrays.
!******************************************************************************


!*******************************************************************************
!  QSE_init 10/22/09
!  The subroutines in this file initalize QSE network by building the qse groups
!  sorting the reaction, and allocating variable that only need reset when the 
!  QSE groups change
!
!  The data.f file contains the nuclear and reaction data structures and the
!  subroutines which read in the data and allocate the arrays.
!*******************************************************************************

      module qse_data
!-----------------------------------------------------------------------------
!  This module contains the essential data for each included specie.
!c  aa, zz, & nn are the mass, proton and neutron numbers, and be is the
!c  binding energy.
!c  nname is the nuclear name (nname(0) is placeholder for non-nuclei).
!c  Their array sizes are set in the routine read_nuclear_data.
!c regrp1-3 Tells group membership of each reactant.
!c gp(1,2,3 or 0)number counts inside each group(1-3) or the set of singlenuc (0)
!c gp(1,2,3 or 0)i holds the number of species in each group or the # of
!c singlenuc.
!C The reactions are sorted by species first - all reactions involving n go in
!c the first set, all involving p in the second
!c all involving si28 in the 77th etc.
!c Later I need to know the N and Z of the
!c nuclues that defines each set; qnnum1-3 and qznum1-3 contian this info.
!c ,n##'s, tn##'s, mu's tmu's, la's le's-these are all tracking numbers for the
!c species involved in each reaction. They help me determine if a reaction
!c spans a group boarder and also help me sort the reactions based on this.
!c part1 dYGdot/dYR
!c part2 dYG/dYR solved 4 vertors.
!c part2i building blocks of dYG/dYR inital matrix.
!c Cqse and Cqseb Qse coeff. See qsecoeff subroutine in net.f for more info.
!c zilch - this zeros out the qsecoeff where necessey in the Jacobian build part1b.
!c f1-5 are the main list number of each of the focal nuclei.
!c-----------------------------------------------------------------------------
      use nuc_number


      integer :: numberofgroups,gn,f1,f2,f3,f4,f5
      integer :: numberofgroups_old
      integer :: adaptive
      integer, dimension(:), allocatable  :: gp1number,gp2number
      integer, dimension(:), allocatable  :: gp3number,gp0number
      integer, dimension(:), allocatable  :: gp4number
      integer  :: gp1i,gp2i,gp3i,gp0i,gp4i
      integer, dimension(:), allocatable  :: drvN,drvZ,ID
      integer, dimension(:), allocatable  :: singlenuc
      integer, dimension(:,:), allocatable  :: regrp1,regrp2,regrp3
      integer,  dimension(:),   allocatable :: qnnum,qznum
      integer,  dimension(:),   allocatable :: qnnum2,qznum2
      integer,  dimension(:),   allocatable :: qnnum3,qznum3
      integer, dimension(:), allocatable :: aaq
! Temp variables.
      integer, dimension(:,:), allocatable :: laq,leq
      integer, dimension(:,:), allocatable :: laq1,leq1,laq2,leq2
      integer, dimension(:,:), allocatable :: laq3,leq3
      integer, dimension(:), allocatable :: tmu1,tmu2,tmu3
      integer, dimension(:), allocatable :: tn11,tn21,tn22
      integer, dimension(:), allocatable :: tn31,tn32,tn33
      integer, dimension(:), allocatable :: tn1p1,tn1p2,tn1p3
      integer, dimension(:), allocatable :: tn2p1,tn2p2,tn2p3
      integer, dimension(:), allocatable :: tn3p1,tn3p2,tn3p3
      real(8),dimension(:), allocatable   ::  reldy
      real(8),  dimension(:), allocatable :: ta1,ta2,ta3,tb1,tb2,tb3 ! dim(nan(i))
      real(8), dimension(:,:), allocatable :: part1,part1b,part2b
      real(8), dimension(:,:), allocatable :: part2i,part2,part2ib
      real(8), dimension(:), allocatable :: dlCqse
!     Jacobian stuff Note in xnet 5 these are in jacobi_---.f
      real(8), dimension(:), allocatable :: f,ydot0,um
      real(8), dimension(:,:), allocatable :: am ! the Jacobian Matrix
      real(8), dimension(:,:), allocatable :: dydot ! Part 1 ofJ Matrix
      integer, dimension(:), allocatable :: indx,indx4


      end module qse_data



      subroutine group_number(data_dir,t9t)
!c-----------------------------------------------------
! This subroutine sets up the nubmer of groups.
!c-----------------------------------------------------
      use qse_data
      use controls
      use abundances
      use nuclear_data
      character (LEN=*) :: data_dir
      integer :: nsi,zsi,nfe,zfe,ncr,zcr,dummy2

         numberofgroups_old= numberofgroups
      If (adaptive>0)then

      	 If(t9t.ge.6)then
             dummy2=1
         Elseif(t9t<6.and.t9t>4)then
            Write(*,*) "group2"
             dummy2=2
       	 Elseif(t9t<4)then 
            dummy2=3
         Endif

         If (numberofgroups_old.ne.dummy2)then
            
               gp0i=0
               gp1i=0
               gp2i=0
               gp3i=0
            Deallocate(indx,f,ydot0,am,dydot,indx4)
            Deallocate(ydot,dy,reldy,yqo,yq,yqt)
            If(numberofgroups==1)then
               Deallocate(aaq,gp0number,gp1number)
            ElseIf(numberofgroups==2)then
               Deallocate(aaq,gp0number,gp1number,gp2number)
       	    Elseif(numberofgroups==3)then
               Deallocate(aaq,gp0number,gp1number)
               Deallocate(gp2number,gp3number)
             Endif
          Endif
            numberofgroups=dummy2

       Else
       Endif
       
      gn=numberofgroups+1
!c Read in singlenuc
       If(adaptive==0.or.numberofgroups_old.ne.dummy2)then  
          If (numberofgroups==1)then
            Write(*,*) "group1"
             call group1
             call qse_sort
           ElseIf (numberofgroups==2)then
            Write(*,*) "group2"
             call group2
             write(*,*) "1"
             call qse_sort
             write(*,*) "2"
!            call read_match_data(data_dir)
             write(*,*) "3"
!            call flux_init
             write(*,*) "4"
           ElseIf (numberofgroups==3)then
             call group3
            Write(*,*) "group3"
            ! write(*,*) " I call group 3",gp0i,gp1i,gp2i,gp3i
             call qse_sort
!            call read_match_data(data_dir)
!            call flux_init
           Endif
        
            If(.not.allocated(f)) Allocate(indx(gp0i), &
     &         f(gp0i),ydot0(ny))      
            If(.not.allocated(am)) Allocate(am(gp0i,gp0i) &
     &         ,dydot(gp0i,gp0i)) 
            Allocate(indx4(gn:gn))
            Allocate (ydot(gp0i),dy(gp0i),reldy(gp0i))
            Allocate (yqo(gp0i),yq(gp0i),yqt(gp0i))
             write(*,*) "5"
        Endif 

      If (adaptive>0.and.numberofgroups_old.ne.dummy2)then
! Build QSE coeff. from 4 focal
          call qse_coeff !
       write(*,*) "6"
!repalce initial abund. with group members.
          call qse_make(yt)
       write(*,*) "7"
!Solve for QSE this is step 0.
          call update(yt) ! y is getting screwed up
       write(*,*) "8"
          ye=sum(zz*yt)
       write(*,*) "9"
           yq =0.0
       write(*,*) "9"

           do i=1,gp1i
              yq(1) =yq(1) + nn(i)*yt(i)
              yq(2) =yq(2)+  zz(i)*yt(i)
            enddo
       write(*,*) "9"
            do ig=1,gp2i
               i=gp2number(ig)
               nsi=nn(i)-nn(f3)
               zsi=zz(i)-zz(f3)
               yq(1)=yq(1)+nsi*yt(i)
               yq(2)=yq(2)+zsi*yt(i)
               yq(3)=yq(3)+yt(i)
             enddo
       write(*,*) "9", gp3i
             do ig=1,gp3i
               i=gp3number(ig)
               nfe=nn(i)-nn(f4)
               zfe=zz(i)-zz(f4)
               yq(1)=yq(1)+nfe*yt(i)
               yq(2)=yq(2)+zfe*yt(i)
               yq(4)=yq(4)+yt(i)
             enddo
       write(*,*) "9"
             do ig=gn+1,gp0i
                i=gp0number(ig)
                yq(ig)=yt(i)
             enddo
       write(*,*) "9"

       write(*,*) "9"
                yqt=yq
          call qse_coeff ! must be called every time rho or T  changes
! Updates groups and puts
          call qse_make(yt)
        Endif

       write(*,*) "6"

       Return
       End




      subroutine group1
!c----------------------------------------------------c
!c This routine selects the group members for qse with
!c two groups.
!c----------------------------------------------------c
      use nuclear_data
      use qse_data
      use ffn_data
      use cross_sect_data
      integer :: i,j,k,n,l,m, dummysngl(ny)
      integer ::  dummysi(ny),dummyfe(ny)
        j=0
        k=3
        l=0
        dummysi=0
! Set the focal nuclei.
!Perhaps this could be done automaticaly by
! looking through the group lists to find the largest abundance. I need to do
! a systematic study of what makes a good focal nuc. The follwoing are
! chosen based on tradition from the literature.
       f1=1 !p
       f2=2 !n

       gp1i=ny
        gp0i=2
        gp2i=0
        gp3i=0
        gp4i=0
       Allocate(aaq(gp0i))
       Allocate (gp0number(0:gp0i))
       Allocate (gp1number(0:gp1i))

        Do i=1,ny
!       Set up the light nuc.
!        if (i<=gp1i)then
           gp1number(i)=i
!        Endif
        Enddo
! The focal nuc.

!--------------------------
         gp0number(1)=f1 !p
         gp0number(2)=f2 !n
!---------------------------

        drvN=0
        drvZ=0
        ID=0

        Do i=1,gp1i
           i1=gp1number(i)
            drvN(i1)=nn(i1)
            drvZ(i1)=zz(i1)
            ID(i1)=2
        Enddo
! sets atomic mass for testc2 and xtot in the qse_net convergence test
        aaq(1)=aa(1)
        aaq(2)=aa(2)
!        zzq(1)=zz(1)
!        zzq(2)=zz(2)
!        nnq(1)=nn(1)
!        nnq(2)=nn(2)

       Return
       End !subroutine group1

      subroutine group2_alpharich
!c----------------------------------------------------c
!c This routine selects the group members for qse with
!c two groups.
!c----------------------------------------------------c
      use nuclear_data
      use qse_data
      use ffn_data
      use cross_sect_data
      integer :: i,j,k,n,l,m, dummysngl(ny)
      integer ::  dummysi(ny),dummyfe(ny)
        j=0
        k=3
        l=0
        dummyfe=0
! Set the focal nuclei.
!Perhaps this could be done automaticaly by
! looking through the group lists to find the largest abundance. I need to do
! a systematic study of what makes a good focal nuc. The follwoing are
! chosen based on tradition from the literature.
       f1=1 !p
       f2=2 !n
       f3=186 !Cr52

       gp1i=6
       Allocate (gp1number(0:gp1i))

        Do i=1,ny
!       Set up the light nuc.
        if (i<=gp1i)then
           gp1number(i)=i
        Endif
!Setup the single nuc
        if(i>6.and.zz(i)<=22)then ! uses "i>6" because the 1st 6 are lt .
             k=k+1
             dummysngl(k)=i
         Elseif(zz(i)==23.and.nn(i)<=28)then
            k=k+1
            dummysngl(k)=i
         Elseif(zz(i)==24.and.nn(i)<=27)then
            k=k+1
            dummysngl(k)=i
         Elseif(zz(i)==25.and.nn(i)<=26)then
            k=k+1
            dummysngl(k)=i
         Elseif(zz(i)==26.and.nn(i)<=26)then
            k=k+1
            dummysngl(k)=i
         Endif
! Setup the Fe group
         If(zz(i)==23.and.nn(i)>=29)then !Sc
            l=l+1
            dummyfe(l)=i
         Elseif(zz(i)==24.and.nn(i)>=28)then !Ti-Ge
            l=l+1
            dummyfe(l)=i
         Elseif(zz(i)==25.and.nn(i)>=27)then !Ti-Ge
            l=l+1
            dummyfe(l)=i
         Elseif(zz(i)==26.and.nn(i)>=27)then !Ti-Ge
            l=l+1
            dummyfe(l)=i
         Elseif(zz(i)>=27)then !Ti-Ge
            l=l+1
            dummyfe(l)=i
         Endif
         Enddo !i
       gp0i=k
       gp2i=l
       Allocate (aaq(gp0i))
       Allocate (gp0number(0:gp0i))
       Allocate (gp2number(0:gp2i))

! The focal nuc.
!--------------------------
         gp0number(1)=f1 !p
         gp0number(2)=f2 !n
         gp0number(3)=f3 !si2
!---------------------------

        drvN=0
        drvZ=0
        ID=0

      Do i=gn+1,gp0i
          gp0number(i)=dummysngl(i)
       enddo
      Do i=1,gp2i
          gp2number(i)=dummyfe(i)
         i3=gp2number(i)
           drvN(i3)=nn(i3)-nn(f3)
           drvZ(i3)=zz(i3)-zz(f3)
           ID(i3)=3
       enddo
        Do i=1,gp1i
           i1=gp1number(i)
            drvN(i1)=nn(i1)
            drvZ(i1)=zz(i1)
            ID(i1)=2
        Enddo

        Do i=1,gp0i
           Do j=1,ny
             ig=gp0number(i)
              If(j==ig)then
                 aaq(i)=aa(j)
              Exit
             Endif
           Enddo
        Enddo


       Return
       End

!-------------------------
      subroutine group2_ext
!c----------------------------------------------------c
!c This routine selects the group members for qse with
!c two groups.
!c----------------------------------------------------c
      use nuclear_data
      use qse_data
      use ffn_data
      use cross_sect_data
      integer :: i,j,k,n,l,m, dummysngl(ny)
      integer ::  dummysi(ny),dummyfe(ny)
        j=0
        k=3
        l=0
        dummysi=0
! Set the focal nuclei.
!Perhaps this could be done automaticaly by
! looking through the group lists to find the largest abundance. I need to do
! a systematic study of what makes a good focal nuc. The follwoing are
! chosen based on tradition from the literature.
       f1=1 !p
       f2=2 !n
       f3=19 ! si28

       gp1i=3
       Allocate (gp1number(0:gp1i))

        Do i=1,ny
!       Set up the light nuc.
        if (i<=gp1i)then
           gp1number(i)=i
        Endif
!Setup the single nuc
        if(i>3.and.zz(i)<10)then ! uses "i>6" because the 1st 6 are lt .
             k=k+1
             dummysngl(k)=i
!           write(502,*) k,nname(dummysngl(k))
        Elseif(zz(i)==12.and.nn(i)<=10)then
             k=k+1
             dummysngl(k)=i
        Elseif(zz(i)==11.and.nn(i)<=12)then
             k=k+1
             dummysngl(k)=i
        Elseif(zz(i)==10.and.nn(i)<=13)then
             k=k+1
             dummysngl(k)=i
        endif

! Set up the Si group From mg to ge
        If(zz(i)==10.and.nn(i)>=14)then
            j=j+1
           dummysi(j)=i
        ElseIf(zz(i)==11.and.nn(i)>=13)then
            j=j+1
           dummysi(j)=i
        ElseIf(zz(i)==12.and.nn(i)>=11)then !mg
           j=j+1
           dummysi(j)=i
         ElseIf(zz(i)>=13)then
            j=j+1
            dummysi(j)=i
         Endif
         Enddo !i
       gp0i=k
       gp2i=j
       Allocate (gp0number(0:gp0i))
       Allocate (gp2number(0:gp2i))

! The focal nuc.
!--------------------------
         gp0number(1)=f1 !p
         gp0number(2)=f2 !n
         gp0number(3)=f3 !si2
!---------------------------

        drvN=0
        drvZ=0
        ID=0

      Do i=gn+1,gp0i
          gp0number(i)=dummysngl(i)
       enddo
      Do i=1,gp2i
          gp2number(i)=dummysi(i)
         i3=gp2number(i)
           drvN(i3)=nn(i3)-nn(f3)
           drvZ(i3)=zz(i3)-zz(f3)
           ID(i3)=3
       enddo
        Do i=1,gp1i
           i1=gp1number(i)
            drvN(i1)=nn(i1)
            drvZ(i1)=zz(i1)
            ID(i1)=2
        Enddo

       Allocate(aaq(gp0i))
       Do i=1,gp0i
           Do j=1,ny
             ig=gp0number(i)
              If(j==ig)then
                 aaq(i)=aa(j)
              Exit
             Endif
           Enddo
        Enddo

       Return
       End
      subroutine group2
!c----------------------------------------------------c
!c This routine selects the group members for qse with
!c two groups.
!c----------------------------------------------------c
      use nuclear_data
      use qse_data
      use ffn_data
      use cross_sect_data
      integer :: i,j,k,n,l,m, dummysngl(ny)
      integer ::  dummysi(ny),dummyfe(ny)
        j=0
        k=3
        l=0
        dummysi=0
! Set the focal nuclei.
!Perhaps this could be done automaticaly by
! looking through the group lists to find the largest abundance. I need to do
! a systematic study of what makes a good focal nuc. The follwoing are
! chosen based on tradition from the literature.
       f1=1 !p
       f2=2 !n
       f3=74 ! si28


       gp1i=6
       Allocate (gp1number(0:gp1i))
       gp1number=0
        Do i=1,ny
!       Set up the light nuc.
        if (i<=gp1i)then
           gp1number(i)=i
        Endif
!Setup the single nuc
        if(i>6.and.zz(i)<12)then ! uses "i>6" because the 1st 6 are lt .
             k=k+1
             dummysngl(k)=i
!           write(502,*) k,nname(dummysngl(k))
        Elseif(zz(i)==12.and.nn(i)<=10)then
             k=k+1
             dummysngl(k)=i
        endif

! Set up the Si group From mg to ge

        If(zz(i)==12.and.nn(i)>=11)then !mg
           j=j+1
           dummysi(j)=i
         ElseIf(zz(i)>=13)then
            j=j+1
            dummysi(j)=i
         Endif
         Enddo !i
       gp0i=k
       gp2i=j
       Allocate (gp0number(0:gp0i))
       Allocate (gp2number(0:gp2i))
! The focal nuc.
!--------------------------
         gp0number(1)=f1 !p
         gp0number(2)=f2 !n
         gp0number(3)=f3 !si28
!---------------------------

        drvN=0
        drvZ=0
        ID=0

      Do i=gn+1,gp0i
          gp0number(i)=dummysngl(i)
       enddo
      Do i=1,gp2i
          gp2number(i)=dummysi(i)
         i3=gp2number(i)
           drvN(i3)=nn(i3)-nn(f3)
           drvZ(i3)=zz(i3)-zz(f3)
           ID(i3)=3
       enddo
        Do i=1,gp1i
           i1=gp1number(i)
            drvN(i1)=nn(i1)
            drvZ(i1)=zz(i1)
            ID(i1)=2
        Enddo
       Allocate(aaq(gp0i))
       Do i=1,gp0i
           Do j=1,ny
             ig=gp0number(i)
              If(j==ig)then
                aaq(i)=aa(j)
              Exit
             Endif
           Enddo
        Enddo
       Return
       End !subroutine group2
      subroutine group3_smallest
!c----------------------------------------------------c
!c This routine selects the group members for qse with
!c three groups.
!c----------------------------------------------------c
      use nuclear_data
      use qse_data
      use ffn_data
      use cross_sect_data
      integer :: i,j,k,n,l,m, dummysngl(ny)
      integer ::  dummysi(ny),dummyfe(ny)
        j=0
        k=4
        l=0
        dummysi=0
! Set the focal nuclei.
!Perhaps this could be done automaticaly by
! looking through the group lists to find the largest abundance. I need to do
! a systematic study of what makes a good focal nuc. The follwoing are
! chosen based on tradition from the literature.
       f1=1 !p
       f2=2 !n
       f3=19 !si28
       f4=74  !fe56

       gp1i=3
       Allocate (gp1number(0:gp1i))

        Do i=1,ny
!       Set up the light nuc.
        if (i<=gp1i)then
           gp1number(i)=i
        Endif

!Setup the single nuc
        if(i>3.and.zz(i)<12)then ! uses "i>6" because the 1st 6 are lt .
             k=k+1
             dummysngl(k)=i
!           write(502,*) k,nname(dummysngl(k))
        Elseif(zz(i)==12.and.nn(i)<=10)then
             k=k+1
             dummysngl(k)=i

! 3 or 4 groups
         Elseif(zz(i)==18.and.nn(i)==26)then
            k=k+1
            dummysngl(k)=i
         Elseif(zz(i)==19.and.nn(i)>=24)then
            k=k+1
            dummysngl(k)=i
         Elseif(zz(i)==20.and.nn(i)>=22)then !.and.nn(i)<=27)then
            k=k+1
            dummysngl(k)=i
         Elseif(zz(i)==21.and.nn(i)>=21.and.nn(i)<=26)then
            k=k+1
            dummysngl(k)=i
         Elseif(zz(i)==22.and.nn(i)>=21.and.nn(i)<=26)then !ti
            k=k+1
            dummysngl(k)=i
         Elseif(zz(i)==23.and.nn(i)<=26)then !ti
            k=k+1
            dummysngl(k)=i
         Elseif(zz(i)>=24.and.zz(i)<=25.and.nn(i)>=21.
     & and.nn(i)<=23)then
            k=k+1
            dummysngl(k)=i
         Endif

! Set up the Si group

        If(zz(i)==12.and.nn(i)>=11)then
           j=j+1
           dummysi(j)=i
!           write(501,*) j,nname(dummysi(j)),'begin'
         ElseIf(zz(i)>=13.and.nn(i)>=10.and.zz(i)<=18.and.
     &      nn(i)<=25)then
            j=j+1
            dummysi(j)=i
         Elseif(zz(i)==19.and.nn(i)<=23)then
            j=j+1
            dummysi(j)=i
         Elseif(zz(i)==20.and.nn(i)<=21)then
            j=j+1
            dummysi(j)=i
         Elseif(zz(i)==21.and.nn(i)<=20)then
            j=j+1
            dummysi(j)=i
         Elseif(zz(i)==22.and.nn(i)==20)then
            j=j+1
            dummysi(j)=i
         Endif
! Setup the Fe group
!         if(zz(i)==20.and.nn(i)>=28)then !Ca
!            l=l+1
!            dummyfe(l)=i
         if(zz(i)==21.and.nn(i)>=27)then !Sc
            l=l+1
            dummyfe(l)=i
         Elseif(zz(i)==22.and.nn(i)>=27)then !Ti
            l=l+1
            dummyfe(l)=i
         Elseif(zz(i)==23.and.nn(i)>=27)then !V
            l=l+1
            dummyfe(l)=i
         Elseif(zz(i)>=24.and.nn(i)>=24)then !cr-Ge
            l=l+1
            dummyfe(l)=i
         Endif
         Enddo !i
       gp0i=k
       gp2i=j
       gp3i=l
       Allocate (gp0number(0:gp0i))
       Allocate (gp2number(0:gp2i))
       Allocate (gp3number(0:gp3i))
! The four focal nuc.
!--------------------------
         gp0number(1)=f1
         gp0number(2)=f2
         gp0number(3)=f3
         gp0number(4)=f4
!        qn(1)=0
!        qn(2)=0
!        qn(3)=nn(f3)
!        qn(4)=nn(f4)
!        qz(1)=0
!        qz(2)=0
!        qz(3)=zz(f3)
!        qz(4)=zz(f4)
!---------------------------

        drvN=0
        drvZ=0
        ID=0
       Do i=gn+1,gp0i
          gp0number(i)=dummysngl(i)
!            write(501,*) i,nname(gp0number(i))
       enddo
      Do i=1,gp2i
          gp2number(i)=dummysi(i)
           i3=gp2number(i)
           drvN(i3)=nn(i3)-nn(f3)
           drvZ(i3)=zz(i3)-zz(f3)
           ID(i3)=3
       enddo
      Do i=1,gp3i
          gp3number(i)=dummyfe(i)
           i4=gp3number(i)
         drvN(i4)=nn(i4)-nn(f4)
         drvZ(i4)=zz(i4)-zz(f4)
         ID(i4)=4
       enddo
        Do i=1,gp1i
           i1=gp1number(i)
            drvN(i1)=nn(i1)
            drvZ(i1)=zz(i1)
            ID(i1)=2
        Enddo

       Allocate(aaq(gp0i))
       Do i=1,gp0i
           Do j=1,ny
             ig=gp0number(i)
              If(j==ig)then
                 aaq(i)=aa(j)
              Exit
             Endif
           Enddo
        Enddo

!        Do i=1,ny
!         write(510,*) i, drvN(i),drvZ(i)
!        Enddo


       Return
       End ! subroutine group3

!----------------------
      subroutine group3_small
!c----------------------------------------------------c
!c This routine selects the group members for qse with
!c three groups.
!c----------------------------------------------------c
      use nuclear_data
      use qse_data
      use ffn_data
      use cross_sect_data
      integer :: i,j,k,n,l,m, dummysngl(ny)
      integer ::  dummysi(ny),dummyfe(ny)
        j=0
        k=4
        l=0
        dummysi=0
! Set the focal nuclei.
!Perhaps this could be done automaticaly by
! looking through the group lists to find the largest abundance. I need to do
! a systematic study of what makes a good focal nuc. The follwoing are
! chosen based on tradition from the literature.
       f1=1 !p
       f2=2 !n
       f3=19 !si28
       f4=74  !fe56

       gp1i=3
       Allocate (gp1number(0:gp1i))

        Do i=1,ny
!       Set up the light nuc.
        if (i<=gp1i)then
           gp1number(i)=i
        Endif

!Setup the single nuc
        if(i>3.and.zz(i)<12)then ! uses "i>6" because the 1st 6 are lt .
             k=k+1
             dummysngl(k)=i
!           write(502,*) k,nname(dummysngl(k))
        Elseif(zz(i)==12.and.nn(i)<=10)then
             k=k+1
             dummysngl(k)=i

! 3 or 4 groups
         Elseif(zz(i)==18.and.nn(i)==26)then
            k=k+1
            dummysngl(k)=i
         Elseif(zz(i)==19.and.nn(i)>=24)then
            k=k+1
            dummysngl(k)=i
         Elseif(zz(i)==20.and.nn(i)>=22)then !.and.nn(i)<=27)then
            k=k+1
            dummysngl(k)=i
         Elseif(zz(i)==21.and.nn(i)>=21.and.nn(i)<=26)then
            k=k+1
            dummysngl(k)=i
         Elseif(zz(i)==22.and.nn(i)>=21.and.nn(i)<=26)then !ti
            k=k+1
            dummysngl(k)=i
         Elseif(zz(i)==23.and.nn(i)<=26)then !ti
            k=k+1
            dummysngl(k)=i
         Elseif(zz(i)>=24.and.zz(i)<=25.and.nn(i)>=21.
     & and.nn(i)<=23)then
            k=k+1
            dummysngl(k)=i
         Endif

! Set up the Si group

        If(zz(i)==12.and.nn(i)>=11)then
           j=j+1
           dummysi(j)=i
!           write(501,*) j,nname(dummysi(j)),'begin'
         ElseIf(zz(i)>=13.and.nn(i)>=10.and.zz(i)<=18.and.
     &      nn(i)<=25)then
            j=j+1
            dummysi(j)=i
         Elseif(zz(i)==19.and.nn(i)<=23)then
            j=j+1
            dummysi(j)=i
         Elseif(zz(i)==20.and.nn(i)<=21)then
            j=j+1
            dummysi(j)=i
         Elseif(zz(i)==21.and.nn(i)<=20)then
            j=j+1
            dummysi(j)=i
         Elseif(zz(i)==22.and.nn(i)==20)then
            j=j+1
            dummysi(j)=i
         Endif
! Setup the Fe group
!         if(zz(i)==20.and.nn(i)>=28)then !Ca
!            l=l+1
!            dummyfe(l)=i
         if(zz(i)==21.and.nn(i)>=27)then !Sc
            l=l+1
            dummyfe(l)=i
         Elseif(zz(i)==22.and.nn(i)>=27)then !Ti
            l=l+1
            dummyfe(l)=i
         Elseif(zz(i)==23.and.nn(i)>=27)then !V
            l=l+1
            dummyfe(l)=i
         Elseif(zz(i)>=24.and.nn(i)>=24)then !cr-Ge
            l=l+1
            dummyfe(l)=i
         Endif
         Enddo !i
       gp0i=k
       gp2i=j
       gp3i=l
       Allocate (gp0number(0:gp0i))
       Allocate (gp2number(0:gp2i))
       Allocate (gp3number(0:gp3i))
! The four focal nuc.
!--------------------------
         gp0number(1)=f1
         gp0number(2)=f2
         gp0number(3)=f3
         gp0number(4)=f4
!---------------------------

        drvN=0
        drvZ=0
        ID=0
       Do i=gn+1,gp0i
          gp0number(i)=dummysngl(i)
!            write(501,*) i,nname(gp0number(i))
       enddo
      Do i=1,gp2i
          gp2number(i)=dummysi(i)
           i3=gp2number(i)
           drvN(i3)=nn(i3)-nn(f3)
           drvZ(i3)=zz(i3)-zz(f3)
           ID(i3)=3
       enddo
      Do i=1,gp3i
          gp3number(i)=dummyfe(i)
           i4=gp3number(i)
         drvN(i4)=nn(i4)-nn(f4)
         drvZ(i4)=zz(i4)-zz(f4)
         ID(i4)=4
       enddo
        Do i=1,gp1i
           i1=gp1number(i)
            drvN(i1)=nn(i1)
            drvZ(i1)=zz(i1)
            ID(i1)=2
        Enddo

       Allocate(aaq(gp0i))
       Do i=1,gp0i
           Do j=1,ny
             ig=gp0number(i)
              If(j==ig)then
                 aaq(i)=aa(j)
              Exit
             Endif
           Enddo
        Enddo

!        Do i=1,ny
!         write(510,*) i, drvN(i),drvZ(i)
!        Enddo


       Return
       End ! subroutine group3

      subroutine group3
!c----------------------------------------------------c
!c This routine selects the group members for qse with
!c three groups.
!c----------------------------------------------------c
      use nuclear_data
      use qse_data
      use ffn_data
      use cross_sect_data
      integer :: i,j,k,n,l,m, dummysngl(ny)
      integer ::  dummysi(ny),dummyfe(ny)
        j=0
        k=4
        l=0
        dummysi=0
! Set the focal nuclei.
!Perhaps this could be done automaticaly by
! looking through the group lists to find the largest abundance. I need to do
! a systematic study of what makes a good focal nuc. The follwoing are
! chosen based on tradition from the literature.
       f1=1 !p
       f2=2 !n
       f3=77 !si28
       f4=186 !fe56

       gp1i=6
       Allocate (gp1number(0:gp1i))

        Do i=1,ny
!       Set up the light nuc.
        if (i<=gp1i)then
           gp1number(i)=i
        Endif
!Setup the single nuc
        if(i>6.and.zz(i)<12)then ! uses "i>6" because the 1st 6 are lt .
             k=k+1
             dummysngl(k)=i
!           write(502,*) k,nname(dummysngl(k))
        Elseif(zz(i)==12.and.nn(i)<=10)then
             k=k+1
             dummysngl(k)=i

! 3 or 4 groups
         Elseif(zz(i)==18.and.nn(i)==26)then
            k=k+1
            dummysngl(k)=i
         Elseif(zz(i)==19.and.nn(i)>=24)then
            k=k+1
            dummysngl(k)=i
         Elseif(zz(i)==20.and.nn(i)>=22.and.nn(i)<=27)then
            k=k+1
            dummysngl(k)=i
         Elseif(zz(i)==21.and.nn(i)>=21.and.nn(i)<=25)then
            k=k+1
            dummysngl(k)=i
         Elseif(zz(i)>=22.and.zz(i)<=25.and.nn(i)>=21.
     & and.nn(i)<=23)then
            k=k+1
            dummysngl(k)=i
         Endif

! Set up the Si group

        If(zz(i)==12.and.nn(i)>=11)then
           j=j+1
           dummysi(j)=i
!           write(501,*) j,nname(dummysi(j)),'begin'
         ElseIf(zz(i)>=13.and.nn(i)>=10.and.zz(i)<=18.and.
     &      nn(i)<=25)then
            j=j+1
            dummysi(j)=i
         Elseif(zz(i)==19.and.nn(i)<=23)then
            j=j+1
            dummysi(j)=i
         Elseif(zz(i)==20.and.nn(i)<=21)then
            j=j+1
            dummysi(j)=i
         Elseif(zz(i)==21.and.nn(i)<=20)then
            j=j+1
            dummysi(j)=i
         Elseif(zz(i)==22.and.nn(i)==20)then
            j=j+1
            dummysi(j)=i
         Endif
! Setup the Fe group
         if(zz(i)==20.and.nn(i)>=28)then !Ca
            l=l+1
            dummyfe(l)=i
         Elseif(zz(i)==21.and.nn(i)>=26)then !Sc
            l=l+1
            dummyfe(l)=i
         Elseif(zz(i)>=22.and.nn(i)>=24)then !Ti-Ge
            l=l+1
            dummyfe(l)=i
         Endif

         Enddo !i
       gp0i=k
       gp2i=j
       gp3i=l
       Allocate (gp0number(0:gp0i))
       Allocate (gp2number(0:gp2i))
       Allocate (gp3number(0:gp3i))
! The four focal nuc.
!--------------------------
         gp0number(1)=f1
         gp0number(2)=f2
         gp0number(3)=f3
         gp0number(4)=f4
!---------------------------

        drvN=0
        drvZ=0
        ID=0
       Do i=gn+1,gp0i
          gp0number(i)=dummysngl(i)
!            write(501,*) i,nname(gp0number(i))
       enddo
      Do i=1,gp2i
          gp2number(i)=dummysi(i)
           i3=gp2number(i)
           drvN(i3)=nn(i3)-nn(f3)
           drvZ(i3)=zz(i3)-zz(f3)
           ID(i3)=3
       enddo
      Do i=1,gp3i
          gp3number(i)=dummyfe(i)
           i4=gp3number(i)
         drvN(i4)=nn(i4)-nn(f4)
         drvZ(i4)=zz(i4)-zz(f4)
         ID(i4)=4
       enddo
        Do i=1,gp1i
           i1=gp1number(i)
            drvN(i1)=nn(i1)
            drvZ(i1)=zz(i1)
            ID(i1)=2
        Enddo
       Allocate(aaq(gp0i))
       Do i=1,gp0i
           Do j=1,ny
             ig=gp0number(i)
              If(j==ig)then
                 aaq(i)=aa(j)
              Exit
             Endif
           Enddo
        Enddo

       Return
       End ! subroutine group3

      Subroutine group4
!c----------------------------------------------------c
!c This routine selects the group members for qse with
!c four groups.
!c----------------------------------------------------c
      use nuclear_data
      use qse_data
      use ffn_data
      use cross_sect_data
      integer :: i,j,k,n,l,m, dummysngl(ny)
      integer ::  dummysi(ny),dummyfe(ny),dummycr(ny)
        j=0
        k=5
        l=0
        m=0
        dummysi=0
! Set the focal nuclei.
!Perhaps this could be done automaticaly by
! looking through the group lists to find the largest abundance. I need to do
! a systematic study of what makes a good focal nuc. The follwoing are
! chosen based on tradition from the literature.
       f1=1 ! n
       f2=2 !p
       f3=77 !si28
       f4=210 !fe56
       f5=184 !Cr48
       gp1i=6
       Allocate (gp1number(0:gp1i))

        Do i=1,ny
!       Set up the light nuc.
        if (i<=gp1i)then
           gp1number(i)=i
        Endif
!Setup the single nuc
        if(i>6.and.zz(i)<12)then ! uses "i>6" because the 1st 6 are lt .
             k=k+1
             dummysngl(k)=i
         Elseif(zz(i)==12.and.nn(i)<=10)then
             k=k+1
             dummysngl(k)=i
! 3 or 4 groups
         Elseif(zz(i)==18.and.nn(i)==26)then
            k=k+1
            dummysngl(k)=i
         Elseif(zz(i)==19.and.nn(i)>=24)then
            k=k+1
            dummysngl(k)=i
         Elseif(zz(i)==20.and.nn(i)>=22.and.nn(i)<=27)then
            k=k+1
            dummysngl(k)=i
         Elseif(zz(i)==21.and.nn(i)>=21.and.nn(i)<=25)then
            k=k+1
            dummysngl(k)=i
         Elseif(zz(i)>=22.and.zz(i)<=25.and.nn(i)>=21.
     & and.nn(i)<=23)then
            k=k+1
            dummysngl(k)=i
         Endif

! Set up the Si group

        If(zz(i)==12.and.nn(i)>=11)then
           j=j+1
           dummysi(j)=i
         ElseIf(zz(i)>=13.and.nn(i)>=10.and.zz(i)<=18.and.
     &      nn(i)<=25)then
            j=j+1
            dummysi(j)=i
!            write(502,*) j,nname(dummysi(j))
         Elseif(zz(i)==19.and.nn(i)<=23)then
            j=j+1
            dummysi(j)=i
         Elseif(zz(i)==20.and.nn(i)<=21)then
            j=j+1
            dummysi(j)=i
         Elseif(zz(i)==21.and.nn(i)<=20)then
            j=j+1
            dummysi(j)=i
         Elseif(zz(i)==22.and.nn(i)==20)then
            j=j+1
            dummysi(j)=i
         Endif
! Setup the Cr group
         if(zz(i)==20.and.nn(i)>=28)then !Ca
            l=l+1
            dummycr(l)=i
         Elseif(zz(i)==21.and.nn(i)>=26)then !Sc
            l=l+1
         Elseif(zz(i)==22.and.nn(i)>=24)then !Ti
            l=l+1
            dummycr(l)=i
         Elseif(zz(i)==23.and.nn(i)>=24)then !Ti
            l=l+1
            dummycr(l)=i
         Elseif(zz(i)>=24.and.zz(i)<=25.and.
     & nn(i)>=24.and.nn(i)<=27)then !Cr and Mn
            l=l+1
            dummycr(l)=i
         Elseif(zz(i)>=26.and.zz(i)<=28.and.nn(i)<=26)then
            l=l+1
            dummycr(l)=i
         Endif
!Setup fe groups
         if(zz(i)>=24.and.zz(i)<=25.and.nn(i)>=28)then
           m=m+1
           dummyfe(m)=i
         elseif(zz(i)>=26.and.nn(i)>=27)then
           m=m+1
           dummyfe(m)=i
         Endif


         Enddo !i
       gp0i=k
       gp2i=j
       gp3i=m
       gp4i=l
       Allocate (gp0number(0:gp0i))
       Allocate (gp2number(0:gp2i))
       Allocate (gp3number(0:gp3i))
       Allocate (gp4number(0:gp4i))
! The four focal nuc.
!--------------------------
         gp0number(1)=f1
         gp0number(2)=f2
         gp0number(3)=f3
         gp0number(4)=f4
         gp0number(5)=f5
!---------------------------

        drvN=0
        drvZ=0
        ID=0

      Do i=gn+1,gp0i
          gp0number(i)=dummysngl(i)
       enddo
      Do i=1,gp2i
          gp2number(i)=dummysi(i)
           i3=gp2number(i)
           drvN(i3)=nn(i3)-nn(f3)
           drvZ(i3)=zz(i3)-zz(f3)
           ID(i3)=3

       enddo
      Do i=1,gp3i
          gp3number(i)=dummyfe(i)
          i4=gp3number(i)
         drvN(i4)=nn(i4)-nn(f4)
         drvZ(i4)=zz(i4)-zz(f4)
         ID(i4)=4
         Write(993,*) i,gp3number(i),nname(gp3number(i))
       enddo
      Do i=1,gp4i
          gp4number(i)=dummycr(i)
          i5=gp4number(i)
         drvN(i5)=nn(i5)-nn(f5)
         drvZ(i5)=zz(i5)-zz(f5)
         ID(i5)=5
         Write(994,*) i,gp4number(i),nname(gp4number(i))
       Enddo
       Do i=1,gp1i
           i1=gp1number(i)
            drvN(i1)=nn(i1)
            drvZ(i1)=zz(i1)
            ID(i1)=2
        Enddo

       Return
       End ! subroutine group4

      subroutine qread_reaction_data(data_dir)

!!****************************************************************************
!! This routine is nearly identical to xnet's read_reaction_data. The
! differnce comes from keeping track of all of the reactants so the 
! reactions can be later sorted. xnet does not keep tral of all reactancts.
!c-----------------------------------------------------------------------------
      use nuclear_data
      use qse_data
      use reac_rate_data
      use ffn_data
      use cross_sect_data
!      use reac_rate_data
      integer :: i,j,n,l,ki,ig,cn,cm,cn2,cn3
      integer :: i0,la1,le1,la2,le2,la3,le3,j1,l1,l2,l3
       character(len=5)  :: nnname
      character (LEN=1) :: null
      character (LEN=*) :: data_dir
      integer :: nr,idel,jdel,ijd,ndi,nds,plotindex(ny)
      common /ceign/ nds,ndi,idel,jdel,ijd
      j=0

!c-----
!c  Read in nuclear set and numbers of reactions
      Open(4,file=trim(data_dir)//"/nets4",form='unformatted',
     &       status='old')
      Read(4) ny
      Read(4) (nname(i),i=1,ny)
      Read(4) nffn
      Read(4) (nreac(i),i=1,3)
!c  If there are FFN rates, read in the FFN data and set FFN array sizes
      If(nffn>0) Then
          Open(3,FILE=trim(data_dir)//"/netweak",status='old')
          Allocate (rf(nffn),r1(nffn),r2(nffn))
          Allocate (ffnsum(nffn,143),ffnenu(nffn,143))
          Do i=1,nffn
            Read(3,"(a1)") null
            Read(3,"(9(f8.3))") (ffnsum(i,j),ffnenu(i,j),j=1,143)
          Enddo
      Endif
!-----
!  Read in the reaction cross section data
      Open(2,file=trim(data_dir)//"/nets3",form='unformatted',
     &       status='old')
!-----
!  Allocate and read in reaction arrays for 1 reactant reactions
      nr=nreac(1)
      Allocate (csect1(nr),rc1(7,nr),q1(nr),rpf1(nr))
      Allocate (n1i(4,nr),iwk1(nr),ires1(nr),irev1(nr),iffn(nr))
      Do j=1,nr
        Read(2) n,(n1i(l,j),l=1,4),iwk1(j),ires1(j),irev1(j),
     &          (rc1(l,j),l=1,7),q1(j)
        If(n/=j) Then
            Write(6,*) 'Error in nets3, 1',j,n
            Exit
        Endif
      Enddo
!-----
!  Allocate and read in reaction arrays for 2 reactant reactions
      nr=nreac(2)
      Allocate (csect2(nr),rc2(7,nr),q2(nr),rpf2(nr),h2(nr))
      Allocate (n2i(5,nr),iwk2(nr),ires2(nr),irev2(nr))
      Do j=1,nr
        Read(2)  n,(n2i(l,j),l=1,5),iwk2(j),ires2(j),irev2(j),
     &         (rc2(l,j),l=1,7),q2(j)
        If(n/=j) Then
            Write(6,*) 'Error in nets3, 2',j,n
            Exit
        Endif
      Enddo
!-----
!  Allocate and read in reaction arrays for 3 reactant reactions
      nr=nreac(3)
      Allocate (csect3(nr),rc3(7,nr),q3(nr))
      Allocate (n3i(6,nr),iwk3(nr),ires3(nr),irev3(nr))
      Do j=1,nr
        Read(2)  n,(n3i(l,j),l=1,6),iwk3(j),ires3(j),irev3(j),
     &         (rc3(l,j),l=1,7),q3(j)
        If(n/=j) Then
            Write(6,*) 'Error in nets3, 3',j,n
            Exit
        Endif
      Enddo
!-----
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
!-----
!  Create and fill extended reaction->nuclei arrays
      nan(1)=le(1,ny)
      Allocate (tmu1(nan(1)),ta1(nan(1)),tb1(nan(1)),tn11(nan(1)))
      Allocate (tn1p1(nan(1)),tn1p2(nan(1)),tn1p3(nan(1)))
      Do j=1,nan(1)
        Read( 2)  ta1(j),tmu1(j)
        tn11(j)=n1i(1,tmu1(j))
        tn1p1(j)=n1i(2,tmu1(j))
        tn1p2(j)=n1i(3,tmu1(j))
        tn1p3(j)=n1i(4,tmu1(j))
      Enddo
      nan(2)=le(2,ny)
      Allocate (tmu2(nan(2)),ta2(nan(2)),tb2(nan(2)))
      Allocate (tn21(nan(2)),tn22(nan(2)))
      Allocate (tn2p1(nan(2)),tn2p2(nan(2)),tn2p3(nan(2)))
      Do j=1,nan(2)
        Read( 2)  ta2(j),tmu2(j)
        tn21(j)=n2i(1,tmu2(j))
        tn22(j)=n2i(2,tmu2(j))
        tn2p1(j)=n2i(3,tmu2(j))
        tn2p2(j)=n2i(4,tmu2(j))
        tn2p3(j)=n2i(5,tmu2(j))
      Enddo
      nan(3)=le(3,ny)
      Allocate (tmu3(nan(3)),ta3(nan(3)),tb3(nan(3)))
      Allocate (tn31(nan(3)),tn32(nan(3)),tn33(nan(3)))
      Allocate (tn3p1(nan(3)),tn3p2(nan(3)),tn3p3(nan(3)))
      Do j=1,nan(3)
        Read( 2)  ta3(j),tmu3(j)
        tn31(j)=n3i(1,tmu3(j))
        tn32(j)=n3i(2,tmu3(j))
        tn33(j)=n3i(3,tmu3(j))
        tn3p1(j)=n3i(4,tmu3(j))
        tn3p2(j)=n3i(5,tmu3(j))
        tn3p3(j)=n3i(6,tmu3(j))
      Enddo
      WRITE(*,*) 'QSE2'
      Return
      end ! qread_reaction_data

      subroutine qse_sort() 
!****************************************************************************
! The nubmer of groups and gp0i is determined in the subroutine
! group_number which must be called before this routine. 
! This will sort the reactions keeping only weak reactions and 
! those that cross into the independent nuclei. 
!----------------------------------------------------------------------------
      use nuclear_data
      use qse_data
      use reac_rate_data
      use ffn_data
      use cross_sect_data
      integer :: i,j,n,l,ki,ig,cn,cm,cn2,cn3
      integer :: i0,la1,le1,la2,le2,la3,le3,j1,l1,l2,l3
       character(len=5)  :: nnname
      character (LEN=1) :: null
      integer :: nr,idel,jdel,ijd,ndi,nds,plotindex(ny)
      common /ceign/ nds,ndi,idel,jdel,ijd
      j=0

! When the nubmer of groups changes in the middle of a run,
! the reactions must be resorted to suit the new groups
! Thus the massive deallocation below

      If (adaptive>0)then
        DeAllocate (singlenuc)
     	DeAllocate (laq,leq,laq1,leq1,laq2,leq2)
      	DeAllocate (laq3,leq3,part1,part2i,part2)
      	DeAllocate (part1b,part2ib,part2b)

     	DeAllocate (mu1,a1,b1,qnnum,qznum)
      	DeAllocate (n11,n1p1,n1p2,n1p3,regrp1)

     	DeAllocate (mu2,a2,b2,qnnum2,qznum2)
      	DeAllocate (n21,n22,n2p1,n2p2,n2p3,regrp2)

        DeAllocate (mu3,a3,b3,qnnum3,qznum3)
       	DeAllocate (n31,n32,n33,n3p1,n3p2,n3p3)
      	DeAllocate(regrp3)
          
       Endif
!Allocate all the sorting arrays
      Allocate (singlenuc(0:ny))
      Allocate (laq(3,gp0i),leq(3,gp0i))
      Allocate (laq1(1,gp0i),leq1(1,gp0i))
      Allocate (laq2(2,gp0i),leq2(2,gp0i))
      Allocate (laq3(2,gp0i),leq3(2,gp0i))
      Allocate (part1(gp0i,gp0i))
      Allocate (part2i(gn,gn))
      Allocate (part2(gn,gn))
      Allocate (part1b(gp0i,gp0i))
      Allocate (part2ib(gn,gn))
      Allocate (part2b(gn,gn))
          singlenuc=0
          Do k=1,gp0i
             singlenuc(gp0number(k))=k
          Enddo !k
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!3
        plotindex=1
        DO i=1,ny
          DO  j=1,gp1i
           If(i==gp1number(j)) then
              plotindex(i)=2
           Endif
         Enddo !j
        Do j=1,gp2i
           If(i==gp2number(j)) then
              plotindex(i)=3
           Endif
         Enddo !j
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!3
      If(numberofgroups.ge.3)then
        Do j=1,gp3i
           If(i==gp3number(j)) then
              plotindex(i)=4
            Endif
         Enddo !j
       Endif !numberofgroups
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!3
        Enddo !i
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C Now we grab all the reaction in group1 and establish a laq
C and leq for the whole group for 1 particle reactions.
Cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! We do this loop twice. The first time we count
! how many reactions are needed for
! the light group. Then we allocate that many.
! qmu1-contains the rxn  number
! qa1-contains the doubble counting corrections
! qnnum-contains the number of neutrons in the lightgroup nucleus
!  that is related to the reactions qmu1.
! qznum-contains the number of protons in the light group
! nucleus that is related to the reaction.
! The last two quantites are needed for eq. 16 in HPFT07.
! The second pass theough the loop assigns values to the varibales.
! There may be a better way to do this, but the routine is only called once
! at the begining so its speed is not that important.

! Trip 1 count the reactions of the light group.
! We keep reactions that cross the light group boundry:
! Since all these reaction have at least 1 light member, I search for reactions
! that also have singlenuc.
! We also keep weak reactions that occur within the lightgroup. Xnet already
! sets the flag iwk1 > 0 if it is a weak reaction, so I use this falg to sort.
! I have to start the loop at 2+ the number of groups becase the first several
! merbers of the singlenuc are the focal nuclei + protons.
      cn=0
      cn2=0
      cn3=0
      leq=0
      laq=0
          laq(1,2)=1
          laq(2,2)=1
          laq(3,2)=1
      Do i=1,gp1i
        ig=gp1number(i)
        la1=la(1,ig)
        le1=le(1,ig)
        Do k=la1,le1
         Do j=1+numberofgroups+1,gp0i
           If (tn11(k)==gp0number(j).or.tn1p1(k)==gp0number(j).or.
     &        tn1p2(k)==gp0number(j).or.tn1p3(k)==gp0number(j).and.
     &        iwk1(tmu1(k)).le.0)then
              cn=cn+1
!             If(ID(tn11(k))==2)then
!             rcn2=rcn2+1
!           Elseif(ID(tn11(k))==3)then
!               rcn3=rcn3+1
!            Elseif(ID(tn11(k))==4)then
!              rcn4=rcn4+1
!            Elseif(ID(tn11(k))==5)then
!              rcn5=rcn5+1
!            Elseif(ID(tn11(k))==0)then
!                rcn0=rcn0+1
!             Endif
!           Write(993,*) rcn2,rcn3,rcn4,rcn5,rcn0

           Exit !so that you don't count twice for rxns with 2 singlenuc.
           Endif
        Enddo ! j
       Enddo ! k
! now for 2p rxns
        la2=la(2,ig)
        le2=le(2,ig)
        Do k=la2,le2
         Do j=1+numberofgroups+1,gp0i
           If (tn21(k)==gp0number(j).or.tn22(k)==gp0number(j).or.
     &        tn2p1(k)==gp0number(j).or.
     &        tn2p2(k)==gp0number(j).or.tn2p3(k)==gp0number(j).and.
     &        iwk2(tmu2(k)).le.0)then
              cn2=cn2+1
           Exit !so that you don't count twice for rxns with 2 singlenuc.
           Endif
        Enddo ! j
       Enddo ! k
! Now for 3p reactions
        la3=la(3,ig)
        le3=le(3,ig)
        Do k=la3,le3
         Do j=1+numberofgroups+1,gp0i
           If (tn31(k)==gp0number(j).or.tn32(k)==gp0number(j).or.
     &        tn33(k)==gp0number(j).or.tn3p1(k)==gp0number(j).or.
     &        tn3p2(k)==gp0number(j).or.tn3p3(k)==gp0number(j).and.
     &        iwk3(tmu3(k)).le.0)then
              cn3=cn3+1
           Exit !so that you don't count twice for rxns with 2 singlenuc.
           Endif
        Enddo ! j
       Enddo ! k

      Enddo !i
       leq(1,2)=cn ! endpoint  of rxns that involve the lightgroup
       leq(2,2)=cn2! endpoint  of rxns that involve the lightgroup
       leq(3,2)=cn3 ! endpooint of rxns that involve the lightgroup
! Now count the 2nd group
          laq(1,3)=leq(1,2)+1
          laq(2,3)=leq(2,2)+1
          laq(3,3)=leq(3,2)+1
!          write(916,*) laq(1,2), leq(1,1)
      Do i=1,gp2i
        ig=gp2number(i)
        la1=la(1,ig)
        le1=le(1,ig)
        Do k=la1,le1
         Do j=1+numberofgroups+1,gp0i
           If (tn11(k)==gp0number(j).or.tn1p1(k)==gp0number(j).or.
     &        tn1p2(k)==gp0number(j).or.tn1p3(k)==gp0number(j).and.
     &        iwk1(tmu1(k)).le.0)then
              cn=cn+1
           Exit ! so that you don't count twice for rxns with 2 singlenuc.
           Endif
        Enddo ! j
       Enddo ! k
        la2=la(2,ig)
        le2=le(2,ig)
        Do k=la2,le2
         Do j=1+numberofgroups+1,gp0i
           If (tn21(k)==gp0number(j).or.tn22(k)==gp0number(j).or.
     &        tn2p1(k)==gp0number(j).or.
     &        tn2p2(k)==gp0number(j).or.tn2p3(k)==gp0number(j).and.
     &        iwk2(tmu2(k)).le.0)then
              cn2=cn2+1
           Exit ! so that you don't count twice for rxns with 2 singlenuc.
           Endif
        Enddo ! j
       Enddo ! k

! Now for 3p reactions
        la3=la(3,ig)
        le3=le(3,ig)
        Do k=la3,le3
         Do j=1+numberofgroups+1,gp0i
           If (tn31(k)==gp0number(j).or.tn32(k)==gp0number(j).or.
     &        tn33(k)==gp0number(j).or.tn3p1(k)==gp0number(j).or.
     &        tn3p2(k)==gp0number(j).or.tn3p3(k)==gp0number(j).and.
     &        iwk3(tmu3(k)).le.0)then
              cn3=cn3+1
           Exit !so that you don't count twice for rxns with 2 singlenuc.
           Endif
        Enddo ! j
       Enddo ! k

       leq(1,3)=cn !end point of number of rxns that involve the 2nd group
       leq(2,3)=cn2 !end point of number of rxns that involve the 2nd group
       leq(3,3)=cn3 !end point of number of rxns that involve the 2nd group

      Enddo !i
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!3
       If(numberofgroups.ge.3)then
! Now count the 3nd group
          laq(1,4)=leq(1,3)+1
          laq(2,4)=leq(2,3)+1
          laq(3,4)=leq(3,3)+1
      Do i=1,gp3i
        ig=gp3number(i)
        la1=la(1,ig)
        le1=le(1,ig)
        Do k=la1,le1
         Do j=1+numberofgroups+1,gp0i
           If (tn11(k)==gp0number(j).or.tn1p1(k)==gp0number(j).or.
     &        tn1p2(k)==gp0number(j).or.tn1p3(k)==gp0number(j).and.
     &        iwk1(tmu1(k)).le.0)then
              cn=cn+1
           Exit !so that you don't count twice for rxns with 2 singlenuc.
           Endif
        Enddo ! j
       Enddo ! k
  
        la2=la(2,ig)
        le2=le(2,ig)
        Do k=la2,le2
         Do j=1+numberofgroups+1,gp0i
           If (tn21(k)==gp0number(j).or.tn22(k)==gp0number(j).or.
     &        tn2p1(k)==gp0number(j).or.
     &        tn2p2(k)==gp0number(j).or.tn2p3(k)==gp0number(j).and.
     &        iwk2(tmu2(k)).le.0)then
              cn2=cn2+1
           Exit ! so that you don't count twice for rxns with 2 singlenuc.
           Endif
        Enddo ! j
       Enddo ! k
! Now for 3p reactions
        la3=la(3,ig)
        le3=le(3,ig)
        Do k=la3,le3
         Do j=1+numberofgroups+1,gp0i
           If (tn31(k)==gp0number(j).or.tn32(k)==gp0number(j).or.
     &        tn33(k)==gp0number(j).or.tn3p1(k)==gp0number(j).or.
     &        tn3p2(k)==gp0number(j).or.tn3p3(k)==gp0number(j).and.
     &        iwk3(tmu3(k)).le.0)then
              cn3=cn3+1
           Exit !so that you don't count twice for rxns with 2 singlenuc.
           Endif
        Enddo ! j
       Enddo ! k

       leq(1,4)=cn !end point of number of rxns that involve the 3rd group
       leq(2,4)=cn2 !end point of number of rxns that involve the 3rd group
       leq(3,4)=cn3!end point of number of rxns that involve the 3rd group


      Enddo ! i
      Endif !number of groups
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!3
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!3
       If(numberofgroups.ge.4)then
! Now count the 4nd group
          laq(1,5)=leq(1,4)+1
          laq(2,5)=leq(2,4)+1
          laq(3,5)=leq(3,4)+1
      Do i=1,gp4i
        ig=gp4number(i)
        la1=la(1,ig)
        le1=le(1,ig)
        Do k=la1,le1
         Do j=1+numberofgroups+1,gp0i
           If (tn11(k)==gp0number(j).or.tn1p1(k)==gp0number(j).or.
     &        tn1p2(k)==gp0number(j).or.tn1p3(k)==gp0number(j).and.
     &        iwk1(tmu1(k)).le.0)then
              cn=cn+1
           Exit !so that you don't count twice for rxns with 2 singlenuc.
           Endif
        Enddo ! j
       Enddo ! k

        la2=la(2,ig)
        le2=le(2,ig)
        Do k=la2,le2
         Do j=1+numberofgroups+1,gp0i
           If (tn21(k)==gp0number(j).or.tn22(k)==gp0number(j).or.
     &        tn2p1(k)==gp0number(j).or.
     &        tn2p2(k)==gp0number(j).or.tn2p3(k)==gp0number(j).and.
     &        iwk2(tmu2(k)).le.0)then
              cn2=cn2+1
           Exit ! so that you don't count twice for rxns with 2 singlenuc.
           Endif
        Enddo ! j
       Enddo ! k
! Now for 3p reactions
        la3=la(3,ig)
        le3=le(3,ig)
        Do k=la3,le3
         Do j=1+numberofgroups+1,gp0i
           If (tn31(k)==gp0number(j).or.tn32(k)==gp0number(j).or.
     &        tn33(k)==gp0number(j).or.tn3p1(k)==gp0number(j).or.
     &        tn3p2(k)==gp0number(j).or.tn3p3(k)==gp0number(j).and.
     &        iwk3(tmu3(k)).le.0)then
              cn3=cn3+1
           Exit !so that you don't count twice for rxns with 2 singlenuc.
           Endif
        Enddo ! j
       Enddo ! k

       leq(1,4)=cn !end point of number of rxns that involve the 3rd group
       leq(2,4)=cn2 !end point of number of rxns that involve the 3rd group
       leq(3,4)=cn3!end point of number of rxns that involve the 3rd group


!grab 1
      Enddo ! i
      Endif !number of groups

! now we count all the single nuc reactions
! now we count all the single nuc reactions
       Do i=gn+1,gp0i
          laq(1,i)=cn+1
          ig=gp0number(i)
        la1=la(1,ig)
        le1=le(1,ig)
        Do k=la1,le1
           If(iwk1(tmu1(k)).le.0)then
!        Do j=1+numberofgroups+1,gp0i
!           If (tn11(k)==gp0number(j))then
             cn=cn+1
!          Exit !so that you don't count twice for rxns with 2 singlenuc.
           Endif
!        Enddo ! j
        Enddo !k

         leq(1,i)=cn
         laq(2,i)=cn2+1 !!!!!!!!!!!!!!!!!!!!!!!!!
          ig=gp0number(i)
        la2=la(2,ig)
        le2=le(2,ig)
       Do k=la2,le2
           If(iwk2(tmu2(k)).le.0)then
             cn2=cn2+1
        Endif
        Enddo !k
         leq(2,i)=cn2 !!!!!!!!!!!!!!!!!!!!!!!!
         laq(3,i)=cn3+1 !!!!!!!!!!!!!!!!!!!!!!!!!
        la3=la(3,ig)
        le3=le(3,ig)
        Do k=la3,le3
           If(iwk3(tmu3(k)).le.0)then
             cn3=cn3+1
         Endif
        Enddo !k
         leq(3,i)=cn3 !!!!!!!!!!!!!!!!!!!!!!!!
       Enddo !i

!Allocate the varaibles
      Allocate (mu1(cn))
      Allocate (a1(cn),b1(cn))
      Allocate (qnnum(cn))
      Allocate (qznum(cn))
      Allocate (n11(cn),n1p1(cn),n1p2(cn),n1p3(cn))
      Allocate(regrp1(cn,3))

      Allocate (mu2(cn2))
      Allocate (a2(cn2),b2(cn2))
      Allocate (qnnum2(cn2))
      Allocate (qznum2(cn2))
      Allocate (n21(cn2),n22(cn2),n2p1(cn2),n2p2(cn2),n2p3(cn2))
      Allocate(regrp2(cn2,3))

      Allocate (mu3(cn3))
      Allocate (a3(cn3),b3(cn3))
      Allocate (qnnum3(cn3))
      Allocate (qznum3(cn3))
      Allocate (n31(cn3),n32(cn3),n33(cn3),n3p1(cn3),
     &         n3p2(cn3),n3p3(cn3))
      Allocate(regrp3(cn3,3))
      regrp1=0
      regrp2=0
      regrp3=0
!Trip 2 Assign and reorganize the variables.
       cn=0
       cn2=0
       cn3=0
      Do i=1,gp1i
        ig=gp1number(i)
        la1=la(1,ig)
        le1=le(1,ig)
        Do k=la1,le1

          Do j=1+numberofgroups+1,gp0i
             If (tn11(k)==gp0number(j).or.tn1p1(k)==gp0number(j).or.
     &          tn1p2(k)==gp0number(j).or.tn1p3(k)==gp0number(j).and.
     &          iwk1(tmu1(k)).le.0)then
                cn=cn+1
                mu1(cn)=tmu1(k)
                a1(cn)=ta1(k)
                b1(cn)=tb1(k)
                qnnum(cn)=int(nn(ig))
                qznum(cn)=int(zz(ig))
                n11(cn)=tn11(k)
                n1p1(cn)=tn1p1(k)
                n1p2(cn)=tn1p2(k)
                n1p3(cn)=tn1p3(k)
          exit
            Endif
          Enddo ! j
       Enddo ! k
! now for 2p rxns
        la2=la(2,ig)
        le2=le(2,ig)
        Do k=la2,le2
         Do j=1+numberofgroups+1,gp0i
           If (tn21(k)==gp0number(j).or.tn22(k)==gp0number(j).or.
     &        tn2p1(k)==gp0number(j).or.
     &        tn2p2(k)==gp0number(j).or.tn2p3(k)==gp0number(j).and.
     &        iwk2(tmu2(k)).le.0)then
              cn2=cn2+1
                mu2(cn2)=tmu2(k)
                a2(cn2)=ta2(k)
                b2(cn2)=tb2(k)
                qnnum2(cn2)=int(nn(ig))
                qznum2(cn2)=int(zz(ig))
                n21(cn2)=tn21(k)
                n22(cn2)=tn22(k)
                n2p1(cn2)=tn2p2(k)
                n2p2(cn2)=tn2p2(k)
                n2p3(cn2)=tn2p3(k)
           Exit !so that you don't count twice for rxns with 2 singlenuc.
           Endif
        Enddo ! j
       Enddo ! k
! Now for 3p reactions
        la3=la(3,ig)
        le3=le(3,ig)
        Do k=la3,le3
         Do j=1+numberofgroups+1,gp0i
           If (tn31(k)==gp0number(j).or.tn32(k)==gp0number(j).or.
     &        tn33(k)==gp0number(j).or.tn3p1(k)==gp0number(j).or.
     &        tn3p2(k)==gp0number(j).or.tn3p3(k)==gp0number(j).and.
     &        iwk3(tmu3(k)).le.0)then
              cn3=cn3+1
               mu3(cn3)=tmu3(k)
                a3(cn3)=ta3(k)
                b3(cn3)=tb3(k)
                qnnum3(cn3)=int(nn(ig))
                qznum3(cn3)=int(zz(ig))
                n31(cn3)=tn31(k)
                n32(cn3)=tn32(k)
                n33(cn3)=tn33(k)
                n3p1(cn3)=tn3p2(k)
                n3p2(cn3)=tn3p2(k)
                n3p3(cn3)=tn3p3(k)

           Exit !so that you don't count twice for rxns with 2 singlenuc.
           Endif
        Enddo ! j
       Enddo ! k

      Enddo !i
!2nd Group
      Do i=1,gp2i
        ig=gp2number(i)
        la1=la(1,ig)
        le1=le(1,ig)
        Do k=la1,le1
          Do j=1+numberofgroups+1,gp0i
             If (tn11(k)==gp0number(j).or.tn1p1(k)==gp0number(j).or.
     &          tn1p2(k)==gp0number(j).or.tn1p3(k)==gp0number(j).and.
     &          iwk1(tmu1(k)).le.0)then
                cn=cn+1
                mu1(cn)=tmu1(k)
                a1(cn)=ta1(k)
                b1(cn)=tb1(k)
                qnnum(cn)=int(nn(ig)) !nn-14 to simplify!
                qznum(cn)=int(zz(ig))
                n11(cn)=tn11(k)
                n1p1(cn)=tn1p1(k)
                n1p2(cn)=tn1p2(k)
                n1p3(cn)=tn1p3(k)
          exit
            Endif

          Enddo ! j
       Enddo ! k
! now for 2p rxns
        la2=la(2,ig)
        le2=le(2,ig)
        Do k=la2,le2
         Do j=1+numberofgroups+1,gp0i
           If (tn21(k)==gp0number(j).or.tn22(k)==gp0number(j).or.
     &        tn2p1(k)==gp0number(j).or.
     &        tn2p2(k)==gp0number(j).or.tn2p3(k)==gp0number(j).and.
     &        iwk2(tmu2(k)).le.0)then
              cn2=cn2+1
                mu2(cn2)=tmu2(k)
                a2(cn2)=ta2(k)
                b2(cn2)=tb2(k)
                qnnum2(cn2)=int(nn(ig))
                qznum2(cn2)=int(zz(ig))
                n21(cn2)=tn21(k)
                n22(cn2)=tn22(k)
                n2p1(cn2)=tn2p2(k)
                n2p2(cn2)=tn2p2(k)
                n2p3(cn2)=tn2p3(k)
           Exit !so that you don't count twice for rxns with 2 singlenuc.
           Endif
        Enddo ! j
       Enddo !k
! Now for 3p reactions
        la3=la(3,ig)
        le3=le(3,ig)
        Do k=la3,le3
         Do j=1+numberofgroups+1,gp0i
           If (tn31(k)==gp0number(j).or.tn32(k)==gp0number(j).or.
     &        tn33(k)==gp0number(j).or.tn3p1(k)==gp0number(j).or.
     &        tn3p2(k)==gp0number(j).or.tn3p3(k)==gp0number(j).and.
     &        iwk2(tmu3(k)).le.0)then
              cn3=cn3+1
                mu3(cn3)=tmu3(k)
                a3(cn3)=ta3(k)
                b3(cn3)=tb3(k)
                qnnum3(cn3)=int(nn(ig))
                qznum3(cn3)=int(zz(ig))
                n31(cn3)=tn31(k)
                n32(cn3)=tn32(k)
                n33(cn3)=tn33(k)
                n3p1(cn3)=tn3p2(k)
                n3p2(cn3)=tn3p2(k)
                n3p3(cn3)=tn3p3(k)

           Exit !so that you don't count twice for rxns with 2 singlenuc.
           Endif
        Enddo ! j
       Enddo ! k
      Enddo !i
! Third group !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!3
      If(numberofgroups.ge.3)then
      Do i=1,gp3i
        ig=gp3number(i)
        la1=la(1,ig)
        le1=le(1,ig)
        Do k=la1,le1
          Do j=1+numberofgroups+1,gp0i
             If (tn11(k)==gp0number(j).or.tn1p1(k)==gp0number(j).or.
     &          tn1p2(k)==gp0number(j).or.tn1p3(k)==gp0number(j).and.
     &          iwk1(tmu1(k)).le.0)then
                cn=cn+1
                mu1(cn)=tmu1(k)
                a1(cn)=ta1(k)
                qnnum(cn)=int(nn(ig))
                b1(cn)=tb1(k)
                qznum(cn)=int(zz(ig))
                n11(cn)=tn11(k)
                n1p1(cn)=tn1p1(k)
                n1p2(cn)=tn1p2(k)
                n1p3(cn)=tn1p3(k)
          exit
            Endif
          Enddo ! j
       Enddo ! k
! now for 2p rxns
        la2=la(2,ig)
        le2=le(2,ig)
        Do k=la2,le2
         Do j=1+numberofgroups+1,gp0i
           If (tn21(k)==gp0number(j).or.tn22(k)==gp0number(j).or.
     &        tn2p1(k)==gp0number(j).or.
     &        tn2p2(k)==gp0number(j).or.tn2p3(k)==gp0number(j).and.
     &        iwk2(tmu2(k)).le.0)then
              cn2=cn2+1
                mu2(cn2)=tmu2(k)
                a2(cn2)=ta2(k)
                b2(cn2)=tb2(k)
                qnnum2(cn2)=int(nn(ig))
                qznum2(cn2)=int(zz(ig))
                n21(cn2)=tn21(k)
                n22(cn2)=tn22(k)
                n2p1(cn2)=tn2p2(k)
                n2p2(cn2)=tn2p2(k)
                n2p3(cn2)=tn2p3(k)
           Exit !so that you don't count twice for rxns with 2 singlenuc.
           Endif
        Enddo ! j
        Enddo !k
! Now for 3p reactions
        la3=la(3,ig)
        le3=le(3,ig)
        Do k=la3,le3
         Do j=1+numberofgroups+1,gp0i
           If (tn31(k)==gp0number(j).or.tn32(k)==gp0number(j).or.
     &        tn33(k)==gp0number(j).or.tn3p1(k)==gp0number(j).or.
     &        tn3p2(k)==gp0number(j).or.tn3p3(k)==gp0number(j).and.
     &        iwk3(tmu3(k)).le.0)then
              cn3=cn3+1
                mu3(cn3)=tmu3(k)
                a3(cn3)=ta3(k)
                b3(cn3)=tb3(k)
                qnnum3(cn3)=int(nn(ig))
                qznum3(cn3)=int(zz(ig))
                n31(cn3)=tn31(k)
                n32(cn3)=tn32(k)
                n33(cn3)=tn33(k)
                n3p1(cn3)=tn3p2(k)
                n3p2(cn3)=tn3p2(k)
                n3p3(cn3)=tn3p3(k)

           Exit !so that you don't count twice for rxns with 2 singlenuc.
           Endif
!grab2 299
        Enddo ! j
       Enddo ! k
      Enddo !i
       Endif ! number of groups!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!3
! 4th group !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!3
      If(numberofgroups.ge.4)then
      Do i=1,gp4i
        ig=gp4number(i)
        la1=la(1,ig)
        le1=le(1,ig)
        Do k=la1,le1
          Do j=1+numberofgroups+1,gp0i
             If (tn11(k)==gp0number(j).or.tn1p1(k)==gp0number(j).or.
     &          tn1p2(k)==gp0number(j).or.tn1p3(k)==gp0number(j).and.
     &          iwk1(tmu1(k)).le.0)then
                cn=cn+1
                mu1(cn)=tmu1(k)
                a1(cn)=ta1(k)
                qnnum(cn)=int(nn(ig))
                b1(cn)=tb1(k)
                qznum(cn)=int(zz(ig))
                n11(cn)=tn11(k)
                n1p1(cn)=tn1p1(k)
                n1p2(cn)=tn1p2(k)
                n1p3(cn)=tn1p3(k)
          exit
            Endif
          Enddo ! j
       Enddo ! k
! now for 2p rxns
        la2=la(2,ig)
        le2=le(2,ig)
        Do k=la2,le2
         Do j=1+numberofgroups+1,gp0i
           If (tn21(k)==gp0number(j).or.tn22(k)==gp0number(j).or.
     &        tn2p1(k)==gp0number(j).or.
     &        tn2p2(k)==gp0number(j).or.tn2p3(k)==gp0number(j).and.
     &        iwk2(tmu2(k)).le.0)then
              cn2=cn2+1
                mu2(cn2)=tmu2(k)
                a2(cn2)=ta2(k)
                b2(cn2)=tb2(k)
                qnnum2(cn2)=int(nn(ig))
                qznum2(cn2)=int(zz(ig))
                n21(cn2)=tn21(k)
                n22(cn2)=tn22(k)
                n2p1(cn2)=tn2p2(k)
                n2p2(cn2)=tn2p2(k)
                n2p3(cn2)=tn2p3(k)
           Exit !so that you don't count twice for rxns with 2 singlenuc.
           Endif
        Enddo ! j
        Enddo !k
! Now for 3p reactions
        la3=la(3,ig)
        le3=le(3,ig)
        Do k=la3,le3
         Do j=1+numberofgroups+1,gp0i
           If (tn31(k)==gp0number(j).or.tn32(k)==gp0number(j).or.
     &        tn33(k)==gp0number(j).or.tn3p1(k)==gp0number(j).or.
     &        tn3p2(k)==gp0number(j).or.tn3p3(k)==gp0number(j).and.
     &        iwk3(tmu3(k)).le.0)then
              cn3=cn3+1
                mu3(cn3)=tmu3(k)
                a3(cn3)=ta3(k)
                b3(cn3)=tb3(k)
                qnnum3(cn3)=int(nn(ig))
                qznum3(cn3)=int(zz(ig))
                n31(cn3)=tn31(k)
                n32(cn3)=tn32(k)
                n33(cn3)=tn33(k)
                n3p1(cn3)=tn3p2(k)
                n3p2(cn3)=tn3p2(k)
                n3p3(cn3)=tn3p3(k)

           Exit !so that you don't count twice for rxns with 2 singlenuc.
           Endif
        Enddo ! j
       Enddo ! k
      Enddo !i
       Endif ! number of groups!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!3

! now we reorder all the single nuc reactions
       Do i=gn+1,gp0i
          ig=gp0number(i)
        la1=la(1,ig)
        le1=le(1,ig)
        Do k=la1,le1
           If(iwk1(tmu1(k)).le.0)then
             cn=cn+1
                mu1(cn)=tmu1(k)
                a1(cn)=ta1(k)
                b1(cn)=tb1(k)
                qnnum(cn)=int(nn(ig))
                qznum(cn)=int(zz(ig))
                n11(cn)=tn11(k)
                n1p1(cn)=tn1p1(k)
                n1p2(cn)=tn1p2(k)
                n1p3(cn)=tn1p3(k)
            Endif
        Enddo !j
! now for 2p rxns
        la2=la(2,ig)
        le2=le(2,ig)
        Do k=la2,le2
           If(iwk2(tmu2(k)).le.0)then
              cn2=cn2+1
                mu2(cn2)=tmu2(k)
                a2(cn2)=ta2(k)
                b2(cn2)=tb2(k)
                qnnum2(cn2)=int(nn(ig))
                qznum2(cn2)=int(zz(ig))
                n21(cn2)=tn21(k)
                n22(cn2)=tn22(k)
                n2p1(cn2)=tn2p2(k)
                n2p2(cn2)=tn2p2(k)
                n2p3(cn2)=tn2p3(k)
             Endif
        Enddo !k
!And 3p rxns
        la3=la(3,ig)
        le3=le(3,ig)
        Do k=la3,le3
           If(iwk3(tmu3(k)).le.0)then
             cn3=cn3+1
                mu3(cn3)=tmu3(k)
                a3(cn3)=ta3(k)
                b3(cn3)=tb3(k)
                qnnum3(cn3)=int(nn(ig))
                qznum3(cn3)=int(zz(ig))
                n31(cn3)=tn31(k)
                n32(cn3)=tn32(k)
                n33(cn3)=tn33(k)
                n3p1(cn3)=tn3p2(k)
                n3p2(cn3)=tn3p2(k)
                n3p3(cn3)=tn3p3(k)
            endif
        Enddo !k

      Enddo !i


      Return
      End
!-------------------------------------------------------------------------
!*************************************************************************
!*****************************************************************************
!*****************************************************************************
!*****************************************************************************
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine qse_coeff
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! Makes the qse coefficients for all the group members.

! I have run 1-299 instead of just handling the group members
! in order to aviod indirect addressing. I'm not sure if it
! is the best idea.
! I have made a second array CQSEb. This will eventually be used for
! all the coeff but right now it is just used for the Fe group.
! I have done this because dexp(dlCqse(ig)-dlCqse(186)) becomes too large
! for T9<3.1, i.e. it's larger than 711; Exp(x>711) = NaN. Thus changes
! were also made in the Jacobian so that Cqseb was multiplied
! by something small, such as log(y(p))*(Zi-ZFe56), before
! putting it in the exponent.
! I will make this change for the other groups when I have figured out a
! better way to build the Jacobian.
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      use nuclear_data
      use qse_data
      use controls
      use conditions
      use abundances
      use match_data
      use flux_data
      use part_funct_data
      real(8) :: xg(7),enuc,edot,tmp,bkt,dummy
      integer :: i,j,k,ig,kstep,kout,zz1(ny),aa1(ny),nn1(ny),is
      parameter (pi=3.14159,hbr=6.58217D-22,amu=1.036435E-18)
      parameter (bok=8.6174D-02,avn=6.02205D+23)
      real(8) ::  tmpa,Yqse(ny),Yqse2(ny)
     & ,Fesum,t9oldy,c(ny)
!         call t9rhofind(kstep)
         call partf(t9t)
         bkt=bok*t9t
         tmp=dlog(((avn*rhot)**(2./3.)*2.*pi*hbr**2)/(bkt*amu))
!cc      write(928,*) 'Const',rhot,avn,hbr,pi,bkt,tmp
         dlCqse(1)=0.
         dlCqse(2)=0.
         DO 220 l=3,ny
            tmpa=dlog(angm(l)*gg(l)*(aa(l)**1.5)*(0.5**aa(l)))
            dlCqse(l)=tmpa+(1.5*(aa(l)-1))*tmp+(be(l)/bkt)

220    CONTINUE

! Get the ratios of the qse coeff for each of the groups.
!      Cqse1=0.0
!      Cqse3=0.0
!      Cqse4=0.0
!      Do i=1,gp1i
!           Cqse1(gp1number(i))=dlCqse(gp1number(i))
!      Enddo
!      Do i=1,gp2i
!           ig=gp2number(i)
!           Cqse3(ig)=(dlCqse(ig)-dlCqse(f3))
!      Enddo
!      If(numberofgroups>=3)then
!        Do i=1,gp3i
!              ig=gp3number(i)
!               Cqse4(ig)=(dlCqse(ig)-dlCqse(f4))
!        Enddo
!      Endif

        Return
        End

      subroutine qse_make(yf)
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c This routine used the focal abundances of from the
!c full set to build the groups.
!c 1. reads in 299 full set
!c 2. uses 4 focal from full set to build the groups.
!c 3. Substitutes the groups back in to the full set leaving
!c  the single nuc therein unchanged.
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      use nuclear_data
      use qse_data
      use controls
      use conditions
      use abundances
      use match_data
      use flux_data
      use part_funct_data
      real(8) :: xg(7),enuc,edot,tmp,bkt,dummy
      integer :: i,j,k,ig,n1,n2,n3,z1,z2,z3
      Real(8) :: tmpa,yf(ny)
     & ,Fesum,t9oldy
!         n3=nn(f3)
!         z3=zz(f3)
         y1=yf(f1)
         y2=yf(f2)
!         y3=yf(f3)


         Do i=1,gp1i
            ig=gp1number(i)
            yf(ig)=dexp(dlCqse(ig)+dlog(y2)*(zz(ig))
     &             +dlog(y1)*(nn(ig)))
         Enddo
         If(numberofgroups>=2)then
           n3=nn(f3)
           z3=zz(f3)
           y3=yf(f3)
           Do i=1,gp2i
             ig=gp2number(i)
             yf(ig)=y3*dexp(dlCqse(ig)-dlCqse(f3)+dlog(y2)*(zz(ig)-z3)
     &              +dlog(y1)*(nn(ig)-n3))
           Enddo
         Endif

         If(numberofgroups>=3)then
           n4=nn(f4)
           z4=zz(f4)
           y4=yf(f4)
            Do i=1,gp3i
              ig=gp3number(i)
              yf(ig)=y4*dexp(dlCqse(ig)-dlCqse(f4)+dlog(y2)*(zz(ig)-z4)
     &               +dlog(y1)*(nn(ig)-n4))
            Enddo
          Endif

        Return
        End

      subroutine createf(yfin,x28,x56,x184,ye,ff,fmax1)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c This makes the vector f that will be used to get yr back out of yg.
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      use nuclear_data
      use qse_data
      use controls
      use conditions
      use abundances
      use match_data
      use flux_data
      use part_funct_data



      real(8) yfin,x28,x56,x184,ff,fmax1,ye
      integer i,indxx,j,ig

      dimension yfin(ny),ff(gn)

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c Generation of Vector ff                                                      c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      ff(1)=-1.
      ff(2)=-ye
      fmax1=0.
      Do i=1,ny
        ff(1)=ff(1)+aa(i)*yfin(i)
        ff(2)=ff(2)+zz(i)*yfin(i)
      Enddo
! Si group
      If(numberofgroups>=2)then
        ff(3)=-x28
        Do ig =1,gp2i
          i=gp2number(ig)
          ff(3)=ff(3)+aa(i)*yfin(i)
        Enddo
      Endif
      If(numberofgroups>=3)then
        ff(4)=-x56
! Fe group
        Do ig =1,gp3i
          i=gp3number(ig)
          ff(4)=ff(4)+aa(i)*yfin(i)
        Enddo
      Endif
      Do i=1,gn
          if (dabs(ff(i)).gt.fmax1) fmax1=dabs(ff(i))
      Enddo
      Return
      End
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine update(ys)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c In: The initial guess at the qse yf abundaces, this is called ys
c as in y used in the subroutine.
c Out: The proper starting qse abundances for yf.
c The goal is to use the derivative to get a normalized set of
c all abudnaces from  the 4 focal abundance in the full set.
c the single nuc remain untouched but the groups nuc get redone.
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      use nuclear_data
      use qse_data
      use controls
      use conditions
      use abundances
      use match_data
      use flux_data
      use part_funct_data
        integer :: i, iter_count1
        Real(8) :: yf(ny),y1,y2,y1n,y2n,ys(ny)
        Real(8) :: df(gn,gn),fin(gn),invdf(gn,gn),dummy(gn)
        Real(8) :: ye,x28,x56,x184
        integer :: n77,z77,n186,z186,nsi,zsi,nfe,zfe,indxx(gn),info
! Get the initial yf and then use it with the QSE coeffients to rebuild the
! all the group nuclei. Put the new nuc in Y.
        x28=0.0
        ye=0.0
        x56=0.0
        x184=0.0
        ye =sum(zz*y)
      If(numberofgroups>=2)then
        Do ig=1,gp2i
           i=gp2number(ig)
           x28=x28+y(i)*aa(i)
        Enddo
      Endif
      If(numberofgroups>=3)then
        Do ig=1,gp3i
          i=gp3number(ig)
          x56=x56+y(i)*aa(i)
        Enddo
      Endif
         call qse_make(ys) ! this is the first time that you solve for qse.
!         call norm(ys) !  Normalizes the abundances
         iter_count1=0
c This is the start of the loop that iterates until df converges

20       continue
ccccccccccccccccccccccccccccccccccccccccccccccc
        iter_count1=iter_count1+1
        If(iter_count1>20) then
      Write(*,*) "Initial abundances do not converge."
           stop
        EndIf

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c Make f f(1)=-1+ff(1),-ye+ff(2),-x28+ff(3),-x56+ff(4)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
          call createf(ys,x28,x56,x184,ye,fin,fmax)
            If (fmax.lt.1d-13) goto 21
            If (fmax.gt.1D12) stop 'keine Konvergenz'

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c make df
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
          ryn=1/ys(1)
          ryp=1/ys(2)
          y1=ys(1)
          y2=ys(2)
          df=0.0
          Do i=1,gp1i
            df(1,1) = df(1,1)+aa(i)*nn(i)*ys(i)*ryn
            df(1,2) = df(1,2)+aa(i)*zz(i)*ys(i)*ryp
            df(2,1) = df(2,1)+zz(i)*nn(i)*ys(i)*ryn
            df(2,2) = df(2,2)+zz(i)*zz(i)*ys(i)*ryp
          Enddo
        If(numberofgroups>=2)then
          Do ig =1,gp2i
            i=gp2number(ig)
            nsi=nn(i)-nn(f3)
            zsi=zz(i)-zz(f3)
             Cpart=dexp(dlCqse(i)-dlCqse(f3)+dlog(y1)*nsi+
     &             dlog(y2)*zsi)
! build Jacobiean
            df(1,1)=df(1,1)+aa(i)*nsi*ys(i)*ryn
            df(1,2)=df(1,2)+aa(i)*zsi*ys(i)*ryp
            df(1,3)=df(1,3)+aa(i)*Cpart
            df(2,1)=df(2,1)+zz(i)*nsi*ys(i)*ryn
            df(2,2)=df(2,2)+zz(i)*zsi*ys(i)*ryp
            df(2,3)=df(2,3)+zz(i)*Cpart
            df(3,1)=df(3,1)+aa(i)*nsi*ys(i)*ryn
            df(3,2)=df(3,2)+aa(i)*zsi*ys(i)*ryp
            df(3,3)=df(3,3)+aa(i)*Cpart
          Enddo
        Endif
        If(numberofgroups>=3)then
          Do ig =1,gp3i
            i=gp3number(ig)
            nfe=nn(i)-nn(f4)
            zfe=zz(i)-zz(f4)
            Cpart=dexp(dlCqse(i)-dlCqse(f4)+dlog(y1)*nfe+dlog(y2)*zfe)
            df(1,1)=df(1,1)+aa(i)*nfe*ys(i)*ryn
            df(1,2)=df(1,2)+aa(i)*zfe*ys(i)*ryp
            df(1,4)=df(1,4)+aa(i)*Cpart
            df(2,1)=df(2,1)+zz(i)*nfe*ys(i)*ryn
            df(2,2)=df(2,2)+zz(i)*zfe*ys(i)*ryp
            df(2,4)=df(2,4)+zz(i)*Cpart
            df(4,1)=df(4,1)+aa(i)*nfe*ys(i)*ryn
            df(4,2)=df(4,2)+aa(i)*zfe*ys(i)*ryp
            df(4,4)=df(4,4)+aa(i)*Cpart
          Enddo
         Endif
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c find inverse of df
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
           invdf=0.0
           do i=1,gn
              invdf(i,i)=1
           enddo
! If you are without LAPACK use these.
!            call ludcmp(df,gn,gn,indxx,info)
!            call lubksb(df,gn,gn,indxx,invdf)
! LAPACK LU decomposition.
           call dgesv(gn,gn,df,gn,indxx,invdf,gn,info)
           dummy=0
! multiply it times f.
           Do i=1,gn
             Do j=1,gn
               dummy(i)=dummy(i)+invdf(i,j)*fin(j)
             Enddo
           Enddo

!           Do i=1,gn
!            ys(gp0number(i))=ys(gp0number(i))-dummy(i)
!            If (ys(i).lt.0) ys(i)=-ys(i)
!           Enddo
           Do i=1,gn
              j=gp0number(i)
            ys(j)=ys(j)-dummy(i)
            If (ys(j).lt.0) ys(j)=0.1*(ys(j)+dummy(i))
           Enddo
       call qse_make(ys)
       goto 20

21        Return
       End
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
       subroutine update2(yg,yf)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c turns yg in to yf. This involves a Jacobian because yg(1-4) c
c are group abudances, not n, p, si and fe.                    c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      use nuclear_data
      use qse_data
      use controls
      use conditions
      use abundances
      use match_data
      use flux_data
      use part_funct_data

       real(8) :: yg(gp0i),ff(gn),df(gn,gn),delta(gn),d
       real(8) :: y1,y2,y1n,y2n,fmax,yf(ny),dummy12
       real(8) :: sumng,sumpg,sumsig,sumfeg,dummy(gn,gn),sumcrg
       integer :: nsi,zsi,nfe,zfe,ncr,zcr,indxx(gn)
       integer :: i,counter,info
        counter=0
20       call qse_make(yf)
         counter=counter+1
         df=0.0
         fmax=0
         sumsig=0
         sumfeg=0
         sumng=0
         sumpg=0
         sumcrg=0
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c Calculate the vector f.
c The first time through sumng=Yg(1) sumpg=Yg(2) etc. so f=0
c After that the sums use the adjusted yfs so f n.e. 0
c This porcess cranks though the Jacobian solve with each new yf
c until  fmax=dabs(ff(i)/yg(i)) is less than 1e-8.
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
         Do i=1, gp1i
            sumng=sumng+nn(i)*yf(i)
            sumpg=sumpg+zz(i)*yf(i)
          Enddo
      If(numberofgroups>=2)then
        Do ig=1,gp2i
          i=gp2number(ig)
          sumng=sumng+(nn(i)-nn(f3))*yf(i)
          sumpg=sumpg+(zz(i)-zz(f3))*yf(i)
          sumsig=sumsig+yf(i)
        Enddo
         ff(3)=yg(3)-sumsig
      Endif
      If(numberofgroups>=3)then
        Do ig=1,gp3i
          i=gp3number(ig)
          sumng=sumng+(nn(i)-nn(f4))*yf(i)
          sumpg=sumpg+(zz(i)-zz(f4))*yf(i)
          sumfeg=sumfeg+yf(i)
        Enddo
       ff(4)=yg(4)-sumfeg
      Endif

      ff(1)=yg(1)-sumng
      ff(2)=yg(2)-sumpg
          Do i=1,gn
            If (yg(i).ne.0.) then
              If (dabs(ff(i)/yg(i)).gt.fmax) then
                 fmax=dabs(ff(i)/yg(i))
              Endif
            Endif
          Enddo
      If (counter.gt.30) then
         write(*,*) fmax, 'Mehr als 20 Iterationen'
         goto 21
      Endif
      If (fmax.lt.1d-8) goto 21

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c Calculate Jacobian df                                                       c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
          ryn=1/yf(1)
          ryp=1/yf(2)
          y1=yf(1)
          y2=yf(2)
          df=0.0
       Do i=1,gp1i
          df(1,1)=df(1,1)+nn(i)*nn(i)*yf(i)/yf(1)
          df(1,2)=df(1,2)+nn(i)*zz(i)*yf(i)/yf(2)
          df(2,1)=df(2,1)+nn(i)*zz(i)*yf(i)/yf(1)
          df(2,2)=df(2,2)+zz(i)*zz(i)*yf(i)/yf(2)
       Enddo
       If(numberofgroups>=2)then
        Do ig =1,gp2i
           i=gp2number(ig)
           nsi=nn(i)-nn(f3)
           zsi=zz(i)-zz(f3)
           Cpart=dexp(dlCqse(i)-dlCqse(f3)+dlog(y1)*nsi+ dlog(y2)*zsi)
! build Jacobiean
           df(1,1)=df(1,1)+nsi*nsi*yf(i)*ryn
           df(1,2)=df(1,2)+nsi*zsi*yf(i)*ryp
           df(1,3)=df(1,3)+nsi*Cpart
           df(2,1)=df(2,1)+nsi*zsi*yf(i)*ryn
           df(2,2)=df(2,2)+zsi*zsi*yf(i)*ryp
           df(2,3)=df(2,3)+zsi*Cpart
           df(3,1)=df(3,1)+nsi*yf(i)*ryn
           df(3,2)=df(3,2)+zsi*yf(i)*ryp
           df(3,3)=df(3,3)+Cpart
         Enddo
        Endif
        If(numberofgroups>=3)then
         Do ig =1,gp3i
           i=gp3number(ig)
           nfe=nn(i)-nn(f4)
           zfe=zz(i)-zz(f4)
           Cpart=dexp(dlCqse(i)-dlCqse(f4)+dlog(y1)*nfe+dlog(y2)*zfe)
           df(1,1)=df(1,1)+nfe*nfe*yf(i)*ryn
           df(1,2)=df(1,2)+nfe*zfe*yf(i)*ryp
           df(1,4)=df(1,4)+nfe*Cpart
           df(2,1)=df(2,1)+zfe*nfe*yf(i)*ryn
           df(2,2)=df(2,2)+zfe*zfe*yf(i)*ryp
           df(2,4)=df(2,4)+zfe*Cpart
           df(4,1)=df(4,1)+nfe*yf(i)*ryn
           df(4,2)=df(4,2)+zfe*yf(i)*ryp
           df(4,4)=df(4,4)+Cpart
          Enddo
         Endif !numberofgroups
          Do i=1,gn
             delta(i)=ff(i)
          Enddo
!In case you don't have LAPACK uses these two.
!       call ludcmp(df,gn,gn,indxx,d)
!       call lubksb(df,gn,gn,indxx,delta)
! LAPACK LU decompostion
          call dgesv(gn,1,df,gn,indxx,delta,gn,info)
            Do ig=1,gn
                i=gp0number(ig)
                yf(i) = yf(i)+delta(ig)
                if(yf(i)<0)then
!                yf(i)=yf(i)-delta(ig)+0.1*delta(ig)
                 goto 21
                  write(*,*) i, "less than 0", Cpart
                endif
            Enddo
        goto 20
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c update the other abundances given in yg(5) until yg(ni)                     c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

21      continue
           Do ig=gn+1,gp0i
              i=gp0number(ig)

              yf(i)=yg(ig)
           Enddo

       Return
        End
!******************************************************************************


      subroutine qse_net(izone,data_dir)
!===============================================================================
!  The abundances are evolved using a Newton-Raphson iteration scheme to solve
!  the equation yt(i)-y(i)/tdel = ydot(i), where y(i) is the abundance at the
!  beginning of the iteration, yt(i) is the trial abundance at the end of the
!  timestep, ydot(i) is the time derivative of the trial abundance  calculated
!  from reaction rates and tdel is the timestep.
!===============================================================================
      use controls
      use qse_data
      use nuclear_data
      use conditions
      use abundances
      use thermo_data
!      real(8) :: reldy(ny)
      real(8) :: toln      !  Network convergence condition
      real(8) :: t9tplus,t9tminus
      real(8) :: testc,testc2,testm,testn    ! Convergence tests
      real(8) :: ta,tdela,xtot,xtoto,enm,enb,enold,en0,edot,ye
      real(8) :: ytot,ztot,atot,Ygroup(ny),ydummy(ny),xsum
      integer irdymx(1),idymx(1),is,j
      integer :: izone,k,kstep,ktsmx,kout,idiag0,countiter
      integer :: nsi,zsi,nfe,zfe,ncr,zcr
      character (LEN=80) :: data_dir
      INTEGER :: clock_start,clock_end,clock_rate
      double precision  elapsed_time,accum2


! Set up QSE
      Allocate(dlCqse(ny)) ! here is where we set the size of the QSE_coeff arra
      Allocate(drvN(ny),drvZ(ny),ID(ny))
         numberofgroups =3
         call qread_reaction_data(data_dir) 
         call read_match_data(data_dir)
         call flux_init
	 t9plus=t9t+1
	 t9minus=t9t-1
         adaptive=0      
        call group_number(data_dir,t9t)




!  Set reaction controls not read in from control
      idiag0=idiag
      ktsmx=30
      kstep=0
      kout=0
      countiter=0
!  Normalize initial abundances and change units if necessary
      yt=y
      call norm(yt)
!  Calculate the total energy of the nuclei
      call benuc(yt,enb,enm,ytot,ztot,atot,t)
      xtot=atot-1.0
      en0=enm
      edot=0.0

!  Start evolution
      Write(6,*) 'Max Step',kstmx,'IDiag=',idiag
      t=tstart
      tt=tstart
      call t9rhofind(kstep,tt,j)
! Build QSE coeff. from 4 focal
      call qse_coeff !
!repalce initial abund. with group members.
        call qse_make(yt)
!Solve for QSE this is step 0.
        call update(yt) ! y is getting screwed up
        ye=sum(zz*yt)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cBuild yg. This is the set with the qse group abundances. Note: it is
c the same size as yr and has the same single nuc except for the focal
c abundances.
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
           yq =0.0

           do i=1,gp1i
              yq(1) =yq(1) + nn(i)*yt(i)
              yq(2) =yq(2)+  zz(i)*yt(i)
            enddo
            do ig=1,gp2i
               i=gp2number(ig)
               nsi=nn(i)-nn(f3)
               zsi=zz(i)-zz(f3)
               yq(1)=yq(1)+nsi*yt(i)
               yq(2)=yq(2)+zsi*yt(i)
               yq(3)=yq(3)+yt(i)
             enddo
             do ig=1,gp3i
               i=gp3number(ig)
               nfe=nn(i)-nn(f4)
               zfe=zz(i)-zz(f4)
               yq(1)=yq(1)+nfe*yt(i)
               yq(2)=yq(2)+zfe*yt(i)
               yq(4)=yq(4)+yt(i)
             enddo
             do ig=gn+1,gp0i
                i=gp0number(ig)
                yq(ig)=yt(i)
             enddo

                yqt=yq
!           
!-----------------------------------------------------------------------------
!  For each step, an initial guess for the timestep is calculated, based on
!  the abundance variations of the previous step if available.
!-----------------------------------------------------------------------------
      Step: Do kstep=1,kstmx
          call timestep(kstep) ! called t9rhofind,which changes rho and T

          call qse_coeff ! must be called every time rho or T  changes
! Updates groups and puts
          call qse_make(yt)
        If(idiag>=1) Write(50,*) 'TDel',tt,tdel
        xtoto=xtot

!  Determine if this is an output step
        idiag=idiag0
!       If(mod(kstep,10).eq.0) idiag=2
!       If(kstep==518.or.kstep==718.or.kstep==803) idiag=5
!-----------------------------------------------------------------------------
!  For each trial timestep, tdel, the Newton-Raphson iteration is attempted.
!  The possible results for a timestep are a converged set of yt or a
!  failure to converge.  If convergence fails, iteration is retried with
!  the trial timestep reduced by tdelmm, up to ktsmx times.
!-----------------------------------------------------------------------------
        TS: Do kts=1,ktsmx
!  Calculate the thermodynamic factors necessary for reaction rates,
!  including screening, and the reaction rates.
         If(gn>2)then
          call cross_sect
        Endif
!  The Newton-Raphson iteration occurs for at most knrmx iterations.
          NR: Do knr=1,knrmx

!  Calculate the changes in abundaces, dy
!            call qyderiv
!Build Jacobian
            call netmatr(kstep,elapsed_time,clock_rate)
         accum2= elapsed_time/clock_rate
         accum1= accum1+ accum2

!  Evolve the abundances and calculate convergence tests
!            yt=yt+dy
!            Where(yt<ymin)
!              yt=0.0
!              reldy=0.0
!            ElseWhere
!              reldy=abs(dy/yt)
!            EndWhere
!            If(idiag>=3) Then
!              irdymx=maxloc(reldy)
!              idymx=maxloc(dy)
!              Write(50,"(a2,i5,2i3,2(a5,2es12.4))")
!     &          'dY',kstep,kts,knr,nname(idymx(1)),
!     &          dy(idymx(1)),y(idymx(1)),nname(irdymx(1)),
!     &          reldy(irdymx(1)),y(irdymx(1))
!              If(idiag>=4) Write(50,"(a5,5es12.4)")
!     &          (nname(k),yt(k),dy(k),reldy(k),(aa(k)*dy(k)),
!     &          (aa(k)*yt(k)),k=1,ny)
!            Endif

!  Evolve the abundances and calculate convergence tests
            testc=0.0
!            testc2=0.0
            yqt=yqt+dy
            call update2(yqt,yt)
            Where(yqt<ymin)
              yqt=0.0
              reldy=0.0
            ElseWhere
              reldy=abs(dy/yqt)
            EndWhere
            If(idiag>=3) Then
              irdymx=maxloc(reldy)
              idymx=maxloc(dy)
              Write(50,"(a2,i5,2i3,2(a5,2es12.4))") !Still 299
     &          'dY',kstep,kts,knr,nname(idymx(1)),
     &          dy(idymx(1)),y(idymx(1)),nname(irdymx(1)),
     &          reldy(irdymx(1)),y(irdymx(1))
              If(idiag>=4) Write(50,"(a5,5es12.4)")
     &          (nname(k),yq(k),dy(k),reldy(k),(aa(k)*dy(k)),
     &          (aa(k)*yt(k)),k=1,gp0i)
            Endif
!-----------------------------------------------------------------------------
!  There are 3 included convergence tests: testc, which measures relative
!  changes, testc2 which measures total abundance changes, and testm
!  which tests mass conservation.
            testc=sum(reldy)
!            testc2=sum(aaq*reldy)
            xtot=sum(aaq*yqt)-1.0
            testm=xtot-xtoto ! 299 needs to be 82
            If(idiag>=2) Write(50,"(a3,i5,i3,3es14.6)")
     &        'KNR',kstep,knr,testm,testc,testc2
!-----------------------------------------------------------------------------
!  testc is the most stringent test, and requires the most iterations.
!  testm is the most lax, and therefore the fastest, often requiring only one
!  iteration.  Considering the uncertainties of the reaction rates, it is
!  doubtful that the increased precision of testc is truly increased
!  accuracy.
!-----------------------------------------------------------------------------
!  Ordinarily, test for true convergence
            If (iconvc/=0.or.tt>=tstop) Then
              testn=testc
!              testn=testc2
              toln=tolc

!  Otherwise, use mass conservation for convergence condition
            Else
              testn=testm
              toln=tolm
            Endif
            If(abs(testn)<=toln) Exit TS
          Enddo NR

!  If convergence is not achieved in knrmx iterations, reset abundances
!  and try again with the timestep reduced.
          If(idiag>=1) Write(50,*) 'TS Failure',knr,kts,xtot,testn
          tdel=tdel/tdelmm
          tt=t+tdel
          yqt=yq
          yt=y
          call t9rhofind(kstep,tt,j)
          call qse_coeff
          call update2(yqt,yt)
           WRITE(*,*) "convergenve not achieved testc=", testc
        Enddo TS

!  If convergence is successful, update time and abundances
        If(kts<ktsmx) Then
          If(idiag>=1) Write(50,"(a4,i5,i3,3es12.4)")
     &      'Conv',kstep,knr,xtot,testn,toln
          ta=t
          tdela=tdel
          t=tt
!          call norm(y)
           yo=y ! 299
           yqo=yq

           yq=yqt
 
!            Write(50,"(a)") 'delta Y'
!            Write(50,"(a5,4es12.4)") (nname(gp0number(k)),yq(k),yqo(k),
!     &        (yqt(k)-yq(k)),(tdel*ydot(k)),k=1,gp0i)

! The next bit is involved in updating y. It makes the groups and then
! puts yq in to the proper places for the single nuc.
           call update2(yq,yt)
! This is an experiment to elimated some errors in the consatnat case.
           Do i=1,ny
              if (yt(i)<1e-27)then
                  yt(i)=0
              Endif
           Enddo
           y=yt
          ye=sum(zz*y)
          xsum=sum(y*aa)
!          write(*,*) ye,xsum
          enold=enm
          call benuc(yt,enb,enm,ytot,ztot,atot,t)
          edot=-(enm-enold)/tdel
          If(itso==3)then
!          write(888,*) t, edot
          endif
          If(idiag>=1) Then
            Write(50,"(i5,5es14.7)") kstep,t,tdel,t9t,rhot,ye
            Write(50,"(5(a5,es11.4))") (nname(k),y(k),k=1,ny)
          Endif
          If(itso>0) call ts_output(kstep,(enm-en0),edot,kts,knr,is)
          countiter=countiter + knr

! check to see if you are running the appropraite number of groups
            If (t9t.le.t9minus.or.t9t.ge.t9plus)then
!            write(*,*) "I'm checking groups", numberofgroups
            adaptive=1
            call group_number(data_dir,t9t)
            t9plus=t9t+1
            t9minus=t9t-1
            adaptive=0
            endif
         If(t>=tstop) Then

            write(*,*) t,t9t,rhot,kstep,countiter,knr
            write(*,*) t, t9t,'end',gp0i
            If(itso==3)then
              write(80,*) t,t9t,rhot
              write(80,"(4(a5,es14.7,1x))") (nname(i),y(i),i=1,ny)
            endif

            Exit STEP
          Endif

!  If reduced timesteps fail to yield convergence, warn and exit
        Else
          Write(6,"(i5,4es12.4,2i3)") kstep,t,tdel,t9t,rhot,knr,kts
          Write(6,*) 'Timestep retrys fail after ',kts,' attempts'
          Exit STEP
        Endif
      Enddo STEP

!  Test that the stop time is reached
      If(t<tstop) Then
        Write(50,"(a,es12.4,a,es12.4,a)") 'Evolution incomplete!!!'
        Write(6,"(a,es12.4,a,es12.4,a)") 'Evolution stopped at time=',t,
     &  '.  Stop time (',tstop,') not reached!'
        Write(6,"(a,i6)")
     &    'Approximately',int((tstop-t)/tdel),'more steps needed'
      Endif
!  End Post Processing cycle
      call final_output(kstep,accum1)
         
      Return
      End subroutine qse_net

      subroutine timestep(kstep)
!===============================================================================
!  This routine calculates the trial timestep, based on the prior timestep.
!  For tdel >0, this calculation is based on the relative changes of the
!               abundances during the previous step in the evolution.
!           =0, the timestep is calculated using the time derivatives of
!               the abundances based on reaction rates.
!           <0, the timestep is held constant, tdel=-tdelo.
!  The timestep is further constrained to be no more than a factor tdelmm
!  larger than the timestep in the previous step.  There is also the
!  provision for limiting the timestep in the event that the thermodynamic
!  conditions are changing too rapidly.
!===============================================================================
      use nuclear_data
      use controls
      use conditions
      use qse_data
      use abundances
      use thermo_data
      integer, save :: ints(1),intso(1)   ! nucleus governing timestep
      integer :: i,kstep

      real(8), dimension(:), allocatable :: ydotoy,tdt
      real(8) :: changeth,changest,tdelo,t9old,rhold,dt,dtherm
      real(8) :: tdels,tdelm,tdeln,t_lim
      If(.not.allocated(ydotoy)) Allocate(ydotoy(gp0i),tdt(ny))
!  Retain old values of timestep and thermo and calculate remaining time
!      changeth=.01
      changeth=.012
      changest=100 ! involved only if inital timestep
      tdelo=tdel
      t9old=t9t
      rhold=rhot
      tdels=tstop-t
      dt=1.0e+20

!  If this is not the initial timestep, calculate timestep from changes
!  in last timestep.
      If(tdelo>0.0) Then
        tdelm=tdelmm*tdelo ! xnet
        Where(yq>ytime)
           ydotoy=abs((yq-yqo)/yq)!
        ElseWhere
          ydotoy=0.0
        EndWhere
        ints=maxloc(ydotoy)
        tdeln=changemx*tdelo/ydotoy(ints(1))
        tdel=min(tdeln,tdels,tdelm)

!  If this is an initial timestep, yo does not exist, so calculate
!  timestep from derivatives. If derivatives are zero, take a small step.
      Elseif(tdelo==0.0) Then
      If(gn>2)then
        tdelm=1.0e-2*tdels
        intso(1)=0
        call cross_sect
        call qyderiv
        Where(yq>ytime)
          ydotoy=abs(ydot/yq)
        ElseWhere
          ydotoy=0.0
        EndWhere
        ints=maxloc(ydotoy)
        If(ydotoy(ints(1)).ne.0) Then
          tdeln=changest*changemx/ydotoy(ints(1))
        Else
          tdeln=tdelm
        Endif
        tdel=min(tdels,tdeln,tdelm)
        write(*,*) "tdels,tdeln,tdelm",tdels,tdeln,tdelm

        Else
        tdel=1.0e-4*tdels
        Endif
!  Keep timestep constant
      Else
        tdel=-tdelo
      Endif

!  Diagnostic Output
      If(idiag>=1) Write(50,"(a4,i5,4es12.4)")
     &  'tdel',kstep,tdel,tdeln,tdelo,tdels
!     If(idiag>=2) Then
!       Write(50,"(a5,i4,2es12.4)")
!    &    (nname(k),k,y(k),ydotoy(k),k=1,ny)
!     Endif

!  Retain the index of the species setting the timestep
      If(ints(1)/=intso(1)) Then
        If(idiag>=1) Write(50,*) 'ITC ',nname(ints(1)),t,tdel
        intso=ints
      Endif

!  Limit timestep if Thermodynamic variation is too large
      Do i=1,10
        tt=t+tdel
        call t9rhofind(kstep,tt,j)
        If(t9old>0) Then
!o           dtherm=abs(rhot-rhold)/rhold !qse version
          dtherm=abs(t9t-t9old)/t9old+abs(rhot-rhold)/rhold
        Endif
        If(dtherm<changeth) Exit
        tdel=.5*tdel
        If(i==10) Write(6,*) 'Error in Thermo variations after ',i,
     &    'reductions',tdel,t9t,rhot
      Enddo

! Limit timestep to 1/10 the current time.
!        t_lim=.1*t
!        If (tdel> t_lim) then
!            tdel= t_lim
!        Endif


!     If(idiag>=1) Write(50,"(a5,i5,2es12.4)")
!    &  'T9del',kstep,tdel,dtherm
      Return
      End subroutine timestep

       subroutine qyderiv
!-----------------------------------------------------------------------------
!  This routine calculates time derivatives for each nuclear species.
!  This calculation is performed by looping over nuclei, and summing the
!  reaction rates for each reaction which effects that nucleus.
! This version makes 82 Ydots.
!-----------------------------------------------------------------------------
      use controls
      use qse_data
      use reac_rate_data
      use nuclear_data
      use abundances
      use cross_sect_data
!      use reac_rate_data
      integer :: i,j,i0,i1,la1,le1,la2,le2,la3,le3
      integer :: qn(4),qz(4)
      real(8) :: s1,s11,s2,s22,s3,s33,ratepart,summand
!  From the cross sections and the counting array, calculate the reaction
!  rates
      b1=a1*csect1(mu1)
      b2=a2*csect2(mu2)
      b3=a3*csect3(mu3)

!  Calculate Ydot for each nucleus, summing over the reactions which affect it.
          ydot=0.0
          ratepart=0.0
       qn(1)=0
       qn(2)=0
       qn(3)=nn(f3)
       qn(4)=nn(f4)
       qz(1)=0
       qz(2)=0
       qz(3)=zz(f3)
       qz(4)=zz(f4)

        la1=laq(1,2)
        le1=leq(1,2)
        Do i1=la1,le1
          ratepart=b1(i1)*yt(n11(i1))
          ydot(1)=ydot(1)+qnnum(i1)*ratepart
          ydot(2)=ydot(2)+qznum(i1)*ratepart
         Enddo
        Do i=3,gn
          la1=laq(1,i)
          le1=leq(1,i)
          Do i1=la1,le1
             ratepart=b1(i1)*yt(n11(i1))
             ydot(1)=ydot(1)+(qnnum(i1)-qn(i))*ratepart
             ydot(2)=ydot(2)+(qznum(i1)-qz(i))*ratepart
             ydot(i)=ydot(i)+ratepart
          Enddo
         Enddo

!  Sum over the reactions with 2 reactants
        la2=laq(2,2)
        le2=leq(2,2)
        Do i1=la2,le2
          ratepart=b2(i1)*yt(n21(i1))*yt(n22(i1))
          ydot(1)=ydot(1)+(qnnum2(i1))*ratepart
          ydot(2)=ydot(2)+(qznum2(i1))*ratepart
         Enddo
        Do i=3,gn
        la2=laq(2,i)
        le2=leq(2,i)
        Do i1=la2,le2
          ratepart=b2(i1)*yt(n21(i1))*yt(n22(i1))
          ydot(1)=ydot(1)+(qnnum2(i1)-qn(i))*ratepart
          ydot(2)=ydot(2)+(qznum2(i1)-qz(i))*ratepart
          ydot(i)=ydot(i)+ratepart
         Enddo
        Enddo
        la3=laq(3,2)
        le3=leq(3,2)
        Do i1=la3,le3
          ratepart=b3(i1)*yt(n31(i1))*yt(n32(i1))*yt(n33(i1))
          ydot(1)=ydot(1)+(qnnum3(i1))*ratepart
          ydot(2)=ydot(2)+(qznum3(i1))*ratepart
         Enddo
       Do i=3,gn
        la3=laq(3,i)
        le3=leq(3,i)
        Do i1=la3,le3
          ratepart=b3(i1)*yt(n31(i1))*yt(n32(i1))*yt(n33(i1))
          ydot(1)=ydot(1)+(qnnum3(i1)-qn(i))*ratepart
          ydot(2)=ydot(2)+(qznum3(i1)-qz(i))*ratepart
          ydot(i)=ydot(i)+ratepart
         Enddo
        Enddo
        If(idiag>=5) Write(50,"(a3,a6,i4)") 'NUC',nname(i0),i0

         Do i0=gn+1,gp0i
!  Sum over the reactions with 1 reactant
        la1=laq(1,i0)
        le1=leq(1,i0)
        s1=0.0
        Do i1=la1,le1
          s11=b1(i1)*yt(n11(i1))
          s1=s1+s11
          If(idiag>=5) Write(50,"(3a5,'  1  ',4es15.7)")
     &      (nname(n1i(j,mu1(i1))),j=1,3),b1(i1),yt(n11(i1)),
     &      s11,a1(i1)
        Enddo
        If(idiag>=5) Write(50,*) '1->',nname(i0),la1,le1,s1

!  Sum over the reactions with 2 reactants
        la2=laq(2,i0)
        le2=leq(2,i0)
        s2=0.0
        Do i1=la2,le2
          s22=b2(i1)*yt(n21(i1))*yt(n22(i1))
          s2=s2+s22
          If(idiag>=5) Write(50,"(4a5,4es15.7)")
     &      (nname(n2i(i,mu2(i1))),i=1,4),b2(i1),yt(n21(i1)),
     &      yt(n22(i1)),s22
        Enddo
        If(idiag>=5) Write(50,*) '2->',nname(i0),la2,le2,s2

 
!  Sum over the reactions with 3 reactants
        la3=laq(3,i0)
        le3=leq(3,i0)
        s3=0.0
        Do i1=la3,le3
          s33=b3(i1)*yt(n31(i1))*yt(n32(i1))*yt(n33(i1))
          s3=s3+s33
          If(idiag>=5) Write(50,"(3a5,'  3  ',5es12.4)")
     &      (nname(n3i(i,mu3(i1))),i=1,3),b3(i1),yt(n31(i1)),
     &      yt(n32(i1)),yt(n33(i1)),s33
        Enddo
        If(idiag>=5) Write(50,*) '3->',nname(i0),la3,le3,s3
!  Sum the 3 components of Ydot
        ydot(i0)=s1+s2+s3
        If(idiag>=5) Write(50,"(a4,a5,2es24.16)")
     &    'YDOT',nname(i0),yt(i0),ydot(i0)
      Enddo
      Return
      End



