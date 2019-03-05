!tt*******************************************************************************
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
      integer :: firstcaller
      integer, dimension(:), allocatable  :: gp1number,gp2number
      integer, dimension(:), allocatable  :: gp3number,gp0number
      integer, dimension(:), allocatable  :: gp4number
      integer  :: gp1i,gp2i,gp3i,gp0i,gp4i,i2,i3,i4,ig,i1
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
      integer, dimension(:), allocatable :: qmu1,qmu2,qmu3
      integer, dimension(:), allocatable :: qn11,qn21,qn22
      integer, dimension(:), allocatable :: qn31,qn32,qn33
      integer, dimension(:), allocatable :: qn12,qn13,qn14
      integer, dimension(:), allocatable :: qn23,qn24,qn25
      integer, dimension(:), allocatable :: qn34,qn35,qn36
      real(8),  dimension(:), allocatable :: qa1,qa2,qa3,qb1,qb2,qb3 ! dim(nan(i))
      real(8), dimension(:,:), allocatable :: part1,part1b,part2b
      real(8), dimension(:,:), allocatable :: part2i,part2,part2ib
      real(8), dimension(:), allocatable :: dlCqse
!     Jacobian stuff Note in xnet 5 these are in jacobi_---.f
      real(8), dimension(:), allocatable :: f,ydot0,um
      real(8), dimension(:,:), allocatable :: am ! the Jacobian Matrix
      real(8), dimension(:,:), allocatable :: dydot ! Part 1 ofJ Matrix
      integer, dimension(:), allocatable :: indx,indx4
      real suma ! qse faction by mass

      end module qse_data
      module qse_abundances
!===============================================================================
!  This module contains the abundances of the nuclear species, at the previous
!  time (yo), current time (y), and trial time (yt), as well as the time
!  derivatives (ydot), at the trial time, and the abundance change due to the
!  Newton-Raphson iteration.  The size of these arrays is allocated in the
!  external routine, where the initial values are set.
!  82 single nuc abundances at previous time(yqo),current time (yq),
!  and trial time (yqt).
!===============================================================================
      use nuc_number
      real(8), dimension(:), allocatable :: yqo,yq,yqt
!      integer, dimension(:), allocatable :: aaq,zzq
      end module qse_abundances

Module qjacobian_data
!===============================================================================
! The jacobian matrix for the solver.
!===============================================================================
  Use nuclear_data
  Real(8), Dimension(:,:), Allocatable :: jac !,jact ! the Jacobian Matrix
 ! Integer, Dimension(:), Allocatable :: indx
End Module qjacobian_data




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
!            Write(*,*) "group3 called"
       f1=1 !p
       f2=2 !n
       f3=43 !si28
       f4=104 !fe56

       gp1i=6
       If(.not.allocated(gp1number)) Allocate (gp1number(0:gp1i))
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
         Elseif(zz(i)>=22.and.zz(i)<=25.and.nn(i)>=21.and. &
     &           nn(i)<=23)then
            k=k+1
            dummysngl(k)=i
         Endif
! Set up the Si group

        If(zz(i)==12.and.nn(i)>=11)then
           j=j+1
           dummysi(j)=i
!           write(501,*) j,nname(dummysi(j)),'begin'
         ElseIf(zz(i)>=13.and.nn(i)>=10.and.zz(i)<=18.and. &
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
       If(.not.allocated(gp0number))Allocate (gp0number(0:gp0i))
       If(.not.allocated(gp2number))Allocate (gp2number(0:gp2i))
       If(.not.allocated(gp3number))Allocate (gp3number(0:gp3i))
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
!          write(*,*) "gp2",gp2number(i),dummysi(i)
           i3=gp2number(i)
           drvN(i3)=nn(i3)-nn(f3)
           drvZ(i3)=zz(i3)-zz(f3)
!           write(*,*) drvN(i3),nn(i3),nn(f3)
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
    If(.not.allocated(aaq))Allocate(aaq(gp0i))
       Do i=1,gp0i
           Do j=1,ny
             ig=gp0number(i)
              If(j==ig)then
                 aaq(i)=aa(j)
              Exit
          Endif
           Enddo
        Enddo
!      write(*,*) "I got to the end of group 3" 
       Return
       End ! subroutine group3

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
       If(.not.allocated(gp0number))Allocate (gp0number(0:gp0i))
       If(.not.allocated(gp2number))Allocate (gp2number(0:gp2i))

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
      If(.not.allocated(aaq))Allocate(aaq(gp0i))
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
       If(.not.allocated(aaq))Allocate(aaq(gp0i))
       If(.not.allocated(gp0number))Allocate (gp0number(0:gp0i))
       If(.not.allocated(gp2number))Allocate (gp2number(0:gp2i))

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
       If(.not.allocated(aaq))Allocate(aaq(gp0i))
       If(.not.allocated(gp0number))Allocate (gp0number(0:gp0i))
       If(.not.allocated(gp1number))Allocate (gp1number(0:gp1i))

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

       Return
       End !subroutine group1


subroutine group_number()
!c-----------------------------------------------------
! This subroutine sets up the nubmer of groups.
!c-----------------------------------------------------
       use qse_data
       use qse_abundances
       use abundances
!             write(*,*) "group_number was called "
If(.not.allocated(drvN)) Allocate(drvN(ny),drvZ(ny),ID(ny))
numberofgroups=3
gn=numberofgroups+1
            Write(*,*) "group3",numberofgroups
             call group3
             write(*,*) " I call group 3",gp0i,gp1i,gp2i,gp3i
             call qse_sort
  If(.not.allocated(indx4)) Allocate(indx4(gn:gn))
  If(.not.allocated(yqo)) Allocate(yqo(gp0i))
  If(.not.allocated(yq)) Allocate(yq(gp0i))
  If(.not.allocated(yqt)) Allocate(yqt(gp0i))

       Return
       End
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
      integer :: i,j,n,l,k,cn,cm,cn2,cn3
      integer :: i0,la1,le1,la2,le2,la3,le3,j1,l1,l2,l3
       character(len=5)  :: nnname
      character (LEN=1) :: null
      integer :: nr,idel,jdel,ijd,ndi,nds,plotindex(ny)
      common /ceign/ nds,ndi,idel,jdel,ijd
      j=0

! When the nubmer of groups changes in the middle of a run,
! the reactions must be resorted to suit the new groups
! Thus the massive deallocation below

!      If (adaptive>0)then
!        DeAllocate (singlenuc)
!        DeAllocate (laq,leq,laq1,leq1,laq2,leq2)
!        DeAllocate (laq3,leq3,part1,part2i,part2)
!        DeAllocate (part1b,part2ib,part2b)
!
!        DeAllocate (qmu1,qa1,qb1,qnnum,qznum)
!        DeAllocate (qn11,qn12,qn13,qn14,regrp1)
!
!        DeAllocate (qmu2,qa2,qb2,qnnum2,qznum2)
!        DeAllocate (qn21,qn22,qn23,qn24,qn25,regrp2)
!
!        DeAllocate (qmu3,qa3,qb3,qnnum3,qznum3)
!        DeAllocate (qn31,qn32,qn33,qn34,qn35,qn36)
!        DeAllocate(regrp3)

!       Endif
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
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!C Now we grab all the reaction in group1 and establish a laq
!C and leq for the whole group for 1 particle reactions.
!Cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! We do this loop twice. The first time we count
! how many reactions are needed for
! the light group. Then we allocate that many.
! qqmu1-contains the rxn  number
! qqa1-contains the doubble counting corrections
! qnnum-contains the number of neutrons in the lightgroup nucleus
!  that is related to the reactions qqmu1.
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
           If (n11(k)==gp0number(j).or.n12(k)==gp0number(j).or. &
     &        n13(k)==gp0number(j).or.n14(k)==gp0number(j).and. &
     &        iwk1(mu1(k)).le.0)then
              cn=cn+1
     !          write(100,*) n11(k),n12(k),n13(k),n14(k)
           Exit !so that you don't count twice for rxns with 2 singlenuc.
           Endif
        Enddo ! j
       Enddo ! k
! now for 2p rxns
        la2=la(2,ig)
        le2=le(2,ig)
        Do k=la2,le2
         Do j=1+numberofgroups+1,gp0i
           If (n21(k)==gp0number(j).or.n22(k)==gp0number(j).or.   &
     &        n23(k)==gp0number(j).or.                            &
     &        n24(k)==gp0number(j).or.n25(k)==gp0number(j).and. &
     &        iwk2(mu2(k)).le.0)then
              cn2=cn2+1
              ! write(100,*) n21(k),n22(k),n23(k),n24(k),n25(k)
           Exit !so that you don't count twice for rxns with 2 singlenuc.
           Endif
        Enddo ! j
       Enddo ! k
! Now for 3p reactions
        la3=la(3,ig)
        le3=le(3,ig)
        Do k=la3,le3
         Do j=1+numberofgroups+1,gp0i
           If (n31(k)==gp0number(j).or.n32(k)==gp0number(j).or.   &
     &        n33(k)==gp0number(j).or.n34(k)==gp0number(j).or.   &
     &        n35(k)==gp0number(j).or.n36(k)==gp0number(j).and. &
     &        iwk3(mu3(k)).le.0)then
              cn3=cn3+1
!               write(100,*) nname(n31(k)),nname(n32(k)),nname(n33(k)) &
!     & ,nname(n34(k)),nname(n35(k)),nname(n36(k))
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
           If (n11(k)==gp0number(j).or.n12(k)==gp0number(j).or.  &
     &        n13(k)==gp0number(j).or.n14(k)==gp0number(j).and. & 
     &        iwk1(mu1(k)).le.0)then
              cn=cn+1
           Exit ! so that you don't count twice for rxns with 2 singlenuc.
           Endif
        Enddo ! j
       Enddo ! k
        la2=la(2,ig)
        le2=le(2,ig)
        Do k=la2,le2
         Do j=1+numberofgroups+1,gp0i
           If (n21(k)==gp0number(j).or.n22(k)==gp0number(j).or.   &
     &        n23(k)==gp0number(j).or.                            &  
     &        n24(k)==gp0number(j).or.n25(k)==gp0number(j).and. &
     &        iwk2(mu2(k)).le.0)then                               
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
           If (n31(k)==gp0number(j).or.n32(k)==gp0number(j).or.   &
     &        n33(k)==gp0number(j).or.n34(k)==gp0number(j).or.   &
     &        n35(k)==gp0number(j).or.n36(k)==gp0number(j).and. &
     &        iwk3(mu3(k)).le.0)then
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
           If (n11(k)==gp0number(j).or.n12(k)==gp0number(j).or.  &
     &        n13(k)==gp0number(j).or.n14(k)==gp0number(j).and. & 
     &        iwk1(mu1(k)).le.0)then
              cn=cn+1
           Exit !so that you don't count twice for rxns with 2 singlenuc.
           Endif
        Enddo ! j
       Enddo ! k

        la2=la(2,ig)
        le2=le(2,ig)
        Do k=la2,le2
         Do j=1+numberofgroups+1,gp0i
           If (n21(k)==gp0number(j).or.n22(k)==gp0number(j).or.   &
     &        n23(k)==gp0number(j).or.                            &      
     &        n24(k)==gp0number(j).or.n25(k)==gp0number(j).and. &
     &        iwk2(mu2(k)).le.0)then
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
           If (n31(k)==gp0number(j).or.n32(k)==gp0number(j).or.   &
     &        n33(k)==gp0number(j).or.n34(k)==gp0number(j).or.   &
     &        n35(k)==gp0number(j).or.n36(k)==gp0number(j).and. &
     &        iwk3(mu3(k)).le.0)then
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
           If (n11(k)==gp0number(j).or.n12(k)==gp0number(j).or.  &
     &        n13(k)==gp0number(j).or.n14(k)==gp0number(j).and. &
     &        iwk1(mu1(k)).le.0)then
              cn=cn+1
           Exit !so that you don't count twice for rxns with 2 singlenuc.
           Endif
        Enddo ! j
       Enddo ! k

        la2=la(2,ig)
        le2=le(2,ig)
        Do k=la2,le2
         Do j=1+numberofgroups+1,gp0i
           If (n21(k)==gp0number(j).or.n22(k)==gp0number(j).or.   &
     &        n23(k)==gp0number(j).or.                            &
     &        n24(k)==gp0number(j).or.n25(k)==gp0number(j).and. &
     &        iwk2(mu2(k)).le.0)then
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
           If (n31(k)==gp0number(j).or.n32(k)==gp0number(j).or.   &
     &        n33(k)==gp0number(j).or.n34(k)==gp0number(j).or.   &
     &        n35(k)==gp0number(j).or.n36(k)==gp0number(j).and. & 
     &        iwk3(mu3(k)).le.0)then
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
           If(iwk1(mu1(k)).le.0)then
!        Do j=1+numberofgroups+1,gp0i
!           If (n11(k)==gp0number(j))then
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
           If(iwk2(mu2(k)).le.0)then
             cn2=cn2+1
        Endif
        Enddo !k
         leq(2,i)=cn2 !!!!!!!!!!!!!!!!!!!!!!!!
         laq(3,i)=cn3+1 !!!!!!!!!!!!!!!!!!!!!!!!!
        la3=la(3,ig)
        le3=le(3,ig)
        Do k=la3,le3
           If(iwk3(mu3(k)).le.0)then
             cn3=cn3+1
         Endif
        Enddo !k
         leq(3,i)=cn3 !!!!!!!!!!!!!!!!!!!!!!!!
       Enddo !i

!Allocate the varaibles
      Allocate (qmu1(cn))
      Allocate (qa1(cn),qb1(cn))
      Allocate (qnnum(cn))
      Allocate (qznum(cn))
      Allocate (qn11(cn),qn12(cn),qn13(cn),qn14(cn))
      Allocate(regrp1(cn,3))

      Allocate (qmu2(cn2))
      Allocate (qa2(cn2),qb2(cn2))
      Allocate (qnnum2(cn2))
      Allocate (qznum2(cn2))
      Allocate (qn21(cn2),qn22(cn2),qn23(cn2),qn24(cn2),qn25(cn2))
      Allocate(regrp2(cn2,3))

      Allocate (qmu3(cn3))
      Allocate (qa3(cn3),qb3(cn3))
      Allocate (qnnum3(cn3))
      Allocate (qznum3(cn3))
      Allocate (qn31(cn3),qn32(cn3),qn33(cn3),qn34(cn3), & 
     &         qn35(cn3),qn36(cn3))
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
             If (n11(k)==gp0number(j).or.n12(k)==gp0number(j).or.  &
     &          n13(k)==gp0number(j).or.n14(k)==gp0number(j).and. &
     &          iwk1(mu1(k)).le.0)then
                cn=cn+1
                qmu1(cn)=mu1(k)
                qa1(cn)=a1(k)
                qb1(cn)=b1(k)
                qnnum(cn)=int(nn(ig))
                qznum(cn)=int(zz(ig))
                qn11(cn)=n11(k)
                qn12(cn)=n12(k)
                qn13(cn)=n13(k)
                qn14(cn)=n14(k)
          exit
            Endif
          Enddo ! j
       Enddo ! k
! now for 2p rxns
        la2=la(2,ig)
        le2=le(2,ig)
        Do k=la2,le2
         Do j=1+numberofgroups+1,gp0i
           If (n21(k)==gp0number(j).or.n22(k)==gp0number(j).or.  &
     &        n23(k)==gp0number(j).or.                           &  
     &        n24(k)==gp0number(j).or.n25(k)==gp0number(j).and.&
     &        iwk2(mu2(k)).le.0)then
              cn2=cn2+1
                qmu2(cn2)=mu2(k)
                qa2(cn2)=a2(k)
                qb2(cn2)=b2(k)
                qnnum2(cn2)=int(nn(ig))
                qznum2(cn2)=int(zz(ig))
                qn21(cn2)=n21(k)
                qn22(cn2)=n22(k)
                qn23(cn2)=n23(k)
                qn24(cn2)=n24(k)
                qn25(cn2)=n25(k)
           Exit !so that you don't count twice for rxns with 2 singlenuc.
           Endif
        Enddo ! j
       Enddo ! k
! Now for 3p reactions
        la3=la(3,ig)
        le3=le(3,ig)
        Do k=la3,le3
         Do j=1+numberofgroups+1,gp0i
           If (n31(k)==gp0number(j).or.n32(k)==gp0number(j).or. &
     &        n33(k)==gp0number(j).or.n34(k)==gp0number(j).or. &
     &        n35(k)==gp0number(j).or.n36(k)==gp0number(j).and. &
     &        iwk3(mu3(k)).le.0)then
              cn3=cn3+1
               qmu3(cn3)=mu3(k)
                qa3(cn3)=a3(k)
                qb3(cn3)=b3(k)
                qnnum3(cn3)=int(nn(ig))
                qznum3(cn3)=int(zz(ig))
                qn31(cn3)=n31(k)
                qn32(cn3)=n32(k)
                qn33(cn3)=n33(k)
                qn34(cn3)=n34(k)
                qn35(cn3)=n35(k)
                qn36(cn3)=n36(k)

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
             If (n11(k)==gp0number(j).or.n12(k)==gp0number(j).or.  &
     &          n13(k)==gp0number(j).or.n14(k)==gp0number(j).and. &
     &          iwk1(mu1(k)).le.0)then
                cn=cn+1
                qmu1(cn)=mu1(k)
                qa1(cn)=a1(k)
                qb1(cn)=b1(k)
                qnnum(cn)=int(nn(ig)) !nn-14 to simplify!
                qznum(cn)=int(zz(ig))
                qn11(cn)=n11(k)
                qn12(cn)=n12(k)
                qn13(cn)=n13(k)
                qn14(cn)=n14(k)
          exit
            Endif

          Enddo ! j
       Enddo ! k
! now for 2p rxns
        la2=la(2,ig)
        le2=le(2,ig)
        Do k=la2,le2
         Do j=1+numberofgroups+1,gp0i
           If (n21(k)==gp0number(j).or.n22(k)==gp0number(j).or.   &
     &        n23(k)==gp0number(j).or.                            &  
     &        n24(k)==gp0number(j).or.n25(k)==gp0number(j).and. &
     &        iwk2(mu2(k)).le.0)then
              cn2=cn2+1
                qmu2(cn2)=mu2(k)
                qa2(cn2)=a2(k)
                qb2(cn2)=b2(k)
                qnnum2(cn2)=int(nn(ig))
                qznum2(cn2)=int(zz(ig))
                qn21(cn2)=n21(k)
                qn22(cn2)=n22(k)
                qn23(cn2)=n23(k)
                qn24(cn2)=n24(k)
                qn25(cn2)=n25(k)
           Exit !so that you don't count twice for rxns with 2 singlenuc.
           Endif
        Enddo ! j
       Enddo !k
! Now for 3p reactions
        la3=la(3,ig)
        le3=le(3,ig)
        Do k=la3,le3
         Do j=1+numberofgroups+1,gp0i
           If (n31(k)==gp0number(j).or.n32(k)==gp0number(j).or.   &
     &        n33(k)==gp0number(j).or.n34(k)==gp0number(j).or.   &
     &        n35(k)==gp0number(j).or.n36(k)==gp0number(j).and. &
     &        iwk2(mu3(k)).le.0)then
              cn3=cn3+1
                qmu3(cn3)=mu3(k)
                qa3(cn3)=a3(k)
                qb3(cn3)=b3(k)
                qnnum3(cn3)=int(nn(ig))
                qznum3(cn3)=int(zz(ig))
                qn31(cn3)=n31(k)
                qn32(cn3)=n32(k)
                qn33(cn3)=n33(k)
                qn34(cn3)=n34(k)
                qn35(cn3)=n35(k)
                qn36(cn3)=n36(k)

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
             If (n11(k)==gp0number(j).or.n12(k)==gp0number(j).or.  &
     &          n13(k)==gp0number(j).or.n14(k)==gp0number(j).and. &
     &          iwk1(mu1(k)).le.0)then
                cn=cn+1
                qmu1(cn)=mu1(k)
                qa1(cn)=a1(k)
                qnnum(cn)=int(nn(ig))
                qb1(cn)=b1(k)
                qznum(cn)=int(zz(ig))
                qn11(cn)=n11(k)
                qn12(cn)=n12(k)
                qn13(cn)=n13(k)
                qn14(cn)=n14(k)
          exit
            Endif
          Enddo ! j
       Enddo ! k
! now for 2p rxns
        la2=la(2,ig)
        le2=le(2,ig)
        Do k=la2,le2
         Do j=1+numberofgroups+1,gp0i
           If (n21(k)==gp0number(j).or.n22(k)==gp0number(j).or.  &
     &        n23(k)==gp0number(j).or.                           & 
     &        n24(k)==gp0number(j).or.n25(k)==gp0number(j).and. &
     &        iwk2(mu2(k)).le.0)then
              cn2=cn2+1
                qmu2(cn2)=mu2(k)
                qa2(cn2)=a2(k)
                qb2(cn2)=b2(k)
                qnnum2(cn2)=int(nn(ig))
                qznum2(cn2)=int(zz(ig))
                qn21(cn2)=n21(k)
                qn22(cn2)=n22(k)
                qn23(cn2)=n23(k)
                qn24(cn2)=n24(k)
                qn25(cn2)=n25(k)
           Exit !so that you don't count twice for rxns with 2 singlenuc.
           Endif
        Enddo ! j
        Enddo !k
! Now for 3p reactions
        la3=la(3,ig)
        le3=le(3,ig)
        Do k=la3,le3
         Do j=1+numberofgroups+1,gp0i
           If (n31(k)==gp0number(j).or.n32(k)==gp0number(j).or.  &
     &        n33(k)==gp0number(j).or.n34(k)==gp0number(j).or.  &
     &        n35(k)==gp0number(j).or.n36(k)==gp0number(j).and.& 
     &        iwk3(mu3(k)).le.0)then
              cn3=cn3+1
                qmu3(cn3)=mu3(k)
                qa3(cn3)=a3(k)
                qb3(cn3)=b3(k)
                qnnum3(cn3)=int(nn(ig))
                qznum3(cn3)=int(zz(ig))
                qn31(cn3)=n31(k)
                qn32(cn3)=n32(k)
                qn33(cn3)=n33(k)
                qn34(cn3)=n34(k)
                qn35(cn3)=n35(k)
                qn36(cn3)=n36(k)

           Exit !so that you don't count twice for rxns with 2 singlenuc.
           Endif
!graqb2 299
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
             If (n11(k)==gp0number(j).or.n12(k)==gp0number(j).or. &
     &          n13(k)==gp0number(j).or.n14(k)==gp0number(j).and.&
     &          iwk1(mu1(k)).le.0)then
                cn=cn+1
                qmu1(cn)=mu1(k)
                qa1(cn)=a1(k)
                qnnum(cn)=int(nn(ig))
                qb1(cn)=b1(k)
                qznum(cn)=int(zz(ig))
                qn11(cn)=n11(k)
                qn12(cn)=n12(k)
                qn13(cn)=n13(k)
                qn14(cn)=n14(k)
          exit
            Endif
          Enddo ! j
       Enddo ! k
! now for 2p rxns
        la2=la(2,ig)
        le2=le(2,ig)
        Do k=la2,le2
         Do j=1+numberofgroups+1,gp0i
           If (n21(k)==gp0number(j).or.n22(k)==gp0number(j).or.   &
     &        n23(k)==gp0number(j).or.                            &
     &        n24(k)==gp0number(j).or.n25(k)==gp0number(j).and. &
     &        iwk2(mu2(k)).le.0)then
              cn2=cn2+1
                qmu2(cn2)=mu2(k)
                qa2(cn2)=a2(k)
                qb2(cn2)=b2(k)
                qnnum2(cn2)=int(nn(ig))
                qznum2(cn2)=int(zz(ig))
                qn21(cn2)=n21(k)
                qn22(cn2)=n22(k)
                qn23(cn2)=n23(k)
                qn24(cn2)=n24(k)
                qn25(cn2)=n25(k)
           Exit !so that you don't count twice for rxns with 2 singlenuc.
           Endif
        Enddo ! j
        Enddo !k
! Now for 3p reactions
        la3=la(3,ig)
        le3=le(3,ig)
        Do k=la3,le3
         Do j=1+numberofgroups+1,gp0i
           If (n31(k)==gp0number(j).or.n32(k)==gp0number(j).or.   &
     &        n33(k)==gp0number(j).or.n34(k)==gp0number(j).or.   &
     &        n35(k)==gp0number(j).or.n36(k)==gp0number(j).and. &
     &        iwk3(mu3(k)).le.0)then
              cn3=cn3+1
                qmu3(cn3)=mu3(k)
                qa3(cn3)=a3(k)
                qb3(cn3)=b3(k)
                qnnum3(cn3)=int(nn(ig))
                qznum3(cn3)=int(zz(ig))
                qn31(cn3)=n31(k)
                qn32(cn3)=n32(k)
                qn33(cn3)=n33(k)
                qn34(cn3)=n34(k)
                qn35(cn3)=n35(k)
                qn36(cn3)=n36(k)

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
           If(iwk1(mu1(k)).le.0)then
             cn=cn+1
                qmu1(cn)=mu1(k)
                qa1(cn)=a1(k)
                qb1(cn)=b1(k)
                qnnum(cn)=int(nn(ig))
                qznum(cn)=int(zz(ig))
                qn11(cn)=n11(k)
                qn12(cn)=n12(k)
                qn13(cn)=n13(k)
                qn14(cn)=n14(k)
            Endif
        Enddo !j
! now for 2p rxns
        la2=la(2,ig)
        le2=le(2,ig)
        Do k=la2,le2
           If(iwk2(mu2(k)).le.0)then
              cn2=cn2+1
                qmu2(cn2)=mu2(k)
                qa2(cn2)=a2(k)
                qb2(cn2)=b2(k)
                qnnum2(cn2)=int(nn(ig))
                qznum2(cn2)=int(zz(ig))
                qn21(cn2)=n21(k)
                qn22(cn2)=n22(k)
                qn23(cn2)=n23(k)
                qn24(cn2)=n24(k)
                qn25(cn2)=n25(k)
             Endif
        Enddo !k
!And 3p rxns
        la3=la(3,ig)
        le3=le(3,ig)
        Do k=la3,le3
           If(iwk3(mu3(k)).le.0)then
             cn3=cn3+1
                qmu3(cn3)=mu3(k)
                qa3(cn3)=a3(k)
                qb3(cn3)=b3(k)
                qnnum3(cn3)=int(nn(ig))
                qznum3(cn3)=int(zz(ig))
                qn31(cn3)=n31(k)
                qn32(cn3)=n32(k)
                qn33(cn3)=n33(k)
                qn34(cn3)=n34(k)
                qn35(cn3)=n35(k)
                qn36(cn3)=n36(k)
            endif
        Enddo !k

      Enddo !i
!       write(300,*) qn11,qn12,qn13,qn14,qn21,qn22,qn23,qn24, &
!      &             qn25,qn31,qn32,qn33,qn34,qn35,qn36
        
      Return
      End
!-------------------------------------------------------------------------

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine update(ys)
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! In: The initial guess at the qse yf abundaces, this is called ys
! as in y used in the subroutine.
! Out: The proper starting qse abundances for yf.
! The goal is to use the derivative to get a normalized set of
! all abudnaces from  the 4 focal abundance in the full set.
! the single nuc remain untouched but the groups nuc get redone.
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      use nuclear_data
      use qse_data
      use controls
      use conditions
      use abundances
      use match_data
      use flux_data
      use part_funct_data
      use screen_qse
        integer :: i, iter_count1,j
        Real(8) :: yf(ny),y1,y2,y1n,y2n,ys(ny)
        Real(8) :: df(gn,gn),fin(gn),invdf(gn,gn),dummy(gn)
        Real(8) :: Cpart,fmax
        Real(8) :: ye,x28,x56,x184,ryn,ryp
        integer :: n77,z77,n186,z186,nsi,zsi,nfe,zfe,indxx(gn),info
! Get the initial yf and then use it with the QSE coeffients to rebuild the
! all the group nuclei. Put the new nuc in Y.
!       write(808,*) "ys",ys,"ys"
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
! This is the start of the loop that iterates until df converges

20       continue
!cccccccccccccccccccccccccccccccccccccccccccccc
        iter_count1=iter_count1+1
        If(iter_count1>20) then
      Write(*,*) "Initial abundances do not converge."
           stop
        EndIf

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! Make f f(1)=-1+ff(1),-ye+ff(2),-x28+ff(3),-x56+ff(4)
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
          call createf(ys,x28,x56,x184,ye,fin,fmax)
            If (fmax.lt.1d-13) goto 21
            If (fmax.gt.1D12) stop 'keine Konvergenz'

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! make df
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
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
             Cpart=dexp(dlCqse(i)-dlCqse(f3)+dlog(y1)*nsi+ &
     &             dlog(y2)*zsi)
   !    write(*,*) "he", i, he(int(zz(i)))
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
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! find inverse of df
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
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
      subroutine createf(yfin,x28,x56,x184,ye,ff,fmax1)
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c This makes the vector f that will be used to get yr back out of yg.
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      use nuclear_data
      use qse_data
      use controls
      use conditions
      use abundances
      use match_data
      use flux_data
      use part_funct_data



      real(8) yfin,x28,x56,x184,ff,fmax1,ye
      integer i,indxx,j

      dimension yfin(ny),ff(gn)

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c Generation of Vector ff                                                      c
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

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

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
       subroutine update2(yg,yf)
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c turns yg in to yf. This involves a Jacobian because yg(1-4) c
!c are group abudances, not n, p, si and fe.                    c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
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
       real(8) :: ryn, ryp,Cpart
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
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c Calculate the vector f.
!c The first time through sumng=Yg(1) sumpg=Yg(2) etc. so f=0
!c After that the sums use the adjusted yfs so f n.e. 0
!c This porcess cranks though the Jacobian solve with each new yf
!c until  fmax=dabs(ff(i)/yg(i)) is less than 1e-8.
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
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

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c Calculate Jacobian df                                                       c
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
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
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c update the other abundances given in yg(5) until yg(ni)                     c
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

21      continue
           Do ig=gn+1,gp0i
              i=gp0number(ig)

              yf(i)=yg(ig)
           Enddo

       Return
        End
!******************************************************************************

















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
      real(8) :: y1,y2,y3,y4                  
      integer :: i,j,k,n1,n2,n3,z1,z2,z3,z4,n4
      Real(8) :: tmpa,yf(ny),Fesum,t9oldy
!write(*,*) "in qse_make"
! write(809,*) "1",yf
         y1=yf(f1)
         y2=yf(f2)
!         write(809,*) y1,y2, "y1,y2"
!         write(809,*) "dlCqse",t9t,rhot, dlCqse

         Do i=1,gp1i
            ig=gp1number(i)
            yf(ig)=dexp(dlCqse(ig)+dlog(y2)*(zz(ig)) &
     &             +dlog(y1)*(nn(ig)))
         Enddo
! write(809,*) "2",yf
         If(numberofgroups>=2)then
           n3=nn(f3)
           z3=zz(f3)
           y3=yf(f3)
           Do i=1,gp2i
             ig=gp2number(i)
             yf(ig)=y3*dexp(dlCqse(ig)-dlCqse(f3)+dlog(y2)*(zz(ig)-z3) &
     &              +dlog(y1)*(nn(ig)-n3))
           Enddo
         Endif

         If(numberofgroups>=3)then
           n4=nn(f4)
           z4=zz(f4)
           y4=yf(f4)
            Do i=1,gp3i
              ig=gp3number(i)
              yf(ig)=y4*dexp(dlCqse(ig)-dlCqse(f4)+dlog(y2)*(zz(ig)-z4) &
     &               +dlog(y1)*(nn(ig)-n4))
            Enddo
          Endif

!write(*,*) "end qse_make"
        Return
        End

      subroutine qse_coeff
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
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
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      use nuclear_data
      use qse_data
      use controls
      use conditions
      use abundances
      use match_data
      use flux_data
      use part_funct_data
      use screen_qse

      real(8) :: xg(7),enuc,edot,tmp,bkt,dummy
      integer :: i,j,k,kstep,kout,zz1(ny),aa1(ny),nn1(ny),is,l,iout,ieos
      real, parameter :: pi=3.14159,hbr=6.58217D-22,amu=1.036435E-18
      real, parameter :: bok=8.6174D-02,avn=6.02205D+23
      real(8) ::  tmpa,Yqse(ny),Yqse2(ny),Fesum,t9oldy,c(ny)
      Real(8) :: ene,ye,yps,y2e,v,pa,ea,dlma,pta,pva,eta,eva,emas,beta

! Call EOS, necessary for screening
! write(300,*) "en from QSE"
! Call en(ny,y,ye,yps,y2e)
!        write(*,*) ye,yps,y2e,t9t
!  If(iscrn>=1) Then
!    v=1./rhot
!    Call state(t9t,v,y,pa,ea,dlma,ye,pta,pva,eta,eva,ieos,ny,emas,yps,beta)
!  EndIf
!      write(300,*) t9t,rhot,ye
!!       ye=sum(zz*y)
!      If(iscrn>0) Then
!        If(t9t>1.0) Then
!          write(*,*) "screen_called"
!          call qse_screen(t9t,rhot,ye,iout)
!          write(*,*) he(int(zz(77))), nname(77), int(zz(77)), "hesi"
!        Else
!          call qse_screen2(t9t,rhot,ye,iout)
!        Endif
!      Else
!        he=0.0
!      Endif

         call partf(t9t)



         bkt=bok*t9t
         tmp=dlog(((avn*rhot)**(2./3.)*2.*pi*hbr**2)/(bkt*amu))
!      write(*,*) 'Const',rhot,avn,hbr,pi,bkt,tmp
         dlCqse(1)=0.
         dlCqse(2)=0.
         DO  l=3,ny
            tmpa=dlog(angm(l)*gg(l)*(aa(l)**1.5)*(0.5**aa(l)))
            dlCqse(l)=tmpa+(1.5*(aa(l)-1))*tmp+(be(l)/bkt) + he(int(zz(l)))
         ENDDO
        Return
        End

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
      use qse_abundances
      use cross_sect_data
!      use reac_rate_data
      integer :: i,j,i0,la1,le1,la2,le2,la3,le3
      integer :: qn(4),qz(4)
      real(8) :: s1,s11,s2,s22,s3,s33,ratepart,summand
!  From the cross sections and the counting array, calculate the reaction
!  rates
      qb1=qa1*csect1(qmu1)
      qb2=qa2*csect2(qmu2)
      qb3=qa3*csect3(qmu3)

!      write(801,*) qa1
!      write(802,*) qmu1
!      write(803,*) csect1
!      write(804,*) yt 
!      write(805,*) qb1
!  Calculate idot for each nucleus, summing over the reactions which affect it.
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
          ratepart=qb1(i1)*yt(qn11(i1))
          ydot(1)=ydot(1)+qnnum(i1)*ratepart
          ydot(2)=ydot(2)+qznum(i1)*ratepart
         Enddo
!      write(800,*) qb1
        Do i=3,gn
          la1=laq(1,i)
          le1=leq(1,i)
          Do i1=la1,le1
             ratepart=qb1(i1)*yt(qn11(i1))
             ydot(1)=ydot(1)+(qnnum(i1)-qn(i))*ratepart
             ydot(2)=ydot(2)+(qznum(i1)-qz(i))*ratepart
             ydot(i)=ydot(i)+ratepart
          Enddo
         Enddo

!  Sum over the reactions with 2 reactants
        la2=laq(2,2)
        le2=leq(2,2)
        Do i1=la2,le2
          ratepart=qb2(i1)*yt(qn21(i1))*yt(qn22(i1))
          ydot(1)=ydot(1)+(qnnum2(i1))*ratepart
          ydot(2)=ydot(2)+(qznum2(i1))*ratepart
         Enddo
        Do i=3,gn
        la2=laq(2,i)
        le2=leq(2,i)
        Do i1=la2,le2
          ratepart=qb2(i1)*yt(qn21(i1))*yt(qn22(i1))
          ydot(1)=ydot(1)+(qnnum2(i1)-qn(i))*ratepart
          ydot(2)=ydot(2)+(qznum2(i1)-qz(i))*ratepart
          ydot(i)=ydot(i)+ratepart
         Enddo
        Enddo
        la3=laq(3,2)
        le3=leq(3,2)
        Do i1=la3,le3
          ratepart=qb3(i1)*yt(qn31(i1))*yt(qn32(i1))*yt(qn33(i1))
          ydot(1)=ydot(1)+(qnnum3(i1))*ratepart
          ydot(2)=ydot(2)+(qznum3(i1))*ratepart
         Enddo
       Do i=3,gn
        la3=laq(3,i)
        le3=leq(3,i)
        Do i1=la3,le3
          ratepart=qb3(i1)*yt(qn31(i1))*yt(qn32(i1))*yt(qn33(i1))
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
          s11=qb1(i1)*yt(qn11(i1))
          s1=s1+s11
          If(idiag>=5) Write(50,"(3a5,'  1  ',4es15.7)") &
     &      (nname(n1i(j,qmu1(i1))),j=1,3),qb1(i1),yt(qn11(i1)), &
     &      s11,qa1(i1)
        Enddo
        If(idiag>=5) Write(50,*) '1->',nname(i0),la1,le1,s1

!  Sum over the reactions with 2 reactants
        la2=laq(2,i0)
        le2=leq(2,i0)
        s2=0.0
        Do i1=la2,le2
          s22=qb2(i1)*yt(qn21(i1))*yt(qn22(i1))
          s2=s2+s22
          If(idiag>=5) Write(50,"(4a5,4es15.7)") &
     &      (nname(n2i(i,qmu2(i1))),i=1,4),qb2(i1),yt(qn21(i1)), &
     &      yt(qn22(i1)),s22
        Enddo
        If(idiag>=5) Write(50,*) '2->',nname(i0),la2,le2,s2


!  Sum over the reactions with 3 reactants
        la3=laq(3,i0)
        le3=leq(3,i0)
        s3=0.0
        Do i1=la3,le3
          s33=qb3(i1)*yt(qn31(i1))*yt(qn32(i1))*yt(qn33(i1))
          s3=s3+s33
          If(idiag>=5) Write(50,"(3a5,'  3  ',5es12.4)") &
     &      (nname(n3i(i,qmu3(i1))),i=1,3),qb3(i1),yt(qn31(i1)), &
     &      yt(qn32(i1)),yt(qn33(i1)),s33 
        Enddo
        If(idiag>=5) Write(50,*) '3->',nname(i0),la3,le3,s3
!  Sum the 3 components of Ydot
        ydot(i0)=s1+s2+s3
        If(idiag>=5) Write(50,"(a4,a5,2es24.16)") &
     &    'YDOT',nname(i0),yt(i0),ydot(i0)
      Enddo
      Return
      End

Subroutine the_great_deallocater()
      
!****************************************************************************
! It deallocates so the groups can switch
!----------------------------------------------------------------------------
      use nuclear_data
      use qse_data
      use qjacobian_data
      use reac_rate_data
      use ffn_data 
      use cross_sect_data
      use controls
      use abundances
      Use conditions
      write(*,*) "the great deallocater"

      If (isolv==3)then
        DeAllocate (singlenuc)
        DeAllocate (laq,leq,laq1,leq1,laq2,leq2)
        DeAllocate (laq3,leq3,part1,part2i,part2)
        DeAllocate (part1b,part2ib,part2b)

        DeAllocate (qmu1,qa1,qb1,qnnum,qznum)
        DeAllocate (qn11,qn12,qn13,qn14,regrp1)

        DeAllocate (qmu2,qa2,qb2,qnnum2,qznum2)
        DeAllocate (qn21,qn22,qn23,qn24,qn25,regrp2)

        DeAllocate (qmu3,qa3,qb3,qnnum3,qznum3)
        DeAllocate (qn31,qn32,qn33,qn34,qn35,qn36)
        DeAllocate(regrp3)
        DeAllocate(jac,indx)
      else


      endif
    End

Subroutine switcher(kstep)

!****************************************************************************
! It deallocates so the groups can switch
!----------------------------------------------------------------------------
      Use conditions
      use controls
      use nuclear_data
      use qse_data
      use qjacobian_data
      use reac_rate_data
      use ffn_data
      use cross_sect_data
      use abundances
     Integer :: izone,kstep,idiag0,its
     Real sumsi, sumfe
     integer nsolv,i
     logical :: firstCall
       Do i =43,89
          sumsi=sumsi+ y(i)*aa(i)
       Enddo
       Do i =89,ny
          sumfe=sumfe+ y(i)*aa(i)
       Enddo
      If (t9t>=3.5.and.sumfe.gt.0.3)then
         nsolv=3
      Else
        nsolv=1
      Endif
!        write(*,*) "isolv,nsolv",isolv,nsolv ,nname(43),y(43)
     If(nsolv.ne.isolv)then
 !     write(*,*) "switch si28=", y(43),nsolv,isolv
      Call the_great_deallocater()
        firstcaller=0
        isolv=nsolv
! Take Integration step
    Select Case (isolv)
      Case(3)
        firstcaller=1
        Call solve_qse(kstep,its)
      Case(2)
        Call solve_bd(kstep,its)
      Case Default
!        write(*,*) "full called"
        Call solve_be(kstep,its)
    End Select
    
     endif
     Return
     END

Subroutine group_ab_test(kstep)
      use nuclear_data
      use qse_data
      use controls
      use conditions
      use abundances
      use match_data
      use flux_data
      use part_funct_data

       real(8) :: y1,y2,y1n,y2n,fmax,yf(ny),dummy12
       real(8) :: ryn, ryp,Cpart
       real(8) :: sumng,sumpg,sumsig,sumfeg,dummy(gn,gn),sumcrg
       integer :: nsi,zsi,nfe,zfe,ncr,zcr,indxx(gn)
       integer :: i,counter,info,kstep
        counter=0
         counter=counter+1
         fmax=0
         sumsig=0
         sumfeg=0
         sumng=0
         sumpg=0
         sumcrg=0

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c Calculate the vector f.
!c The first time through sumng=Yg(1) sumpg=Yg(2) etc. so f=0
!c After that the sums use the adjusted yfs so f n.e. 0
!c This porcess cranks though the Jacobian solve with each new yf
!c until  fmax=dabs(ff(i)/yg(i)) is less than 1e-8.
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
         Do i=1, gp1i
            sumng=sumng+nn(i)*y(i)
            sumpg=sumpg+zz(i)*y(i)
          Enddo
        Do ig=1,gp2i
          i=gp2number(ig)
          sumng=sumng+(nn(i)-nn(f3))*y(i)
          sumpg=sumpg+(zz(i)-zz(f3))*y(i)
          sumsig=sumsig+y(i)
        Enddo
        Do ig=1,gp3i
          i=gp3number(ig)
          sumng=sumng+(nn(i)-nn(f4))*y(i)
          sumpg=sumpg+(zz(i)-zz(f4))*y(i)
          sumfeg=sumfeg+y(i)
        Enddo
!        write(807,*) "step", kstep
!        write(807,*) sumng,sumpg,sumsig,sumfeg

  Return
  End
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine qse_fraction(ys)
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! In: The initial guess at the qse yf abundaces, this is called ys
! as in y used in the subroutine.
! Out: The proper starting qse abundances for yf.
! The goal is to use the derivative to get a normalized set of
! all abudnaces from  the 4 focal abundance in the full set.
! the single nuc remain untouched but the groups nuc get redone.
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      use nuclear_data
      use qse_data
      use controls
      use conditions
      use abundances
      use match_data
      use flux_data
      use part_funct_data
      use screen_qse
        integer :: i, iter_count1,j
        Real(8) :: yf(ny),y1,y2,y1n,y2n,ys(ny)
        Real(8) :: df(gn,gn),fin(gn),invdf(gn,gn),dummy(gn)
        Real(8) :: Cpart,fmax
        Real(8) :: xpn,x28,x56,x184,ryn,ryp
        integer :: n77,z77,n186,z186,nsi,zsi,nfe,zfe,indxx(gn),info
! Get the initial yf and then use it with the QSE coeffients to rebuild the
! all the group nuclei. Put the new nuc in Y.
        x28=0.0
        xpn=0.0
        x56=0.0
        x184=0.0
         suma=0.000
   
        Do ig=1,gp1i
           i=gp1number(ig)
           xpn=xpn+ys(i)*aa(i)
           write(999,*) gp1number(ig)
        Enddo
        write(999,*) "si"
      If(numberofgroups>=2)then
        Do ig=1,gp2i
           i=gp2number(ig)
           x28=x28+ys(i)*aa(i)
           write(999,*) gp2number(ig)
        Enddo
        write(999,*) "fe"
 
      Endif
      If(numberofgroups>=3)then
        Do ig=1,gp3i
          i=gp3number(ig)
          x56=x56+ys(i)*aa(i)
         
           write(999,*) gp3number(ig)
        Enddo
      Endif
         suma=xpn+x28+x56
        write(*,*) "total mass fration in QSE groups =",suma
        write(*,*) xpn,x28,x56
      Return
     end
