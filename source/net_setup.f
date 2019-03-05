!===============================================================================
!  Net Setup 4.10 2/26/07
!  This program takes the human readable, ASCII data files sunet, netsu, 
!  netweak and netwinv and prepares binary versions for the reaction network. 
!  This involves translating characters to indicies and organizing the data  
!  for more efficient computation.  Several additional tests are performed.
!===============================================================================

      module reaction_data
!-----------------------------------------------------------------------------
!  This module contains the reaction rate data, including the indices
!  which type the reaction and the cross-references.
!-----------------------------------------------------------------------------
      use nuc_number
      integer, parameter :: nrm1=20000,nrm2=50000,nrm3=100
      integer            :: nnucr(8)=(/2,3,4,3,4,5,6,5/)
      integer            :: nreactant(8)=(/1,1,1,2,2,2,2,3/)
      integer            :: la(8),le(8)
      real(8)            ::    q1(nrm1),   q2(nrm2),   q3(nrm3)
      integer            ::  iwk1(nrm1), iwk2(nrm2), iwk3(nrm3)
      integer            :: ires1(nrm1),ires2(nrm2),ires3(nrm3)
      integer            :: irev1(nrm1),irev2(nrm2),irev3(nrm3)
      character(4)       :: desc1(nrm1),desc2(nrm2),desc3(nrm3)
      integer            ::  n1(4,nrm1), n2(6,nrm2), n3(5,nrm3)
      real(8)            :: an1(4,nrm1),an2(6,nrm2),an3(5,nrm3)
      integer            :: mu1(4*nrm1),mu2(6*nrm2),mu3(5*nrm3)
      real(8)            ::  a1(4*nrm1), a2(6*nrm2),a3(5*nrm3)
      integer, dimension(:), allocatable :: l1a,l2a,l3a
      integer, dimension(:), allocatable :: l1e,l2e,l3e
      end module reaction_data

      module flux_data
!-----------------------------------------------------------------------------
!  This module contains the 
!-----------------------------------------------------------------------------
      integer,parameter :: mxflx=50000
      integer      :: nflx(7,mxflx),iwflx(mxflx)
      real(8)      :: qflx(mxflx)
      character(4) :: descx(mxflx)
      end module flux_data

      program net_setup
!------------------------------------------------------------------------
!  This is the outer shell of Net Setup.  This routine also reformats 
!  the reaction data.
!------------------------------------------------------------------------
      use nuclear_data
      use reaction_data
!-------------------------------------------------------------------            
      integer n,i,j,k,l,ni(6),ki(6)
      integer, parameter :: ndesc=43
      character(LEN=5) :: name(6),blank5='     '
      character(LEN=4) desc
      character(LEN=4), dimension(ndesc) :: desc_known =
     &  (/'  ec','bet+','betp','bet-','betm','btyk',' bec','bkmo',
     &    'mo92',' ffn',' lmp','fkth','cf88','wies','laur','bb92',
     &    'baka','rolf','wawo','mafo','ms01',' wag','bb93','ca85',
     &    'ornl',' wfh','rath','rpsm','hg94','ro95','wien','iucf',
     &    'smwa','smmo','wi2p','aexp','pexp','bex+','bex-','bqa+',
     &    'bhi+','bec ','bklq'/)
      character(LEN=1) :: nr,vw,nrr='n',vvw='v',blank=' '
      character(80) data_desc
      real(8) p(7),q
      Open(unit=2,file='netsu')
      Open(3,file='nets3',form='unformatted')
      Open(4,file='nets4',form='unformatted')
      Open(14,file='net_diag')

!  Read and Reformat the nuclear data
      call format_nuclear_data
      Write(6,"(10a6)") (nname(i),i=1,ny)
      Write(6,"(1x,'n=',i4)") ny

!  Read and Reformat the reaction data
      Read(2,"(i5)") nm
      n=1
      Do jj=1,100000

!  Read in each entry from the reaction data file
        Read(2,"(i1,4x,6a5,8x,a4,a1,a1,3x,1pe12.5)",end=120) 
     &       k,(name(j),j=1,6),desc,nr,vw,q
        Read(2,"(4e13.6)",end=120) (p(j),j=1,7)

!  If the entry is a reaction
        If(k==0) Then

!  Convert nuclear names to indicies
          Do i=1,6
            call index_from_name(name(i),ni(i))
          Enddo
          
          if (any(ni>ny)) then
             write(6,*)'Dropping reaction ',name,
     &         'Some nuclei not included in netsu.'
            cycle
          endif
          
! Identify the reaction source against the list of known sources.
          Do i=1,ndesc
            If(desc==desc_known(i)) Exit
          Enddo
          If(i>ndesc) Then
            Write(6,'(3a)') 'Reaction Desc. ',desc,' not found!'
!           Stop
            idesc=i
          Else
            idesc=i
          Endif           
!-------------------------------------------------------------------            
!  The reaction rates are calculated via a number of templates, depending
!  on the type of reaction and its source.  In this section, flags are 
!  set for correct treatment of these different types of reactions
!-------------------------------------------------------------------  
!  For weak interactions; ,  (iec=1), Fuller, Fowler, Neuman Rates (iec=2,3),
!  Other weak rates not covered above (iec=4)
          Select Case (idesc)
            Case Default 
              iec=0 ! not a weak rate 
            Case (1) ! electron capture rates from Fowler and Co.
              iec=1
            Case (10,11) ! FFN type tabulated rate
              If(zz(ni(1))>zz(ni(2)).and.zz(ni(1))/=1) Then
                iec=3
              Else
                iec=2
              Endif
            Case(2:9,38:43) ! beta decay and non-tabulated ec
               iec=4
          End Select 

!  Check resonant nature of rate, resonant rate (ier=2), nonresonant (ier=1),
!  or none of the above (ier=0)
          If(nr==nrr) Then
            ier=1
          Elseif(nr==blank) Then
            ier=0
          Else
            ier=2
          Endif

!  Check for inverse rate (iev=1). Inverse rates must be multiplied with 
!  Temperature dependant ratio of partition functions. iev=0 is a forward rate.
          If(vw==vvw) Then
            iev=1
          Else
            iev=0
          Endif

!  Divide nuclei name list into reactants (ki=-1), products (ki=+1) and
!  blanks (ki=0)
          Do i=1,6
            If(i>nreactant(krt)) Then
              If(name(i)==blank5) Then
                ki(i)=0
              Else
                ki(i)=1
              Endif
            Else
              ki(i)=-1
            Endif
          Enddo

!  Count identical reactants
          kl=nreactant(krt)
          j=1
          Do i=2,kl
            If(name(i-1)/=name(i)) Then
              j=i
            Else
              ki(j)=ki(j)+ki(i)
              ki(i)=0
            Endif
          Enddo

!  Count identical products
          kl=nreactant(krt)+2
          j=kl-1
          Do i=kl,6
            If(name(i-1)/=name(i)) Then
              j=i
            Else
              ki(j)=ki(j)+ki(i)
              ki(i)=0
            Endif
          Enddo

!  Calculate double counting corrections and output reaction data to nets3
          If(nreactant(krt)==1) Then
            r=1.0
            Do i=1,nnucr(krt)
              n1(i,n)=ni(i)
              an1(i,n)=dble(ki(i))*r
            Enddo
            q1(n)=q
            iwk1(n)=iec
            ires1(n)=ier
            irev1(n)=iev
            desc1(n)=desc
            write(3) n,(ni(i),i=1,3),iec,ier,iev,p,q
!           write(21,*) n,(nname(ni(i)),i=1,3),iec,ier,iev,p,q
          Else
            lo=0
            Do i=1,nreactant(krt)
              If(ki(i)<lo) lo=ki(i)
            Enddo
            r=1.0/float(ifactorial(-lo))
            If(nreactant(krt)==2) Then
              Do i=1,nnucr(krt)
                n2(i,n)=ni(i)
                an2(i,n)=dble(ki(i))*r
              Enddo
              q2(n)=q
              iwk2(n)=iec
              ires2(n)=ier
              irev2(n)=iev
              desc2(n)=desc
              write(3) n,(ni(i),i=1,4),iec,ier,iev,p,q
!             write(21,*) n,(ni(i),i=1,4),iec,ier,iev,p,q
            Elseif(nreactant(krt)==3) Then
              Do i=1,nnucr(krt)
                n3(i,n)=ni(i)
                an3(i,n)=dble(ki(i))*r
              Enddo
              q3(n)=q
              iwk3(n)=iec
              ires3(n)=ier
              irev3(n)=iev
              desc3(n)=desc
              write(3) n,(ni(i),i=1,3),iec,ier,iev,p,q
!             write(21,*) n,(ni(i),i=1,3),iec,ier,iev,p,q
            Endif
          Endif
        n=n+1

!  If data entry is not a reaction, it is a type or sub-type marker
!  For type marker, reset type and sub-type counters
        Elseif(k==1) Then
          krt=k
          la(k)=1
          n=1
        Elseif(k==4.or.k==8) Then
          krt=k
          la(k)=1
          le(k-1)=n-1
          n=1

!  For sub-type marker, restart sub-type counters
        Else
          krt=k
          la(k)=n
          le(k-1)=n-1
        Endif
!       Write(14,"(1x,i2,6i4,6i2,3i2,2e10.4/1x,7e10.4)") 
!    &        k,ni(:),ki(:),iec,ier,iev,r,q,p(:)
      Enddo
  120 continue
      le(8)=n-1
      Close(1)
      Close(2)
!----------------------------------------------------------------------------
!  In order to build the Jacobian efficiently (row or column wise rather
!  than scattered), information about which reactions affect each nucleus 
!  is necessary.  
!----------------------------------------------------------------------------
!  For each nucleus, register each single reactant reaction which involves it.
      Allocate(l1a(ny),l1e(ny),l2a(ny),l2e(ny),l3a(ny),l3e(ny))
      n=0
      Do i=1,ny
        l1a(i)=n+1
        Do j=1,3 ; Do k=la(j),le(j) ;  Do l=1,nnucr(j)
          If(n1(l,k)==i.and.an1(l,k)/=0.) Then
            n=n+1
            a1(n)=an1(l,k)
            mu1(n)=k
            write(3) a1(n),mu1(n)
!           write(21,*) a1(n),mu1(n)
          Endif
        Enddo ; Enddo ; Enddo
        l1e(i)=n
      Enddo

!  For each nucleus, register each two reactant reaction which involves it.
      n=0
      Do i=1,ny
        l2a(i)=n+1
        Do j=4,7 ;  Do k=la(j),le(j) ; Do l=1,nnucr(j)
          If(n2(l,k)==i.and.an2(l,k)/=0.) Then
              n=n+1
              a2(n)=an2(l,k)
              mu2(n)=k
              write(3) a2(n),mu2(n)
!             write(21,*) a2(n),mu2(n)
          Endif
        Enddo ; Enddo ; Enddo
        l2e(i)=n
      Enddo

!  For each nucleus, register each three reactant reaction which involves it.
      n=0
      Do i=1,ny
        l3a(i)=n+1
        Do k=la(8),le(8) ; Do l=1,nnucr(8)
          If(n3(l,k)==i.and.an3(l,k)/=0.) Then
              n=n+1
              a3(n)=an3(l,k)
              mu3(n)=k
              write(3) a3(n),mu3(n)
!             write(21,*) a3(n),mu3(n)
          Endif
        Enddo ; Enddo
        l3e(i)=n
      Enddo
      Close(3)

!  Output the nucleus to reaction registration to nets4
      write(4) ny
      write(4) (nname(i),i=1,ny)
      write(4) nm
      write(4) le(3),le(7),le(8)
      Do i=1,ny
        write(4) i,l1a(i),l1e(i),l2a(i),l2e(i),l3a(i),l3e(i)
      Enddo
      Write(6,"(10x,'number of reactions of different types')")
      Write(6,"(3(10x,i8))") (j,la(j),le(j),j=1,8)
      Write(6,"(10x,'necessary dimensions')")
      Write(6,"(3(10x,i8))") le(3),le(7),le(8),l1e(ny),l2e(ny),l3e(ny)
      Write(6,"(3(10x,i8))") nm

!  Prepare network description file
      Write(6,*) "Provide a one line description of this data set"
      Read(5,"(a80)") data_desc
      Open(11,file='net_desc')
      Write(11,"(a80)") data_desc
      Write(11,"('Number of Nuclear Species=',i5)") ny
      Write(11,"('Reactions Numbers for the 8 different types')")
      Write(11,"(3i10)") (j,la(j),le(j),j=1,8)
      Write(11,"('Necessary dimensions')")
      Write(11,"(a15,3i8)") 'nreac(1,2,3)=',le(3),le(7),le(8)
      Write(11,"(a15,3i8)") 'nan(1,2,3)=',l1e(ny),l2e(ny),l3e(ny)
      Write(11,"(a15,3i8)") 'nffn=',nm

!  Match reactions
      call match_react

!  Calculate sparseness
      call sparse_check
      Close(4)
      Stop
      End

      function ifactorial(n)
!-----------------------------------------------------------------------------
!  This function computes factorials
!-----------------------------------------------------------------------------
      integer ifactorial
      ifactorial=1
      Do i=1,n
        ifactorial=ifactorial*i
      Enddo
      Return
      End

      subroutine match_react
!------------------------------------------------------------------------
!  In this data set, the forward and reverse reactions are separate.
!  For the purposes of this routine, forward reactions are defined as
!  those with positive Q values.  Additional, some reaction have multiple 
!  components.  This routine matches these multiple components of forward 
!  and reverse rates.  While this data is unnecessary for the basic 
!  functioning of the network, it is useful for some analysis.  For 
!  example, to calculate the net flux through a reaction channel, one must 
!  sum the components.
!------------------------------------------------------------------------
      use nuclear_data
      use reaction_data
      use flux_data
      integer               :: itest(7)
      integer               :: ifl1(nrm1),ifl2(nrm2),ifl3(nrm3)
      real(8)               :: q,rvar
      character(5)          :: blank,arrow
      character(4)          :: desc
      mflx=0

!  Match the one reactant reactions
      Do j=1,3
        nprod=nnucr(j)-nreactant(j)
        nreact=nreactant(j)
        nnuc=nnucr(j)
        Do k=la(j),le(j)
          q=q1(k)
          iw=iwk1(k)
          desc=desc1(k)
          itest=0
          If(q>0.0) Then
            itest(1:nreact)=n1(1:nreact,k)
            itest(4:nprod+3)=n1(nreact+1:nnuc,k)
            call flux_search(iflx,iw,itest,q,desc,mflx)
            ifl1(k)=iflx
          Else
            itest(1:nprod)=n1(nreact+1:nnuc,k)
            itest(4:nreact+3)=n1(1:nreact,k)
            call flux_search(iflx,iw,itest,q,desc,mflx)
            ifl1(k)=-iflx
          Endif
        Enddo
      Enddo

!  Match two reactant reactions
      Do j=4,7
        nprod=nnucr(j)-nreactant(j)
        nreact=nreactant(j)
        nnuc=nnucr(j)
        Do k=la(j),le(j)
          q=q2(k)
          iw=iwk2(k)
          desc=desc2(k)
          itest=0
          If(q>0.0) Then
            itest(1:nreact)=n2(1:nreact,k)
            itest(4:nprod+3)=n2(nreact+1:nnuc,k)
            call flux_search(iflx,iw,itest,q,desc,mflx)
            ifl2(k)=iflx
          Else
            itest(1:nprod)=n2(nreact+1:nnuc,k)
            itest(4:nreact+3)=n2(1:nreact,k)
            call flux_search(iflx,iw,itest,q,desc,mflx)
            ifl2(k)=-iflx
          Endif
        Enddo
      Enddo

!  Match three reactant reactions
      nprod=nnucr(8)-nreactant(8)
      nreact=nreactant(8)
      nnuc=nnucr(8)
      Do k=la(8),le(8)
        q=q3(k)
        iw=iwk3(k)
        desc=desc3(k)
        itest=0
        If(q>0.0) Then
          itest(1:nreact)=n3(1:nreact,k)
          itest(4:nprod+3)=n3(nreact+1:nnuc,k)
          call flux_search(iflx,iw,itest,q,desc,mflx)
          ifl3(k)=iflx
        Else
          itest(1:nprod)=n3(nreact+1:nnuc,k)
          itest(4:nreact+3)=n3(1:nreact,k)
          call flux_search(iflx,iw,itest,q,desc,mflx)
          ifl3(k)=-iflx
        Endif
      Enddo

!  Output reaction matching
      Open(10,file='match_data',FORM='unformatted')
      Write(10) mflx,le(3),le(7),le(8)
      Write(10) ifl1(1:le(3)),ifl2(1:le(7)),ifl3(1:le(8))
      Write(10) nflx(:,1:mflx),qflx(1:mflx),iwflx(1:mflx),descx(1:mflx)

!  Output ASCII match info
      arrow=' --> '
      rvar=0.0
      Open(12,file='match_read')
      Write(12,"(9a5,1es10.3)") (nname(nflx(1:3,i)),arrow,
     &      nname(nflx(4:7,i)),descx(i),rvar,i=1,mflx)

!  Test output
      blank='     '
      Do ii=1,mflx
        Write(14,"(a2)") '--'
        Do k=1,le(3)
          If(abs(ifl1(k))==ii) write(14,"(9a5,i1,i5,1es12.4)") 
     &        nname(n1(1,k)),blank,blank,arrow,nname(n1(2:4,k)),blank,
     &        desc1(k),ires1(k),ifl1(k),q1(k)
        Enddo
        Do k=1,le(7)
          If(abs(ifl2(k))==ii) write(14,"(9a5,i1,i5,1es12.4)") 
     &      nname(n2(1:2,k)),blank,arrow,nname(n2(3:6,k)),
     &      desc2(k),ires2(k),ifl2(k),q2(k)
        Enddo
        Do k=1,le(8)
          If(abs(ifl3(k))==ii) write(14,"(9a5,i1,i5,1es12.4)") 
     &      nname(n3(1:3,k)),arrow,nname(n3(4:5,k)),blank,blank,
     &      desc3(k),ires3(k),ifl3(k),q3(k)
        Enddo
! Write flux output formatted for Matlab
        ifl_orig=nflx(count(nflx(1:3,ii)/=0),ii)
        ifl_term=nflx(count(nflx(4:7,ii)/=0)+3,ii)
        Write(12,'(i5,4f6.1,es13.5,a5)') 
     &    ii,zz(ifl_orig),nn(ifl_orig),zz(ifl_term),
     &    nn(ifl_term),1.0,descx(ii)
      Enddo

      Return
      End

      subroutine flux_search(iflx,iw,i7,q,desc,mflx)
!------------------------------------------------------------------------
!  This routine searches for a given reaction among the matched reaction
!  pairs (fluxes), returning the flux index
!------------------------------------------------------------------------
      use reaction_data
      use flux_data
      integer            :: i7(7)
      real(8)             :: q
      character(4)       :: desc

!  Check reaction against existing fluxes
      iflx=0
!     Write(14,"(a6,8i5)") 'Find',iflx,i7
      If(iw==0.and.mflx/=0) Then
!     If(mflx/=0)
        Do m=1,mflx
          If(All(i7==nflx(:,m))) Then
            iflx=m
            Exit
          Endif
        Enddo
      Endif

!  If flux does not match existing flux, create new flux
      If(iflx==0) Then
        mflx=mflx+1
        nflx(:,mflx)=i7
        qflx(mflx)=abs(q)
        descx(mflx)=desc
        iwflx(mflx)=iw
!       Write(14,"(a6,8i5)") 'NFLX+',mflx,nflx(:,mflx)
        iflx=mflx
      Else
!       Write(14,"(a6,8i5)") 'Found',iflx,nflx(:,iflx)
      Endif
      Return
      End

      subroutine sparse_check
!------------------------------------------------------------------------
!  This routine maps the network Jacobian and calculates parameters 
!  for the sparseness.  The network Jacobian is best described as a doubly 
!  bordered band diagonal, though even this structure is sparse.  To examine 
!  further, the file matr_shape lists the coordinates of the non-zero 
!  elements in a form suitable for plotting.
!------------------------------------------------------------------------
      use nuclear_data
      use reaction_data
      integer, dimension(:),   allocatable :: pb,pe
      integer, dimension(:),   allocatable :: col,row
      integer,dimension(:),allocatable :: ns11,ns21,ns22
      integer,dimension(:),allocatable :: ns31,ns32,ns33 
      real(8),  dimension(ny)    :: um,llum
      real(8),  dimension(ny,ny) :: am
      integer :: vl,ll,i0,i1
      
      ll=ny*ny
      
       
      allocate(ns11(l1e(ny)),ns21(l2e(ny)))
      allocate(ns22(l2e(ny)),ns31(l3e(ny)))
      allocate(ns32(l3e(ny)),ns33(l3e(ny)))
      
      allocate(pb(ny+1),pe(ny))
      allocate(col(ll),row(ll))
      
      col = 0
      row = 0
      pb = 0
      pe = 0
      vl = 0
      am = 0.0
      
!  Allocate Matrix
      Open(9,file='matr_shape')
      Open(22,file='sparse_ind',form='unformatted')

!  Output the matrix dimension and trial right hand side
      Write(9,"(a3,i4)") 'NY=',ny
      um=1.0 
      Write(9,"(8es10.3)") um

!  Build the test matrix
      Do i0=1,ny
        um=0.0                                                                
        um(i0)=10.0                                                   
        Do j1=l1a(i0),l1e(i0)
          l1=n1(1,mu1(j1))
          um(l1)=um(l1)+1.0
        Enddo
        Do j1=l2a(i0),l2e(i0)
            l1=n2(1,mu2(j1))
            l2=n2(2,mu2(j1))
            um(l1)=um(l1)+1.0                                            
            um(l2)=um(l2)+1.0                                            
        Enddo
        Do j1=l3a(i0),l3e(i0)
            l1=n3(1,mu3(j1))
            l2=n3(2,mu3(j1))
            l3=n3(3,mu3(j1))
            um(l1)=um(l1)+1.0                                            
            um(l2)=um(l2)+1.0                                            
            um(l3)=um(l3)+1.0                                            
        Enddo
        
        write (14,*) i0
        am(i0,:)= um

!  Check for non-zero elements
        nlum=0
        Do i1=1,ny
          If(um(i1)/=0.0) Then
            nlum=nlum+1
            llum(nlum)=i1
            write(9,"(i4,2x,i4,es11.3)") i0,i1,um(i1)
          Endif
        Enddo
        Write(14,"(18i4)") (llum(i),i=1,nlum)
      Enddo
      
      do i0=1,ny
        pb(i0) = vl + 1
        do i1=1,ny
          If(am(i0,i1)/=0.0) Then
            vl = vl + 1
            col(vl) = i1
            row(vl) = i0
          Endif
        enddo
        pe(i0) = vl
      enddo
      
      pb(ny+1) = vl + 1
      
      i0 = 0
      i1 = 0
      do i0=1,ny
          do j1=l1a(i0),l1e(i0)
            l1=n1(1,mu1(j1))
            do i1=1,vl
              if(col(i1)==l1.and.row(i1)==i0)ns11(j1)=i1
            enddo
          enddo
          do j1=l2a(i0),l2e(i0)
            l1=n2(1,mu2(j1))
            l2=n2(2,mu2(j1))
            do i1=1,vl
              if(col(i1)==l1.and.row(i1)==i0)ns21(j1)=i1
              if(col(i1)==l2.and.row(i1)==i0)ns22(j1)=i1
            enddo
          enddo
          do j1=l3a(i0),l3e(i0)
            l1=n3(1,mu3(j1))
            l2=n3(2,mu3(j1))
            l3=n3(3,mu3(j1))
            do i1=1,vl
              if(col(i1)==l1.and.row(i1)==i0)ns31(j1)=i1
              if(col(i1)==l2.and.row(i1)==i0)ns32(j1)=i1
              if(col(i1)==l3.and.row(i1)==i0)ns33(j1)=i1
            enddo
          enddo
      enddo
      
      write(22)vl
      write(22)row(1:vl),col(1:vl),pb
      write(22)l1e(ny),l2e(ny),l3e(ny)
      write(22)ns11,ns21,ns22
      write(22)ns31
      write(22)ns32
      write(22)ns33

      close(22)
      
!  Find heaviest incident nucleus, which sets border width
      Do i=1,ny
        If(nname(i)=='    n'.or.nname(i)=='    p'.or. 
     &      nname(i)=='    d'.or.nname(i)=='    t'.or. 
     &      nname(i)=='  he3'.or.nname(i)=='  he4'.or. 
     &      nname(i)=='  c12'.or.nname(i)=='  n14'.or. 
     &      nname(i)=='  o16'.or.nname(i)==' ne20') leftm=i
      Enddo

!  Find band widths
      imax=0
      jmax=0
      imin=0
      jmin=0
      ijd=leftm+1

!  From the diagonal, for each row beyond the border
      Do i=ijd,ny

!  Search rightmost (lowermost) elements of the band diagonal
        Do j=1,ny-i
          If(am(i,i+j)/=0.0.and.j>jmax) Then
            jmax=j
            jwhere=i
          Endif
          If(am(i+j,i)/=0.0.and.j>imax) Then
            imax=j
            iwhere=i
          Endif
        Enddo

!  Search for the leftmost (uppermost) element of the band diagonal
        Do j=1,i-1
          If(i-j>leftm.and.am(i,i-j)/=0.0.and.j>jmin) Then
            jmin=j
            If(i-j.le.leftm) jmin=i-leftm-1
          Endif
          If(i-j>leftm.and.am(i-j,i)/=0.0.and.j>imin) Then
            imin=j
            If(i-j.le.leftm) imin=i-leftm-1
          Endif
        Enddo
      Enddo

!  Output sparseness parameters
      write(4) leftm,leftm,ijd,imin,imax
      Write(11,"('Matrix Sparseness parameters')")
      Write(11,"('Border Widths   ',2i5)") leftm,leftm
      Write(11,"('Diagonal Widths ',2i5,'(j<i,j>i)')") jmin,jmax
      
      deallocate(ns11,ns21,ns22)
      deallocate(ns31,ns32,ns33)
      
      deallocate(pb,pe)
      deallocate(col,row)
      Return
  210 Format(3(10x,i8))
      End

      subroutine format_nuclear_data
!-----------------------------------------------------------------------------  
!  This routine reads from the file netwinv, the nuclei included along 
!  with the nuclear data which will be needed for later calculations.  This 
!  data includes the atomic number, the number of protons and neutrons, and 
!  the binding energy (calculated from the tabulated mass excess).  Also the 
!  tabulations of the  partition functions, g, are read in for later 
!  interpolation.  A binary data file is then prepared for faster input.
!-----------------------------------------------------------------------------  
      use nuclear_data
      use part_funct_data
      integer, parameter :: nymx=10000
      real(8)  me                                      
      integer it9i(24)
      character(5) nam,ntest(nymx)                                 

!  Read in sunet
      Open(7,FILE='sunet',STATUS='old')                      
      Do i=1,nymx
        Read(7,"(2a5)",END=100) ntest(i)
      Enddo  
  100 ny=i-1
      Close(7)
      Write(6,*) ny

!  Read nuclear data from netwinv
      Open(8,FILE='netwinv',STATUS='old')                      
      Read(8,"(a5)") nam   

!  Read in the partition function iteration grid, and fix endpoints
      Read(8,"(24i3)") (it9i(i),i=1,24)                   
      Do i=1,24                                                    
        t9i(i)=it9i(i)*0.01     
      Enddo                                             
      t9i(24)=t9i(24)*10.                                                       

!  Allocate the nuclear name array.  nname(0) exists since blanks in the
!  input file get indexed to 0.
      Allocate (nname(0:ny))
      nname(0)='     '

!  Read the nuclear names from netwinv and check against sunet
      Do n=1,ny
        Read(8,"(a5)") nname(n)
        If(nname(n)/=ntest(n)) Then
          Write(6,*) 'Netwinv /= sunet',n,nname(n),ntest(n)
          Stop
        Endif
      Enddo                     

!  Set size of nuclear parameter arrays and read in nuclear parameters 
!  and partition function interpolation table. 
      Allocate (aa(ny),zz(ny),nn(ny),be(ny),g(24,ny),angm(ny))
      Do l=1,ny                                                  
        Read(8,"(a5,f12.3,2i4,f6.1,f10.3)") nam,a,na,nb,sp,me  
        Read(8,"(8f9.2)") (g(m,l),m=1,24)     
        aa(l)=a                                                         
        zz(l)=dble(float(na))                                                
        nn(l)=dble(float(nb))                                                
        be(l)=8.07144*nn(l)+7.28899*zz(l)-me 
        angm(l)=2.*sp+1.
        nname(l)=nam
      Enddo                                                           
      Close(8)
!     Write(g,"(a5,4es12.4)") (nname(i),aa(i),zz(i),nn(i),be(i),i=1,ny)        

!  Write binary data file
      Open(8,FILE='nuc_data',FORM='unformatted')
      Write(8) ny
      Write(8) t9i
      Write(8) (nname(i),aa(i),zz(i),nn(i),be(i),g(:,i),angm(i),i=1,ny)
      Close(8)
      Return                                                                    
      End                                                                       
