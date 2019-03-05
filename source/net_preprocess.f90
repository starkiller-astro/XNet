!*******************************************************************************
! net_preprocess.f90, 2/21/12
! Reaction Network Data Preprocessing
!
! These routines take the human readable, ASCII data files sunet, netsu, 
! netweak, netneutrino and netwinv and prepares binary versions for the reaction network. 
! This involves translating Characters to indicies and organizing the data  
! for more efficient computation.  Several additional tests are performed.
!*******************************************************************************

Module process_reaction_data
!===============================================================================
! This module contains the reaction rate data, including the indices
! which type the reaction and the cross-references.
!===============================================================================
  Use nuc_number
  Integer, parameter :: nrm1=50000,nrm2=60000,nrm3=100
  Integer            :: nnucr(8)=(/2,3,4,3,4,5,6,5/)
  Integer            :: nreactant(8)=(/1,1,1,2,2,2,2,3/)
  Integer            :: la(8),le(8)
  Real(8)            ::    q1(nrm1),   q2(nrm2),   q3(nrm3)
  Integer            ::  iwk1(nrm1), iwk2(nrm2), iwk3(nrm3)
  Integer            :: ires1(nrm1),ires2(nrm2),ires3(nrm3)
  Integer            :: irev1(nrm1),irev2(nrm2),irev3(nrm3)
  Character(LEN=4)   :: desc1(nrm1),desc2(nrm2),desc3(nrm3)
  Integer            ::  n1(4,nrm1), n2(6,nrm2), n3(5,nrm3)
  Real(8)            :: an1(4,nrm1),an2(6,nrm2),an3(5,nrm3)
  Integer            :: mu1(4*nrm1),mu2(6*nrm2),mu3(5*nrm3)
  Real(8)            ::  a1(4*nrm1), a2(6*nrm2),a3(5*nrm3)
  Integer, dimension(:), allocatable :: l1a,l2a,l3a
  Integer, dimension(:), allocatable :: l1e,l2e,l3e
End module process_reaction_data

Module process_flux_data
!===============================================================================
! This module contains the 
!===============================================================================
  Integer,parameter :: mxflx=60000
  Integer           :: nflx(7,mxflx),iwflx(mxflx)
  Real(8)           :: qflx(mxflx)
  Character(LEN=4)  :: descx(mxflx)
End module process_flux_data

Function ifactorial(n)
!===============================================================================
! This function computes factorials
!===============================================================================
  Integer, Intent(in) :: n
  Integer             :: ifactorial,i
  ifactorial=1
  Do i=1,n
    ifactorial=ifactorial*i
  Enddo
  Return
End Function ifactorial

Subroutine net_preprocess(lun_out,data_dir,data_desc)
!===============================================================================
! This is the outer shell of Net Setup.  This routine also reformats 
! the reaction data.
!===============================================================================
  Use nuclear_data
  Use process_reaction_data
  Integer, Intent(in)            :: lun_out   ! The logical unit for progress output
  Character(LEN=80), Intent(in)  :: data_dir  ! The location of the data directory
  Character(LEN=80), Intent(in)  :: data_desc ! Description of data directory
!-------------------------------------------------------------------            
  Integer, parameter :: ndesc=152
  Integer          :: n,i,j,k,l,jj ! Loop indicies
  Integer          :: ni(6),ki(6),ierr
  Integer          :: nffn,nnnu,krt,kl,lo,ifactorial
  Integer          :: iec,ier,iev  ! Reaction flags
  Integer          :: idesc,nnew_desc ! Reaction descriptor indices
  Character(LEN=5) :: name(6),blank5='     '
  Character(LEN=4) :: desc,desc_new(100)
  Character(LEN=1) :: nr,vw,nrr='r',nrn='n',nrw='w',vvw='v',blank=' '
  Character(LEN=4), dimension(ndesc) :: desc_known = &
  & (/'  ec','bet+','betp','bet-','betm','btyk',' bec','bkmo', &
  &   'mo92',' ffn',' lmp','fkth','cf88','wies','laur','bb92', &
  &   'baka','rolf','wawo','mafo','ms01',' wag','bb93','ca85', &
  &   'ornl',' wfh','rath','rpsm','hg94','ro95','wien','iucf', &
  &   'smwa','smmo','wi2p','aexp','pexp','bex+','bex-','bqa+', &
  &   'bhi+','bec ','bklq','nuba','bede','bqa+','gamm','ja01', &
  &   'nosm','nacr','tang','ta04','ch04','ku02','il04','go00', &
  &   'most','il03','de04','da03','vanc','il01','ha00','ku92', &
  &   'ro04','sc83','he98','po00','il00','he95','sc98','vo00', &
  &   'il99','beau','da04','ca00','da77','bl03','ch01','fynb', &
  &   'De99','bu96','ba00','Ha96','sh03','tr04','sc05','nfis', &
  &   'GW95','ch05','re98','nb03','mo03','wc07','an06','ka02', &
  &   'gl07','GG90','lg06','im05','GG92','ct07','ga00','ww02', &
  &   'wh87','dh03','ha01','cb01','HF99','hg95','pg06','jm06', &
  &   'sb05','SM93','SM86','thra','fy05','ha04','pt05','bk91', &
  &   'bk92','kg91','SM91','iwam',' nen','nebp','lznu','lzan', &
  &   'ths8','wc12','nk06','mb11','jz10','cd08','nw00','il10', &
  &   'li10','chw0','co10','kd02','st08','pkrF','br10','hi09', &
  &   'ld11','hd08','ua08','ol11','ls09','ia08','dc11','mb07'/) 
  Real(8) p(7),q,qtemp,r
  
  nnew_desc = 1
  desc_new(1) = '    '
  
  Open(2,file=trim(data_dir)//'/netsu')
  Open(3,file=trim(data_dir)//'/nets3',form='unformatted')
  Open(4,file=trim(data_dir)//'/nets4',form='unformatted')
  Open(14,file=trim(data_dir)//'/net_diag')

! Read and Reformat the nuclear data
  Call format_nuclear_data(lun_out,data_dir)
  Write(lun_out,"('Listing Nuclei')") 
  Write(lun_out,"(10a6)") (nname(i),i=1,ny)
  Write(lun_out,"(1x,'Number of species=',i4)") ny

! Read and Reformat the reaction data
  Write(lun_out,"('Reading and reformating reactions data')") 
  Read(2,"(2i5)") nffn,nnnu
  n=1
  Do jj=1,100000

! Read in each entry from the reaction data file
    Read(2,"(i1,4x,6a5,8x,a4,a1,a1,3x,1pe12.5)",Iostat=ierr) k,(name(j),j=1,6),desc,nr,vw,q
    If(ierr < 0) Then
      Exit
    ElseIf(ierr > 0) Then
      Write(lun_out,*) 'Problem reading reaction ',j
      Stop
    Endif
    Read(2,"(4e13.6)",Iostat=ierr) (p(j),j=1,7)
    If(ierr < 0) Then
      Exit
    ElseIf(ierr > 0) Then
      Write(lun_out,*) 'Problem reading reaction ',j
      Stop
    Endif

! If the entry is a reaction
    If(k==0) Then

! Convert nuclear names to indicies
      Do i=1,6
        Call index_from_name(name(i),ni(i))
      Enddo
      
      If(any(ni>ny)) Then
        Write(lun_out,*)'Dropping reaction ',name,'Some nuclei not included in netsu.'
        Cycle
      Endif
      
! Identify the reaction source against the list of known sources.
      Do i=1,ndesc
        If(desc==desc_known(i)) Exit
      EndDo
      If(i>ndesc) Then
! The reaction source is unknown, treat it as a normal strong rate.
        idesc=i
! Check if this is the first time this new reaction source has showed up in this particular netsu.
        Do i=1,nnew_desc
          If(desc==desc_new(i)) Exit
        EndDo
! If this is the first time, remember the name so the User can be informed of unknown codes at the end.
        If (i>nnew_desc) Then
          nnew_desc = nnew_desc + 1
          desc_new(i) = desc
          Write(lun_out,'(3a,i6,a)') 'Reaction Desc. ',desc, &
          & ' not found for entry ',jj,', treating as strong rate.'
        EndIf
      Else
        idesc=i
      Endif           
!-------------------------------------------------------------------------------
! The reaction rates are calculated via a number of templates, depending
! on the type of reaction and its source.  In this section, flags are 
! set for correct treatment of these different types of reactions
!-------------------------------------------------------------------------------
! For weak interactions; ,  (iec=1), Fuller, Fowler, Neuman Rates (iec=2,3),
! Other weak rates not covered above (iec=4)
      Select Case (idesc)
        Case Default 
          iec=0               ! not a weak rate 
        Case (1)              ! electron capture rates from Fowler and Co.
          iec=1
        Case (10,11)          ! FFN type tabulated rate
          If(zz(ni(1))>zz(ni(2)).and.zz(ni(1))/=1) Then
            iec=3
          Else
            iec=2
          Endif
        Case(2:9,38:45,92:94,130) ! beta decay and non-tabulated ec
          iec=4
        Case(125,127)         ! electron neutrino capture 
          iec=7                                  
        Case(126,128)         ! electron antineutrino capture 
          iec=8                                  
      End Select 

! Check resonant nature of rate, resonant rate (ier=2), nonresonant (ier=1),
! or none of the above (ier=0).  JINA REALCIB also includes a weak flag.
      If(nr==blank) Then
        ier=0
      ElseIf(nr==nrw) Then ! weak reaction
        ier=0
        If(iec==0) Then
          Write(lun_out,'(5a,i2)') 'Reaction with description ',desc,' has resonant flag ',nr,' but iec=',iec
        Endif
      Elseif(nr==nrn) Then
        ier=1
      Elseif(nr==nrr) Then
        ier=2
      Else
        ier=2
      Endif

! Check for inverse rate (iev=1). Inverse rates must be multiplied with 
! Temperature dependant ratio of partition functions. iev=0 is a forward rate.
      If(vw==vvw) Then
        iev=1
      Else
        iev=0
      Endif

! Divide nuclei name list into reactants (ki=-1), products (ki=+1) and
! blanks (ki=0)
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
      
! Don't trust the Q-values in the netsu file
      qtemp = 0.0
      Do i=1,6
        If(i>nreactant(krt)) Then
          If(name(i)==blank5) Then
            qtemp = qtemp
          Else
            qtemp = qtemp + be(ni(i))-8.07144*nn(ni(i))-7.28899*zz(ni(i))
          Endif
        Else
          qtemp = qtemp - be(ni(i))+8.07144*nn(ni(i))+7.28899*zz(ni(i))
        Endif
      Enddo
      If (abs(qtemp-q)>abs(.1*q)) Then
        Write(14,'(a,7a6,a,es9.2,a,es9.2)') &
        & 'Inconsistent q-value for ',desc,name,'  netsu: ',q,' netwinv: ',qtemp
      Endif
      q=qtemp
      
! Count identical reactants
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

! Count identical products
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

! Calculate double counting corrections and output reaction data to nets3
      If(nreactant(krt)==1) Then
        r=1.0
        Do i=1,nnucr(krt)
          n1(i,n)=ni(i)
          an1(i,n)=real(ki(i))*r
        Enddo
        q1(n)=q
        iwk1(n)=iec
        ires1(n)=ier
        irev1(n)=iev
        desc1(n)=desc
        Write(3) n,(ni(i),i=1,4),iec,ier,iev,p,q
!       Write(21,*) n,(nname(ni(i)),i=1,4),iec,ier,iev,p,q
      Else
        lo=0
        Do i=1,nreactant(krt)
          If(ki(i)<lo) lo=ki(i)
        Enddo
        r=1.0/float(ifactorial(-lo))
        If(nreactant(krt)==2) Then
          Do i=1,nnucr(krt)
            n2(i,n)=ni(i)
            an2(i,n)=real(ki(i))*r
          Enddo
          q2(n)=q
          iwk2(n)=iec
          ires2(n)=ier
          irev2(n)=iev
          desc2(n)=desc
          Write(3) n,(ni(i),i=1,5),iec,ier,iev,p,q
!         Write(21,*) n,(ni(i),i=1,5),iec,ier,iev,p,q
        Elseif(nreactant(krt)==3) Then
          Do i=1,nnucr(krt)
            n3(i,n)=ni(i)
            an3(i,n)=real(ki(i))*r
          Enddo
          q3(n)=q
          iwk3(n)=iec
          ires3(n)=ier
          irev3(n)=iev
          desc3(n)=desc
          Write(3) n,(ni(i),i=1,6),iec,ier,iev,p,q
!         Write(21,*) n,(ni(i),i=1,6),iec,ier,iev,p,q
        Endif
      Endif
      n=n+1

! If data entry is not a reaction, it is a type or sub-type marker
! For type marker, reset type and sub-type counters
    Elseif(k==1) Then
      krt=k
      la(k)=1
      n=1
    Elseif(k==4.or.k==8) Then
      krt=k
      la(k)=1
      le(k-1)=n-1
      n=1

! For sub-type marker, restart sub-type counters
    Else
      krt=k
      la(k)=n
      le(k-1)=n-1
    Endif
!   Write(14,"(1x,i2,6i4,6i2,3i2,2e10.4/1x,7e10.4)") k,ni(:),ki(:),iec,ier,iev,r,q,p(:)
  Enddo
  le(8)=n-1
  Close(1)
  Close(2)

!-------------------------------------------------------------------------------
! In order to build the Jacobian efficiently (row or column wise rather
! than scattered), information about which reactions affect each nucleus 
! is necessary.  
!-------------------------------------------------------------------------------
  Write(lun_out,"('Registering reactions to species')") 
! For each nucleus, register each single reactant reaction which involves it.
  Allocate(l1a(ny),l1e(ny),l2a(ny),l2e(ny),l3a(ny),l3e(ny))
  n=0
  Do i=1,ny
    l1a(i)=n+1
    Do j=1,3 ; Do k=la(j),le(j) ;  Do l=1,nnucr(j)
      If(n1(l,k)==i.and.an1(l,k)/=0.) Then
        n=n+1
        a1(n)=an1(l,k)
        mu1(n)=k
        Write(3) a1(n),mu1(n)
!       Write(21,*) a1(n),mu1(n)
      Endif
    Enddo ; Enddo ; Enddo
    l1e(i)=n
  Enddo

! For each nucleus, register each two reactant reaction which involves it.
  n=0
  Do i=1,ny
    l2a(i)=n+1
    Do j=4,7 ;  Do k=la(j),le(j) ; Do l=1,nnucr(j)
      If(n2(l,k)==i.and.an2(l,k)/=0.) Then
        n=n+1
        a2(n)=an2(l,k)
        mu2(n)=k
        Write(3) a2(n),mu2(n)
!       Write(21,*) a2(n),mu2(n)
      Endif
    Enddo ; Enddo ; Enddo
    l2e(i)=n
  Enddo

! For each nucleus, register each three reactant reaction which involves it.
  n=0
  Do i=1,ny
    l3a(i)=n+1
    Do k=la(8),le(8) ; Do l=1,nnucr(8)
      If(n3(l,k)==i.and.an3(l,k)/=0.) Then
        n=n+1
        a3(n)=an3(l,k)
        mu3(n)=k
        Write(3) a3(n),mu3(n)
!       Write(21,*) a3(n),mu3(n)
      Endif
    Enddo ; Enddo
    l3e(i)=n
  Enddo
  Close(3)

! Output the nucleus to reaction registration to nets4
  Write(4) ny
  Write(4) (nname(i),i=1,ny)
  Write(4) nffn,nnnu
  Write(4) le(3),le(7),le(8)
  Do i=1,ny
    Write(4) i,l1a(i),l1e(i),l2a(i),l2e(i),l3a(i),l3e(i)
  Enddo
  Write(lun_out,"(10x,'number of reactions of different types')")
  Write(lun_out,"(3(10x,i8))") (j,la(j),le(j),j=1,8)
  Write(lun_out,"(10x,'necessary dimensions')")
  Write(lun_out,"(3(10x,i8))") le(3),le(7),le(8),l1e(ny),l2e(ny),l3e(ny)
  Write(lun_out,"(3(10x,i8))") nffn, nnnu

! Prepare network description file
  Open(11,file=trim(data_dir)//'/net_desc')
  Write(11,"(a80)") data_desc
  Write(11,"('Number of Nuclear Species=',i5)") ny
  Write(11,"('Reaction Count for the 8 different types')")
  Write(11,"(3i10)") (j,la(j),le(j),j=1,8)
  Write(11,"('Necessary dimensions')")
  Write(11,"(a15,3i8)") 'nreac(1,2,3)=',le(3),le(7),le(8)
  Write(11,"(a15,3i8)") 'nan(1,2,3)=',l1e(ny),l2e(ny),l3e(ny)
  Write(11,"('Reaction Count for non-REACLIB rates')")
  Write(11,"(a7,i8,a7,i8)") 'nffn=',nffn,' nnu=',nnnu
  
! Match reactions
  Write(lun_out,*) "Matching forward and reverse reactions"
  Call match_react(data_dir)

! Calculate sparseness
  Write(lun_out,*) "Determining network sparsity pattern"
  Call sparse_check(data_dir)
  Close(4)

! Create inab template
  Open(13,file=trim(data_dir)//'/ab_blank')
  Write(13,"(a)") 'Abundance Description'
  Write(13,"(4(a5,es14.7,1x))") (nname(i),0.0,i=1,ny)
  Close(13)
  
! Tell User to add new descriptions to net_setup
  If(nnew_desc>1) Then
    Write(lun_out,*) "These descriptions need to be added to net_setup:"
    Write(lun_out,"(8(a4))")(desc_new(i),i=2,nnew_desc)
  Endif

! Deallocate Nuclear Data arrays
  Deallocate(nname,aa,zz,nn,be)
  
  Close(11)
  Close(14)
  Return
End Subroutine net_preprocess

Subroutine match_react(data_dir)
!===============================================================================
! In this data set, the forward and reverse reactions are separate.
! For the purposes of this routine, forward reactions are defined as
! those with positive Q values.  Additional, some reaction have multiple 
! components.  This routine matches these multiple components of forward 
! and reverse rates.  While this data is unnecessary for the basic 
! functioning of the network, it is Useful for some analysis.  For 
! example, to calculate the net flux through a reaction channel, one must 
! sum the components.
!===============================================================================
  Use nuclear_data
  Use process_reaction_data
  Use process_flux_data
  Character (LEN=*),  INTENT(in)  :: data_dir
  Integer          :: mflx ! the number of reaction sets
  Integer          :: iflx ! the reaction set index
  Integer          :: i,j,k,ii ! Loop variables
  Integer          :: itest(7),nprod,nreact,nnuc,iw
  Integer          :: ifl1(nrm1),ifl2(nrm2),ifl3(nrm3)
  Integer          :: revtest,comtest
  Integer          :: ifl_orig,ifl_term ! indicies of originating and termination species
  Real(8)          :: q,rvar
  Character(LEN=5) :: blank,arrow
  Character(LEN=4) :: desc
  mflx=0

! Match the one reactant reactions
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
        Call flux_search(iflx,iw,itest,q,desc,mflx)
        ifl1(k)=iflx
      Else
        itest(1:nprod)=n1(nreact+1:nnuc,k)
        itest(4:nreact+3)=n1(1:nreact,k)
        Call flux_search(iflx,iw,itest,q,desc,mflx)
        ifl1(k)=-iflx
      Endif
    Enddo
  Enddo

! Match two reactant reactions
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
        Call flux_search(iflx,iw,itest,q,desc,mflx)
        ifl2(k)=iflx
      Else
        itest(1:nprod)=n2(nreact+1:nnuc,k)
        itest(4:nreact+3)=n2(1:nreact,k)
        Call flux_search(iflx,iw,itest,q,desc,mflx)
        ifl2(k)=-iflx
      Endif
    Enddo
  Enddo

! Match three reactant reactions
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
      Call flux_search(iflx,iw,itest,q,desc,mflx)
      ifl3(k)=iflx
    Else
      itest(1:nprod)=n3(nreact+1:nnuc,k)
      itest(4:nreact+3)=n3(1:nreact,k)
      Call flux_search(iflx,iw,itest,q,desc,mflx)
      ifl3(k)=-iflx
    Endif
  Enddo

! Output reaction matching
  Open(10,file=trim(data_dir)//'/match_data',FORM='unformatted')
  Write(10) mflx,le(3),le(7),le(8)
  Write(10) ifl1(1:le(3)),ifl2(1:le(7)),ifl3(1:le(8))
  Write(10) nflx(:,1:mflx),qflx(1:mflx),iwflx(1:mflx),descx(1:mflx)
  Close(10)

! Output ASCII match info
  arrow=' --> '
  rvar=0.0
  Open(12,file=trim(data_dir)//'/match_read')
  Write(12,"(9a5,1es10.3)") (nname(nflx(1:3,i)),arrow,nname(nflx(4:7,i)),descx(i),rvar,i=1,mflx)

! Test for matching components, reverse reactions and output
  blank='     '
  Do ii=1,mflx
    revtest=0;comtest=0
    Write(14,"(a2)") '--'
    Do k=1,le(3)
      If(abs(ifl1(k))==ii) Then 
        Write(14,"(9a5,2i2,i6,1es12.4)") & 
        & nname(n1(1,k)),blank,blank,arrow,nname(n1(2:4,k)),blank,desc1(k),ires1(k),irev1(k),ifl1(k),q1(k)
        comtest=comtest+1
        revtest=revtest+irev1(k)
      Endif
    Enddo
    Do k=1,le(7)
      If(abs(ifl2(k))==ii) Then 
        Write(14,"(9a5,2i2,i6,1es12.4)") &
        & nname(n2(1:2,k)),blank,arrow,nname(n2(3:6,k)),desc2(k),ires2(k),irev2(k),ifl2(k),q2(k)
        comtest=comtest+1
        revtest=revtest+irev2(k)
      Endif
    Enddo
    Do k=1,le(8)
      If(abs(ifl3(k))==ii) Then
        Write(14,"(9a5,2i2,i6,1es12.4)") &
        & nname(n3(1:3,k)),arrow,nname(n3(4:5,k)),blank,blank,desc3(k),ires3(k),irev3(k),ifl3(k),q3(k)
        comtest=comtest+1
        revtest=revtest+irev3(k)
      Endif
    Enddo

! Write test results
    If(mod(comtest,2)/=0.and.iwflx(ii)==0) Write(14,"(a,i4)") 'Unbalanced reaction',comtest
    If(comtest>1.and.revtest==0) Write(14,"(a,i4)") 'No reverse defined',revtest   
    
! Write flux output formatted for Matlab
    ifl_orig=nflx(count(nflx(1:3,ii)/=0),ii)
    ifl_term=nflx(count(nflx(4:7,ii)/=0)+3,ii)
    Write(12,'(i5,4f6.1,es13.5,a5)') &
    & ii,zz(ifl_orig),nn(ifl_orig),zz(ifl_term),nn(ifl_term),1.0,descx(ii)
  Enddo
  Close(12)

  Return
End Subroutine match_react

Subroutine flux_search(iflx,iw,i7,q,desc,mflx)
!===============================================================================
! This routine searches for a given reaction among the matched reaction
! pairs (fluxes), returning the flux index
!===============================================================================
  Use process_reaction_data
  Use process_flux_data
  Integer, Intent(out)         :: iflx  ! Search result
  Integer, Intent(in)          :: iw    ! weak reaction flag
  Integer, Intent(in)          :: i7(7) ! nuclear indicies to be matched
  Real(8), Intent(in)          :: q     ! Reaction Q value
  Character(LEN=4), Intent(in) :: desc  ! Reaction source descriptor
  Integer, Intent(inout)       :: mflx  ! Number of reaction pairs
  Integer                      :: m     ! loop variable

! Check reaction against existing fluxes
  iflx=0
! Write(14,"(a6,8i5)") 'Find',iflx,i7
  If(iw==0.and.mflx/=0) Then
! If(mflx/=0)
    Do m=1,mflx
      If(All(i7==nflx(:,m))) Then
        iflx=m
        Exit
      Endif
    Enddo
  Endif

! If flux does not match existing flux, create new flux
  If(iflx==0) Then
    mflx=mflx+1
    nflx(:,mflx)=i7
    qflx(mflx)=abs(q)
    descx(mflx)=desc
    iwflx(mflx)=iw
!   Write(14,"(a6,8i5)") 'NFLX+',mflx,nflx(:,mflx)
    iflx=mflx
  Else
!   Write(14,"(a6,8i5)") 'Found',iflx,nflx(:,iflx)
  Endif
  Return
End Subroutine flux_search

Subroutine sparse_check(data_dir)
!===============================================================================
! This routine maps the network Jacobian and calculates parameters for the 
! sparseness.  
! The network Jacobian is reasonably described as a doubly bordered band diagonal, 
! though even this structure contains considerable zero elements.  
!===============================================================================
  Use nuclear_data
  Use process_reaction_data
  Character (LEN=*),  INTENT(in)       :: data_dir
  Real(8), dimension(ny)               :: um       ! Jacobian Row
  Integer, dimension(ny)               :: llum     ! List of non-zero row elements
  Real(8), dimension(:,:), allocatable :: am       ! Dense Jacobian
  Integer, dimension(:), allocatable   :: col,row  ! Column and Row positions of non-zero elements
  Integer, dimension(ny+1)             :: pb       ! First column index for each row in row vector
  Integer, dimension(ny)               :: pe       ! Last column index for each row in row vector
  Integer :: ns11(l1e(ny)),ns21(l2e(ny)),ns22(l2e(ny))  ! Map of sparse matrix location 
  Integer :: ns31(l3e(ny)),ns32(l3e(ny)),ns33(l3e(ny))  ! for each reaction
  Integer :: vl,ll,i0,i1
  Integer :: i,j,j1,l1,l2,l3 ! Loop Indicies
  Integer :: nlum,leftm,imin,imax,iwhere,jmin,jmax,jwhere,ijd ! sparseness temp variables
!-------------------------------------------------------------------------------
! To exploit the double bordered band form, we must indentify the band and 
! diagonal widths.
!-------------------------------------------------------------------------------

! The file matr_shape lists the coordinates of the non-zero elements in a form suitable for plotting.
! Output the matrix dimension and trial right hand side
  Open(9,file=trim(data_dir)//'/matr_shape')
  Write(9,"(a3,i4)") 'NY=',ny
  um=1.0
  Write(9,"(8es10.3)") um

! Build the test matrix
  Allocate(am(ny,ny))
  am = 0.0
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
    
    Write (14,*) i0
    am(i0,:)= um

! Check for non-zero elements
    nlum=0
    Do i1=1,ny
      If(um(i1)/=0.0) Then
        nlum=nlum+1
        llum(nlum)=i1
        Write(9,"(i4,2x,i4,es11.3)") i0,i1,um(i1)
      Endif
    Enddo
    Write(14,"(18i4)") (llum(i),i=1,nlum)
  Enddo

  Close(9)
! Find heaviest incident nucleus, which sets border width
  Do i=1,ny
    If(nname(i)=='    n'.or.nname(i)=='    p'.or. nname(i)=='    d'.or.nname(i)=='    t'.or. & 
    &  nname(i)=='  he3'.or.nname(i)=='  he4'.or. nname(i)=='  c12'.or.nname(i)=='  n14'.or. &
    &  nname(i)=='  o16'.or.nname(i)==' ne20') leftm=i
  Enddo

! Find band widths
  imax=0
  jmax=0
  imin=0
  jmin=0
  ijd=leftm+1

! From the diagonal, for each row beyond the border
  Do i=ijd,ny

! Search rightmost (lowermost) elements of the band diagonal
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

! Search for the leftmost (uppermost) element of the band diagonal
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

! Output parameters describing doubly-bordered band diagonal form.
  Write(4) leftm,leftm,ijd,imin,imax
  Write(11,"('Matrix Sparseness parameters')")
  Write(11,"('Border Widths   ',2i5)") leftm,leftm
  Write(11,"('Diagonal Widths ',2i5,'(j<i,j>i)')") jmin,jmax

!-------------------------------------------------------------------------------
! Sparse solver methods like PARDISO and MA28 use Compressed Row storage (CRS) 
! to efficiently story the matrix to be solved
!-------------------------------------------------------------------------------

! Build the Jacobian in CRS
  Allocate(col(ny*ny),row(ny*ny))
  col = 0 ; row = 0 ; pb = 0 ; pe = 0 ; vl = 0
  Do i0=1,ny
    pb(i0) = vl + 1
    Do i1=1,ny
      If(am(i0,i1)/=0.0) Then
        vl = vl + 1
        col(vl) = i1
        row(vl) = i0
      Endif
    Enddo
    pe(i0) = vl
  Enddo
  pb(ny+1) = vl + 1
  
! Build arrays to map reaction rates into locations in the CRS Jacobian 
  i0 = 0
  i1 = 0
  Do i0=1,ny
    Do j1=l1a(i0),l1e(i0)
      l1=n1(1,mu1(j1))
      Do i1=1,vl
        If(col(i1)==l1.and.row(i1)==i0) ns11(j1)=i1
      Enddo
    Enddo
    Do j1=l2a(i0),l2e(i0)
      l1=n2(1,mu2(j1))
      l2=n2(2,mu2(j1))
      Do i1=1,vl
        If(col(i1)==l1.and.row(i1)==i0) ns21(j1)=i1
        If(col(i1)==l2.and.row(i1)==i0) ns22(j1)=i1
      Enddo
    Enddo
    Do j1=l3a(i0),l3e(i0)
      l1=n3(1,mu3(j1))
      l2=n3(2,mu3(j1))
      l3=n3(3,mu3(j1))
      Do i1=1,vl
        if(col(i1)==l1.and.row(i1)==i0)ns31(j1)=i1
        if(col(i1)==l2.and.row(i1)==i0)ns32(j1)=i1
        if(col(i1)==l3.and.row(i1)==i0)ns33(j1)=i1
      Enddo
    Enddo
  Enddo
  
! Output the CRS jacobian format data
  Open(22,file=trim(data_dir)//'/sparse_ind',form='unformatted')
  Write(22)vl
  Write(22)row(1:vl),col(1:vl),pb
  Write(22)l1e(ny),l2e(ny),l3e(ny)
  Write(22)ns11,ns21,ns22
  Write(22)ns31
  Write(22)ns32
  Write(22)ns33
  Close(22)

! Recover matrix memory
  Deallocate(am,col,row)
  
  Return
End Subroutine sparse_check

Subroutine format_nuclear_data(lun_out,data_dir)
!===============================================================================
! This routine reads from the file netwinv, the nuclei included along 
! with the nuclear data which will be needed for later calculations.  This 
! data includes the atomic number, the number of protons and neutrons, and 
! the binding energy (calculated from the tabulated mass excess).  Also the 
! tabulations of the  partition functions, g, are read in for later 
! interpolation.  A binary data file is Then prepared for faster input.
!===============================================================================
  Use nuclear_data
  Use part_funct_data
  Character (LEN=*),  INTENT(in)  :: data_dir
  Integer, Intent(in) :: lun_out
  Integer, parameter :: nymx=10000
  Real(8)  :: a,sp,me  ! Temp read variables
  Integer  :: na,nb    ! Temp read variables                                    
  Integer :: it9i(24)  ! Temp read partition function grid
  Character(LEN=5) :: nam,ntest(nymx)                                 
  Integer i,n,l,m   ! Loop Indicies
  
  
! Read in sunet
  Open(7,FILE=trim(data_dir)//'/sunet',STATUS='old')                      
  Do i=1,nymx
    Read(7,"(2a5)",END=100) ntest(i)
  Enddo  
  100 ny=i-1
  Close(7)
  Write(6,*) ny

! Read nuclear data from netwinv
  Open(8,FILE=trim(data_dir)//'/netwinv',STATUS='old')                      
  Read(8,"(a5)") nam   

! Read in the partition function iteration grid, and fix endpoints
  Read(8,"(24i3)") (it9i(i),i=1,24)                   
  Do i=1,24                                                    
    t9i(i)=it9i(i)*0.01     
  Enddo                                             
  t9i(24)=t9i(24)*10.                                                       

! Allocate the nuclear name array.  nname(0) exists since blanks in the
! input file get indexed to 0.
  Allocate (nname(0:ny))
  nname(0)='     '

! Read the nuclear names from netwinv and check against sunet
  Do n=1,ny
    Read(8,"(a5)") nname(n)
    If(nname(n)/=ntest(n)) Then
      Write(lun_out,*) 'Netwinv /= sunet',n,nname(n),ntest(n)
      Stop
    Endif
  Enddo                     

! Set size of nuclear parameter arrays and read in nuclear parameters 
! and partition function interpolation table. 
  Allocate (aa(ny),zz(ny),nn(ny),be(ny),g(24,ny),angm(ny))
  Do l=1,ny  
    Write(6,*) l                                                
    Read(8,"(a5,f12.3,2i4,f6.1,f15.8)") nam,a,na,nb,sp,me  
    Read(8,*) (g(m,l),m=1,24)     
    aa(l)=a                                                         
    zz(l)=float(na)                                                
    nn(l)=float(nb)                                                
    be(l)=8.07144*nn(l)+7.28899*zz(l)-me 
    angm(l)=2.*sp+1.
    nname(l)=nam
  Enddo                                                           
  Close(8)
!    Write(g,"(a5,4es12.4)") (nname(i),aa(i),zz(i),nn(i),be(i),i=1,ny)        

! Write binary data file
  Open(8,FILE=trim(data_dir)//'/nuc_data',FORM='unformatted')
  Write(8) ny
  Write(8) t9i
  Write(8) (nname(i),aa(i),zz(i),nn(i),be(i),g(:,i),angm(i),i=1,ny)
  Close(8)

! Deallocate partition function data
  Deallocate (g,angm)

  Return                                                                    
End Subroutine format_nuclear_data
