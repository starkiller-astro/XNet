!*******************************************************************************
! Jacobian Dense, part of XNet 7, 5/25/10 
!
! The routines in this file assume a dense Jacobian and use a variety of linear 
! algebra packages.  Choice of matrix solver is determined by uncommented lines.
!
! The bulk of the computational cost of the network (60-95%) is the solving 
! of the matrix equation.  Careful selection of the matrix solver is therefore 
! very important to fast computation.  For networks from a few dozen up to a 
! couple hundred species, hand tuned dense solvers such as those supplied by the 
! hardware manufacturer (often LAPACK) or third-parties like NAG, IMSL, etc. are 
! the fastest. However for larger matrices, sparse solvers are faster.  
!*******************************************************************************
  
Subroutine qread_jacobian_data(data_dir)
!===============================================================================
! Initializes the Jacobian data. For the dense solver, the only jacobian data is 
! the matrix dimension, ny.
!===============================================================================
  Use qjacobian_data
  Use qse_data
  Character (LEN=*),  Intent(in)  :: data_dir

! Allocate (jac(gp0i,gp0i),indx(gp0i)))
End Subroutine qread_jacobian_data
  
Subroutine qjacobian_build(diag,mult)
!===============================================================================
! This routine calculates the reaction Jacobian matrix, dYdot/dY, and augments 
! by multiplying all elements by mult and adding diag to the diagonal elements.
!===============================================================================
  Use controls
  Use qjacobian_data
  Use nuclear_data
  Use qse_data
  Use reac_rate_data
  Use conditions
  Use abundances
  Use qse_abundances
  Use part_funct_data
  integer :: i,j,kstep,i0,la1,le1,la2,le2,la3,le3,j1,l1,l2,l3
  real(8) :: alph,beta,rdt,d,Cpart,dummy,hold1,hold2,hold3,hold4
  real(8) :: yn,yz,yt1,yt2,ylog1,ylog2,yp
  real(8) :: cross1(gn),same1,Qse3(ny),Qse4(ny)
  integer :: qn,qz
  integer :: n3,z3,n4,z4,k
  integer :: nsi,zsi,nfe,zfe,ncr,zcr,info
  Real(8), Intent(IN) :: diag, mult
  integer over, clock_start,clock_end,clock_rate
  double precision  elapsed_time,accum2,accum1
! write(*,*) "in jacobian build"
!       If(.not.allocated(f)) Allocate(indx(ny),f(ny),ydot0(ny),um(ny))
!       If(.not.allocated(f)) Allocate(indx(gp0i),
!     & f(gp0i),ydot0(ny))
!       If(.not.allocated(am)) Allocate(am(ny,ny),dydot(gpi0,gpi0))
!       If(.not.allocated(am)) Allocate(am(gp0i,gp0i),dydot(gp0i,gp0i))
!       Allocate(indx4(gn:gn))

 If(.not.allocated(jac)) Allocate (jac(gp0i,gp0i))
!         write(*,*) "tag gp0i", gp0i
      jac=0.0
      beta=1.0
       If(gn>2)then   ! This needs to be  commneted out when weak rxns are used.
      alph=0.0
      accum1=0
      accum2=0
      part1=0.0
      part2i=0.0
      part2=0.0
      cross1=0.0

! Calculate the reaction rates and abundance time derivatives

      call qyderiv
!  Build the Jacobian
! regrp1-3 (reactant group) tells the group membership, lt,si of Fe,
! for each reactant so that the proper
! derivative is used (dY/dYn,dY/dSi28, etc).
! The regrp's are set up in data.f.
! part1 is dYgdot/dYR.
! part2 is dYR/dYG.
 !     diag rdt=1.0/tdel
        dydot=0.0
        yn=1/yt(1)
        yp=1/yt(2)
        yt1=yt(1)
        yt2=yt(2)
        ylog1 = dlog(yt(1))
        ylog2 = dlog(yt(2))

        Do i=1,ny
         If(ID(i)==3)then
           Qse3(i)=dexp(dlCqse(i)-dlCqse(f3)+ylog1*(drvN(i)) &
     &     +ylog2*(drvZ(i)))
           Qse4(i)=0
         ElseIf(ID(i)==4)then
           Qse4(i)=dexp(dlCqse(i)-dlCqse(f4)+ylog1*(drvN(i)) &
     &     +ylog2*(drvZ(i)))
           Qse3(i)=0
         ElseIf(ID(i)<3)then
          Qse3(i)=0
          Qse4(i)=0
         Endif
!      write(*,*) "tag",kstep,i,yt1,yt2, Qse3(i),Qse4(i)
        Enddo

! Start building the light groups pieces for eq 12.
       Do j1=laq(1,2),leq(1,2)
          l1= qn11(j1)
          cross1=0
          same1=qb1(j1)*yt(l1)
        If(ID(l1)>0)then
          cross1(1)=same1*drvN(l1)*yn
          cross1(2)=same1*drvZ(l1)*yp
           If(ID(l1)==3)then
             cross1(3)=qb1(j1)*Qse3(l1)
           ElseIf(ID(l1)==4)then
              cross1(4)=qb1(j1)*Qse4(l1)
           Else
              cross1(3)=0
              cross1(4)=0
           Endif
            Do j=1,gn
              part1(1,j)=part1(1,j)-qnnum(j1)*cross1(j)
              part1(2,j)=part1(2,j)-qznum(j1)*cross1(j)
            Enddo
            cross1=0
        ElseIf(ID(l1)==0) then
            part1(1,singlenuc(l1))=part1(1,singlenuc(l1)) &
     &       -qnnum(j1)*qb1(j1)
            part1(2,singlenuc(l1))=part1(2,singlenuc(l1)) &
     &       -qznum(j1)*qb1(j1)
        Endif
       Enddo ! j1

! now build the Si or, more generally, the second group 1p pieces for eq 12
!      If(gn>2)then ! Sould be uncommented for weak rxns
       Do i=3,gn
           qn=nn(gp0number(i))
           qz=zz(gp0number(i))
         Do j1=laq(1,i),leq(1,i)
            cross1=0
            l1= qn11(j1)
            same1=qb1(j1)*yt(l1)
            If(ID(l1)>0)then
              cross1(1)=same1*drvN(l1)*yn
              cross1(2)=same1*drvZ(l1)*yp
           If(ID(l1)==3)then
             cross1(3)=qb1(j1)*Qse3(l1)
            ElseIf(ID(l1)==4)then
             cross1(4)=qb1(j1)*Qse4(l1)
           Else
              cross1(3)=0
              cross1(4)=0
          Endif
            Do j=1,gn
              part1(1,j)=part1(1,j)-(qnnum(j1)-qn)*cross1(j)
              part1(2,j)=part1(2,j)-(qznum(j1)-qz)*cross1(j)
              part1(i,j)=part1(i,j)-cross1(j)
            Enddo
         ElseIf(ID(l1)==0) then
            part1(1,singlenuc(l1))=part1(1,singlenuc(l1)) &
     &      -(qnnum(j1)-qn)*qb1(j1)
            part1(2,singlenuc(l1))=part1(2,singlenuc(l1)) &
     &      -(qznum(j1)-qz)*qb1(j1)
            part1(i,singlenuc(l1))=part1(i,singlenuc(l1))-qb1(j1)
          Endif
        Enddo !j1
       Enddo !i
!      Endif !gn
! Now Build the single nuc square of Jacobian (larger lower right block).
       Do i=gn+1,gp0i
            cross1=0
           Do j1=laq(1,i),leq(1,i)
              l1= qn11(j1)
              same1=qb1(j1)*yt(l1)
               If(ID(l1)>0)then
                  cross1(1)=same1*drvN(l1)*yn
                  cross1(2)=same1*drvZ(l1)*yp
                If(ID(l1)==3)then
                   cross1(3)=qb1(j1)*Qse3(l1)
                ElseIf(ID(l1)==4)then
                  cross1(4)=qb1(j1)*Qse4(l1)
                Else
                 cross1(3)=0
                 cross1(4)=0
                Endif
                Do j=1,gn
                  part1(i,j)=part1(i,j)-cross1(j)
               Enddo
              ElseIf(ID(l1)==0) then
               part1(i,singlenuc(l1))=part1(i,singlenuc(l1))-qb1(j1)
              Endif
           Enddo !j1
        Enddo !i
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Now for the 2 particle reactions
!regrp(i,1) - first reactant
!regrp(i,2) - second reactant
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      Do j1=laq(2,2),leq(2,2)
          l1= qn21(j1)
          l2= qn22(j1)
           cross1=0.0 ! Fix this !!!!
          same1=qb2(j1)*yt(l1)*yt(l2)

          If(ID(l1)>0)then
            cross1(1)=same1*yn*(drvN(l1))
            cross1(2)= same1*yp*(drvZ(l1))
             If(ID(l1)==3)then
                cross1(3)=qb2(j1)*yt(l2)*Qse3(l1)
             ElseIf(ID(l1)==4)then
                cross1(4)=qb2(j1)*yt(l2)*Qse4(l1)
             Endif

          Elseif(ID(l1)==0) then
                part1(1,singlenuc(l1))=part1(1,singlenuc(l1)) &
     &                -qnnum2(j1)*qb2(j1)*yt(l2)
                part1(2,singlenuc(l1))=part1(2,singlenuc(l1)) &
     &                -qznum2(j1)*qb2(j1)*yt(l2)
          Endif

          If(ID(l2)>0)then
            cross1(1)=cross1(1)+ same1*yn*(drvN(l2))
            cross1(2)=cross1(2)+ same1*yp*(drvZ(l2))
              If(ID(l2)==3)then
                 cross1(3)=cross1(3)+qb2(j1)*yt(l1)*Qse3(l2)
              ElseIf(ID(l2)==4)then
                 cross1(4)=cross1(4)+qb2(j1)*yt(l1)*Qse4(l2)
              Endif

          ElseIf(ID(l2)==0) then
            part1(1,singlenuc(l2))=part1(1,singlenuc(l2)) &
     &           -qnnum2(j1)*qb2(j1)*yt(l1)
            part1(2,singlenuc(l2))=part1(2,singlenuc(l2)) &
     &            -qznum2(j1)*qb2(j1)*yt(l1)
           Endif !l2

            Do j=1,gn
              part1(1,j)=part1(1,j)-qnnum2(j1)*cross1(j)
              part1(2,j)=part1(2,j)-qznum2(j1)*cross1(j)
            Enddo
      Enddo !j1

!      If(gn>2)then ! Should be uncommented for weak rxns
         Do i=3,gn
           qn=nn(gp0number(i))
           qz=zz(gp0number(i))

           Do j1=laq(2,i),leq(2,i)
              l1= qn21(j1)
              l2= qn22(j1)
              cross1=0.0
              same1=qb2(j1)*yt(l1)*yt(l2)
              If(ID(l1)>0)then
                cross1(1)=same1*yn*(drvN(l1))
                cross1(2)=same1*yp*(drvZ(l1))
                 If(ID(l1)==3)then
                    cross1(3)=qb2(j1)*yt(l2)*Qse3(l1)
                 ElseIf(ID(l1)==4)then
                   cross1(4)=qb2(j1)*yt(l2)*Qse4(l1)
                 Endif

              Elseif(ID(l1)==0) then
                     part1(1,singlenuc(l1))=part1(1,singlenuc(l1))  &
     &                      -(qnnum2(j1)-qn)*qb2(j1)*yt(l2)
                      part1(2,singlenuc(l1))=part1(2,singlenuc(l1)) &
     &                      -(qznum2(j1)-qz)*qb2(j1)*yt(l2)
                      part1(i,singlenuc(l1))=part1(i,singlenuc(l1)) &
     &                      -qb2(j1)*yt(l2)
              Endif ! l1
              If(ID(l2)>0)then
                cross1(1)=cross1(1)+same1*yn*(drvN(l2))
                cross1(2)=cross1(2)+same1*yp*(drvZ(l2))
                    If(ID(l2)==3)then
                       cross1(3)=cross1(3)+qb2(j1)*yt(l1)*Qse3(l2)
                 ElseIf(ID(l2)==4)then
                   cross1(4)=cross1(4)+qb2(j1)*yt(l1)*Qse4(l2)
                 Endif
              ElseIf(ID(l2)==0) then
                  part1(1,singlenuc(l2))=part1(1,singlenuc(l2)) &
     &                  -(qnnum2(j1)-qn)*qb2(j1)*yt(l1)
                  part1(2,singlenuc(l2))=part1(2,singlenuc(l2)) &
     &                  -(qznum2(j1)-qz)*qb2(j1)*yt(l1)
                  part1(i,singlenuc(l2))=part1(i,singlenuc(l2)) &
     &                  -qb2(j1)*yt(l1)
              Endif !j1

            Do j=1,gn
              part1(1,j)=part1(1,j)-(qnnum2(j1)-qn)*cross1(j)
              part1(2,j)=part1(2,j)-(qznum2(j1)-qz)*cross1(j)
              part1(i,j)=part1(i,j)-cross1(j)
!              write(*,*) "tag2",i,j, cross1(j),gn
         Enddo
      Enddo
      Enddo !i
!23456789012345678988888888888888888888888888888888888888888888888888888888888888

      Do i=gn+1,gp0i
        Do j1=laq(2,i),leq(2,i)
           l1= qn21(j1)
           l2= qn22(j1)
           cross1=0
           same1=qb2(j1)*yt(l1)*yt(l2)
           If(ID(l1)>0)then
             cross1(1)=same1*yn*(drvN(l1))
             cross1(2)=same1*yp*(drvZ(l1))

             If(ID(l1)==3)then
               cross1(3)=qb2(j1)*yt(l2)*Qse3(l1)
             ElseIf(ID(l1)==4)then
               cross1(4)=qb2(j1)*yt(l2)*Qse4(l2)
             Endif

           ElseIf (ID(l1)==0)then
              part1(i,singlenuc(l1))=part1(i,singlenuc(l1)) &
     &        -qb2(j1)*yt(l2)
           Endif

           If(ID(l2)>0)then
              cross1(1)=cross1(1)+same1*yn*(drvN(l2))
              cross1(2)=cross1(2)+same1*yp*(drvZ(l2))
              If(ID(l2)==3)then
                cross1(3)=cross1(3)+qb2(j1)*yt(l1)*Qse3(l2)
              ElseIf(ID(l2)==4)then
                cross1(4)=cross1(4)+qb2(j1)*yt(l1)*Qse4(l2)
              Endif
            ElseIf(ID(l2)==0) then
              part1(i,singlenuc(l2))=part1(i,singlenuc(l2)) &
     &        -qb2(j1)*yt(l1)
            Endif

            Do j=1,gn
              part1(i,j)=part1(i,j)-cross1(j)
            Enddo
           Enddo !j1
        Enddo !i
!      Endif !gn
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Now for the 3 particle reactions
! Many of the following are not necessary since
! the only 3-p rxns in the current network involve reactants
! (lt lt single), and (lt lt lt) .
! This is more genral for any configuration of the three possible group
! reactants.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!2345678901234567890
      Do j1=laq(3,2),leq(3,2)
         l1= qn31(j1)
         l2= qn32(j1)
         l3= qn33(j1)
         same1=qb3(j1)*yt(l1)*yt(l2)*yt(l3)
         cross1=0
         If(ID(l1)>0)then
            cross1(1)=same1*yn*(drvN(l1))
            cross1(2)=same1*yp*(drvZ(l1))
            If(ID(l1)==3)then
             cross1(3)=qb3(j1)*yt(l2)*yt(l3)*(dlCqse(l1)-dlCqse(f3))
            ElseIf(ID(l1)==4)then
             cross1(4)=qb3(j1)*yt(l2)*yt(l3)*(dlCqse(l1)-dlCqse(f4))
            Endif
         Elseif(ID(l1)==0) then
           part1(1,singlenuc(l1))=part1(1,singlenuc(l1)) &
     &                 -qnnum3(j1)*qb3(j1)*yt(l2)*yt(l3)
           part1(2,singlenuc(l1))=part1(2,singlenuc(l1)) &
     &                  -qznum3(j1)*qb3(j1)*yt(l2)*yt(l3)
         Endif
         If(ID(l2)>0)then
            cross1(1)=cross1(1)+same1*yn*(drvN(l2))
            cross1(2)=cross1(2)+same1*yp*(drvZ(l2))
            If(ID(l2)==3)then
              cross1(3)=cross1(3)+qb3(j1)*yt(l1)*yt(l3)*Qse3(l2)
            ElseIf(ID(l2)==4)then
              cross1(4)=cross1(4)+qb3(j1)*yt(l1)*yt(l3)*Qse4(l2)
            Endif
         Elseif(ID(l2)==0)then
            part1(1,singlenuc(l2))=part1(1,singlenuc(l2)) &
     &                 -qnnum3(j1)*qb3(j1)*yt(l1)*yt(l3)
            part1(2,singlenuc(l2))=part1(2,singlenuc(l2)) &
     &                  -qznum3(j1)*qb3(j1)*yt(l1)*yt(l3)
         Endif

         If(ID(l3)>0)then
            cross1(1)=cross1(1)+same1*yn*(drvN(l3))
            cross1(2)=cross1(2)+same1*yp*(drvZ(l3))
            If(ID(l3)==3)then
               cross1(3)=cross1(3)+qb3(j1)*yt(l1)*yt(l2)*Qse3(l3)
            ElseIf(ID(l3)==4)then
               cross1(4)=cross1(4)+qb3(j1)*yt(l1)*yt(l2)*Qse4(l3)
            Endif
          Elseif(ID(l3)==0)then
             part1(1,singlenuc(l3))=part1(1,singlenuc(l3)) &
     &                 -qnnum3(j1)*qb3(j1)*yt(l1)*yt(l2)
             part1(2,singlenuc(l3))=part1(2,singlenuc(l3)) &
     &                  -qznum3(j1)*qb3(j1)*yt(l1)*yt(l2)
           Endif

            Do j=1,gn
              part1(1,j)=part1(1,j)-qnnum3(j1)*cross1(j)
              part1(2,j)=part1(2,j)-qznum3(j1)*cross1(j)
            Enddo
       Enddo !j1

!      If(gn>2) then shoud be uncommented for weak rxns.
      Do i=3,gn
           qn=nn(gp0number(i))
           qz=zz(gp0number(i))
       Do j1=laq(3,i),leq(3,i)
          l1= qn31(j1)
          l2= qn32(j1)
          l3= qn33(j1)
          cross1=0
          same1=qb3(j1)*yt(l1)*yt(l2)*yt(l3)
          If(ID(l1)>0)then
            cross1(1)=same1*yn*(drvN(l1))
            cross1(2)=same1*yp*(drvZ(l1))
            If(ID(l1)==3)then
              cross1(3)=qb3(j1)*yt(l2)*yt(l3)*Qse3(l1)
            ElseIf(ID(l1)==4)then
              cross1(4)=qb3(j1)*yt(l2)*yt(l3)*Qse4(l1)
            Endif
          Elseif(ID(l1)==0) then
            part1(1,singlenuc(l1))=part1(1,singlenuc(l1))    &
     &                 -(qnnum3(j1)-qn)*qb3(j1)*yt(l2)*yt(l3)
            part1(2,singlenuc(l1))=part1(2,singlenuc(l1))    &
     &                  -(qznum3(j1)-qz)*qb3(j1)*yt(l2)*yt(l3)
            part1(i,singlenuc(l1))=part1(i,singlenuc(l1))    &
     &                     -qb3(j1)*yt(l2)*yt(l3)
           Endif
           If(ID(l2)>0)then
              cross1(1)=cross1(1)+same1*yn*(drvN(l2))
              cross1(2)=cross1(2)+same1*yp*(drvZ(l2))
              If(ID(l2)==3)then
                 cross1(3)=cross1(3)*qb3(j1)*yt(l1)*yt(l3)*Qse3(l2)
              ElseIf(ID(l2)==4)then
                 cross1(4)=cross1(3)*qb3(j1)*yt(l1)*yt(l3)*Qse4(l2)
              Endif
            Elseif(ID(l2)==0)then
              part1(1,singlenuc(l2))=part1(1,singlenuc(l2))  &   
     &                 -(qnnum3(j1)-qn)*qb3(j1)*yt(l1)*yt(l3)
              part1(2,singlenuc(l2))=part1(2,singlenuc(l2))  &
     &                  -(qznum3(j1)-qz)*qb3(j1)*yt(l1)*yt(l3)
              part1(i,singlenuc(l2))=part1(i,singlenuc(l2))  &
     &                     -qb3(j1)*yt(l1)*yt(l3)
              Endif

       If(ID(l3)>0)then
         cross1(1)=cross1(1)+same1*yn*(drvN(l3))
         cross1(2)=cross1(2)+same1*yp*(drvZ(l3))
         If(ID(l3)==3)then
           cross1(3)=cross1(3)*qb3(j1)*yt(l1)*yt(l2)*Qse3(l3)
         ElseIf(ID(l3)==4)then
           cross1(4)=cross1(4)*qb3(j1)*yt(l1)*yt(l2)*Qse4(l3)
         Endif
       Elseif(ID(l3)==0)then
          part1(1,singlenuc(l3))=part1(1,singlenuc(l3))       &
     &                 -(qnnum3(j1)-qn)*qb3(j1)*yt(l1)*yt(l2)
          part1(2,singlenuc(l3))=part1(2,singlenuc(l3))       &
     &                  -(qznum3(j1)-qz)*qb3(j1)*yt(l1)*yt(l2)
          part1(i,singlenuc(l3))=part1(i,singlenuc(l3))       &
     &                     -qb3(j1)*yt(l2)*yt(l1)
       Endif

            Do j=1,gn
              part1(1,j)=part1(1,j)-(qnnum3(j1)-qn)*cross1(j)
              part1(2,j)=part1(2,j)-(qznum3(j1)-qz)*cross1(j)
              part1(i,j)=part1(i,j)-cross1(j)
            Enddo
       Enddo !j1
      Enddo !i
!!large block
!2345678901234567890999999999999999999999999999999999999999999999999
      Do i=gn+1,gp0i
        Do j1=laq(3,i),leq(3,i)
           l1= qn31(j1)
           l2= qn32(j1)
           l3= qn33(j1)
           cross1=0.0
           same1=qb3(j1)*yt(l1)*yt(l2)*yt(l3)
           If(ID(l1)>0)then
             cross1(1)=same1*yn*(drvN(l1))
             cross1(2)=same1*yp*(drvZ(l1))
             If(ID(l1)==3)then
               cross1(3)=qb3(j1)*yt(l2)*yt(l3)*Qse3(l1)
             Endif
             If(ID(l1)==4)then
               cross1(4)=qb3(j1)*yt(l2)*yt(l3)*Qse4(l1)
             Endif
          Elseif(ID(l1)==0) then
            part1(i,singlenuc(l1))=part1(i,singlenuc(l1)) &
     &           -qb3(j1)*yt(l2)*yt(l3)
          Endif
          If(ID(l2)>0)then
             cross1(1)=cross1(1)+same1*yn*(drvN(l2))
             cross1(2)=cross1(2)+same1*yp*(drvZ(l2))
             If(ID(l2)==3)then
               cross1(3)=cross1(3)+qb3(j1)*yt(l1)*yt(l3)*Qse3(l2)
             Endif
             If(ID(l2)==4)then
               cross1(4)=cross1(4)+qb3(j1)*yt(l1)*yt(l3)*Qse4(l2)
             Endif
           Elseif(ID(l2)==0)then
             part1(i,singlenuc(l2))=part1(i,singlenuc(l2)) &
     &           -qb3(j1)*yt(l1)*yt(l3)
           Endif
           If(ID(l3)>0)then
              cross1(1)=cross1(1)+same1*yn*(drvN(l3))
              cross1(2)=cross1(2)+same1*yp*(drvZ(l3))
              If(ID(l3)==3)then
                cross1(3)=cross1(3)+qb3(j1)*yt(l1)*yt(l2)*Qse3(l3)
              Endif
              If(ID(l3)==4)then
                cross1(4)=cross1(4)+qb3(j1)*yt(l1)*yt(l2)*Qse4(l3)
              Endif
           Elseif(ID(l3)==0)then
             part1(i,singlenuc(l3))=part1(i,singlenuc(l3)) &
     &           -qb3(j1)*yt(l2)*yt(l1)
          Endif

            Do j=1,gn
               part1(i,j)=part1(i,j)-cross1(j)
            Enddo
         Enddo !j1
        Enddo !i
!       Endif !gn
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Now build Part2 the pieces of the Derivative of YG with Respect to YG.
! Part2i contains the building blocks for each element of 4 X 4 martix.
! Part2 contains the 4-vectors.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!write(*,*) "tag building part2"
! Light group pieces.
           n3=nn(f3)
           z3=zz(f3)
           Part2i=0.0
           Do i=1,gp1i
       part2i(1,1)=part2i(1,1)+nn(i)*nn(i)*yt(i)*yn
       part2i(1,2)=part2i(1,2)+nn(i)*zz(i)*yt(i)*yp
       part2i(2,1)=part2i(2,1)+nn(i)*zz(i)*yt(i)*yn
       part2i(2,2)=part2i(2,2)+zz(i)*zz(i)*yt(i)*yp
              Enddo
!Si group pieces
         Do ig =1,gp2i
           i=gp2number(ig)
           nsi=nn(i)-n3
           zsi=zz(i)-z3
           Cpart=dexp((dlCqse(i)-dlCqse(f3))+ylog1*nsi+ylog2*zsi)

           part2i(1,1)=part2i(1,1)+nsi*nsi*yt(i)*yn
           part2i(1,2)=part2i(1,2)+nsi*zsi*yt(i)*yp
           part2i(1,3)=part2i(1,3)+nsi*Cpart
           part2i(2,1)=part2i(2,1)+nsi*zsi*yt(i)*yn
           part2i(2,2)=part2i(2,2)+zsi*zsi*yt(i)*yp
           part2i(2,3)=part2i(2,3)+zsi*Cpart
           part2i(3,1)=part2i(3,1)+nsi*yt(i)*yn
           part2i(3,2)=part2i(3,2)+zsi*yt(i)*yp
           part2i(3,3)=part2i(3,3)+Cpart

         Enddo
         If(numberofgroups>=3)then
           n4=nn(f4)
          z4=zz(f4)
!Fe group pieces
          Do ig =1,gp3i
           i=gp3number(ig)
           nfe=nn(i)-n4
           zfe=zz(i)-z4
           Cpart=dexp(dlCqse(i)-dlCqse(f4)+ylog1*nfe+ylog2*zfe)
           part2i(1,1)=part2i(1,1)+nfe*nfe*yt(i)*yn
           part2i(1,2)=part2i(1,2)+nfe*zfe*yt(i)*yp
           part2i(1,4)=part2i(1,4)+nfe*Cpart
           part2i(2,1)=part2i(2,1)+nfe*zfe*yt(i)*yn
           part2i(2,2)=part2i(2,2)+zfe*zfe*yt(i)*yp
           part2i(2,4)=part2i(2,4)+zfe*Cpart
           part2i(4,1)=part2i(4,1)+nfe*yt(i)*yn
           part2i(4,2)=part2i(4,2)+zfe*yt(i)*yp
           part2i(4,4)=part2i(4,4)+Cpart
         Enddo
        Endif
! Initialize variables. The four vectors will go in the 4x4 part2 matrix.
         part2=0.0
!write(*,*) "tag about to do small sove part2"
          DO j=1,gn
             part2(j,j) =1
          Enddo
!          write(*,*) "gn", gn
          call dgesv(gn,gn,part2i,gn,indx4,part2,gn,info)
!write(*,*) "tag after solve part2"

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Make Matrix
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!         write(*,*) "tag2 gp0i", gp0i
! diag = rdt=1.0/tdel/beta
      Do i=1,gp0i
           Do j=1,gn
             Do k=1,gn
              jac(i,j) = jac(i,j)+part1(i,k)*part2(k,j)
             Enddo
            Enddo
           Do j=gn+1,gp0i
              jac(i,j)=part1(i,j)
            Enddo
        Enddo
      Endif ! gn
! Augment matrix with externally provided factors  
      Do i=1,gp0i
         jac(i,i)=jac(i,i)+diag
      Enddo
!write(*,*) "tag after matrix make with jac"



  If(idiag>=5) Then
    Write(lun_diag,"(a9,2es14.7)") 'JAC_build',diag,mult
    Do i=1,ny
      Write(lun_diag,"(14es9.1)") (jac(i,j),j=1,ny)
    EndDo
  EndIf
  
! Test the eigenvalues 
! If(idiag>=6) Then
!   Call eigen_test(kstep,jac,rdt)
! EndIf
! write(*,*) 'at the end of jacobain build gp01', gp0i
!write (*,*) 'jacobian', jac
  Return   
End Subroutine qjacobian_build
  
Subroutine qjacobian_solve(kstep,rhs,dy) 
!===============================================================================
! This routine solves the system of abundance equations composed of the jacobian
! matrix and rhs vector.
!===============================================================================
  Use controls
  Use qjacobian_data
  Use qse_data
  Integer, Intent(in)  :: kstep
  Real(8), Intent(in)  :: rhs(gp0i)
  Real(8), Intent(out) ::  dy(gp0i)
  Real(8) :: d 
  Integer :: info,i,j
!  write(*,*) "in solve"
  If(.not.allocated(indx)) Allocate(indx(gp0i)) 
!  write(*,*) "in solve"
  
! To Use Num Rec LU Decomp, uncomment 3 lines below
! Call ludcmp(jac,gp0i,gp0i,indx,d)
! Call lubksb(jac,gp0i,gp0i,indx,rhs)
! dy=rhs
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Debug
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Do i= 1,gp0i
!   do j= 1,gp0i
!        write(900,*) i,j, jac(i,j)
!   enddo
!enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!










! To Use LAPACK solver,  uncomment one of the first two lines below and the third line 
! Call sgesv(ny,1,jac,ny,indx,rhs,ny,info) ! Single precision version
!  write(*,*) "before dgesv",gp0i
  Call dgesv(gp0i,1,jac,gp0i,indx,rhs,gp0i,info) ! Double precision version
!  write(*,*) "after  dgesv"
  dy=rhs
  
!  write(*,*) "after dy    "
! Diagnostic output
  If(idiag>=4) Then
    Write(lun_diag,"(a)") 'JAC_SOLV'
    Write(lun_diag,"(14es10.3)") dy
  EndIf

  Return                                                                    
End Subroutine qjacobian_solve                                                                       
  
Subroutine qjacobian_decomp(kstep) 
!===============================================================================
! This routine performs a matrix decomposition for the jacobian
!===============================================================================
  Use controls
  Use qjacobian_data
  Use qse_data
  Integer, Intent(in)  :: kstep
  Real(8) :: d 
  Integer :: i,j,info
  
! To Use Num Rec LU Decomp, uncomment line below
! Call ludcmp(jac,gp0i,gp0i,indx,d)
  
! To Use LAPACK ,  uncomment one of the first two lines below 
! Call sgetrf(gp0i,gp0i,jac,gp0i,indx,info)! Single precision version
  Call dgetrf(gp0i,gp0i,jac,gp0i,indx,info)
  
! Diagnostic output
  If(idiag>=4) Then
    Write(lun_diag,"(a3,i4)") 'LUD',info
    Do i=1,ny
      Write(lun_diag,"(14es9.1)") (jac(i,j),j=1,gp0i)
    EndDo
  EndIf
  
  Return                                                                    
End Subroutine qjacobian_decomp                                                                       
  
Subroutine qjacobian_bksub(rhs,dy) 
!===============================================================================
! This routine performs back-substitution for a previously factored matrix and 
! the vector rhs.   
!===============================================================================
  Use controls
  Use qjacobian_data
  Use qse_data
  Real(8), Intent(IN) :: rhs(gp0i)
  Real(8), Intent(out) ::  dy(gp0i)
  Real(8) :: d 
  Integer :: i,info
  
! To Use Num Rec LU Decomp, uncomment 2 lines below
! Call lubksb(jac,gp0i,gp0i,indx,rhs)
! dy=rhs
  
! To use LAPACK solver, uncomment one of the first two lines below and the third line 
! Call sgetrs('No transpose',gp0i,1,jac,gp0i,indx,rhs,gp0i,info) ! Single precision version
  Call dgetrs('No transpose',gp0i,1,jac,gp0i,indx,rhs,gp0i,info) ! Double precision version
  dy=rhs
  
! Diagnostic output
  If(idiag>=4) Then
    Write(lun_diag,"(a)") 'BKSUB'
    Write(lun_diag,"(14es10.3)") dy
  EndIf

  Return                                                                    
End Subroutine qjacobian_bksub                                                                       
  
