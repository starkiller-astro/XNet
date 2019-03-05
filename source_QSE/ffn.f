!*******************************************************************************
!  FFN.f for XNet
!  This file contains reaction data structures for Fuller, Fowler, Neuman 
!  (FFN; 1982,1985) formated weak reactions and routines to read the data 
!  and calculate the reaction rates.
!*******************************************************************************

      module ffn_data
!-----------------------------------------------------------------------------
!  This module contains the data to calculate FFN formatted weak reactions.
!-----------------------------------------------------------------------------
      integer :: nffn ! The number of FFN reactions
      real(8), dimension(:),   allocatable :: rf,r1,r2       ! dim(nffn)
      real(8), dimension(:,:), allocatable :: ffnsum,ffnenu  ! dim(nffn,143)

! Threading Scope
!$OMP THREADPRIVATE(rf,r1,r2)

      end module ffn_data

      subroutine read_ffn_data(nin,data_dir)
!-----------------------------------------------------------------------------  
!  This routine allocates and loads the data structures for FFN reaction rates 
!-----------------------------------------------------------------------------  
      use ffn_data
      integer nin
      character (LEN=*) :: data_dir
      Open(3,FILE=trim(data_dir)//"/netweak",status='old')
      nffn=nin
!$OMP PARALLEL
      Allocate (rf(nffn),r1(nffn),r2(nffn))
!$OMP END PARALLEL
      Allocate (ffnsum(nffn,143),ffnenu(nffn,143))
      Do i=1,nffn
        Read(3,"(a1)") null  
        Read(3,"(9(f8.3))") (ffnsum(i,j),ffnenu(i,j),j=1,143)  
      Enddo
      End subroutine read_ffn_data

      subroutine ffn_rate(t9,ene) 
!-----------------------------------------------------------------------------  
!  This routine calculates the reaction rates for FFN weak rates
!-----------------------------------------------------------------------------  
      use ffn_data
      real(8) :: t9,ene,tg(13),egl(11),enl,dt,de,ddt
      integer i1(4)
      integer i,le1,lt1
      data tg/0.01,0.1,0.2,0.4,0.7,1.0,1.5,2.,3.,5.,10.,30.,100./               
      data egl/1.,2.,3.,4.,5.,6.,7.,8.,9.,10.,11./                              
      Do i=1,13
        If(t9<=tg(i)) Exit
      Enddo
      lt1=i-1
      If(lt1<1) lt1=1
      enl=log10(ene)                                                           
      le1=int(enl)
      If(le1<1) le1=1
      dt=tg(lt1+1)-tg(lt1)                                                      
      de=egl(le1+1)-egl(le1)                                                    
      ddt=t9-tg(lt1)                                                            
      i1(1)=13*(le1-1)+lt1                                                      
      i1(2)=i1(1)+1                                                             
      i1(3)=13*le1+lt1                                                          
      i1(4)=i1(3)+1                                                             
      Do i=1,nffn                                                              
        r1(i)=ffnsum(i,i1(1))+(ffnsum(i,i1(2))-ffnsum(i,i1(1)))/dt*ddt
        r2(i)=ffnsum(i,i1(3))+(ffnsum(i,i1(4))-ffnsum(i,i1(3)))/dt*ddt
        rf(i)=r1(i)+(r2(i)-r1(i))/de*(enl-egl(le1))    
        If(rf(i).lt.-30.) Then
          rf(i)=0.
        Else
          rf(i)=10.**rf(i)     
        Endif
      Enddo                                                                     
      Return                                                                    
      End subroutine ffn_rate

