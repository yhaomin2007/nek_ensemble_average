c-----------------------------------------------------------------------
      subroutine savg(ua,u,n,hndl)

c     Compute the average of quantity u() across sessions

      include 'SIZE'
      include 'TOTAL'

      real u,ua,weight
      integer hndl

      call copy(ua,u,n)
      call msgsop(ua,'sum',hndl)   ! Sum ua() across sessions

      weight = 1./nsessions
c      weight = 1.0

      call cmult(ua,weight,n)   ! Scale with 1./nsession for sess avg

      return
      end
c-----------------------------------------------------------------------
      subroutine msgsop_get_hndl(hndl,nel,lx,ly,lz)

c     Get a multi-session gsop handle

      include 'SIZE'
      include 'TOTAL'
      include 'mpif.h'

      integer e,eg

      common /c_is1/ glo_num(lx1*ly1*lz1*lelt)
      integer*8 glo_num
      integer hndl

      n = nel*lx*ly*lz

      do e=1,nel
         eg = lglel(e)
         do k=1,lz      
         do j=1,ly
         do i=1,lx
            ii = i + lx*(j-1) + lx*ly*(k-1) + lx*ly*lz*(e-1)
            glo_num(ii) = i + lx*(j-1) + lx*ly*(k-1) + 
     $                                   lx*ly*lz*(eg-1)
         enddo
         enddo
         enddo
      enddo


      call mpi_comm_size(mpi_comm_world,np_global,ierr)
      call fgslib_gs_setup(hndl,glo_num,n,mpi_comm_world,np_global)

      return
      end
c-----------------------------------------------------------------------
      subroutine msgsop(u,op,hndl)

c     multi-session version of gsop

      include 'SIZE'
      include 'INPUT'
      include 'TSTEP'
      include  'CTIMER'

      real u(1)
      character*3 op

      if(ifsync) call nekgsync()

      if (op.eq.'+  ' .or. op.eq.'sum' .or. op.eq.'SUM')
     &   call fgslib_gs_op(hndl,u,1,1,0)

      if (op.eq.'*  ' .or. op.eq.'mul' .or. op.eq.'MUL')
     &   call fgslib_gs_op(hndl,u,1,2,0)

      if (op.eq.'m  ' .or. op.eq.'min' .or. op.eq.'mna' 
     &                .or. op.eq.'MIN' .or. op.eq.'MNA')
     &   call fgslib_gs_op(hndl,u,1,3,0)

      if (op.eq.'M  ' .or. op.eq.'max' .or. op.eq.'mxa'
     &                .or. op.eq.'MAX' .or. op.eq.'MXA')
     &   call fgslib_gs_op(hndl,u,1,4,0)

      return
      end
c------------------------------------------------------------------------
      subroutine multimesh_create

c     Dummy for ensemble

      return
      end
c------------------------------------------------------------------------
      subroutine userchk_set_xfer

c     Dummy for ensemble

      return
      end
c------------------------------------------------------------------------
      subroutine bcopy

c     Dummy for ensemble

      return
      end
C---------------------------------------------------------------------
      subroutine setintercomm(nekcommtrue,nptrue) 

c     Dummy for ensemble

      return
      end
c-----------------------------------------------------------------------
      subroutine unsetintercomm(nekcommtrue,nptrue)

c     Dummy for ensemble

      return
      end
c-----------------------------------------------------------------------
      function uglmin(a,n)
      real a(1)

      call happy_check(1)
      call setintercomm(nekcommtrue,nptrue)    ! nekcomm=iglobalcomml
      uglmin=glmin(a,n)
      call unsetintercomm(nekcommtrue,nptrue)  ! nekcomm=nekcomm_original

      return
      end
c-----------------------------------------------------------------------
      function uglamax(a,n)
      real a(1)

      call happy_check(1)
      call setintercomm(nekcommtrue,nptrue)    ! nekcomm=iglobalcomml
      uglamax=glamax(a,n)
      call unsetintercomm(nekcommtrue,nptrue)  ! nekcomm=nekcomm_original

      return
      end
c------------------------------------------------------------------------
      function uglmax(a,n)
      real a(1)

      call happy_check(1)
      call setintercomm(nekcommtrue,nptrue)    ! nekcomm=iglobalcomml
      uglmax=glmax(a,n)
      call unsetintercomm(nekcommtrue,nptrue)  ! nekcomm=nekcomm_original

      return
      end
c-----------------------------------------------------------------------
      subroutine happy_check(ihappy)

c     Dummy for ensemble

      return
      end
c-----------------------------------------------------------------------
      subroutine chk_outflow ! Assign neighbor velocity to outflow if
                             ! characteristics are going the wrong way.

c     Dummy for ensemble

      return
      end
c-----------------------------------------------------------------------
