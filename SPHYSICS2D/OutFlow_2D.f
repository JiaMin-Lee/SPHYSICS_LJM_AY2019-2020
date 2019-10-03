      subroutine outflow_2D
      include 'common.2D'

      double precision dxp_dble

      zmax=0
      ncases=0
      totalmass=0
      tmassout=0
      tmassremain=0
      volumerate=0
      totalmassini=0

      do i=nbp1,np
         totalmassini=totalmassini+pm(i)
         tmassremain=tmassremain+pm(i)*iout(i)
         pmout=0
      end do

      do i=nbp1,np ! Check if any fluid particle is over the box

          if(tmassremain.lt.totalmassini) then
             totalmass=tmassremain
          else
             totalmass=totalmassini
          end if

        if(zp(i).gt.zmax) zmax=zp(i)
          if(iout(i).eq.1.and.(xp(i).gt.0.0.and.
     +      xp(i).lt.0.0015.and.zp(i).lt.0.0)) then
!              i_out=1
!            if(i_periodicOBs(1).eq.1.and.i.gt.nb
!     +              .and.zp(i).gt.zmin_container)then  !Replace Particle on other side of domain
!              if(xp(i).gt.xmax_container)then
!                dxp_dble = dble(xp(i)) - xmax_container_double
!                xp(i) = real(dxp_dble +xmin_container_double)
!              else if(xp(i).lt.xmin_container)then
!                dxp_dble = xmin_container_double - dble(xp(i))
!                xp(i) = real(xmax_container_double - dxp_dble)
!              end if
!            else
            write(80,*) 'Particle Flowing Out: ',i,
     +      '  last X position: ',xp(i), '  last Z-position: ',zp(i)
            write(*,*) 'Particle Flowing Out: ',i,
     +      '  last X position: ',xp(i), '  last Z-position: ',zp(i)


             iout(i)=0
             ncases=ncases+1

             if((time-timeprev).ne.0) then
               tmassremain=totalmass-pm(i)
               tmassout=totalmass-tmassremain
               flowratemass=tmassout/(time-timeprev)
              !totalmassin-tmassremain
               volumerate=tmassout/(1000*(time-timeprev))
               pmprev=pm(i)
               deltime=time-timeprev
               nout=1
               nprev=nout
             else if((time-timeprev).eq.0)then
               tmassremain=totalmass-pm(i)
               tmassout=totalmass-tmassremain+pmprev
               flowratemass=tmassout/(deltime)
               volumerate=tmassout/(1000*deltime)
               pmout=pmprev+pm(i)
               nout=nprev+1
               nprev=nout
               pmprev=pmout
             end if

             if(i_out.ne.1)then
               open(26,file='OutFlow')
               write(26,*) i, xp(i), zp(i), nout, time, tmassremain,
     +                     tmassout, volumerate, flowratemass
               i_out=1
             else
              open(26,file='OutFlow',status='old',POSITION='append')
              write(26,*) i, xp(i), zp(i), nout, time, tmassremain,
     +                    tmassout, volumerate, flowratemass
             end if
             close(26)

             up(i)=0.
             wp(i)=0.
             um1(i)=0.
             uo(i)=0.
             wm1(i)=0.
             wo(i)=0.

             timeprev=time

!            endif
         endif
        enddo

c        if(ipoute.eq.1) then
c        write(80,*) ' '
!        write(80,*) 'Total volume = ',totalvolume
c        write(80,*) 'Total mass = ',totalmassini
c        write(80,*) 'Mass remain = ',tmassremain
c        flowout=totalmassini-tmassremain
!        write(80,*) 'Total volume flow out = ',volumeout
c        write(80,*) 'Total mass flow out = ',flowout
c        write(80,*) 'Time = ',time
c        write(80,*) 'Volume flow rate =',flowout/(1000*time)
!        write(80,*) 'Volume remain = ',volumeremain
c        write(80,*) ' '
c        write(*,*) ' '
!        write(*,*) 'Total volume = ',totalvolume
c        write(*,*) 'Total mass = ',totalmassini
c        write(*,*) 'Mass remain = ',tmassremain
!        write(*,*) 'Total volume flow out  = ',volumeout
c        write(*,*) 'Total mass flow out = ',flowout
c        write(*,*) 'Time = ',time
c        write(*,*) 'Volume flow rate = ',flowout/(1000*time)
!        write(*,*) 'Volume remain = ',volumeremain
c        write(*,*) ' '

!        prevtotalmass = tmassremain
c      end if

      return
        end

