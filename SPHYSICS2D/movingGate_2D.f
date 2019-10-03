c    "Copyright 2009 Prof. Robert Dalrymple, Prof. M. Gomez Gesteira, Dr Benedict Rogers, 
c     Dr Alejandro Crespo, Dr Muthukumar Narayanaswamy, Dr Shan Zou, Dr Andrea Panizzo "
c
c    This file is part of SPHYSICS.
c
c    SPHYSICS is free software; you can redistribute it and/or modify
c    it under the terms of the GNU General Public License as published by
c    the Free Software Foundation; either version 3 of the License, or
c    (at your option) any later version.
c
c    SPHYSICS is distributed in the hope that it will be useful,
c    but WITHOUT ANY WARRANTY; without even the implied warranty of
c    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
c    GNU General Public License for more details.
c
c    You should have received a copy of the GNU General Public License
c    along with this program.  If not, see <http://www.gnu.org/licenses/>.

c
c*********************************************************************
c
      subroutine movingGate

      include 'common.2D'         

        do i=ngate_ini,ngate_end
          if (time.ge.tgate) then 
             
!              if ((xp(i)+VXgate*dt).gt.xmax_ini)then
!                  xp(i)= xmin
!              end if      
             
             xp(i)= xp(i)+VXgate*dt
             zp(i)= min(zp(i)+VZgate*dt,zmax_ini)

             up(i)=VXgate
             wp(i)=VZgate
          endif
        enddo
      
	return      
      end
c
c*********************************************************************
