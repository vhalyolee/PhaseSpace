      DOUBLE PRECISION FUNCTION  circe_ee(zz1,zz2)
      IMPLICIT NONE
      DOUBLE PRECISION z1,z2,gamma,t1,t2,zz1,zz2,x1,x2
      DOUBLE PRECISION cir2dn
      INTEGER ierror,imap
      INTEGER init
      DATA    init /0/
      DATA    imap /1/
*------------------------------------------
      IF( init.EQ.0) THEN
         CALL cir2ld('./teslaee_500_polavg.circe',  '*', 500d0, ierror) !
         init=1
      ENDIF
      IF( imap.EQ.0) THEN
         circe_ee =  cir2dn(11, 0, -11, 0, zz1, zz2)
      ELSE
         gamma=0.4; ! this is optimal
*        gamma=0.3; ! it also OK
*        gamma=0.2; ! crazy
         gamma=0.10; ! this is optimal
         t1= zz1
         t2= zz2
*         z1= 1-exp(1/gamma*log(t1))
*         z2= 1-exp(1/gamma*log(t2))
         x1= t1**(1d0/gamma)
         x2= t2**(1d0/gamma)
         circe_ee = cir2dn(11, 0, -11, 0, 1d0-x1, 1d0-x2)
         circe_ee = circe_ee/(gamma*t1/x1);
         circe_ee = circe_ee/(gamma*t2/x2);
      ENDIF
*
      IF( ierror.LT.0) THEN
         WRITE(*,*) ' STOP in circe_init, ierror=', ierror
         STOP
      ENDIF
      END
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     circe2.f -- beam spectra for linear colliders and photon colliders
c     $Id: circe2.f,v 1.1.1.1 2001/12/14 09:58:20 jadach Exp $
c     Copyright (C) 2001 by Thorsten Ohl <ohl@hep.tu-darmstadt.de>
c     
c     Circe2 is free software; you can redistribute it and/or modify it
c     under the terms of the GNU General Public License as published by
c     the Free Software Foundation; either version 2, or (at your option)
c     any later version.
c     
c     Circe2 is distributed in the hope that it will be useful, but
c     WITHOUT ANY WARRANTY; without even the implied warranty of
c     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
c     GNU General Public License for more details.
c     
c     You should have received a copy of the GNU General Public License
c     along with this program; if not, write to the Free Software
c     Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine cir2gn (p1, h1, p2, h2, y1, y2, rng)
      implicit none
      integer p1, h1, p2, h2
      double precision y1, y2
      external rng
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      integer NBMAX, NCMAX
      parameter (NBMAX = 100, NCMAX = 36)
      integer POLAVG, POLHEL, POLGEN
      parameter (POLAVG = 1)
      parameter (POLHEL = 2)
      parameter (POLGEN = 3)
      double precision wgt(0:NBMAX*NBMAX,NCMAX)
      common /cir2cm/ wgt
      double precision val(NBMAX,NBMAX,NCMAX)
      common /cir2cm/ val
      double precision xb1(0:NBMAX,NCMAX), xb2(0:NBMAX,NCMAX)
      common /cir2cm/ xb1, xb2
      double precision lumi(NCMAX)
      common /cir2cm/ lumi
      double precision cwgt(0:NCMAX)
      common /cir2cm/ cwgt
      double precision yb1(0:NBMAX,NCMAX), yb2(0:NBMAX,NCMAX)
      common /cir2cm/ yb1, yb2
      double precision alpha1(NBMAX,NCMAX), alpha2(NBMAX,NCMAX)
      double precision xi1(NBMAX,NCMAX), xi2(NBMAX,NCMAX)
      double precision eta1(NBMAX,NCMAX), eta2(NBMAX,NCMAX)
      double precision a1(NBMAX,NCMAX), a2(NBMAX,NCMAX)
      double precision b1(NBMAX,NCMAX), b2(NBMAX,NCMAX)
      common /cir2cm/ alpha1, xi1, eta1, a1, b1
      common /cir2cm/ alpha2, xi2, eta2, a2, b2
      integer nb1(NCMAX), nb2(NCMAX)
      common /cir2cm/ nb1, nb2
      logical triang(NCMAX)
      common /cir2cm/ triang
      integer nc
      common /cir2cm/ nc
      integer pid1(NCMAX), pid2(NCMAX)
      integer pol1(NCMAX), pol2(NCMAX)
      common /cir2cm/ pid1, pol1, pid2, pol2
      integer map1(NBMAX,NCMAX), map2(NBMAX,NCMAX)
      common /cir2cm/ map1, map2
      integer polspt
      common /cir2cm/ polspt
      save /cir2cm/
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      integer i, ic, i1, i2, ibot, itop
      double precision x1, x2
      double precision u, tmp
      ic = 0
      if (((polspt .eq. POLAVG) .or. (polspt .eq. POLGEN))
     $     .and. ((h1 .ne. 0) .or. (h2 .ne. 0))) then
         write (*, '(2A)') 'circe2: current beam description ',
     $        'supports only polarization averages'
      else if ((polspt .eq. POLHEL)
     $        .and. ((h1 .eq. 0) .or. (h2 .eq. 0))) then
         write (*, '(2A)') 'circe2: polarization averages ',
     $        'not supported by current beam description'
      else
         do 10 i = 1, nc
            if (       (p1 .eq. pid1(i)) .and. (h1 .eq. pol1(i))
     $           .and. (p2 .eq. pid2(i)) .and. (h2 .eq. pol2(i))) then
               ic = i
            end if
 10      continue
      end if
      if (ic .le. 0) then
         write (*, '(A,2I4,A,2I3)')
     $        'circe2: no channel for particles', p1, p2,
     $        ' and polarizations', h1, h2
         y1 = -3.4E+38
         y2 = -3.4E+38
         return
      end if
      call rng (u)
      ibot = 0
      itop = nb1(ic) * nb2(ic)
 20   continue
      if (itop .le. (ibot + 1)) then
         i = ibot + 1
      else
         i = (ibot + itop) / 2
         if (u .lt. wgt(i,ic)) then
            itop = i
         else
            ibot = i
         end if
         goto 20
      end if
      i2 = 1 + (i - 1) / nb1(ic)
      i1 = i - (i2 - 1) * nb1(ic)
      call rng (u)
      x1 = xb1(i1,ic)*u + xb1(i1-1,ic)*(1-u)
      call rng (u)
      x2 = xb2(i2,ic)*u + xb2(i2-1,ic)*(1-u)
      if (map1(i1,ic) .eq. 0) then
         y1 = x1
      else if (map1(i1,ic) .eq. 1) then
         y1 = (a1(i1,ic)*(x1-xi1(i1,ic)))**alpha1(i1,ic) / b1(i1,ic)
     $        + eta1(i1,ic)
      else if (map1(i1,ic) .eq. 2) then
         y1 = a1(i1,ic) * tan(a1(i1,ic)*(x1-xi1(i1,ic))/b1(i1,ic)**2)
     $        + eta1(i1,ic)
      else
         write (*, '(A,I3)')
     $        'circe2: internal error: invalid map: ', map1(i1,ic)
      end if
      if (map2(i2,ic) .eq. 0) then
         y2 = x2
      else if (map2(i2,ic) .eq. 1) then
         y2 = (a2(i2,ic)*(x2-xi2(i2,ic)))**alpha2(i2,ic) / b2(i2,ic)
     $        + eta2(i2,ic)
      else if (map2(i2,ic) .eq. 2) then
         y2 = a2(i2,ic) * tan(a2(i2,ic)*(x2-xi2(i2,ic))/b2(i2,ic)**2)
     $        + eta2(i2,ic)
      else
         write (*, '(A,I3)')
     $        'circe2: internal error: invalid map: ', map2(i2,ic)
      end if
      if (triang(ic)) then
         y2 = y1 * y2
         call rng (u)
         if (2*u .ge. 1) then
            tmp = y1
            y1 = y2
            y2 = tmp
         end if
      end if
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine cir2ch (p1, h1, p2, h2, rng)
      implicit none
      integer p1, h1, p2, h2
      external rng
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      integer NBMAX, NCMAX
      parameter (NBMAX = 100, NCMAX = 36)
      integer POLAVG, POLHEL, POLGEN
      parameter (POLAVG = 1)
      parameter (POLHEL = 2)
      parameter (POLGEN = 3)
      double precision wgt(0:NBMAX*NBMAX,NCMAX)
      common /cir2cm/ wgt
      double precision val(NBMAX,NBMAX,NCMAX)
      common /cir2cm/ val
      double precision xb1(0:NBMAX,NCMAX), xb2(0:NBMAX,NCMAX)
      common /cir2cm/ xb1, xb2
      double precision lumi(NCMAX)
      common /cir2cm/ lumi
      double precision cwgt(0:NCMAX)
      common /cir2cm/ cwgt
      double precision yb1(0:NBMAX,NCMAX), yb2(0:NBMAX,NCMAX)
      common /cir2cm/ yb1, yb2
      double precision alpha1(NBMAX,NCMAX), alpha2(NBMAX,NCMAX)
      double precision xi1(NBMAX,NCMAX), xi2(NBMAX,NCMAX)
      double precision eta1(NBMAX,NCMAX), eta2(NBMAX,NCMAX)
      double precision a1(NBMAX,NCMAX), a2(NBMAX,NCMAX)
      double precision b1(NBMAX,NCMAX), b2(NBMAX,NCMAX)
      common /cir2cm/ alpha1, xi1, eta1, a1, b1
      common /cir2cm/ alpha2, xi2, eta2, a2, b2
      integer nb1(NCMAX), nb2(NCMAX)
      common /cir2cm/ nb1, nb2
      logical triang(NCMAX)
      common /cir2cm/ triang
      integer nc
      common /cir2cm/ nc
      integer pid1(NCMAX), pid2(NCMAX)
      integer pol1(NCMAX), pol2(NCMAX)
      common /cir2cm/ pid1, pol1, pid2, pol2
      integer map1(NBMAX,NCMAX), map2(NBMAX,NCMAX)
      common /cir2cm/ map1, map2
      integer polspt
      common /cir2cm/ polspt
      save /cir2cm/
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      integer ic, ibot, itop
      double precision u
      call rng (u)
      ibot = 0
      itop = nc
 10   continue
      if (itop .le. (ibot + 1)) then
         ic = ibot + 1
         p1 = pid1(ic)
         h1 = pol1(ic)
         p2 = pid2(ic)
         h2 = pol2(ic)
         return
      else
         ic = (ibot + itop) / 2
         if (u .lt. cwgt(ic)) then
            itop = ic
         else
            ibot = ic
         end if
         goto 10
      end if
      write (*, '(A)') 'circe2: internal error'
      stop
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine cir2gp (p1, p2, x1, x2, pol, rng)
      implicit none
      integer p1, p2
      double precision x1, x2
      double precision pol(0:3,0:3)
      external rng
      integer h1, h2, i1, i2
      double precision pol00
      call cir2ch (p1, h1, p2, h2, rng)
      call cir2gn (p1, h1, p2, h2, x1, x2, rng)
      call cir2dm (p1, p2, x1, x2, pol)
      pol00 = pol(0,0)
      do 10 i1 = 0, 4
         do 11 i2 = 0, 4
            pol(i1,i2) = pol(i1,i2) / pol00
 11      continue
 10   continue
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      double precision function cir2lm (p1, h1, p2, h2)
      implicit none
      integer p1, h1, p2, h2
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      integer NBMAX, NCMAX
      parameter (NBMAX = 100, NCMAX = 36)
      integer POLAVG, POLHEL, POLGEN
      parameter (POLAVG = 1)
      parameter (POLHEL = 2)
      parameter (POLGEN = 3)
      double precision wgt(0:NBMAX*NBMAX,NCMAX)
      common /cir2cm/ wgt
      double precision val(NBMAX,NBMAX,NCMAX)
      common /cir2cm/ val
      double precision xb1(0:NBMAX,NCMAX), xb2(0:NBMAX,NCMAX)
      common /cir2cm/ xb1, xb2
      double precision lumi(NCMAX)
      common /cir2cm/ lumi
      double precision cwgt(0:NCMAX)
      common /cir2cm/ cwgt
      double precision yb1(0:NBMAX,NCMAX), yb2(0:NBMAX,NCMAX)
      common /cir2cm/ yb1, yb2
      double precision alpha1(NBMAX,NCMAX), alpha2(NBMAX,NCMAX)
      double precision xi1(NBMAX,NCMAX), xi2(NBMAX,NCMAX)
      double precision eta1(NBMAX,NCMAX), eta2(NBMAX,NCMAX)
      double precision a1(NBMAX,NCMAX), a2(NBMAX,NCMAX)
      double precision b1(NBMAX,NCMAX), b2(NBMAX,NCMAX)
      common /cir2cm/ alpha1, xi1, eta1, a1, b1
      common /cir2cm/ alpha2, xi2, eta2, a2, b2
      integer nb1(NCMAX), nb2(NCMAX)
      common /cir2cm/ nb1, nb2
      logical triang(NCMAX)
      common /cir2cm/ triang
      integer nc
      common /cir2cm/ nc
      integer pid1(NCMAX), pid2(NCMAX)
      integer pol1(NCMAX), pol2(NCMAX)
      common /cir2cm/ pid1, pol1, pid2, pol2
      integer map1(NBMAX,NCMAX), map2(NBMAX,NCMAX)
      common /cir2cm/ map1, map2
      integer polspt
      common /cir2cm/ polspt
      save /cir2cm/
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      integer ic
      cir2lm = 0
      do 10 ic = 1, nc
         if (       ((p1 .eq. pid1(ic)) .or. (p1 .eq. 0))
     $        .and. ((h1 .eq. pol1(ic)) .or. (h1 .eq. 0))
     $        .and. ((p2 .eq. pid2(ic)) .or. (p2 .eq. 0))
     $        .and. ((h2 .eq. pol2(ic)) .or. (h2 .eq. 0))) then
            cir2lm = cir2lm + lumi(ic)
         end if
 10   continue
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      double precision function cir2dn (p1, h1, p2, h2, yy1, yy2)
      implicit none
      integer p1, h1, p2, h2
      double precision yy1, yy2
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      integer NBMAX, NCMAX
      parameter (NBMAX = 100, NCMAX = 36)
      integer POLAVG, POLHEL, POLGEN
      parameter (POLAVG = 1)
      parameter (POLHEL = 2)
      parameter (POLGEN = 3)
      double precision wgt(0:NBMAX*NBMAX,NCMAX)
      common /cir2cm/ wgt
      double precision val(NBMAX,NBMAX,NCMAX)
      common /cir2cm/ val
      double precision xb1(0:NBMAX,NCMAX), xb2(0:NBMAX,NCMAX)
      common /cir2cm/ xb1, xb2
      double precision lumi(NCMAX)
      common /cir2cm/ lumi
      double precision cwgt(0:NCMAX)
      common /cir2cm/ cwgt
      double precision yb1(0:NBMAX,NCMAX), yb2(0:NBMAX,NCMAX)
      common /cir2cm/ yb1, yb2
      double precision alpha1(NBMAX,NCMAX), alpha2(NBMAX,NCMAX)
      double precision xi1(NBMAX,NCMAX), xi2(NBMAX,NCMAX)
      double precision eta1(NBMAX,NCMAX), eta2(NBMAX,NCMAX)
      double precision a1(NBMAX,NCMAX), a2(NBMAX,NCMAX)
      double precision b1(NBMAX,NCMAX), b2(NBMAX,NCMAX)
      common /cir2cm/ alpha1, xi1, eta1, a1, b1
      common /cir2cm/ alpha2, xi2, eta2, a2, b2
      integer nb1(NCMAX), nb2(NCMAX)
      common /cir2cm/ nb1, nb2
      logical triang(NCMAX)
      common /cir2cm/ triang
      integer nc
      common /cir2cm/ nc
      integer pid1(NCMAX), pid2(NCMAX)
      integer pol1(NCMAX), pol2(NCMAX)
      common /cir2cm/ pid1, pol1, pid2, pol2
      integer map1(NBMAX,NCMAX), map2(NBMAX,NCMAX)
      common /cir2cm/ map1, map2
      integer polspt
      common /cir2cm/ polspt
      save /cir2cm/
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      double precision y1, y2
      integer i, ic, i1, i2, ibot, itop
      ic = 0
      if (((polspt .eq. POLAVG) .or. (polspt .eq. POLGEN))
     $     .and. ((h1 .ne. 0) .or. (h2 .ne. 0))) then
         write (*, '(2A)') 'circe2: current beam description ',
     $        'supports only polarization averages'
      else if ((polspt .eq. POLHEL)
     $        .and. ((h1 .eq. 0) .or. (h2 .eq. 0))) then
         write (*, '(2A)') 'circe2: polarization averages ',
     $        'not supported by current beam description'
      else
         do 10 i = 1, nc
            if (       (p1 .eq. pid1(i)) .and. (h1 .eq. pol1(i))
     $           .and. (p2 .eq. pid2(i)) .and. (h2 .eq. pol2(i))) then
               ic = i
            end if
 10      continue
      end if
      if (ic .le. 0) then
         cir2dn = 0
         return
      end if
      if (triang(ic)) then
         y1 = max (yy1, yy2)
         y2 = min (yy1, yy2) / y1
      else
         y1 = yy1
         y2 = yy2
      end if
      if (     (y1 .lt. yb1(0,ic)) .or. (y1 .gt. yb1(nb1(ic),ic))
     $     .or. (y2 .lt. yb2(0,ic)) .or. (y2 .gt. yb2(nb2(ic),ic))) then
         cir2dn = 0
         return
      end if
      ibot = 0
      itop = nb1(ic)
 20   continue
      if (itop .le. (ibot + 1)) then
         i1 = ibot + 1
      else
         i1 = (ibot + itop) / 2
         if (y1 .lt. yb1(i1,ic)) then
            itop = i1
         else
            ibot = i1
         end if
         goto 20
      end if
      ibot = 0
      itop = nb2(ic)
 30   continue
      if (itop .le. (ibot + 1)) then
         i2 = ibot + 1
      else
         i2 = (ibot + itop) / 2
         if (y2 .lt. yb2(i2,ic)) then
            itop = i2
         else
            ibot = i2
         end if
         goto 30
      end if
      cir2dn = val(i1,i2,ic)
      if (map1(i1,ic) .eq. 0) then
      else if (map1(i1,ic) .eq. 1) then
         cir2dn = cir2dn * b1(i1,ic) / (a1(i1,ic)*alpha1(i1,ic))
     $        * (b1(i1,ic)*(y1-eta1(i1,ic)))**(1/alpha1(i1,ic)-1)
      else if (map1(i1,ic) .eq. 2) then
         cir2dn = cir2dn * b1(i1,ic)**2
     $        / ((y1-eta1(i1,ic))**2 + a1(i1,ic)**2)
      else
         write (*, '(A,I3)')
     $        'circe2: internal error: invalid map: ', map1(i1,ic)
         stop
      end if
      if (map2(i2,ic) .eq. 0) then
      else if (map2(i2,ic) .eq. 1) then
         cir2dn = cir2dn * b2(i2,ic) / (a2(i2,ic)*alpha2(i2,ic))
     $        * (b2(i2,ic)*(y2-eta2(i2,ic)))**(1/alpha2(i2,ic)-1)
      else if (map2(i2,ic) .eq. 2) then
         cir2dn = cir2dn * b2(i2,ic)**2
     $        / ((y2-eta2(i2,ic))**2 + a2(i2,ic)**2)
      else
         write (*, '(A,I3)')
     $        'circe2: internal error: invalid map: ', map2(i2,ic)
         stop
      end if
      if (triang(ic)) then
         cir2dn = cir2dn / y1
      end if
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine cir2dm (p1, p2, x1, x2, pol)
      implicit none
      integer p1, p2
      double precision x1, x2
      double precision pol(0:3,0:3)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      integer NBMAX, NCMAX
      parameter (NBMAX = 100, NCMAX = 36)
      integer POLAVG, POLHEL, POLGEN
      parameter (POLAVG = 1)
      parameter (POLHEL = 2)
      parameter (POLGEN = 3)
      double precision wgt(0:NBMAX*NBMAX,NCMAX)
      common /cir2cm/ wgt
      double precision val(NBMAX,NBMAX,NCMAX)
      common /cir2cm/ val
      double precision xb1(0:NBMAX,NCMAX), xb2(0:NBMAX,NCMAX)
      common /cir2cm/ xb1, xb2
      double precision lumi(NCMAX)
      common /cir2cm/ lumi
      double precision cwgt(0:NCMAX)
      common /cir2cm/ cwgt
      double precision yb1(0:NBMAX,NCMAX), yb2(0:NBMAX,NCMAX)
      common /cir2cm/ yb1, yb2
      double precision alpha1(NBMAX,NCMAX), alpha2(NBMAX,NCMAX)
      double precision xi1(NBMAX,NCMAX), xi2(NBMAX,NCMAX)
      double precision eta1(NBMAX,NCMAX), eta2(NBMAX,NCMAX)
      double precision a1(NBMAX,NCMAX), a2(NBMAX,NCMAX)
      double precision b1(NBMAX,NCMAX), b2(NBMAX,NCMAX)
      common /cir2cm/ alpha1, xi1, eta1, a1, b1
      common /cir2cm/ alpha2, xi2, eta2, a2, b2
      integer nb1(NCMAX), nb2(NCMAX)
      common /cir2cm/ nb1, nb2
      logical triang(NCMAX)
      common /cir2cm/ triang
      integer nc
      common /cir2cm/ nc
      integer pid1(NCMAX), pid2(NCMAX)
      integer pol1(NCMAX), pol2(NCMAX)
      common /cir2cm/ pid1, pol1, pid2, pol2
      integer map1(NBMAX,NCMAX), map2(NBMAX,NCMAX)
      common /cir2cm/ map1, map2
      integer polspt
      common /cir2cm/ polspt
      save /cir2cm/
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      if (polspt .ne. POLGEN) then
         write (*, '(2A)') 'circe2: current beam ',
     $        'description supports no density matrices'
         return
      end if
      print *, 'circe2: cir2dm not implemented yet!'
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine cir2ld (file, design, roots, ierror)
      implicit none
      character*(*) file, design
      double precision roots
      integer ierror
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      integer NBMAX, NCMAX
      parameter (NBMAX = 100, NCMAX = 36)
      integer POLAVG, POLHEL, POLGEN
      parameter (POLAVG = 1)
      parameter (POLHEL = 2)
      parameter (POLGEN = 3)
      double precision wgt(0:NBMAX*NBMAX,NCMAX)
      common /cir2cm/ wgt
      double precision val(NBMAX,NBMAX,NCMAX)
      common /cir2cm/ val
      double precision xb1(0:NBMAX,NCMAX), xb2(0:NBMAX,NCMAX)
      common /cir2cm/ xb1, xb2
      double precision lumi(NCMAX)
      common /cir2cm/ lumi
      double precision cwgt(0:NCMAX)
      common /cir2cm/ cwgt
      double precision yb1(0:NBMAX,NCMAX), yb2(0:NBMAX,NCMAX)
      common /cir2cm/ yb1, yb2
      double precision alpha1(NBMAX,NCMAX), alpha2(NBMAX,NCMAX)
      double precision xi1(NBMAX,NCMAX), xi2(NBMAX,NCMAX)
      double precision eta1(NBMAX,NCMAX), eta2(NBMAX,NCMAX)
      double precision a1(NBMAX,NCMAX), a2(NBMAX,NCMAX)
      double precision b1(NBMAX,NCMAX), b2(NBMAX,NCMAX)
      common /cir2cm/ alpha1, xi1, eta1, a1, b1
      common /cir2cm/ alpha2, xi2, eta2, a2, b2
      integer nb1(NCMAX), nb2(NCMAX)
      common /cir2cm/ nb1, nb2
      logical triang(NCMAX)
      common /cir2cm/ triang
      integer nc
      common /cir2cm/ nc
      integer pid1(NCMAX), pid2(NCMAX)
      integer pol1(NCMAX), pol2(NCMAX)
      common /cir2cm/ pid1, pol1, pid2, pol2
      integer map1(NBMAX,NCMAX), map2(NBMAX,NCMAX)
      common /cir2cm/ map1, map2
      integer polspt
      common /cir2cm/ polspt
      save /cir2cm/
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      character*(72) buffer
      character*(72) fdesgn
      character*(72) fpolsp
      double precision froots
      integer lun, loaded, prefix
      logical match
      integer i, ic
      integer i1, i2
      double precision w
      integer status
      logical exists, isopen
      integer EOK, EFILE, EMATCH, EFORMT, ESIZE
      parameter (EOK = 0)
      parameter (EFILE = -1)
      parameter (EMATCH = -2)
      parameter (EFORMT = -3)
      parameter (ESIZE = -4)
      do 10 lun = 10, 99
         inquire (unit = lun, exist = exists,
     $        opened = isopen, iostat = status)
         if ((status .eq. 0) .and. exists .and. .not.isopen) then
            goto 11
         end if
 10   continue
      write (*, '(A)') 'cir2ld: no free unit'
      ierror = ESIZE
      return
 11   continue
      loaded = 0
      open (unit = lun, file = file, status = 'old', iostat = status)
      if (status .ne. 0) then
         write (*, '(2A)') 'cir2ld: can''t open ', file
         ierror = EFILE
         return
      end if
      if (ierror .gt. 0) then
         write (*, '(2A)') 'cir2ld: ',
     $        '$Id: circe2.f,v 1.1.1.1 2001/12/14 09:58:20 jadach Exp $'
      end if
      prefix = index (design, '*') - 1
 100  continue
 20   continue
      read (lun, '(A)', end = 29) buffer
      if (buffer(1:6) .eq. 'CIRCE2') then
         goto 21
      else if (buffer(1:1) .eq. '!') then
         if (ierror .gt. 0) then
            write (*, '(A)') buffer
         end if
         goto 20
      end if
      write (*, '(A)') 'cir2ld: invalid file'
      ierror = EFORMT
      return
 29   continue
      if (loaded .gt. 0) then           
         close (unit = lun)
         ierror = EOK
      else
         ierror = EMATCH
      end if
      return
 21   continue
      if (buffer(8:15) .eq. 'FORMAT#1') then
         read (lun, *)
         read (lun, *) fdesgn, froots
         match = .false.
         if (fdesgn .eq. design) then
            match = .true.
         else if (prefix .eq. 0) then
            match = .true.
         else if (prefix .gt. 0) then
            if (fdesgn(1:min(prefix,len(fdesgn)))
     $           .eq. design(1:min(prefix,len(design)))) then
               match = .true.
            end if
         end if
         if (match .and. (abs (froots - roots) .le. 1d0)) then
            read (lun, *) 
            read (lun, *) nc, fpolsp
            if (nc .gt. NCMAX) then
               write (*, '(A)') 'cir2ld: too many channels'
               ierror = ESIZE
               return
            end if
            if (      (fpolsp(1:1).eq.'a')
     $           .or. (fpolsp(1:1).eq.'A')) then
               polspt = POLAVG
            else if (      (fpolsp(1:1).eq.'h')
     $              .or. (fpolsp(1:1).eq.'H')) then
               polspt = POLHEL
            else if (      (fpolsp(1:1).eq.'d')
     $              .or. (fpolsp(1:1).eq.'D')) then
               polspt = POLGEN
            else
               write (*, '(A,I5)')
     $              'cir2ld: invalid polarization support: ', fpolsp
               ierror = EFORMT
               return
            end if
            cwgt(0) = 0
            do 30 ic = 1, nc
               read (lun, *)
               read (lun, *)
     $              pid1(ic), pol1(ic), pid2(ic), pol2(ic), lumi(ic)
               cwgt(ic) = cwgt(ic-1) + lumi(ic)
               if (polspt .eq. POLAVG
     $              .and. (      (pol1(ic) .ne. 0)
     $              .or. (pol2(ic) .ne. 0))) then
                  write (*, '(A)')
     $                 'cir2ld: expecting averaged polarization'
                  ierror = EFORMT
                  return
               else if (polspt .eq. POLHEL
     $                 .and. (      (pol1(ic) .eq. 0)
     $                 .or. (pol2(ic) .eq. 0))) then
                  write (*, '(A)')
     $                 'cir2ld: expecting helicities'
                  ierror = EFORMT
                  return
               else if (polspt .eq. POLGEN) then
                  write (*, '(A)')
     $                 'cir2ld: general polarizations not supported yet'
                  ierror = EFORMT
                  return
               else if (polspt .eq. POLGEN
     $                 .and. (      (pol1(ic) .ne. 0)
     $                 .or. (pol2(ic) .ne. 0))) then
                  write (*, '(A)') 'cir2ld: expecting pol = 0'
                  ierror = EFORMT
                  return
               end if
               read (lun, *)
               read (lun, *) nb1(ic), nb2(ic), triang(ic)
               if ((nb1(ic) .gt. NBMAX) .or. (nb2(ic) .gt. NBMAX)) then
                  write (*, '(A)') 'cir2ld: too many bins'
                  ierror = ESIZE
                  return
               end if
               read (lun, *)
               read (lun, *) xb1(0,ic)
               do 31 i1 = 1, nb1(ic)
                  read (lun, *) xb1(i1,ic), map1(i1,ic), alpha1(i1,ic),
     $                 xi1(i1,ic), eta1(i1,ic), a1(i1,ic), b1(i1,ic)
 31            continue
               read (lun, *)
               read (lun, *) xb2(0,ic)
               do 32 i2 = 1, nb2(ic)
                  read (lun, *) xb2(i2,ic), map2(i2,ic), alpha2(i2,ic),
     $                 xi2(i2,ic), eta2(i2,ic), a2(i2,ic), b2(i2,ic)
 32            continue
               do 33 i = 0, nb1(ic)
                  i1 = max (i, 1)
                  if (map1(i1,ic) .eq. 0) then
                     yb1(i,ic) = xb1(i,ic)
                  else if (map1(i1,ic) .eq. 1) then
                     yb1(i,ic) =
     $                    (a1(i1,ic)
     $                    * (xb1(i,ic)-xi1(i1,ic)))**alpha1(i1,ic)
     $                    / b1(i1,ic) + eta1(i1,ic)
                  else if (map1(i1,ic) .eq. 2) then
                     yb1(i,ic) = a1(i1,ic)
     $                    * tan(a1(i1,ic)/b1(i1,ic)**2
     $                    * (xb1(i,ic)-xi1(i1,ic)))
     $                    + eta1(i1,ic)
                  else
                     write (*, '(A,I3)')
     $                    'cir2ld: invalid map: ', map1(i1,ic)
                     ierror = EFORMT
                     return
                  end if
 33            continue
               do 34 i = 0, nb2(ic)
                  i2 = max (i, 1)
                  if (map2(i2,ic) .eq. 0) then
                     yb2(i,ic) = xb2(i,ic)
                  else if (map2(i2,ic) .eq. 1) then
                     yb2(i,ic)
     $                    = (a2(i2,ic)
     $                    * (xb2(i,ic)-xi2(i2,ic)))**alpha2(i2,ic)
     $                    / b2(i2,ic) + eta2(i2,ic)
                  else if (map2(i2,ic) .eq. 2) then
                     yb2(i,ic) = a2(i2,ic)
     $                    * tan(a2(i2,ic)/b2(i2,ic)**2
     $                    * (xb2(i,ic)-xi2(i2,ic)))
     $                    + eta2(i2,ic)
                  else
                     write (*, '(A,I3)')
     $                    'cir2ld: invalid map: ', map2(i2,ic)
                     ierror = EFORMT
                     return
                  end if
 34            continue
               read (lun, *)
               wgt(0,ic) = 0
               do 35 i = 1, nb1(ic)*nb2(ic)
                  read (lun, *) w
                  wgt(i,ic) = wgt(i-1,ic) + w
                  i2 = 1 + (i - 1) / nb1(ic)
                  i1 = i - (i2 - 1) * nb1(ic)
                  val(i1,i2,ic) = w
     $                 / (  (xb1(i1,ic) - xb1(i1-1,ic))
     $                 * (xb2(i2,ic) - xb2(i2-1,ic)))
 35            continue
               wgt(nb1(ic)*nb2(ic),ic) = 1
 30         continue
            do 40 ic = 1, nc
               cwgt(ic) = cwgt(ic) / cwgt(nc)
 40         continue
            loaded = loaded + 1
         else
 101        continue
            read (lun, *) buffer
            if (buffer(1:6) .ne. 'ECRIC2') then
               goto 101
            end if
            goto 100
         end if
      else
         write (*, '(2A)') 'cir2ld: invalid format: ', buffer(8:72)
         ierror = EFORMT
         return
      end if
      read (lun, '(A)') buffer
      if (buffer(1:6) .ne. 'ECRIC2') then
         write (*, '(A)') 'cir2ld: invalid file'
         ierror = EFORMT
         return
      end if
      goto 100
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
