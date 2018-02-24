      Program hplaq

* Two dimensional hexagonal plane!! 2 bands!!
* Compute Chern number C for band ibb using plaquette method!!


      implicit real*8 (a-h,o-z)
      parameter (maxp=1500,kdim=2,nq0=25000)

      double precision eA,eB,t1,t2,phi
      double precision rtemp1,rtemp2,rtemp3 


      double precision pi,twopi 
      double precision uxm,uym,uzm
      double precision kx,ky,kx0,ky0,dkx,dky
      double precision ux,uy,uz,dux,duy,duz, kmaxGK,kmaxKKP

      double precision nos(maxp),dos(maxp)
      double precision tnos(maxp),tdos(maxp),x1(2),f1(2)
      double precision emin,emax,ebot,etop,e(2,8),ee(4),bb(4)
      double precision freq(maxp),rdos(maxp),cdos(maxp),
     .                 v2(maxp),v2d(maxp),rgamma(maxp)

      double precision rsp(maxp),csp(maxp),sumr(maxp),sumc(maxp)
      double precision h1(2,2,2), o1(2,2,2), z1(2,2,2), w1(2,11), 
     .                 e1(2)

      double precision z1xp(2,2,2)  !kx+d
      double precision z1xm(2,2,2)  !kx-d
      double precision z1yp(2,2,2)  !ky+d
      double precision z1ym(2,2,2)  !ky-d
      double precision phiC(maxp) 
      complex*16 cdiff
      double precision  F(1000,1000), Ftilde(1000,1000) 
      double precision numbk, ekn(2), eqn(2), G1(2), G2(2)
      double precision gam1(2), gam2(2), 
     .                 phi12, phi23, phi34, phi41

      complex*16 b(2,2,8), xi, csum, csum1, csum2, csum3, csum4
      complex*16 csum12, csum23, csum34, csum41
      complex*16 csumAx, csumAy
      complex*16 utilde(nq0,2,2)
      complex*16 derutilde(nq0,2,2), arg, arg1
      complex*16 ukxp(2)
      complex*16 ukxm(2)
      complex*16 ukyp(2)
      complex*16 ukym(2)
      complex*16 derxuk(2)
      complex*16 deryuk(2)
      complex*16 H0, Hx, Hy, Hz 
      complex*16 ukn(2,2), uqn(2,2), A(2,2)
      complex*16 ukn1(2,2)
      complex*16 ukn2(2,2)
      complex*16 ukn3(2,2)
      complex*16 ukn4(2,2)
      complex*16 Gphase, expphi, ctemp1, ctemp2, expphiD 
      complex*16 expphiD1, expphiD2 
      complex*16 expphi0, expphi2pi,delphix(1500), delphiy(1500)
      complex*16 ca1, cb1, ca2, cb2, cnorm 1, cnorm2

      integer nkabc(3), p(4,6), iw1(2)
      integer itest 



      data p/1,2,4,5,4,5,7,8,2,4,5,7,1,3,4,5,4,5,6,8,3,4,5,6/

      pi    = 4.d0*atan(1.d0)
      twopi = 8.d0*atan(1.d0)
      xi    = dcmplx(0.d0,1.d0)

* Open files.
      open(unit=1,file='haldane2.d',status='OLD')
      open(unit=2,file='ftilde.d',status='unknown')

      read(1,*)ipp,itotal
      read(1,*)nkabc
      read(1,*)npts,emin,emax
      read(1,*)t1,t2
      read(1,*)phi
      read(1,*)dx,dy
      close(1,status='keep')
      if(npts.gt.maxp)then
      print*,' Dimension too small, increase maxp!!'
      stop
      endif
*      if(nkabc(1)*nkabc(2).gt.nq0)stop 'Increase nq0 !!'

* Control parameters
      write(6,*)'kdim:',kdim
      write(6,*)'nkabc(i)',nkabc(1),nkabc(2),nkabc(3)
      write(6,'(a20,f10.3)')'t1:',t1
      write(6,'(a20,f10.3)')'t2:',t2
      write(6,'(a20,f10.3)')'phi (units of pi):',phi
      write(6,'(a20,f10.6)')'dx:',dx
      write(6,'(a20,f10.6)')'dy:',dy

      phi = phi*pi
      write(6,'(a20,f10.3)')'phi (in radians):',phi


* Calculate Chern number C.
      volwgt = ( (8.d0*pi*pi)/(3.d0*dsqrt(3.d0)) )*
     .         ( 1.d0/((nkabc(1)+1)*(nkabc(2)+1)) )

      numbk = (nkabc(1)+1)*(nkabc(2)+1)
      write(6,*)'Number k:',numbk

      uxm=(2.d0/3.d0)*2.d0*pi
      uym=uxm
      dux=dble(uxm/nkabc(1))
      duy=dble(uym/nkabc(2))

* Choose band.
      ibb = 1
*      ibb = 2
      if(ibb.eq.1) write(6,'(a20,i3)')'Chern number band:',ibb
      if(ibb.eq.2) write(6,'(a20,i3)')'Chern number band:',ibb

* Choose yy = Delta/t2 !!
*      yy = 6.00d0
*      yy = 4.00d0
*      yy = 3.69d0
*      yy = 3.68d0
*      yy = 3.67d0
*      yy = 3.d0*dsqrt(3.d0)/dsqrt(2.d0)   !Critical value
*      yy = 3.65d0
      yy = 3.60d0
*      yy = 3.00d0
*      yy = 1.00d0
*      yy = 0.10d0

      eA = -t2*yy
      eB = -eA
      write(6,'(a20,f10.3)')'eA:',eA
      write(6,'(a20,f10.3)')'eB:',eB
      write(6,'(a20,f10.3)')'Delta/t2:',dabs(eA)/t2
      write(6,*)

      nphi =   1

      dphi = twopi/dble(nphi-1)
      do i = 1,nphi
      phiC(i) = -pi +dphi*dble(i-1)
      enddo

      do ip = 1, nphi 

*      phi = phiC(ip) 


* For Chern number!
      csum1 = dcmplx(0.d0,0.d0)

      ux=-dux

      icount = 0

      do ix = 1, nkabc(1) +1 
*      do ix = 1, 1

      ux=ux+dux
      uy=-duy

      do iy = 1, nkabc(2) +1 
*      do iy = 1, 1 

      uy=uy+duy

      icount = icount +1

* For fix (ux,uy) solve for (kx,ky).
      kx  =  ux*dsqrt(3.d0)/2.d0
      ky  = -ux/2.d0 + uy
*      write(6,'(2i4,2f10.6)')ix,iy,kx/twopi, ky/twopi



* Fix diagonal elements.
      do i = 1,kdim
      do j = 1,kdim
      h1(i,j,1)  = 0.d0
      h1(i,j,2)  = 0.d0
      o1(i,j,1)  = 0.d0
      o1(i,j,2)  = 0.d0
      enddo
      enddo

      do i = 1,kdim
      o1(i,i,1)  = 1.d0
      enddo

* Diagonal elements H11 and H22 first!
      rtemp1     = dcos(dsqrt(3.d0)*kx - phi)
      rtemp2     = dcos(-dsqrt(3.d0)*kx/2.d0 + 3.d0*ky/2.d0 - phi)
      rtemp3     = dcos(-dsqrt(3.d0)*kx/2.d0 - 3.d0*ky/2.d0 - phi)

      h1(1,1,1)  =  -2.d0*t2*(rtemp1 + rtemp2 + rtemp3) + eA

      rtemp1     = dcos(dsqrt(3.d0)*kx + phi)
      rtemp2     = dcos(-dsqrt(3.d0)*kx/2.d0 + 3.d0*ky/2.d0 + phi)
      rtemp3     = dcos(-dsqrt(3.d0)*kx/2.d0 - 3.d0*ky/2.d0 + phi)

      h1(2,2,1)  =  -2.d0*t2*(rtemp1 + rtemp2 + rtemp3) + eB

* Hopping H12.
      rtemp1     = dcos(ky)
      rtemp2     = dcos(-dsqrt(3.d0)*kx/2.d0 - ky/2.d0)
      rtemp3     = dcos( dsqrt(3.d0)*kx/2.d0 - ky/2.d0)

      h1(1,2,1)  =  -t1*(rtemp1 + rtemp2 + rtemp3)

      rtemp1     = dsin(ky)
      rtemp2     = dsin(-dsqrt(3.d0)*kx/2.d0 - ky/2.d0)
      rtemp3     = dsin( dsqrt(3.d0)*kx/2.d0 - ky/2.d0)

      h1(1,2,2)  =   t1*(rtemp1 + rtemp2 + rtemp3)

* H21.
      h1(2,1,1)  =    h1(1,2,1)
      h1(2,1,2)  =  - h1(1,2,2)

* For kx,ky!!
      call diagno(kdim,h1,o1,w1,iw1,z1,e1)

      do ib      = 1,kdim  !Band all
      do i1      = 1,kdim
      ukn(i1,ib) = dcmplx(z1(i1,ib,1),z1(i1,ib,2))
      enddo
      enddo


* For kx+dx,ky!!
      kx  =  ux*dsqrt(3.d0)/2.d0 + dx
      ky  = -ux/2.d0 + uy
* Fix diagonal elements.
      do i = 1,kdim
      do j = 1,kdim
      h1(i,j,1)  = 0.d0
      h1(i,j,2)  = 0.d0
      o1(i,j,1)  = 0.d0
      o1(i,j,2)  = 0.d0
      enddo
      enddo

      do i = 1,kdim
      o1(i,i,1)  = 1.d0
      enddo

* Diagonal elements H11 and H22 first!
      rtemp1     = dcos(dsqrt(3.d0)*kx - phi)
      rtemp2     = dcos(-dsqrt(3.d0)*kx/2.d0 + 3.d0*ky/2.d0 - phi)
      rtemp3     = dcos(-dsqrt(3.d0)*kx/2.d0 - 3.d0*ky/2.d0 - phi)

      h1(1,1,1)  =  -2.d0*t2*(rtemp1 + rtemp2 + rtemp3) + eA

      rtemp1     = dcos(dsqrt(3.d0)*kx + phi)
      rtemp2     = dcos(-dsqrt(3.d0)*kx/2.d0 + 3.d0*ky/2.d0 + phi)
      rtemp3     = dcos(-dsqrt(3.d0)*kx/2.d0 - 3.d0*ky/2.d0 + phi)

      h1(2,2,1)  =  -2.d0*t2*(rtemp1 + rtemp2 + rtemp3) + eB

* Hopping H12.
      rtemp1     = dcos(ky)
      rtemp2     = dcos(-dsqrt(3.d0)*kx/2.d0 - ky/2.d0)
      rtemp3     = dcos( dsqrt(3.d0)*kx/2.d0 - ky/2.d0)

      h1(1,2,1)  =  -t1*(rtemp1 + rtemp2 + rtemp3)

      rtemp1     = dsin(ky)
      rtemp2     = dsin(-dsqrt(3.d0)*kx/2.d0 - ky/2.d0)
      rtemp3     = dsin( dsqrt(3.d0)*kx/2.d0 - ky/2.d0)

      h1(1,2,2)  =   t1*(rtemp1 + rtemp2 + rtemp3)

* H21.
      h1(2,1,1)  =    h1(1,2,1)
      h1(2,1,2)  =  - h1(1,2,2)


* For kx+dx,ky!!
      call diagno(kdim,h1,o1,w1,iw1,z1xp,e1)

* Determine dual state! Band ibb only occupied.
      csum  = dcmplx(0.d0,0.d0)
      do i1 = 1,kdim
      csum  = csum + dcmplx(z1  (i1,ibb,1),-z1  (i1,ibb,2))*
     .               dcmplx(z1xp(i1,ibb,1), z1xp(i1,ibb,2))
      enddo 

      do i1    = 1,kdim
      ukxp(i1) = dcmplx(z1xp(i1,ibb,1),z1xp(i1,ibb,2))/csum
      enddo 


* For kx-dx,ky!!
      kx  =  ux*dsqrt(3.d0)/2.d0 - dx
      ky  = -ux/2.d0 + uy
* Fix diagonal elements.
      do i = 1,kdim
      do j = 1,kdim
      h1(i,j,1)  = 0.d0
      h1(i,j,2)  = 0.d0
      o1(i,j,1)  = 0.d0
      o1(i,j,2)  = 0.d0
      enddo
      enddo

      do i = 1,kdim
      o1(i,i,1)  = 1.d0
      enddo

* Diagonal elements H11 and H22 first!
      rtemp1     = dcos(dsqrt(3.d0)*kx - phi)
      rtemp2     = dcos(-dsqrt(3.d0)*kx/2.d0 + 3.d0*ky/2.d0 - phi)
      rtemp3     = dcos(-dsqrt(3.d0)*kx/2.d0 - 3.d0*ky/2.d0 - phi)

      h1(1,1,1)  =  -2.d0*t2*(rtemp1 + rtemp2 + rtemp3) + eA

      rtemp1     = dcos(dsqrt(3.d0)*kx + phi)
      rtemp2     = dcos(-dsqrt(3.d0)*kx/2.d0 + 3.d0*ky/2.d0 + phi)
      rtemp3     = dcos(-dsqrt(3.d0)*kx/2.d0 - 3.d0*ky/2.d0 + phi)

      h1(2,2,1)  =  -2.d0*t2*(rtemp1 + rtemp2 + rtemp3) + eB

* Hopping H12.
      rtemp1     = dcos(ky)
      rtemp2     = dcos(-dsqrt(3.d0)*kx/2.d0 - ky/2.d0)
      rtemp3     = dcos( dsqrt(3.d0)*kx/2.d0 - ky/2.d0)

      h1(1,2,1)  =  -t1*(rtemp1 + rtemp2 + rtemp3)

      rtemp1     = dsin(ky)
      rtemp2     = dsin(-dsqrt(3.d0)*kx/2.d0 - ky/2.d0)
      rtemp3     = dsin( dsqrt(3.d0)*kx/2.d0 - ky/2.d0)

      h1(1,2,2)  =   t1*(rtemp1 + rtemp2 + rtemp3)

* H21.
      h1(2,1,1)  =    h1(1,2,1)
      h1(2,1,2)  =  - h1(1,2,2)


* For kx-dx,ky!!
      call diagno(kdim,h1,o1,w1,iw1,z1xm,e1)

* Determine dual state! Band ibb only occupied.
      csum  = dcmplx(0.d0,0.d0)
      do i1 = 1,kdim
      csum  = csum + dcmplx(z1  (i1,ibb,1),-z1  (i1,ibb,2))*
     .               dcmplx(z1xm(i1,ibb,1), z1xm(i1,ibb,2))
      enddo


      do i1    = 1,kdim
      ukxm(i1) = dcmplx(z1xm(i1,ibb,1),z1xm(i1,ibb,2))/csum
      enddo


* For kx,ky+dy!!
      kx  =  ux*dsqrt(3.d0)/2.d0 
      ky  = -ux/2.d0 + uy + dy
* Fix diagonal elements.
      do i = 1,kdim
      do j = 1,kdim
      h1(i,j,1)  = 0.d0
      h1(i,j,2)  = 0.d0
      o1(i,j,1)  = 0.d0
      o1(i,j,2)  = 0.d0
      enddo
      enddo

      do i = 1,kdim
      o1(i,i,1)  = 1.d0
      enddo

* Diagonal elements H11 and H22 first!
      rtemp1     = dcos(dsqrt(3.d0)*kx - phi)
      rtemp2     = dcos(-dsqrt(3.d0)*kx/2.d0 + 3.d0*ky/2.d0 - phi)
      rtemp3     = dcos(-dsqrt(3.d0)*kx/2.d0 - 3.d0*ky/2.d0 - phi)

      h1(1,1,1)  =  -2.d0*t2*(rtemp1 + rtemp2 + rtemp3) + eA

      rtemp1     = dcos(dsqrt(3.d0)*kx + phi)
      rtemp2     = dcos(-dsqrt(3.d0)*kx/2.d0 + 3.d0*ky/2.d0 + phi)
      rtemp3     = dcos(-dsqrt(3.d0)*kx/2.d0 - 3.d0*ky/2.d0 + phi)

      h1(2,2,1)  =  -2.d0*t2*(rtemp1 + rtemp2 + rtemp3) + eB

* Hopping H12.
      rtemp1     = dcos(ky)
      rtemp2     = dcos(-dsqrt(3.d0)*kx/2.d0 - ky/2.d0)
      rtemp3     = dcos( dsqrt(3.d0)*kx/2.d0 - ky/2.d0)

      h1(1,2,1)  =  -t1*(rtemp1 + rtemp2 + rtemp3)

      rtemp1     = dsin(ky)
      rtemp2     = dsin(-dsqrt(3.d0)*kx/2.d0 - ky/2.d0)
      rtemp3     = dsin( dsqrt(3.d0)*kx/2.d0 - ky/2.d0)

      h1(1,2,2)  =   t1*(rtemp1 + rtemp2 + rtemp3)

* H21.
      h1(2,1,1)  =    h1(1,2,1)
      h1(2,1,2)  =  - h1(1,2,2)


* For kx,ky+dy!!
      call diagno(kdim,h1,o1,w1,iw1,z1yp,e1)

* Determine dual state! Band ibb only occupied.
      csum  = dcmplx(0.d0,0.d0)
      do i1 = 1,kdim
      csum  = csum + dcmplx(z1  (i1,ibb,1),-z1  (i1,ibb,2))*
     .               dcmplx(z1yp(i1,ibb,1), z1yp(i1,ibb,2))
      enddo


      do i1    = 1,kdim
      ukyp(i1) = dcmplx(z1yp(i1,ibb,1),z1yp(i1,ibb,2))/csum
      enddo

* For kx,ky-dy!!
      kx  =  ux*dsqrt(3.d0)/2.d0
      ky  = -ux/2.d0 + uy - dy
* Fix diagonal elements.
      do i = 1,kdim
      do j = 1,kdim
      h1(i,j,1)  = 0.d0
      h1(i,j,2)  = 0.d0
      o1(i,j,1)  = 0.d0
      o1(i,j,2)  = 0.d0
      enddo
      enddo

      do i = 1,kdim
      o1(i,i,1)  = 1.d0
      enddo

* Diagonal elements H11 and H22 first!
      rtemp1     = dcos(dsqrt(3.d0)*kx - phi)
      rtemp2     = dcos(-dsqrt(3.d0)*kx/2.d0 + 3.d0*ky/2.d0 - phi)
      rtemp3     = dcos(-dsqrt(3.d0)*kx/2.d0 - 3.d0*ky/2.d0 - phi)

      h1(1,1,1)  =  -2.d0*t2*(rtemp1 + rtemp2 + rtemp3) + eA

      rtemp1     = dcos(dsqrt(3.d0)*kx + phi)
      rtemp2     = dcos(-dsqrt(3.d0)*kx/2.d0 + 3.d0*ky/2.d0 + phi)
      rtemp3     = dcos(-dsqrt(3.d0)*kx/2.d0 - 3.d0*ky/2.d0 + phi)

      h1(2,2,1)  =  -2.d0*t2*(rtemp1 + rtemp2 + rtemp3) + eB

* Hopping H12.
      rtemp1     = dcos(ky)
      rtemp2     = dcos(-dsqrt(3.d0)*kx/2.d0 - ky/2.d0)
      rtemp3     = dcos( dsqrt(3.d0)*kx/2.d0 - ky/2.d0)
      
      h1(1,2,1)  =  -t1*(rtemp1 + rtemp2 + rtemp3)

      rtemp1     = dsin(ky)
      rtemp2     = dsin(-dsqrt(3.d0)*kx/2.d0 - ky/2.d0)
      rtemp3     = dsin( dsqrt(3.d0)*kx/2.d0 - ky/2.d0)

      h1(1,2,2)  =   t1*(rtemp1 + rtemp2 + rtemp3)

* Hopping H12.
      rtemp1     = dcos(ky)
      rtemp2     = dcos(-dsqrt(3.d0)*kx/2.d0 - ky/2.d0)
      rtemp3     = dcos( dsqrt(3.d0)*kx/2.d0 - ky/2.d0)

      h1(1,2,1)  =  -t1*(rtemp1 + rtemp2 + rtemp3)

      rtemp1     = dsin(ky)
      rtemp2     = dsin(-dsqrt(3.d0)*kx/2.d0 - ky/2.d0)
      rtemp3     = dsin( dsqrt(3.d0)*kx/2.d0 - ky/2.d0)

      h1(1,2,2)  =   t1*(rtemp1 + rtemp2 + rtemp3)

* H21.
      h1(2,1,1)  =    h1(1,2,1)
      h1(2,1,2)  =  - h1(1,2,2)


* For kx,ky-dy!!
      call diagno(kdim,h1,o1,w1,iw1,z1ym,e1)

* Determine dual state! Band ibb only occupied.
      csum  = dcmplx(0.d0,0.d0)
      do i1 = 1,kdim
      csum  = csum + dcmplx(z1  (i1,ibb,1),-z1  (i1,ibb,2))*
     .               dcmplx(z1ym(i1,ibb,1), z1ym(i1,ibb,2))
      enddo


      do i1    = 1,kdim
      ukym(i1) = dcmplx(z1ym(i1,ibb,1),z1ym(i1,ibb,2))/csum
      enddo


* Determine the derivatives!
      do i2      = 1,kdim
      derxuk(i2) = ( ukxp(i2) - ukxm(i2) )/(2.d0*dx) 
      deryuk(i2) = ( ukyp(i2) - ukym(i2) )/(2.d0*dy) 
      enddo

*      csumAx  = dcmplx(0.d0,0.d0)
*      csumAy  = dcmplx(0.d0,0.d0)
*      do i2      = 1,kdim
*      csumAx = csumAx + derxuk(i2)*dconjg(ukn(i2,1)) 
*      csumAy = csumAy + deryuk(i2)*dconjg(ukn(i2,1)) 
*      csumAx = csumAx + derxuk(i2)*dconjg(ukn(i2,1)) 
*      csumAy = csumAy + deryuk(i2)*dconjg(ukn(i2,1)) 
*      csumAx = csumAx + ukxp(i2)*dconjg(ukn(i2,1)) 
*      csumAy = csumAy + ukyp(i2)*dconjg(ukn(i2,1)) 
*      enddo
*      write(6,'(a10,2f8.4,8x,6f14.8)')'Ax Ay ukn:',kx/twopi,ky/twopi,
*     .           csumAx, csumAy 
*      if(cdabs(csumAx).gt.0.00000001d0)then
*      write(6,'(a10,2f8.4,8x,6f10.8)')'Ax:',kx/twopi,ky/twopi,
*     .           csumAx
*      endif
*      if(cdabs(csumAy).gt.0.00000001d0)then
*      write(6,'(a10,2f8.4,8x,6f10.8)')'Ay ukn:',kx/twopi,ky/twopi,
*     .           csumAy
*      endif

*      write(6,*)'derxuk(1):', derxuk(1)*2.d0*dx
*      write(6,*)'derxuk(2):', derxuk(2)*2.d0*dx
*      write(6,*)'deryuk(1):', deryuk(1)*2.d0*dy
*      write(6,*)'deryuk(2):', deryuk(2)*2.d0*dy

      csum1 = csum1  +  dconjg(derxuk(1))*deryuk(1)
     .               +  dconjg(derxuk(2))*deryuk(2) 
     .               -  dconjg(deryuk(1))*derxuk(1) 
     .               -  dconjg(deryuk(2))*derxuk(2) 

      ctemp1 =  dconjg(derxuk(1))*deryuk(1)
     .       +  dconjg(derxuk(2))*deryuk(2)
     .       -  dconjg(deryuk(1))*derxuk(1)
     .       -  dconjg(deryuk(2))*derxuk(2)


      enddo
      enddo     !End k-sums

      write(6,*)
      write(6,*)'icount:',icount
      write(6,*)'Number k:',numbk
      write(6,'(a30,f10.2,4x,2f20.4)')'Chern number (sum over BZ):',
     .                                 phi/pi,
     .                                 csum1*volwgt*xi/twopi

      enddo     !End phi-loop
*      stop





      write(6,*)
      write(6,*)
      write(6,'(a40,2f12.4)')'Chern number (summing plaquettes):'
      write(6,*)
      write(6,*)

      do ibb1 = 1,2  !Band loop
*      do ibb1 = 1,1  !Band loop
*      do ibb1 = 2,2  !Band loop

      write(6,*)
      write(6,*)'Band:',ibb1
      write(6,*)

      rsumF      = 0.d0
      rsumFtilde = 0.d0
      ux=-dux
* Sum plaquettes !!

      icount = 0
* First all states !!
      do ix = 1, nkabc(1)+1      !No boundary

      ux=ux+dux
      uy=-duy

      do iy = 1, nkabc(2)+1       !No boundary

      uy=uy+duy

      icount = icount + 1

* For fix (ux,uy) solve for (kx,ky).
* Left corner label 1.
      kx  =  ux*dsqrt(3.d0)/2.d0
      ky  = -ux/2.d0 + uy

*      write(6,'(a20,2f10.4)')'u v:',ux/twopi, uy/twopi
*      write(6,'(a20,2f10.4)')'1:',kx/twopi,ky/twopi

* Fix diagonal elements.
      do i = 1,kdim
      do j = 1,kdim
      h1(i,j,1)  = 0.d0
      h1(i,j,2)  = 0.d0
      o1(i,j,1)  = 0.d0
      o1(i,j,2)  = 0.d0
      enddo
      enddo

      do i = 1,kdim
      o1(i,i,1)  = 1.d0
      enddo


* Diagonal elements H11 and H22 first!
      rtemp1     = dcos(dsqrt(3.d0)*kx - phi)
      rtemp2     = dcos(-dsqrt(3.d0)*kx/2.d0 + 3.d0*ky/2.d0 - phi)
      rtemp3     = dcos(-dsqrt(3.d0)*kx/2.d0 - 3.d0*ky/2.d0 - phi)

      h1(1,1,1)  =  -2.d0*t2*(rtemp1 + rtemp2 + rtemp3) + eA

      rtemp1     = dcos(dsqrt(3.d0)*kx + phi)
      rtemp2     = dcos(-dsqrt(3.d0)*kx/2.d0 + 3.d0*ky/2.d0 + phi)
      rtemp3     = dcos(-dsqrt(3.d0)*kx/2.d0 - 3.d0*ky/2.d0 + phi)

      h1(2,2,1)  =  -2.d0*t2*(rtemp1 + rtemp2 + rtemp3) + eB

* Hopping H12.
      rtemp1     = dcos(ky)
      rtemp2     = dcos(-dsqrt(3.d0)*kx/2.d0 - ky/2.d0)
      rtemp3     = dcos( dsqrt(3.d0)*kx/2.d0 - ky/2.d0)

      h1(1,2,1)  =  -t1*(rtemp1 + rtemp2 + rtemp3)

      rtemp1     = dsin(ky)
      rtemp2     = dsin(-dsqrt(3.d0)*kx/2.d0 - ky/2.d0)
      rtemp3     = dsin( dsqrt(3.d0)*kx/2.d0 - ky/2.d0)

      h1(1,2,2)  =   t1*(rtemp1 + rtemp2 + rtemp3)

* H21.
      h1(2,1,1)  =   h1(1,2,1)
      h1(2,1,2)  =  -h1(1,2,2)

* For kx,ky!!
      call diagno(kdim,h1,o1,w1,iw1,z1,e1)

      do ib      = 1,kdim  !Band
      do i1      = 1,kdim
      utilde(icount,i1,ib) = dcmplx(z1(i1,ib,1),z1(i1,ib,2)) 
      enddo
      enddo

      enddo  !All states
      enddo

      write(6,*)'All states:',icount
*      stop

      if(nkabc(1).gt.1000) stop 'Increase dimension Ftilde and F!!'

      ux=-dux

      do ix = 1, nkabc(1)       !No boundary

      ux=ux+dux
      uy=-duy

      do iy = 1, nkabc(1)       !No boundary

      uy=uy+duy

      icount = icount +1

      kx  =  ux*dsqrt(3.d0)/2.d0
      ky  = -ux/2.d0 + uy


* For this plaquette!!
* Band ibb.
      csum1  = dcmplx(0.d0,0.d0)
      ip1    = 1 + (nkabc(1)+1)*(ix-1) + (iy-1)
      ip2    = ip1 + (nkabc(1) + 1)  
      do i1  = 1,kdim
      csum1  = csum1 + dconjg(utilde(ip1,i1,ibb1))*utilde(ip2,i1,ibb1)   !1 to 2
      enddo
      phi12  = cdlog(csum1/cdabs(csum1))/xi
*      write(6,'(a30,2i4,3f10.4)')'phi12:',ip1,ip2,phi12

      csum2  = dcmplx(0.d0,0.d0)
      ip1    = ip2
      ip2    = ip2 + 1 
      do i1  = 1,kdim
      csum2  = csum2 + dconjg(utilde(ip1,i1,ibb1))*utilde(ip2,i1,ibb1)   !2 to 3
      enddo
      phi23  = cdlog(csum2/cdabs(csum2))/xi
*      write(6,'(a30,2i4,3f10.4)')'phi23:',ip1,ip2,phi23

      csum3  = dcmplx(0.d0,0.d0)
      ip1    = ip2
      ip2    = ip2 - (nkabc(1)+1)
      do i1  = 1,kdim
      csum3  = csum3 + dconjg(utilde(ip1,i1,ibb1))*utilde(ip2,i1,ibb1)   !3 to 4
      enddo
      phi34  = cdlog(csum3/cdabs(csum3))/xi
*      write(6,'(a30,2i4,3f10.4)')'phi34:',ip1,ip2,phi34

      csum4  = dcmplx(0.d0,0.d0)
      ip1    = ip2
      ip2    = ip2 - 1
      do i1  = 1,kdim
      csum4  = csum4 + dconjg(utilde(ip1,i1,ibb1))*utilde(ip2,i1,ibb1)   !4 to 1
      enddo
      phi41 = cdlog(csum4/cdabs(csum4))/xi
*      write(6,'(a30,2i4,3f10.4)')'phi41:',ip1,ip2,phi41




      F(ix,iy)      = cdlog( 
     .              (csum1/cdabs(csum1))*(csum2/cdabs(csum2))*
     .              (csum3/cdabs(csum3))*(csum4/cdabs(csum4)) )/xi 
      Ftilde(ix,iy) = phi12 + phi23 + phi34 + phi41

      rdiff = dabs(F(ix,iy) - Ftilde(ix,iy))
      if(rdiff.gt.0.001d0)then
      write(6,*)
      write(6,'(a20,2f10.4)')'kx ky:',kx/twopi,ky/twopi 
      write(6,'(a20,2f10.4)')'(F-Ftilde):',
     .(F(ix,iy) - Ftilde(ix,iy))/twopi 
      write(6,'(a20,2f10.4)')'F:',F(ix,iy)/twopi 
      write(6,'(a20,2f10.4)')'Ftilde:',Ftilde(ix,iy)/twopi 
      write(6,'(a30,3f10.4)')'phi12:',phi12/twopi
      write(6,'(a30,3f10.4)')'phi23:',phi23/twopi
      write(6,'(a30,3f10.4)')'phi34:',phi34/twopi
      write(6,'(a30,3f10.4)')'phi41:',phi41/twopi
      write(6,*)
      endif
     
      rsumF      = rsumF + cdlog( 
     .              (csum1/cdabs(csum1))*(csum2/cdabs(csum2))*
     .              (csum3/cdabs(csum3))*(csum4/cdabs(csum4)) )/xi 
      rsumFtilde = rsumFtilde + phi12 + phi23 + phi34 + phi41
     

      enddo
      enddo

* Sign due to direction of curve-integral along edges.
      write(6,'(a20,3f10.4)') 'Chern number:', -rsumF/twopi
      write(6,*) 
      write(6,'(a20,3f10.4)') 'Sum F (Chern):', rsumF/twopi
      write(6,'(a20,3f10.4)') 'Sum Ftilde:', rsumFtilde
      write(6,'(a20,3f10.4)') 'Sum F-Ftilde:',phi/pi, 
     .(rsumF-rsumFtilde)/twopi

      rdiff = dabs(csumFtilde)
      if(rdiff.gt.0.001d0)then
      write(6,*)'Error sum Ftilde not zero !!'
      endif

      enddo   !ibb1 loop

      do ibb1 = 1,2
* Edges !!
      arg    = dcmplx(1.d0,0.d0)   !Gauge-invariant
      arg1   = dcmplx(0.d0,0.d0)
      do  ix = 1, nkabc(1) 
      csum1  = dcmplx(0.d0,0.d0)
      ip1    = 1 + (nkabc(1)+1)*(ix-1) 
      ip2    = ip1 + (nkabc(1) + 1)  
*      write(6,'(a30,2i4,3f10.4)')'Index:',ip1,ip2
      do i1  = 1,kdim
      csum1  = csum1 + dconjg(utilde(ip1,i1,ibb1))*utilde(ip2,i1,ibb1)   !1 to 2
      enddo
      arg    = arg*(csum1/cdabs(csum1))
      arg1   = arg1 + cdlog(csum1/cdabs(csum1))/xi
      enddo

      do  ix = 1, nkabc(1)
      csum1  = dcmplx(0.d0,0.d0)
      ip1    = ip2 
      ip2    = ip1 + 1
*      write(6,'(a30,2i4,3f10.4)')'Index:',ip1,ip2
      do i1  = 1,kdim
      csum1  = csum1 + dconjg(utilde(ip1,i1,ibb1))*utilde(ip2,i1,ibb1)   !2 to 3
      enddo
      arg    = arg*(csum1/cdabs(csum1))
      arg1   = arg1 + cdlog(csum1/cdabs(csum1))/xi
      enddo

      do  ix = 1, nkabc(1)
      csum1  = dcmplx(0.d0,0.d0)
      ip1    = ip2 
      ip2    = ip1 - (nkabc(1) + 1)
*      write(6,'(a30,2i4,3f10.4)')'Index:',ip1,ip2
      do i1  = 1,kdim
      csum1  = csum1 + dconjg(utilde(ip1,i1,ibb1))*utilde(ip2,i1,ibb1)   !3 to 4
      enddo
      arg    = arg*(csum1/cdabs(csum1))
      arg1   = arg1 + cdlog(csum1/cdabs(csum1))/xi
      enddo

      do  ix = 1, nkabc(1)
      csum1  = dcmplx(0.d0,0.d0)
      ip1    = ip2
      ip2    = ip1 - 1
*      write(6,'(a30,2i4,3f10.4)')'Index:',ip1,ip2
      do i1  = 1,kdim
      csum1  = csum1 + dconjg(utilde(ip1,i1,ibb1))*utilde(ip2,i1,ibb1)   !4 to 1
      enddo
      arg    = arg*(csum1/cdabs(csum1))
      arg1   = arg1 + cdlog(csum1/cdabs(csum1))/xi
      enddo

      write(6,'(a20,i4,4f12.8)')'Along edges for band:',ibb1,
     .                           cdlog(arg)/xi, arg1
      enddo




      write(6,*)'Done!'

      stop
      end
