      subroutine gibbs(iseed)
      implicit real*8 (a-h, o-z)
	integer, parameter :: ni = 1180, nm = 14, np = 91
	integer, parameter :: nc = 8, nb = 3, nl = 5, nk = 5
			
	real*8 yi(ni)
      real*8 wi(ni), taui(ni), vnu 
	real*8 zi(ni,nc), betaa(nc)
      
      real*8 xim(ni,nm,nb), alpham(nb,nl)  
      real*8 zetam(nl), gammam
      real*8 gm(nm), wgtm(nl), thetam(nl)

      real*8 xip(ni,np,nb), alphap(nb,nk)  
      real*8 zetap(nk), gammap
      real*8 gp(np), wgtp(nk), thetap(nk)
      
      real*8 sigma_betaa
      real*8 sigma_thetam, sigma_thetap
      real*8 a0, b0, a1, b1
      
      real*8 zb(ni)
      real*8 fxim(ni), deltam(nm,nb)
      real*8 fxip(ni), deltap(np,nb)

      real*8 probm(nm,nl), probp(np,nk)
      
      common /vyi/yi
      common /vwi/wi
      common /vtaui/taui
      common /vvnu/vnu
      
      common /vzi/zi
      common /vbetaa/betaa
      
      common /vxim/xim
      common /valpham/alpham
      common /vzetam/zetam
      common /vgammam/gammam
      
      common /vgm/gm
      common /vwgtm/wgtm
      common /vthetam/thetam
            
      common /vxip/xip
      common /valphap/alphap
      common /vzetap/zetap
      common /vgammap/gammap
      
      common /vgp/gp
      common /vwgtp/wgtp
      common /vthetap/thetap
      
      common /vsigma_betaa/sigma_betaa
      common /vsigma_thetam/sigma_thetam
      common /vsigma_thetap/sigma_thetap
      common /va0/a0
      common /vb0/b0
      common /va1/a1
      common /vb1/b1
            
      common /vzb/zb
      common /vfxim/fxim
      common /vdeltam/deltam
      common /vfxip/fxip
      common /vdeltap/deltap
      
      common /vprobm/probm
      common /vprobp/probp
      
      call gibbs_gm(iseed) 
      call gibbs_thetam(iseed)          
      
      call gibbs_gp(iseed) 
      call gibbs_thetap(iseed)          
      
      call gibbs_wi(iseed)          
      call gibbs_taui(iseed)          
      
      call gibbs_alpham(iseed)     
      call gibbs_alphap(iseed)           
      call gibbs_betaa(iseed)          
                  
      call gibbs_gammam(iseed)              
      call gibbs_zetam(iseed)         

      call gibbs_gammap(iseed)              
      call gibbs_zetap(iseed)         
      
      end subroutine
          
      
	subroutine gibbs_gp(iseed)  
      implicit real*8 (a-h, o-z)
	integer, parameter :: ni = 1180, nm = 14, np = 91
	integer, parameter :: nc = 8, nb = 3, nl = 5, nk = 5
		
	real*8 yi(ni), vnu
      real*8 xip(ni,np,nb), alphap(nb,nk)  
      real*8 gp(np), wgtp(nk)
     
      real*8 zb(ni)
      real*8 fxim(ni), fxip(ni), deltap(np,nb)
      
      real*8 probp(np,nk)
      
      real*8 xa, cdf, cprobp(nk)
      real*8 amax, u, summ
      
      common /vyi/yi
      
      common /vvnu/vnu
      
      common /vxip/xip
      common /valphap/alphap
      
      common /vgp/gp
      common /vwgtp/wgtp
            
      common /vzb/zb
      common /vfxim/fxim
      common /vfxip/fxip
      common /vdeltap/deltap
      
      common /vprobp/probp
      
      external drnunf, dtdf
                  
      do j = 1, np
          
          do l = 1, nk
              probp(j,l) = dlog(wgtp(l))
          enddo
          do i = 1, ni
                 
              if (xip(i,j,1) .ne. 0.d0) then
              
                  xa = 0.d0
                  do jj = 1, nb
                      xa = xa + xip(i,j,jj)*deltap(j,jj)
                  enddo
              
                  fxip(i) = fxip(i) - xa                            
              
                  do l = 1, nk
                            
                      xa = 0.d0
                      do jj = 1, nb
                          xa = xa + xip(i,j,jj)*alphap(jj,l)
                      enddo
                      
                      cdf = dtdf(zb(i) + fxim(i) + fxip(i) + xa,vnu)
                  
                      if (yi(i) .eq. 1.d0) then
                                            
                          probp(j,l) = probp(j,l) + dlog(cdf)
                      
                      else

                          probp(j,l) = probp(j,l) + dlog(1.d0 - cdf)
                      
                      endif
              
                  enddo

              endif
                  
          enddo
         	
          amax = probp(j,1)
          do l = 2, nk
              if (amax .lt. probp(j,l)) then
                  amax = probp(j,l)
              endif
          enddo       
                
          summ = 0.d0
          do l = 1, nk
              probp(j,l) = dexp(probp(j,l) - amax)
              summ = summ + probp(j,l)
          enddo
          do l = 1, nk
              probp(j,l) = probp(j,l)/summ
          enddo
        
          cprobp(1) = probp(j,1)
          do l = 2, nk
              cprobp(l) = cprobp(l-1) + probp(j,l)
          enddo

          call rnset(iseed)
          u = drnunf()
          call rnget(iseed)
                		                        
          if (u .le. cprobp(1)) then
              gp(j) = 1.d0
          endif
          do l = 2, nk
              if ((u .gt. cprobp(l-1)) 
     +             .and. (u .le. cprobp(l))) then
                  gp(j) = dfloat(l)
              endif
          enddo
          
          l = nint(gp(j))
          do jj = 1, nb
              deltap(j,jj) = alphap(jj,l)               
          enddo

          do i = 1, ni
              
              xa = 0.d0
              do jj = 1, nb
                  xa = xa + xip(i,j,jj)*deltap(j,jj)
              enddo
              
              fxip(i) = fxip(i) + xa
              
          enddo
          
      enddo
                  
      end subroutine	   
                         

      subroutine gibbs_thetap(iseed) 
      implicit real*8 (a-h, o-z)
	integer, parameter :: ni = 1180, nm = 14, np = 91
	integer, parameter :: nc = 8, nb = 3, nl = 5, nk = 5

      real*8 gp(np), wgtp(nk), thetap(nk)
      real*8 sigma_thetap
      
	real*8 start(1), xmin(1), ynewlo, reqmin, step(1)
	integer konvge, kcount, icount, numres, ifault, nopt

	real*8 bold, anew, amean, asigma
	real*8 rv, u, bpdf, apdf, ratio
      
      integer ldum
      real*8 thetak(nk), sumw
      real*8 suma, sumb, sumc, sumd, summ
            
      common /vgp/gp
      common /vwgtp/wgtp
      common /vthetap/thetap
                  
      common /vsigma_thetap/sigma_thetap
      
	common /vldum/ldum
	common /vthetak/thetak
      
      external fthetap, drnnof, drnunf
            
      thetak(1) = 0.d0      
      do l = 2, nk
          if (l .eq. 2) then
              thetak(l) = dlog(-thetap(l))
          else
              thetak(l) = dlog(thetap(l-1) - thetap(l))
          endif
      enddo
            
      do l = 2, nk
          
          ldum = l
          
          bold = thetak(l)
                                  
          nopt = 1
	    reqmin = 1.0d-10 ; konvge = 5 ; kcouns = 1000
	    step(1) = 0.2d0 ; start(1) = bold
          call nelmin(fthetap, nopt, start, xmin, ynewlo,  
     +                reqmin, step, konvge, kcount, 
     +                icount, numres, ifault)
	    amean = xmin(1)
          thetak(l) = amean
                             
          thetap(1) = 0.d0      
          do ll = 2, nk
              if (ll .eq. 2) then
                  thetap(ll) = -dexp(thetak(ll))
              else
                  thetap(ll) = thetap(ll-1) - dexp(thetak(ll))
              endif
          enddo
          
          summ = 0.d0
          do j = 1, np          

              suma = 1.d0
              do ll = 2, nk
                  suma = suma + dexp(thetap(ll))                  
              enddo
              
              sumb = 0.d0; sumc = 0.d0 
              do ll = l, nk
                  
                  if (nint(gp(j)) .eq. ll) then
                      sumb = sumb + 1.d0
                  endif
                      
                  sumc = sumc + dexp(thetap(ll))
                                    
              enddo

              summ = summ 
     +             - sumb*dexp(thetak(l))
     +             + ( dexp(thetak(l))
     +                 *(1.d0 - dexp(thetak(l)))
     +                 *sumc*suma 
     +                 + (dexp(thetak(l))*sumc)**2 )
     +               /suma**2
                            
          enddo
          
          sumd = 0.d0 
          do ll = l, nk
                                    
              sumd = sumd
     +             + dexp(thetak(l))*thetap(ll)/sigma_thetap
     +             - dexp(2.d0*thetak(l))/sigma_thetap
                  
          enddo  
          
          asigma = summ + sumd
          
          asigma = -1.d0/asigma*2.d0
          if (asigma .lt. 0.0d0) asigma = -asigma      
                                                                        
          bpdf = -fthetap(bold) 
     +         + (bold - amean)**2/(2.d0*asigma)
        
          do ii = 1, 100
	  
              call rnset(iseed)
              rv = drnnof()
              call rnget(iseed)
              anew = amean + rv*dsqrt(asigma)

              apdf = -fthetap(anew) 
     +             + (anew - amean)**2/(2.d0*asigma)
              ratio = apdf - bpdf 
     
              if (ratio .ge. 0.0d0) then
                  bpdf = apdf
                  bold = anew
              else
                  call rnset(iseed)
                  u = drnunf()
                  call rnget(iseed)
                  if (dlog(u) .le. ratio) then
                      bpdf = apdf
                      bold = anew
                  endif
              endif
		
          enddo     
                    
          thetak(l) = bold
                                   
      enddo
      
      thetap(1) = 0.d0      
      do l = 2, nk
          if (l .eq. 2) then
              thetap(l) = -dexp(thetak(l))
          else
              thetap(l) = thetap(l-1) - dexp(thetak(l))
          endif
      enddo
                    
      sumw = 0.d0
      do l = 1, nk
          sumw = sumw + dexp(thetap(l))
      enddo
      do l = 1, nk
          wgtp(l) = dexp(thetap(l))/sumw
      enddo
      
      end subroutine                   
             
      real*8 function fthetap(star)
      implicit real*8 (a-h, o-z)
	integer, parameter :: ni = 1180, nm = 14, np = 91
	integer, parameter :: nc = 8, nb = 3, nl = 5, nk = 5

      real*8 gp(np), wgtp(nk), thetap(nk)
      real*8 sigma_thetap

      integer ldum
      real*8 thetak(nk)

      real*8 star, sumw, sum1, sum2, pdf
      
      common /vgp/gp
                  
      common /vsigma_thetap/sigma_thetap
      
	common /vldum/ldum
	common /vthetak/thetak
            
      l = ldum
                   
      thetak(l) = star
                   
      thetap(1) = 0.d0      
      do ll = 2, nk
          if (ll .eq. 2) then
              thetap(ll) = -dexp(thetak(ll))
          else
              thetap(ll) = thetap(ll-1) - dexp(thetak(ll))
          endif
      enddo
                
      sumw = 0.d0
      do ll = 1, nk
          sumw = sumw + dexp(thetap(ll))
      enddo
      do ll = 1, nk
          wgtp(ll) = dexp(thetap(ll))/sumw
      enddo
                             
      sum1 = 0.d0
      do j = 1, np         
          ll = nint(gp(j))
          sum1 = sum1 + dlog(wgtp(ll))
      enddo
            
      sum2 = 0.d0
      do ll = 2, nk
          sum2 = sum2 
     +         + thetak(ll) 
     +         - thetap(ll)**2/(2.d0*sigma_thetap)
      enddo
      
      pdf = sum1 + sum2 
                          
      fthetap = -pdf
                        
      end function 
      
      
      subroutine gibbs_gammap(iseed)
      implicit real*8 (a-h, o-z)
	integer, parameter :: ni = 1180, nm = 14, np = 91
	integer, parameter :: nc = 8, nb = 3, nl = 5, nk = 5

	real*8 alphap(nb,nk), gammap
      real*8 a1, b1
      
	real*8 start(1), xmin(1), ynewlo, reqmin, step(1)
	integer konvge, kcount, icount, numres, ifault, nopt

	real*8 bold, anew, amean, asigma
	real*8 rv, u, bpdf, apdf, ratio
      
      real*8 summ, temp
      
      common /valphap/alphap
      common /vgammap/gammap
                        
      common /va1/a1
      common /vb1/b1
            
      external fgammap, drnnof, drnunf
       
      bold = dlog(gammap)

      nopt = 1
      reqmin = 1.0d-10 ; konvge = 5 ; kcouns = 1000
      step(1) = 0.2d0 ; start(1) = bold
      call nelmin(fgammap, nopt, start, xmin, ynewlo,  
     +            reqmin, step, konvge, kcount, 
     +            icount, numres, ifault)
      amean = xmin(1)
      gammap = dexp(amean)
      
      summ = 0.d0
      do l = 1, nk                             
          
          temp = 0.d0
          do jj = 1, nb
              temp = temp + alphap(jj,l)**2
          enddo          
          
          summ = summ 
     +         - dsqrt(temp)*dexp(amean/2.d0)/4.d0
          
      enddo
            
      asigma = summ - b1*dexp(amean)
                
      asigma = -1.d0/asigma*2.d0       
      if (asigma .lt. 0.0d0) asigma = -asigma      
                 
      bpdf = -fgammap(bold) 
     +     + (bold - amean)**2/(2.d0*asigma)
        
      do ii = 1, 100
	  
          call rnset(iseed)
          rv = drnnof()
          call rnget(iseed)
          anew = amean + rv*dsqrt(asigma)

          apdf = -fgammap(anew) 
     +         + (anew - amean)**2/(2.d0*asigma)
          ratio = apdf - bpdf 
     
          if (ratio .ge. 0.0d0) then
              bpdf = apdf
              bold = anew
          else
              call rnset(iseed)
              u = drnunf()
              call rnget(iseed)
              if (dlog(u) .le. ratio) then
                  bpdf = apdf
                  bold = anew
              endif
          endif
		
      enddo     

      gammap = dexp(bold)      
          
      end subroutine                   
   
      real*8 function fgammap(star)
      implicit real*8 (a-h, o-z)
	integer, parameter :: ni = 1180, nm = 14, np = 91
	integer, parameter :: nc = 8, nb = 3, nl = 5, nk = 5

	real*8 alphap(nb,nk)
      real*8 a1, b1
      
      real*8 star, summ, temp, pdf
      
      common /valphap/alphap
                        
      common /va1/a1
      common /vb1/b1
                  
      summ = 0.d0
      do l = 1, nk                            
          
          temp = 0.d0
          do jj = 1, nb
              temp = temp + alphap(jj,l)**2
          enddo          
          
          summ = summ 
     +         + dfloat(nb)*star/2.d0
     +         - dsqrt(temp)*dexp(star/2.d0)
          
      enddo
      
      pdf = summ + a1*star - b1*dexp(star)

      fgammap = -pdf
    
      end function
      
      
      subroutine gibbs_zetap(iseed)
      implicit real*8 (a-h, o-z)
	integer, parameter :: ni = 1180, nm = 14, np = 91
	integer, parameter :: nc = 8, nb = 3, nl = 5, nk = 5

	real*8 alphap(nb,nk)
      real*8 zetap(nk), gammap
      
      real*8 pp, aa, bb, temp, rv
      
      common /valphap/alphap
      common /vzetap/zetap
      common /vgammap/gammap
                                    
      external rgenGIG
         
      do l = 1, nk
                    
          temp = 0.d0
          do jj = 1, nb
              temp = temp + alphap(jj,l)**2
          enddo
          
          pp = -1.d0/2.d0
          aa = gammap*temp
          bb = 1.d0
          
          rv = rgenGIG(pp,aa,bb,iseed)
          zetap(l) = 1.d0/rv
          
      enddo
                    
      end subroutine                   
        
      
      subroutine gibbs_alphap(iseed)
      implicit real*8 (a-h, o-z)
	integer, parameter :: ni = 1180, nm = 14, np = 91
	integer, parameter :: nc = 8, nb = 3, nl = 5, nk = 5

      real*8 wi(ni), taui(ni)
	real*8 zi(ni,nc)
      real*8 xip(ni,np,nb), alphap(nb,nk)  
      real*8 zetap(nk), gammap
      real*8 gp(np)
                
      real*8 sigma_betaa
      
      real*8 fxim(ni), fxip(ni), deltap(np,nb)
      
      real*8 xstar(nb), wstar, xa, ffxip
      real*8 sumxx(nb,nb), sumxw(nb)
      real*8 sumzz(nc,nc), sumzx(nc,nb), sumzw(nc)
      real*8 AA(nc,nc), AAi(nc,nc)
      real*8 BB(nb,nc), CC(nb,nb), DD(nb)
      
      real*8 hmean(nb), hsigma(nb,nb)
      real*8 amean(nb), asigma(nb,nb)
	real*8 tol, rsig(nb,nb), rv(nb)
      
      real*8 xxip(np*nb), ddelta(np*nb)
            
      common /vwi/wi
      common /vtaui/taui

      common /vzi/zi
      
      common /vxip/xip
      common /valphap/alphap
      common /vzetap/zetap
      common /vgammap/gammap
      
      common /vgp/gp
            
      common /vsigma_betaa/sigma_betaa
      
      common /vfxim/fxim
      common /vfxip/fxip
      common /vdeltap/deltap
      
      external dlinrg, dblinf, dmach, dchfac, drnmvn
            
      do l = 1, nk
                                                      
          do j1 = 1, nb
              
              sumxw(j1) = 0.d0
              
              do j2 = 1, nb
                  sumxx(j1,j2) = 0.d0
              enddo
                      
          enddo
              
          do j1 = 1, nc      
        
              sumzw(j1) = 0.d0
              
              do j2 = 1, nc           
                  sumzz(j1,j2) = 0.d0
              enddo

              do j2 = 1, nb           
                  sumzx(j1,j2) = 0.d0
              enddo
                  
          enddo             
          
          do i = 1, ni
                                  
              do jj = 1, nb
                  xstar(jj) = 0.d0
              enddo              
              ffxip = 0.d0               
              do j = 1, np
                  if (xip(i,j,1) .ne. 0.d0) then
                  
                      ll = nint(gp(j))
                
                      if (ll .eq. l) then
                      
                          do jj = 1, nb
                              xstar(jj) = xstar(jj) + xip(i,j,jj)
                          enddo
                      
                      else                          
                          
                          xa = 0.d0
                          do jj = 1, nb
                              xa = xa + xip(i,j,jj)*alphap(jj,ll)
                          enddo
                      
                          ffxip = ffxip + xa
                          
                      endif  
                      
                  endif  
              enddo
                            
              wstar = wi(i) - fxim(i) - ffxip 
                                                    
              do j1 = 1, nb
                      
                  sumxw(j1) = sumxw(j1) + taui(i)*xstar(j1)*wstar
                      
                  do j2 = 1, nb
                  
                      sumxx(j1,j2) = sumxx(j1,j2) 
     +                             + taui(i)*xstar(j1)*xstar(j2)
                      
                  enddo
                      
              enddo
              
              do j1 = 1, nc      
        
                  sumzw(j1) = sumzw(j1) + taui(i)*zi(i,j1)*wstar
            
                  do j2 = 1, nc           
                
                      sumzz(j1,j2) = sumzz(j1,j2) 
     +                             + taui(i)*zi(i,j1)*zi(i,j2)

                  enddo

                  do j2 = 1, nb           
                
                      sumzx(j1,j2) = sumzx(j1,j2) 
     +                             + taui(i)*zi(i,j1)*xstar(j2)

                  enddo
                  
              enddo             
                                                  
          enddo
            
          do j1 = 1, nc                           
              do j2 = 1, nc           
                  AA(j1,j2) = sumzz(j1,j2)
              enddo                  
              AA(j1,j1) = AA(j1,j1) + 1.d0/sigma_betaa
          enddo      
                    
          call dlinrg(nc, AA, nc, AAi, nc)

          do j1 = 1, nb                           
              do j2 = 1, nc
                  BB(j1,j2) = 0.d0
                  do jj = 1, nc
                      BB(j1,j2) = BB(j1,j2) 
     +                          + sumzx(jj,j1)*AAi(jj,j2)
                  enddo
              enddo                  
          enddo             

          do j1 = 1, nb 
              
              do j2 = 1, nb
                  CC(j1,j2) = 0.d0
                  do jj = 1, nc
                      CC(j1,j2) = CC(j1,j2) 
     +                          + BB(j1,jj)*sumzx(jj,j2)
                  enddo
              enddo        
              
              DD(j1) = 0.d0
              do j2 = 1, nc
                  DD(j1) = DD(j1) + BB(j1,j2)*sumzw(j2)
              enddo
              
          enddo             
          
          do j1 = 1, nb     
        
              hmean(j1) = sumxw(j1) - DD(j1)
            
              do j2 = 1, nb
                
                  hsigma(j1,j2) = sumxx(j1,j2) 
     +                          -  CC(j1,j2)

              enddo
              
              hsigma(j1,j1) = hsigma(j1,j1) 
     +                      + gammap/zetap(l)
              
          enddo             
          
          call dlinrg(nb, hsigma, nb, asigma, nb)
	
          do j1 = 1, nb
              amean(j1) = 0.d0
              do j2 = 1, nb
                  amean(j1) = amean(j1) 
     +                      + asigma(j1,j2)*hmean(j2)
              enddo
          enddo

          tol = 100.d0*dmach(4)
          call dchfac(nb, asigma, nb, tol, irank, rsig, nb)
          call rnset(iseed)
          call drnmvn(1, nb, rsig, nb, rv, 1)
          call rnget(iseed)      
	
          do jj = 1, nb
              alphap(jj,l)= amean(jj) + rv(jj)
          enddo  
          
      enddo
       
      do j = 1, np
          l = nint(gp(j))
          do jj = 1, nb
              deltap(j,jj) = alphap(jj,l)               
          enddo
      enddo
      do i = 1, ni
          xxip = reshape(xip(i,:,:),(/np*nb/))
          ddelta = reshape(deltap,(/np*nb/))
          fxip(i) = dot_product(xxip,ddelta)
      enddo
      
      end subroutine  
            
      
	subroutine gibbs_gm(iseed)  
      implicit real*8 (a-h, o-z)
	integer, parameter :: ni = 1180, nm = 14, np = 91
	integer, parameter :: nc = 8, nb = 3, nl = 5, nk = 5
		
	real*8 yi(ni), vnu
      real*8 xim(ni,nm,nb), alpham(nb,nl)  
      real*8 gm(nm), wgtm(nl)
     
      real*8 zb(ni)
      real*8 fxim(ni), deltam(nm,nb), fxip(ni)
      
      real*8 probm(nm,nl)
      
      real*8 xa, cdf, cprobm(nl)
      real*8 amax, u, summ
      
      common /vyi/yi
      
      common /vvnu/vnu
      
      common /vxim/xim
      common /valpham/alpham
      
      common /vgm/gm
      common /vwgtm/wgtm
            
      common /vzb/zb
      common /vfxim/fxim
      common /vdeltam/deltam
      common /vfxip/fxip
      
      common /vprobm/probm
      
      external drnunf, dtdf
                  
      do j = 1, nm
          
          do l = 1, nl
              probm(j,l) = dlog(wgtm(l))
          enddo
          do i = 1, ni
                 
              if (xim(i,j,1) .ne. 0.d0) then
              
                  xa = 0.d0
                  do jj = 1, nb
                      xa = xa + xim(i,j,jj)*deltam(j,jj)
                  enddo
              
                  fxim(i) = fxim(i) - xa                            
              
                  do l = 1, nl
                            
                      xa = 0.d0
                      do jj = 1, nb
                          xa = xa + xim(i,j,jj)*alpham(jj,l)
                      enddo
                      
                      cdf = dtdf(zb(i) + fxim(i) + xa + fxip(i),vnu)
                  
                      if (yi(i) .eq. 1.d0) then
                                            
                          probm(j,l) = probm(j,l) + dlog(cdf)
                      
                      else

                          probm(j,l) = probm(j,l) + dlog(1.d0 - cdf)
                      
                      endif
              
                  enddo

              endif
                  
          enddo
         	
          amax = probm(j,1)
          do l = 2, nl
              if (amax .lt. probm(j,l)) then
                  amax = probm(j,l)
              endif
          enddo       
                
          summ = 0.d0
          do l = 1, nl
              probm(j,l) = dexp(probm(j,l) - amax)
              summ = summ + probm(j,l)
          enddo
          do l = 1, nl
              probm(j,l) = probm(j,l)/summ
          enddo
        
          cprobm(1) = probm(j,1)
          do l = 2, nl
              cprobm(l) = cprobm(l-1) + probm(j,l)
          enddo

          call rnset(iseed)
          u = drnunf()
          call rnget(iseed)
                		                        
          if (u .le. cprobm(1)) then
              gm(j) = 1.d0
          endif
          do l = 2, nl
              if ((u .gt. cprobm(l-1)) 
     +             .and. (u .le. cprobm(l))) then
                  gm(j) = dfloat(l)
              endif
          enddo
          
          l = nint(gm(j))
          do jj = 1, nb
              deltam(j,jj) = alpham(jj,l)               
          enddo

          do i = 1, ni
              
              xa = 0.d0
              do jj = 1, nb
                  xa = xa + xim(i,j,jj)*deltam(j,jj)
              enddo
              
              fxim(i) = fxim(i) + xa
              
          enddo
          
      enddo
                  
      end subroutine	   
                   

      subroutine gibbs_thetam(iseed) 
      implicit real*8 (a-h, o-z)
	integer, parameter :: ni = 1180, nm = 14, np = 91
	integer, parameter :: nc = 8, nb = 3, nl = 5, nk = 5

      real*8 gm(nm), wgtm(nl), thetam(nl)
      real*8 sigma_thetam
      
	real*8 start(1), xmin(1), ynewlo, reqmin, step(1)
	integer konvge, kcount, icount, numres, ifault, nopt

	real*8 bold, anew, amean, asigma
	real*8 rv, u, bpdf, apdf, ratio
      
      integer ldum
      real*8 thetal(nl), sumw
      real*8 suma, sumb, sumc, sumd, summ
            
      common /vgm/gm
      common /vwgtm/wgtm
      common /vthetam/thetam
                  
      common /vsigma_thetam/sigma_thetam
      
	common /vldum/ldum
	common /vthetal/thetal
      
      external fthetam, drnnof, drnunf
            
      thetal(1) = 0.d0      
      do l = 2, nl
          if (l .eq. 2) then
              thetal(l) = dlog(-thetam(l))
          else
              thetal(l) = dlog(thetam(l-1) - thetam(l))
          endif
      enddo
            
      do l = 2, nl
          
          ldum = l
          
          bold = thetal(l)
                                  
          nopt = 1
	    reqmin = 1.0d-10 ; konvge = 5 ; kcouns = 1000
	    step(1) = 0.2d0 ; start(1) = bold
          call nelmin(fthetam, nopt, start, xmin, ynewlo,  
     +                reqmin, step, konvge, kcount, 
     +                icount, numres, ifault)
	    amean = xmin(1)
          thetal(l) = amean
                             
          thetam(1) = 0.d0      
          do ll = 2, nl
              if (ll .eq. 2) then
                  thetam(ll) = -dexp(thetal(ll))
              else
                  thetam(ll) = thetam(ll-1) - dexp(thetal(ll))
              endif
          enddo
          
          summ = 0.d0
          do j = 1, nm        

              suma = 1.d0
              do ll = 2, nl
                  suma = suma + dexp(thetam(ll))                  
              enddo
              
              sumb = 0.d0; sumc = 0.d0 
              do ll = l, nl
                  
                  if (nint(gm(j)) .eq. ll) then
                      sumb = sumb + 1.d0
                  endif
                      
                  sumc = sumc + dexp(thetam(ll))
                                    
              enddo

              summ = summ 
     +             - sumb*dexp(thetal(l))
     +             + ( dexp(thetal(l))
     +                 *(1.d0 - dexp(thetal(l)))
     +                 *sumc*suma 
     +                 + (dexp(thetal(l))*sumc)**2 )
     +               /suma**2
                            
          enddo
          
          sumd = 0.d0 
          do ll = l, nl
                                    
              sumd = sumd
     +             + dexp(thetal(l))*thetam(ll)/sigma_thetam
     +             - dexp(2.d0*thetal(l))/sigma_thetam
                  
          enddo  
          
          asigma = summ + sumd
          
          asigma = -1.d0/asigma*2.d0
          if (asigma .lt. 0.0d0) asigma = -asigma      
                                                                  
          bpdf = -fthetam(bold) 
     +         + (bold - amean)**2/(2.d0*asigma)
        
          do ii = 1, 100
	  
              call rnset(iseed)
              rv = drnnof()
              call rnget(iseed)
              anew = amean + rv*dsqrt(asigma)

              apdf = -fthetam(anew) 
     +             + (anew - amean)**2/(2.d0*asigma)
              ratio = apdf - bpdf 
     
              if (ratio .ge. 0.0d0) then
                  bpdf = apdf
                  bold = anew
              else
                  call rnset(iseed)
                  u = drnunf()
                  call rnget(iseed)
                  if (dlog(u) .le. ratio) then
                      bpdf = apdf
                      bold = anew
                  endif
              endif
		
          enddo     
                    
          thetal(l) = bold
                                   
      enddo
      
      thetam(1) = 0.d0      
      do l = 2, nl
          if (l .eq. 2) then
              thetam(l) = -dexp(thetal(l))
          else
              thetam(l) = thetam(l-1) - dexp(thetal(l))
          endif
      enddo
                    
      sumw = 0.d0
      do l = 1, nl
          sumw = sumw + dexp(thetam(l))
      enddo
      do l = 1, nl
          wgtm(l) = dexp(thetam(l))/sumw
      enddo
      
      end subroutine                   
             
      real*8 function fthetam(star)
      implicit real*8 (a-h, o-z)
	integer, parameter :: ni = 1180, nm = 14, np = 91
	integer, parameter :: nc = 8, nb = 3, nl = 5, nk = 5

      real*8 gm(nm), wgtm(nl), thetam(nl)
      real*8 sigma_thetam

      integer ldum
      real*8 thetal(nl)

      real*8 star, sumw, sum1, sum2, pdf
      
      common /vgm/gm
                  
      common /vsigma_thetam/sigma_thetam
      
	common /vldum/ldum
	common /vthetal/thetal
            
      l = ldum
                   
      thetal(l) = star
                   
      thetam(1) = 0.d0      
      do ll = 2, nl
          if (ll .eq. 2) then
              thetam(ll) = -dexp(thetal(ll))
          else
              thetam(ll) = thetam(ll-1) - dexp(thetal(ll))
          endif
      enddo
                
      sumw = 0.d0
      do ll = 1, nl
          sumw = sumw + dexp(thetam(ll))
      enddo
      do ll = 1, nl
          wgtm(ll) = dexp(thetam(ll))/sumw
      enddo
                             
      sum1 = 0.d0
      do j = 1, nm          
          ll = nint(gm(j))
          sum1 = sum1 + dlog(wgtm(ll))
      enddo
            
      sum2 = 0.d0
      do ll = 2, nl
          sum2 = sum2 
     +         + thetal(ll) 
     +         - thetam(ll)**2/(2.d0*sigma_thetam)
      enddo
      
      pdf = sum1 + sum2 
                          
      fthetam = -pdf
                        
      end function 
         
      
      subroutine gibbs_gammam(iseed)
      implicit real*8 (a-h, o-z)
	integer, parameter :: ni = 1180, nm = 14, np = 91
	integer, parameter :: nc = 8, nb = 3, nl = 5, nk = 5

	real*8 alpham(nb,nl), gammam
      real*8 a0, b0
      
	real*8 start(1), xmin(1), ynewlo, reqmin, step(1)
	integer konvge, kcount, icount, numres, ifault, nopt

	real*8 bold, anew, amean, asigma
	real*8 rv, u, bpdf, apdf, ratio
      
      real*8 summ, temp
      
      common /valpham/alpham
      common /vgammam/gammam
                        
      common /va0/a0
      common /vb0/b0
            
      external fgammam, drnnof, drnunf
       
      bold = dlog(gammam)

      nopt = 1
      reqmin = 1.0d-10 ; konvge = 5 ; kcouns = 1000
      step(1) = 0.2d0 ; start(1) = bold
      call nelmin(fgammam, nopt, start, xmin, ynewlo,  
     +            reqmin, step, konvge, kcount, 
     +            icount, numres, ifault)
      amean = xmin(1)
      gammam = dexp(amean)
      
      summ = 0.d0
      do l = 1, nl                              
          
          temp = 0.d0
          do jj = 1, nb
              temp = temp + alpham(jj,l)**2
          enddo          
          
          summ = summ 
     +         - dsqrt(temp)*dexp(amean/2.d0)/4.d0
          
      enddo
            
      asigma = summ - b0*dexp(amean)
                
      asigma = -1.d0/asigma*2.d0       
      if (asigma .lt. 0.0d0) asigma = -asigma      
                 
      bpdf = -fgammam(bold) 
     +     + (bold - amean)**2/(2.d0*asigma)
        
      do ii = 1, 100
	  
          call rnset(iseed)
          rv = drnnof()
          call rnget(iseed)
          anew = amean + rv*dsqrt(asigma)

          apdf = -fgammam(anew) 
     +         + (anew - amean)**2/(2.d0*asigma)
          ratio = apdf - bpdf 
     
          if (ratio .ge. 0.0d0) then
              bpdf = apdf
              bold = anew
          else
              call rnset(iseed)
              u = drnunf()
              call rnget(iseed)
              if (dlog(u) .le. ratio) then
                  bpdf = apdf
                  bold = anew
              endif
          endif
		
      enddo     

      gammam = dexp(bold)      
          
      end subroutine                   
   
      real*8 function fgammam(star)
      implicit real*8 (a-h, o-z)
	integer, parameter :: ni = 1180, nm = 14, np = 91
	integer, parameter :: nc = 8, nb = 3, nl = 5, nk = 5

	real*8 alpham(nb,nl)
      real*8 a0, b0
      
      real*8 star, summ, temp, pdf
      
      common /valpham/alpham
                        
      common /va0/a0
      common /vb0/b0
                  
      summ = 0.d0
      do l = 1, nl                              
          
          temp = 0.d0
          do jj = 1, nb
              temp = temp + alpham(jj,l)**2
          enddo          
          
          summ = summ 
     +         + dfloat(nb)*star/2.d0
     +         - dsqrt(temp)*dexp(star/2.d0)
          
      enddo
      
      pdf = summ + a0*star - b0*dexp(star)

      fgammam = -pdf
    
      end function

      
      subroutine gibbs_zetam(iseed)
      implicit real*8 (a-h, o-z)
	integer, parameter :: ni = 1180, nm = 14, np = 91
	integer, parameter :: nc = 8, nb = 3, nl = 5, nk = 5

	real*8 alpham(nb,nl)
      real*8 zetam(nl), gammam
      
      real*8 pp, aa, bb, temp, rv
      
      common /valpham/alpham
      common /vzetam/zetam
      common /vgammam/gammam
                                    
      external rgenGIG
         
      do l = 1, nl
                    
          temp = 0.d0
          do jj = 1, nb
              temp = temp + alpham(jj,l)**2
          enddo
          
          pp = -1.d0/2.d0
          aa = gammam*temp
          bb = 1.d0
          
          rv = rgenGIG(pp,aa,bb,iseed)
          zetam(l) = 1.d0/rv
          
      enddo
                    
      end subroutine                   
        
   
      subroutine gibbs_alpham(iseed)
      implicit real*8 (a-h, o-z)
	integer, parameter :: ni = 1180, nm = 14, np = 91
	integer, parameter :: nc = 8, nb = 3, nl = 5, nk = 5

      real*8 wi(ni), taui(ni)
	real*8 zi(ni,nc)
      real*8 xim(ni,nm,nb), alpham(nb,nl)  
      real*8 zetam(nl), gammam
      real*8 gm(nm)
                
      real*8 sigma_betaa
      
      real*8 fxim(ni), deltam(nm,nb), fxip(ni)
      
      real*8 xstar(nb), wstar, xa, ffxim
      real*8 sumxx(nb,nb), sumxw(nb)
      real*8 sumzz(nc,nc), sumzx(nc,nb), sumzw(nc)
      real*8 AA(nc,nc), AAi(nc,nc)
      real*8 BB(nb,nc), CC(nb,nb), DD(nb)
      
      real*8 hmean(nb), hsigma(nb,nb)
      real*8 amean(nb), asigma(nb,nb)
	real*8 tol, rsig(nb,nb), rv(nb)
      
      real*8 xxim(nm*nb), ddelta(nm*nb)
            
      common /vwi/wi
      common /vtaui/taui

      common /vzi/zi
      
      common /vxim/xim
      common /valpham/alpham
      common /vzetam/zetam
      common /vgammam/gammam
      
      common /vgm/gm
            
      common /vsigma_betaa/sigma_betaa
      
      common /vfxim/fxim
      common /vdeltam/deltam
      common /vfxip/fxip
      
      external dlinrg, dblinf, dmach, dchfac, drnmvn
            
      do l = 1, nl
                                                      
          do j1 = 1, nb
              
              sumxw(j1) = 0.d0
              
              do j2 = 1, nb
                  sumxx(j1,j2) = 0.d0
              enddo
                      
          enddo
              
          do j1 = 1, nc      
        
              sumzw(j1) = 0.d0
              
              do j2 = 1, nc           
                  sumzz(j1,j2) = 0.d0
              enddo

              do j2 = 1, nb           
                  sumzx(j1,j2) = 0.d0
              enddo
                  
          enddo             
          
          do i = 1, ni
                                  
              do jj = 1, nb
                  xstar(jj) = 0.d0
              enddo              
              ffxim = 0.d0               
              do j = 1, nm              
                  if (xim(i,j,1) .ne. 0.d0) then
                  
                      ll = nint(gm(j))
                
                      if (ll .eq. l) then
                      
                          do jj = 1, nb
                              xstar(jj) = xstar(jj) + xim(i,j,jj)
                          enddo
                      
                      else                          
                          
                          xa = 0.d0
                          do jj = 1, nb
                              xa = xa + xim(i,j,jj)*alpham(jj,ll)
                          enddo
                      
                          ffxim = ffxim + xa
                          
                      endif  
                      
                  endif  
              enddo
                            
              wstar = wi(i) - ffxim - fxip(i) 
                                                    
              do j1 = 1, nb
                      
                  sumxw(j1) = sumxw(j1) + taui(i)*xstar(j1)*wstar
                      
                  do j2 = 1, nb
                  
                      sumxx(j1,j2) = sumxx(j1,j2) 
     +                             + taui(i)*xstar(j1)*xstar(j2)
                      
                  enddo
                      
              enddo
              
              do j1 = 1, nc      
        
                  sumzw(j1) = sumzw(j1) + taui(i)*zi(i,j1)*wstar
            
                  do j2 = 1, nc           
                
                      sumzz(j1,j2) = sumzz(j1,j2) 
     +                             + taui(i)*zi(i,j1)*zi(i,j2)

                  enddo

                  do j2 = 1, nb           
                
                      sumzx(j1,j2) = sumzx(j1,j2) 
     +                             + taui(i)*zi(i,j1)*xstar(j2)

                  enddo
                  
              enddo             
                                                  
          enddo
            
          do j1 = 1, nc                           
              do j2 = 1, nc           
                  AA(j1,j2) = sumzz(j1,j2)
              enddo                  
              AA(j1,j1) = AA(j1,j1) + 1.d0/sigma_betaa
          enddo      
                    
          call dlinrg(nc, AA, nc, AAi, nc)

          do j1 = 1, nb                           
              do j2 = 1, nc
                  BB(j1,j2) = 0.d0
                  do jj = 1, nc
                      BB(j1,j2) = BB(j1,j2) 
     +                          + sumzx(jj,j1)*AAi(jj,j2)
                  enddo
              enddo                  
          enddo             

          do j1 = 1, nb 
              
              do j2 = 1, nb
                  CC(j1,j2) = 0.d0
                  do jj = 1, nc
                      CC(j1,j2) = CC(j1,j2) 
     +                          + BB(j1,jj)*sumzx(jj,j2)
                  enddo
              enddo        
              
              DD(j1) = 0.d0
              do j2 = 1, nc
                  DD(j1) = DD(j1) + BB(j1,j2)*sumzw(j2)
              enddo
              
          enddo             
          
          do j1 = 1, nb     
        
              hmean(j1) = sumxw(j1) - DD(j1)
            
              do j2 = 1, nb
                
                  hsigma(j1,j2) = sumxx(j1,j2) 
     +                          -  CC(j1,j2)

              enddo
              
              hsigma(j1,j1) = hsigma(j1,j1) 
     +                      + gammam/zetam(l)
              
          enddo             
          
          call dlinrg(nb, hsigma, nb, asigma, nb)
	
          do j1 = 1, nb
              amean(j1) = 0.d0
              do j2 = 1, nb
                  amean(j1) = amean(j1) 
     +                      + asigma(j1,j2)*hmean(j2)
              enddo
          enddo

          tol = 100.d0*dmach(4)
          call dchfac(nb, asigma, nb, tol, irank, rsig, nb)
          call rnset(iseed)
          call drnmvn(1, nb, rsig, nb, rv, 1)
          call rnget(iseed)      
	
          do jj = 1, nb
              alpham(jj,l)= amean(jj) + rv(jj)
          enddo  
          
      enddo
       
      do j = 1, nm
          l = nint(gm(j))
          do jj = 1, nb
              deltam(j,jj) = alpham(jj,l)               
          enddo
      enddo
      do i = 1, ni
          xxim = reshape(xim(i,:,:),(/nm*nb/))
          ddelta = reshape(deltam,(/nm*nb/))
          fxim(i) = dot_product(xxim,ddelta)
      enddo
      
      end subroutine                   
                 

      subroutine gibbs_betaa(iseed)
      implicit real*8 (a-h, o-z)
	integer, parameter :: ni = 1180, nm = 14, np = 91
	integer, parameter :: nc = 8, nb = 3, nl = 5, nk = 5
		
      real*8 wi(ni), taui(ni)
	real*8 zi(ni,nc), betaa(nc)
            
      real*8 sigma_betaa
      
      real*8 zb(ni), fxim(ni), fxip(ni)
      
      real*8 hmean(nc), hsigma(nc,nc)
      real*8 amean(nc), asigma(nc,nc)
	real*8 tol, rsig(nc,nc), rv(nc)
      
      real*8 wstar
            
      common /vwi/wi
      common /vtaui/taui

      common /vzi/zi      
      common /vbetaa/betaa
                  
      common /vsigma_betaa/sigma_betaa
      
      common /vzb/zb
      common /vfxim/fxim
      common /vfxip/fxip
      
      external dmach, dchfac, drnmvn
      
      do j1 = 1, nc
          
          hmean(j1) = 0.d0
          
          do j2 = 1, nc
              hsigma(j1,j2) = 0.d0
          enddo
          hsigma(j1,j1) = 1.d0/sigma_betaa
                        
      enddo
      
      do i = 1, ni
                          
          wstar = wi(i) - fxim(i) - fxip(i)
          
          do j1 = 1, nc      
        
              hmean(j1) = hmean(j1) + taui(i)*zi(i,j1)*wstar
            
              do j2 = 1, nc           
                
                  hsigma(j1,j2) = hsigma(j1,j2) 
     +                          + taui(i)*zi(i,j1)*zi(i,j2)

              enddo

          enddo             
                              		                        
      enddo  
                           
	call dlinrg(nc, hsigma, nc, asigma, nc)
	
	do j1 = 1, nc
	    amean(j1) = 0.d0
		do j2 = 1, nc
		    amean(j1) = amean(j1) + asigma(j1,j2)*hmean(j2)
		enddo
      enddo

	tol = 100.d0*dmach(4)
	call dchfac(nc, asigma, nc, tol, irank, rsig, nc)
	call rnset(iseed)
	call drnmvn(1, nc, rsig, nc, rv, 1)
	call rnget(iseed)   
		
	do jj = 1, nc
          betaa(jj) = amean(jj) + rv(jj)
      enddo      
      
      zb = matmul(zi,betaa)
      
      end subroutine            
   

      subroutine gibbs_taui(iseed)
      implicit real*8 (a-h, o-z)
	integer, parameter :: ni = 1180, nm = 14, np = 91
	integer, parameter :: nc = 8, nb = 3, nl = 5, nk = 5
		
      real*8 wi(ni), taui(ni), vnu            
      real*8 zb(ni), fxim(ni), fxip(ni)
      
      real*8 wstar, shape, scale, rv
            
      common /vwi/wi
      common /vtaui/taui
      common /vvnu/vnu
               
      common /vzb/zb
      common /vfxim/fxim
      common /vfxip/fxip
      
      external drngam
      
      do i = 1, ni
                                                                
          wstar = wi(i) - zb(i) - fxim(i) - fxip(i)
                    
		shape = (vnu + 1.d0)/2.d0
		scale = (vnu + wstar**2)/2.d0
          
          call rnset(iseed)
          call drngam(1,shape,rv)
          call rnget(iseed)

          taui(i) = rv/scale      
                    		                        
      enddo  
                         
      end subroutine      
        
      
      subroutine gibbs_wi(iseed)
      implicit real*8 (a-h, o-z)
	integer, parameter :: ni = 1180, nm = 14, np = 91
	integer, parameter :: nc = 8, nb = 3, nl = 5, nk = 5
		
	real*8 yi(ni)  
      real*8 wi(ni), taui(ni)      
      real*8 zb(ni), fxim(ni), fxip(ni)
      
      real*8 wmean, wsigma, trim, rv       
      
      logical la, lb
            
      common /vyi/yi
      
      common /vwi/wi
      common /vtaui/taui
                       
      common /vzb/zb
      common /vfxim/fxim
      common /vfxip/fxip
      
      external ytuvn
      
      do i = 1, ni
                                                  
          wmean = zb(i) + fxim(i) + fxip(i)
          wsigma = 1.d0/taui(i)

		trim = -wmean/dsqrt(wsigma)
          
		if (yi(i) .eq. 0.d0) then
		    la = .true. ; lb = .false.
          else
		    la = .false. ; lb = .true.
		endif

		call rnset(iseed)
		rv = ytuvn(trim, trim, la, lb, iseed)
		call rnget(iseed)
          
		wi(i) = wmean + rv*dsqrt(wsigma)                             
          	          
      enddo  
                
      end subroutine  