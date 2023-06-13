	program main
      implicit real*8 (a-h, o-z)
	integer, parameter :: ni = 1180, nm = 14, np = 91
	integer, parameter :: nc = 8, nb = 3, nl = 5, nk = 5
      integer, parameter :: initer = 20000
      integer, parameter :: niter = 100000, niter1 = 20000
	      
	real*8 cancer, base(nc-1)
      real*8 duratm(nm), duratp(np)
      real*8 zmean(nc), zstd(nc)
      real*8 pmax(np)
	
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
      
      integer ngm(nl), ngp(nk)
      real*8 zb(ni), temp, summ, sump
      
      real*8 fxim(ni), deltam(nm,nb)
      real*8 xxim(nm*nb), ddeltam(nm*nb)
      
      real*8 fxip(ni), deltap(np,nb)
      real*8 xxip(np*nb), ddeltap(np*nb)
      
      real*8 probm(nm,nl), probp(np,nk)
                 
	real*8 mean_betaa(nc)
      real*8 mean_deltam(nm,nb)     
      real*8 mean_deltap(np,nb)   
	real*8 bardic, dicbar, pd, dic 
      real*8 pp, mean_pp(ni)
      
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
      
      external gen_dic
      
    	open(unit = 5, file = 'nDustData.txt') 
     	open(unit = 6, file = 'nPesticideF_Initial.txt') 
     
    	open(unit = 11, file = 'nPesticideF_Output1.txt') 
    	open(unit = 12, file = 'nPesticideF_Output2.txt') 
    	open(unit = 13, file = 'nPesticideF_Output3.txt') 
    	open(unit = 14, file = 'nPesticideF_Output4.txt') 
    	open(unit = 15, file = 'nPesticideF_Output5.txt') 
    	open(unit = 16, file = 'nPesticideF_Output6.txt') 
    	open(unit = 17, file = 'nPesticideF_Output7.txt') 
    	open(unit = 18, file = 'nPesticideF_Output8.txt') 
    	open(unit = 19, file = 'nPesticideF_Output9.txt') 
    	open(unit = 20, file = 'nPesticideF_Output10.txt') 
            
      iseed = 9999999
                
      do jj = 1, nc
          zmean(jj) = 0.d0; zstd(jj) = 0.d0
      enddo   
	do ii = 1, ni
		             
	    read(5,*) i, cancer, base, duratm
                    
          jcount = 0
          do j1 = 1, nm-1
              do j2 = j1+1,nm
                  jcount = jcount + 1
                  duratp(jcount) = duratm(j1)*duratm(j2)
              enddo
          enddo
          
	    yi(i) = cancer
          
          zi(i,1) = 1.d0
          do jj = 2, nc
              zi(i,jj) = base(jj-1)
              zmean(jj) = zmean(jj) + zi(i,jj)/dfloat(ni)
              zstd(jj) = zstd(jj) + zi(i,jj)**2
          enddo
          
          do j = 1, nm
              xim(i,j,1) = duratm(j)
              xim(i,j,2) = duratm(j)**2
              xim(i,j,3) = duratm(j)**3
          enddo

          do j = 1, np
              xip(i,j,1) = duratp(j)
              xip(i,j,2) = duratp(j)**2
              xip(i,j,3) = duratp(j)**3
          enddo
          
      enddo 
               
      do jj = 2, nc  
          temp = (zstd(jj) - dfloat(ni)*zmean(jj)**2)
     +            /dfloat(ni-1)
          zstd(jj) = dsqrt(temp)
      enddo

      do i = 1, ni
          do jj = 2, nc  
              zi(i,jj) = (zi(i,jj) - zmean(jj))/zstd(jj)
          enddo
      enddo
      
      do j = 1, nm
          pmax(j) = -10000.d0
          do i = 1, ni
              if (pmax(j) .le. xim(i,j,1)) then
                  pmax(j) = xim(i,j,1)
              endif      
          enddo
      enddo
      
      do i = 1, ni
          do j = 1, nm
              duratm(j) = xim(i,j,1)/pmax(j)
              xim(i,j,1) = duratm(j)
              xim(i,j,2) = duratm(j)**2
              xim(i,j,3) = duratm(j)**3
          enddo
      enddo
      
      do j = 1, np
          pmax(j) = -10000.d0
          do i = 1, ni
              if (pmax(j) .le. xip(i,j,1)) then
                  pmax(j) = xip(i,j,1)
              endif      
          enddo
      enddo
      
      do i = 1, ni
          do j = 1, np
              duratp(j) = xip(i,j,1)/pmax(j)
              xip(i,j,1) = duratp(j)
              xip(i,j,2) = duratp(j)**2
              xip(i,j,3) = duratp(j)**3
          enddo
      enddo
      
c     set hyper-parameters
      sigma_betaa = 1000.d0
      sigma_thetam = 1000.d0
      sigma_thetap = 1000.d0
      a0 = 1.0d0 ; b0 = 0.1d0
      a1 = 1.0d0 ; b1 = 0.1d0
      vnu = 7.d0 
                       
c     set initial values    
                       
      read(6,*) betaa
      read(6,*) alpham
      read(6,*) zetam
      read(6,*) gammam
      read(6,*) thetam
      read(6,*) alphap
      read(6,*) zetap
      read(6,*) gammap
      read(6,*) thetap                 
                        
      do i = 1, ni
          wi(i) = 0.d0
          taui(i) = 1.d0
      enddo                     
      
c     Run Gibbs
      call rnset(iseed)
      
      do j = 1, nm
          call rnset(iseed)
          call rnund(1,nl,igm)
          call rnget(iseed)              
          gm(j) = dfloat(igm)
      enddo
      
      do j = 1, np
          call rnset(iseed)
          call rnund(1,nk,igp)
          call rnget(iseed)              
          gp(j) = dfloat(igp)
      enddo
                  
      zb = matmul(zi,betaa)
      
      do j = 1, nm
          l = nint(gm(j))
          do jj = 1, nb
              deltam(j,jj) = alpham(jj,l)               
          enddo
      enddo
      do i = 1, ni
          xxim = reshape(xim(i,:,:),(/nm*nb/))
          ddeltam = reshape(deltam,(/nm*nb/))
          fxim(i) = dot_product(xxim,ddeltam)
      enddo

      do j = 1, np
          l = nint(gp(j))
          do jj = 1, nb
              deltap(j,jj) = alphap(jj,l)               
          enddo
      enddo
      do i = 1, ni
          xxip = reshape(xip(i,:,:),(/np*nb/))
          ddeltap = reshape(deltap,(/np*nb/))
          fxip(i) = dot_product(xxip,ddeltap)
      enddo
      
      icount = 0
      do ir = 1, initer
              
          call gibbs(iseed)

          do l = 1, nl
              ngm(l) = 0
              do j = 1, nm
                  if (nint(gm(j)) .eq. l) then
                      ngm(l) = ngm(l) + 1
                  endif
              enddo
          enddo

          do l = 1, nk
              ngp(l) = 0
              do j = 1, np
                  if (nint(gp(j)) .eq. l) then
                      ngp(l) = ngp(l) + 1
                  endif
              enddo
          enddo
          
          write(*,1) ir, ngm, ngp
          write(*,2) ir, betaa
          do l = 1, nl
              write(*,2) ir, (alpham(jj,l),jj=1,nb)
          enddo
          write(*,2) ir, wgtm
          write(*,2) ir, zetam, gammam
          do l = 1, nk
              write(*,2) ir, (alphap(jj,l),jj=1,nb)
          enddo
          write(*,2) ir, wgtp
          write(*,2) ir, zetap, gammap
          write(*,*)
                    
          if (ir - ir/5*5 .eq. 1) then
        
              icount = icount + 1
              
              write(11,1) icount, ngm, ngp
              write(12,3) icount, betaa
              write(13,3) icount, alpham, zetam, gammam
              write(14,3) icount, wgtm, thetam
              write(15,3) icount, gm
              write(16,3) icount, alphap, zetap, gammap
              write(17,3) icount, wgtp, thetap
              write(18,3) icount, gp
                                                              
          endif
                                   
      enddo    	  
               
      do jj = 1, nc
	    mean_betaa(jj) = 0.d0          
      enddo
            
      do j = 1, nm
          do jj = 1, nb
	        mean_deltam(j,jj) = 0.d0          
          enddo
      enddo

      do j = 1, np
          do jj = 1, nb
	        mean_deltap(j,jj) = 0.d0          
          enddo
      enddo
      
      do i = 1, ni
          mean_pp(i) = 0.d0
      enddo                     
      
      icount = 0
      do ir = 1, niter
              
          call gibbs(iseed)

          do l = 1, nl
              ngm(l) = 0
              do j = 1, nm
                  if (nint(gm(j)) .eq. l) then
                      ngm(l) = ngm(l) + 1
                  endif
              enddo
          enddo

          do l = 1, nk
              ngp(l) = 0
              do j = 1, np
                  if (nint(gp(j)) .eq. l) then
                      ngp(l) = ngp(l) + 1
                  endif
              enddo
          enddo
                            
          write(*,1) ir, ngm, ngp
          write(*,2) ir, betaa
          do l = 1, nl
              write(*,2) ir, (alpham(jj,l),jj=1,nb)
          enddo
          write(*,2) ir, wgtm
          write(*,2) ir, zetam, gammam
          do l = 1, nk
              write(*,2) ir, (alphap(jj,l),jj=1,nb)
          enddo
          write(*,2) ir, wgtp
          write(*,2) ir, zetap, gammap
          write(*,*)
                    
          if (ir - ir/5*5 .eq. 1) then
        
              icount = icount + 1
              
              write(11,1) icount, ngm, ngp
              write(12,3) icount, betaa
              write(13,3) icount, alpham, zetam, gammam
              write(14,3) icount, wgtm, thetam
              write(15,3) icount, gm
              write(16,3) icount, alphap, zetap, gammap
              write(17,3) icount, wgtp, thetap
              write(18,3) icount, gp
                         
              do jj = 1, nc
	            mean_betaa(jj) = mean_betaa(jj) 
     +                           + betaa(jj)/dfloat(niter1)
              enddo
                            
              do j = 1, nm
                  do jj = 1, nb
                      mean_deltam(j,jj) = mean_deltam(j,jj) 
     +                                  + deltam(j,jj)
     +                                    /dfloat(niter1)
                  enddo
              enddo

              do j = 1, np
                  do jj = 1, nb
                      mean_deltap(j,jj) = mean_deltap(j,jj) 
     +                                  + deltap(j,jj)
     +                                    /dfloat(niter1)
                  enddo
              enddo
              
              bardic = bardic 
     +               - 2.d0*gen_dic(betaa,deltam,deltap)
     +                 /dfloat(niter1)
                          

              do i = 1, ni
              
                  pp = dtdf(zb(i) + fxim(i) + fxip(i),vnu)

                  mean_pp(i) = mean_pp(i) + pp/dfloat(niter1)
                  
              enddo
                      
          endif
                                   
      enddo    	  
      
      call rnget(iseed)
                     
      dicbar = -2.d0*gen_dic(mean_betaa,mean_deltam,mean_deltap)
      pd = bardic - dicbar
      dic = dicbar + 2.0d0*pd 

      write(19,4) dicbar, pd, dic

      do i = 1, ni
          write(20,3) i, yi(i), mean_pp(i) 
      enddo
            
    1 format(1000i5)
    2 format(i5,1000f10.5)
    3 format(i5,1000f20.10)
    4 format(1000f20.10)
            
      end program
                     
      real*8 function gen_dic(betaa,deltam,deltap)
      implicit real*8 (a-h, o-z)
	integer, parameter :: ni = 1180, nm = 14, np = 91
	integer, parameter :: nc = 8, nb = 3, nl = 5, nk = 5
		
	real*8 yi(ni), zi(ni,nc)    
	real*8 xim(ni,nm,nb), xip(ni,np,nb)    
      real*8 vnu, betaa(nc)
      real*8 deltam(nm,nb), deltap(np,nb)     
      
      real*8 zb, fxim, fxip, cdf, dic
      real*8 xxim(nm*nb), ddeltam(nm*nb)
      real*8 xxip(np*nb), ddeltap(np*nb)
                        
      common /vyi/yi
      common /vzi/zi
      common /vxim/xim
      common /vxip/xip
      
      common /vvnu/vnu
                  
      external dtdf
            
      dic = 0.d0
      do i = 1, ni
                   
          zb = dot_product(zi(i,:),betaa)
                    
          xxim = reshape(xim(i,:,:),(/nm*nb/))
          ddeltam = reshape(deltam,(/nm*nb/))
          fxim = dot_product(xxim,ddeltam)

          xxip = reshape(xip(i,:,:),(/np*nb/))
          ddeltap = reshape(deltap,(/np*nb/))
          fxip = dot_product(xxip,ddeltap)
          
          cdf = dtdf(zb + fxim + fxip, vnu)
                    
          if (yi(i) .eq. 1.d0) then
                                            
              dic = dic + dlog(cdf)
                      
          else

              dic = dic + dlog(1.d0 - cdf)
                      
          endif
              
      enddo
                        
	gen_dic = dic

	end function          
      
 	include 'nPesticideF_Gibbs.f'
	include 'tnorm.f'
	include 'optim1.f'
      include 'rGiGDist.f'                        
      include 'hpd.f'                        
            