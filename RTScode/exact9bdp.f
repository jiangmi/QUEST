
c  12/15/04 PUT IN D-WAVE WITHOUT SIGN.
c  (x) AGREES WITH SIGN AT HIGH T WHEN <sgn>=1.0
c  (x) CHECKS AGAINST OLD RESULTS AT LOW T.

c	DOUBLE PRECISION
c	rand-->ran2
c	MIXED REAL/INTEGER COMMONS ELIMINATED
c	INCLUDE GLOBAL FLIP MOVE.
c********************************EXACT9.FOR**********************
 
c        The main program for the two-dimensional Hubbard model with
c        positive U.
 
        program exact9dp
        implicit none
        integer n,l,volume,toff
        parameter(n=4,l=48)
        parameter(volume=n*n*l,toff=n*n)
 
        integer warms,sweeps,msr,nwrap
        integer iran,numtry
        common/integers/warms,sweeps,msr,nwrap, iran
 
         real*8 t,u,mu,dtau,expmu,gam,lambda,tdtau,dens
         common/couple/t,u,mu,dtau,expmu,gam,lambda,tdtau,dens
 
c        index vectors for indirect addressing and phase vector
        integer xplus(0:toff-1),xminus(0:toff-1)
        integer yplus(0:toff-1),yminus(0:toff-1)
 
        common/ivectors/xplus,xminus,yplus,yminus
 
c        vectors:
 
        real*8 hub(0:volume-1)
        real*8 vup(0:volume-1),vdn(0:volume-1)
        real*8 mphase(0:toff-1)
 
c        common block for vectors
        common/vectors/hub,vup,vdn, mphase
        real*8 twopi
        parameter(twopi=6.283185307179586)
 
        integer i,k,accept2,reject2,accept,reject,wraps,tausk
        parameter (tausk=20)
        integer j1,j2,mx,my,lx,ly,ti,mpx,mpy,kx,ky,kmx,kmy
        real*8 time,aveval,errval,aveval2,errval2
        complex*16 cveval,crrval,cveval2,crrval2
 
        real*8 bpairv(10,9),bpairl(10,9,0:l),gpairl(10,9,0:l)
        complex*16 bpairw(10,9,0:l),gpairw(10,9,0:l)
        real*8 waves(-1:1,-1:1,9)
        data waves/0.0, 0.0, 0.0, 0.0,1.0,0.0,0.0, 0.0, 0.0,
     1             0.0, 0.5, 0.0, 0.5,0.0,0.5,0.0, 0.5, 0.0,
     2             0.0,-0.5, 0.0, 0.5,0.0,0.5,0.0,-0.5, 0.0,
     3             0.5, 0.0, 0.5, 0.0,0.0,0.0,0.5, 0.0, 0.5,
     4            -0.5, 0.0, 0.5, 0.0,0.0,0.0,0.5, 0.0,-0.5,
     5             0.0, 0.0, 0.0,-0.5,0.0,0.5,0.0, 0.0, 0.0,
     6             0.0,-0.5, 0.0, 0.0,0.0,0.0,0.0, 0.5, 0.0,
     7            -0.5, 0.0, 0.0, 0.0,0.0,0.0,0.0, 0.0, 0.5,
     8             0.0, 0.0,-0.5, 0.0,0.0,0.0,0.5, 0.0, 0.0/
 
c       *****
 
        real*8 sgnup,sgndn,sgn,nup,ndn,nud,ke
        real*8 saf,saf2,sferro,sfer2,grfun(0:n/2,0:n/2),safsq,safsq2
        real*8 spinxx(0:n/2,0:n/2),spinzz(0:n/2,0:n/2)
        real*8 pairmat(-1:1,-1:1,-1:1,-1:1),errrat,asgnp
        integer nmeas0,redo,noredo,dopair,nmeasp,numpair
        real*8 asgnup,asgndn,asgn,anup,andn,anud,ake,asafsq,asafsq2
        real*8 asaf,asaf2,asferro,asfer2,agrfun(0:n/2,0:n/2)
        real*8 aspinxx(0:n/2,0:n/2),aspinzz(0:n/2,0:n/2)
        real*8 apairmat(-1:1,-1:1,-1:1,-1:1)
        real*8 den(0:n/2,0:n/2,0:1),aden(0:n/2,0:n/2,0:1)
        common/m0var/ sgnup,sgndn,sgn,nup,ndn,nud,ke,safsq,safsq2,
     1       saf,saf2,sferro,sfer2,grfun,
     1       spinxx,spinzz,aspinxx,aspinzz,asafsq,asafsq2,
     2       asgnup,asgndn,asgn,anup,andn,anud,ake,den,aden,
     3       asaf,asaf2,asferro,asfer2,agrfun,errrat,asgnp
     4       ,pairmat,apairmat
        common/m0vari/dopair,nmeasp,numpair,nmeas0,redo,noredo
        real*8 gnl(0:n/2,0:n/2,0:l),agnl(0:n/2,0:n/2,0:l)
        real*8 chinl(0:n/2,0:n/2,0:l),achinl(0:n/2,0:n/2,0:l)
        real*8 pairsus(-1:1,-1:1,-1:1,-1:1,0:l)
        real*8 apairsus(-1:1,-1:1,-1:1,-1:1,0:l),asgnt
        integer nmeast
        common/mtauvar/gnl,agnl,chinl,achinl,pairsus,
     1       apairsus,asgnt,nmeast
 
        real*8 gql(0:n/2,0:n/2,0:l)
        complex*16 gnw(0:n/2,0:n/2,0:l),gqw(0:n/2,0:n/2,0:l)
        complex*16 sigma(0:n/2,0:n/2,0:l)
        real*8 chiql(0:n/2,0:n/2,0:l)
        complex*16 chinw(0:n/2,0:n/2,0:l),chiqw(0:n/2,0:n/2,0:l)
 
        real*8 bgnl(0:n/2,0:n/2,0:l),bgql(0:n/2,0:n/2,0:l)
        complex*16 bgnw(0:n/2,0:n/2,0:l),bgqw(0:n/2,0:n/2,0:l)
        complex*16 bsigma(0:n/2,0:n/2,0:l)
        real*8 bchinl(0:n/2,0:n/2,0:l),bchiql(0:n/2,0:n/2,0:l)
        complex*16 bchinw(0:n/2,0:n/2,0:l),bchiqw(0:n/2,0:n/2,0:l)
 
        real*8 sgnl(0:n/2,0:n/2,0:l),sgql(0:n/2,0:n/2,0:l)
        complex*16 sgnw(0:n/2,0:n/2,0:l),sgqw(0:n/2,0:n/2,0:l)
        complex*16 ssigma(0:n/2,0:n/2,0:l)
        real*8 schinl(0:n/2,0:n/2,0:l),schiql(0:n/2,0:n/2,0:l)
        complex*16 schinw(0:n/2,0:n/2,0:l),schiqw(0:n/2,0:n/2,0:l)
 
        real*8 bsgnup(10),bsgndn(10),bsgn(10),bnup(10),bndn(10)
        real*8 bntot(10),bnud(10),bke(10),benergy(10)
        real*8 bsaf(10),bsaf2(10),bsferro(10),bsfer2(10)
        real*8 bgrfun(10,0:n/2,0:n/2),bsafsq(10),bsafsq2(10)
        real*8 bden(10,0:n/2,0:n/2,0:1)
        real*8 bspinxx(10,0:n/2,0:n/2), bspinzz(10,0:n/2,0:n/2)
 
        real*8 pairgsus(-1:1,-1:1,-1:1,-1:1,0:l)
        real*8 gmatup(0:n*n,0:n*n),gmatdn(0:n*n,0:n*n),epsk
	real*8 detup,detdn
        integer nmax
 
c        write header
        write(76,*) 'Version exact9bdp'
 
        write(76,*) ' '
        write(76,*) 'n=',n,'  l=',l
        write(76,*) ' '
        write(76,*) 'afeps = ',0.0
 
c        set up the phase vectors
        call phaseset()
 
        call readin()
 
c        set up the index vectors
        call indexset()
 
          mu = 0.d0
        expmu=dexp(dtau*mu)
        call setvup()
 
        call autoset()
 
        call ranlat()
        write (6,*) 'INSTEAD, READ IN MU ITSELF'
        read  (5,*)  mu
        expmu=dexp(dtau*mu)
        write (76,*) 'Using mu = ',mu
        write (6,*) 'Using mu = ',mu
        write (6,*) 'ENTER NUMTRY'
        read  (5,*)  numtry
        write (76,*) 'Using numtry = ',numtry
 
        nmax = l / 2
        accept = 0
        reject = 0
        accept2= 0
        reject2= 0
        do 141 ti = 0, l
         do 141 j1 = 0, n/2
          do 141 j2 = 0, n/2
           bgnl(j1,j2,ti) = 0.d0
           bgnw(j1,j2,ti) = 0.d0
           bgql(j1,j2,ti) = 0.d0
           bgqw(j1,j2,ti) = 0.d0
           bsigma(j1,j2,ti) = 0.d0
           bchinl(j1,j2,ti) = 0.d0
           bchinw(j1,j2,ti) = 0.d0
           bchiql(j1,j2,ti) = 0.d0
           bchiqw(j1,j2,ti) = 0.d0
           sgnl(j1,j2,ti) = 0.d0
           sgnw(j1,j2,ti) = 0.d0
           sgql(j1,j2,ti) = 0.d0
           sgqw(j1,j2,ti) = 0.d0
           ssigma(j1,j2,ti) = 0.d0
           schinl(j1,j2,ti) = 0.d0
           schinw(j1,j2,ti) = 0.d0
           schiql(j1,j2,ti) = 0.d0
           schiqw(j1,j2,ti) = 0.d0
141     continue
        call setvup()
 
c        perform warmup sweeps
        wraps = nwrap
        redo = 0
        noredo = 0
        call getgp(vup,0,gmatup,sgnup,detup)
        if(mu .ne. 0.d0) then
          call getgp(vdn,0,gmatdn,sgndn,detdn)
          sgndn = sgnup
        endif
        do 10 i=1,warms
          write(6,*)'Starting warmup sweep ',i
          call sweep(gmatup,gmatdn,accept,reject,wraps)
          call sweep2(gmatup,gmatdn,accept2,reject2,wraps,numtry)
10      continue
        write(76,*)'after warmups, accept ratio is ',
     1                float(accept)/(accept+reject)
        if (numtry.ne.0) then
        write(76,*)'after warmups, accept2 ratio is ',
     1                float(accept2)/(accept2+reject2)
        endif
        write(76,*)'gamma is ',gam
        write(76,*)'redo ratio is ',
     1                float(redo)/(redo+noredo)
 
c        perform measurement sweeps
        call setvup()
        call zeroas()
        do 20 i=1,sweeps
          write(6,*)'Starting measurement sweep ',i
          if(mod(i,10) .eq. 0) then
            write(6,*)'accept, redo ratios are ',
     1                float(accept)/(accept+reject),
     2                float(redo)/(noredo+redo)
            if (numtry.ne.0) then
                 write(6,*)'accept2 ratio is ',
     1                float(accept2)/(accept2+reject2)
             endif
          endif
          call sweep(gmatup,gmatdn,accept,reject,wraps)
          call sweep2(gmatup,gmatdn,accept2,reject2,wraps,numtry)
          if(mod(i,tausk) .eq. 0)then
            call meastau
          endif
          if(mod(i,sweeps/10) .eq. 0) then
            write(76,*)'Finished measurement sweep ',i
            if(dopair .eq. 0)nmeasp = 1
            write(76,9045) asgn/nmeas0,asgnp/nmeasp,
     1            float(accept)/(accept+reject),
     2            float(redo)/(noredo+redo)
9045        format('asgn, asgnp: ',2(f8.3,' '),
     1          ';accept,redo ratios: ',2('  ',f9.4))
            if (numtry.ne.0) then
                 write(76,*)'accept2 ratio is ',
     1                float(accept2)/(accept2+reject2)
             endif
            k = (i*10)/sweeps
c In the following lines, we have replaced nmeas0 with asgn
c for a different error estimate.
            if(asgn .eq. 0.d0) then
              write(76,*)'adding .01 to asgn'
              asgn = asgn + .01
            endif
            bnup(k)= anup / asgn
            bndn(k)= andn / asgn
            bntot(k)= (anup+andn) / asgn
            bsaf(k)= asaf / asgn
            bsafsq(k)= dsqrt(dabs(asafsq / asgn))
            bsferro(k)= asferro / asgn
            bsfer2(k)= asfer2 / asgn
            bsaf2(k)= asaf2 / asgn
            bke(k)= ake / asgn
            benergy(k)= (ake+u*anud) / asgn
            bnud(k)= anud / asgn
c
            bsafsq2(k)= dsqrt(dabs(asafsq2) / nmeas0)
            bsgnup(k) = asgnup / nmeas0
            bsgndn(k) = asgndn / nmeas0
            bsgn(k) = asgn / nmeas0
            do 41 j1 = 0, n/2
             do 41 j2 = 0, n/2
              bgrfun(k,j1,j2) = agrfun(j1,j2) / asgn
              bden(k,j1,j2,0) = aden(j1,j2,0) / asgn
              bden(k,j1,j2,1) = aden(j1,j2,1) / asgn
              bspinxx(k,j1,j2) = aspinxx(j1,j2) / asgn
              bspinzz(k,j1,j2) = aspinzz(j1,j2) / asgn
c Set gnl and chinl for this tenth of the run.
              do 41 ti = 0, l
               gnl(j1,j2,ti) = agnl(j1,j2,ti)/asgnt
               gql(j1,j2,ti) = 0.d0
               chinl(j1,j2,ti) = achinl(j1,j2,ti)/asgnt
               chiql(j1,j2,ti) = 0.d0
41          continue
            call ftntok(gnl,gql,n,l)
            call ftntok(chinl,chiql,n,l)
            call ftltow(gnl,gnw,n,l,dtau,0,nmax)
            call ftltow(gql,gqw,n,l,dtau,0,nmax)
            call ftltow(chinl,chinw,n,l,dtau,1,nmax)
            call ftltow(chiql,chiqw,n,l,dtau,1,nmax)
          
            do 586 mx = 0, n/2
             do 586 my = 0, n/2
              do 586 ti = 0, l
               bgnl(mx,my,ti) = bgnl(mx,my,ti) + gnl(mx,my,ti)
               bgql(mx,my,ti) = bgql(mx,my,ti) + gql(mx,my,ti)
               bchinl(mx,my,ti) = bchinl(mx,my,ti)+chinl(mx,my,ti)
               bchiql(mx,my,ti) = bchiql(mx,my,ti)+chiql(mx,my,ti)
               sgnl(mx,my,ti) = sgnl(mx,my,ti) + gnl(mx,my,ti)**2
               sgql(mx,my,ti) = sgql(mx,my,ti) + gql(mx,my,ti)**2
               schinl(mx,my,ti)=schinl(mx,my,ti)+chinl(mx,my,ti)**2
               schiql(mx,my,ti)=schiql(mx,my,ti)+chiql(mx,my,ti)**2
              if(ti .le. nmax)then
               bgnw(mx,my,ti) = bgnw(mx,my,ti) + gnw(mx,my,ti)
               bgqw(mx,my,ti) = bgqw(mx,my,ti) + gqw(mx,my,ti)
               bchinw(mx,my,ti) = bchinw(mx,my,ti)+chinw(mx,my,ti)
               bchiqw(mx,my,ti) = bchiqw(mx,my,ti)+chiqw(mx,my,ti)
               sgnw(mx,my,ti) = sgnw(mx,my,ti) + cmplx(dreal
     1                 (gnw(mx,my,ti))**2,dimag(gnw(mx,my,ti))**2)
               sgqw(mx,my,ti) = sgqw(mx,my,ti) + cmplx(dreal
     1                 (gqw(mx,my,ti))**2,dimag(gqw(mx,my,ti))**2)
               schinw(mx,my,ti) = schinw(mx,my,ti) + cmplx(dreal
     1              (chinw(mx,my,ti))**2,dimag(chinw(mx,my,ti))**2)
               schiqw(mx,my,ti) = schiqw(mx,my,ti) + cmplx(dreal
     1              (chiqw(mx,my,ti))**2,dimag(chiqw(mx,my,ti))**2)
               epsk = -2.d0*t*(cos(twopi*mx/n)+cos(twopi*my/n))
               sigma(mx,my,ti)=cmplx(-epsk+mu,(ti+0.5)*twopi/l/dtau)
     1                             + 1.d0/gqw(mx,my,ti) 
               bsigma(mx,my,ti) = bsigma(mx,my,ti) + sigma(mx,my,ti)
               ssigma(mx,my,ti) = ssigma(mx,my,ti) + cmplx(dreal
     1              (sigma(mx,my,ti))**2,dimag(sigma(mx,my,ti))**2)
              endif
586         continue
            if(dopair .eq. 1) then
              if(asgnp .eq. 0.d0) then
                write(76,*)'adding .01 to asgnp'
                asgnp = asgnp + .01
              endif
            do 8379 mx = -1, 1
            do 8379 my = -1, 1
             do 8379 mpx = -1, 1
             do 8379 mpy = -1, 1
              do 8356 ti = 0,l
8356           pairgsus(mx,my,mpx,mpy,ti) = 0.d0
              do 8379 lx = 0, n-1
              do 8379 ly = 0, n-1
               kx = min(lx,n-lx)
               ky = min(ly,n-ly)
               kmx = mod(n+n+lx-(mx-mpx),n)
               kmx = min(kmx,n-kmx)
               kmy = mod(n+n+ly-(my-mpy),n)
               kmy = min(kmy,n-kmy)
               do 8379 ti = 0,l
              pairgsus(mx,my,mpx,mpy,ti)=pairgsus(mx,my,mpx,mpy,ti)
     1          +gnl(kmx,kmy,ti)*gnl(kx,ky,ti)
8379        continue
            do 2378 j1 = 1,9
             bpairv(k,j1) = 0.d0
             do 2338 ti = 0, l
               gpairl(k,j1,ti) = 0.d0
2338           bpairl(k,j1,ti) = 0.d0
             do 2378 mx = -1, 1
              do 2378 my = -1, 1
               do 2378 mpx = -1, 1
                do 2378 mpy = -1, 1
                  bpairv(k,j1)=bpairv(k,j1)+waves(mx,my,j1)*
     1            apairmat(mx,my,mpx,mpy)*waves(mpx,mpy,j1)/asgnp
                 do 2378 ti = 0, l
                bpairl(k,j1,ti)=bpairl(k,j1,ti)+waves(mx,my,j1)*
     1           apairsus(mx,my,mpx,mpy,ti)*waves(mpx,mpy,j1)/asgnt
                gpairl(k,j1,ti)=gpairl(k,j1,ti)+waves(mx,my,j1)*
     1            pairgsus(mx,my,mpx,mpy,ti)*waves(mpx,mpy,j1)
2378          continue
              call ftltowp(bpairl,bpairw,l,k,dtau,nmax)
              call ftltowp(gpairl,gpairw,l,k,dtau,nmax)
            endif                                  
            call zeroas()
          endif
20      continue
 
        write(76,*)'At end, redo ratio is ',float(redo)/(redo+noredo)
        write(76,*)'gamma is ',gam
       write(76,*)'Acceptance ratio = ',float(accept)/(accept+reject)
        call geterr(bsgnup,aveval,errval)
        write(76,*) 'Average up sign =',aveval,' +- ',errval
        call geterr(bsgndn,aveval,errval)
        write(76,*) 'Average dn sign =',aveval,' +- ',errval
        call geterr(bsgn,aveval,errval)
        write(76,*) 'Average total sign =',aveval,' +- ',errval
        call geterr(bntot,aveval,errval)
        write(76,*) 'Average density =',
     1         aveval,' +- ',errval
        call geterr(bnup,aveval,errval)
        write(76,*) 'Average up occupancy =',
     1         aveval,' +- ',errval
        call geterr(bndn,aveval,errval)
        write(76,*) 'Average dn occupancy =',
     1         aveval,' +- ',errval
        call geterr(benergy,aveval,errval)
        write(76,*) 'Average Energy =',
     1         aveval,' +- ',errval
        call geterr(bke,aveval,errval)
        write(76,*) 'Average Kinetic Energy =',
     1         aveval,' +- ',errval
        call geterr(bnud,aveval,errval)
        write(76,*) 'Average Nup*Ndn =',
     1         aveval,' +- ',errval
        call geterr(bsaf,aveval,errval)
        write(76,*) 'AF correlation function (xx) = ',
     1         aveval,' +- ',errval
        call geterr(bsaf2,aveval,errval)
        write(76,*) 'AF correlation function (zz)= ',
     1         aveval,' +- ',errval
        call geterr(bsferro,aveval,errval)
        write(76,*) 'Ferro correlation function(xx)= ',
     1         aveval,' +- ',errval
        call geterr(bsfer2,aveval,errval)
        write(76,*) 'Ferro correlation function(zz)= ',
     1         aveval,' +- ',errval
        write(76,*)'Green''s function:'
        do 677 j1 = 0, n/2
          do 677 j2 = 0, n/2
            call geterr(bgrfun(1,j1,j2),aveval,errval)
            if(j2 .ge. j1)
     1        write(76,*)j1,j2,aveval,' +- ',errval
            grfun(j1,j2) = aveval
677     continue
        write(76,*)'density-density correlation fn: (up-up,up-dn)'
        do 977 j1 = 0, n/2
          do 977 j2 = j1, n/2
            call geterr(bden(1,j1,j2,0),aveval,errval)
            call geterr(bden(1,j1,j2,1),aveval2,errval2)
            write(76,1987)j1,j2,aveval,errval,aveval2,errval2
1987        format(2i4,2('    ',f12.6,' +- ',f12.6))
977     continue
        write(76,*)'zz Spin correlation function:'
        do 687 j1 = 0, n/2
          do 687 j2 = j1, n/2
            call geterr(bspinzz(1,j1,j2),aveval,errval)
            write(76,*)j1,j2,aveval,' +- ',errval
687     continue
        write(76,*)'xx Spin correlation function:'
        do 697 j1 = 0, n/2
          do 697 j2 = j1, n/2
            call geterr(bspinxx(1,j1,j2),aveval,errval)
            write(76,*)j1,j2,aveval,' +- ',errval
697     continue
        call geterr(bsafsq,aveval,errval)
        write(76,*) 'RMS AF correlation function (xx) = ',
     1         aveval,' +- ',errval
        call geterr(bsafsq2,aveval,errval)
        write(76,*) 'RMS AF correlation function (zz) = ',
     1         aveval,' +- ',errval
        write(76,*)' '
 
        write(78,*)'G(nx,ny,ti):'
1234      format(i5,f14.6,' +- ',f14.6)
        do 987 j1 = 0, n/2
         do 987 j2 = j1, n/2
          write(78,*)'nx = ',j1,' ny = ',j2
          do 987 ti = 0,l
          call geterr2(bgnl(j1,j2,ti),sgnl(j1,j2,ti),aveval,errval)
          write(78,1234)ti,-aveval,errval
987     continue
 
        write(78,*)'G(qx,qy,ti):'
        do 9871 j1 = 0, n/2
         do 9871 j2 = j1, n/2
          write(78,*)'qx = ',j1,' qy = ',j2
          do 9871 ti = 0,l
          call geterr2(bgql(j1,j2,ti),sgql(j1,j2,ti),aveval,errval)
          write(78,1234)ti,-aveval,errval
9871     continue
 
        write(78,*)'G(nx,ny,omega), omega = (n+.5) 2 pi T :'
        do 9872 j1 = 0, n/2
         do 9872 j2 = j1, n/2
          write(78,*)'nx = ',j1,' ny = ',j2
          do 9872 ti = 0,nmax
          call geterr2(dreal(bgnw(j1,j2,ti)),dreal(sgnw(j1,j2,ti)),
     1                      aveval,errval)
          call geterr2(dimag(bgnw(j1,j2,ti)),dimag(sgnw(j1,j2,ti)),
     1                      aveval2,errval2)
          write(78,1235)ti,-aveval,errval,-aveval2,errval2
1235      format(i5,'(',f12.6,' +- ',f12.6,') + i * ('
     1                  ,f12.6,' +- ',f12.6,')')
9872     continue
 
        write(78,*)'G(qx,qy,omega):'
        do 9873 j1 = 0, n/2
         do 9873 j2 = j1, n/2
          write(78,*)'qx = ',j1,' qy = ',j2
          do 9873 ti = 0,nmax
          call geterr2(dreal(bgqw(j1,j2,ti)),dreal(sgqw(j1,j2,ti)),
     1                      aveval,errval)
          call geterr2(dimag(bgqw(j1,j2,ti)),dimag(sgqw(j1,j2,ti)),
     1                      aveval2,errval2)
          write(78,1235)ti,-aveval,errval,-aveval2,errval2
9873     continue
 
        write(78,*)'SIGMA(qx,qy,omega):'
        do 9874 j1 = 0, n/2
         do 9874 j2 = j1, n/2
          write(78,*)'qx = ',j1,' qy = ',j2
          do 9875 ti = nmax,0,-1
           call geterr2(dreal(bsigma(j1,j2,ti)),dreal(ssigma(j1,j2,ti)),
     1                      aveval,errval)
         call geterr2(dimag(bsigma(j1,j2,ti)),dimag(ssigma(j1,j2,ti))
     1                      ,aveval2,errval2)
           write(78,1235)-(2*ti+1), aveval,errval,-aveval2,errval2
9875      continue
          do 9874 ti = 0,nmax
           call geterr2(dreal(bsigma(j1,j2,ti)),dreal(ssigma(j1,j2,ti)),
     1                      aveval,errval)
         call geterr2(dimag(bsigma(j1,j2,ti)),dimag(ssigma(j1,j2,ti))
     1                      ,aveval2,errval2)
           write(78,1235)(2*ti+1), aveval,errval,aveval2,errval2
9874     continue
 
        write(78,*)'chi(nx,ny,ti):'
        do 417 j1 = 0, n/2
         do 417 j2 = j1, n/2
          write(78,*)'nx = ',j1,' ny = ',j2
          do 417 ti = 0,l
        call geterr2(bchinl(j1,j2,ti),schinl(j1,j2,ti),aveval,errval)
          write(78,1234)ti, aveval,errval
417     continue
 
        write(78,*)'chi(qx,qy,ti):'
        do 4171 j1 = 0, n/2
         do 4171 j2 = j1, n/2
          write(78,*)'qx = ',j1,' qy = ',j2
          do 4171 ti = 0,l
        call geterr2(bchiql(j1,j2,ti),schiql(j1,j2,ti),aveval,errval)
          write(78,1234)ti, aveval,errval
4171     continue
 
        write(78,*)'chi(nx,ny,omega), omega = 2 n pi T :'
        do 4172 j1 = 0, n/2
         do 4172 j2 = j1, n/2
          write(78,*)'nx = ',j1,' ny = ',j2
          do 4172 ti = 0,nmax
       call geterr2(dreal(bchinw(j1,j2,ti)),dreal(schinw(j1,j2,ti)),
     1                      aveval,errval)
       call geterr2(dimag(bchinw(j1,j2,ti)),dimag(schinw(j1,j2,ti)),
     1                      aveval2,errval2)
          write(78,1235)ti, aveval,errval,aveval2,errval2
4172     continue
 
        write(78,*)'chi(qx,qy,omega):'
        do 4173 j1 = 0, n/2
         do 4173 j2 = j1, n/2
          write(78,*)'qx = ',j1,' qy = ',j2
          do 4173 ti = 0,nmax
       call geterr2(dreal(bchiqw(j1,j2,ti)),dreal(schiqw(j1,j2,ti)),
     1                      aveval,errval)
       call geterr2(dimag(bchiqw(j1,j2,ti)),dimag(schiqw(j1,j2,ti)),
     1                      aveval2,errval2)
          write(78,1235)ti, aveval,errval,aveval2,errval2
4173     continue
 
       if(dopair .eq. 1) then
1236    format(2(f12.5,' +- ',f12.5,'   ')) 
        write(76,*) 's-wave: (corr. fn, no vertex)'
        write(76,*) '        (suscept., no vertex)'
        call geterr(bpairv(1,1),aveval,errval)
        call geterr(gpairl(1,1,0),aveval2,errval2)
        write(76,1236)aveval,errval,aveval2,errval2
        call geterrc(bpairw(1,1,0),cveval,crrval)
        call geterrc(gpairw(1,1,0),cveval2,crrval2)
        write(76,1236)dreal(cveval),dreal(crrval),dreal(cveval2)
     1               ,dreal(crrval2)
 
        write(76,*) 'sx-wave: '
        call geterr(bpairv(1,2),aveval,errval)
        call geterr(gpairl(1,2,0),aveval2,errval2)
        write(76,1236)aveval,errval,aveval2,errval2
        call geterrc(bpairw(1,2,0),cveval,crrval)
        call geterrc(gpairw(1,2,0),cveval2,crrval2)
        write(76,1236)dreal(cveval),dreal(crrval),dreal(cveval2)
     1               ,dreal(crrval2)
 
        write(76,*) 'd-wave: '
        call geterr(bpairv(1,3),aveval,errval)
        call geterr(gpairl(1,3,0),aveval2,errval2)
        write(76,1236)aveval,errval,aveval2,errval2
        call geterrc(bpairw(1,3,0),cveval,crrval)
        call geterrc(gpairw(1,3,0),cveval2,crrval2)
        write(76,1236)dreal(cveval),dreal(crrval),dreal(cveval2)
     1               ,dreal(crrval2)
 
        write(76,*) 'sxx-wave: '
        call geterr(bpairv(1,4),aveval,errval)
        call geterr(gpairl(1,4,0),aveval2,errval2)
        write(76,1236)aveval,errval,aveval2,errval2
        call geterrc(bpairw(1,4,0),cveval,crrval)
        call geterrc(gpairw(1,4,0),cveval2,crrval2)
        write(76,1236)dreal(cveval),dreal(crrval),dreal(cveval2)
     1               ,dreal(crrval2)
 
        write(76,*) 'dxx-wave: '
        call geterr(bpairv(1,5),aveval,errval)
        call geterr(gpairl(1,5,0),aveval2,errval2)
        write(76,1236)aveval,errval,aveval2,errval2
        call geterrc(bpairw(1,5,0),cveval,crrval)
        call geterrc(gpairw(1,5,0),cveval2,crrval2)
        write(76,1236)dreal(cveval),dreal(crrval),dreal(cveval2)
     1               ,dreal(crrval2)
 
        write(76,*) 'px-wave: '
        call geterr(bpairv(1,6),aveval,errval)
        call geterr(gpairl(1,6,0),aveval2,errval2)
        write(76,1236)aveval,errval,aveval2,errval2
        call geterrc(bpairw(1,6,0),cveval,crrval)
        call geterrc(gpairw(1,6,0),cveval2,crrval2)
        write(76,1236)dreal(cveval),dreal(crrval),dreal(cveval2)
     1               ,dreal(crrval2)
 
        write(76,*) 'py-wave: '
        call geterr(bpairv(1,7),aveval,errval)
        call geterr(gpairl(1,7,0),aveval2,errval2)
        write(76,1236)aveval,errval,aveval2,errval2
        call geterrc(bpairw(1,7,0),cveval,crrval)
        call geterrc(gpairw(1,7,0),cveval2,crrval2)
        write(76,1236)dreal(cveval),dreal(crrval),dreal(cveval2)
     1               ,dreal(crrval2)
 
        write(76,*) 'pxy-wave: '
        call geterr(bpairv(1,8),aveval,errval)
        call geterr(gpairl(1,8,0),aveval2,errval2)
        write(76,1236)aveval,errval,aveval2,errval2
        call geterrc(bpairw(1,8,0),cveval,crrval)
        call geterrc(gpairw(1,8,0),cveval2,crrval2)
        write(76,1236)dreal(cveval),dreal(crrval),dreal(cveval2)
     1               ,dreal(crrval2)
 
        write(76,*) 'pyx-wave: '
        call geterr(bpairv(1,9),aveval,errval)
        call geterr(gpairl(1,9,0),aveval2,errval2)
        write(76,1236)aveval,errval,aveval2,errval2
        call geterrc(bpairw(1,9,0),cveval,crrval)
        call geterrc(gpairw(1,9,0),cveval2,crrval2)
        write(76,1236)dreal(cveval),dreal(crrval),dreal(cveval2)
     1               ,dreal(crrval2)
 
1237    format(i5,2('    ',f12.5,' +- ',f12.5)) 
        write(78,*)'P_d(l): (vertex,novertex)'
        do 7754 ti = 0,l
          call geterr(bpairl(1,3,ti),aveval,errval)
          call geterr(gpairl(1,3,ti),aveval2,errval2)
          write(78,1237)ti,aveval,errval,aveval2,errval2
7754    continue
 
        write(78,*)'P_d(w):(vertex)'
        do 7756 ti = 0,nmax
          call geterrc(bpairw(1,3,ti),cveval,crrval)
          write(78,1235)ti,dreal(cveval),dreal(crrval),dimag(cveval)
     1               ,dimag(crrval)
7756    continue
        write(78,*)'P_d(w):(novertex)'
        do 9756 ti = 0,nmax
          call geterrc(gpairw(1,3,ti),cveval,crrval)
          write(78,1235)ti,dreal(cveval),dreal(crrval),dimag(cveval)
     1               ,dimag(crrval)
9756    continue
          
        write(78,*)'P_sx(l):'
        do 7759 ti = 0,l
          call geterr(bpairl(1,2,ti),aveval,errval)
          call geterr(gpairl(1,2,ti),aveval2,errval2)
          write(78,1237)ti,aveval,errval,aveval2,errval2
7759    continue
 
        write(78,*)'P_sx(w):(vertex)'
        do 7856 ti = 0,nmax
          call geterrc(bpairw(1,2,ti),cveval,crrval)
          write(78,1235)ti,dreal(cveval),dreal(crrval),dimag(cveval)
     1               ,dimag(crrval)
7856    continue
        write(78,*)'P_sx(w):(novertex)'
        do 9856 ti = 0,nmax
          call geterrc(gpairw(1,2,ti),cveval,crrval)
          write(78,1235)ti,dreal(cveval),dreal(crrval),dimag(cveval)
     1               ,dimag(crrval)
9856    continue

        write(78,*)'P_s(l):'
        do 7760 ti = 0,l
          call geterr(bpairl(1,1,ti),aveval,errval)
          call geterr(gpairl(1,1,ti),aveval2,errval2)
          write(78,1237)ti,aveval,errval,aveval2,errval2
7760    continue   
 
        write(78,*)'P_s(w):(vertex)'
        do 7857 ti = 0,nmax
          call geterrc(bpairw(1,1,ti),cveval,crrval)
          write(78,1235)ti,dreal(cveval),dreal(crrval),dimag(cveval)
     1               ,dimag(crrval)
7857    continue
        write(78,*)'P_s(w):(novertex)'
        do 9857 ti = 0,nmax
          call geterrc(gpairw(1,1,ti),cveval,crrval)
          write(78,1235)ti,dreal(cveval),dreal(crrval),dimag(cveval)
     1               ,dimag(crrval)
9857    continue

      
       endif
 
        saf = 2.d0*grfun(0,0)
        do 3879 lx = 0, n-1
          do 3879 ly = 0, n-1
            kx = min(lx,n-lx)
            ky = min(ly,n-ly)
3879          saf = saf -(-1)**(lx+ly)*2*grfun(kx,ky)*grfun(kx,ky)
        write(76,*)'saf with no vertex is ',saf
 
c        call second(time)
        time=time/60.d0
        write(76,*)' '
        write(76,*) 'cpu time in minutes =',time
        stop
        end
c*******************************************
        subroutine sweep2(gmatup,gmatdn,accept,reject,wraps,numtry)
	implicit none
        integer n,l,volume,toff
        parameter(n=4,l=48)
        parameter(volume=n*n*l,toff=n*n)
        integer warms,sweeps,msr,nwrap
        integer iran
        common/integers/warms,sweeps,msr,nwrap, iran
        integer orthlen,doauto,torth
        real*8 eorth,difflim,errpam
        common/getgpam/eorth,difflim,errpam,orthlen,doauto,
     1                torth
 
        real*8 t,u,mu,dtau,expmu,gam,lambda,tdtau,dens
        common/couple/t,u,mu,dtau,expmu,gam,lambda,tdtau,dens
        real*8 ch,soc,ch4
        common/kinetic/ch,soc,ch4
        real*8 hub(0:volume-1)
        real*8 vup(0:volume-1),vdn(0:volume-1)
        real*8 mphase(0:toff-1)
 
c        common block for vectors
        common/vectors/hub,vup,vdn, mphase
 
        integer accept,reject,wraps
        integer ti,numtry,isite,itry
        real*8 gmatup(0:n*n,0:n*n),gmatdn(0:n*n,0:n*n)
        real*8 gmatupp(0:n*n,0:n*n),gmatdnp(0:n*n,0:n*n)
	real*8 detup,detdn,detupp,detdnp,sgnupp,sgndnp,ranf,ran2
 
        real*8 sgnup,sgndn,sgn,nup,ndn,nud,ke
        real*8 saf,saf2,sferro,sfer2,grfun(0:n/2,0:n/2),safsq,safsq2
        real*8 spinxx(0:n/2,0:n/2),spinzz(0:n/2,0:n/2)
        real*8 pairmat(-1:1,-1:1,-1:1,-1:1),errrat,asgnp
        integer nmeas0,redo,noredo,dopair,nmeasp,numpair
        real*8 asgnup,asgndn,asgn,anup,andn,anud,ake,asafsq,asafsq2
        real*8 asaf,asaf2,asferro,asfer2,agrfun(0:n/2,0:n/2)
        real*8 aspinxx(0:n/2,0:n/2),aspinzz(0:n/2,0:n/2)
        real*8 apairmat(-1:1,-1:1,-1:1,-1:1)
        real*8 den(0:n/2,0:n/2,0:1),aden(0:n/2,0:n/2,0:1)
        common/m0var/ sgnup,sgndn,sgn,nup,ndn,nud,ke,safsq,safsq2,
     1       saf,saf2,sferro,sfer2,grfun,
     1       spinxx,spinzz,aspinxx,aspinzz,asafsq,asafsq2,
     2       asgnup,asgndn,asgn,anup,andn,anud,ake,den,aden,
     3       asaf,asaf2,asferro,asfer2,agrfun,errrat,asgnp
     4       ,pairmat,apairmat
        common/m0vari/dopair,nmeasp,numpair,nmeas0,redo,noredo

        call getgp(vup,0,gmatup,sgnup,detup)
        call getgp(vdn,0,gmatdn,sgndn,detdn)
        sgn=sgnup*sgndn

        do 1000 itry=1,numtry
            isite=int(ran2(iran)*toff)
            do 100 ti=0,L-1
                hub(isite+ti*toff)=-hub(isite+ti*toff)
100         continue
            call setvup()
            call getgp(vup,0,gmatupp,sgnupp,detupp)
            call getgp(vdn,0,gmatdnp,sgndnp,detdnp)
            ranf=ran2(iran)
            if (ranf.le.abs(detupp*detdnp/detup/detdn) ) then
               accept=accept+1
               call matcop(gmatupp,gmatup)
               call matcop(gmatdnp,gmatdn)
               sgnup=sgnupp
               sgndn=sgndnp
               detup=detupp
               detdn=detdnp
               sgn=sgnupp*sgndnp
               else
               do 200 ti=0,L-1
                   hub(isite+ti*toff)=-hub(isite+ti*toff)
200            continue
               call setvup()
               reject=reject+1
            endif
1000    continue

	return
	end
c*******************************************
        subroutine sweep(gmatup,gmatdn,accept,reject,wraps)
        implicit none
        integer n,l,volume,toff
        parameter(n=4,l=48)
        parameter(volume=n*n*l,toff=n*n)
        integer warms,sweeps,msr,nwrap
        integer iran
        common/integers/warms,sweeps,msr,nwrap, iran
        integer orthlen,doauto,torth
        real*8 eorth,difflim,errpam
        common/getgpam/eorth,difflim,errpam,orthlen,doauto,
     1                torth
 
        real*8 t,u,mu,dtau,expmu,gam,lambda,tdtau,dens
        common/couple/t,u,mu,dtau,expmu,gam,lambda,tdtau,dens
        real*8 ch,soc,ch4
        common/kinetic/ch,soc,ch4
        real*8 hub(0:volume-1)
        real*8 vup(0:volume-1),vdn(0:volume-1)
        real*8 mphase(0:toff-1)
 
c        common block for vectors
        common/vectors/hub,vup,vdn, mphase
 
        integer accept,reject,wraps
        integer ti,ti1
        real*8 gmatup(0:n*n,0:n*n),ogmat(0:n*n,0:n*n),diffup
        real*8 gmatdn(0:n*n,0:n*n),diffdn
        real*8 accrat,redorat
 
        real*8 sgnup,sgndn,sgn,nup,ndn,nud,ke,detup,detdn
        real*8 saf,saf2,sferro,sfer2,grfun(0:n/2,0:n/2),safsq,safsq2
        real*8 spinxx(0:n/2,0:n/2),spinzz(0:n/2,0:n/2)
        real*8 pairmat(-1:1,-1:1,-1:1,-1:1),errrat,asgnp
        integer nmeas0,redo,noredo,dopair,nmeasp,numpair
        real*8 asgnup,asgndn,asgn,anup,andn,anud,ake,asafsq,asafsq2
        real*8 asaf,asaf2,asferro,asfer2,agrfun(0:n/2,0:n/2)
        real*8 aspinxx(0:n/2,0:n/2),aspinzz(0:n/2,0:n/2)
        real*8 apairmat(-1:1,-1:1,-1:1,-1:1)
        real*8 den(0:n/2,0:n/2,0:1),aden(0:n/2,0:n/2,0:1)
        common/m0var/ sgnup,sgndn,sgn,nup,ndn,nud,ke,safsq,safsq2,
     1       saf,saf2,sferro,sfer2,grfun,
     1       spinxx,spinzz,aspinxx,aspinzz,asafsq,asafsq2,
     2       asgnup,asgndn,asgn,anup,andn,anud,ake,den,aden,
     3       asaf,asaf2,asferro,asfer2,agrfun,errrat,asgnp
     4       ,pairmat,apairmat
        common/m0vari/dopair,nmeasp,numpair,nmeas0,redo,noredo
        integer lastwr,skipair
 
        do 10 ti = 0, l-1
          call multb(gmatup,vup,ti,-1)
          call multbi(gmatup,vup,ti,1)
          if(mu .ne. 0.d0)then
            call multb(gmatdn,vdn,ti,-1)
            call multbi(gmatdn,vdn,ti,1)
          endif
          if(wraps .gt. 0) then
            wraps = wraps - 1
          else
            wraps = nwrap
            ti1 = mod(ti+1,l)
            call matcop(gmatup,ogmat)
            call getgp(vup,ti1,gmatup,sgnup,detup)
            call matdif(gmatup,ogmat,diffup)
            if(mu .ne. 0.d0) then
              call matcop(gmatdn,ogmat)
              call getgp(vdn,ti1,gmatdn,sgndn,detdn)
              call matdif(gmatdn,ogmat,diffdn)
            else
              diffdn = diffup
              sgndn = sgnup
            endif
            if(diffup.gt.difflim .or. diffdn.gt.difflim)then
              redo = redo+1
            else
              noredo=noredo+1
            endif
          endif
          call swpslice(gmatup,gmatdn,sgnup,sgndn,accept,
     1                   reject,ti)
c Now do the tau=0 measurements with the free of charge Green's functions.
c Only do them every third time slice.
          if(mod(ti,12) .eq. 0)call meas0(gmatup,gmatdn)
c Measure pair correlations every skipair slices.
c Set skipair so numpair meas. are done every sweep.
          skipair = l/numpair
          if(mod(l,numpair) .ne. 0)skipair = skipair+1
          if(mod(ti,skipair).eq.0)call measpair(gmatup,gmatdn)
10      continue
 
        if(accept+reject .gt. 10000) then
          accrat = float(accept)/float(accept+reject)
          if(accrat .gt. 0.52 .or. accrat .lt. 0.48)then
            gam = gam + (accrat - 0.5)
            gam = max(0.0,real(gam))
            gam = min(1.0,real(gam))
c            write(76,*)'accrat is ',accrat,' changing gam to ',gam
            accept = 100*accrat
            reject = 100*(1.d0-accrat)
          endif
        endif
 
        if(redo .gt. 20) then
          redorat = float(redo)/float(redo+noredo)
          if(redorat.gt. errrat)then
            nwrap = nwrap - 1
            write(76,*)'reducing nwrap to ',nwrap
            redo = 0
            noredo = 1
          endif
        endif
 
        if(noredo .gt. 500) then
          redorat = float(redo)/float(redo+noredo)
          if(redorat .lt. 0.2*errrat) then
            if(nwrap .ge. lastwr)then
c Only increase wrap if nwrap has not been reduced since last increase.
              nwrap = nwrap + 2
              write(76,*)'increasing nwrap to ',nwrap
              redo = 0
              noredo = 1
              lastwr = nwrap
            endif
          endif
        endif
 
        return    
        end
c*******************************************
        subroutine swpslice(gmatup,gmatdn,sgnup,sgndn,accept,
     1                   reject,ti)
        implicit none
        integer n,l,volume,toff
        parameter(n=4,l=48)
        parameter(volume=n*n*l,toff=n*n)
        integer warms,sweeps,msr,nwrap
        integer iran
        common/integers/warms,sweeps,msr,nwrap, iran
 
         real*8 t,u,mu,dtau,expmu,gam,lambda,tdtau,dens
         common/couple/t,u,mu,dtau,expmu,gam,lambda,tdtau,dens
        real*8 hub(0:volume-1)
        real*8 vup(0:volume-1),vdn(0:volume-1)
        real*8 mphase(0:toff-1)
 
c        common block for vectors
        common/vectors/hub,vup,vdn, mphase
 
        integer accept,reject
        integer ti,nbar,nbari,j1,j2
        real*8 gmatup(0:n*n,0:n*n)
        real*8 gmatdn(0:n*n,0:n*n),sgnup,sgndn
        real*8 bvecup(0:n*n-1),bvecup2(0:n*n-1)
        real*8 bvecdn(0:n*n-1),bvecdn2(0:n*n-1)
        real*8 vnbarup,vnbardn,ratup,ratdn,rat,p,ran2,ranf
    
        do 20 nbar = 0, toff-1
          nbari = nbar + ti * toff
          vnbarup = dexp(-2.d0*lambda*hub(nbari))
          vnbardn = 1.d0 / vnbarup - 1.d0
          vnbarup = vnbarup - 1.d0
          ratup = 1.d0 + (1.d0 - gmatup(nbar,nbar))*vnbarup
          if(mu .ne. 0.d0) then
            ratdn = 1.d0 + (1.d0 - gmatdn(nbar,nbar))*vnbardn
          else
            ratdn = 1.d0 + gmatup(nbar,nbar)*vnbardn
          endif
          ranf=ran2(iran)
          rat = dabs(ratup*ratdn)
          if(rat .le. 1.d0) p = rat/(1.d0+gam*rat)
          if(rat .gt. 1.d0) p = rat/(gam+rat)
           if(p .gt. ranf) then
            accept = accept + 1
            if(ratup .lt. 0.d0)sgnup = -sgnup
            if(ratdn .lt. 0.d0)sgndn = -sgndn
            if(mu .ne. 0.d0) then
              do 30 j1 = 0, toff-1
                bvecup(j1) = -vnbarup*gmatup(nbar,j1)
30              bvecdn(j1) = -vnbardn*gmatdn(nbar,j1)
              bvecup(nbar) = bvecup(nbar) + vnbarup
              bvecdn(nbar) = bvecdn(nbar) + vnbardn
              do 40 j1 = 0, toff-1
                bvecup2(j1) = gmatup(j1,nbar)/(1.d0+bvecup(nbar))
40              bvecdn2(j1) = gmatdn(j1,nbar)/(1.d0+bvecdn(nbar))
              do 50 j1 = 0, toff-1
               do 50 j2 = 0, toff-1
                gmatup(j1,j2)=gmatup(j1,j2)-bvecup2(j1)*bvecup(j2)
50              gmatdn(j1,j2)=gmatdn(j1,j2)-bvecdn2(j1)*bvecdn(j2)
            else
              do 35 j1 = 0, toff-1
35              bvecup(j1) = -vnbarup*gmatup(nbar,j1)
              bvecup(nbar) = bvecup(nbar) + vnbarup
              do 45 j1 = 0, toff-1
45              bvecup2(j1) = gmatup(j1,nbar)/(1.d0+bvecup(nbar))
              do 55 j1 = 0, toff-1
               do 55 j2 = 0, toff-1
55              gmatup(j1,j2)=gmatup(j1,j2)-bvecup2(j1)*bvecup(j2)
            endif  

c                  ti1=mod(ti+1,l)
c                  call debug(vup,ti1,dbg,detdbg)
c                  aaa1=detdbg
c                  call debug(vup,ti1,dbg,detdbg)
c                  aaa2=detdbg

            hub(nbari) = -hub(nbari)
            vup(nbari) = vup(nbari) * (vnbarup+1.d0)
            vdn(nbari) = vdn(nbari) * (vnbardn+1.d0)
      
c                  call debug(vup,ti1,dbg,detdbg)
c                  maxdif=0.d0
c                  do 9000 ixx=0,toff-1
c                  do 9000 iyy=0,toff-1   
c                  dif=dabs(dbg(ixx,iyy)-gmatup(ixx,iyy))
c9000              if (dif.gt.maxdif) maxdif=dif
c                  write (6,*) maxdif,detdbg/aaa1/ratup
c                  call debug(vdn,ti1,dbg,detdbg)
c                  maxdif=0.d0
c                  do 9001 ixx=0,toff-1
c                  do 9001 iyy=0,toff-1   
c                  dif=dabs(dbg(ixx,iyy)-gmatdn(ixx,iyy))
c9001              if (dif.gt.maxdif) maxdif=dif
c                  write (6,*) maxdif,detdbg/aaa2/ratdn

          else
            reject = reject + 1
          endif
20      continue
 
        return
        end
 
c********* phaseset() - set and load phase vectors **********
        subroutine phaseset()
        implicit none
        integer n,l,volume,toff
        parameter(n=4,l=48)
        parameter(volume=n*n*l,toff=n*n)
        real*8 hub(0:volume-1)
        real*8 vup(0:volume-1),vdn(0:volume-1)
        real*8 mphase(0:toff-1)
        common/vectors/hub,vup,vdn, mphase
        integer ix,iy,addr
        real*8 ph
 
c        Calculate mphase - the phase vector needed in the measurement process
c        mphase = +1 on odd spatial sites
c        mphase = -1 on even spatial sites
 
        addr=0
        ph=1.d0
        do 30 iy=0,n-1
        do 25 ix=0,n-1
        mphase(addr)=ph
        addr=addr+1
25        ph=-ph
30        ph=-ph
 
        return
        end
 
c**************************indexset()************************
c        this subroutine sets and loads the index vectors
        subroutine indexset()
        implicit none
        integer n,l,volume,toff
        parameter(n=4,l=48)
        parameter(volume=n*n*l,toff=n*n)
c        index vectors for indirect addressing and phase vector
        integer xplus(0:toff-1),xminus(0:toff-1)
        integer yplus(0:toff-1),yminus(0:toff-1)
 
        common/ivectors/xplus,xminus,yplus,yminus
 
        integer i,neighbor
 
c        calculate the index vectors
        do 10 i=0,toff-1
            xplus(i) =neighbor(i,1,0,0)
            xminus(i)=neighbor(i,-1,0,0)
            yplus(i) =neighbor(i,0,1,0)
            yminus(i)=neighbor(i,0,-1,0)
10        continue
 
        return
        end
 
c*******************************************
        subroutine readin()
        implicit none
        integer n,l,volume,toff
        parameter(n=4,l=48)
        parameter(volume=n*n*l,toff=n*n)
 
        integer warms,sweeps,msr,nwrap
        integer iran
        common/integers/warms,sweeps,msr,nwrap, iran
        integer orthlen,doauto,torth
        real*8 eorth,difflim,errpam
        common/getgpam/eorth,difflim,errpam,orthlen,doauto,
     1                torth
 
        real*8 sgnup,sgndn,sgn,nup,ndn,nud,ke
        real*8 saf,saf2,sferro,sfer2,grfun(0:n/2,0:n/2),safsq,safsq2
        real*8 spinxx(0:n/2,0:n/2),spinzz(0:n/2,0:n/2)
        real*8 pairmat(-1:1,-1:1,-1:1,-1:1),errrat,asgnp
        integer nmeas0,redo,noredo,dopair,nmeasp,numpair
        real*8 asgnup,asgndn,asgn,anup,andn,anud,ake,asafsq,asafsq2
        real*8 asaf,asaf2,asferro,asfer2,agrfun(0:n/2,0:n/2)
        real*8 aspinxx(0:n/2,0:n/2),aspinzz(0:n/2,0:n/2)
        real*8 apairmat(-1:1,-1:1,-1:1,-1:1)
        real*8 den(0:n/2,0:n/2,0:1),aden(0:n/2,0:n/2,0:1)
        common/m0var/ sgnup,sgndn,sgn,nup,ndn,nud,ke,safsq,safsq2,
     1       saf,saf2,sferro,sfer2,grfun,
     1       spinxx,spinzz,aspinxx,aspinzz,asafsq,asafsq2,
     2       asgnup,asgndn,asgn,anup,andn,anud,ake,den,aden,
     3       asaf,asaf2,asferro,asfer2,agrfun,errrat,asgnp
     4       ,pairmat,apairmat
        common/m0vari/dopair,nmeasp,numpair,nmeas0,redo,noredo
 
         real*8 t,u,mu,dtau,expmu,gam,lambda,tdtau,dens
         common/couple/t,u,mu,dtau,expmu,gam,lambda,tdtau,dens
        real*8 ch,soc,ch4
        common/kinetic/ch,soc,ch4
 
        real*8 temp
 
        write(6,*) 'enter warms and sweeps'
        read(5,*) warms,sweeps
        write(76,*) 'warms=',warms,'  sweeps=',sweeps
 
        write(6,*)'enter t,u,dens '
        read(5,*) t,u,dens
        write(76,*) 't=',t,'  u=',u,' dens=',dens
 
        write(6,*) 'enter dtau, nwrap, difflim, errrat'
        read(5,*)dtau,nwrap,difflim,errrat
        write(76,*)'dtau=',dtau,' nwrap = ',nwrap,
     1         ' difflim= ', difflim,' errrat= ',errrat
 
        write(6,*) 'enter doauto, orthlen, eorth, dopair, numpair'
        read(5,*)doauto,orthlen,eorth,dopair, numpair
        write(76,*)' doauto= ',doauto,' orthlen= ',orthlen,
     1   ' eorth= ',eorth,' dopair= ',dopair,' numpair= ',numpair
 
        write(6,*) 'enter torth'
        read(5,*)torth
        write(76,*)' torth= ',torth
 
        gam = 0.5
        errpam = 1.0e13
        write(76,*)'errpam is ',errpam
 
        tdtau=t*dtau
        temp = dexp(dtau*u*0.5)
        lambda = dlog(temp+dsqrt(temp**2-1.d0))
        write(76,*)'lambda is ',lambda
        ch = dcosh(tdtau)
        soc = dsinh(tdtau)/ch
        ch4 = ch**4
 
        write(6,*)'enter random number seed'
        read(5,*)iran
        write(76,*) 'iran=',iran
 
        return
        end
 
c*****************************ranlat()*****************************
 
c        This subroutine makes a lattice with random values for the
c        Hubbard-Stratanovich variables.
 
        subroutine ranlat()
        implicit none
        integer n,l,volume,toff
        parameter(n=4,l=48)
        parameter(volume=n*n*l,toff=n*n)
        integer warms,sweeps,msr,nwrap
        integer iran
        common/integers/warms,sweeps,msr,nwrap, iran
 
         real*8 t,u,mu,dtau,expmu,gam,lambda,tdtau,dens
         common/couple/t,u,mu,dtau,expmu,gam,lambda,tdtau,dens
        real*8 hub(0:volume-1)
        real*8 vup(0:volume-1),vdn(0:volume-1)
        real*8 mphase(0:toff-1),ran2
        common/vectors/hub,vup,vdn, mphase
 
        integer i
        real*8 ranf
 
        do 10 i=0,volume-1
          hub(i)=1.d0
          ranf = ran2(iran)
          if(ranf .gt. 0.5) hub(i) = -1.d0
10      continue
 
        call setvup()
 
        return
        end
 
c****************************siteindx**********************************
 
c        this function finds the index of a site given its three coordinates
 
        integer function siteindx(x,y,ti)
        implicit none
        integer n,l,volume,toff
        parameter(n=4,l=48)
        parameter(volume=n*n*l,toff=n*n)
        integer x,y,ti
 
        siteindx = x + n*(y + n*(ti))
 
        return
        end
 
c****************************neighbor()*****************************
 
c        this function finds the index of a neighboring site
 
        integer function neighbor(site,delx,dely,delt)
        implicit none
        integer n,l,volume,toff
        parameter(n=4,l=48)
        parameter(volume=n*n*l,toff=n*n)
        integer site,delx,dely,delt
        integer x,y,ti,siteindx
 
c        find the coordinates of site
        x = mod(site,n)
        y = mod(site/n,n)
        ti = mod(site/(n*n),l)
 
c        find the coordinates of the neighbor
        x = mod(x+delx+n,n)
        y = mod(y+dely+n,n)
        ti = mod(ti+delt+l,l)
 
c        find the index of the neighboring site
        neighbor = siteindx(x,y,ti)
 
        return
        end
 
c************************GETERR*****************************
        subroutine geterr(bin,ave,err)
        real*8 bin(10),err,sum,sumsq,ave
        integer i
 
        sum = 0.d0
        do 10 i = 1, 10
10        sum = sum + bin(i)
        ave = sum * 0.1
        sumsq = 0.d0
        do 20 i = 1, 10
20          sumsq = sumsq + (bin(i)-ave)**2
        sumsq = 0.1 * sumsq
        err = dsqrt(sumsq)/3.d0
        if(err*1.0e12 .lt. dabs(ave)) err = 0.d0
        return
        end
 
c************************GETERRC*****************************
        subroutine geterrc(bin,ave,err)
        complex*16 bin(10),err,ave
        real*8 rsumsq,isumsq
        integer i
 
        ave = 0.d0
        do 10 i = 1, 10
10        ave = ave + bin(i)
        ave = ave * 0.1
        rsumsq = 0.d0
        do 20 i = 1, 10
          rsumsq = rsumsq + (dreal(bin(i)-ave))**2
20        isumsq = isumsq + (dimag(bin(i)-ave))**2
        rsumsq = 0.1 * rsumsq
        isumsq = 0.1 * isumsq
        err = cmplx(dsqrt(rsumsq)/3.d0,dsqrt(isumsq)/3.d0)
        return
        end
 
c************************GETERR2*****************************
        subroutine geterr2(sum,sumsq,ave,err)
        real*8 err,sum,sumsq,ave,ssq
 
        ave = sum * 0.1
        ssq = 0.1*sumsq - ave**2
        err = dsqrt(dabs(ssq))/3.d0
        if(err*1.0e12 .lt. dabs(ave)) err = 0.d0
        return
        end
 
ccccccccccccccccccccc Subroutine addaf() cccccccc
c        This subroutine adds in the small antiferromagnetic field to
c        vup and vdn.
 
        subroutine addaf()
        implicit none
        integer n,l,volume,toff
        parameter(n=4,l=48)
        parameter(volume=n*n*l,toff=n*n)
         real*8 t,u,mu,dtau,expmu,gam,lambda,tdtau,dens
         common/couple/t,u,mu,dtau,expmu,gam,lambda,tdtau,dens
        real*8 hub(0:volume-1)
        real*8 vup(0:volume-1),vdn(0:volume-1)
        real*8 mphase(0:toff-1)
        common/vectors/hub,vup,vdn, mphase
 
        integer x,y,ti,siteindx,s,i
        real*8 afeps,epsfacup,epsfacdn
        parameter (afeps=0.d0)
c
        if(afeps .ne. 0.d0) then
          do 10 x = 0, n-1
            do 10 y = 0, n-1
c            s = (-1)**(x+y)
              i = siteindx(x,y,0)
              s = mphase(i)
              epsfacup = dexp(-dtau * afeps * s)
              epsfacdn = 1.d0 / epsfacup
              do 10 ti = 0, l-1
                i = siteindx(x,y,ti)
                vup(i)=vup(i) * epsfacup
                vdn(i)=vdn(i) * epsfacdn
10        continue
        endif
        return
        end
c*****************************setvup()*****************************
 
c        This subroutine sets vup and vdn given hubs.
 
        subroutine setvup()
        implicit none
        integer n,l,volume,toff
        parameter(n=4,l=48)
        parameter(volume=n*n*l,toff=n*n)
        integer warms,sweeps,msr,nwrap
        integer iran
        common/integers/warms,sweeps,msr,nwrap, iran
 
 
         real*8 t,u,mu,dtau,expmu,gam,lambda,tdtau,dens
         common/couple/t,u,mu,dtau,expmu,gam,lambda,tdtau,dens
 
c        vectors:
 
        real*8 hub(0:volume-1)
        real*8 vup(0:volume-1),vdn(0:volume-1)
        real*8 mphase(0:toff-1)
 
c        common block for vectors
        common/vectors/hub,vup,vdn, mphase
 
        integer i
 
        do 20 i=0,volume-1
          vup(i)=dexp(lambda*hub(i))
          vdn(i)=dexp(-lambda*hub(i))
20      continue
        call addaf()
 
        return
        end
c*******************************multt(m,lr)********************
 
c        This subroutine multiplies m by the matrix T on the
c        left or right depending on whether lr = -1 or 1.
c        The factor expmu is also included.
c        This is the real*8 version.
 
        subroutine multt(mat,lr)
        implicit none
        integer n,l,volume,toff
        parameter(n=4,l=48)
        parameter(volume=n*n*l,toff=n*n)
 
         real*8 t,u,mu,dtau,expmu,gam,lambda,tdtau,dens
         common/couple/t,u,mu,dtau,expmu,gam,lambda,tdtau,dens
        real*8 ch,soc,ch4
        common/kinetic/ch,soc,ch4
        real*8 rux(0:n*n,0:n*n)
        common/rtmats/rux
 
        real*8 mat(0:n*n,0:n*n)
        integer i,j,k,lr
 
        if(lr .eq. -1) then
c Do odd x first
          do 40 i = 1, n*n-3,2
            do 40 j = 0, n*n-1
              rux(i,j) = mat(i,j) + soc*mat(i+1,j)
40            rux(i+1,j) = mat(i+1,j) + soc*mat(i,j)
          do 50 i = n-1, n*n-1,n
            do 50 j = 0, n*n-1
              rux(i,j) = mat(i,j) + soc*mat(i+1-n,j)
50            rux(i+1-n,j) = mat(i+1-n,j) + soc*mat(i,j)
c Do odd y next
          do 60 k = 0, n-1
           do 60 i = n, n*n-3*n,2*n
            do 60 j = 0, n*n-1
              mat(k+i,j) = rux(k+i,j) + soc*rux(k+i+n,j)
60            mat(k+i+n,j) = rux(k+i+n,j) + soc*rux(k+i,j)
          do 70 k = 0, n-1
            do 70 j = 0, n*n-1
              mat(k+n*n-n,j) = rux(k+n*n-n,j) + soc*rux(k,j)
70            mat(k,j) = rux(k,j) + soc*rux(k+n*n-n,j)
c Do even y next
          do 80 k = 0, n-1
           do 80 i = n, n*n-n,2*n
            do 80 j = 0, n*n-1
              rux(k+i,j) = mat(k+i,j) + soc*mat(k+i-n,j)
80            rux(k+i-n,j) = mat(k+i-n,j) + soc*mat(k+i,j)
c Do even x next
          do 90 i = 1, n*n-1,2
            do 90 j = 0, n*n-1
              mat(i,j)=(rux(i,j)+soc*rux(i-1,j)) * expmu*ch4
90            mat(i-1,j)=(rux(i-1,j)+soc*rux(i,j)) * expmu*ch4
        else
c Do even y first, using transposed mat.
          do 180 k = 0, n-1
           do 180 i = n, n*n-n,2*n
            do 180 j = 0, n*n-1
              rux(k+i,j) = mat(j,k+i) + soc*mat(j,k+i-n)
180           rux(k+i-n,j) = mat(j,k+i-n) + soc*mat(j,k+i)
c Do even x next
          do 190 i = 1, n*n-1,2
            do 190 j = 0, n*n-1
              mat(i,j)=(rux(i,j)+soc*rux(i-1,j))*expmu*ch4
190           mat(i-1,j)=(rux(i-1,j)+soc*rux(i,j))*expmu*ch4
c Do odd x next
          do 140 i = 1, n*n-3,2
            do 140 j = 0, n*n-1
              rux(i,j) = mat(i,j) + soc*mat(i+1,j)
140           rux(i+1,j) = mat(i+1,j) + soc*mat(i,j)
          do 150 i = n-1, n*n-1,n
            do 150 j = 0, n*n-1
              rux(i,j) = mat(i,j) + soc*mat(i+1-n,j)
150           rux(i+1-n,j) = mat(i+1-n,j) + soc*mat(i,j)
c Do odd y next
          do 160 k = 0, n-1
           do 160 i = n, n*n-3*n,2*n
            do 160 j = 0, n*n-1
c The following 6 lines includes the transposing of mat(k,j)
              mat(j,k+i) = rux(k+i,j) + soc*rux(k+i+n,j)
160           mat(j,k+i+n) = rux(k+i+n,j) + soc*rux(k+i,j)
          do 170 k = 0, n-1
            do 170 j = 0, n*n-1
              mat(j,k+n*n-n) = rux(k+n*n-n,j) + soc*rux(k,j)
170           mat(j,k) = rux(k,j) + soc*rux(k+n*n-n,j)
        endif
 
        return
        end
c*******************************multti(m,lr)********************
 
c        This subroutine multiplies m by the matrix T^-1 on the
c        left or right depending on whether lr = -1 or 1.
c        The factor expmu^-1 is also included.
c        This is the real*8 version.
 
        subroutine multti(mat,lr)
        implicit none
        integer n,l,volume,toff
        parameter(n=4,l=48)
        parameter(volume=n*n*l,toff=n*n)
 
         real*8 t,u,mu,dtau,expmu,gam,lambda,tdtau,dens
         common/couple/t,u,mu,dtau,expmu,gam,lambda,tdtau,dens
        real*8 ch,soc,ch4
        common/kinetic/ch,soc,ch4
 
        real*8 rux(0:n*n,0:n*n)
        common/rtmats/rux
 
        real*8 mat(0:n*n,0:n*n)
        integer i,j,k,lr
 
        if(lr .eq. 1) then
c Do odd x first, using transposed mat
          do 40 i = 1, n*n-3,2
            do 40 j = 0, n*n-1
              rux(i,j) = mat(j,i) - soc*mat(j,i+1)
40            rux(i+1,j) = mat(j,i+1) - soc*mat(j,i)
          do 50 i = n-1, n*n-1,n
            do 50 j = 0, n*n-1
              rux(i,j) = mat(j,i) - soc*mat(j,i+1-n)
50            rux(i+1-n,j) = mat(j,i+1-n) - soc*mat(j,i)
c Do odd y next
          do 60 k = 0, n-1
           do 60 i = n, n*n-3*n,2*n
            do 60 j = 0, n*n-1
              mat(k+i,j) = rux(k+i,j) - soc*rux(k+i+n,j)
60            mat(k+i+n,j) = rux(k+i+n,j) - soc*rux(k+i,j)
          do 70 k = 0, n-1
            do 70 j = 0, n*n-1
              mat(k+n*n-n,j) = rux(k+n*n-n,j) - soc*rux(k,j)
70            mat(k,j) = rux(k,j) - soc*rux(k+n*n-n,j)
c Do even y next
          do 80 k = 0, n-1
           do 80 i = n, n*n-n,2*n
            do 80 j = 0, n*n-1
              rux(k+i,j) = mat(k+i,j) - soc*mat(k+i-n,j)
80            rux(k+i-n,j) = mat(k+i-n,j) - soc*mat(k+i,j)
c Do even x next, and transpose back mat.
          do 90 i = 1, n*n-1,2
            do 90 j = 0, n*n-1
              mat(j,i)=(rux(i,j)-soc*rux(i-1,j)) /expmu *ch4
90            mat(j,i-1)=(rux(i-1,j)-soc*rux(i,j)) /expmu *ch4
        else
c for multiplying on left, do even y first.
          do 180 k = 0, n-1
           do 180 i = n, n*n-n,2*n
            do 180 j = 0, n*n-1
              rux(k+i,j) = mat(k+i,j) - soc*mat(k+i-n,j)
180           rux(k+i-n,j) = mat(k+i-n,j) - soc*mat(k+i,j)
c Do even x next
          do 190 i = 1, n*n-1,2
            do 190 j = 0, n*n-1
              mat(i,j)=(rux(i,j)-soc*rux(i-1,j))/expmu* ch4
190           mat(i-1,j)=(rux(i-1,j)-soc*rux(i,j))/expmu* ch4
c Do odd x next
          do 140 i = 1, n*n-3,2
            do 140 j = 0, n*n-1
              rux(i,j) = mat(i,j) - soc*mat(i+1,j)
140           rux(i+1,j) = mat(i+1,j) - soc*mat(i,j)
          do 150 i = n-1, n*n-1,n
            do 150 j = 0, n*n-1
              rux(i,j) = mat(i,j) - soc*mat(i+1-n,j)
150           rux(i+1-n,j) = mat(i+1-n,j) - soc*mat(i,j)
c Do odd y next
          do 160 k = 0, n-1
           do 160 i = n, n*n-3*n,2*n
            do 160 j = 0, n*n-1
              mat(k+i,j) = rux(k+i,j) - soc*rux(k+i+n,j)
160           mat(k+i+n,j) = rux(k+i+n,j) - soc*rux(k+i,j)
          do 170 k = 0, n-1
            do 170 j = 0, n*n-1
              mat(k+n*n-n,j) = rux(k+n*n-n,j) - soc*rux(k,j)
170           mat(k,j) = rux(k,j) - soc*rux(k+n*n-n,j)
        endif
 
        return
        end
c*******************************multb(m,vvv,ti,lr)********************
 
c        This subroutine multiplies m by the matrix B(ti) on the
c        left or right depending on whether lr = -1 or 1.
c        This is the real*8 version.
 
        subroutine multb(mat,vvv,ti,lr)
        implicit none
        integer n,l,volume,toff
        parameter(n=4,l=48)
        parameter(volume=n*n*l,toff=n*n)
        real*8 vvv(0:volume-1)
        real*8 mat(0:n*n,0:n*n)
        integer i,j,ti,lr
 
        if(lr .eq. -1) then
          call multt(mat,lr)
          do 10 i = 0, n*n-1
            do 10 j = 0, n*n-1
10            mat(i,j) = mat(i,j) * vvv(i+toff*ti)
        else
          do 20 j = 0, n*n-1
            do 20 i = 0, n*n-1
20            mat(i,j) = mat(i,j) * vvv(j+toff*ti)
          call multt(mat,lr)
        endif
 
        return
        end
c*******************************multbi(m,vvv,ti,lr)********************
 
c        This subroutine multiplies m by the matrix B(ti)^-1 on the
c        left or right depending on whether lr = -1 or 1.
c        This is the real*8 version.
 
        subroutine multbi(mat,vvv,ti,lr)
        implicit none
        integer n,l,volume,toff
        parameter(n=4,l=48)
        parameter(volume=n*n*l,toff=n*n)
        real*8 vvv(0:volume-1)
        real*8 mat(0:n*n,0:n*n)
        integer i,j,ti,lr
 
        if(lr .eq. 1) then
          call multti(mat,lr)
          do 10 j = 0, n*n-1
            do 10 i = 0, n*n-1
10            mat(i,j) = mat(i,j) / vvv(j+toff*ti)
        else
          do 20 i = 0, n*n-1
            do 20 j = 0, n*n-1
20            mat(i,j) = mat(i,j) / vvv(i+toff*ti)
          call multti(mat,lr)
        endif
 
        return
        end
c       ******************unit(mat)********************
c        This subroutine puts the unit matrix in aux, which is real*8.
 
        subroutine unit(mat)
        implicit none
        integer n,l,volume,toff
        parameter(n=4,l=48)
        parameter(volume=n*n*l,toff=n*n)
        real*8 mat(0:n*n,0:n*n)
        integer i,j
        do 10 i = 0, n*n-1
          do 20 j = 0, n*n-1
20          mat(i,j) = 0.d0
10        mat(i,i) = 1.d0
        return
        end
c       ******************zeromat(mat)********************
c        This subroutine zeroes mat.
 
        subroutine zeromat(mat)
        implicit none
        integer n,l,volume,toff
        parameter(n=4,l=48)
        parameter(volume=n*n*l,toff=n*n)
        real*8 mat(0:n*n,0:n*n)
        integer i,j
        do 20 i = 0, n*n-1
          do 20 j = 0, n*n-1
20          mat(i,j) = 0.d0
        return
        end
c **************************TRANSP(MAT,AUX)**********************
c        This subroutine multiplies mat by the even-x matrix
 
        subroutine transp(mat,aux)
        implicit none
        integer n,l
        parameter(n=4,l=48)
 
        integer i,j
        real*8 mat(0:n*n,0:n*n),aux(0:n*n,0:n*n)
 
        do 20 i = 0, n*n-1, 4
          do 20 j = 0, n*n-1
            aux(i,j)   = mat(j,i)
            aux(i+1,j) = mat(j,i+1)
            aux(i+2,j) = mat(j,i+2)
            aux(i+3,j) = mat(j,i+3)
20      continue
        return
        end
 
c**********************************
       subroutine matmult(mat1,mat2,mat3)
       implicit none
       integer n,l
       parameter(n=4,l=48)
       integer i,j,k
       real*8 mat1(0:n*n,0:n*n),mat2(0:n*n,0:n*n)
       real*8 mat3(0:n*n,0:n*n)
c
       do 10 i = 0, n*n-1
         do 10 j = 0, n*n-1
10         mat3(i,j) = 0.d0
       do 20 j = 0, n*n-1
         do 20 k = 0, n*n-1
           do 20 i = 0, n*n-1
20           mat3(i,j) = mat3(i,j) + mat1(i,k)*mat2(k,j)
       return
       end
c**********************************
       subroutine matdif(mat1,mat2,diff)
       implicit none
       integer n,l
       parameter(n=4,l=48)
       integer i,j
       real*8 mat1(0:n*n,0:n*n),mat2(0:n*n,0:n*n),diff
c
       diff = 0.d0
       do 10 i = 0, n*n-1
         do 10 j = 0, n*n-1
10         diff = diff + (mat1(i,j) - mat2(i,j))**2
       diff = dsqrt(diff) / (n*n)
       return
       end
c********************************
       subroutine matcop(mat1,mat2)
       implicit none
       integer n,l
       parameter(n=4,l=48)
       integer i,j
       real*8 mat1(0:n*n,0:n*n),mat2(0:n*n,0:n*n)
c
       do 10 i = 0, n*n-1
         do 10 j = 0, n*n-1
10         mat2(i,j) = mat1(i,j)
       return
       end
c**************************************
       subroutine donorm(vec,vnorm)
       implicit none
       integer n,l
       parameter(n=4,l=48)
       integer i
       real*8 vec(0:n*n),vnorm,temp
c
       temp = 0.d0
       do 10 i = 0, n*n-1
10       temp = temp + vec(i)**2
       vnorm = dsqrt(temp)
       if(vnorm .ne. 0.d0) then
         temp = 1.d0 / vnorm
         do 20 i = 0, n*n-1
20         vec(i) = vec(i) * temp
       endif
       return
       end
c*******************************meas0(lots of variables)********************
 
c        This subroutine does the measurements.
 
        subroutine meas0(gmatup,gmatdn)
        implicit none
        integer n,l,volume,toff
        parameter(n=4,l=48)
        parameter(volume=n*n*l,toff=n*n)
 
        integer warms,sweeps,msr,nwrap
        integer iran
        common/integers/warms,sweeps,msr,nwrap, iran
 
         real*8 t,u,mu,dtau,expmu,gam,lambda,tdtau,dens
         common/couple/t,u,mu,dtau,expmu,gam,lambda,tdtau,dens
 
c        index vectors for indirect addressing and phase vector
        integer xplus(0:toff-1),xminus(0:toff-1)
        integer yplus(0:toff-1),yminus(0:toff-1)
        common/ivectors/xplus,xminus,yplus,yminus
        real*8 hub(0:volume-1)
        real*8 vup(0:volume-1),vdn(0:volume-1)
        real*8 mphase(0:toff-1)
 
        real*8 sgnup,sgndn,sgn,nup,ndn,nud,ke
        real*8 saf,saf2,sferro,sfer2,grfun(0:n/2,0:n/2),safsq,safsq2
        real*8 spinxx(0:n/2,0:n/2),spinzz(0:n/2,0:n/2)
        real*8 pairmat(-1:1,-1:1,-1:1,-1:1),errrat,asgnp
        integer nmeas0,redo,noredo,dopair,nmeasp,numpair
        real*8 asgnup,asgndn,asgn,anup,andn,anud,ake,asafsq,asafsq2
        real*8 asaf,asaf2,asferro,asfer2,agrfun(0:n/2,0:n/2)
        real*8 aspinxx(0:n/2,0:n/2),aspinzz(0:n/2,0:n/2)
        real*8 apairmat(-1:1,-1:1,-1:1,-1:1)
        real*8 den(0:n/2,0:n/2,0:1),aden(0:n/2,0:n/2,0:1)
        common/m0var/ sgnup,sgndn,sgn,nup,ndn,nud,ke,safsq,safsq2,
     1       saf,saf2,sferro,sfer2,grfun,
     1       spinxx,spinzz,aspinxx,aspinzz,asafsq,asafsq2,
     2       asgnup,asgndn,asgn,anup,andn,anud,ake,den,aden,
     3       asaf,asaf2,asferro,asfer2,agrfun,errrat,asgnp
     4       ,pairmat,apairmat
        common/m0vari/dopair,nmeasp,numpair,nmeas0,redo,noredo
 
c        common block for vectors
        common/vectors/hub,vup,vdn, mphase
 
        real*8 gmatup(0:n*n,0:n*n),gmatdn(0:n*n,0:n*n)
        real*8 safch,gtemp,xxtemp,zztemp,xfac,yfac
 
        integer i,j,ix,jx,iy,jy,kx,ky
 
        saf = 0.d0
        safsq = 0.d0
        sferro = 0.d0
        sfer2 = 0.d0
        saf2 = 0.d0
        safsq2 = 0.d0
        safch = 0.d0
        nup = 0.d0
        ndn = 0.d0
        ke = 0.d0
        nud = 0.d0
        if(mu .eq. 0.d0) then
          do 10 i = 0,n*n-1
           do 20 j = 0,n*n-1
20           gmatdn(i,j) = -mphase(i)*mphase(j)*gmatup(j,i)
10        gmatdn(i,i) = gmatdn(i,i) + 1.d0
        endif
 
        do 30 kx = 0, n/2
          do 30 ky = 0, n/2
            grfun(kx,ky) = 0.d0
            spinxx(kx,ky) = 0.d0
            spinzz(kx,ky) = 0.d0
            den(kx,ky,0) = 0.d0
            den(kx,ky,1) = 0.d0
30      continue
 
        do 200 i = 0,n*n-1
          nup = nup + 1.d0 - gmatup(i,i)
          ndn = ndn + 1.d0 - gmatdn(i,i)
          ke=ke+gmatup(i,xplus(i))+gmatup(i,xminus(i))
     1              +gmatup(i,yplus(i))+gmatup(i,yminus(i))
     2              +gmatdn(i,xplus(i))+gmatdn(i,xminus(i))
     3              +gmatdn(i,yplus(i))+gmatdn(i,yminus(i))
          nud=nud+(1.d0-gmatup(i,i))*(1.d0-gmatdn(i,i))
          sferro = sferro + gmatup(i,i)+gmatdn(i,i)
          sfer2 = sfer2 + gmatup(i,i)+gmatdn(i,i)
          saf = saf + gmatup(i,i)+gmatdn(i,i)
200       saf2 = saf2 + gmatup(i,i)+gmatdn(i,i)
 
        spinxx(0,0) = saf
        spinzz(0,0) = saf
        do 210 ix = 0,n-1
         do 210 iy = 0,n-1
          i = ix+n*iy
          do 210 jx = 0,n-1
           do 210 jy = 0,n-1
            j = jx+n*jy
            kx = abs(ix-jx)
            kx = min(kx,n-kx)
            ky = abs(iy-jy)
            ky = min(ky,n-ky)
            gtemp = gmatup(i,j)+ gmatdn(i,j)
            grfun(kx,ky) = grfun(kx,ky) + gtemp
c The following line is not correct for kx=ky=0, but it will be fixed later.
            den(kx,ky,0) = den(kx,ky,0)+
     1             (1.d0-gmatup(i,i))*(1.d0-gmatup(j,j)) +
     2             (1.d0-gmatdn(i,i))*(1.d0-gmatdn(j,j)) -
     3           gmatup(i,j)*gmatup(j,i)-gmatdn(i,j)*gmatdn(j,i)
            den(kx,ky,1) = den(kx,ky,1)+
     1             (1.d0-gmatup(i,i))*(1.d0-gmatdn(j,j))
            xxtemp = -2.d0*gmatup(i,j)*gmatdn(j,i)
            zztemp = (gmatup(i,i)*gmatup(j,j)
     1        +gmatdn(i,i)*gmatdn(j,j)-2.d0*gmatup(i,i)*gmatdn(j,j)
     2        -gmatup(j,i)*gmatup(i,j) - gmatdn(j,i)*gmatdn(i,j))
            spinxx(kx,ky) = spinxx(kx,ky) + xxtemp
            spinzz(kx,ky) = spinzz(kx,ky) + zztemp
            sferro=sferro + xxtemp
            sfer2=sfer2 + zztemp
            saf=saf + mphase(i)*mphase(j)*xxtemp
210         saf2=saf2 + mphase(i)*mphase(j)*zztemp
 
        saf = saf / toff
        sferro = sferro / toff
        sfer2 = sfer2 / toff
        saf2 = saf2 / toff
        ke = ke*t / toff
        nud = nud / toff
        nup = nup / toff
        ndn = ndn / toff
        safsq = saf**2
        safsq2 = saf2**2
 
        do 41 kx = 0, n/2
          if(kx .eq. 0 .or. kx .eq. n/2)then
            xfac = 1.d0
          else
            xfac = 2.d0
          endif
          do 41 ky = 0, n/2
            if(ky .eq. 0 .or. ky .eq. n/2)then
              yfac = 1.d0
            else
              yfac = 2.d0
            endif
            spinxx(kx,ky) = spinxx(kx,ky) / (n*n*xfac*yfac)
            spinzz(kx,ky) = spinzz(kx,ky) / (n*n*xfac*yfac)
            grfun(kx,ky) = grfun(kx,ky) / (2*n*n*xfac*yfac)
            den(kx,ky,0) = den(kx,ky,0) / (2*n*n*xfac*yfac)
            den(kx,ky,1) = den(kx,ky,1) / (n*n*xfac*yfac)
41      continue
 
c See note above:
        den(0,0,0) = (nup + ndn) * 0.5
 
        do 31 kx = 0, n/2
          do 31 ky = 0, n/2
            grfun(kx,ky) = (grfun(kx,ky)+grfun(ky,kx))*0.5
            grfun(ky,kx) = grfun(kx,ky)
            den(kx,ky,0) = (den(kx,ky,0)+den(ky,kx,0))*0.5
            den(ky,kx,0) = den(kx,ky,0)
            den(kx,ky,1) = (den(kx,ky,1)+den(ky,kx,1))*0.5
            den(ky,kx,1) = den(kx,ky,1)
            spinxx(kx,ky) = (spinxx(kx,ky) + spinxx(ky,kx))*0.5
            spinxx(ky,kx) = spinxx(kx,ky)
            spinzz(kx,ky) = (spinzz(kx,ky) + spinzz(ky,kx))*0.5
            spinzz(ky,kx) = spinzz(kx,ky)
31      continue
 
        nmeas0 = nmeas0 + 1
        sgn = sgnup*sgndn
        anup = anup + nup*sgn
        andn = andn + ndn*sgn
        asaf = asaf + saf*sgn
        asafsq = asafsq + safsq*sgn
        asferro = asferro + sferro*sgn
        asfer2 = asfer2 + sfer2*sgn
        asaf2 = asaf2 + saf2*sgn
        asafsq2 = asafsq2 + safsq2*sgn
        anud = anud + nud*sgn
        ake = ake + ke*sgn
        asgn = asgn + sgn
        asgnup = asgnup + sgnup
        asgndn = asgndn + sgndn
        do 51 kx = 0, n/2
          do 51 ky = 0, n/2
            agrfun(kx,ky) = agrfun(kx,ky) + grfun(kx,ky)*sgn
            aden(kx,ky,0) = aden(kx,ky,0) + den(kx,ky,0)*sgn
            aden(kx,ky,1) = aden(kx,ky,1) + den(kx,ky,1)*sgn
            aspinxx(kx,ky) = aspinxx(kx,ky) + spinxx(kx,ky)*sgn
51          aspinzz(kx,ky) = aspinzz(kx,ky) + spinzz(kx,ky)*sgn
 
        return
        end
c*******************************measpair(gmatup,gmatdn)********************
 
c        This subroutine does the measurements.
 
        subroutine measpair(gmatup,gmatdn)
        implicit none
        integer n,l,volume,toff
        parameter(n=4,l=48)
        parameter(volume=n*n*l,toff=n*n)
 
        integer warms,sweeps,msr,nwrap
        integer iran
        common/integers/warms,sweeps,msr,nwrap, iran
 
         real*8 t,u,mu,dtau,expmu,gam,lambda,tdtau,dens
         common/couple/t,u,mu,dtau,expmu,gam,lambda,tdtau,dens
 
c        index vectors for indirect addressing and phase vector
        integer xplus(0:toff-1),xminus(0:toff-1)
        integer yplus(0:toff-1),yminus(0:toff-1)
        common/ivectors/xplus,xminus,yplus,yminus
        real*8 hub(0:volume-1)
        real*8 vup(0:volume-1),vdn(0:volume-1)
        real*8 mphase(0:toff-1)
 
        real*8 sgnup,sgndn,sgn,nup,ndn,nud,ke
        real*8 saf,saf2,sferro,sfer2,grfun(0:n/2,0:n/2),safsq,safsq2
        real*8 spinxx(0:n/2,0:n/2),spinzz(0:n/2,0:n/2)
        real*8 pairmat(-1:1,-1:1,-1:1,-1:1),errrat,asgnp
        integer nmeas0,redo,noredo,dopair,nmeasp,numpair
        real*8 asgnup,asgndn,asgn,anup,andn,anud,ake,asafsq,asafsq2
        real*8 asaf,asaf2,asferro,asfer2,agrfun(0:n/2,0:n/2)
        real*8 aspinxx(0:n/2,0:n/2),aspinzz(0:n/2,0:n/2)
        real*8 apairmat(-1:1,-1:1,-1:1,-1:1)
        real*8 den(0:n/2,0:n/2,0:1),aden(0:n/2,0:n/2,0:1)
        common/m0var/ sgnup,sgndn,sgn,nup,ndn,nud,ke,safsq,safsq2,
     1       saf,saf2,sferro,sfer2,grfun,
     1       spinxx,spinzz,aspinxx,aspinzz,asafsq,asafsq2,
     2       asgnup,asgndn,asgn,anup,andn,anud,ake,den,aden,
     3       asaf,asaf2,asferro,asfer2,agrfun,errrat,asgnp
     4       ,pairmat,apairmat
        common/m0vari/dopair,nmeasp,numpair,nmeas0,redo,noredo
 
c        common block for vectors
        common/vectors/hub,vup,vdn, mphase
 
        real*8 gmatup(0:n*n,0:n*n),gmatdn(0:n*n,0:n*n)
 
        integer i,j,mx,my,mpx,mpy
        integer inn,in0,in1,i0n,i01,i1n,i10,i11,jp(-1:1),jpp
 
        if(mu .eq. 0.d0) then
          do 10 i = 0,n*n-1
           do 20 j = 0,n*n-1
20           gmatdn(i,j) = -mphase(i)*mphase(j)*gmatup(j,i)
10        gmatdn(i,i) = gmatdn(i,i) + 1.d0
        endif
 
        do 294 mx = -1,1
         do 294 my = -1,1
          do 294 mpx = -1,1
           do 294 mpy = -1,1
294          pairmat(mx,my,mpx,mpy) = 0.d0
        do 522 i = 0,n*n-1
          inn = xminus(yminus(i))
          in0 = xminus(i)
          in1 = xminus(yplus(i))
          i0n = yminus(i)
          i01 = yplus(i)
          i1n = xplus(yminus(i))
          i10 = xplus(i)
          i11 = xplus(yplus(i))
          do 522 j = 0,n*n-1
            jpp = xplus(xplus(j))
            do 500 mx = -1,1
             jpp = xminus(jpp)
             jp(-1) = yplus(jpp)
             jp(0) = jpp
             jp(1) = yminus(jpp)
             do 500 my = -1,1
               pairmat(mx,my,-1,-1) = pairmat(mx,my,-1,-1) +
     1             gmatup(jp(my),i11)*gmatdn(j,i)
               pairmat(mx,my,-1,0) = pairmat(mx,my,-1,0) +
     1             gmatup(jp(my),i10)*gmatdn(j,i)
               pairmat(mx,my,-1,1) = pairmat(mx,my,-1,1) +
     1             gmatup(jp(my),i1n)*gmatdn(j,i)
               pairmat(mx,my,0,-1) = pairmat(mx,my,0,-1) +
     1             gmatup(jp(my),i01)*gmatdn(j,i)
               pairmat(mx,my,0,0) = pairmat(mx,my,0,0) +
     1             gmatup(jp(my),i)*gmatdn(j,i)
               pairmat(mx,my,0,1) = pairmat(mx,my,0,1) +
     1             gmatup(jp(my),i0n)*gmatdn(j,i)
               pairmat(mx,my,1,-1) = pairmat(mx,my,1,-1) +
     1             gmatup(jp(my),in1)*gmatdn(j,i)
               pairmat(mx,my,1,0) = pairmat(mx,my,1,0) +
     1             gmatup(jp(my),in0)*gmatdn(j,i)
               pairmat(mx,my,1,1) = pairmat(mx,my,1,1) +
     1             gmatup(jp(my),inn)*gmatdn(j,i)
500         continue
522     continue
 
        nmeasp = nmeasp + 1
        do 494 mx = -1,1
         do 494 my = -1,1
          do 494 mpx = -1,1
           do 494 mpy = -1,1
            pairmat(mx,my,mpx,mpy)=pairmat(mx,my,mpx,mpy)/toff
494         apairmat(mx,my,mpx,mpy)=apairmat(mx,my,mpx,mpy) +
     1                  sgnup*sgndn*pairmat(mx,my,mpx,mpy)
        asgnp = asgnp + sgnup*sgndn
 
        return
        end
c*******************************meastau********************
 
c        This subroutine does the measurements.
 
        subroutine meastau
        implicit none
        integer n,l,volume,toff
        parameter(n=4,l=48)
        parameter(volume=n*n*l,toff=n*n)
 
        integer warms,sweeps,msr,nwrap
        integer iran
        common/integers/warms,sweeps,msr,nwrap, iran
        integer orthlen,doauto,torth
        real*8 eorth,difflim,errpam
        common/getgpam/eorth,difflim,errpam,orthlen,doauto,
     1                torth
 
 
         real*8 t,u,mu,dtau,expmu,gam,lambda,tdtau,dens
         common/couple/t,u,mu,dtau,expmu,gam,lambda,tdtau,dens
 
c        index vectors for indirect addressing and phase vector
        integer xplus(0:toff-1),xminus(0:toff-1)
        integer yplus(0:toff-1),yminus(0:toff-1)
        common/ivectors/xplus,xminus,yplus,yminus
        real*8 hub(0:volume-1)
        real*8 vup(0:volume-1),vdn(0:volume-1)
        real*8 mphase(0:toff-1)
 
        real*8 sgnup,sgndn,sgn,nup,ndn,nud,ke
        real*8 saf,saf2,sferro,sfer2,grfun(0:n/2,0:n/2),safsq,safsq2
        real*8 spinxx(0:n/2,0:n/2),spinzz(0:n/2,0:n/2)
        real*8 pairmat(-1:1,-1:1,-1:1,-1:1),errrat,asgnp
        integer nmeas0,redo,noredo,dopair,nmeasp,numpair
        real*8 asgnup,asgndn,asgn,anup,andn,anud,ake,asafsq,asafsq2
        real*8 asaf,asaf2,asferro,asfer2,agrfun(0:n/2,0:n/2)
        real*8 aspinxx(0:n/2,0:n/2),aspinzz(0:n/2,0:n/2)
        real*8 apairmat(-1:1,-1:1,-1:1,-1:1)
        real*8 den(0:n/2,0:n/2,0:1),aden(0:n/2,0:n/2,0:1)
        common/m0var/ sgnup,sgndn,sgn,nup,ndn,nud,ke,safsq,safsq2,
     1       saf,saf2,sferro,sfer2,grfun,
     1       spinxx,spinzz,aspinxx,aspinzz,asafsq,asafsq2,
     2       asgnup,asgndn,asgn,anup,andn,anud,ake,den,aden,
     3       asaf,asaf2,asferro,asfer2,agrfun,errrat,asgnp
     4       ,pairmat,apairmat
        common/m0vari/dopair,nmeasp,numpair,nmeas0,redo,noredo
        real*8 gnl(0:n/2,0:n/2,0:l),agnl(0:n/2,0:n/2,0:l)
        real*8 chinl(0:n/2,0:n/2,0:l),achinl(0:n/2,0:n/2,0:l)
        real*8 pairsus(-1:1,-1:1,-1:1,-1:1,0:l)
        real*8 apairsus(-1:1,-1:1,-1:1,-1:1,0:l),asgnt
        integer nmeast
        common/mtauvar/gnl,agnl,chinl,achinl,pairsus,
     1       apairsus,asgnt,nmeast
 
c        common block for vectors
        common/vectors/hub,vup,vdn, mphase
 
        integer ti,i,j,ix,jx,iy,jy,kx,ky,sx,sy,mx,my,mpx,mpy
        integer inn,in0,in1,i0n,i01,i1n,i10,i11,jp(-1:1),jpp
        real*8 gupt0(0:n*n,0:n*n),gup0t(0:n*n,0:n*n),chitemp
        real*8 gdnt0(0:n*n,0:n*n),gdn0t(0:n*n,0:n*n),gtemp
        integer inacc
 
100     do 15 ti = 0, l
         do 15 mx = -1,1
          do 15 my = -1,1
           do 15 mpx = -1,1
            do 15 mpy = -1,1
15           pairsus(mx,my,mpx,mpy,ti) = 0.d0
        do 10 ti = 0, l-1
         call makegt(0,ti,gupt0,gup0t,gdnt0,gdn0t,inacc)
         if(inacc .eq. 1) then
           torth = max(1,torth-1)
           write(6,*)'reducing torth to ',torth
           write(76,*)'reducing torth to ',torth
           goto 100
         endif
         do 40 kx = 0,n/2
          do 40 ky = 0, n/2
           gtemp = 0.d0
           chitemp = 0.d0
           do 30 sx = -1, 1, 2
            do 30 sy = -1, 1, 2
             do 30 ix = 0,n-1
              do 30 iy = 0,n-1
               jx = mod(ix+sx*kx+n,n)
               jy = mod(iy+sy*ky+n,n)
               i = ix+n*iy
               j = jx+n*jy
               gtemp = gtemp + gupt0(i,j)+ gdnt0(i,j)
               chitemp = chitemp - gup0t(j,i)*gdnt0(i,j)
     1                           - gdn0t(j,i)*gupt0(i,j)
30         continue
           gnl(kx,ky,ti) = gtemp / (8*n*n)
           chinl(kx,ky,ti) = chitemp / (4*n*n)
           if(ti .eq. 0) then
             gnl(kx,ky,l) = -gnl(kx,ky,0)
             chinl(kx,ky,l) = chinl(kx,ky,0)
           endif
40       continue
         gnl(0,0,l) = 1.d0 - gnl(0,0,0)
         
        if(dopair .eq. 1) then
        do 522 i = 0,n*n-1
          inn = xminus(yminus(i))
          in0 = xminus(i)
          in1 = xminus(yplus(i))
          i0n = yminus(i)
          i01 = yplus(i)
          i1n = xplus(yminus(i))
          i10 = xplus(i)
          i11 = xplus(yplus(i))
          do 522 j = 0,n*n-1
            jpp = xplus(xplus(j))
            do 500 mx = -1,1
             jpp = xminus(jpp)
             jp(-1) = yplus(jpp)
             jp(0) = jpp
             jp(1) = yminus(jpp)
             do 510 my = -1,1
               pairsus(mx,my,-1,-1,ti) = pairsus(mx,my,-1,-1,ti) +
     1             gupt0(jp(my),i11)*gdnt0(j,i)
               pairsus(mx,my,-1,0,ti) = pairsus(mx,my,-1,0,ti) +
     1             gupt0(jp(my),i10)*gdnt0(j,i)
               pairsus(mx,my,-1,1,ti) = pairsus(mx,my,-1,1,ti) +
     1             gupt0(jp(my),i1n)*gdnt0(j,i)
               pairsus(mx,my,0,-1,ti) = pairsus(mx,my,0,-1,ti) +
     1             gupt0(jp(my),i01)*gdnt0(j,i)
               pairsus(mx,my,0,0,ti) = pairsus(mx,my,0,0,ti) +
     1             gupt0(jp(my),i)*gdnt0(j,i)
               pairsus(mx,my,0,1,ti) = pairsus(mx,my,0,1,ti) +
     1             gupt0(jp(my),i0n)*gdnt0(j,i)
               pairsus(mx,my,1,-1,ti) = pairsus(mx,my,1,-1,ti) +
     1             gupt0(jp(my),in1)*gdnt0(j,i)
               pairsus(mx,my,1,0,ti) = pairsus(mx,my,1,0,ti) +
     1             gupt0(jp(my),in0)*gdnt0(j,i)
               pairsus(mx,my,1,1,ti) = pairsus(mx,my,1,1,ti) +
     1             gupt0(jp(my),inn)*gdnt0(j,i)
510         continue
            if(ti .eq. 0)then
c g(beta) = 1-g(0); but for ti=0, g0t(0) = g(0)-1, minus signs cancel.
             do 515 my = -1,1
               pairsus(mx,my,-1,-1,l) = pairsus(mx,my,-1,-1,l) +
     1             gup0t(jp(my),i11)*gdn0t(j,i)
               pairsus(mx,my,-1,0,l) = pairsus(mx,my,-1,0,l) +
     1             gup0t(jp(my),i10)*gdn0t(j,i)
               pairsus(mx,my,-1,1,l) = pairsus(mx,my,-1,1,l) +
     1             gup0t(jp(my),i1n)*gdn0t(j,i)
               pairsus(mx,my,0,-1,l) = pairsus(mx,my,0,-1,l) +
     1             gup0t(jp(my),i01)*gdn0t(j,i)
               pairsus(mx,my,0,0,l) = pairsus(mx,my,0,0,l) +
     1             gup0t(jp(my),i)*gdn0t(j,i)
               pairsus(mx,my,0,1,l) = pairsus(mx,my,0,1,l) +
     1             gup0t(jp(my),i0n)*gdn0t(j,i)
               pairsus(mx,my,1,-1,l) = pairsus(mx,my,1,-1,l) +
     1             gup0t(jp(my),in1)*gdn0t(j,i)
               pairsus(mx,my,1,0,l) = pairsus(mx,my,1,0,l) +
     1             gup0t(jp(my),in0)*gdn0t(j,i)
               pairsus(mx,my,1,1,l) = pairsus(mx,my,1,1,l) +
     1             gup0t(jp(my),inn)*gdn0t(j,i)
515          continue
            endif
500         continue
522     continue
        endif
 
10     continue
 
        nmeast = nmeast + 1
        asgnt = asgnt + sgnup*sgndn
        do 594 mx = 0, n/2
         do 594 my = 0, n/2
          do 594 ti = 0,l
           agnl(mx,my,ti)=agnl(mx,my,ti)+gnl(mx,my,ti)*sgnup*sgndn
594        achinl(mx,my,ti)=achinl(mx,my,ti)+
     1                          chinl(mx,my,ti)*sgnup*sgndn
        do 494 mx = -1,1
        do 494 my = -1,1
         do 494 mpx = -1,1
         do 494 mpy = -1,1
          do 494 ti = 0, l
          pairsus(mx,my,mpx,mpy,ti)=pairsus(mx,my,mpx,mpy,ti)/toff
494       apairsus(mx,my,mpx,mpy,ti)=apairsus(mx,my,mpx,mpy,ti) +
     1                  sgnup*sgndn*pairsus(mx,my,mpx,mpy,ti)
 
        return
        end
c***************************************
        subroutine zeroas()
        implicit none
        integer n,l,volume,toff
        parameter(n=4,l=48)
        parameter(volume=n*n*l,toff=n*n)
 
        real*8 sgnup,sgndn,sgn,nup,ndn,nud,ke
        real*8 saf,saf2,sferro,sfer2,grfun(0:n/2,0:n/2),safsq,safsq2
        real*8 spinxx(0:n/2,0:n/2),spinzz(0:n/2,0:n/2)
        real*8 pairmat(-1:1,-1:1,-1:1,-1:1),errrat,asgnp
        integer nmeas0,redo,noredo,dopair,nmeasp,numpair
        real*8 asgnup,asgndn,asgn,anup,andn,anud,ake,asafsq,asafsq2
        real*8 asaf,asaf2,asferro,asfer2,agrfun(0:n/2,0:n/2)
        real*8 aspinxx(0:n/2,0:n/2),aspinzz(0:n/2,0:n/2)
        real*8 apairmat(-1:1,-1:1,-1:1,-1:1)
        real*8 den(0:n/2,0:n/2,0:1),aden(0:n/2,0:n/2,0:1)
        common/m0var/ sgnup,sgndn,sgn,nup,ndn,nud,ke,safsq,safsq2,
     1       saf,saf2,sferro,sfer2,grfun,
     1       spinxx,spinzz,aspinxx,aspinzz,asafsq,asafsq2,
     2       asgnup,asgndn,asgn,anup,andn,anud,ake,den,aden,
     3       asaf,asaf2,asferro,asfer2,agrfun,errrat,asgnp
     4       ,pairmat,apairmat
        common/m0vari/dopair,nmeasp,numpair,nmeas0,redo,noredo
        real*8 gnl(0:n/2,0:n/2,0:l),agnl(0:n/2,0:n/2,0:l)
        real*8 chinl(0:n/2,0:n/2,0:l),achinl(0:n/2,0:n/2,0:l)
        real*8 pairsus(-1:1,-1:1,-1:1,-1:1,0:l)
        real*8 apairsus(-1:1,-1:1,-1:1,-1:1,0:l),asgnt
        integer nmeast
        common/mtauvar/gnl,agnl,chinl,achinl,pairsus,
     1       apairsus,asgnt,nmeast
 
        integer i,j,i2,j2,ti
 
        nmeas0 = 0
        nmeasp = 0
        nmeast = 0
        call setto0(asgnup,asgndn,asgn,anup,andn)
        call setto0(anud,ake,asaf,asaf2,asferro)
        call setto0(asfer2,asafsq,asafsq2,asgnp,asgnt)
        do 10 i = 0, n/2
          do 10 j = 0, n/2
            agrfun(i,j) = 0.d0
            aspinxx(i,j) = 0.d0
            aspinzz(i,j) = 0.d0
            aden(i,j,0) = 0.d0
            aden(i,j,1) = 0.d0
            do 10 ti = 0, l
              agnl(i,j,ti) = 0.d0
              achinl(i,j,ti) = 0.d0
10      continue
        if(dopair .eq. 1)then
          do 20 i = -1,1
            do 20 j = -1,1
              do 20 i2 = -1,1
                do 20 j2 = -1,1
                  apairmat(i,j,i2,j2) = 0.d0
                  do 20 ti = 0, l
                    apairsus(i,j,i2,j2,ti) = 0.d0
20        continue
        endif
 
        return
        end
c**************************************************
        subroutine autoset()
        implicit none
        integer n,l,volume,toff
        parameter(n=4,l=48)
        parameter(volume=n*n*l,toff=n*n)
        integer warms,sweeps,msr,nwrap
        integer iran
        common/integers/warms,sweeps,msr,nwrap, iran
        integer orthlen,doauto,torth
        real*8 eorth,difflim,errpam
        common/getgpam/eorth,difflim,errpam,orthlen,doauto,
     1                torth
 
         real*8 t,u,mu,dtau,expmu,gam,lambda,tdtau,dens
         common/couple/t,u,mu,dtau,expmu,gam,lambda,tdtau,dens
        real*8 hub(0:volume-1)
        real*8 vup(0:volume-1),vdn(0:volume-1)
        real*8 mphase(0:toff-1)
 
c        common block for vectors
        common/vectors/hub,vup,vdn, mphase
 
        integer start
        real*8 gmatup(0:n*n,0:n*n),diffup,sgnup,detup
        real*8 gmatacc(0:n*n,0:n*n)
        
        write(76,*)'eorth is ',eorth
        start=orthlen
        orthlen = 1
        call ranlat()
        call getgp(vup,0,gmatacc,sgnup,detup)
5       continue
        do 10 orthlen = start, 2, -1
          call getgp(vup,0,gmatup,sgnup,detup)
          call matdif(gmatup,gmatacc,diffup)
          if(doauto .eq. 1)
     1      write(6,*)'diffup for orthlen= ',orthlen,' is ',diffup
          if(doauto .ne. 1 .or. diffup .le. eorth)goto 40
10      continue
        eorth = eorth * 100.d0
        write(76,*)'resetting eorth to ',eorth
        goto 5
40      continue
        write(6,*)'Using orthlen= ',orthlen,' diffup is ',diffup
        write(76,*)'Using orthlen= ',orthlen,' diffup is ',diffup
        return
        end
c
       subroutine setto0(a,b,c,d,e)
       real*8 a,b,c,d,e
       a=0.d0
       b=0.d0
       c=0.d0
       d=0.d0
       e=0.d0
       return
       end
c******************makegt(ti,dti,gtup,gsup,gtdn,gsdn)*************
c This subroutine returns the unequal time Green's functions
c      gtup=Gup(ti+dti,ti)   gsup=Gup(ti,ti+dti)
c      gtdn=Gdn(ti+dti,ti)   gsdn=Gdn(ti,ti+dti)
 
c In the measurement routines gtup,... should not be changed since unless
c dti=0, it is assumed that they contain Gup(ti+dti-1,ti),...
c Note that for dti=0 gtup and gtdn return the equal time Green's functions,
c while gsup and gsdn return gtup-1 and gtdn-1, which are not needed for
c measurements, but are needed for future calculations of gsup and gsdn.
c Note that dti should range from 0 to l-1.
 
       subroutine makegt(ti,dti,gtup,gsup,gtdn,gsdn,inacc)
       implicit none
        integer n,l,volume,toff
        parameter(n=4,l=48)
        parameter(volume=n*n*l,toff=n*n)
       integer orthlen,doauto,torth
       real*8 eorth,difflim,errpam
       common/getgpam/eorth,difflim,errpam,orthlen,doauto,
     1                torth
       integer ti,dti,i,j,kk,inacc
       real*8 hub(0:volume-1),vup(0:volume-1),vdn(0:volume-1)
       real*8 gtup(0:n*n,0:n*n),gsup(0:n*n,0:n*n)
        real*8 mphase(0:toff-1)
        common/vectors/hub,vup,vdn, mphase
       real*8 gtdn(0:n*n,0:n*n),gsdn(0:n*n,0:n*n)
       real*8 gt2(0:n*n,0:n*n),gs2(0:n*n,0:n*n),difft,diffs
 
c Calculate the Green's functions on the next time slice.
       inacc = 0
       kk=mod(dti,torth)
c Multiply gt on the left by B(ti+dti) and gs on the right by B^-1(ti+dti)
       i=mod(ti+dti-1,l)
       call multb(gtup,vup,i,-1)
       call multb(gtdn,vdn,i,-1)
       call multbi(gsup,vup,i,1)
       call multbi(gsdn,vdn,i,1)
       if (kk.eq.0) then
c If kk=0 it is time to calculate the Green's functions from scratch.
        call matcop(gtup,gt2)
        call matcop(gsup,gs2)
        call getgtau(vup,ti,dti,gtup,gsup)
c getgtau returns gsup and gsdn with the wrong sign. Correct this.
        do 20 j=0,n*n-1
         do 20 i=0,n*n-1
          gsup(i,j)=-gsup(i,j)
20      continue
        call matdif(gtup,gt2,difft)
        call matdif(gsup,gs2,diffs)
        if(difft .gt. 1.0e-4 .or. diffs .gt. 1.0e-4)then
          if(dti .ne. 0)then
            inacc = 1
            write(6,*)'difft,diffs are ',difft,diffs
            return
          endif
        endif
        call matcop(gtdn,gt2)
        call matcop(gsdn,gs2)
        call getgtau(vdn,ti,dti,gtdn,gsdn)
        do 25 j=0,n*n-1
         do 25 i=0,n*n-1
          gsdn(i,j)=-gsdn(i,j)
25      continue
        call matdif(gtdn,gt2,difft)
        call matdif(gsdn,gs2,diffs)
        if(difft .gt. 1.0e-4 .or. diffs .gt. 1.0e-4)then
          if(dti .ne. 0)then
            inacc = 1
            write(6,*)'difft,diffs are ',difft,diffs
            return
          endif
        endif
       endif
 
        return
        end
c***********************************
c This getudr includes pivoting.
       subroutine getudr(vvv,ti1,ti2,dodepth,orthlen,omat,bvec,rmat)
       implicit none
        integer n,l,volume,toff
        parameter(n=4,l=48)
        parameter(volume=n*n*l,toff=n*n)
       integer i,mx,orthlen,my,kk,ti1,ti2,dodepth
       real*8 vvv(0:volume-1)
       real*8 omat(0:n*n,0:n*n)
       real*8 rmat(0:n*n,0:n*n),rmat3(0:n*n,0:n*n)
       real*8 rmat2(0:n*n,0:n*n),bvec(0:n*n)
       integer depth
c
c Assume ti2 > ti1
        call unit(omat)
        call unit(rmat)
        do 10 mx = 0, n*n-1
10        bvec(mx) = 1.d0
        do 40 i = ti1, ti2
          kk = mod(i+l,l)
          call multb(omat,vvv,kk,-1)
          if(mod(i+1-ti1,orthlen) .eq. 0 .or. i .eq. ti2)then
            do 30 my = 0, n*n-1
              do 30 mx = 0, n*n-1
30              omat(mx,my) = omat(mx,my) * bvec(my)
            if(i .eq. ti2 .or. dodepth .eq. 0)then
              depth = n*n
            else if(mod(i+1-ti1,2*orthlen) .eq. 0)then
              depth = 2*n*n/3
            else
              depth = n*n/3
            endif
            if(n .eq. 2)then
              depth = 4
            endif
            call orthfacp(omat,bvec,rmat2,depth)
            call matmult(rmat2,rmat,rmat3)
            call matcop(rmat3,rmat)
          endif
40      continue
       do 35 my = 0, n*n-1
         do 35 mx = 0, n*n-1
35         omat(mx,my) = omat(mx,my) * bvec(my)
       call orthfacp(omat,bvec,rmat2,n*n)
       call matmult(rmat2,rmat,rmat3)
       call matcop(rmat3,rmat)
       return
       end
c
       subroutine getudri(vvv,ti1,ti2,dodepth,orthlen,omat,bvec,rmat)
       implicit none
        integer n,l,volume,toff
        parameter(n=4,l=48)
        parameter(volume=n*n*l,toff=n*n)
       integer i,mx,orthlen,my,kk,ti1,ti2,dodepth
       real*8 vvv(0:volume-1)
       real*8 omat(0:n*n,0:n*n)
       real*8 rmat(0:n*n,0:n*n),rmat3(0:n*n,0:n*n)
       real*8 rmat2(0:n*n,0:n*n),bvec(0:n*n)
       integer depth
c
c Assume ti2 > ti1
        call unit(omat)
        call unit(rmat)
        do 10 mx = 0, n*n-1
10        bvec(mx) = 1.d0
        do 40 i = ti2, ti1, -1
          kk = mod(i+l,l)
          call multbi(omat,vvv,kk,-1)
          if(mod(i+1-ti1,orthlen) .eq. 0 .or. i .eq. ti2)then
            do 30 my = 0, n*n-1
              do 30 mx = 0, n*n-1
30              omat(mx,my) = omat(mx,my) * bvec(my)
            if(i .eq. ti2 .or. dodepth .eq. 0)then
              depth = n*n
            else if(mod(i+1-ti1,2*orthlen) .eq. 0)then
              depth = 2*n*n/3
            else
              depth = n*n/3
            endif
            if(n .eq. 2)then
              depth = 4
            endif
            call orthfacp(omat,bvec,rmat2,depth)
            call matmult(rmat2,rmat,rmat3)
            call matcop(rmat3,rmat)
          endif
40      continue
       do 35 my = 0, n*n-1
         do 35 mx = 0, n*n-1
35         omat(mx,my) = omat(mx,my) * bvec(my)
       call orthfacp(omat,bvec,rmat2,n*n)
       call matmult(rmat2,rmat,rmat3)
       call matcop(rmat3,rmat)
       return
       end
c
       subroutine ftntok(gn,gq,ndim,maxl)
c Fourier transform g(n,l) to get g(q,l).
       implicit none
       real*8 twopi
       parameter(twopi=6.283185307179586)
       integer ndim,maxl
       real*8 gn(0:ndim/2,0:ndim/2,0:maxl)
       real*8 gq(0:ndim/2,0:ndim/2,0:maxl),cfac,cs(0:400)
       integer mx,my,lx,ly,ti,lxp,lyp,first
 
       if(first .ne. 1273) then
         if(ndim .gt. 20) write(76,*)'Help#1 in ftntok!,',ndim
         first = 1273
         do 100 lx = 0,ndim**2
100        cs(lx) = cos(twopi/ndim*lx) 
       endif
       do 10 mx = 0, ndim/2
        do 10 my = 0, ndim/2
         do 15 ti = 0, maxl
15        gq(mx,my,ti) = 0.d0
         do 10 lx = 0, ndim-1
          lxp = min(lx,ndim-lx)
          do 10 ly = 0, ndim-1
           lyp = min(ly,ndim-ly)
           cfac = cs(mx*lx) * cs(my*ly)
           do 10 ti = 0, maxl
            gq(mx,my,ti) = gq(mx,my,ti)+ cfac * gn(lxp,lyp,ti)
10     continue
       return
       end
c
       subroutine ftltow(gl,gw,ndim,maxl,dtau,bose,nmax)
c Fourier transform g(n,l) to get g(n,w) (Fermi frequencies).
       implicit none
       real*8 twopi
       parameter(twopi=6.283185307179586)
       integer ndim,maxl,bose,nmax
       real*8 gl(0:ndim/2,0:ndim/2,0:maxl),dtau
       real omega
       real*8 larray(1000),terpray(1000),tauray(1000)
       real*8 temp,temp2,temp3,rti2
       complex*16 gw(0:ndim/2,0:ndim/2,0:maxl)
       integer ti,ti2,intrat,mx,my
       parameter (intrat=4)
c dtau/intrat is the spacing used for the tau integration.
 
       do 20 mx = 0, ndim/2
        do 20 my = 0, ndim/2
         do 10 ti = 0, maxl
           larray(ti+1) = gl(mx,my,ti)
10         tauray(ti+1) = ti*dtau
         call spline(tauray,larray,maxl+1,2e30,2e30,terpray)
c Fourier transform over time to get g(n,w).
         do 20 ti = 0, nmax
          if(bose .eq. 0)then
            omega = (ti+0.5)*twopi/(maxl)
          else
            omega = ti*twopi/(maxl)
          endif
          call splint(tauray,larray,terpray,maxl+1,0.d0,temp)
          call splint(tauray,larray,terpray,maxl+1,maxl*dtau,temp2)
          rti2 = (intrat*maxl-1.d0)/intrat*dtau
          call splint(tauray,larray,terpray,maxl+1,rti2,temp3)
          gw(mx,my,ti) = temp*dtau / intrat / 3.d0
     1       + temp2*exp(cmplx(0.0,1.0)*omega*maxl)*dtau/intrat/3.d0
     2       + temp3*exp(cmplx(0.0,1.0)*omega*rti2/dtau)*
     3                                 dtau/intrat*(4.d0/3.d0)
          do 20 ti2 = 1,intrat*maxl-3,2
           call splint(tauray,larray,terpray,maxl+1,
     1                                     ti2*dtau/intrat,temp)
           gw(mx,my,ti) = gw(mx,my,ti) + temp*(4.d0/3.d0) *
     1        exp(cmplx(0.0,1.0)*omega*ti2/intrat)*dtau/intrat
           call splint(tauray,larray,terpray,maxl+1,
     1                                (ti2+1)*dtau/intrat,temp)
           gw(mx,my,ti) = gw(mx,my,ti) + temp*(2.d0/3.d0) *
     1       exp(cmplx(0.0,1.0)*omega*(ti2+1)/intrat)*dtau/intrat
20      continue
       return
       end
c
       subroutine ftltowp(pairl,pairw,maxl,k,dtau,nmax)
c Fourier transform pair(n,l) to get pair(n,w) (Fermi frequencies).
       implicit none
       real*8 twopi
       parameter(twopi=6.283185307179586)
       integer maxl,k,j,nmax
       real*8 pairl(10,9,0:maxl),dtau
       real omega
       complex*16 pairw(10,9,0:maxl)
       integer ti,ti2,intrat
       real*8 larray(1000),terpray(1000),tauray(1000)
       real*8 temp,temp2,temp3,rti2
       parameter (intrat=4)
c dtau/intrat is the spacing used for the tau integration.
 
       do 20 j = 1,9
        do 10 ti = 0, maxl
         larray(ti+1) = pairl(k,j,ti)
10       tauray(ti+1) = ti*dtau
        call spline(tauray,larray,maxl+1,2e30,2e30,terpray)
        do 20 ti = 0, nmax
         omega = ti*twopi/(maxl)
         call splint(tauray,larray,terpray,maxl+1,0.d0,temp)
         call splint(tauray,larray,terpray,maxl+1,maxl*dtau,temp2)
         rti2 = (intrat*maxl-1.d0)/intrat*dtau
         call splint(tauray,larray,terpray,maxl+1,rti2,temp3)
          pairw(k,j,ti) = temp*dtau / intrat / 3.d0
     1       + temp2*exp(cmplx(0.0,1.0)*omega*maxl)*dtau/intrat/3.d0
     2       + temp3*exp(cmplx(0.0,1.0)*omega*rti2/dtau)*
     3                                 dtau/intrat*(4.d0/3.d0)
          do 20 ti2 = 1,intrat*maxl-3,2
           call splint(tauray,larray,terpray,maxl+1,
     1                                     ti2*dtau/intrat,temp)
           pairw(k,j,ti) = pairw(k,j,ti) + temp*(4.d0/3.d0) *
     1        exp(cmplx(0.0,1.0)*omega*ti2/intrat)*dtau/intrat
           call splint(tauray,larray,terpray,maxl+1,
     1                                (ti2+1)*dtau/intrat,temp)
           pairw(k,j,ti) = pairw(k,j,ti) + temp*(2.d0/3.d0) *
     1       exp(cmplx(0.0,1.0)*omega*(ti2+1)/intrat)*dtau/intrat
20      continue
       return
       end
c


       SUBROUTINE SPLINE(X,Y,N,YP1,YPN,Y2)
       PARAMETER (NMAX=400)
       REAL*8 X(N),Y(N),Y2(N),U(NMAX)
c       DIMENSION X(N),Y(N),Y2(N),U(NMAX)
       IF (YP1.GT..99E30) THEN
         Y2(1)=0.
         U(1)=0.
       ELSE
         Y2(1)=-0.5
         U(1)=(3./(X(2)-X(1)))*((Y(2)-Y(1))/(X(2)-X(1))-YP1)
       ENDIF
       DO 10 I=2,N-1
         SIG=(X(I)-X(I-1))/(X(I+1)-X(I-1))
         P=SIG*Y2(I-1)+2.
         Y2(I)=(SIG-1.)/P
         U(I)=(6.*((Y(I+1)-Y(I))/(X(I+1)-X(I))-(Y(I)-Y(I-1))
     *         /(X(I)-X(I-1)))/(X(I+1)-X(I-1))-SIG*U(I-1))/P
10      continue
       IF (YPN.GT..99E30) THEN
         QN=0.
         UN=0.
       ELSE
         QN=0.5
         UN=(3./(X(N)-X(N-1)))*(YPN-(Y(N)-Y(N-1))/(X(N)-X(N-1)))
       ENDIF
       Y2(N)=(UN-QN*U(N-1))/(QN*Y2(N-1)+1.)
       DO 20 K=N-1,1,-1
              Y2(K)=Y2(K)*Y2(K+1)+U(K)
20     continue
       RETURN
       END
c
       SUBROUTINE SPLINT(XA,YA,Y2A,N,X,Y)
	REAL*8 XA(N),YA(N),Y2A(N),X,Y
c       DIMENSION XA(N),YA(N),Y2A(N)
       KLO=1
       KHI=N
1       IF (KHI-KLO.GT.1) THEN
         K=(KHI+KLO)/2
         IF(XA(K).GT.X)THEN
           KHI=K
         ELSE
           KLO=K
         ENDIF
       GOTO 1
       ENDIF
       H=XA(KHI)-XA(KLO)
       IF (H.EQ.0.) then
         write(6,*) 'Bad XA input.'
         stop
       endif
       A=(XA(KHI)-X)/H
       B=(X-XA(KLO))/H
       Y=A*YA(KLO)+B*YA(KHI)+
     *     ((A**3-A)*Y2A(KLO)+(B**3-B)*Y2A(KHI))*(H**2)/6.
       RETURN
       END
c**********************************
c getgp -- routine to include pivoting in the orthogonalization. 1/30/89
       subroutine getgp(vvv,ti,gmat,sgndet,deta)
       implicit none
        integer n,l,volume,toff
        parameter(n=4,l=48)
        parameter(volume=n*n*l,toff=n*n)
       integer orthlen,doauto,torth
       real*8 eorth,difflim,errpam
        common/getgpam/eorth,difflim,errpam,orthlen,doauto,
     1                torth
       integer i,j,mx,k,my,kk,ti,depth
       real*8 vvv(0:volume-1)
       real*8 gmat(0:n*n,0:n*n)
 
       real*8 rmat(0:n*n,0:n*n),rmat2(0:n*n,0:n*n),rmat3(0:n*n,0:n*n)
       real*8 oimat(0:n*n,0:n*n),omat(0:n*n,0:n*n)
       real*8 omat2(0:n*n,0:n*n),omat3(0:n*n,0:n*n)
       common/getgmats/rmat,rmat2,rmat3,oimat,omat,omat2,omat3
 
       real*8 bvec(0:n*n),sgndet,det,deta
       integer ipvt(toff),nbig
 
c
        nbig=0
        call unit(omat)
        call unit(rmat)
        do 10 mx = 0, toff-1
10        bvec(mx) = 1.0
        do 40 i = 0, l-1
          kk = mod(i+ti+l,l)
          call multb(omat,vvv,kk,-1)
          if(mod(i+1,orthlen) .eq. 0 .or. i .eq. l-1)then
            do 30 my = 0, toff-1
              do 30 mx = 0, toff-1
30              omat(mx,my) = omat(mx,my) * bvec(my)
              depth = n*n
            call orthfacp(omat,bvec,rmat2,depth)
            call matmult(rmat2,rmat,rmat3)
            call matcop(rmat3,rmat)
          endif
40      continue
c One last orthogonalization to make sure omat is orthogonal.
        do 35 my = 0, toff-1
          do 35 mx = 0, toff-1
35          omat(mx,my) = omat(mx,my) * bvec(my)
        call orthfacp(omat,bvec,rmat2,depth)
        call matmult(rmat2,rmat,rmat3)
        call matcop(rmat3,rmat)
c
c oimat = omat^-1
        call invertr(omat,oimat,det,toff+1,toff,2,ipvt)
        deta=det
        sgndet = 1.0
        if(det .lt. 0.0) sgndet = -1.0
c rmat2 = rmat^-1
        call invertr(rmat,rmat2,det,toff+1,toff,2,ipvt)
        deta=deta*det

c Calculate omat2 = D + omat^-1 R^-1
c Dimension of omat2 is toff-nbig.
        call zeromat(omat2)
        do 65 i = 0, toff-1
65        omat2(i,i) = bvec(i)
        do 60 k = 0, toff-1
         do 60 j = 0, toff-1
          do 60 i = 0, toff-1
60         omat2(i,j)=omat2(i,j)+oimat(i,k)*rmat2(k,j)
c calculate inverse of omat2, put in omat3
        call invertr(omat2,omat3,det,toff+1,toff-0,2,ipvt)
        deta=deta*det
        if(det .lt. 0.0) sgndet = -sgndet

        call zeromat(omat)
        do 80 i = 0,toff-1
         do 80 j = 0,toff-1
          do 80 k = 0,toff-1
80         omat(i,k) = omat(i,k) + omat3(i-nbig,j-nbig)*oimat(j,k)

c Multiply to get gmat=rmat2*omat.
c Remember:omat(j,k) = 0 for j < 0
        call zeromat(gmat)
        do 95 i = 0,toff-1
          do 95 j = nbig,toff-1
            do 95 k = 0,toff-1
95            gmat(i,k) = gmat(i,k) + rmat2(i,j)*omat(j,k)

       return
       end
 
c**********************************
c getgp -- routine to include pivoting in the orthogonalization. 1/30/89
       subroutine debug(vvv,ti,gmat,   det)
       implicit none
        integer n,l,volume,toff
        parameter(n=4,l=48)
        parameter(volume=n*n*l,toff=n*n)
       integer orthlen,doauto,torth
       real*8 eorth,difflim,errpam
        common/getgpam/eorth,difflim,errpam,orthlen,doauto,
     1                torth
       integer i,j,mx,k,my,kk,ti,depth
       real*8 vvv(0:volume-1)
       real*8 gmat(0:n*n,0:n*n)
 
       real*8 rmat(0:n*n,0:n*n),rmat2(0:n*n,0:n*n),rmat3(0:n*n,0:n*n)
       real*8 oimat(0:n*n,0:n*n),omat(0:n*n,0:n*n)
       real*8 omat2(0:n*n,0:n*n),omat3(0:n*n,0:n*n)
       common/getgmats/rmat,rmat2,rmat3,oimat,omat,omat2,omat3
 
       real*8 bvec(0:n*n),sgndet,det
       integer ipvt(toff),nbig
              
        call unit(omat)
        call unit(rmat)
        do 10 mx = 0, toff-1
10        bvec(mx) = 1.d0
        do 40 i = 0, l-1
          kk = mod(i+ti+l,l)
          call multb(omat,vvv,kk,-1)
          if(mod(i+1,orthlen) .eq. 0 .or. i .eq. l-1)then
            do 30 my = 0, toff-1
              do 30 mx = 0, toff-1
30              omat(mx,my) = omat(mx,my) * bvec(my)
            if(i .eq. l-1 )then
              depth = n*n
            else if(mod(i+1,2*orthlen) .eq. 0)then
              depth = 2*n*n/3
            else
              depth = n*n/3
            endif
            if(n .eq. 2)then
              depth = 4
            endif
c 
c	    depth = toff
            call orthfacp(omat,bvec,rmat2,depth)
            call matmult(rmat2,rmat,rmat3)
            call matcop(rmat3,rmat)
          endif
40      continue
c One last orthogonalization to make sure omat is orthogonal.
        do 35 my = 0, toff-1
          do 35 mx = 0, toff-1
35          omat(mx,my) = omat(mx,my) * bvec(my)
        call orthfacp(omat,bvec,rmat2,depth)
        call matmult(rmat2,rmat,rmat3)
        call matcop(rmat3,rmat)
c
c oimat = omat^-1
        call invertr(omat,oimat,det,toff+1,toff,2,ipvt)
        sgndet = 1.d0
        if(det .lt. 0.d0) sgndet = -1.d0
c rmat2 = rmat^-1
        call invertr(rmat,rmat2,det,toff+1,toff,2,ipvt)
 
c Calculate omat2 = D + omat^-1 R^-1
c Dimension of omat2 is toff-nbig.
        call zeromat(omat2)
        do 65 i = 0, toff-1
65        omat2(i,i) = bvec(i)
        do 60 k = 0, toff-1
         do 60 j = 0, toff-1
          do 60 i = 0, toff-1
60         omat2(i,j)=omat2(i,j)+oimat(i,k)*rmat2(k,j)
c calculate inverse of omat2, put in omat3
        call invertr(omat2,omat3,det,toff+1,toff-0,2,ipvt)
        if(det .lt. 0.d0) sgndet = -sgndet
 
        call zeromat(omat)
        do 80 i = 0,toff-1
         do 80 j = 0,toff-1
          do 80 k = 0,toff-1
80         omat(i,k) = omat(i,k) + omat3(i-nbig,j-nbig)*oimat(j,k)
 
c Multiply to get gmat=rmat2*omat.
c Remember:omat(j,k) = 0 for j < 0
        call zeromat(gmat)
        do 95 i = 0,toff-1
          do 95 j = nbig,toff-1
            do 95 k = 0,toff-1
95            gmat(i,k) = gmat(i,k) + rmat2(i,j)*omat(j,k)
 
       return
       end
 
c****** NEW VERSION ******* (with pivoting)
       subroutine orthogp(mat,rmat,depth)
       implicit none
       integer n,l
       parameter(n=4,l=48)
 
       integer i,j,k,depth, done(0:n*n-1),ib
       real*8 mat(0:n*n,0:n*n),rmat(0:n*n,0:n*n),maxsize,colsize
c
       call unit(rmat)
       do 10 i = 0, n*n-1
10       done(i) = 0
c done(i)=0 means the other columns havent been orthogonalized to col i yet.
       do 20 i = 0, min(depth,n*n)-1
c First find the column which is largest and hasn't been done yet.
c That col. will be "ib".
         maxsize = -10.d0
         ib = -1
         do 30 j = 0, n*n-1
           if(done(j) .eq. 0)then
             colsize = 0.d0
             do 40 k=0, n*n-1
40             colsize = colsize + mat(k,j)**2
             if(colsize .gt. maxsize) then
               maxsize = colsize
               ib = j
             endif
           endif
30       continue
         done(ib) = 1
         call donorm(mat(0,ib),rmat(ib,ib))
c Now orthogonalize all the other columns to ib.
         do 50 j = 0, n*n-1
           if(done(j) .eq. 0) then
             do 60 k=0, n*n-1
60             rmat(ib,j) = rmat(ib,j) + mat(k,ib)*mat(k,j)
             do 70 k=0, n*n-1
70             mat(k,j) = mat(k,j) - rmat(ib,j)*mat(k,ib)
           endif
50       continue
 
20     continue
 
       return
       end
 
c********************************
 
       subroutine orthfacp(mat,dvec,rmat,depth)
       implicit none
       integer n,l
       parameter(n=4,l=48)
       integer i,j,depth
       real*8 mat(0:n*n,0:n*n),rmat(0:n*n,0:n*n),dvec(0:n*n-1)
       real*8 temp
c
       call orthogp(mat,rmat,depth)
       do 10 i = 0, n*n-1
         dvec(i) = rmat(i,i)
         if(dvec(i) .eq. 0.d0) then
           rmat(i,i) = 1.d0
         else if(dvec(i) .ne. 1.d0) then
             temp = 1.d0 / dvec(i)
             do 20 j = 0, n*n-1
20             rmat(i,j) = rmat(i,j) * temp
         endif
10     continue
       return
       end
c**************************subroutine getgtau*********************
c This subroutine uses Eugene's idea for unequal time Green's functions
c with Steve's factorizing out diagonal matrices trick.
       subroutine getgtau(vvv,ti,dti,gmat,gmat2)
       implicit none
        integer n,l,volume,toff
        parameter(n=4,l=48)
        parameter(volume=n*n*l,toff=n*n)
       integer orthlen,doauto,torth
       real*8 eorth,difflim,errpam
       common/getgpam/eorth,difflim,errpam,orthlen,doauto,
     1                torth
       integer i,j,ti,dti,depth
       real*8 vvv(0:volume-1)
       real*8 gmat(0:n*n,0:n*n), gmat2(0:n*n,0:n*n)
 
       real*8 r1(0:n*n,0:n*n),r1i(0:n*n,0:n*n),r21i(0:n*n,0:n*n)
       real*8 r12i(0:n*n,0:n*n)
       real*8 oimat(0:n*n,0:n*n),u1(0:n*n,0:n*n)
       real*8 u2i1(0:n*n,0:n*n),u2i(0:n*n,0:n*n)
       real*8 u1i2(0:n*n,0:n*n)
       real*8 u2(0:n*n,0:n*n),r2(0:n*n,0:n*n)
       real*8 imat(0:n*n,0:n*n)
 
       real*8 d1(0:n*n),d2(0:n*n)
       real*8 dbar1i(0:n*n),dbar2i(0:n*n),dtil1(0:n*n),dtil2(0:n*n)
       real*8 det
       integer ipvt(toff)
 
        depth = n*n
        call getudr(vvv,ti+dti,ti+l-1,0,orthlen,u1,d1,r1)
        call getudri(vvv,ti,ti+dti-1,0,orthlen,u2,d2,r2)
c r1i = r1^-1 ; r21i = r2*r1^-1; r12i = r1*r2^-1
        call matcop(r1,imat)
        call invertr(imat,r1i,det,toff+1,toff,2,ipvt)
        call matmult(r2,r1i,r21i)
        call matcop(r21i,imat)
        call invertr(imat,r12i,det,toff+1,toff,2,ipvt)
c Calculate u2i=u2^-1; u2i1 = u2^-1*u1; u1i2 = u1^-1*u2
        call transp(u2,u2i)
        call matmult(u2i,u1,u2i1)
        call matcop(u2i1,imat)
        call invertr(imat,u1i2,det,toff+1,toff,2,ipvt)
c To make the last inversion stable, factor out some diagonal matrices.
c Define  dbar1i(i) = 1.0 / max(d1(i),1); dtil1(i)=dbar1i(i)*d1(i)=min(1,d1(i))
c same for dbar2i, etc.
c The formula for G(ti+dti,ti) is:
c  G = r1i dbar1i (dtil2 r21i dbar1i + dbar2i u2i1 dtil1)^-1 dbar2i u2i
c
        do 76 j = 0, n*n-1
          dbar1i(j) = 1.d0 / max(d1(j),1.d0)
          dtil1(j) = d1(j) * dbar1i(j)
          dbar2i(j) = 1.d0 / max(d2(j),1.d0)
76        dtil2(j) = d2(j) * dbar2i(j)
        do 75 j = 0, n*n-1
          do 75 i = 0, n*n-1         
            r21i(i,j) = dtil2(i)*r21i(i,j)*dbar1i(j)
75          u2i1(i,j) = dbar2i(i)*u2i1(i,j)*dtil1(j)
c Now add the two matrices and invert
        do 65 i = 0, n*n-1
          do 65 j = 0, n*n-1       
            gmat(i,j)=u2i1(i,j)+r21i(i,j)
65      continue    
        call invertr(gmat,oimat,det,toff+1,toff,2,ipvt)
        do 85 i = 0, n*n-1
          do 85 j = 0, n*n-1       
            oimat(i,j)=oimat(i,j)*dbar1i(i)*dbar2i(j)
85      continue    
c  Multiply on the left by r1i and the right by u2i
        call matmult(r1i,oimat,r21i)
        call matmult(r21i,u2i,gmat)
c
c Now do G(ti,ti+dti)
        do 176 j = 0, n*n-1
          d1(j) = 1.d0 / d1(j)
          d2(j) = 1.d0 / d2(j)
          dbar1i(j) = 1.d0 / max(d1(j),1.d0)
          dtil1(j) = d1(j) * dbar1i(j)
          dbar2i(j) = 1.d0 / max(d2(j),1.d0)
176       dtil2(j) = d2(j) * dbar2i(j)
        do 175 j = 0, n*n-1
          do 175 i = 0, n*n-1         
            r12i(i,j) = dtil2(j)*r12i(i,j)*dbar1i(i)
175         u1i2(i,j) = dbar2i(j)*u1i2(i,j)*dtil1(i)
        do 165 i = 0, n*n-1
          do 165 j = 0, n*n-1       
            gmat2(i,j)=u1i2(i,j)+r12i(i,j)
165      continue    
        call invertr(gmat2,oimat,det,toff+1,toff,2,ipvt)
        do 185 i = 0, n*n-1
          do 185 j = 0, n*n-1       
            oimat(i,j)=oimat(i,j)*dbar1i(j)*dbar2i(i)
185      continue    
c  Multiply on the left by r1i and the right by u2i
        call matmult(u2,oimat,r21i)
        call matmult(r21i,r1,gmat2)
 
       return
       end
           
c**************************INVERTR(BBB,BBBIN,D,IDET,IPVT)**********************

        subroutine invertr(b,c,d,dim,toff,idet,ipvt)   
                   
c       This subroutine finds the inverse of the matrix b
c       using Gaussian elimination with partial pivoting.
c       b is the input matrix.
c       c is the inverse of b.
c       d is the determinent of b.
c       toff is the size of the matrix
c       If idet=0 the inverse of b is returned but not the determinent.
c       If idet=1 the determinent is returned,but not the inverse.
c       If idet=2 both the inverse and the determinent are returned.
c                                     
        implicit none
        integer dim,toff,ipvt(dim) 
        real*8 b(dim,dim),c(dim,dim),d,tb            
        integer iwn,nm1,k,kp1,mb,i,j,idet,kb,km1
 
        ipvt(toff)=1
       iwn=1
       nm1=toff-1
c                                   
c       Finding the pivot
c
       do 35 k=1,nm1
       kp1=k+1
       mb=k
       do 15 i=kp1,toff
       if (dabs(b(i,k)).gt.dabs(b(mb,k))) mb=i
15       continue
       ipvt(k)=mb                                          
       tb=b(mb,k)
       if (mb.ne.k) then
       ipvt(toff)=-ipvt(toff)
       b(mb,k)=b(k,k)
       b(k,k)=tb
       endif
       if (dabs(tb).lt.1.d-10) then
       iwn=0
       go to 35
       endif
c
c       Computing multipliers
c
       do 20 i=kp1,toff
20       b(i,k)=-b(i,k)/tb
c
c       Interchange and elimination by columns
c
       do 30 j=kp1,toff
       tb=b(mb,j)
       if (mb.ne.k) then
       b(mb,j)=b(k,j)
        b(k,j)=tb
       endif
       do 25 i=kp1,toff
25       b(i,j)=b(i,j)+b(i,k)*tb
cray       call saxpy(toff-k,tb,b(kp1,k),1,b(kp1,j),1)
30       continue
35       continue
       if (iwn.eq.0) then
       write (6,*) 'Warning: Pivot=0'
       d=0.d0    
       go to 150
       endif
       if (idet.gt.0) then
c      
c       Compute the determinent
c
       d=ipvt(toff)
       do 40 i=1,toff
40       d=d*b(i,i)
       endif
       if (idet.ne.1) then
c
c       Compute the inverse
c       Ax=b goes to A~x=Lb~
c
       do 50 i=1,toff
       do 45 j=1,toff
45       c(i,j)=0.d0
50       c(i,i)=1.d0
       do 100 j=1,toff
       do 60 k=1,nm1
       kp1=k+1
       mb=ipvt(k)
       tb=c(mb,j)
       if (mb.ne.k) then
       c(mb,j)=c(k,j)
       c(k,j)=tb
       endif
       do 55 i=kp1,toff
55       c(i,j)=c(i,j)+b(i,k)*tb
cray       call saxpy(toff-k,tb,b(kp1,k),1,c(kp1,j),1)
60       continue
c
c       Inverting A~
c
       do 70 kb=1,nm1
       km1=toff-kb
       k=km1+1
       c(k,j)=c(k,j)/b(k,k)
       tb=-c(k,j)
       do 65 i=1,km1    
65       c(i,j)=c(i,j)+b(i,k)*tb
cray       call saxpy(km1,tb,b(1,k),1,c(1,j),1)
70       continue
       c(1,j)=c(1,j)/b(1,1)
100       continue
       endif
150       continue
       return
       end
     

c       USE THESE COMMENTED OUT LINES IF REAL*8 DESIRED.
      REAL*8 FUNCTION RAN2(IDUM)
      IMPLICIT REAL*8(A-H,O-Z)
c      FUNCTION RAN2(IDUM)
      save
      PARAMETER (M=714025,IA=1366,IC=150889,RM=1.4005112E-6)
      DIMENSION IR(97)
      DATA IFF /0/
      IF(IDUM.LT.0.OR.IFF.EQ.0)THEN
        IFF=1
        IDUM=MOD(IC-IDUM,M)
        DO 11 J=1,97
          IDUM=MOD(IA*IDUM+IC,M)
          IR(J)=IDUM
11      CONTINUE
        IDUM=MOD(IA*IDUM+IC,M)
        IY=IDUM
      ENDIF
      J=1+(97*IY)/M
      IF(J.GT.97.OR.J.LT.1)PAUSE
      IY=IR(J)
      RAN2=IY*RM
      IDUM=MOD(IA*IDUM+IC,M)
      IR(J)=IDUM
      RETURN
      END



       
