
c	MODIFIED JULY 1997 TO INCLUDE t(ij) IN MEASUREMENTS
c       MODIFIED FEB. 2000 TO INCLUDE PH SYMMETRIC SITE DISORDER
c       MODIFIED JUNE 2001 TO INCLUDE NN HOPPING DISORDER IN A
c                          PARALLEL (ZEEMAN) MAGNETIC FIELD

c	MODIFICATION OF CAREY'S negu96b.f BACK TO +U, AND THENCE
c	     TO GENERAL LATTICE STRUCTURES.
c	(*) CHECK CAREY'S CODE AGAINST MY OLD negu94c6.f AT U=0.
c	(*) CHECK CAREY'S CODE AGAINST MY OLD negu94c6.f AT U=4.
c	(*) MODIFY BACK TO +U AND CHECK AGAINST exact9adp.f AT U=0.
c	(*) MODIFY BACK TO +U AND CHECK AGAINST exact9adp.f AT U=4.
c	(*) CHECK AGAINST MARTIN'S VERSION.
c       (*) PUT IN DIFFERENT MULTT.
c	    NB:  MANY DEAD SECTIONNS OF CODE REMOVED AND CODE REARRANGED.
c	(*) RECOMPILE AND LOOK FOR ERRORS:  IBM.
c	(*) RECOMPILE AND LOOK FOR ERRORS:  SUN.
c	(*) RECOMPILE AND LOOK FOR ERRORS:  DEC.
c	(*) RUN FTNCHEK.
c	(*) RUN EXECUTION TIME BOUNDS CHECKER.
c	(*) EXAMINE GETEK ETC IN DEBUGGER.
c	(*) CHECK AGAINST exact9adp.f AT U=0. [INCLUDE TROTTER ANALYSIS.]
c	(*) CHECK AGAINST exact9adp.f AT U=4. [INCLUDE TROTTER ANALYSIS.]
c	(*) CHECK AGAINST MARTIN'S VERSION. 
c       (*) TIMING TRIALS
c	(*) CHECK BOND/SITE RANDOMNESS AT U=0 WITH INDEPENDENT DIAG.
c       ++  version posu96c ++
c       ( ) PRINT OUT 4L BINS OF TIME DEPENDENT DATA FOR ANALYTIC CONTINUATION.
c	( ) PUT IN SUBROUTINES FOR 2D ANDERSON LATTICE GEOMETRY.
c	( ) PUT IN SUBROUTINES FOR BCC ANDERSON LATTICE GEOMETRY.
c	( ) PUT IN SUBROUTINES FOR CUO2 GEOMETRY.

      program dishbpar
cpd  A previous version of this program is called posu96c.f in Davis 
cpd  and was used on Cray SDSC.

      implicit none
      include 'paramp.dat'
      integer warms,sweeps,nwrap,iran,iran0
      common/integers/warms,sweeps,nwrap,iran,iran0
 
      double precision u,mu,bpar,delt,dtau,gam,lambda
      common/couple/u,mu,bpar,delt,dtau,gam,lambda
 
c     common block for vectors
      double precision hub(0:volume-1)

      double precision vup(0:volume-1),vdn(0:volume-1)
      double precision mphase(0:toff-1)
      common/vectors/hub,vup,vdn,mphase

      double precision sgnup,sgndn,sgn,nup,ndn,nud,ke
      double precision saf,saf2,sferro,sfer2,grfun(0:n/2,0:n/2),
     1     safsq,safsq2
      double precision spinxx(0:n/2,0:n/2),spinzz(0:n/2,0:n/2)
      double precision pairmat(-1:1,-1:1,-1:1,-1:1),errrat,asgnp
      double precision fpairmat(-1:1,-1:1,-1:1,-1:1,0:n/2,0:n/2)
      double precision afpairmat(-1:1,-1:1,-1:1,-1:1,0:n/2,0:n/2)
      double precision asgnup,asgndn,asgn,anup,andn,anud,ake,
     1     asafsq,asafsq2
      double precision asaf,asaf2,asferro,asfer2,agrfun(0:n/2,0:n/2)
      double precision aspinxx(0:n/2,0:n/2),aspinzz(0:n/2,0:n/2)
      double precision apairmat(-1:1,-1:1,-1:1,-1:1)

      double precision den(0:n/2,0:n/2,0:1),aden(0:n/2,0:n/2,0:1)
      integer nmeas0,redo,noredo,dopair,nmeasp,numpair,nmeaspf
      common/m0var/ sgnup,sgndn,sgn,nup,ndn,nud,ke,safsq,safsq2,
     1       saf,saf2,sferro,sfer2,grfun,spinxx,spinzz,aspinxx,
     2       aspinzz,asafsq,asafsq2,asgnup,asgndn,asgn,anup,andn,
     3       anud,ake,den,aden,asaf,asaf2,asferro,asfer2,agrfun,
     4       errrat,asgnp,pairmat,apairmat,
     6       fpairmat,afpairmat
      common/m0vari/dopair,nmeasp,numpair,nmeas0,redo,noredo,nmeaspf
 
      double precision gnl(0:n/2,0:n/2,0:l),agnl(0:n/2,0:n/2,0:l)
      double precision chinl(0:n/2,0:n/2,0:l),achinl(0:n/2,0:n/2,0:l)
      double precision dent(0:n/2,0:n/2,0:l),adent(0:n/2,0:n/2,0:l)
      double precision chinlz(0:n/2,0:n/2,0:l),achinlz(0:n/2,0:n/2,0:l)
      double complex chinwz(0:n/2,0:n/2,0:l)
      double complex chiqwz(0:n/2,0:n/2,0:l)
      double precision cnl(0:n/2,0:n/2,0:l),acnl(0:n/2,0:n/2,0:l)
      double precision pstau(0:l),pdtau(0:l),psxtau(0:l)
      double complex lam(0:n*n,0:n*n),alam(0:n*n,0:n*n)
      double complex ag1w(0:n*n,0:n*n)
      double precision apstau(0:l),apdtau(0:l),apsxtau(0:l),asgnt
      integer nmeast
      common/mtauvar/gnl,agnl,chinl,achinl,chinlz,achinlz,cnl,
     1   acnl,pstau,pdtau,psxtau,apstau,apdtau,apsxtau,lam,alam,
     2	 ag1w,asgnt,dent,adent
      common/mtauvari/nmeast
 
      integer i,k,accept,reject,wraps,ione
      integer j1,j2,mx,my,lx,ly,ti,mpx,mpy,kx,ky,kmx,kmy
      double precision aveval,errval,aveval2,errval2,avedens,
     1                 sigmadc 
      double precision bpairv(10,9),gpairl(10,9,0:l),
     1     bfpairv(10,9,0:n/2,0:n/2)
 
      double precision gql(0:n/2,0:n/2,0:l)
      double complex gnw(0:n/2,0:n/2,0:l),gqw(0:n/2,0:n/2,0:l)
c      double complex sigma(0:n/2,0:n/2,0:l)
      double complex cnw(0:n/2,0:n/2,0:l),cqw(0:n/2,0:n/2,0:l)
      double precision chiql(0:n/2,0:n/2,0:l)
      double precision dentq(0:n/2,0:n/2,0:l)
      double precision chiqlz(0:n/2,0:n/2,0:l),cql(0:n/2,0:n/2,0:l)
      double complex chinw(0:n/2,0:n/2,0:l),chiqw(0:n/2,0:n/2,0:l)
 
      double precision bgnl(10,0:n/2,0:n/2,0:l),bgql(10,0:n/2,0:n/2,0:l)
      double complex bgnw(0:n/2,0:n/2,0:l),bgqw(0:n/2,0:n/2,0:l)
      double complex bsigma(0:n/2,0:n/2,0:l)
      double complex bcnw(0:n/2,0:n/2,0:l),bcqw(0:n/2,0:n/2,0:l)
      double precision bchi(0:n/2,0:n/2,0:l),
     1     bchiz(0:n/2,0:n/2,0:l)
      double precision bchinl(10,0:n/2,0:n/2,0:l),
     1     bchiql(10,0:n/2,0:n/2,0:l)
      double precision bdent(10,0:n/2,0:n/2,0:l),
     1     bdentq(10,0:n/2,0:n/2,0:l)
      double precision bchinlz(10,0:n/2,0:n/2,0:l),
     1     bchiqlz(10,0:n/2,0:n/2,0:l)
      double precision bcnl(10,0:n/2,0:n/2,0:l),
     1     bcql(10,0:n/2,0:n/2,0:l)
      double complex bchinw(0:n/2,0:n/2,0:l),bchiqw(0:n/2,0:n/2,0:l)
      double precision bpstau(10,0:l),bpdtau(10,0:l),bpsxtau(10,0:l)
 
      double complex sgnw(0:n/2,0:n/2,0:l),sgqw(0:n/2,0:n/2,0:l)
      double complex ssigma(0:n/2,0:n/2,0:l)
      double complex scnw(0:n/2,0:n/2,0:l),scqw(0:n/2,0:n/2,0:l)
      double complex schinw(0:n/2,0:n/2,0:l),schiqw(0:n/2,0:n/2,0:l)
 
      double precision bsgnup(10),bsgndn(10),bsgn(10),bnup(10),bndn(10)
      double precision bsgnp(10),bsgnt(10),bone(10)
      double precision bntot(10),bnud(10),bke(10),benergy(10)
      double precision bsaf(10),bsaf2(10),bsferro(10),bsfer2(10)
      double precision bgrfun(10,0:n/2,0:n/2),bsafsq(10),bsafsq2(10)
      double precision bden(10,0:n/2,0:n/2,0:1)
      double precision bspinxx(10,0:n/2,0:n/2), bspinzz(10,0:n/2,0:n/2)
 
      double precision pairgsus(-1:1,-1:1,-1:1,-1:1,0:l)
      double precision gmatup(0:n*n,0:n*n),gmatdn(0:n*n,0:n*n)
c      double precision epsk
      integer nmax
      character*79 outname
      character*80 rstring,gstring
      double precision alpha1,alpha2
      integer start
      common/startup/start
      double precision waves(-1:1,-1:1,9)

c     initialization flags
      integer ifbigt,ifftnt,ifgetc,ifgetr,ifmeas,ifgetd
      common/iflags/ifbigt,ifftnt,ifgetc,ifgetr,ifmeas,ifgetd

c     occupancy histogram variables (uses parameter histno)
      integer hist1(histno*3+2),hist2(histno*3+2),hist3(histno*3+2)
      common/histogram/hist1,hist2,hist3

c     variable so we can avoid divide by zero errors when tausk=0
      integer tauskp
      common/tauskdiv/tauskp

      data waves/0.0, 0.0, 0.0, 0.0,1.0,0.0,0.0, 0.0, 0.0,
     1             0.0, 0.5, 0.0, 0.5,0.0,0.5,0.0, 0.5, 0.0,
     2             0.0,-0.5, 0.0, 0.5,0.0,0.5,0.0,-0.5, 0.0,
     3             0.5, 0.0, 0.5, 0.0,0.0,0.0,0.5, 0.0, 0.5,
     4            -0.5, 0.0, 0.5, 0.0,0.0,0.0,0.5, 0.0,-0.5,
     5             0.0, 0.0, 0.0,-0.5,0.0,0.5,0.0, 0.0, 0.0,
     6             0.0,-0.5, 0.0, 0.0,0.0,0.0,0.0, 0.5, 0.0,
     7            -0.5, 0.0, 0.0, 0.0,0.0,0.0,0.0, 0.0, 0.5,
     8             0.0, 0.0,-0.5, 0.0,0.0,0.0,0.5, 0.0, 0.0/

c     the work starts here

c     set the initialization flags
      ifbigt = 0
      ifftnt = 0
      ifgetc = 0
      ifgetr = 0
      ifmeas = 0
      ifgetd = 0

c     initialize
      if(tausk.eq.0) then
         tauskp = 1
      else
         tauskp = tausk
      endif

c     get the input
      write (6,*) 'enter outname'
      read(5,1789)  outname
 1789 format(A40)
c
c     open output files
c                 (except unit 72; opened only at very end)
      rstring = 'r'//outname
      gstring = 'g'//outname
      open(unit=66, status='new',file=rstring)
c 
      if(tausk.ne.0) then
         open(unit=68, status='new',file=gstring)
      endif
 
c     write header
c     first to "r" file
      write(66,*) 'Version dishbpar'
      write(66,*) ' '
      write(66,235) '    n = ',n
      write(66,235) '    l = ',l
 235  format(a10,i4)
c      write(66,*) ' '
      write(66,235) 'tausk = ',tausk
      write(66,235) 'doall = ',doall
      write(66,235) 'denswp= ',denswp
      write(66,235) 'histno= ',histno
c     then to "g" file
      if(tausk.ne.0) then
         write(68,*) 'Version dishbpar'
         write(68,*) ' '
         write(68,235) '    n = ',n
         write(68,235) '    l = ',l
c     write(68,*) ' '
         write(68,235) 'tausk = ',tausk
         write(68,235) 'doall = ',doall
         write(68,235) 'denswp= ',denswp
         write(68,235) 'histno= ',histno
c     write (68,*) '  '
      endif

c     set up the one body part of H and its exponential.
      call getek

      call readin()

c     set random number generator
      write (66,*) 'after disorder set, iran = ',iran
      write (66,*) '  '
      if(tausk.ne.0) then
         write (68,*) 'after disorder set, iran = ',iran
         write (68,*) '  '
      endif

      if (start.ne.999) call autoset()
      if (start.ne.999) call ranlat()

 234  format(a12,f16.12)
      write (66,234) 'Using mu = ',mu
      write (66,*) ' '
      if(tausk.ne.0) then
         write (68,234) 'Using mu = ',mu
         write (68,*) ' '
      endif

c     DEBUGGING:
c     lr=1
c     call   unit(gmatup)
c     call  multt(gmatup,lr)
c     call multti(gmatup,lr)
c     lr=-1
c     call   unit(gmatup)
c     call  multt(gmatup,lr)
c     call multti(gmatup,lr)
c     lr=1
c     call   unit(gmatup)
c     call  multt(gmatup,lr)
c     call multti(gmatup,-lr)
c     lr=1
c     call   unit(gmatup)
c     call  multt(gmatup,-lr)
c     call multti(gmatup,lr)

 
      nmax = l / 2
      accept = 1
      reject = 0
      do 141 ti = 0, l
         do 141 j1 = 0, n/2
            do 141 j2 = 0, n/2
               bchi(j1,j2,ti) = 0.d0
               bchiz(j1,j2,ti) = 0.d0
               bgnw(j1,j2,ti) = 0.d0
               bgqw(j1,j2,ti) = 0.d0
               bsigma(j1,j2,ti) = 0.d0
               bcnw(j1,j2,ti) = 0.d0
               bcqw(j1,j2,ti) = 0.d0
               bchinw(j1,j2,ti) = 0.d0
               bchiqw(j1,j2,ti) = 0.d0
               sgnw(j1,j2,ti) = 0.d0
               sgqw(j1,j2,ti) = 0.d0
               ssigma(j1,j2,ti) = 0.d0
               scnw(j1,j2,ti) = 0.d0
               scqw(j1,j2,ti) = 0.d0
               schinw(j1,j2,ti) = 0.d0
               schiqw(j1,j2,ti) = 0.d0
 141  continue

c     initialize the histogram bins
      do j1 = 1,histno*3+2
         hist1(j1) = 0
         hist2(j1) = 0
         hist3(j1) = 0
      enddo

      call setvup()
 
c     perform warmup sweeps
      wraps = 0
      redo = 0
      noredo = 1
      call getgp(vup,0,gmatup,sgnup,1)
c      if(mu .ne. 0.d0.or.bpar.ne.0.d0) then
         call getgp(vdn,0,gmatdn,sgndn,-1)
c         sgndn = sgnup
c      endif
      do 10 i=1,warms
         write(6,*)'Starting warmup sweep ',i
         call sweep(gmatup,gmatdn,accept,reject,wraps)
 10   continue
      write(66,*)'after warmups, accept ratio is ',
     1                float(accept)/(accept+reject)
      write(66,*)'gamma is ',gam
      write(66,*)'redo ratio is ',
     1                float(redo)/(redo+noredo)
      write(66,*) ' '
 
c     perform measurement sweeps
      call setvup()
      call zeroas()
      do 20 i=1,sweeps
         write(6,*)'Starting measurement sweep ',i
         if(mod(i,10) .eq. 0) then
            write(6,*)'accept, redo ratios are ',
     1         float(accept)/(accept+reject),
     2         float(redo)/(noredo+redo)
         endif
         call sweep(gmatup,gmatdn,accept,reject,wraps)
         if(tausk.ne.0) then
            if(mod(i,tauskp) .eq. 0) then
               call meastausq
            endif
         endif

         if(mod(i,sweeps/10) .eq. 0) then
            write(66,*)'Finished measurement sweep ',i
               if(nmeasp .eq. 0)nmeasp = 1
               if(nmeast .eq. 0)nmeast = 1
            write(66,9045) asgn/nmeas0,asgnp/nmeasp,
     1         float(accept)/(accept+reject),
     2         float(redo)/(noredo+redo)
 9045       format('  asgn, asgnp: ',2(f8.3,' '),
     1          ';accept,redo ratios: ',2('  ',f9.4))
            k = (i*10)/sweeps
            bsgnup(k) = asgnup / nmeas0
            bsgndn(k) = asgndn / nmeas0
            bsgn(k) = asgn / nmeas0
	    bsgnp(k) = asgnp / nmeasp
	    bsgnt(k) = asgnt / nmeast
	    bone(k) = 1.d0
            bnup(k)= anup / nmeas0
            bndn(k)= andn / nmeas0
            bntot(k)= (anup+andn) / nmeas0
            bsaf(k)= asaf / nmeas0
            bsafsq(k)= dsqrt(dabs(asafsq / nmeas0))
            bsferro(k)= asferro / nmeas0
            bsfer2(k)= asfer2 / nmeas0
            bsaf2(k)= asaf2 / nmeas0
            bke(k)= ake / nmeas0
            benergy(k)= (ake+u*anud) / nmeas0
            bnud(k)= anud / nmeas0

            bsafsq2(k)= dsqrt(dabs(asafsq2) / nmeas0)
            do 41 j1 = 0, n/2
               do 41 j2 = 0, n/2
                  bgrfun(k,j1,j2) = agrfun(j1,j2) / nmeas0
                  bden(k,j1,j2,0) = aden(j1,j2,0) / nmeas0
                  bden(k,j1,j2,1) = aden(j1,j2,1) / nmeas0
                  bspinxx(k,j1,j2) = aspinxx(j1,j2) / nmeas0
                  bspinzz(k,j1,j2) = aspinzz(j1,j2) / nmeas0
c      Set gnl and chinl for this tenth of the run.
                  do 41 ti = 0, l
                     gnl(j1,j2,ti) = agnl(j1,j2,ti)/nmeast
                     gql(j1,j2,ti) = 0.d0
                     chinl(j1,j2,ti) = achinl(j1,j2,ti)/nmeast
                     dent(j1,j2,ti) = adent(j1,j2,ti)/nmeast
                     chinlz(j1,j2,ti) = achinlz(j1,j2,ti)/nmeast
                     cnl(j1,j2,ti) = acnl(j1,j2,ti)/nmeast
                     chiql(j1,j2,ti) = 0.d0
                     dentq(j1,j2,ti) = 0.d0
                     chiqlz(j1,j2,ti) = 0.d0
                     cql(j1,j2,ti) = 0.d0
 41         continue
	    if (nmeast.ne.1) then
               call ftntok(chinlz,chiqlz,n,l)
               call ftntok(chinl,chiql,n,l)
               call ftntok(dent,dentq,n,l)
               call ftntok(gnl,gql,n,l)
               call ftntok(cnl,cql,n,l)
	    endif
	    if (nmeast.ne.1) then
               call ftltow(gnl,gnw,n,l,dtau,0,nmax)
               call ftltow(gql,gqw,n,l,dtau,0,nmax)
               call ftltow(chinl,chinw,n,l,dtau,1,nmax)
               call ftltow(chinlz,chinwz,n,l,dtau,1,nmax)
               call ftltow(chiql,chiqw,n,l,dtau,1,nmax)
               call ftltow(chiqlz,chiqwz,n,l,dtau,1,nmax)
               call ftltow(cnl,cnw,n,l,dtau,1,nmax)
               call ftltow(cql,cqw,n,l,dtau,1,nmax)
	    endif

            do 586 mx = 0, n/2
              do 586 my = 0, n/2
                do 586 ti = 0, l
                  bgnl(k,mx,my,ti) = gnl(mx,my,ti)
                  bgql(k,mx,my,ti) = gql(mx,my,ti)
                  bchinl(k,mx,my,ti) = chinl(mx,my,ti)
                  bdent(k,mx,my,ti) = dent(mx,my,ti)
                  bchinlz(k,mx,my,ti) = chinlz(mx,my,ti)
                  bchiql(k,mx,my,ti) = chiql(mx,my,ti)
                  bdentq(k,mx,my,ti) = dentq(mx,my,ti)
                  bchiqlz(k,mx,my,ti) = chiqlz(mx,my,ti)
                  bcnl(k,mx,my,ti) = cnl(mx,my,ti)
                  bcql(k,mx,my,ti) = cql(mx,my,ti)
                    if(ti .le. nmax)then
                      bgnw(mx,my,ti) = bgnw(mx,my,ti) + gnw(mx,my,ti)
                      bgqw(mx,my,ti) = bgqw(mx,my,ti) + gqw(mx,my,ti)
                      bcnw(mx,my,ti) = bcnw(mx,my,ti) + cnw(mx,my,ti)
                      bcqw(mx,my,ti) = bcqw(mx,my,ti) + cqw(mx,my,ti)
                      bchinw(mx,my,ti) = bchinw(mx,my,ti)+
     1                     chinw(mx,my,ti)
                      bchiqw(mx,my,ti) = bchiqw(mx,my,ti)+
     1                     chiqw(mx,my,ti)

                      alpha1= dreal(gnw(mx,my,ti))**2
                      alpha2=dimag(gnw(mx,my,ti))**2
                      sgnw(mx,my,ti) = sgnw(mx,my,ti)+
     1                     dcmplx(alpha1,alpha2)
 
                      alpha1= dreal(gqw(mx,my,ti))**2
                      alpha2=dimag(gqw(mx,my,ti))**2
                      sgqw(mx,my,ti) = sgqw(mx,my,ti)+
     1                     dcmplx(alpha1,alpha2)
 
                      alpha1= dreal(cnw(mx,my,ti))**2
                      alpha2=dimag(cnw(mx,my,ti))**2
                      scnw(mx,my,ti) = scnw(mx,my,ti)+
     1                     dcmplx(alpha1,alpha2)
 
                      alpha1= dreal(cqw(mx,my,ti))**2
                      alpha2=dimag(cqw(mx,my,ti))**2
                      scqw(mx,my,ti) = scqw(mx,my,ti)+
     1                     dcmplx(alpha1,alpha2)

                      alpha1= dreal(chinw(mx,my,ti))**2
                      alpha2=dimag(chinw(mx,my,ti))**2
                      schinw(mx,my,ti) = schinw(mx,my,ti)+
     1                     dcmplx(alpha1,alpha2)
 
                      alpha1= dreal(chiqw(mx,my,ti))**2
                      alpha2=dimag(chiqw(mx,my,ti))**2
                      schiqw(mx,my,ti) = schiqw(mx,my,ti)+
     1                     dcmplx(alpha1,alpha2)


c	comment this out since it needs the free particle dispersion
c                 epsk = (-2.d0)*t*(dcos(twopi*mx/n)+dcos(twopi*my/n))
c
c     Since the following line is commented you are getting sigma=0.
carey  the next lines were the offending lines which prevented
c     only time-independent measurements on the cray
c     and the division: 1.d0/gqw(mx,my,ti) was the offending statement
c                 if(tausk.ne.0) then
c                    sigma(mx,my,ti)=
c     1                   dcmplx(-epsk+mu,(ti+0.5d0)*twopi/l/dtau)
c     2                   + 1.d0/gqw(mx,my,ti)
c                    bsigma(mx,my,ti) = bsigma(mx,my,ti) + 
c     1                   sigma(mx,my,ti)
c                    alpha1= dreal(sigma(mx,my,ti))**2
c                    alpha2=dimag(sigma(mx,my,ti))**2
c                    ssigma(mx,my,ti) = ssigma(mx,my,ti)+
c     1                   dcmplx(alpha1,alpha2)
c                 endif
              endif
 586       continue

            if(dopair .eq. 1) then
               do 8372 ti = 0,l
                  bpstau(k,ti) = 0.d0
                  bpdtau(k,ti) = 0.d0
                  bpsxtau(k,ti) = 0.d0
 8372          continue
               do 8379 mx = -1, 1
                  do 8379 my = -1, 1
                     do 8379 mpx = -1, 1
                        do 8379 mpy = -1, 1
                           do 8356 ti = 0,l
 8356                         pairgsus(mx,my,mpx,mpy,ti) = 0.d0
                            do 8379 lx = 0, n-1
                               do 8379 ly = 0, n-1
                                  kx = min(lx,n-lx)
                                  ky = min(ly,n-ly)
                                  kmx = mod(n+n+lx-(mx-mpx),n)
                                  kmx = min(kmx,n-kmx)
                                  kmy = mod(n+n+ly-(my-mpy),n)
                                  kmy = min(kmy,n-kmy)
                                     do 8379 ti = 0,l
c     the next line just won't fit
      pairgsus(mx,my,mpx,mpy,ti)=pairgsus(mx,my,mpx,mpy,ti)
     1            +gnl(kmx,kmy,ti)*gnl(kx,ky,ti)
 8379         continue
	      do 2374 ti = 0, l
                 bpstau(k,ti) = bpstau(k,ti) + apstau(ti)/nmeast
                 bpdtau(k,ti) = bpdtau(k,ti) + apdtau(ti)/nmeast
                 bpsxtau(k,ti) = bpsxtau(k,ti) + apsxtau(ti)/nmeast
 2374         continue
              do 2378 j1 = 1,9
                 bpairv(k,j1) = 0.d0
                    do 2338 ti = 0, l
                       gpairl(k,j1,ti) = 0.d0
 2338               continue
               do 2378 mx = -1, 1
                do 2378 my = -1, 1
                 do 2378 mpx = -1, 1
                  do 2378 mpy = -1, 1
                    bpairv(k,j1)=bpairv(k,j1)+waves(mx,my,j1)*
     1              apairmat(mx,my,mpx,mpy)*waves(mpx,mpy,j1)/nmeasp
                     do 2378 ti = 0, l
                       gpairl(k,j1,ti)=gpairl(k,j1,ti)+waves(mx,my,j1)*
     1                     pairgsus(mx,my,mpx,mpy,ti)*waves(mpx,mpy,j1)
 2378         continue
              do 2379 j1 = 1,9
               do 2379 kx=0,n/2
                do 2379 ky=0,n/2
                 bfpairv(k,j1,kx,ky) = 0.d0
                  do 2379 mx = -1, 1
                   do 2379 my = -1, 1
                    do 2379 mpx = -1, 1
                     do 2379 mpy = -1, 1
               bfpairv(k,j1,kx,ky)=bfpairv(k,j1,kx,ky)+waves(mx,my,j1)*
     1         afpairmat(mx,my,mpx,mpy,kx,ky)*waves(mpx,mpy,j1)/nmeaspf
 2379          continue

            endif
            call zeroas()
         endif
 20   continue

      write(66,*) ' '
      write(66,*) ' '
      write(66,*)'At end, redo ratio is ',float(redo)/(redo+noredo)
      write(66,*)'gamma is ',gam
      write(66,*)'Acceptance ratio = ',float(accept)/(accept+reject)
      call geterr(bsgnup,bone,aveval,errval)
 8877 format(a32,f9.6,' +- ',f9.6)
      write(66,8877)'Average up sign = ',aveval,errval
      call geterr(bsgndn,bone,aveval,errval)
      write(66,8877)'Average dn sign = ',aveval,errval
      call geterr(bsgn,bone,aveval,errval)
      write(66,8877)'Average total sign = ',aveval,errval
      asgn=aveval
      call geterr(bntot,bsgn,aveval,errval)
      write(66,8877)'Average density = ',aveval,errval
cpd99 define avedens for compact output
      avedens=aveval
      call geterr(bnup,bsgn,aveval,errval)
      write(66,8877)'Average up occupancy = ',aveval,errval
      call geterr(bndn,bsgn,aveval,errval)
      write(66,8877)'Average dn occupancy = ',aveval,errval
      call geterr(benergy,bsgn,aveval,errval)
      write (66,*) '**please note average energy and kinetic energy**'
      write (66,*) '**  include mu terms, both uniform and random  **'
      write(66,8877)'Average Energy = ',aveval,errval
      call geterr(bke,bsgn,aveval,errval)
      write(66,8877)'Average Kinetic Energy = ',aveval,errval
      call geterr(bnud,bsgn,aveval,errval)
      write(66,8877)'Average Nup*Ndn = ',aveval,errval
      call geterr(bsaf,bsgn,aveval,errval)
      write(66,8877)'AF correlation function (xx) = ',aveval,errval
      call geterr(bsaf2,bsgn,aveval,errval)
      write(66,8877)'AF correlation function (zz) = ',aveval,errval
      call geterr(bsferro,bsgn,aveval,errval)
      write(66,8877)'Ferro corr. func. (xx) = ',aveval,errval
      call geterr(bsfer2,bsgn,aveval,errval)
      write(66,8877)'Ferro corr. func. (zz) = ',aveval,errval
      write(66,*)'Green''s function:'
 8878 format(i4,i4,f10.6,' +- ',f10.6)
      do 677 j1 = 0, n/2
         do 677 j2 = 0, n/2
            call geterr(bgrfun(1,j1,j2),bsgn,aveval,errval)
            if(j2 .ge. j1)
     1        write(66,8878)j1,j2,-aveval,errval
            grfun(j1,j2) = aveval
 677  continue
      write(66,*)'density-density correlation fn: (up-up,up-dn)'
      do 977 j1 = 0, n/2
         do 977 j2 = j1, n/2
            call geterr(bden(1,j1,j2,0),bsgn,aveval,errval)
            call geterr(bden(1,j1,j2,1),bsgn,aveval2,errval2)
            write(66,1987)j1,j2,aveval,errval,aveval2,errval2
 1987       format(2i4,2('    ',f12.6,' +- ',f12.6))
 977  continue
      write(66,*)'zz Spin correlation function:'
      do 687 j1 = 0, n/2
         do 687 j2 = j1, n/2
            call geterr(bspinzz(1,j1,j2),bsgn,aveval,errval)
            write(66,8878)j1,j2,aveval,errval
 687  continue
      write(66,*)'xx Spin correlation function:'
      do 697 j1 = 0, n/2
         do 697 j2 = j1, n/2
            call geterr(bspinxx(1,j1,j2),bsgn,aveval,errval)
            write(66,8878)j1,j2,aveval,errval
 697  continue
      if(doall .eq. 1)then
         write(66,*)'lambda(p=(pi,0),pprime):'
      endif
      call geterr(bsafsq,bsgn,aveval,errval)
      write(66,8877)'RMS AF corr. func. (xx) = ',aveval,errval
      call geterr(bsafsq2,bsgn,aveval,errval)
      write(66,8877)'RMS AF corr. func. (zz) = ',aveval,errval
      write(66,*)' '
 
      if(tausk.ne.0) then
         write(68,*)'signs: '
         write(68,101)(bsgnt(j1),j1=1,10)
 101     format(5(f10.6))
         write(68,*) ' '
         write(68,*)'G(nx,ny,ti):'
 1234    format(i4,10f7.4)
 102     format(i4,'  ',f14.10,' +- ',f14.10)
         do 987 j1 = 0, n/2
            do 987 j2 = j1, n/2
               write(68,*)'nx = ',j1,' ny = ',j2
               do 987 ti = 0,l
                  call geterr(bgnl(1,j1,j2,ti),bsgn,aveval,errval)
                  write(68,102)ti, -aveval, errval
 987           continue
 
      write(68,*)'G(qx,qy,ti):'
      do 9871 j1 = 0, n/2
         do 9871 j2 = j1, n/2
            write(68,*)'qx = ',j1,' qy = ',j2
            do 9871 ti = 0,l
               call geterr(bgql(1,j1,j2,ti),bsgn,aveval,errval)
               write(68,102)ti, -aveval, errval
 9871 continue
 
      write(68,*)'G(nx,ny,omega), omega = (n+.5d0) 2 pi T :'
      do 9872 j1 = 0, n/2
         do 9872 j2 = j1, n/2
            write(68,*)'nx = ',j1,' ny = ',j2
            do 9872 ti = 0,nmax
               call geterr2( dreal(bgnw(j1,j2,ti)),
     1              dreal(sgnw(j1,j2,ti)),aveval,errval)
               call geterr2(dimag(bgnw(j1,j2,ti)),
     1              dimag(sgnw(j1,j2,ti)),aveval2,errval2)
               write(68,1235)ti,-aveval,errval,-aveval2,errval2
 1235          format(i5,' ( ',f12.7,' +- ',f12.7,' ) + i * ( '
     1              ,f12.7,' +- ',f12.7,' ) ')
 9872 continue
 
      write(68,*)'G(qx,qy,omega):'
      do 9873 j1 = 0, n/2
         do 9873 j2 = j1, n/2
            write(68,*)'qx = ',j1,' qy = ',j2
            do 9873 ti = 0,nmax
               call geterr2( dreal(bgqw(j1,j2,ti)),
     1              dreal(sgqw(j1,j2,ti)),aveval,errval)
               call geterr2(dimag(bgqw(j1,j2,ti)),
     1              dimag(sgqw(j1,j2,ti)),aveval2,errval2)
               write(68,1235)ti,-aveval,errval,-aveval2,errval2
 9873 continue

c      write(68,*)'SIGMA(qx,qy,omega):'
c      do 9874 j1 = 0, n/2
c         do 9874 j2 = j1, n/2
c            write(68,*)'qx = ',j1,' qy = ',j2
c            do 9875 ti = nmax,0,-1
c               call geterr2( dreal(bsigma(j1,j2,ti)), 
c     1              dreal(ssigma(j1,j2,ti)),aveval,errval)
c               call geterr2(dimag(bsigma(j1,j2,ti)),
c     1              dimag(ssigma(j1,j2,ti)),aveval2,errval2)
c               write(68,1235)-(2*ti+1), aveval,errval,
c     1              -aveval2,errval2
c 9875       continue
c            do 9874 ti = 0,nmax
c               call geterr2( dreal(bsigma(j1,j2,ti)),
c     1              dreal(ssigma(j1,j2,ti)),aveval,errval)
c               call geterr2(dimag(bsigma(j1,j2,ti)),
c     1              dimag(ssigma(j1,j2,ti)),aveval2,errval2)
c               write(68,1235)(2*ti+1), aveval,errval,aveval2,
c     1              errval2
c 9874 continue

      write(68,*)'chi_xx(nx,ny,ti):'
      do 417 j1 = 0, n/2
         do 417 j2 = j1, n/2
            write(68,*)'nx = ',j1,' ny = ',j2
            do 417 ti = 0,l
               call geterr(bchinl(1,j1,j2,ti),bsgn,aveval,errval)
               write(68,102)ti, aveval, errval
 417  continue

      write(68,*)'den(nx,ny,ti):'
      do 4117 j1 = 0, n/2
         do 4117 j2 = j1, n/2
            write(68,*)'nx = ',j1,' ny = ',j2
            do 4117 ti = 0,l
               call geterr(bdent(1,j1,j2,ti),bsgn,aveval,errval)
               write(68,102)ti, aveval, errval
 4117 continue
 
      write(68,*)'chi_zz(nx,ny,ti):'
      do 418 j1 = 0, n/2
         do 418 j2 = j1, n/2
            write(68,*)'nx = ',j1,' ny = ',j2
            do 418 ti = 0,l
               call geterr(bchinlz(1,j1,j2,ti),bsgn,aveval,errval)
               write(68,102)ti, aveval, errval
 418  continue

      write(68,*)'den(qx,qy,ti):'
      do 41711 j1 = 0, n/2
         do 41711 j2 = j1, n/2
            write(68,*)'qx = ',j1,' qy = ',j2
            do 41711 ti = 0,l
               call geterr(bdentq(1,j1,j2,ti),bsgn,aveval,errval)
               write(68,102)ti, aveval, errval
41711 continue
 
      write(68,*)'chi_xx(qx,qy,ti):'
      do 4171 j1 = 0, n/2
         do 4171 j2 = j1, n/2
            write(68,*)'qx = ',j1,' qy = ',j2
            do 4171 ti = 0,l
               call geterr(bchiql(1,j1,j2,ti),bsgn,aveval,errval)
               write(68,102)ti, aveval, errval
 4171 continue
 
      write(68,*)'chi_zz(qx,qy,ti):'
      do 4174 j1 = 0, n/2
         do 4174 j2 = j1, n/2
            write(68,*)'qx = ',j1,' qy = ',j2
            do 4174 ti = 0,l
               call geterr(bchiqlz(1,j1,j2,ti),bsgn,aveval,errval)
               write(68,102)ti, aveval, errval
 4174 continue

      write(68,*)'chi,z:'
 103  format('chi,z(',i4,',',i4,')=',f14.10,f14.10)
      do 9991 j1=0,n/2
         do 9991 j2=j1,n/2
            do 9992 k=1,10
               bchi(j1,j2,l/2)=bchi(j1,j2,l/2)+bchiql(k,j1,j2,l/2)/10.d0
               bchiz(j1,j2,l/2)=bchiz(j1,j2,l/2)+
     1              bchiqlz(k,j1,j2,l/2)/10.d0
 9992       continue
            write(68,103)j1,j2,bchi(j1,j2,l/2),bchiz(j1,j2,l/2)
 9991 continue
      write(68,*)'current(nx,ny,ti):'
      do 419 j1 = 0, n/2
         do 419 j2 = 0, n/2
            write(68,*)'nx = ',j1,' ny = ',j2
            do 419 ti = 0,l
               call geterr(bcnl(1,j1,j2,ti),bsgn,aveval,errval)
               write(68,102)ti, aveval, errval
 419  continue
 
      write(68,*)'current(qx,qy,ti):'
      do 4175 j1 = 0, n/2
         do 4175 j2 = 0, n/2
            write(68,*)'qx = ',j1,' qy = ',j2
            do 4175 ti = 0,l
               call geterr(bcql(1,j1,j2,ti),bsgn,aveval,errval)
               write(68,102)ti, aveval, errval
cpd99 added to output sigmadc
               if (ti.EQ.l/2.AND.j1.EQ.0.AND.j2.EQ.0) then 
               sigmadc =  -aveval*dtau*l*dtau*l/3.141592654
               endif
 4175 continue

      write(68,*)'current(nx,ny,omega), omega=2n*pi T:'
      do j1 = 0, n/2
         do j2 = 0, n/2
            write(68,*)'nx = ',j1,' ny = ',j2
            do ti = 0,nmax
               call geterr2( dreal(bcnw(j1,j2,ti)),
     1              dreal(scnw(j1,j2,ti)),aveval,errval)
               call geterr2(dimag(bcnw(j1,j2,ti)),
     1              dimag(scnw(j1,j2,ti)),aveval2,errval2)
               write (68,1235) ti,aveval, errval, aveval2, errval2
            enddo
         enddo
      enddo
      write(68,*)'current(qx,qy,omega), omega=2n*pi T:'
      do j1 = 0, n/2
         do j2 = 0, n/2
            write(68,*)'qx = ',j1,' qy = ',j2
            do ti = 0,nmax
               call geterr2( dreal(bcqw(j1,j2,ti)),
     1              dreal(scqw(j1,j2,ti)),aveval,errval)
               call geterr2(dimag(bcqw(j1,j2,ti)),
     1              dimag(scqw(j1,j2,ti)),aveval2,errval2)
               write (68,1235) ti,aveval, errval, aveval2, errval2
            enddo
         enddo
      enddo
 
      write(68,*)'chi(nx,ny,omega), omega = 2 n pi T :'
      do 4172 j1 = 0, n/2
         do 4172 j2 = j1, n/2
            write(68,*)'nx = ',j1,' ny = ',j2
            do 4172 ti = 0,nmax
               call geterr2( dreal(bchinw(j1,j2,ti)),
     1              dreal(schinw(j1,j2,ti)),aveval,errval)
               call geterr2(dimag(bchinw(j1,j2,ti)),
     1              dimag(schinw(j1,j2,ti)),aveval2,errval2)
               write(68,1235)ti, aveval,errval,aveval2,errval2
 4172 continue

      write(68,*)'chi(qx,qy,omega):'
      do 4173 j1 = 0, n/2
         do 4173 j2 = j1, n/2
            write(68,*)'qx = ',j1,' qy = ',j2
            do 4173 ti = 0,nmax
               call geterr2( dreal(bchiqw(j1,j2,ti)), 
     1              dreal(schiqw(j1,j2,ti)),aveval,errval)
               call geterr2(dimag(bchiqw(j1,j2,ti)),
     1              dimag(schiqw(j1,j2,ti)),aveval2,errval2)
               write(68,1235)ti, aveval,errval,aveval2,errval2
 4173 continue
      endif
c     this endif ends the if(tausk.ne.0) statement above
        
      if(dopair .eq. 1) then
 1236    format(2(f12.7,' +- ',f12.7,'   '))
         write(66,*) 's-wave: (corr. fn, no vertex)'
         call geterr(bpairv(1,1),bsgnp,aveval,errval)
         call geterr(gpairl(1,1,0),bsgnp,aveval2,errval2)
         write(66,1236)aveval,errval,aveval2,errval2
         
         write(66,*) 'sx-wave: '
         call geterr(bpairv(1,2),bsgnp,aveval,errval)
         call geterr(gpairl(1,2,0),bsgnp,aveval2,errval2)
         write(66,1236)aveval,errval,aveval2,errval2
         
         write(66,*) 'd-wave: '
         call geterr(bpairv(1,3),bsgnp,aveval,errval)
         call geterr(gpairl(1,3,0),bsgnp,aveval2,errval2)
         write(66,1236)aveval,errval,aveval2,errval2

 106     format(i4,i4,f14.6,' +- ',f14.6)
         write(66,*) 's-wave pair correlation function = '
         do 2233 j1=0,n/2
            do 2233 j2=j1,n/2
               call geterrp(bfpairv(1,1,j1,j2),aveval,errval)
               write(66,106)j1,j2,aveval/asgn,errval/asgn
 2233 continue
         write(66,*) 'Extended s-wave pair correlation function = '
         do 2234 j1=0,n/2
            do 2234 j2=j1,n/2
               call geterrp(bfpairv(1,2,j1,j2),aveval,errval)
               write(66,106)j1,j2,aveval/asgn,errval/asgn
 2234   continue
c     possible kjr output point
        write(66,*) 'd-wave pair correlation function = '
        do 2235 j1=0,n/2
           do 2235 j2=j1,n/2
              call geterrp(bfpairv(1,3,j1,j2),aveval,errval)
               write(66,106)j1,j2,aveval/asgn,errval/asgn
 2235   continue
        write(66,*) 'sxx-wave pair correlation function = '
        do 2236 j1=0,n/2
           do 2236 j2=j1,n/2
              call geterrp(bfpairv(1,4,j1,j2),aveval,errval)
               write(66,106)j1,j2,aveval/asgn,errval/asgn
 2236   continue
        write(66,*) 'dxx-wave pair correlation function = '
        do 2237 j1=0,n/2
           do 2237 j2=j1,n/2
              call geterrp(bfpairv(1,5,j1,j2),aveval,errval)
               write(66,106)j1,j2,aveval/asgn,errval/asgn
 2237   continue
        write(66,*) 'px-wave pair correlation function = '
        do 2238 j1=0,n/2
           do 2238 j2=j1,n/2
              call geterrp(bfpairv(1,6,j1,j2),aveval,errval)
               write(66,106)j1,j2,aveval/asgn,errval/asgn
 2238   continue
        write(66,*) 'py-wave pair correlation function = '
        do 2239 j1=0,n/2
           do 2239 j2=j1,n/2
              call geterrp(bfpairv(1,7,j1,j2),aveval,errval)
               write(66,106)j1,j2,aveval/asgn,errval/asgn
 2239   continue
        write(66,*) 'pxy-wave pair correlation function = '
        do 2240 j1=0,n/2
           do 2240 j2=j1,n/2
              call geterrp(bfpairv(1,8,j1,j2),aveval,errval)
               write(66,106)j1,j2,aveval/asgn,errval/asgn
 2240   continue
        write(66,*) 'pyx-wave pair correlation function = '
        do 2241 j1=0,n/2
           do 2241 j2=j1,n/2
              call geterrp(bfpairv(1,9,j1,j2),aveval,errval)
               write(66,106)j1,j2,aveval/asgn,errval/asgn
 2241   continue

        if(tausk.ne.0) then
        write(68,*)'P_d(l):'
        do 7754 ti = 0,l
           call geterrp(bpdtau(1,ti),aveval,errval)
           write(68,102) ti,aveval/asgn,errval/asgn
 7754   continue
        write(68,*)'P_d(l): (novertex)'
        do 9754 ti = 0,l
           call geterrp(gpairl(1,3,ti),aveval,errval)
           write(68,102) ti,aveval/asgn,errval/asgn
9754    continue
 
        write(68,*)'P_sx(l):'
        do 7759 ti = 0,l
           call geterrp(bpsxtau(1,ti),aveval,errval)
           write(68,102) ti,aveval/asgn,errval/asgn
7759    continue
        write(68,*)'P_sx(l): (novertex)'
        do 9759 ti = 0,l
           call geterrp(gpairl(1,2,ti),aveval,errval)
           write(68,102) ti,aveval/asgn,errval/asgn
9759    continue
 
        write(68,*)'P_s(l):'
        do 7755 ti = 0,l
           call geterrp(bpstau(1,ti),aveval,errval)
           write(68,102) ti,aveval/asgn,errval/asgn
7755    continue
        write(68,*)'P_s(l): (novertex)'
        do 9755 ti = 0,l
           call geterrp(gpairl(1,1,ti),aveval,errval)
           write(68,102) ti,aveval/asgn,errval/asgn
9755    continue
      endif
      endif
 
      saf = 2.d0*grfun(0,0)
      do 3879 lx = 0, n-1
         do 3879 ly = 0, n-1
            kx = min(lx,n-lx)
            ky = min(ly,n-ly)
 3879       saf = saf -(-1)**(lx+ly)*2*grfun(kx,ky)*grfun(kx,ky)
      write(66,*)'saf with no vertex is ',saf
 
c     write out the histograms
cpd99 Following commented out PD 9/99
c      write(66,*) ' '
c      write (66,*) 'nup, single occupancy histogram data'
c      write(66,*) ' '
c 334  format(f10.6,'  ',i12)
c 335  format('counts past lower bounds: ',i12)
c 336  format('counts past upper bounds: ',i12)
c      do i=1,histno*3
c         write(66,334) dble(i-1-histno)/dble(histno),hist1(i)
c      enddo
c      write(66,335) hist1(histno*3+1)
c      write(66,336) hist1(histno*3+2)
c      write(66,*) ' '
c      write (66,*) 'ndn, single occupancy histogram data'
c      write(66,*) ' '
c      do i=1,histno*3
c         write(66,334) dble(i-1-histno)/dble(histno),hist2(i)
c      enddo
c      write(66,335) hist2(histno*3+1)
c      write(66,336) hist2(histno*3+2)
c      write(66,*) ' '
c      write (66,*) 'double occupancy histogram data'
c      write(66,*) ' '
c       do i=1,histno*3
c         write(66,334) dble(i-1-histno)/dble(histno),hist3(i)
c      enddo
c      write(66,335) hist3(histno*3+1)
c      write(66,336) hist3(histno*3+2)

c     AT VERY END WRITE OUT H-S FIELD AND EXPONENTIAL OF SITE ENERGIES:
cpd99 Following commented out PD 9/99
c      write (66,*) '  '
c      write (66,*) ' random number seed at very end '
c      write (66,*) iran
c      write (66,*) '  '
c         write (66,*) ' final h-s field '
c         do 880 i=0,volume-1
c 880        write (66,*) hub(i)
c
cpd99 Extra file for compact output (opened at end to allow 
c     different jobs to write into same file !?)
c
      open(unit=72, status='replace',file='results.dat')
cpd99 write compact output to file results.dat
 8548 format(f7.3, i4, f7.3, f7.3, f7.3, f9.6, i9, f12.6, i3)
      ione=1
      write(72,8548)bpar,l,dtau,dtau*l,mu,avedens,iran0,sigmadc,
     #      ione
 
      stop
      end
 
c********************************************
      subroutine autoset()
      implicit none
      include 'paramp.dat'

      integer orthlen,doauto
      double precision eorth,difflim
      common/getgpam/eorth,difflim,orthlen,doauto
      save /getgpam/
 
c     common block for vectors
      double precision hub(0:volume-1)
      double precision vup(0:volume-1),vdn(0:volume-1)
      double precision mphase(0:toff-1)
      common/vectors/hub,vup,vdn,mphase

      integer start
      double precision gmatup(0:n*n,0:n*n),diffup,sgnup
      double precision gmatacc(0:n*n,0:n*n)

      write(66,*)'eorth is ',eorth
      start=orthlen
      orthlen = 1
      call ranlat()
      call getgp(vup,0,gmatacc,sgnup,1)
c     oh no!! jump statements!! see the goto below
 5    continue
      do 10 orthlen = start, 2, -1
         call getgp(vup,0,gmatup,sgnup,1)
         call matdif(gmatup,gmatacc,diffup)
         if(doauto .eq. 1)
     1        write(6,*)'diffup for orthlen= ',orthlen,' is ',diffup
              if(doauto .ne. 1 .or. diffup .le. eorth)goto 40
 10   continue
      eorth = eorth * 100.d0
      write(66,*)'resetting eorth to ',eorth,' #NO_AVE'
      goto 5
c     but what is this 40 doing here?????
 40   continue
 220  format(a16,i4)
      write(6,220)'Using orthlen= ',orthlen
      write(6,*)'diffup is ',diffup,' #NO_AVE'
      write(66,220)'Using orthlen= ',orthlen
      write(66,*)'diffup is ',diffup,' #NO_AVE'
      return
      end
c*******************************biggtau********************
 
c     This subroutine gets the big g(tau,tau prime) matrix.
cpd00 Feb.2000: extra variable sigma to discrimate for spin up/down
cpd00          (see also multb,...)
 
      subroutine biggtau(vvv,bigmat,sigma)
      implicit none
      include 'paramp.dat'
 
      double precision vvv(0:volume-1),
     1     bigmat(0:toff*ndiv-1,0:toff*ndiv-1)
      double precision inmat(0:toff*ndiv-1,0:toff*ndiv-1)
      integer llim(0:ndiv)
 
      integer ti,i,j,k,sigma
      double precision bmat(0:n*n,0:n*n),det

c     initialization flags
      integer ifbigt,ifftnt,ifgetc,ifgetr,ifmeas,ifgetd
      common/iflags/ifbigt,ifftnt,ifgetc,ifgetr,ifmeas,ifgetd

      save llim
 
c     initialization of llim using the "ifbigt" variable
c     if we are going to use this one-time initialization, then we should
c     make sure to use the save llim statement above

      if(ifbigt .ne. 373) then
         ifbigt = 373
         llim(0) = -1
         do 1 k = 1,ndiv
 1          llim(k) = (l * k)/ ndiv - 1
      endif
      do 10 i = 0, toff*ndiv-1
         do 10 j = 0,toff*ndiv-1
 10         inmat(i,j) = 0.d0
      do 30 k = 1,ndiv
         call unit(bmat)
         do 40 ti = llim(k-1)+1,llim(k)
 40         call multb(bmat,vvv,ti,-1,sigma)
            if(k .ne. ndiv) then
               do 50 i = 0,toff-1
                  do 50 j = 0,toff-1
 50                  inmat(k*toff+i,(k-1)*toff+j) = -bmat(i,j)
            else
               do 60 i = 0,toff-1
                  do 60 j = 0,toff-1
 60                  inmat(i,(k-1)*toff+j) = bmat(i,j)
            endif
 30   continue
c     We add in the 1's in case ndiv = 1 .
      do 20 i = 0, toff*ndiv-1
 20      inmat(i,i) = inmat(i,i) + 1.d0
      call invertr(inmat,bigmat,det,toff*ndiv,toff*ndiv,2)
      return
      end
c**************************************
       subroutine donorm(vec,vnorm)
       implicit none
       include 'paramp.dat'
       integer i
       double precision vec(0:n*n),vnorm,temp
c
       temp = 0.d0
       do 10 i = 0, n*n-1
 10       temp = temp + vec(i)**2
       vnorm = dsqrt(temp)
       if(vnorm .ne. 0.d0) then
          temp = 1.d0 / vnorm
          do 20 i = 0, n*n-1
 20          vec(i) = vec(i) * temp
       endif
       return
       end
c
c *****************************************************************
      subroutine ftltow(gl,gw,ndim,maxl,dtau,bose,nmax)
c     Fourier transform g(n,l) to get g(n,w) (Fermi or Bose frequencies).

c     12/30/2015:
c     Different from QUEST, here first do spline interpolation on tau data
c     then perform FT to iwn
    
      implicit none
      include 'paramp.dat'
      integer maxl,bose,nmax,ndim
      double precision gl(0:ndim/2,0:ndim/2,0:maxl),dtau,omega
      double precision larray(1000),terpray(1000),tauray(1000)
      double precision temp,temp2,temp3,rti2
      double precision splarg,spzero
      double complex gw(0:ndim/2,0:ndim/2,0:maxl)
      integer mx,my,ti,ti2
c  dtau/intrat is the spacing used for the tau integration.

c  splarg -- we want to pass 2d30 to spline, pass it via splarg
      splarg = 2.d30
c  spzero -- we also want to pass 0.d0 to splint....
      spzero = 0.d0

      do 20 mx = 0, ndim/2
       do 20 my = 0, ndim/2
        do 10 ti = 0, maxl
         larray(ti+1) = gl(mx,my,ti)
 10      tauray(ti+1) = ti*dtau
        call spline(tauray,larray,maxl+1,splarg,splarg,terpray)
c Fourier transform over time to get g(n,w).
        do 20 ti = 0, nmax
         if(bose .eq. 0)then
          omega = (ti+0.5d0)*twopi/maxl
         else
          omega = ti*twopi/maxl
          endif
          call splint(tauray,larray,terpray,maxl+1,spzero,temp)
          call splint(tauray,larray,terpray,maxl+1,maxl*dtau,temp2)
          rti2 = (intrat*maxl-1.d0)/intrat*dtau
          call splint(tauray,larray,terpray,maxl+1,rti2,temp3)
          gw(mx,my,ti) = temp*dtau / intrat / 3.d0
     1       + temp2*cdexp(dcmplx(0.d0,1.d0)*omega*maxl)*
     2                                 dtau/intrat/3.d0
     3       + temp3*cdexp(dcmplx(0.d0,1.d0)*omega*rti2/dtau)*
     4                                 dtau/intrat*(4.d0/3.d0)
          do 20 ti2 = 1,intrat*maxl-3,2
           call splint(tauray,larray,terpray,maxl+1,
     1                                     ti2*dtau/intrat,temp)
           gw(mx,my,ti) = gw(mx,my,ti) + temp*(4.d0/3.d0) *
     1        cdexp(dcmplx(0.d0,1.d0)*omega*ti2/intrat)*dtau/intrat
           call splint(tauray,larray,terpray,maxl+1,
     1                                (ti2+1)*dtau/intrat,temp)
           gw(mx,my,ti) = gw(mx,my,ti) + temp*(2.d0/3.d0) *
     1       cdexp(dcmplx(0.d0,1.d0)*omega*(ti2+1)/intrat)*dtau/intrat
20      continue
       return
       end
c**************************************************************************
      subroutine ftntok(gn,gq,ndim,maxl)
c Fourier transform g(n,l) to get g(q,l).

c 12/30/2015:
c Different from QUEST, only consider cos real part of exp(i*q*l)
c since l can be +/- [1,n/2]

      implicit none
      include 'paramp.dat'
       integer ndim,maxl
       double precision gn(0:ndim/2,0:ndim/2,0:maxl)
       double precision gq(0:ndim/2,0:ndim/2,0:maxl),cfac,cs(0:400)
       integer mx,my,lx,ly,ti,lxp,lyp

c     initialization flags
      integer ifbigt,ifftnt,ifgetc,ifgetr,ifmeas,ifgetd
      common/iflags/ifbigt,ifftnt,ifgetc,ifgetr,ifmeas,ifgetd

      save cs

       if(ifftnt .ne. 1273) then
	 write(6,*)'ndim,maxl are ',ndim,maxl
         if(ndim .gt. 20) write(66,*)'Help#1 in ftntok!,',ndim
     &	,' #NO_AVE'
         ifftnt = 1273
         do 100 lx = 0,ndim**2
100        cs(lx) = dcos(twopi/ndim*lx)
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
c*******************************getcol********************
 
      subroutine getcol(bigmat,colg,vvv,k,sigma)
      implicit none
	include 'paramp.dat'
        integer llim(0:ndiv)
        double precision vvv(0:volume-1)
 
        integer ti,i,j,nn,k,kp,dk,sigma
        double precision bigmat(0:toff*ndiv-1,0:toff*ndiv-1)
        double precision colg(0:n*n,0:n*n,0:l),fermi

c     initialization flags
      integer ifbigt,ifftnt,ifgetc,ifgetr,ifmeas,ifgetd
      common/iflags/ifbigt,ifftnt,ifgetc,ifgetr,ifmeas,ifgetd

      save llim

        if(ifgetc .ne. 373) then
          ifgetc = 373
          do 1 kp = 0,ndiv
1           llim(kp) = (l * kp)/ ndiv - 1
        endif

        nn = toff
        do 134 kp = 0, ndiv-1
         dk = mod(l+llim(kp)-llim(k),l)
         fermi=1.d0
         if(kp .lt. k) fermi = -1.d0
         do 134 i = 0, n*n-1
          do 134 j = 0, n*n-1
           colg(i,j,dk) = bigmat(kp*nn+i,k*nn+j)*fermi
134     continue
        call onesub(colg(0,0,0),colg(0,0,l))
        do 138 kp = 0, ndiv-1
         do 138 ti = llim(kp)+2,(llim(kp+1)+llim(kp)+1)/2
          dk = mod(l+ti-llim(k)-1,l)
          call matcop(colg(0,0,dk-1),colg(0,0,dk))
          call multb(colg(0,0,dk),vvv,ti-1,-1,sigma)
138     continue
        do 137 kp = 1, ndiv
         do 137 ti = llim(kp),(llim(kp)+llim(kp-1)+3)/2,-1
          dk = mod(l+ti-llim(k)-1,l)
          call matcop(colg(0,0,dk+1),colg(0,0,dk))
          call multbi(colg(0,0,dk),vvv,ti,-1,sigma)
137     continue
        return
        end

c************************GETERR*****************************
        subroutine geterr(bin,signs,ave,err)
        implicit none
        double precision bin(10),signs(10),err,sumsq,ave,denom,num
        integer i
 
        denom = 0.d0
        num = 0.d0
        do 10 i = 1, 10
          num = num + bin(i)
10        denom = denom + signs(i)
        ave = 0.d0
        do 20 i = 1, 10
20        ave = ave + (num-bin(i))/(denom-signs(i))/10.0d0
	sumsq = 0.d0
        do 30 i = 1, 10
30          sumsq = sumsq + ((num-bin(i))/(denom-signs(i))-ave)**2
        err = dsqrt(0.9d0*sumsq)
	ave = num / denom
        if(err*1.0d12 .lt. dabs(ave)) err = 0.d0
        return
        end
 
c************************GETERR2*****************************
        subroutine geterr2(sum,sumsq,ave,err)
        implicit none
        double precision err,sum,sumsq,ave,ssq
 
        ave = sum * 0.1d0
        ssq = 0.1d0*sumsq - ave**2
        err = dsqrt(dabs(ssq))/3.d0
        if(err*1.0d12 .lt. dabs(ave)) err = 0.d0
        return
        end
 
c************************GETERRP*****************************
        subroutine geterrp(bin,ave,err)
        implicit none
        double precision bin(10),err,sum,sumsq,ave
        integer i

        sum = 0.d0
        do 10 i = 1, 10
10        sum = sum + bin(i)
        ave = sum * 0.1d0
        sumsq = 0.d0
        do 20 i = 1, 10
20          sumsq = sumsq + (bin(i)-ave)**2
        sumsq = 0.1d0 * sumsq
        err = dsqrt(sumsq)/3.d0
        if(err*1.0d12 .lt. dabs(ave)) err = 0.d0
        return
        end

c**********************************
c getgp -- routine to include pivoting in the orthogonalization. 1/30/89
       subroutine getgp(vvv,ti,gmat,sgndet,sigma)
       implicit none
	include 'paramp.dat'
 
        integer orthlen,doauto
        double precision eorth,difflim
        common/getgpam/eorth,difflim,orthlen,doauto
        save /getgpam/

       integer i,j,mx,k,my,kk,ti,depth,sigma
       double precision vvv(0:volume-1)
       double precision gmat(0:n*n,0:n*n)
 
       double precision rmat(0:n*n,0:n*n),rmat2(0:n*n,0:n*n),
     1      rmat3(0:n*n,0:n*n)
       double precision oimat(0:n*n,0:n*n),omat(0:n*n,0:n*n)
       double precision omat2(0:n*n,0:n*n),omat3(0:n*n,0:n*n)
 
       double precision bvec(0:n*n),sgndet,det
       integer nbig
 
	nbig = 0
        call unit(omat)
        call unit(rmat)
        do 10 mx = 0, toff-1
10        bvec(mx) = 1.d0
        do 40 i = 0, l-1
          kk = mod(i+ti+l,l)
          call multb(omat,vvv,kk,-1,sigma)
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
        call invertr(omat,oimat,det,toff+1,toff,2)
        sgndet = 1.d0
        if(det .lt. 0.d0) sgndet = -1.d0
c rmat2 = rmat^-1
        call invertr(rmat,rmat2,det,toff+1,toff,2)
 
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
        call invertr(omat2,omat3,det,toff+1,toff-0,2)
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
c*******************************getrow********************
 
        subroutine getrow(bigmat,rowg,vvv,k,sigma)
        implicit none
	include 'paramp.dat'
        integer llim(0:ndiv)
        double precision vvv(0:volume-1)
 
        integer ti,i,j,nn,k,kp,dk,sigma
        double precision bigmat(0:toff*ndiv-1,0:toff*ndiv-1)
        double precision rowg(0:n*n,0:n*n,0:l),fermi

c     initialization flags
      integer ifbigt,ifftnt,ifgetc,ifgetr,ifmeas,ifgetd
      common/iflags/ifbigt,ifftnt,ifgetc,ifgetr,ifmeas,ifgetd

      save llim

      if(ifgetr .ne. 373) then
         ifgetr = 373
         do 1 kp = 0,ndiv
 1          llim(kp) = (l * kp)/ ndiv - 1
      endif

      nn = toff
      do 134 kp = 0, ndiv-1
         dk = mod(l+llim(k)-llim(kp),l)
         fermi=1.d0
         if(kp .gt. k) fermi = -1.d0
         do 134 i = 0, n*n-1
            do 134 j = 0, n*n-1
               rowg(i,j,dk) = bigmat(k*nn+i,kp*nn+j)*fermi
 134  continue
      call onesub(rowg(0,0,0),rowg(0,0,l))
      do 138 kp = 0, ndiv-1
         do 138 ti = llim(kp)+2,(llim(kp+1)+llim(kp)+1)/2
            dk = mod(l-ti+llim(k)+1,l)
            call matcop(rowg(0,0,dk+1),rowg(0,0,dk))
            call multbi(rowg(0,0,dk),vvv,ti-1,1,sigma)
 138  continue
      do 137 kp = 1, ndiv
         do 137 ti = llim(kp),(llim(kp)+llim(kp-1)+3)/2,-1
            dk = mod(l-ti+llim(k)+1,l)
            call matcop(rowg(0,0,dk-1),rowg(0,0,dk))
            call multb(rowg(0,0,dk),vvv,ti,1,sigma)
 137  continue
      return
      end
c**************************INVERTR(b,c,d,dim,toff,idet)********************
 
         subroutine invertr(b,c,d,dim,toff,idet)
c
c       This subroutine finds the inverse of the matrix b
c       using Gaussian elimination with partial pivoting.
c       b is the input matrix.
c       c is the inverse of b.
c       d is the determinant of b.
c       toff is the size of the matrix
c       dim is the dimensioned size of the matrix
c       If idet=0 the inverse of b is returned but not the determinent.
c       If idet=1 the determinent is returned,but not the inverse.
c       If idet=2 both the inverse and the determinent are returned, but
c            only the sign of the determinant is necessarily right.
c            (Takes care of overflow problems.)
c       If idet=3 both the inverse and the determinent are returned.
c
         implicit none
        integer dim,toff,ipvt(2000)
        double precision b(dim,dim),c(dim,dim),d,tb
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
20       b(i,k)=(-b(i,k))/tb
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
30       continue
35       continue
       if (iwn.eq.0.d0) then
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
         d=d*b(i,i)
	 if(idet .eq. 2) then
c Take care of overflow probs.
	   if(dabs(d) .gt. 1.0d100) d = d * 1.0d-50
	   if(dabs(d) .lt. 1.0d-100) d = d * 1.0d50
         endif
40     continue
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
70       continue
       c(1,j)=c(1,j)/b(1,1)
100       continue
       endif
150       continue
       return
       end
c********************************
      subroutine matcop(mat1,mat2)
      implicit none
      include 'paramp.dat'
      integer i,j
      double precision mat1(0:n*n,0:n*n),mat2(0:n*n,0:n*n)
c
      do 10 i = 0, n*n-1
         do 10 j = 0, n*n-1
 10         mat2(i,j) = mat1(i,j)
      return
      end
c**********************************
      subroutine matdif(mat1,mat2,diff)
      implicit none
      include 'paramp.dat'
      integer i,j
      double precision mat1(0:n*n,0:n*n),mat2(0:n*n,0:n*n),diff
c
       diff = 0.d0
       do 10 i = 0, n*n-1
         do 10 j = 0, n*n-1
10         diff = diff + (mat1(i,j) - mat2(i,j))**2
       diff = dsqrt(diff) / (n*n)
       return
       end
 
c**********************************
      subroutine matmult(mat1,mat2,mat3)
      implicit none
      include 'paramp.dat'
      integer i,j,k
      double precision mat1(0:n*n,0:n*n),mat2(0:n*n,0:n*n)
      double precision mat3(0:n*n,0:n*n)
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
c*******************************multb(m,vvv,ti,lr,sigma)********************
 
c        This subroutine multiplies m by the matrix B(ti) on the
c        left or right depending on whether lr = -1 or 1.
c        This is the double precision version.
cpd00    Put in extra argument sigma to discriminate up/down spin
 
        subroutine multb(mat,vvv,ti,lr,sigma)
        implicit none
	include 'paramp.dat'
        double precision vvv(0:volume-1)
        double precision mat(0:n*n,0:n*n)
        integer i,j,ti,lr,sigma
 
        if(lr .eq. (-1)) then
          call multt(mat,lr,sigma)
          do 10 i = 0, n*n-1
            do 10 j = 0, n*n-1
10            mat(i,j) = mat(i,j) * vvv(i+toff*ti)
        else
          do 20 j = 0, n*n-1
            do 20 i = 0, n*n-1
20            mat(i,j) = mat(i,j) * vvv(j+toff*ti)
          call multt(mat,lr,sigma)
        endif
 
        return
        end
c*******************************multbi(m,vvv,ti,lr,sigma)********************
 
c        This subroutine multiplies m by the matrix B(ti)^-1 on the
c        left or right depending on whether lr = -1 or 1.
c        This is the double precision version.
cpd00    Put in extra argument sigma to discriminate up/down spin
 
        subroutine multbi(mat,vvv,ti,lr,sigma)
        implicit none
	include 'paramp.dat'
        double precision vvv(0:volume-1)
        double precision mat(0:n*n,0:n*n)
        integer i,j,ti,lr,sigma
 
        if(lr .eq. 1) then
          call multti(mat,lr,sigma)
          do 10 j = 0, n*n-1
            do 10 i = 0, n*n-1
10            mat(i,j) = mat(i,j) / vvv(j+toff*ti)
        else
          do 20 i = 0, n*n-1
            do 20 j = 0, n*n-1
20            mat(i,j) = mat(i,j) / vvv(i+toff*ti)
          call multti(mat,lr,sigma)
        endif
 
        return
        end

c*******************************multt(m,lr,sigma)********************
 
c        This subroutine multiplies m by the matrix T on the
c        left or right depending on whether lr = -1 or 1.
c        This is the double precision version.
cpd00    Put in extra argument sigma to discriminate up/down spin.
cpd00    Variables tt,ek,eki replaced by ttup,ttdn,.....
cpd00    Common block: newke
 
        subroutine multt(mat,lr,sigma)
        implicit none
        include 'paramp.dat'

        double precision mat(0:n*n,0:n*n),temp(0:n*n,0:n*n)
        double precision ttup(0:toff,0:toff),ekup(0:toff,0:toff)
        double precision ekupi(0:toff,0:toff)
        double precision ttdn(0:toff,0:toff),ekdn(0:toff,0:toff)
        double precision ekdni(0:toff,0:toff)
        common/newke/ttup,ekup,ekupi,ttdn,ekdn,ekdni
      save /newke/
        integer i,j,k,lr,sigma

cpd00  Discriminate between spin up/down (if/then/else with doubled loops)
 
	if (sigma.eq.1) then

        if(lr .eq. (-1)) then
	  do 65 j=0,toff-1
	  do 65 i=0,toff-1
	  temp(i,j)=0.d0
	  do 65 k=0,toff-1
	       temp(i,j)=temp(i,j)+ekup(i,k)*mat(k,j)
65	  continue
	  do 68 j=0,toff-1
	  do 68 i=0,toff-1
	       mat(i,j)=temp(i,j)
68	  continue

          else 

	  do 75 j=0,toff-1
	  do 75 i=0,toff-1
	  temp(i,j)=0.d0
	  do 75 k=0,toff-1
	       temp(i,j)=temp(i,j)+mat(i,k)*ekup(k,j)
75	  continue
	  do 78 j=0,toff-1
	  do 78 i=0,toff-1
	       mat(i,j)=temp(i,j)
78	  continue

        endif

	else

        if(lr .eq. (-1)) then
	  do 165 j=0,toff-1
	  do 165 i=0,toff-1
	  temp(i,j)=0.d0
	  do 165 k=0,toff-1
	       temp(i,j)=temp(i,j)+ekdn(i,k)*mat(k,j)
165	  continue
	  do 168 j=0,toff-1
	  do 168 i=0,toff-1
	       mat(i,j)=temp(i,j)
168	  continue

          else 

	  do 175 j=0,toff-1
	  do 175 i=0,toff-1
	  temp(i,j)=0.d0
	  do 175 k=0,toff-1
	       temp(i,j)=temp(i,j)+mat(i,k)*ekdn(k,j)
175	  continue
	  do 178 j=0,toff-1
	  do 178 i=0,toff-1
	       mat(i,j)=temp(i,j)
178	  continue

        endif

	endif
 
        return
        end
c*******************************multti(m,lr,sigma)********************
 
c        This subroutine multiplies m by the matrix T^-1 on the
c        left or right depending on whether lr = -1 or 1.
c        This is the double precision version.
cpd00    Put in extra argument sigma to discriminate up/down spin.
cpd00    Variables tt,ek,eki replaced by ttup,ttdn,..... 
cpd00    Common block: newke
 
        subroutine multti(mat,lr,sigma)
        implicit none
        include 'paramp.dat'

        double precision mat(0:n*n,0:n*n),temp(0:n*n,0:n*n)
        double precision ttup(0:toff,0:toff),ekup(0:toff,0:toff)
        double precision ekupi(0:toff,0:toff)
        double precision ttdn(0:toff,0:toff),ekdn(0:toff,0:toff)
        double precision ekdni(0:toff,0:toff)
        common/newke/ttup,ekup,ekupi,ttdn,ekdn,ekdni
      save /newke/
        integer i,j,k,lr,sigma

cpd00  Discriminate between spin up/down (if/then/else with doubled loops)

	if (sigma.eq.1) then
 
        if(lr .eq. (-1)) then
	  do 65 j=0,toff-1
	  do 65 i=0,toff-1
	  temp(i,j)=0.d0
	  do 65 k=0,toff-1
	       temp(i,j)=temp(i,j)+ekupi(i,k)*mat(k,j)
65	  continue
	  do 68 j=0,toff-1
	  do 68 i=0,toff-1
	       mat(i,j)=temp(i,j)
68	  continue

          else 

	  do 75 j=0,toff-1
	  do 75 i=0,toff-1
	  temp(i,j)=0.d0
	  do 75 k=0,toff-1
	       temp(i,j)=temp(i,j)+mat(i,k)*ekupi(k,j)
75	  continue
	  do 78 j=0,toff-1
	  do 78 i=0,toff-1
	       mat(i,j)=temp(i,j)
78	  continue

        endif

	else

        if(lr .eq. (-1)) then
	  do 165 j=0,toff-1
	  do 165 i=0,toff-1
	  temp(i,j)=0.d0
	  do 165 k=0,toff-1
	       temp(i,j)=temp(i,j)+ekdni(i,k)*mat(k,j)
165	  continue
	  do 168 j=0,toff-1
	  do 168 i=0,toff-1
	       mat(i,j)=temp(i,j)
168	  continue

          else 

	  do 175 j=0,toff-1
	  do 175 i=0,toff-1
	  temp(i,j)=0.d0
	  do 175 k=0,toff-1
	       temp(i,j)=temp(i,j)+mat(i,k)*ekdni(k,j)
175	  continue
	  do 178 j=0,toff-1
	  do 178 i=0,toff-1
	       mat(i,j)=temp(i,j)
178	  continue

        endif

	endif
 
        return
        end

c*******************************onesub********************
 
        subroutine onesub(mat,submat)
        implicit none
	include 'paramp.dat'
 
        integer i,j
        double precision mat(0:n*n,0:n*n),submat(0:n*n,0:n*n)
 
        do 10 i = 0, n*n-1
         do 20 j = 0, n*n-1
20        submat(i,j) = -mat(i,j)
         submat(i,i) = 1.d0 - mat(i,i)
10      continue
        return
        end
 
c********************************
 
       subroutine orthfacp(mat,dvec,rmat,depth)
       implicit none
	include 'paramp.dat'
       integer i,j,depth
       double precision mat(0:n*n,0:n*n),rmat(0:n*n,0:n*n),
     1      dvec(0:n*n-1)
       double precision temp
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
 
c****** NEW VERSION ******* (with pivoting)
       subroutine orthogp(mat,rmat,depth)
       implicit none
	include 'paramp.dat'
 
       integer i,j,k,depth, done(0:n*n-1),ib
       double precision mat(0:n*n,0:n*n),rmat(0:n*n,0:n*n),
     1      maxsize,colsize
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
 
c*****************************ranlat()*****************************
 
c        This subroutine makes a lattice with random values for the
c        Hubbard-Stratanovich variables.
 
        subroutine ranlat()
        implicit none
	include 'paramp.dat'
 
        integer warms,sweeps,nwrap,iran,iran0
        common/integers/warms,sweeps,nwrap,iran,iran0
 
c        common block for vectors
        double precision hub(0:volume-1)
      double precision vup(0:volume-1),vdn(0:volume-1)
      double precision mphase(0:toff-1)
      common/vectors/hub,vup,vdn,mphase
 
        integer i
        double precision ran,ran2
 
        do 10 i=0,volume-1
          hub(i)=1.d0
          ran = ran2(iran)
          if(ran .gt. 0.5d0) hub(i) = -1.d0
10      continue
 
        call setvup()
 
        return
        end
c*******************************************
        subroutine readin()
        implicit none
	include 'paramp.dat'
 
        integer warms,sweeps,nwrap,iran,iran0
        common/integers/warms,sweeps,nwrap,iran,iran0
 
        integer orthlen,doauto
        double precision eorth,difflim
        common/getgpam/eorth,difflim,orthlen,doauto
        save /getgpam/

        double precision sgnup,sgndn,sgn,nup,ndn,nud,ke
        double precision saf,saf2,sferro,sfer2,grfun(0:n/2,0:n/2),
     1       safsq,safsq2
        double precision spinxx(0:n/2,0:n/2),spinzz(0:n/2,0:n/2)
        double precision pairmat(-1:1,-1:1,-1:1,-1:1),errrat,asgnp
        double precision asgnup,asgndn,asgn,anup,andn,anud,
     1       ake,asafsq,asafsq2
        double precision asaf,asaf2,asferro,asfer2,agrfun(0:n/2,0:n/2)
        double precision aspinxx(0:n/2,0:n/2),aspinzz(0:n/2,0:n/2)
        double precision apairmat(-1:1,-1:1,-1:1,-1:1)
        double precision fpairmat(-1:1,-1:1,-1:1,-1:1,0:n/2,0:n/2)
        double precision afpairmat(-1:1,-1:1,-1:1,-1:1,0:n/2,0:n/2)

        double precision den(0:n/2,0:n/2,0:1),aden(0:n/2,0:n/2,0:1)
        integer nmeas0,redo,noredo,dopair,nmeasp,numpair,nmeaspf
        common/m0var/ sgnup,sgndn,sgn,nup,ndn,nud,ke,safsq,safsq2,
     1       saf,saf2,sferro,sfer2,grfun,spinxx,spinzz,aspinxx,
     2       aspinzz,asafsq,asafsq2,asgnup,asgndn,asgn,anup,andn,
     3       anud,ake,den,aden,asaf,asaf2,asferro,asfer2,agrfun,
     4       errrat,asgnp,pairmat,apairmat,
     6       fpairmat,afpairmat
       common/m0vari/dopair,nmeasp,numpair,nmeas0,redo,noredo,nmeaspf
 
         double precision u,mu,bpar,delt,dtau,gam,lambda
         common/couple/u,mu,bpar,delt,dtau,gam,lambda
 
        double precision hub(0:volume-1)
      double precision vup(0:volume-1),vdn(0:volume-1)
      double precision mphase(0:toff-1)
      common/vectors/hub,vup,vdn,mphase

	integer start
	common/startup/start
 
        double precision temp
	integer i
 
        write(6,*) 'enter warms and sweeps'
        read(5,*) warms,sweeps
 237    format(a9,i4)
        write(66,237)' warms  = ',warms
        write(66,237)' sweeps = ',sweeps
        if(tausk.ne.0) then
           write(68,237)' warms  = ',warms
           write(68,237)' sweeps = ',sweeps
        endif

c	HOPPING WILL NOW BE READ IN WHEN ONE BODY MATRIX CONSTRUCTED,
c	SINCE THE NUMBER AND NATURE OF THESE PARAMETERS WILL VARY
c	FROM PROBLEM TO PROBLEM.

        write(6,*)'enter u'
        read(5,*) u
 238    format(a9,f16.12)
        write(66,238)'      u = ',u
        if(tausk.ne.0) then
           write(68,238)'      u = ',u
        endif

        write(6,*) 'enter nwrap, difflim, errrat'
        read(5,*)nwrap,difflim,errrat
        write(66,237)'  nwrap = ',nwrap
        write(66,238)'difflim = ',difflim
        write(66,238)' errrat = ',errrat
        if(tausk.ne.0) then
           write(68,237)'  nwrap = ',nwrap
           write(68,238)'difflim = ',difflim
           write(68,238)' errrat = ',errrat
        endif

        write(6,*) 'enter doauto, orthlen, eorth, dopair, numpair'
        read(5,*)doauto,orthlen,eorth,dopair, numpair
        write(66,237)' doauto = ',doauto
        write(66,237)'orthlen = ',orthlen
        write(66,238)'  eorth = ',eorth
        write(66,237)' dopair = ',dopair
        write(66,237)'numpair = ',numpair
        if(tausk.ne.0) then
           write(68,237)' doauto = ',doauto
           write(68,237)'orthlen = ',orthlen
           write(68,238)'  eorth = ',eorth
           write(68,237)' dopair = ',dopair
           write(68,237)'numpair = ',numpair
        endif

        gam = 0.5d0
 
        temp = dexp(dtau*u*0.5d0)
        lambda = dlog(temp+dsqrt(temp**2-1))
        write(66,*)'lambda is ',lambda
        if(tausk.ne.0) then
           write(68,*)'lambda is ',lambda
        endif
 
        write(6,*)'enter start type'
        read(5,*)start
	if (iran.gt.0) iran=-iran
 239    format(a9,i12)
        write(66,239) '  start = ',start
        if(tausk.ne.0) then
           write(68,239) '  start = ',start
        endif

	if (start.eq.999) then
	     do 871 i=0,volume-1
871	     read(5,*) hub(i)
	endif 

        return
        end
c
c *******************************************************
       subroutine setto0(a,b,c,d,e)
       implicit none
       double precision a,b,c,d,e
       a=0.d0
       b=0.d0
       c=0.d0
       d=0.d0
       e=0.d0
       return
       end
c*****************************setvup()*****************************
 
c        This subroutine sets vup and vdn given hubs.
 
        subroutine setvup()
        implicit none
	include 'paramp.dat'
 
         double precision u,mu,bpar,delt,dtau,gam,lambda
         common/couple/u,mu,bpar,delt,dtau,gam,lambda
 
c        common block for vectors
        double precision hub(0:volume-1)
      double precision vup(0:volume-1),vdn(0:volume-1)
      double precision mphase(0:toff-1)
      common/vectors/hub,vup,vdn,mphase
 
        integer i
 
        do 20 i=0,volume-1
          vup(i)=dexp(lambda*hub(i))
          vdn(i)=dexp(-lambda*hub(i))
20      continue
 
        return
        end
 
c****************************siteindx**********************************
 
c        this function finds the index of a site given its three coordinates
 
        integer function siteindx(x,y,ti)
 
	include 'paramp.dat'
        integer x,y,ti
 
        siteindx = x + n*(y + n*ti)
 
        return
        end
c ************************************************************************
       SUBROUTINE SPLINE(X,Y,N,YP1,YPN,Y2)
       implicit none
       integer nmax,n
       PARAMETER (NMAX=400)
       double precision x(n),y(n),y2(n),u(nmax)
       double precision yp1,ypn,p,qn,sig,un
       integer i,k
       IF (YP1.GT..99D30) THEN
         Y2(1)=0.d0
         U(1)=0.d0
       ELSE
         Y2(1)=-0.5d0
         U(1)=(3.d0/(X(2)-X(1)))*((Y(2)-Y(1))/(X(2)-X(1))-YP1)
       ENDIF
       DO 10 I=2,N-1
         SIG=(X(I)-X(I-1))/(X(I+1)-X(I-1))
         P=SIG*Y2(I-1)+2.d0
         Y2(I)=(SIG-1.d0)/P
         U(I)=(6.d0*((Y(I+1)-Y(I))/(X(I+1)-X(I))-(Y(I)-Y(I-1))
     *         /(X(I)-X(I-1)))/(X(I+1)-X(I-1))-SIG*U(I-1))/P
10      continue
       IF (YPN.GT..99D30) THEN
         QN=0.d0
         UN=0.d0
       ELSE
         QN=0.5d0
         UN=(3.d0/(X(N)-X(N-1)))*(YPN-(Y(N)-Y(N-1))/(X(N)-X(N-1)))
       ENDIF
       Y2(N)=(UN-QN*U(N-1))/(QN*Y2(N-1)+1.d0)
       DO 20 K=N-1,1,-1
              Y2(K)=Y2(K)*Y2(K+1)+U(K)
20     continue
       RETURN
       END
c
c ********************************************************************
      SUBROUTINE SPLINT(XA,YA,Y2A,N,X,Y)
      implicit none
      integer n
      double precision XA(N),YA(N),Y2A(N)
      integer klo,khi,k
      double precision x,y,h,a,b
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
       IF (H.EQ.0.d0) then
         write(6,*) 'Bad XA input.'
         stop
       endif
       A=(XA(KHI)-X)/H
       B=(X-XA(KLO))/H
       Y=A*YA(KLO)+B*YA(KHI)+
     1     ((A**3-A)*Y2A(KLO)+(B**3-B)*Y2A(KHI))*(H**2)/6.d0
       RETURN
       END
c*******************************************
        subroutine sweep(gmatup,gmatdn,accept,reject,wraps)
        implicit none
	include 'paramp.dat'
 
        integer warms,sweeps,nwrap,iran,iran0
        common/integers/warms,sweeps,nwrap,iran,iran0
 
        integer orthlen,doauto
        double precision eorth,difflim
        common/getgpam/eorth,difflim,orthlen,doauto
        save /getgpam/

         double precision u,mu,bpar,delt,dtau,gam,lambda
         common/couple/u,mu,bpar,delt,dtau,gam,lambda
 
c        common block for vectors
        double precision hub(0:volume-1)
      double precision vup(0:volume-1),vdn(0:volume-1)
      double precision mphase(0:toff-1)
      common/vectors/hub,vup,vdn,mphase
 
        double precision sgnup,sgndn,sgn,nup,ndn,nud,ke
        double precision saf,saf2,sferro,sfer2,grfun(0:n/2,0:n/2),
     1       safsq,safsq2
        double precision spinxx(0:n/2,0:n/2),spinzz(0:n/2,0:n/2)
        double precision pairmat(-1:1,-1:1,-1:1,-1:1),errrat,asgnp
        double precision asgnup,asgndn,asgn,anup,andn,anud,ake,
     1       asafsq,asafsq2
        double precision asaf,asaf2,asferro,asfer2,agrfun(0:n/2,0:n/2)
        double precision aspinxx(0:n/2,0:n/2),aspinzz(0:n/2,0:n/2)
        double precision apairmat(-1:1,-1:1,-1:1,-1:1)
        double precision fpairmat(-1:1,-1:1,-1:1,-1:1,0:n/2,0:n/2)
        double precision afpairmat(-1:1,-1:1,-1:1,-1:1,0:n/2,0:n/2)

        double precision den(0:n/2,0:n/2,0:1),aden(0:n/2,0:n/2,0:1)
        integer nmeas0,redo,noredo,dopair,nmeasp,numpair,nmeaspf
        common/m0var/ sgnup,sgndn,sgn,nup,ndn,nud,ke,safsq,safsq2,
     1       saf,saf2,sferro,sfer2,grfun,spinxx,spinzz,aspinxx,
     2       aspinzz,asafsq,asafsq2,asgnup,asgndn,asgn,anup,andn,
     3       anud,ake,den,aden,asaf,asaf2,asferro,asfer2,agrfun,
     4       errrat,asgnp,pairmat,apairmat,
     6       fpairmat,afpairmat
       common/m0vari/dopair,nmeasp,numpair,nmeas0,redo,noredo,nmeaspf
 
        integer accept,reject,wraps
        integer ti,ti1
        double precision gmatup(0:n*n,0:n*n),ogmat(0:n*n,0:n*n),diffup
        double precision gmatdn(0:n*n,0:n*n),diffdn
        double precision accrat,redorat
 
        integer lastwr,skipair

c     initialize lastwr to 0
        lastwr = 0

        do 10 ti = 0, l-1
          call multb(gmatup,vup,ti,-1,1)
          call multbi(gmatup,vup,ti,1,1)
c          if(mu .ne. 0.d0.or.bpar.ne.0.d0)then
            call multb(gmatdn,vdn,ti,-1,-1)
            call multbi(gmatdn,vdn,ti,1,-1)
c          endif
          if(wraps .gt. 0) then
            wraps = wraps - 1
          else
            wraps = nwrap
            ti1 = mod(ti+1,l)
            call matcop(gmatup,ogmat)
            call getgp(vup,ti1,gmatup,sgnup,1)
            call matdif(gmatup,ogmat,diffup)
c            if(mu .ne. 0.d0.or.bpar.ne.0.d0) then
              call matcop(gmatdn,ogmat)
              call getgp(vdn,ti1,gmatdn,sgndn,-1)
              call matdif(gmatdn,ogmat,diffdn)
c            else
c              diffdn = diffup
c              sgndn = sgnup
c            endif
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
          if(mod(ti,12) .eq. 0)call meas0sq(gmatup,gmatdn)
c Measure pair correlations every skipair slices.
c Set skipair so numpair meas. are done every sweep.
          skipair = l/numpair
          if(mod(l,numpair) .ne. 0)skipair = skipair+1
          if(mod(ti,skipair).eq.0)call measpairsq(gmatup,gmatdn)
          if(mod(ti,skipair).eq.0)call pairfoursq(gmatup,gmatdn)
10      continue
 
        if(accept+reject .gt. 10000) then
          accrat = dble(accept)/dble(accept+reject)
          if(accrat .gt. 0.52d0 .or. accrat .lt. 0.48d0)then
            gam = gam + (accrat - 0.5d0)
            gam = dmax1(zero,gam)
            gam = dmin1(one,gam)
            accept = int(100*accrat)
            reject = int(100*(1.d0-accrat))
          endif
        endif
 
        if(redo .gt. 20) then
          redorat = dble(redo)/dble(redo+noredo)
          if(redorat.gt. errrat)then
            nwrap = nwrap - 1
            write(66,*)'reducing nwrap to ',nwrap,' #NO_AVE'
            redo = 0
            noredo = 1
          endif
        endif
 
        if(noredo .gt. 500) then
          redorat = dble(redo)/dble(redo+noredo)
          if(redorat .lt. 0.2d0*errrat) then
            if(nwrap .ge. lastwr)then
c Only increase wrap if nwrap has not been reduced since last increase.
              nwrap = nwrap + 2
              write(66,*)'increasing nwrap to ',nwrap,' #NO_AVE'
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
	include 'paramp.dat'
 
        integer warms,sweeps,nwrap,iran,iran0
        common/integers/warms,sweeps,nwrap,iran,iran0
 
         double precision u,mu,bpar,delt,dtau,gam,lambda
         common/couple/u,mu,bpar,delt,dtau,gam,lambda
 
c        common block for vectors
        double precision hub(0:volume-1)
      double precision vup(0:volume-1),vdn(0:volume-1)
      double precision mphase(0:toff-1)
      common/vectors/hub,vup,vdn,mphase
 
        integer accept,reject
        integer ti,nbar,nbari,j1,j2
        double precision gmatup(0:n*n,0:n*n),sgnup
        double precision gmatdn(0:n*n,0:n*n),sgndn
        double precision bvecup(0:n*n-1),bvecup2(0:n*n-1)
        double precision bvecdn(0:n*n-1),bvecdn2(0:n*n-1),ratdn,vnbardn
        double precision vnbarup,ratup,ran,rat,p
        double precision ran2
 
        do 20 nbar = 0, toff-1
          nbari = nbar + ti * toff
          vnbarup = dexp(-2.d0*lambda*hub(nbari))
          vnbardn = 1.d0 / vnbarup - 1.d0
          vnbarup = vnbarup - 1.d0
          ratup = 1.d0 + (1.d0 - gmatup(nbar,nbar))*vnbarup
c          if(mu .ne. 0.d0.or.bpar.ne.0.d0) then
            ratdn = 1.d0 + (1.d0 - gmatdn(nbar,nbar))*vnbardn
c          else
c            ratdn = 1.d0 + gmatup(nbar,nbar)*vnbardn
c          endif
          ran = ran2(iran)
          rat = dabs(ratup*ratdn)
          if(rat .le. 1.d0) p = rat/(1.d0+gam*rat)
          if(rat .gt. 1.d0) p = rat/(gam+rat)
          if(p .gt. ran) then
            accept = accept + 1
            if(ratup .lt. 0.d0)sgnup = -sgnup
            if(ratdn .lt. 0.d0)sgndn = -sgndn
c            if(mu .ne. 0.d0.or.bpar.ne.0.d0) then
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
c            else
c              do 35 j1 = 0, toff-1
c35              bvecup(j1) = (-vnbarup)*gmatup(nbar,j1)
c              bvecup(nbar) = bvecup(nbar) + vnbarup
c              do 45 j1 = 0, toff-1
c45              bvecup2(j1) = gmatup(j1,nbar)/(1.d0+bvecup(nbar))
c              do 55 j1 = 0, toff-1
c               do 55 j2 = 0, toff-1
c55              gmatup(j1,j2)=gmatup(j1,j2)-bvecup2(j1)*bvecup(j2)
c            endif
            hub(nbari) = -hub(nbari)
            vup(nbari) = vup(nbari) * (vnbarup+1.d0)
            vdn(nbari) = vdn(nbari) * (vnbardn+1.d0)
          else
            reject = reject + 1
          endif
20      continue
 
        return
        end
c       ******************unit(mat)********************
c        This subroutine puts the unit matrix in aux, which is real.
 
        subroutine unit(mat)
        implicit none
	include 'paramp.dat'
        double precision mat(0:n*n,0:n*n)
        integer i,j
        do 10 i = 0, n*n-1
          do 20 j = 0, n*n-1
20          mat(i,j) = 0.d0
10        mat(i,i) = 1.d0
        return
        end
c
c***************************************
        subroutine zeroas()
        implicit none
	include 'paramp.dat'
 
        double precision sgnup,sgndn,sgn,nup,ndn,nud,ke
        double precision saf,saf2,sferro,sfer2,grfun(0:n/2,0:n/2),
     1       safsq,safsq2
        double precision spinxx(0:n/2,0:n/2),spinzz(0:n/2,0:n/2)
        double precision pairmat(-1:1,-1:1,-1:1,-1:1),errrat,asgnp
        double precision asgnup,asgndn,asgn,anup,andn,anud,ake,
     1       asafsq,asafsq2
        double precision asaf,asaf2,asferro,asfer2,agrfun(0:n/2,0:n/2)
        double precision aspinxx(0:n/2,0:n/2),aspinzz(0:n/2,0:n/2)
        double precision apairmat(-1:1,-1:1,-1:1,-1:1)
        double precision afpairmat(-1:1,-1:1,-1:1,-1:1,0:n/2,0:n/2)
        double precision fpairmat(-1:1,-1:1,-1:1,-1:1,0:n/2,0:n/2)

        double precision den(0:n/2,0:n/2,0:1),aden(0:n/2,0:n/2,0:1)
        integer nmeas0,redo,noredo,dopair,nmeasp,numpair,nmeaspf

        common/m0var/ sgnup,sgndn,sgn,nup,ndn,nud,ke,safsq,safsq2,
     1       saf,saf2,sferro,sfer2,grfun,spinxx,spinzz,aspinxx,
     2       aspinzz,asafsq,asafsq2,asgnup,asgndn,asgn,anup,andn,
     3       anud,ake,den,aden,asaf,asaf2,asferro,asfer2,agrfun,
     4       errrat,asgnp,pairmat,apairmat,
     6       fpairmat,afpairmat
       common/m0vari/dopair,nmeasp,numpair,nmeas0,redo,noredo,nmeaspf
 
        double precision gnl(0:n/2,0:n/2,0:l),agnl(0:n/2,0:n/2,0:l)
        double precision chinl(0:n/2,0:n/2,0:l),achinl(0:n/2,0:n/2,0:l)
        double precision dent(0:n/2,0:n/2,0:l),adent(0:n/2,0:n/2,0:l)
        double precision chinlz(0:n/2,0:n/2,0:l),
     1       achinlz(0:n/2,0:n/2,0:l)
        double precision cnl(0:n/2,0:n/2,0:l),acnl(0:n/2,0:n/2,0:l)
        double precision pstau(0:l),pdtau(0:l),psxtau(0:l)
        double complex lam(0:n*n,0:n*n),alam(0:n*n,0:n*n)
	double complex ag1w(0:n*n,0:n*n)
        double precision apstau(0:l),apdtau(0:l),apsxtau(0:l),asgnt
        integer nmeast
        common/mtauvar/gnl,agnl,chinl,achinl,chinlz,achinlz,cnl,
     1   acnl,pstau,pdtau,psxtau,apstau,apdtau,apsxtau,lam,alam,
     2	 ag1w,asgnt,dent,adent
        common/mtauvari/nmeast
	double precision ak2
 
        integer i,j,i2,j2,ti,kx,ky
 
        nmeas0 = 0
        nmeasp = 0
        nmeaspf = 0
        nmeast = 0
        call setto0(asgnup,asgndn,asgn,anup,andn)
        call setto0(anud,ake,asaf,asaf2,asferro)
        call setto0(asfer2,asafsq,asafsq2,asgnp,asgnt)
        call setto0(ak2,asgndn,asgn,anup,andn)
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
              adent(i,j,ti) = 0.d0
              achinlz(i,j,ti) = 0.d0
              acnl(i,j,ti) = 0.d0
10      continue
        do 40 i = 0, n*n
          do 40 j = 0, n*n
	    alam(i,j) = 0.d0
 40      ag1w(i,j) = 0.d0
        if(dopair .eq. 1)then
          do 20 i = -1,1
            do 20 j = -1,1
              do 20 i2 = -1,1
                do 20 j2 = -1,1
                  apairmat(i,j,i2,j2) = 0.d0
20        continue
        do 30 ti = 0, l
          apstau(ti) = 0.d0
          apdtau(ti) = 0.d0
30        apsxtau(ti) = 0.d0
      do 21 kx=0,n/2
      do 21 ky=0,n/2
          do 21 i = -1,1
            do 21 j = -1,1
              do 21 i2 = -1,1
                do 21 j2 = -1,1
21                afpairmat(i,j,i2,j2,kx,ky) = 0.d0

        endif
         

        return
        end
c       ******************zeromat(mat)********************
c        This subroutine zeroes mat.
 
        subroutine zeromat(mat)
        implicit none
	include 'paramp.dat'
        double precision mat(0:n*n,0:n*n)
        integer i,j
        do 20 i = 0, n*n-1
          do 20 j = 0, n*n-1
20          mat(i,j) = 0.d0
        return
        end
c*******************************getdiag********************
 
        subroutine getdiag(bigmat,diagg,vvv,sigma)
        implicit none
	include 'paramp.dat'
        integer llim(0:ndiv)
        double precision vvv(0:volume-1)
 
        integer ti,i,j,nn,kp,sigma
        double precision bigmat(0:toff*ndiv-1,0:toff*ndiv-1)
        double precision diagg(0:n*n,0:n*n,0:l)

c     initialization flags
      integer ifbigt,ifftnt,ifgetc,ifgetr,ifmeas,ifgetd
      common/iflags/ifbigt,ifftnt,ifgetc,ifgetr,ifmeas,ifgetd

      save llim

      if(ifgetd .ne. 373) then
         ifgetd = 373
         do 1 kp = 0,ndiv
 1          llim(kp) = (l * kp)/ ndiv - 1
      endif

        nn = toff
        do 10 kp = 0, ndiv-1
         ti = llim(kp)+1
         do 10 i = 0, n*n-1
          do 10 j = 0, n*n-1
           diagg(i,j,ti) = bigmat(kp*nn+i,kp*nn+j)
10      continue
        call matcop(diagg(0,0,0),diagg(0,0,l))
        do 20 kp = 0, ndiv-1
         do 20 ti = llim(kp)+2,(llim(kp+1)+llim(kp)+1)/2
          call matcop(diagg(0,0,ti-1),diagg(0,0,ti))
          call multb(diagg(0,0,ti),vvv,ti-1,-1,sigma)
          call multbi(diagg(0,0,ti),vvv,ti-1,1,sigma)
20      continue
        do 30 kp = 1, ndiv
         do 30 ti = llim(kp),(llim(kp)+llim(kp-1)+3)/2,-1
          call matcop(diagg(0,0,ti+1),diagg(0,0,ti))
          call multbi(diagg(0,0,ti),vvv,ti,-1,sigma)
          call multb(diagg(0,0,ti),vvv,ti,1,sigma)
30      continue
        return
        end

c ************************************************************

      DOUBLE PRECISION FUNCTION RAN2(IDUM)
      implicit none
      double precision rm
      integer ia,ic,idum,iff,ir,iy,j,m
      save
      PARAMETER (M=714025,IA=1366,IC=150889,RM=1.4005112d-6)
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


c     THIS SUBROUTINE EXPONENTIATES THE ONE BODY PART OF HAMILTONIAN,
c     IE HOPPING AND SITE ENERGIES.

      subroutine getek
      implicit none
      include 'paramp.dat'
      integer i,j,k,nrot
      double precision ev(toff),evect(toff,toff)
      double precision ekt(toff,toff),temp(toff,toff)
      double precision ttt(toff,toff)
        double precision ttup(0:toff,0:toff),ekup(0:toff,0:toff)
        double precision ekupi(0:toff,0:toff)
        double precision ttdn(0:toff,0:toff),ekdn(0:toff,0:toff)
        double precision ekdni(0:toff,0:toff)
        common/newke/ttup,ekup,ekupi,ttdn,ekdn,ekdni
      save /newke/
      double precision u,mu,bpar,delt,dtau,gam,lambda
      common/couple/u,mu,bpar,delt,dtau,gam,lambda

c     GETS THE MATRIX TO BE EXPONENTIATED
c     DIFFERENT VERSIONS FOR DIFFERENT GEOMETRIES.

      call gett1band2dsq

c     DIAGONALIZE ttup.  (SAVE A COPY.  SHIFT INDICES.  PUT IN dtau.).

	do 110 i=0,toff-1
	do 100 j=0,toff-1
	     ttt(i+1,j+1)=-dtau*ttup(i,j)
100	continue
110	continue

       call jacobi(ttt,toff,toff,ev,evect,nrot)

c     NOW EXPONENTIATE tt.

       do 260 i=1,toff
       do 240 j=1,toff
             temp(i,j)=exp(ev(i))*evect(j,i)
240     continue
260     continue
        do 290 i=1,toff
	do 280 j=1,toff
	ekt(i,j)=0.d0
        do 270 k=1,toff
	     ekt(i,j)=ekt(i,j)+evect(i,k)*temp(k,j)
270     continue
280     continue
290     continue

c	FINALLY WRITE THINGS IN 0,toff-1 CONVENTION

        do 297 i=1,toff
	do 293 j=1,toff
	     ekup(i-1,j-1)=ekt(i,j)
293     continue
297     continue

c     DO INVERSE

       do 1260 i=1,toff
       do 1240 j=1,toff
             temp(i,j)=exp(-ev(i))*evect(j,i)
1240     continue
1260     continue
        do 1290 i=1,toff
	do 1280 j=1,toff
	ekt(i,j)=0.d0
        do 1270 k=1,toff
	     ekt(i,j)=ekt(i,j)+evect(i,k)*temp(k,j)
1270     continue
1280     continue
1290     continue

c	FINALLY WRITE THINGS IN 0,toff-1 CONVENTION

        do 1297 i=1,toff
	do 1293 j=1,toff
	     ekupi(i-1,j-1)=ekt(i,j)
1293     continue
1297     continue

c     DIAGONALIZE ttdn.  (SAVE A COPY.  SHIFT INDICES.  PUT IN dtau.).

	do 5110 i=0,toff-1
	do 5100 j=0,toff-1
	     ttt(i+1,j+1)=-dtau*ttdn(i,j)
5100	continue
5110	continue

       call jacobi(ttt,toff,toff,ev,evect,nrot)

c     NOW EXPONENTIATE tt.

       do 5260 i=1,toff
       do 5240 j=1,toff
             temp(i,j)=exp(ev(i))*evect(j,i)
5240     continue
5260     continue
        do 5290 i=1,toff
	do 5280 j=1,toff
	ekt(i,j)=0.d0
        do 5270 k=1,toff
	     ekt(i,j)=ekt(i,j)+evect(i,k)*temp(k,j)
5270     continue
5280     continue
5290     continue

c	FINALLY WRITE THINGS IN 0,toff-1 CONVENTION

        do 5297 i=1,toff
	do 5293 j=1,toff
	     ekdn(i-1,j-1)=ekt(i,j)
5293     continue
5297     continue

c     DO INVERSE

       do 6260 i=1,toff
       do 6240 j=1,toff
             temp(i,j)=exp(-ev(i))*evect(j,i)
6240     continue
6260     continue
        do 6290 i=1,toff
	do 6280 j=1,toff
	ekt(i,j)=0.d0
        do 6270 k=1,toff
	     ekt(i,j)=ekt(i,j)+evect(i,k)*temp(k,j)
6270     continue
6280     continue
6290     continue

c	FINALLY WRITE THINGS IN 0,toff-1 CONVENTION

        do 6297 i=1,toff
	do 6293 j=1,toff
	     ekdni(i-1,j-1)=ekt(i,j)
6293     continue
6297     continue

      return
      end
 
 
c %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c      following code until next row of percentage signs is all 
c      2d square lattice stuff.
c %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

c	THIS SUBROUTINE CONSTRUCTS THE ONE BODY MATRIX TO BE
c	DIAGONALIZED.  2-D SQUARE LATTICE.

      subroutine gett1band2dsq
      implicit none
      include 'paramp.dat'
      integer i,j
      double precision t,mu,delt,bpar,siteeav,sitee(0:toff),ran2
         double precision u,dtau,gam,lambda
         common/couple/u,mu,bpar,delt,dtau,gam,lambda
c        index vectors for indirect addressing and phase vector
        integer xplus(0:volume-1),xminus(0:volume-1)
        integer yplus(0:volume-1),yminus(0:volume-1)
        common/ivectors/xplus,xminus,yplus,yminus
        save /ivectors/
        double precision ttup(0:toff,0:toff),ekup(0:toff,0:toff)
        double precision ekupi(0:toff,0:toff)
        double precision ttdn(0:toff,0:toff),ekdn(0:toff,0:toff)
        double precision ekdni(0:toff,0:toff)
        common/newke/ttup,ekup,ekupi,ttdn,ekdn,ekdni
      save /newke/
      integer warms,sweeps,nwrap,iran,iran0
      common/integers/warms,sweeps,nwrap,iran,iran0

c	input stuff
237     format(a9,i4)
257     format(a9,f12.8)
        write(6,*)'enter random number seed'
        read(5,*)iran
        iran0=iran
	if (iran.gt.0) iran=-iran
239     format(a9,i12)
        write(66,239) '   iran = ',iran
        if(tausk.ne.0) then
           write(68,239) '   iran = ',iran
        endif

c     set up the index vectors
      call indexsetsq()

c     set up the phase vectors
      call phaseset()
	
        write(6,*) 'enter t,delt,mu,bpar,dtau'
        read (5,*)        t,delt,mu,bpar,dtau
        write(66,257)'     t  = ',t
        write(66,257)'  delt  = ',delt
        write(66,257)'    mu  = ',mu
        write(66,257)'  bpar  = ',bpar
        write(66,257)'  dtau  = ',dtau 
        if(tausk.ne.0) then
             write(68,257)'     t  = ',t
             write(68,257)'  delt  = ',delt
             write(68,257)'    mu  = ',mu
             write(68,257)'  bpar  = ',bpar
             write(68,257)'  dtau  = ',dtau 
        endif

	write (66,*) 'x bond strengths'
	do 100 i=0,toff-1
	     j=xplus(i)
	     ttup(i,j)=-t+delt*(ran2(iran)-0.5d0)
	     ttup(j,i)=ttup(i,j)
	     ttdn(i,j)=ttup(i,j)
	     ttdn(j,i)=ttdn(i,j)
	     write (66,248) i,j,ttup(i,j),ttdn(i,j)
100	continue
	write (66,*) 'y bond strengths'
	do 110 i=0,toff-1
	     j=yplus(i)
	     ttup(i,j)=-t+delt*(ran2(iran)-0.5d0)
	     ttup(j,i)=ttup(i,j)
	     ttdn(i,j)=ttup(i,j)
	     ttdn(j,i)=ttdn(i,j)
	     write (66,248) i,j,ttup(i,j),ttdn(i,j)
110	continue

cpd  Following part changed/commented out to include parallel (Zeeman)
cpd  magnetic field. (bpar replaces old variable delmu)   (PJHD, June 2001)
cpd  Note: This code cannot handle random site energies anymore because
cpd  of these changes!
c	siteeav=0.d0
c	do 200 i=0,toff-1
c	     sitee(i)=delmu*(ran2(iran)-0.5d0)
c	     siteeav=siteeav+sitee(i)
c200	continue
	write (66,*) 'site energies (including mu)'
	do 210 i=0,toff-1
c	     sitee(i)=sitee(i)-siteeav/dfloat(toff)
	     ttup(i,i)=-mu+bpar
	     ttdn(i,i)=-mu-bpar
	     write (66,247) i,ttup(i,i),ttdn(i,i)
210	continue
247     format(i6,2f12.6)
248     format(2i6,2f12.6)

      return
      end


c**************************indexsetsq()************************
c        this subroutine sets and loads the index vectors
        subroutine indexsetsq()
        implicit none
	include 'paramp.dat'
 
c        index vectors for indirect addressing and phase vector
        integer xplus(0:volume-1),xminus(0:volume-1)
        integer yplus(0:volume-1),yminus(0:volume-1)
        common/ivectors/xplus,xminus,yplus,yminus
        save /ivectors/
 
        integer i,neighborsq

c        calculate the index vectors
        do 10 i=0,volume-1
            xplus(i) =neighborsq(i,1,0,0)
            xminus(i)=neighborsq(i,-1,0,0)
            yplus(i) =neighborsq(i,0,1,0)
            yminus(i)=neighborsq(i,0,-1,0)
10        continue
 
        return
        end
c****************************neighborsq()*****************************
 
c        this function finds the index of a neighboring site
 
        integer function neighborsq(site,delx,dely,delt)
 
	include 'paramp.dat'
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
        neighborsq = siteindx(x,y,ti)
 
        return
        end
c********* phaseset() - set and load phase vectors **********
        subroutine phaseset()
        implicit none
	include 'paramp.dat'
 
c        common block for vectors
        double precision hub(0:volume-1)
      double precision vup(0:volume-1),vdn(0:volume-1)
      double precision mphase(0:toff-1)
      common/vectors/hub,vup,vdn,mphase

        integer ix,iy,addr
        double precision ph
 
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
 
c*******************************meas0sq(lots of variables)********************
 
c        This subroutine does the measurements.
 
        subroutine meas0sq(gmatup,gmatdn)
        implicit none
	include 'paramp.dat'
 
         double precision u,mu,bpar,delt,dtau,gam,lambda
         common/couple/u,mu,bpar,delt,dtau,gam,lambda
 
c        common block for vectors
        double precision hub(0:volume-1)
      double precision vup(0:volume-1),vdn(0:volume-1)
      double precision mphase(0:toff-1)
      common/vectors/hub,vup,vdn,mphase

        double precision sgnup,sgndn,sgn,nup,ndn,nud,ke
        double precision saf,saf2,sferro,sfer2,grfun(0:n/2,0:n/2),
     1       safsq,safsq2
        double precision spinxx(0:n/2,0:n/2),spinzz(0:n/2,0:n/2)
        double precision pairmat(-1:1,-1:1,-1:1,-1:1),errrat,asgnp
        double precision fpairmat(-1:1,-1:1,-1:1,-1:1,0:n/2,0:n/2)
        double precision afpairmat(-1:1,-1:1,-1:1,-1:1,0:n/2,0:n/2)

        double precision asgnup,asgndn,asgn,anup,andn,anud,ake,
     1       asafsq,asafsq2
        double precision asaf,asaf2,asferro,asfer2,agrfun(0:n/2,0:n/2)
        double precision aspinxx(0:n/2,0:n/2),aspinzz(0:n/2,0:n/2)
        double precision apairmat(-1:1,-1:1,-1:1,-1:1)
        double precision den(0:n/2,0:n/2,0:1),aden(0:n/2,0:n/2,0:1)
        integer nmeas0,redo,noredo,dopair,nmeasp,numpair,nmeaspf
        common/m0var/ sgnup,sgndn,sgn,nup,ndn,nud,ke,safsq,safsq2,
     1       saf,saf2,sferro,sfer2,grfun,spinxx,spinzz,aspinxx,
     2       aspinzz,asafsq,asafsq2,asgnup,asgndn,asgn,anup,andn,
     3       anud,ake,den,aden,asaf,asaf2,asferro,asfer2,agrfun,
     4       errrat,asgnp,pairmat,apairmat,
     6       fpairmat,afpairmat
       common/m0vari/dopair,nmeasp,numpair,nmeas0,redo,noredo,nmeaspf
 
        double precision gmatup(0:n*n,0:n*n),gmatdn(0:n*n,0:n*n)
        double precision gtemp,xxtemp,zztemp,xfac,yfac

c     occupancy histogram variables (uses parameter histno)
      integer hist1(histno*3+2),hist2(histno*3+2),hist3(histno*3+2)
      common/histogram/hist1,hist2,hist3
      integer histi

      integer i,j,ix,jx,iy,jy,kx,ky

        double precision ttup(0:toff,0:toff),ekup(0:toff,0:toff)
        double precision ekupi(0:toff,0:toff)
        double precision ttdn(0:toff,0:toff),ekdn(0:toff,0:toff)
        double precision ekdni(0:toff,0:toff)
        common/newke/ttup,ekup,ekupi,ttdn,ekdn,ekdni
      save /newke/
 
        saf = 0.d0
        safsq = 0.d0
        sferro = 0.d0
        sfer2 = 0.d0
        saf2 = 0.d0
        safsq2 = 0.d0
        nup = 0.d0
        ndn = 0.d0
        ke = 0.d0
        nud = 0.d0
c         if(mu .eq. 0.d0.and.bpar.eq.0.d0) then
c          do 10 i = 0,n*n-1
c           do 20 j = 0,n*n-1
c20           gmatdn(i,j) = -mphase(i)*mphase(j)*gmatup(j,i)
c10        gmatdn(i,i) = gmatdn(i,i) + 1.d0
c         endif
 
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
c     make the histograms
c     single occupancy,nup
          histi = histno + nint(dble(histno-1) * 
     1         (1.d0 - gmatup(i,i))) + 1
          if(histi.gt.histno*3) then
             hist1(histno*3+2) = hist1(histno*3+2) + 1
          else
             if(histi.le.0) then
                hist1(histno*3+1) = hist1(histno*3+1) + 1
             else
                hist1(histi) = hist1(histi) + 1
             endif
          endif
c     single occupancy,ndn
          histi = histno + nint(dble(histno-1) * 
     1         (1.d0 - gmatdn(i,i))) + 1
          if(histi.gt.histno*3) then
             hist2(histno*3+2) = hist2(histno*3+2) + 1
          else
             if(histi.le.0) then
                hist2(histno*3+1) = hist2(histno*3+1) + 1
             else
                hist2(histi) = hist2(histi) + 1
             endif
          endif
c     double occupancy
          histi = histno + nint(dble(histno-1) * 
     1         (1.d0 - gmatup(i,i)) * (1.d0 - gmatdn(i,i))) + 1
          if(histi.gt.histno*3) then
             hist3(histno*3+2) = hist3(histno*3+2) + 1
          else
             if(histi.le.0) then
                hist3(histno*3+1) = hist3(histno*3+1) + 1
             else
                hist3(histi) = hist3(histi) + 1
             endif
          endif

c	compute ke by multiplying green's function by tt.

        do 198 j=0,toff-1
                if (i.eq.j) then
		   ke=ke+(1.d0-gmatup(i,j))*ttup(i,j)
     1   	        +(1.d0-gmatdn(i,j))*ttdn(i,j)
                   else
		   ke=ke-gmatup(i,j)*ttup(i,j)
     1   	        -gmatdn(i,j)*ttdn(i,j)
                endif
198	continue

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
            kx = iabs(ix-jx)
            kx = min(kx,n-kx)
            ky = iabs(iy-jy)
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
            xxtemp = (-2.d0)*gmatup(i,j)*gmatdn(j,i)
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
        ke = ke / toff
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
        den(0,0,0) = (nup + ndn) * 0.5d0
 
        do 31 kx = 0, n/2
          do 31 ky = 0, n/2
            grfun(kx,ky) = (grfun(kx,ky)+grfun(ky,kx))*0.5d0
            grfun(ky,kx) = grfun(kx,ky)
            den(kx,ky,0) = (den(kx,ky,0)+den(ky,kx,0))*0.5d0
            den(ky,kx,0) = den(kx,ky,0)
            den(kx,ky,1) = (den(kx,ky,1)+den(ky,kx,1))*0.5d0
            den(ky,kx,1) = den(kx,ky,1)
            spinxx(kx,ky) = (spinxx(kx,ky) + spinxx(ky,kx))*0.5d0
            spinxx(ky,kx) = spinxx(kx,ky)
            spinzz(kx,ky) = (spinzz(kx,ky) + spinzz(ky,kx))*0.5d0
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
c*******************************measpairsq(gmatup)********************
 
c        This subroutine does the measurements.
 
        subroutine measpairsq(gmatup,gmatdn)
        implicit none
	include 'paramp.dat'
 
         double precision u,mu,bpar,delt,dtau,gam,lambda
         common/couple/u,mu,bpar,delt,dtau,gam,lambda
 
c        index vectors for indirect addressing and phase vector
        integer xplus(0:volume-1),xminus(0:volume-1)
        integer yplus(0:volume-1),yminus(0:volume-1)
 
        common/ivectors/xplus,xminus,yplus,yminus
        save /ivectors/

c        common block for vectors
        double precision hub(0:volume-1)
        double precision vup(0:volume-1),vdn(0:volume-1)
        double precision mphase(0:toff-1)
        common/vectors/hub,vup,vdn,mphase

        double precision sgnup,sgndn,sgn,nup,ndn,nud,ke
        double precision saf,saf2,sferro,sfer2,grfun(0:n/2,0:n/2),
     1       safsq,safsq2
        double precision spinxx(0:n/2,0:n/2),spinzz(0:n/2,0:n/2)
        double precision pairmat(-1:1,-1:1,-1:1,-1:1),errrat,asgnp
        double precision asgnup,asgndn,asgn,anup,andn,anud,ake,
     1       asafsq,asafsq2
        double precision asaf,asaf2,asferro,asfer2,agrfun(0:n/2,0:n/2)
        double precision aspinxx(0:n/2,0:n/2),aspinzz(0:n/2,0:n/2)
        double precision apairmat(-1:1,-1:1,-1:1,-1:1)
        double precision fpairmat(-1:1,-1:1,-1:1,-1:1,0:n/2,0:n/2)
        double precision afpairmat(-1:1,-1:1,-1:1,-1:1,0:n/2,0:n/2)

        double precision den(0:n/2,0:n/2,0:1),aden(0:n/2,0:n/2,0:1)
        integer nmeas0,redo,noredo,dopair,nmeasp,numpair,nmeaspf
        common/m0var/ sgnup,sgndn,sgn,nup,ndn,nud,ke,safsq,safsq2,
     1       saf,saf2,sferro,sfer2,grfun,spinxx,spinzz,aspinxx,
     2       aspinzz,asafsq,asafsq2,asgnup,asgndn,asgn,anup,andn,
     3       anud,ake,den,aden,asaf,asaf2,asferro,asfer2,agrfun,
     4       errrat,asgnp,pairmat,apairmat,
     6       fpairmat,afpairmat
       common/m0vari/dopair,nmeasp,numpair,nmeas0,redo,noredo,nmeaspf
 
        double precision gmatup(0:n*n,0:n*n),gmatdn(0:n*n,0:n*n)
 
        integer i,j,mx,my,mpx,mpy
        integer inn,in0,in1,i0n,i01,i1n,i10,i11,jp(-1:1),jpp
 
c        if(mu .eq. 0.d0.and.bpar.eq.0.d0) then
c          do 10 i = 0,n*n-1
c           do 20 j = 0,n*n-1
c20           gmatdn(i,j) = -mphase(i)*mphase(j)*gmatup(j,i)
c10        gmatdn(i,i) = gmatdn(i,i) + 1.d0
c        endif
 
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
c*******************************meastausq********************
 
c        This subroutine does the measurements.
 
        subroutine meastausq
        implicit none
	include 'paramp.dat'
 
        integer warms,sweeps,nwrap,iran,iran0
        common/integers/warms,sweeps,nwrap,iran,iran0
 
         double precision u,mu,bpar,delt,dtau,gam,lambda
         common/couple/u,mu,bpar,delt,dtau,gam,lambda
 
c        index vectors for indirect addressing and phase vector
        integer xplus(0:volume-1),xminus(0:volume-1)
        integer yplus(0:volume-1),yminus(0:volume-1)
 
        common/ivectors/xplus,xminus,yplus,yminus
        save /ivectors/

c        common block for vectors
        double precision hub(0:volume-1)
      double precision vup(0:volume-1),vdn(0:volume-1)
      double precision mphase(0:toff-1)
      common/vectors/hub,vup,vdn,mphase

        double precision sgnup,sgndn,sgn,nup,ndn,nud,ke
        double precision saf,saf2,sferro,sfer2,grfun(0:n/2,0:n/2),
     1       safsq,safsq2
        double precision spinxx(0:n/2,0:n/2),spinzz(0:n/2,0:n/2)
        double precision pairmat(-1:1,-1:1,-1:1,-1:1),errrat,asgnp
        double precision asgnup,asgndn,asgn,anup,andn,anud,ake,
     1       asafsq,asafsq2
        double precision asaf,asaf2,asferro,asfer2,agrfun(0:n/2,0:n/2)
        double precision aspinxx(0:n/2,0:n/2),aspinzz(0:n/2,0:n/2)
        double precision apairmat(-1:1,-1:1,-1:1,-1:1)
        double precision fpairmat(-1:1,-1:1,-1:1,-1:1,0:n/2,0:n/2)
        double precision afpairmat(-1:1,-1:1,-1:1,-1:1,0:n/2,0:n/2)

        double precision den(0:n/2,0:n/2,0:1),aden(0:n/2,0:n/2,0:1)
        integer nmeas0,redo,noredo,dopair,nmeasp,numpair,nmeaspf
        common/m0var/ sgnup,sgndn,sgn,nup,ndn,nud,ke,safsq,safsq2,
     1       saf,saf2,sferro,sfer2,grfun,spinxx,spinzz,aspinxx,
     2       aspinzz,asafsq,asafsq2,asgnup,asgndn,asgn,anup,andn,
     3       anud,ake,den,aden,asaf,asaf2,asferro,asfer2,agrfun,
     4       errrat,asgnp,pairmat,apairmat,
     6       fpairmat,afpairmat
       common/m0vari/dopair,nmeasp,numpair,nmeas0,redo,noredo,nmeaspf
 
        double precision gnl(0:n/2,0:n/2,0:l),
     1      agnl(0:n/2,0:n/2,0:l)
        double precision chinl(0:n/2,0:n/2,0:l),
     1       achinl(0:n/2,0:n/2,0:l)
        double precision dent(0:n/2,0:n/2,0:l),
     1       adent(0:n/2,0:n/2,0:l)
        double precision chinlz(0:n/2,0:n/2,0:l),
     1       achinlz(0:n/2,0:n/2,0:l)
        double precision cnl(0:n/2,0:n/2,0:l),
     1       acnl(0:n/2,0:n/2,0:l)
        double precision pstau(0:l),pdtau(0:l),psxtau(0:l)
        double complex lam(0:n*n,0:n*n),alam(0:n*n,0:n*n)
	double complex ag1w(0:n*n,0:n*n)
        double precision apstau(0:l),apdtau(0:l),
     1       apsxtau(0:l),asgnt
        integer nmeast

        common/mtauvar/gnl,agnl,chinl,achinl,chinlz,achinlz,cnl,
     1   acnl,pstau,pdtau,psxtau,apstau,apdtau,apsxtau,lam,alam,
     2	 ag1w,asgnt,dent,adent
        common/mtauvari/nmeast
 
        integer llim(0:ndiv)
 
        integer ti,i,j,ix,jx,iy,jy,kx,ky,sx,sy,mx,my,jj
	integer lx,ly,dti
	integer ipx,ipy,ip1,ip2,jpx,jpy,jp1,jp2
        double precision bigmatup(0:toff*ndiv-1,0:toff*ndiv-1)
        double precision bigmatdn(0:toff*ndiv-1,0:toff*ndiv-1)
        double precision colgup(0:n*n,0:n*n,0:l)
        double precision rowgdn(0:n*n,0:n*n,0:l)
        double precision rowgup(0:n*n,0:n*n,0:l)
        double precision diagup(0:n*n,0:n*n,0:l)
        double precision diagdn(0:n*n,0:n*n,0:l)
        double precision colgdn(0:n*n,0:n*n,0:l)
        double precision fac
        double precision gplus(0:n*n,0:n*n,0:l), gminus(0:n*n,0:n*n,0:l)
c The following matrices are big only if doall is not 0
c	double precision hugup(0:n*n*doall,0:n*n*doall,0:l-1,0:l-1)
c	double precision hugdn(0:n*n*doall,0:n*n*doall,0:l-1,0:l-1)
	double complex giwup(0:n*n,0:n*n)
        double complex giwdn(0:n*n,0:n*n)
        double complex cfac
	double complex giwppup(0:n*n,0:n*n)
        double complex giwppdn(0:n*n,0:n*n)
	double complex cintfac(0:l),tup(0:n*n,0:n*n)
        double complex tdn(0:n*n,0:n*n)
 	double complex cint2(0:l)
        integer k,dk, xpi,xmi,ypi,ymi,xpj,xmj,ypj,ymj,doy
	integer writegi
        integer ip22(0:3)
        integer writea
	double precision arg

c     initialization flags
      integer ifbigt,ifftnt,ifgetc,ifgetr,ifmeas,ifgetd
      common/iflags/ifbigt,ifftnt,ifgetc,ifgetr,ifmeas,ifgetd

c     variable so we can avoid divide by zero errors when tausk=0
      integer tauskp
      common/tauskdiv/tauskp

        double precision ttup(0:toff,0:toff),ekup(0:toff,0:toff)
        double precision ekupi(0:toff,0:toff)
        double precision ttdn(0:toff,0:toff),ekdn(0:toff,0:toff)
        double precision ekdni(0:toff,0:toff)
        common/newke/ttup,ekup,ekupi,ttdn,ekdn,ekdni
      save /newke/

      save llim,cintfac
      save cint2

c Change doall to 1 to have it do all columns of the big tau-tauprime matrix.
      writegi=0
      writea=0
      if(ifmeas .ne. 373) then
         ifmeas = 373
         do 1 k = 0,ndiv
 1          llim(k) = (l * k)/ ndiv - 1
         if(writegi .eq. 1) then
            open(78,status='old',file="gidata")
            write(78,*)sweeps/tauskp,l
         endif

         if(writea .eq. 1) then
            open(79,status='old',file="pddata")
            if(tausk.ne.0) write(79,*)sweeps/tauskp,l
            open(80,status='old',file="psxdata")
            if(tausk.ne.0) write(80,*)sweeps/tauskp,l
            open(81,status='old',file="psdata")
            if(tausk.ne.0) write(81,*)sweeps/tauskp,l
         endif
         
         do 8833 k=0,l
            arg=twopi/l*0.5d0*k
            cintfac(k) = cdexp(dcmplx(0.d0,arg))/l
            if(k .eq. 0 .or. k .eq. l) then
               cintfac(k) = cintfac(k)*dtau/3.0d0
            else if(mod(k,2) .eq. 1) then
               cintfac(k) = cintfac(k)*dtau*4.0d0/3.0d0
            else
               cintfac(k) = cintfac(k)*dtau*2.d0/3.0d0
            endif
         cint2(k) = dconjg(cintfac(k))
 8833    continue
         ip22(0) = 1
         ip22(1) = 3
         ip22(2) = 0
         ip22(3) = 2
      endif

      do 14 ti = 0, l
         pdtau(ti) = 0.d0
         psxtau(ti) = 0.d0
         pstau(ti) = 0.d0
         do 14 kx=0,n/2
           do 14 ky=0,n/2
             gnl(kx,ky,ti) = 0.d0
             chinl(kx,ky,ti) = 0.d0
             dent(kx,ky,ti) = 0.d0
             chinlz(kx,ky,ti) = 0.d0
             cnl(kx,ky,ti) = 0.d0
14       continue
        if(doall .eq. 1) then
	 do 7887 i=0,n*n-1
	  do 7887 j=0,n*n-1
	  giwup(i,j) = 0.d0
7887	  giwdn(i,j) = 0.d0
        endif
        call biggtau(vup,bigmatup,1)
        call biggtau(vdn,bigmatdn,-1)
	call getdiag(bigmatup,diagup,vup,1)
 	call getdiag(bigmatdn,diagdn,vdn,-1)
        do 5 k = 0, ndiv-1
         jj = 0
1234     continue
         if(jj .eq. 10000)goto 5
         ti = mod(l+llim(k) + 1 + jj,l)
         if(jj .eq. 0) then
          call getcol(bigmatup,colgup,vup,k,1)
          call getcol(bigmatdn,colgdn,vdn,k,-1)
          call getrow(bigmatdn,rowgdn,vdn,k,-1)
          call getrow(bigmatup,rowgup,vup,k,1)
          if(doall .eq. 1)jj = -1
         else if(jj .lt. 0) then
          do 134 dk = l, 1, -1
           call matcop(colgup(0,0,dk-1),colgup(0,0,dk))
           call multb(colgup(0,0,dk),vup,ti,1,1)
           call matcop(colgdn(0,0,dk-1),colgdn(0,0,dk))
           call multb(colgdn(0,0,dk),vdn,ti,1,-1)
134       continue
          do 1341 dk = 0, l-1
           call matcop(rowgup(0,0,dk+1),rowgup(0,0,dk))
           call multbi(rowgup(0,0,dk),vup,ti,-1,1)
           call matcop(rowgdn(0,0,dk+1),rowgdn(0,0,dk))
1341       call multbi(rowgdn(0,0,dk),vdn,ti,-1,-1)
          call onesub(colgup(0,0,l),colgup(0,0,0))
          call onesub(colgdn(0,0,l),colgdn(0,0,0))
          call onesub(rowgdn(0,0,0),rowgdn(0,0,l))
          call onesub(rowgup(0,0,0),rowgup(0,0,l))
          jj = jj - 1
          if(k .eq. 0) then
           if(ti-1.le.(llim(ndiv)+llim(ndiv-1)+3)/2)jj=1
          else
           if(ti-1.le.(llim(k)+llim(k-1)+3)/2)jj=1
          endif
         else
          if(jj .eq. 1) then
            call getcol(bigmatup,colgup,vup,k,1)
            call getcol(bigmatdn,colgdn,vdn,k,-1)
            call getrow(bigmatdn,rowgdn,vdn,k,-1)
            call getrow(bigmatup,rowgup,vup,k,1)
          endif
          do 137 dk = 0, l-1
           call matcop(colgup(0,0,dk+1),colgup(0,0,dk))
           call multbi(colgup(0,0,dk),vup,ti-1,1,1)
           call matcop(colgdn(0,0,dk+1),colgdn(0,0,dk))
           call multbi(colgdn(0,0,dk),vdn,ti-1,1,-1)
137       continue
          do 1371 dk = l, 1, -1
          call matcop(rowgdn(0,0,dk-1),rowgdn(0,0,dk))
          call multb(rowgdn(0,0,dk),vdn,ti-1,-1,-1)
          call matcop(rowgup(0,0,dk-1),rowgup(0,0,dk))
1371      call multb(rowgup(0,0,dk),vup,ti-1,-1,1)
          call onesub(colgup(0,0,0),colgup(0,0,l))
          call onesub(colgdn(0,0,0),colgdn(0,0,l))
          call onesub(rowgdn(0,0,l),rowgdn(0,0,0))
c	  this is not in exact11e.f?
          call onesub(rowgup(0,0,l),rowgup(0,0,0))
          jj = jj + 1
          if(ti + 1 .gt. (llim(k)+llim(k+1)+3)/2) jj = 10000
         endif
	 if(doall .eq. 1)then
	    do 7888 dk=0,l
	     do 7888 i=0,n*n-1
	      do 7888 j=0,n*n-1
	       giwup(i,j) = giwup(i,j) + colgup(i,j,dk)*cintfac(dk)
 	       giwdn(i,j) = giwdn(i,j) + colgdn(i,j,dk)*cint2(dk)
7888	    continue
	 endif
	doy = 0
	sx = 1
	sy = 1
         do 10 kx = 0,n/2
         do 10 ky = 0, n/2
          do 10 ix = 0,n-1
          do 10 iy = 0,n-1
           jx = mod(ix+sx*kx+n,n)
           jy = mod(iy+sy*ky+n,n)
           i = ix+n*iy
           j = jx+n*jy
           do 151 dti = 0, l-ti-1
              dent(kx,ky,dti)=dent(kx,ky,dti)+
     1            (2.d0-diagup(i,i,dti+ti)-diagdn(i,i,dti+ti))*
     2            (2.d0-diagup(j,j,ti)-diagdn(j,j,ti))+
     3            colgup(i,j,dti)*rowgup(j,i,l-dti)+
     4            colgdn(i,j,dti)*rowgdn(j,i,l-dti)
151         chinlz(kx,ky,dti) = chinlz(kx,ky,dti)
     1       + (diagup(i,i,dti+ti)- diagdn(i,i,dti+ti))
     2       * (diagup(j,j,ti) - diagdn(j,j,ti))
     3       + colgup(i,j,dti)*rowgup(j,i,l-dti)
     4       + colgdn(i,j,dti)*rowgdn(j,i,l-dti)
           do 161 dti = l-ti, l
              dent(kx,ky,dti)=dent(kx,ky,dti)+
     1            (2.d0-diagup(i,i,dti+ti-l)-diagdn(i,i,dti+ti-l))*
     2            (2.d0-diagup(j,j,ti)-diagdn(j,j,ti))+
     3            colgup(i,j,dti)*rowgup(j,i,l-dti)+
     4            colgdn(i,j,dti)*rowgdn(j,i,l-dti)
161         chinlz(kx,ky,dti) = chinlz(kx,ky,dti)
     1       + (diagup(i,i,dti+ti-l)- diagdn(i,i,dti+ti-l))
     2       * (diagup(j,j,ti) - diagdn(j,j,ti))
     3       + colgup(i,j,dti)*rowgup(j,i,l-dti)
     4       + colgdn(i,j,dti)*rowgdn(j,i,l-dti)
	    ipx = mod(ix+1+n,n)
	    ipy = mod(iy+1+n,n)
	    jpx = mod(jx+1+n,n)
	    jpy = mod(jy+1+n,n)
	    ip1 = ipx + n*iy
	    jp1 = jpx + n*jy
	    ip2 = ix + n*ipy
	    jp2 = jx + n*jpy
	    if(doy .eq. 1)then
	      ip1 = ip2
	      jp1 = jp2
	    endif
	    if(n .eq. 2)then
c The following make the current go around the loop for 2x2.
	      ip1 = ip22(i)
	      jp1 = ip22(j)
	    endif

cpd00  Measurement of current-current correlation changed to discriminate
cpd00  spin up/down in tt arrays (previous version commented out).

            do 15 dti = 0, l-ti-1
c             cnl(kx,ky,dti) = cnl(kx,ky,dti) +
c     1          tt(i,ip1)*tt(j,jp1)*
c     1          (diagup(i,ip1,dti+ti)+diagdn(i,ip1,dti+ti)
c     2           -diagup(ip1,i,dti+ti)-diagdn(ip1,i,dti+ti))
c     3        * (diagup(j,jp1,ti)+diagdn(j,jp1,ti)
c     4          -diagup(jp1,j,ti)-diagdn(jp1,j,ti) )

             cnl(kx,ky,dti) = cnl(kx,ky,dti) +
     1    (ttup(i,ip1)*(diagup(i,ip1,dti+ti)-diagup(ip1,i,dti+ti))
     2  +  ttdn(i,ip1)*(diagdn(i,ip1,dti+ti)-diagdn(ip1,i,dti+ti)) )
     3  * (ttup(j,jp1)*(diagup(j,jp1,ti)-diagup(jp1,j,ti))
     4  +  ttdn(j,jp1)*(diagdn(j,jp1,ti)-diagdn(jp1,j,ti)) )

15	    continue
            do 16 dti = l-ti,l
c             cnl(kx,ky,dti) = cnl(kx,ky,dti) +
c     1          tt(i,ip1)*tt(j,jp1)*
c     1          (diagup(i,ip1,dti+ti-l)+diagdn(i,ip1,dti+ti-l)
c     2           -diagup(ip1,i,dti+ti-l)-diagdn(ip1,i,dti+ti-l))
c     3        * (diagup(j,jp1,ti)+diagdn(j,jp1,ti)
c     4          -diagup(jp1,j,ti)-diagdn(jp1,j,ti) )

             cnl(kx,ky,dti) = cnl(kx,ky,dti) +
     1    (ttup(i,ip1)*(diagup(i,ip1,dti+ti-l)-diagup(ip1,i,dti+ti-l))
     2  +  ttdn(i,ip1)*(diagdn(i,ip1,dti+ti-l)-diagdn(ip1,i,dti+ti-l)) )
     3  * (ttup(j,jp1)*(diagup(j,jp1,ti)-diagup(jp1,j,ti))
     4  +  ttdn(j,jp1)*(diagdn(j,jp1,ti)-diagdn(jp1,j,ti)) )

16	    continue
            do 17 dti = 0, l
c             cnl(kx,ky,dti) = cnl(kx,ky,dti) + 
c     1        tt(i,ip1)*tt(j,jp1)*
c     1        (colgup(i,jp1,dti)*rowgup(j,ip1,l-dti)
c     2        - colgup(i,j,dti)*rowgup(jp1,ip1,l-dti)
c     3        - colgup(ip1,jp1,dti)*rowgup(j,i,l-dti)
c     4        + colgup(ip1,j,dti)*rowgup(jp1,i,l-dti)
c     1        + colgdn(i,jp1,dti)*rowgdn(j,ip1,l-dti)
c     2        - colgdn(i,j,dti)*rowgdn(jp1,ip1,l-dti)
c     3        - colgdn(ip1,jp1,dti)*rowgdn(j,i,l-dti)
c     4        + colgdn(ip1,j,dti)*rowgdn(jp1,i,l-dti) )

             cnl(kx,ky,dti) = cnl(kx,ky,dti) + 
     1        ttup(i,ip1)*ttup(j,jp1)*
     1        (colgup(i,jp1,dti)*rowgup(j,ip1,l-dti)
     2        - colgup(i,j,dti)*rowgup(jp1,ip1,l-dti)
     3        - colgup(ip1,jp1,dti)*rowgup(j,i,l-dti)
     4        + colgup(ip1,j,dti)*rowgup(jp1,i,l-dti) )
     1        + ttdn(i,ip1)*ttdn(j,jp1)*
     1        ( colgdn(i,jp1,dti)*rowgdn(j,ip1,l-dti)
     2        - colgdn(i,j,dti)*rowgdn(jp1,ip1,l-dti)
     3        - colgdn(ip1,jp1,dti)*rowgdn(j,i,l-dti)
     4        + colgdn(ip1,j,dti)*rowgdn(jp1,i,l-dti) )

17	    continue
           do 10 dti = 0, l
             gnl(kx,ky,dti) = gnl(kx,ky,dti) +
     1         colgup(i,j,dti) + colgdn(i,j,dti)
             chinl(kx,ky,dti) = chinl(kx,ky,dti) +
     1         colgup(j,i,dti)*rowgdn(i,j,l-dti) +
     1         colgdn(j,i,dti)*rowgup(i,j,l-dti)
10         continue
           
        if(dopair .eq. 1) then
c New section of code:
         do 522 dti = 0, l
          do 522 i = 0,n*n-1
           do 522 j = 0,n*n-1
             gplus(i,j,dti) = 0.d0
522          gminus(i,j,dti) = 0.d0
         do 510 i = 0,n*n-1
          xpi = xplus(i)
          xmi = xminus(i)
          ypi = yplus(i)
          ymi = yminus(i)
          do 510 j = 0,n*n-1
           xpj = xplus(j)
           xmj = xminus(j)
           ypj = yplus(j)
           ymj = yminus(j)
           do 510 dti = 0, l
             gplus(i,j,dti) = gplus(i,j,dti)
     1              +  colgup(xpi,xpj,dti)
     2              +  colgup(xmi,xmj,dti)
     3              +  colgup(xpi,xmj,dti)
     4              +  colgup(xmi,xpj,dti)
     1              +  colgup(ypi,ypj,dti)
     2              +  colgup(ymi,ymj,dti)
     3              +  colgup(ypi,ymj,dti)
     4              +  colgup(ymi,ypj,dti)
             gminus(i,j,dti) = gminus(i,j,dti)
     1              +  colgup(xpi,ypj,dti)
     2              +  colgup(xmi,ymj,dti)
     3              +  colgup(xpi,ymj,dti)
     4              +  colgup(xmi,ypj,dti)
     1              +  colgup(ypi,xpj,dti)
     2              +  colgup(ymi,xmj,dti)
     3              +  colgup(ypi,xmj,dti)
     4              +  colgup(ymi,xpj,dti)
510      continue
         do 512 i = 0,n*n-1
          do 512 j = 0,n*n-1
           do 512 dti = 0, l
             pdtau(dti) = pdtau(dti) + 
     1             colgdn(i,j,dti) * (gplus(i,j,dti)
     2             - gminus(i,j,dti))
             psxtau(dti) = psxtau(dti) +
     1            colgdn(i,j,dti)*(gplus(i,j,dti)
     2            + gminus(i,j,dti))
             pstau(dti) = pstau(dti) + 
     1            colgdn(i,j,dti) * colgup(i,j,dti)
512      continue
         do 513 dti = 0, l
513      continue
        endif
        if(doall .eq. 1)goto 1234
5       continue
 
        if(doall .eq. 1) then
	 do 9887 i=0,n*n-1
	  do 9887 j=0,n*n-1
	  tup(i,j) = 0.d0
 	  tdn(i,j) = 0.d0
    	  giwppup(i,j) = 0.d0
9887	  giwppdn(i,j) = 0.d0
	  do 9888 mx=0,n-1
	   do 9888 my=0,n-1
	    do 9888 lx=0,n-1
	     do 9888 ly=0,n-1
	      arg  = twopi/n*(mx*lx+my*ly)
	      cfac = cdexp(dcmplx(0.d0,arg))/n
	      do 9888 j=0,n*n-1
            tup(j,mx+n*my) = tup(j,mx+n*my)+
     1                giwup(j,lx+ly*n)*cfac
9888	       tdn(j,mx+n*my) = tdn(j,mx+n*my)+giwdn(j,lx+ly*n)*cfac
	  do 9988 mx=0,n-1
	   do 9988 my=0,n-1
	    do 9988 lx=0,n-1
	     do 9988 ly=0,n-1
	      arg = (-twopi)/n*(mx*lx+my*ly)
	      cfac = cdexp(dcmplx(0.d0,arg))/n
	      do 9988 j=0,n*n-1
	   giwppup(mx+n*my,j)=giwppup(mx+n*my,j)+
     1                tup(lx+ly*n,j)*cfac
9988	   giwppdn(mx+n*my,j)=giwppdn(mx+n*my,j)+tdn(lx+ly*n,j)*cfac
	  do 9989 mx=0,n-1
	   do 9989 my=0,n-1
	    i = mx + n*my
	    kx = mod(n-mx,n)+n*mod(n-my,n)
	    do 9989 lx=0,n-1
	     do 9989 ly=0,n-1
	     j = lx + n*ly
	     ky = mod(n-lx,n)+n*mod(n-ly,n)
9989	     lam(i,j) = giwppup(i,j)*giwppdn(kx,ky)
	  do 9998 i=0,n*n-1
	   do 9998 j=0,n*n-1
	    ag1w(i,j)=ag1w(i,j)+giwppup(i,j)*sgnup*sgndn
9998	    alam(i,j) = alam(i,j) + lam(i,j)*sgnup*sgndn
        endif
        fac = 1.d0 / (n*n*ndiv)
        if(doall .eq. 1) fac = 1.d0 / (n*n*l)
        do 30 kx = 0,n/2
         do 30 ky = 0, n/2
          do 30 dti = 0, l
           gnl(kx,ky,dti) = gnl(kx,ky,dti) * 0.5d0 * fac
           chinl(kx,ky,dti) = chinl(kx,ky,dti) * 1.d0 * fac
           chinlz(kx,ky,dti) = chinlz(kx,ky,dti) * fac
           dent(kx,ky,dti) = dent(kx,ky,dti) * fac
           cnl(kx,ky,dti) = cnl(kx,ky,dti) * fac
30       continue
 
        nmeast = nmeast + 1
        asgnt = asgnt + sgnup*sgndn
        do 594 mx = 0, n/2
         do 594 my = 0, n/2
          do 594 dti = 0,l
           agnl(mx,my,dti)=agnl(mx,my,dti)+
     1            gnl(mx,my,dti)*sgnup*sgndn
           achinl(mx,my,dti)=achinl(mx,my,dti)+
     1                          chinl(mx,my,dti)*sgnup*sgndn
           adent(mx,my,dti)=adent(mx,my,dti)+
     1                          dent(mx,my,dti)*sgnup*sgndn
           acnl(mx,my,dti)=acnl(mx,my,dti)+
     1                          cnl(mx,my,dti)*sgnup*sgndn
594        achinlz(mx,my,dti)=achinlz(mx,my,dti)+
     1                          chinlz(mx,my,dti)*sgnup*sgndn
        fac = 1.d0 / ndiv
        if(doall .eq. 1) fac = 1.d0 / l
        fac = fac * 0.25d0 / (n*n)
	if(writea .eq. 1)write(79,*)sgnup*sgndn
	if(writea .eq. 1)write(80,*)sgnup*sgndn
	if(writea .eq. 1)write(81,*)sgnup*sgndn
	if(writegi .eq. 1)write(78,*)sgnup*sgndn
        do 494 dti = 0, l
          pdtau(dti)=pdtau(dti)*fac
          psxtau(dti)=psxtau(dti)*fac
          pstau(dti)=pstau(dti)*fac * 4.d0
1093	  format(g16.8)
	  if(writegi .eq. 1)write(78,1093)gnl(0,0,dti)
	  if(writea .eq. 1)write(79,1093)pdtau(dti)
	  if(writea .eq. 1)write(80,1093)psxtau(dti)
	  if(writea .eq. 1)write(81,1093)pstau(dti)
          apdtau(dti)=apdtau(dti) + sgnup*sgndn*pdtau(dti)
          apstau(dti)=apstau(dti) + sgnup*sgndn*pstau(dti)
494       apsxtau(dti)=apsxtau(dti) + sgnup*sgndn*psxtau(dti)
         

        return
        end

c*******************************pairfoursq(gmatup)********************
 
c        This subroutine does the measurements of pair functions
c        allowing to Fourier transform the results.
 
        subroutine pairfoursq(gmatup,gmatdn)
        implicit none
        include 'paramp.dat'
 
         double precision u,mu,bpar,delt,dtau,gam,lambda
         common/couple/u,mu,bpar,delt,dtau,gam,lambda
 
c        index vectors for indirect addressing and phase vector
        integer xplus(0:volume-1),xminus(0:volume-1)
        integer yplus(0:volume-1),yminus(0:volume-1)
        common/ivectors/xplus,xminus,yplus,yminus
        save /ivectors/

        double precision sgnup,sgndn,sgn,nup,ndn,nud,ke
        double precision saf,saf2,sferro,sfer2,grfun(0:n/2,0:n/2),
     1       safsq,safsq2
        double precision spinxx(0:n/2,0:n/2),spinzz(0:n/2,0:n/2)
        double precision pairmat(-1:1,-1:1,-1:1,-1:1),errrat,asgnp
        double precision fpairmat(-1:1,-1:1,-1:1,-1:1,0:n/2,0:n/2)
        integer nmeas0,redo,noredo,dopair,nmeasp,numpair,nmeaspf
        double precision asgnup,asgndn,asgn,anup,andn,anud,ake,
     1       asafsq,asafsq2
        double precision asaf,asaf2,asferro,asfer2,agrfun(0:n/2,0:n/2)
        double precision aspinxx(0:n/2,0:n/2),aspinzz(0:n/2,0:n/2)
        double precision apairmat(-1:1,-1:1,-1:1,-1:1)
        double precision afpairmat(-1:1,-1:1,-1:1,-1:1,0:n/2,0:n/2)
        double precision den(0:n/2,0:n/2,0:1),aden(0:n/2,0:n/2,0:1)
        common/m0var/ sgnup,sgndn,sgn,nup,ndn,nud,ke,safsq,safsq2,
     1       saf,saf2,sferro,sfer2,grfun,spinxx,spinzz,aspinxx,
     2       aspinzz,asafsq,asafsq2,asgnup,asgndn,asgn,anup,andn,
     3       anud,ake,den,aden,asaf,asaf2,asferro,asfer2,agrfun,
     4       errrat,asgnp,pairmat,apairmat,
     6       fpairmat,afpairmat
       common/m0vari/dopair,nmeasp,numpair,nmeas0,redo,noredo,nmeaspf
 
        double precision xfac,yfac
 
cc        common block for vectors
        double precision hub(0:volume-1)
       double precision vup(0:volume-1),vdn(0:volume-1)
       double precision mphase(0:toff-1)
       common/vectors/hub,vup,vdn,mphase

        double precision gmatup(0:n*n,0:n*n)
        double precision gmatdn(0:n*n,0:n*n)
 
        integer i,j,ix,jx,iy,jy,kx,ky,mx,my,mpx,mpy
        integer inn,in0,in1,i0n,i01,i1n,i10,i11,jp(-1:1),jpp
 
c        if(mu.eq.0.d0.and.bpar.eq.0.d0) then
c          do 10 i = 0,n*n-1
c           do 20 j = 0,n*n-1
c20           gmatdn(i,j) = -mphase(i)*mphase(j)*gmatup(j,i)
c10        gmatdn(i,i) = gmatdn(i,i) + 1.d0
c        endif
 
        do 294 mx = -1,1
         do 294 my = -1,1
          do 294 mpx = -1,1
           do 294 mpy = -1,1
            do 294 kx =0,n/2
             do 294 ky =0,n/2
294          fpairmat(mx,my,mpx,mpy,kx,ky) = 0.d0
        do 522 ix = 0,n-1
        do 522 iy = 0,n-1
         i=ix+n*iy
          inn = xminus(yminus(i))
          in0 = xminus(i)
          in1 = xminus(yplus(i))
          i0n = yminus(i)
          i01 = yplus(i)
          i1n = xplus(yminus(i))
          i10 = xplus(i)
          i11 = xplus(yplus(i))
          do 522 jx = 0,n-1
          do 522 jy = 0,n-1
           j=jx+n*jy
           kx=iabs(ix-jx)
           kx=min(kx,n-kx)
           ky=iabs(iy-jy)
           ky=min(ky,n-ky)
            jpp = xplus(xplus(j))
            do 500 mx = -1,1
             jpp = xminus(jpp)
             jp(-1) = yplus(jpp)
             jp(0) = jpp
             jp(1) = yminus(jpp)
             do 500 my = -1,1
           fpairmat(mx,my,-1,-1,kx,ky) = fpairmat(mx,my,-1,-1,kx,ky) +
     1             (gmatup(jp(my),i11)*gmatdn(j,i)  +
     1               gmatup(jp(my),i)*gmatdn(j,i11) +
     1               gmatup(j,i11)*gmatdn(jp(my),i) +
     1               gmatup(j,i)*gmatdn(jp(my),i11) )*0.25d0

           fpairmat(mx,my,-1,0,kx,ky) = fpairmat(mx,my,-1,0,kx,ky) +
     1             (gmatup(jp(my),i10)*gmatdn(j,i)  +
     1               gmatup(jp(my),i)*gmatdn(j,i10) +
     1               gmatup(j,i10)*gmatdn(jp(my),i) +
     1               gmatup(j,i)*gmatdn(jp(my),i10) )*0.25d0

           fpairmat(mx,my,-1,1,kx,ky) = fpairmat(mx,my,-1,1,kx,ky) +
     1             (gmatup(jp(my),i1n)*gmatdn(j,i)  +
     1               gmatup(jp(my),i)*gmatdn(j,i1n) +
     1               gmatup(j,i1n)*gmatdn(jp(my),i) +
     1               gmatup(j,i)*gmatdn(jp(my),i1n) )*0.25d0

           fpairmat(mx,my,0,-1,kx,ky) = fpairmat(mx,my,0,-1,kx,ky) +
     1             (gmatup(jp(my),i01)*gmatdn(j,i)  +
     1               gmatup(jp(my),i)*gmatdn(j,i01) +
     1               gmatup(j,i01)*gmatdn(jp(my),i) +
     1               gmatup(j,i)*gmatdn(jp(my),i01) )*0.25d0

           fpairmat(mx,my,0,0,kx,ky) = fpairmat(mx,my,0,0,kx,ky) +
     1             (gmatup(jp(my),i)*gmatdn(j,i) +
     1               gmatup(jp(my),i)*gmatdn(j,i) +
     1               gmatup(j,i)*gmatdn(jp(my),i) +
     1               gmatup(j,i)*gmatdn(jp(my),i))*0.25d0

           fpairmat(mx,my,0,1,kx,ky) = fpairmat(mx,my,0,1,kx,ky) +
     1             (gmatup(jp(my),i0n)*gmatdn(j,i) +
     1               gmatup(jp(my),i)*gmatdn(j,i0n) +
     1               gmatup(j,i0n)*gmatdn(jp(my),i) +
     1               gmatup(j,i)*gmatdn(jp(my),i0n) )*0.25d0

           fpairmat(mx,my,1,-1,kx,ky) = fpairmat(mx,my,1,-1,kx,ky) +
     1             (gmatup(jp(my),in1)*gmatdn(j,i)  +
     1               gmatup(jp(my),i)*gmatdn(j,in1) +
     1               gmatup(j,in1)*gmatdn(jp(my),i) +
     1               gmatup(j,i)*gmatdn(jp(my),in1))*0.25d0


           fpairmat(mx,my,1,0,kx,ky) = fpairmat(mx,my,1,0,kx,ky) +
     1             (gmatup(jp(my),in0)*gmatdn(j,i) +
     1               gmatup(jp(my),i)*gmatdn(j,in0) +
     1               gmatup(j,in0)*gmatdn(jp(my),i) +
     1               gmatup(j,i)*gmatdn(jp(my),in0) )*0.25d0

           fpairmat(mx,my,1,1,kx,ky) = fpairmat(mx,my,1,1,kx,ky) +
     1             (gmatup(jp(my),inn)*gmatdn(j,i)  +
     1               gmatup(jp(my),i)*gmatdn(j,inn) +
     1               gmatup(j,inn)*gmatdn(jp(my),i) +
     1               gmatup(j,i)*gmatdn(jp(my),inn) )*0.25d0

500         continue
522     continue
 
        nmeaspf= nmeaspf+ 1
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
           do 42 mx=-1,1
           do 42 my=-1,1
           do 42 mpx=-1,1
           do 42 mpy=-1,1
42            fpairmat(mx,my,mpx,mpy,kx,ky)=
     1              fpairmat(mx,my,mpx,mpy,kx,ky)/(n*n*xfac*yfac)
41      continue
 
        do 31 kx = 0, n/2
          do 31 ky = 0, n/2
           do 31 mx=-1,1
           do 31 my=-1,1
           do 31 mpx=-1,1
           do 31 mpy=-1,1
            fpairmat(mx,my,mpx,mpy,kx,ky)=(fpairmat(mx,my,mpx
     1      ,mpy,kx,ky)+fpairmat(mx,my,mpx,mpy,ky,kx))*0.5d0
31      continue
 
        do 51 kx = 0, n/2
          do 51 ky = 0, n/2
           do 51 mx=-1,1
           do 51 my=-1,1
           do 51 mpx=-1,1
           do 51 mpy=-1,1
           afpairmat(mx,my,mpx,mpy,kx,ky)=afpairmat(mx,my,mpx
     1      ,mpy,kx,ky)+fpairmat(mx,my,mpx,mpy,ky,kx)*sgnup*sgndn
51       continue

        return
        end
c %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


      SUBROUTINE JACOBI(A,N,NP,D,V,NROT)
      implicit real*8(a-h,o-z)
      PARAMETER (NMAX=200)
      double precision A(NP,NP),D(NP),V(NP,NP),B(NMAX),Z(NMAX)
      DO 12 IP=1,N
        DO 11 IQ=1,N
          V(IP,IQ)=0.d0
11      CONTINUE
        V(IP,IP)=1.d0
12    CONTINUE
      DO 13 IP=1,N
        B(IP)=A(IP,IP)
        D(IP)=B(IP)
        Z(IP)=0.d0
13    CONTINUE
      NROT=0
      DO 24 I=1,50
        SM=0.d0
        DO 15 IP=1,N-1
          DO 14 IQ=IP+1,N
            SM=SM+DABS(A(IP,IQ))
14        CONTINUE
15      CONTINUE
        IF(SM.EQ.0.d0)RETURN
        IF(I.LT.4)THEN
          TRESH=0.2d0*SM/dfloat(N)**2
        ELSE
          TRESH=0.d0
        ENDIF
        DO 22 IP=1,N-1
          DO 21 IQ=IP+1,N
            G=100.d0*DABS(A(IP,IQ))
            IF((I.GT.4).AND.(DABS(D(IP))+G.EQ.DABS(D(IP)))
     *         .AND.(DABS(D(IQ))+G.EQ.DABS(D(IQ))))THEN
              A(IP,IQ)=0.d0
            ELSE IF(DABS(A(IP,IQ)).GT.TRESH)THEN
              H=D(IQ)-D(IP)
              IF(DABS(H)+G.EQ.DABS(H))THEN
                T=A(IP,IQ)/H
              ELSE
                THETA=0.5d0*H/A(IP,IQ)
                T=1.d0/(DABS(THETA)+DSQRT(1.d0+THETA**2.d0))
                IF(THETA.LT.0.d0)T=-T
              ENDIF
              C=1.d0/DSQRT(1.d0+T**2.d0)
              S=T*C
              TAU=S/(1.d0+C)
              H=T*A(IP,IQ)
              Z(IP)=Z(IP)-H
              Z(IQ)=Z(IQ)+H
              D(IP)=D(IP)-H
              D(IQ)=D(IQ)+H
              A(IP,IQ)=0.d0
              DO 16 J=1,IP-1
                G=A(J,IP)
                H=A(J,IQ)
                A(J,IP)=G-S*(H+G*TAU)
                A(J,IQ)=H+S*(G-H*TAU)
16            CONTINUE
              DO 17 J=IP+1,IQ-1
                G=A(IP,J)
                H=A(J,IQ)
                A(IP,J)=G-S*(H+G*TAU)
                A(J,IQ)=H+S*(G-H*TAU)
17            CONTINUE
              DO 18 J=IQ+1,N
                G=A(IP,J)
                H=A(IQ,J)
                A(IP,J)=G-S*(H+G*TAU)
                A(IQ,J)=H+S*(G-H*TAU)
18            CONTINUE
              DO 19 J=1,N
                G=V(J,IP)
                H=V(J,IQ)
                V(J,IP)=G-S*(H+G*TAU)
                V(J,IQ)=H+S*(G-H*TAU)
19            CONTINUE
              NROT=NROT+1
            ENDIF
21        CONTINUE
22      CONTINUE
        DO 23 IP=1,N
          B(IP)=B(IP)+Z(IP)
          D(IP)=B(IP)
          Z(IP)=0.d0
23      CONTINUE
24    CONTINUE
      Write(63,*) '50 iterations should never happen'
      RETURN
      END



