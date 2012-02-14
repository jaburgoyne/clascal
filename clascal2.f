c
c
c      program clascal
c-------------------------------------------------------------------------------
c                     CLASCAL
c               MDS Model With Latent Classes for WTs
c               --------------------------------------------------
c
c  written by Suzanne Winsberg -- Version 7.01 -- may 1993
c				(when changing the version, update
c				1st output format in getpar.f)
c
c  Input sequence:
c  --------------
c
c     1) Descriptive job title (a80)
c     2) Run parameter record (12i5)
c           nsub    = no. of subjects  
c           nstim   = no. of stimuli (> 2)
c           nlo  = Lower no. of latent classes (> 0)
c           nhi  = Upper no. of latent classes (> 0)
c           ns  = 0: no specificities (see A latent class approach to
c   fitting the weighted euclidean model,CLASCAL in Psychometrika june '93
c   by S Winsberg and G de Soete)
c                 1 : specificities (see the restricted Exscal Model in
c   Fitting an extended Indscal model to three way proximity data to
c   appear in Jrnl of Classification  by JD Carroll and S. Winsberg;
c   this model is applied here with a different weight for the specificites
c   for each class rather than each source 
c                 2 : different pattern of specificities(see the Exscal Model
c   in Fitting an extended Indscal model to three way proximity data to
c   appear in Jrnl of Classification  by JD Carroll and S. Winsberg;
c   this model is applied here with a different specificity for each latent
c   class rather than each source)
c           ndimlo    = Lower no. of dimensions
c           ndimhi    = Upper no. of dimensions
c           istart  = 0 rational start (k means) for class sort
c                   = 1 random start for class sort
c           iseed    = seed for random starts and monte carlo tests >0
c           nrep  =  the no. of samples generated for monte carlo tests
c                      if nrep < 20  no monte carlo tests
c           icarlo
c                  = 0  the monte caro is on the null model
c                  = 1  the monte carlo is on the spatail model 
c           iconsc =0 no conscal comparison
c                  =1 conscal comparison this option is only valid for
c   the case when there is no monte carlo testing
c     3) Logical Variables
c        pgword = .true. for extra print info
c         rswrd = .false. for rational start (torgerson) for spatial parameters
c               = .true. for random start
c         indscwd  =  .true. for indscal start
c                  =  .false. for any other start   
c     4) Set Iteration Constants
c          sst = 1.0 start value for specificities
c     5)   itmax = max no of iterations for m step
c           itxmax= max no of iterations for spatial step
c           itwmax= max no of iterations for weight step
c           conv= convergence criterion for spatial step
c           wconv=convergenc criterion for weight step         
c     6) Data format: format for reading one row of the 
c        paired dissimilarity matrix as floating point numbers (a80)
c        Each row has the data for one subject. So the data in a row
c        is the lower triangle of the dissimilarities re~ad by rows
c     7) Paired Dissimilarity Matrix, entry (i,j) indicates the dissimilarity
c        for the ith subject and the jth pair (j =1 for 21 j =2 for 31 j=3 for
c        32 j=4 for 41 j=5 for 42 j= 6 for 43 j=7 for 51 etc.)
c     8a) if iconcl=1 the value of fit1 for conscal analysis
c     8b) if monte carlo testing on the spatial model the values of
c          nclass no of dimensions and specificity parameter ns for
c          the comparison model in format 3i5
c
c
c     Remarks
c     -------
c
c     The program uses subroutines from the port library for
c     error handling, dynamic storage allocation and specifica-
c     tion of machine-dependent constants
c     cf. Fox, P.A., Hall, A.D., & Schryer, N.L.
c           - The port mathematical subroutine libray. ACM Trans.
c             on Math. Softw., 1978, 4, 104-126.
c           - Algorithm 528. Framework for a portable library.
c             ACM Trans. on Math. Softw., 1978, 4, 177-188.
c 
c     To change the amount of memory available change the following
c     statements below:
c          double precision(NITEM)
c          real rstak(2*NITEM)
c          integer istak(2*NITEM)
c
c          call istkin(NITEM,4)
c     where NITEM indicates the no. of double precision storage units
c     that are desired. No other changes are required.
c
c-------------------------------------------------------------------------------
c
      common /cstak/ dstak
c
      double precision dstak(300000000)
      real rstak(600000000)
      integer istak(600000000)
c
      double precision  loglik,loglikp
      common/parblk/nstim,ndim,nsub,nwgt,nmodel,ndat,
     1nsd,nspr,ns,nspace
      common/cnvblk/conv,xsconv,wconv,
     1itmax,itxmax,itwmax
      logical test,wword,abort,pgword,rswrd,indscwd
c
      equivalence (dstak(1),rstak(1))
      equivalence (dstak(1),istak(1))
      data ndimmax / 10 /
c
      call istkin(300000000,4)
      lp=i1mach(2)
c
      call gtpar(nsub,nstim,nspr,ndimlo,ndimhi,ns,
     1    nlo,nhi,istart,iseed,nrep,icarlo,iconsc)
c
c  compute various parameters
      nspr=(nstim*(nstim-1))/2
      ndat=(nstim*(nstim-1)*nsub)/2
c   set logical variables
      read *,pgword,rswrd,indscwd
      abort=.false.
c   set up convergence parameters
	call conset(sst,rswrd,indscwd,nrep)
c  allocate memory for distance data
      jy=istkgt(nsub*nspr, 3)
      call gtdata(rstak(jy),nsub,nspr)
      do 800 ndim=ndimlo,ndimhi
        if((nrep.lt.20).or.(icarlo.eq.1))print 1211,ndim
 1211   format(//'Analysis with No. of dimensions:',i5)
      nsd=nstim*ndim
      nspace =nsd
      nmax= nspace + nstim
c  allocate memory and set start values for coordinates of spacial models
         jx=istkgt(nmax,3)
         jwin=istkgt(nsub*ndim,3)
c  allocate memory for random start
          nspr2=2*nspr
          nstim2=2*nstim
         jaa=istkgt(nspr2,3)
         jbb=istkgt(nspr2,3)
         ja=istkgt(nstim2,4)
         jb=istkgt(nstim2,4)
         jc=istkgt(nstim2,4)
         nsq=nstim2*nstim2
         nsq1=istkgt(nsq,4)
         nsq2=istkgt(nsq,4)
         nsq3=istkgt(nsq,4)  
      if(indscwd) then
        call gtvdata(rstak(jx),nstim,ndim,nmax)
        call gtdata(rstak(jwin),nsub,ndim)
      else
      if(rswrd) then
        do 17 i=1,nsd
   17    rstak(jx+i-1)= -1.0 +2.0*runif(iseed)
      else
	call xstart(nstim2,nspr2,nmax,rstak(jx),rstak(jy),dstak(ja),
     1      dstak(jb),dstak(jc),dstak(nsq1),dstak(nsq2),
     2      dstak(nsq3),rstak(jaa),rstak(jbb),pgword)
        end if
      end if
      call istkrl(8)
      do 500 nclass=nlo,nhi
         jwmean=istkgt(nclass*ndim,3)
       if(ns.eq.1)nspace=nspace+nstim
       if(ns.eq.2)nspace=nspace+nstim*nclass
          if((nclass.eq.1).and.(ns.eq.0))nmodel=0
          if((nclass.eq.1).and.(ns.ge.1))nmodel=2
          if((nclass.gt.1).and.(ns.eq.0))nmodel=1
          if((nclass.gt.1).and.(ns.eq.1))nmodel=3
          if((nclass.gt.1).and.(ns.eq.2))nmodel=4
           wword=.true.
           if(nclass.eq.1)wword=.false.

        test=.false.
        print 1212,nclass
 1212   format(//'No. of classes:',i5)
	nw=nclass*ndim
	if(nmodel.eq.3)nw=nw+nclass
        nwgt=nw
        nmaxw=nw + ndim + 1
c
c  allocate memory for model params
      jlam=istkgt(nclass,3)
      jmean=istkgt(nclass*nspr,3)
c
c
      jz=istkgt(nsub*nclass,3)
        jw=istkgt(nmaxw,3)
        if(indscwd) then
         call wstart(nmodel,rstak(jlam),rstak(jwin),rstak(jw),nsub,
     1ndim,nclass,rstak(jwmean),rstak(jmean),nspr,nmaxw,rstak(jy)
     2       ,sigma2)
        else
         jend=jw+nw-1
          do 27 ijq=jw,jend
   27          rstak(ijq)=1.0
        end if
      if(nsub.eq.nclass)call zinit(nsub,nspr,rstak(jz),
     1    rstak(jy),rstak(jlam),rstak(jmean))
      jtmp=istkgt(nclass,4)
      if(nsub.eq.nclass) go to 57
      if(nclass.eq.1) go to 23
      if(indscwd) go to 43
c  rational start
      if (istart .eq. 0) call ratin(rstak(jlam),rstak(jmean),
     1   sigma2,rstak(jy),nsub,nspr,nclass)
c  random start
      if (istart .gt. 0) call ranin(rstak(jlam),rstak(jmean),
     1   sigma2,rstak(jy),nsub,nspr,nclass,istart)
         go to 43
   23   call onecl(rstak(jlam),rstak(jmean),sigma2,rstak(jy),
     1      nsub,nspr)
   43     continue
c  no of parameters for null model
          npar=nspr*nclass + nclass
      call logl(rstak(jlam),rstak(jmean),sigma2,rstak(jy),
     1   rstak(jz),nclass,nspr,nsub,dstak(jtmp),loglik,loglikp,
     1    npar,aic,bic)
      if((ndim.eq.ndimlo).and.(nrep.gt.19).and.(icarlo.eq.0))  then
         
         ndown=nclass+1
         nup=nhi+1
         do 56 ncli=ndown,nup
               call rundvrp(
     1    rstak(jlam),rstak(jmean),sigma2,rstak(jy),rstak(jz),
     1    nclass,ncli,logli,dstak(jtmp),npar,aic,bic,nrep,iseed,test)
 56            continue
      else
      write(6,50)
 50   format(//' Initial parameter estimates Null Model')
      call prtpar(rstak(jlam),rstak(jmean),sigma2,nclass,nspr,loglik,
     1    aic,bic)
      end if
c
 57   continue
c  compute no of active parameters
c  adjust for model
       go to(10,20,30,40,45)nmodel+1
   10    npar=nsd-ndim*(ndim+1)/2
         go to 33
   20    npar=nsd+ndim*nclass-2*ndim
         go to 33
   30    npar=nsd-ndim*(ndim+1)/2+nstim
         go to 33
   40    npar=nsd+nclass*(ndim+1)-2*ndim-1+nstim
         go to 33
   45      npar=nsd+nclass*ndim-2*ndim+nstim*nclass
   33  continue 
         npar= npar+nclass
c   set some variables and memory for m-step
        if(ns.lt.2) then
           nmaxs=nstim
        else
           nclp1=nclass+1
           nmaxs=nstim*nclp1
        end if
          js=istkgt(nmaxs,3)
c
          if(nmodel.ge.2) then
           if(nmodel.lt.4)nup=nstim
           if(nmodel.eq.4)nup=nstim*nclass
           do 87 i=1,nup
            rstak(js+i-1)=10*runif(iseed)
   87         continue
           else
            continue
           end if 
   66     jend=js+nmaxs-1
          do 26 ijq=js,jend
   26          rstak(ijq)=1.0

	abort=.false.
c     
c  EM iterations
        if((nrep.gt.19).and.(icarlo.eq.0))go to 789
               eps=0.001
      call rundvr(rstak(jlam),rstak(jmean),sigma2,rstak(jy),
     1   rstak(jz),nclass,eps,loglik,dstak(jtmp),npar,aic,bic,
     2   test,wword,abort,pgword,indscwd,rstak(jx),rstak(jw),
     4   rstak(js),nmax,nmaxp,nmaxs,
     9   nmaxw,nmaxpw,ndimmax,nrep,icarlo,iseed,rstak(jwin),iconsc) 
c      if(nclass.ge.2)c=(nsub-1-(npar-nclass+1)/nclass-nhi/2)/nsub
c      if(nclass.ge.2)aic=-2.0*c*loglik+3*npar
      write(6,70)
 70   format(//' Final parameter estimates')
      call prtpar(rstak(jlam),rstak(jmean),sigma2,nclass,nspr,loglik,
     1  aic,bic)
 889  write(6,80)
 80    format(//'Posterior Probabilities')
      call prmat(rstak(jz),nsub,nclass)
c       go to 900
       if(nrep.gt.99)go to 789
       write(6,90)
 90    format(//'Final Distance parameter estimates')
       call outprn(nmax,nmaxw,nstim,ndim,nmodel,nclass,
     1     rstak(jx),rstak(jw),
     1     rstak(js),wword)
 789   continue
c release memory from m-step
c  release s
       call istkrl(1)
c

c  release memory for w, z ,tmp
      call istkrl(3)
c  release memory for lambda, mean
 900  continue
      call istkrl(2)
      if(test)go to 999
 500  continue
c release memory for x win wmean
      call istkrl(3)
 999  continue
 800  continue
c   release memory for y
      call istkrl(1)
      stop
      end
c   ------------------------------------------
c  --------------------------------------------
      subroutine conset(sst,rswrd,indscwd,nrep)
c this subroutine sets up convergence criteria
      common/parblk/nstim,ndim,nsub,nwgt,nmodel,ndat,
     1nsd,nspr,ns,nspace
      common/cnvblk/conv,xsconv,wconv,
     1itmax,itxmax,itwmax
      real conv,wconv,xsconv
      logical rswrd,indscwd
      lp=i1mach(2)
      read * ,sst
      if(ns.gt.0)print 1000,sst
 1000 format(1x,5x,"s start=",1x,f5.1)
      read * ,itmax,itxmax,itwmax,conv,wconv
      if(nrep.lt.100) then
      if(rswrd) then
        print 1002
      else
        if(indscwd) then
          print 1003
        else
          print 1001
        end if
      end if
 1003 format(/'     indscal start')
 1001 format(1x,5x,"torgerson start")
 1002 format(1x,5x,"random start")
       write(lp,2000) itmax,itxmax,itwmax,conv,wconv
 2000  format(/1x,5x,'max it for mstp:',i3,3x,'max it for xstp:',i3,3x,
     1     'max it for xstp:',i3,/1x,5x,'mstp conv crit:',f8.5,
     2      3x,'xwstp conv crit:',f8.5)     
      else
        continue
      end if
c      xsconv=float(ndat)*(exp(conv/float(ndat))**2)
      xsconv=wconv
      return
      end
c ----------------------------------------
c
      subroutine xstart(nstim2,nspr2,nmax,x,d,ahold,bhold,chold,
     1he,z,u,aahold,bbhold,pgword)
c  this subroutine does classical mds to start x
      logical pgword,matu,matv
      common/parblk/nstim,ndim,nsub,nwgt,nmodel,ndat,
     1nsd,nspr,ns,nspace
      double precision he(nstim2,nstim2),z(nstim2,nstim2),
     1u(nstim2,nstim2)
      real x(nmax),d(nsub,nspr),aahold(nspr2),bbhold(nspr2)
      double precision cc,cmax,ahmax,dsum,
     1ahold(nstim2),bhold(nstim2),chold(nstim2)
      lp=i1mach(2)
      xn=float(nstim)
      xm=float(ndat)
		call vecclr(nspr,aahold)
		if(nsub.eq.1) go to 45
c compute mean of squared dissimilarities
      dmean=0
      dkt=0
      do 5 m=1,nspr
       sum=0
	icount=0
	do 10 l=1,nsub
      dt=d(l,m)
       if(dt.le.1e-20)go to 10
       sum=sum+dt*dt
       icount=icount+1
       dmean=dmean+dt*dt
       dkt=dkt+1.0
   10 continue
      aahold(m)=sum
    5 bbhold(m)=float(icount)
      dmean=dmean/dkt
      do 20 m=1,nspr
       if(bbhold(m).eq.0.0)go to 30
       aahold(m)=aahold(m)/bbhold(m)
       go to 20
   30  aahold(m)=dmean
   20 continue
   45      do 40 i=1,nstim
   40 he(i,i)=0.d0
      m=0
      do 50 i=2,nstim
       im=i-1
       do 55 j=1,im
        m=m+1
		if(nsub.eq.1)aahold(m)=d(1,m)*d(1,m)
        he(i,j)=dble(-0.5*aahold(m))
   55 he(j,i)=he(i,j)
   50 continue
	cmax=0.d0
c	do 61 i=2,nstim
c		im=i-1
c		do 62 j=1,im
c			do 63 l=1,im
c				cc=dsqrt(2.d0*he(i,j))-
c     1		dsqrt(2.d0*he(i,l))-dsqrt(2.d0*he(j,l))
c   63			if(cc.gt.cmax)cmax=cc
c   62		continue
c   61	continue
c        if(pgword)write(lp,4327)
        if(pgword)print 4327
 4327   format(//1x,5x,'Results of Torgerson Start')
c	if(pgword)print 6667,cmax
c 6667 format(11x,21hadditive constant is ,f10.3)
	do 77 i=1,nstim
		do 78 j=1,nstim
   78		he(i,j)=-0.5d0*((dsqrt(-2.d0*he(i,j))+cmax)**2)
   77	continue
	      dsum=0.d0
      call vecclr(nstim,chold)
      do 60 i=1,nstim
        do 70 j=1,nstim
   70   chold(i)=chold(i)+he(i,j)/dble(xn)
   60 dsum=dsum+chold(i)/dble(xn )
c  compute scaler product matrix
      do 85 i=1,nstim
       do 80 j=1,nstim
   80 he(i,j)=1.d0*(-chold(i)-chold(j)+dsum+he(i,j))
   85	continue
c   do eigenanalysis
 1069 format(1x,12(1x,f4.1))
      call vecclr(nstim,ahold)
      call vecclr(nstim,bhold)
      call vecclr(nstim,chold)
      matu=.false.
      matv=.true.
       call svd(nstim2,nstim,nstim,he,ahold,matu,u,matv,z,ierr,bhold)
      if(pgword)print 2010
 2010 format(1x,5x,"eigenvalues of t-soln")
 2011 format(1x,5x,5(1x,f10.3))
      do 90 iq=1,ndim
	ahmax=0.d0
	do 96 iqq=1,nstim
		if(bhold(iqq).eq.1.0)go to 96
		if(ahmax.lt.ahold(iqq)) then
			ahmax=ahold(iqq)
			index=iqq
		else
			continue
		end if
   96	continue
	bhold(index)=1.d0
             do 95 i=1,nstim
       indx=(iq-1)*nstim+i
   95 x(indx)=z(i,index)*sqrt(abs(ahold(index)))
   90 continue
      if(pgword) then
        print 4000
 4000 format(1x,5x,"x")
      do 100 i=1,ndim
      ist=(i-1)*nstim+1
      iend=(i*nstim)
  100 print 2001,(x(ij),ij=ist,iend)
       else
        continue
       end if
 2001 format(1x,5(2x,f10.2))
      return
      end
c  ----------------------------------------
c   -------------------------------------------
      subroutine center(x,nstim,ndim)
c   this subroutine puts origin of configuration x at centroid
      real x(1)
      mj=0
      do 10 j=1,ndim
        sum=0
        mi=mj
        do 20 i=1,nstim
         mi=mi+1
   20    sum=sum+x(mi)/float(nstim)
        mi=mj
        do 30 i=1,nstim
         mi=mi+1
   30   x(mi)=x(mi)-sum
   10 mj=mj+nstim
      return
      end
c ------------------------------------
c
      subroutine svd(nm,m,n,a,w,matu,u,matv,v,ierr,rv1)
c
      integer i,j,k,l,m,n,ii,i1,kk,k1,ll,l1,mn,nm,its,ierr
      double precision a(nm,n),w(n),u(nm,n),v(nm,n),rv1(n)
      double precision c,f,g,h,s,x,y,z,tst1,tst2,scale,pythag
      logical matu,matv
c
c     this subroutine is a translation of the algol procedure svd,
c     num. math. 14, 403-420(1970) by golub and reinsch.
c     handbook for auto. comp., vol ii-linear algebra, 134-151(1971).
c
c     this subroutine determines the singular value decomposition
c          t
c     a=usv  of a real m by n rectangular matrix.  householder
c     bidiagonalization and a variant of the qr algorithm are used.
c
c     on input
c
c        nm must be set to the row dimension of two-dimensional
c          array parameters as declared in the calling program
c          dimension statement.  note that nm must be at least
c          as large as the maximum of m and n.
c
c        m is the number of rows of a (and u).
c
c        n is the number of columns of a (and u) and the order of v.
c
c        a contains the rectangular input matrix to be decomposed.
c
c        matu should be set to .true. if the u matrix in the
c          decomposition is desired, and to .false. otherwise.
c
c        matv should be set to .true. if the v matrix in the
c          decomposition is desired, and to .false. otherwise.
c
c     on output
c
c        a is unaltered (unless overwritten by u or v).
c
c        w contains the n (non-negative) singular values of a (the
c          diagonal elements of s).  they are unordered.  if an
c          error exit is made, the singular values should be correct
c          for indices ierr+1,ierr+2,...,n.
c
c        u contains the matrix u (orthogonal column vectors) of the
c          decomposition if matu has been set to .true.  otherwise
c          u is used as a temporary array.  u may coincide with a.
c          if an error exit is made, the columns of u corresponding
c          to indices of correct singular values should be correct.
c
c        v contains the matrix v (orthogonal) of the decomposition if
c          matv has been set to .true.  otherwise v is not referenced.
c          v may also coincide with a if u is not needed.  if an error
c          exit is made, the columns of v corresponding to indices of
c          correct singular values should be correct.
c
c        ierr is set to
c          zero       for normal return,
c          k          if the k-th singular value has not been
c                     determined after 30 iterations.
c
c        rv1 is a temporary storage array.
c
c     calls pythag for  dsqrt(a*a + b*b) .
c
c     questions and comments should be directed to burton s. garbow,
c     mathematics and computer science div, argonne national laboratory
c
c     this version dated august 1983.
c
c     ------------------------------------------------------------------
c
      ierr = 0
c
      do 100 i = 1, m
c
         do 100 j = 1, n
            u(i,j) = a(i,j)
  100 continue
c     .......... householder reduction to bidiagonal form ..........
      g = 0.0d0
      scale = 0.0d0
      x = 0.0d0
c
      do 300 i = 1, n
         l = i + 1
         rv1(i) = scale * g
         g = 0.0d0
         s = 0.0d0
         scale = 0.0d0
         if (i .gt. m) go to 210
c
         do 120 k = i, m
  120    scale = scale + dabs(u(k,i))
c
         if (scale .eq. 0.0d0) go to 210
c
         do 130 k = i, m
            u(k,i) = u(k,i) / scale
            s = s + u(k,i)**2
  130    continue
c
         f = u(i,i)
         g = -dsign(dsqrt(s),f)
         h = f * g - s
         u(i,i) = f - g
         if (i .eq. n) go to 190
c
         do 150 j = l, n
            s = 0.0d0
c
            do 140 k = i, m
  140       s = s + u(k,i) * u(k,j)
c
            f = s / h
c
            do 150 k = i, m
               u(k,j) = u(k,j) + f * u(k,i)
  150    continue
c
  190    do 200 k = i, m
  200    u(k,i) = scale * u(k,i)
c
  210    w(i) = scale * g
         g = 0.0d0
         s = 0.0d0
         scale = 0.0d0
         if (i .gt. m .or. i .eq. n) go to 290
c
         do 220 k = l, n
  220    scale = scale + dabs(u(i,k))
c
         if (scale .eq. 0.0d0) go to 290
c
         do 230 k = l, n
            u(i,k) = u(i,k) / scale
            s = s + u(i,k)**2
  230    continue
c
         f = u(i,l)
         g = -dsign(dsqrt(s),f)
         h = f * g - s
         u(i,l) = f - g
c
         do 240 k = l, n
  240    rv1(k) = u(i,k) / h
c
         if (i .eq. m) go to 270
c
         do 260 j = l, m
            s = 0.0d0
c
            do 250 k = l, n
  250       s = s + u(j,k) * u(i,k)
c
            do 260 k = l, n
               u(j,k) = u(j,k) + s * rv1(k)
  260    continue
c
  270    do 280 k = l, n
  280    u(i,k) = scale * u(i,k)
c
  290    x = dmax1(x,dabs(w(i))+dabs(rv1(i)))
  300 continue
c     .......... accumulation of right-hand transformations ..........
      if (.not. matv) go to 410
c     .......... for i=n step -1 until 1 do -- ..........
      do 400 ii = 1, n
         i = n + 1 - ii
         if (i .eq. n) go to 390
         if (g .eq. 0.0d0) go to 360
c
         do 320 j = l, n
c     .......... double division avoids possible underflow ..........
  320    v(j,i) = (u(i,j) / u(i,l)) / g
c
         do 350 j = l, n
            s = 0.0d0
c
            do 340 k = l, n
  340       s = s + u(i,k) * v(k,j)
c
            do 350 k = l, n
               v(k,j) = v(k,j) + s * v(k,i)
  350    continue
c
  360    do 380 j = l, n
            v(i,j) = 0.0d0
            v(j,i) = 0.0d0
  380    continue
c
  390    v(i,i) = 1.0d0
         g = rv1(i)
         l = i
  400 continue
c     .......... accumulation of left-hand transformations ..........
  410 if (.not. matu) go to 510
c     ..........for i=min(m,n) step -1 until 1 do -- ..........
      mn = n
      if (m .lt. n) mn = m
c
      do 500 ii = 1, mn
         i = mn + 1 - ii
         l = i + 1
         g = w(i)
         if (i .eq. n) go to 430
c
         do 420 j = l, n
  420    u(i,j) = 0.0d0
c
  430    if (g .eq. 0.0d0) go to 475
         if (i .eq. mn) go to 460
c
         do 450 j = l, n
            s = 0.0d0
c
            do 440 k = l, m
  440       s = s + u(k,i) * u(k,j)
c     .......... double division avoids possible underflow ..........
            f = (s / u(i,i)) / g
c
            do 450 k = i, m
               u(k,j) = u(k,j) + f * u(k,i)
  450    continue
c
  460    do 470 j = i, m
  470    u(j,i) = u(j,i) / g
c
         go to 490
c
  475    do 480 j = i, m
  480    u(j,i) = 0.0d0
c
  490    u(i,i) = u(i,i) + 1.0d0
  500 continue
c     .......... diagonalization of the bidiagonal form ..........
  510 tst1 = x
c     .......... for k=n step -1 until 1 do -- ..........
      do 700 kk = 1, n
         k1 = n - kk
         k = k1 + 1
         its = 0
c     .......... test for splitting.
c                for l=k step -1 until 1 do -- ..........
  520    do 530 ll = 1, k
            l1 = k - ll
            l = l1 + 1
            tst2 = tst1 + dabs(rv1(l))
            if (tst2 .eq. tst1) go to 565
c     .......... rv1(1) is always zero, so there is no exit
c                through the bottom of the loop ..........
            tst2 = tst1 + dabs(w(l1))
            if (tst2 .eq. tst1) go to 540
  530    continue
c     .......... cancellation of rv1(l) if l greater than 1 ..........
  540    c = 0.0d0
         s = 1.0d0
c
         do 560 i = l, k
            f = s * rv1(i)
            rv1(i) = c * rv1(i)
            tst2 = tst1 + dabs(f)
            if (tst2 .eq. tst1) go to 565
            g = w(i)
            h = pythag(f,g)
            w(i) = h
            c = g / h
            s = -f / h
            if (.not. matu) go to 560
c
            do 550 j = 1, m
               y = u(j,l1)
               z = u(j,i)
               u(j,l1) = y * c + z * s
               u(j,i) = -y * s + z * c
  550       continue
c
  560    continue
c     .......... test for convergence ..........
  565    z = w(k)
         if (l .eq. k) go to 650
c     .......... shift from bottom 2 by 2 minor ..........
         if (its .eq. 30) go to 1000
         its = its + 1
         x = w(l)
         y = w(k1)
         g = rv1(k1)
         h = rv1(k)
         f = 0.5d0 * (((g + z) / h) * ((g - z) / y) + y / h - h / y)
         g = pythag(f,1.0d0)
         f = x - (z / x) * z + (h / x) * (y / (f + dsign(g,f)) - h)
c     .......... next qr transformation ..........
         c = 1.0d0
         s = 1.0d0
c
         do 600 i1 = l, k1
            i = i1 + 1
            g = rv1(i)
            y = w(i)
            h = s * g
            g = c * g
            z = pythag(f,h)
            rv1(i1) = z
            c = f / z
            s = h / z
            f = x * c + g * s
            g = -x * s + g * c
            h = y * s
            y = y * c
            if (.not. matv) go to 575
c
            do 570 j = 1, n
               x = v(j,i1)
               z = v(j,i)
               v(j,i1) = x * c + z * s
               v(j,i) = -x * s + z * c
  570       continue
c
  575       z = pythag(f,h)
            w(i1) = z
c     .......... rotation can be arbitrary if z is zero ..........
            if (z .eq. 0.0d0) go to 580
            c = f / z
            s = h / z
  580       f = c * g + s * y
            x = -s * g + c * y
            if (.not. matu) go to 600
c
            do 590 j = 1, m
               y = u(j,i1)
               z = u(j,i)
               u(j,i1) = y * c + z * s
               u(j,i) = -y * s + z * c
  590       continue
c
  600    continue
c
         rv1(l) = 0.0d0
         rv1(k) = f
         w(k) = x
         go to 520
c     .......... convergence ..........
  650    if (z .ge. 0.0d0) go to 700
c     .......... w(k) is made non-negative ..........
         w(k) = -z
         if (.not. matv) go to 700
c
         do 690 j = 1, n
  690    v(j,k) = -v(j,k)
c
  700 continue
c
      go to 1001
c     .......... set error -- no convergence to a
c                singular value after 30 iterations ..........
 1000 ierr = k
 1001 return
      end
      double precision function pythag(a,b)
      double precision a,b
c
c     finds dsqrt(a**2+b**2) without overflow or destructive underflow
c
      double precision p,r,s,t,u
      p = dmax1(dabs(a),dabs(b))
      if (p .eq. 0.0d0) go to 20
      r = (dmin1(dabs(a),dabs(b))/p)**2
   10 continue
         t = 4.0d0 + r
         if (t .eq. 4.0d0) go to 20
         s = r/t
         u = 1.0d0 + 2.0d0*s
         p = u*p
         r = (s/u)**2 * r
      go to 10
   20 pythag = p
      return
      end
c
      subroutine vecclr(n,x)
      double precision x(1)
      do 10 i=1,n
   10 x(i)=0.d0
      return
      end
c  --------------------------------------------
      subroutine vecex(n,istx,incx,isty,incy,x,y,con)
      double precision x(1),y(1),con
      mx=istx
      my=isty
      do 10 i=1,n
      y(my)=x(mx)*con
      mx=mx+incx
   10 my=my+incy
      return                                               
      end
c  -----------------------------
      subroutine veccon(n,x,con)
      double precision x(1),con
      do 10 i=1,n
   10 x(i)=con
      return
      end
c  ----------------------------------------------
      subroutine swap(a1,b1,f1,a0,b0,f0)
      double precision a1,b1,f1,a0,b0,f0
      a1=a0
      b1=b0
      f1=f0
      return
      end
c  -------------------------------------------------
      subroutine vecsc(n,x,fact)
      double precision x(1),fact
      do 10 ii=1,n
   10 x(ii)=fact*x(ii)
      return
      end
c  ----------------------------------
      subroutine hsym(h,a,n,in)
c  this subroutine converts a sym matrix from full to sym storage
c   h ... input sym matrix
c   n ... order of h
c   a ... output in sym storage mode
      double precision a(1),h(in,in)
      icount=0
      do 10 ii=1,n
        do 20 ij=1,ii
        icount=icount+1
   20   a(icount)=h(ij,ii)
   10 continue
      return
      end
      subroutine matvec(n,a,ira,ncola,ista,istx,x)
      double precision a(ira,1),x(1)
      index=0
      do 10 i=1,n
       ix=index+istx
       ia=index+ista
       x(ix)=a(ia,ncola)
   10 index=index+1
      return
      end
c
c
      subroutine vecmat(n,x,istx,ista,ncola,ira,a)
      double precision a(ira,1),x(1)
      index=0
      do 10 i=1,n
        ix=istx+index
        ia=ista+index
        a(ia,ncola)=x(ix)
   10 index=index+1
      return
      end
c .......................
c ........................
c ........................
      subroutine rvecclr(n,x)
      real x(1)
      do 10 i=1,n
   10 x(i)=0.0
      return
      end
c  --------------------------------------------
      subroutine rvecex(n,istx,incx,isty,incy,x,y,con)
      real x(1),y(1),con
      mx=istx
      my=isty
      do 10 i=1,n
      y(my)=x(mx)*con
      mx=mx+incx
   10 my=my+incy
      return                                               
      end
c  -----------------------------
      subroutine rveccon(n,x,con)
      real x(1),con
      do 10 i=1,n
   10 x(i)=con
      return
      end
c  ----------------------------------------------
      subroutine rswap(a1,b1,f1,a0,b0,f0)
      real a1,b1,f1,a0,b0,f0 
      a1=a0
      b1=b0
      f1=f0
      return
      end
c  -------------------------------------------------
      function rbilin(n,x,y)
      real x(1),y(1),sum,rbilin
      sum=0.0
      do 10 i=1,n
   10 sum=sum+x(i)*y(i)
      rbilin=sum
      return
      end
c   -------------------------------------------------
      subroutine rvecsc(n,x,fact)
      real x(1),fact
      do 10 ii=1,n
   10 x(ii)=fact*x(ii)
      return
      end
c  ----------------------------------
      subroutine rhsym(h,a,n,in)
c  this subroutine converts a sym matrix from full to sym storage
c   h ... input sym matrix
c   n ... order of h
c   a ... output in sym storage mode
c   in ... row and col dim of h in calling program
      real a(1),h(in,in)
      icount=0
      do 10 ii=1,n
        do 20 ij=1,ii
        icount=icount+1
   20   a(icount)=h(ii,ij)
   10 continue
      return
      end
      subroutine rmatvec(n,a,ira,ncola,ista,istx,x)
      real a(ira,1),x(1)
      index=0
      do 10 i=1,n
       ix=index+istx
       ia=index+ista
       x(ix)=a(ia,ncola)
   10 index=index+1
      return
      end
c
c
      subroutine rvecmat(n,x,istx,ista,ncola,ira,a)
      real a(ira,1),x(1)
      index=0
      do 10 i=1,n
        ix=istx+index
        ia=ista+index
        a(ia,ncola)=x(ix)
   10 index=index+1
      return
      end
c .......................
c ........................
c ........................
	subroutine dcomp(wword,i,j,lc,m,x,w,s,dstar,difsav,nclass)
* this subroutine finds the stimulus pair and subject no.
* given a datum and calculates dstar.
	real w(1),s(1),x(1)
	common/parblk/nstim,ndim,nsub,nwgt,nmodel,ndat,
     *	   nsd,nspr,ns,nspace
        logical wword
        real dstar,difsav(1),ssqs,sw(20)
	i = 0
	j = 0
	dstar = 0.0
	call subone(m,i,j)
	do 100 k = 1,ndim
	   indxi = (i + nstim*(k - 1))
	   indxj = (j + nstim*(k - 1))
	   indw = (lc + nclass*(k - 1))
c            print 8001,indxi,indxj,x(indxi),x(indxj)
c 8001       format(11x,3hijx,2x,i3,2x,i3,2x,f6.3,2x,f6.3)
	   difsav(k) = x(indxi) - x(indxj)
           if(.not.wword)go to 100
	   sw(k) = w(indw)
  100	continue
	   indsi = (i + nstim*(lc - 1))
	   indsj = (j + nstim*(lc - 1))
        go to (20,20,30,40,50),nmodel+1
   50     ssqs= s(indsi)+s(indsj)
          go to 20
   40 ssqs = w(nclass*ndim+lc)*(s(i)+s(j))
      go to 20
   30 ssqs = s(i) + s(j)
   20 if(nmodel.eq.0)call submd0(dstar,difsav)
	if(nmodel.eq.1) call submd1(sw,dstar,difsav)
	if(nmodel.eq.2) call submd2(ssqs,dstar,difsav)
	if(nmodel.ge.3) call submd3(sw,ssqs,dstar,difsav)
	return
	end

	subroutine subone(m1,i1,j1)
* subroutine finds stimulus pair and subject #.
      common/parblk/nstim,ndim,nsub,nwgt,nmodel,ndat,
     1nsd,nspr,ns,nspace
	integer m1,i1,j1,l1
        if(m1.gt.nspr) then
        l1=mod(m1,nspr)
        l1=m1/nspr
          if(l1.eq.0) then
              l1=m1/nspr
          else
              l1=m1/nspr+1
          end if
         else
           l1=1
         end if
         l2=(l1-1)*nspr
	m2 = m1 -l2
         do 100 k = 2,ndat
	   kk = k + 1
	npr = k*(k - 1)/2
	npr1 = kk*(kk - 1)/2
	if(m2.eq.npr) then
	i1 = k
	j1 = i1 - 1
        return
	endif
	if(m2.gt.npr) then
	  if(m2.lt.npr1) then
	  i1 = kk
	  j1 = m2 - npr
          return
	  endif
	endif
  100	continue
	return
	end

	subroutine submd0(dstar,difsav)
* subroutine evaluates dstar for model = 0, when all  w's = 1
* and all s's = 0
	common/parblk/nstim,ndim,nsub,nwgt,nmodel,ndat,
     *	   nsd,nspr,ns,nspace
	integer ndim
	real summ,dstar,difsav(1)
	summ = 0.0
	n = ndim
	do 100 k = 1,n
	   summ = summ + difsav(k)**2
  100	continue
	dstar = sqrt(summ)
	return
	end

	subroutine submd1(sqw,dstar,difsav)
* subroutine evaluates dstar for model = 1, when all s's = 0
	common/parblk/nstim,ndim,nsub,nwgt,nmodel,ndat,
     *	   nsd,nspr,ns,nspace

	integer ndim
	real summ,dstar,sqw(1),difsav(1)
	summ = 0.0
	n = ndim
	do 100 k = 1,n
	if(sqw(k).eq.0)sqw(k)=0.00001
	   summ = summ + (sqw(k)*difsav(k)**2)
  100	continue
	dstar = sqrt(summ)
	return
	end

	subroutine submd2(ssqs,dstar,difsav)
* subroutine evalutes dstar for model = 2, when  all w's = 1
* and s's independent of l1
	common/parblk/nstim,ndim,nsub,nwgt,nmodel,ndat,
     *	   nsd,nspr,ns,nspace
	integer ndim
	real summ,dstar,difsav(1),ssqs
	summ = 0.0
	n = ndim
	do 100 k = 1,n
	   summ = summ + difsav(k)**2
  100	continue
      dstar= ssqs+ summ
      if(dstar.lt.0.0)print 1000,ssqs,summ
 1000 format(1x,5x,"distance neg s+s,d",2(1x,f20.10))
      dstar=sqrt(dstar)
	return
	end

	subroutine submd3(sqw,ssqs,dstar,difsav)
*subroutine evaluates dstar for model = 3,or m,odel 4 full model
	common/parblk/nstim,ndim,nsub,nwgt,nmodel,ndat,
     *	   nsd,nspr,ns,nspace
	integer ndim
c	real sqw(1),summ,dstar,difsav(1)
        real sqw(1),difsav(1)
	summ = 0.0
	n = ndim
	do 100 k = 1,n
	   summ = summ + sqw(k) * difsav(k)**2
  100	continue
	dstar = (ssqs + summ)
        dstar=sqrt(dstar)
	return
	end
        subroutine gtdata(y,nsub,nstim)
c
        real y(nsub,nstim)
        character*80 fmt
c
        in=i1mach(1)
        lp=i1mach(2)
c
        read(in,10) fmt
 10     format(a80)
        write(lp,20) fmt
 20     format(/' Data format: ',a80)
c
        do 30 i=1,nsub
 30        read(in,fmt)(y(i,j),j=1,nstim)
c
        write(lp,40)
 40     format(//' Input data')
        call prmat(y,nsub,nstim)
        return
        end
        subroutine gtvdata(x,nstim,ndim,nmax)
        real x(nmax)
        in=i1mach(1)
        lp=i1mach(2)
        write(lp,40)
 40     format(//' Input data x')
          np=nstim*ndim
          read *,(x(l),l=1,np)
          print *,(x(l),l=1,np)
 10       continue
        return
        end
        subroutine wstart(nmodel,lambda,win,w,nsub,ndim,nclass,
     1      wmean,mean,nspr,nmaxw,y,sig)
        common /cstak/ dstak
        real lambda(nclass),wmean(nclass,ndim),win(nsub,ndim)
        real mean(nclass,nspr),y(nsub,nspr),w(nmaxw)
        double precision dstak(300000000),sum
        real rstak(600000000)
        integer istak(600000000) 
        equivalence(dstak(1),rstak(1))
        equivalence(dstak(1),istak(1))
         jmem=istkgt(2*nsub,2)
        call ratini(lambda,wmean,win,nsub,ndim,nclass,istak(jmem))
        do 10 k=1,ndim
          do 20 lc=1,nclass
            indw = (lc + nclass*(k-1))
            w(indw)= wmean(lc,k)
 20       continue
 10       continue
        if(nmodel.eq.3) then
           do 30 lc=1,nclass
 30           w(lc+ndim*nclass)=1.0
         else
           continue
         end if
c compute initial estimates of means and sigma
         do 50 j=1,nspr
         do 40 k=1,nclass
          mean(k,j)=0.0
          do 45 i=1,nsub
           if(istak(jmem-1+i).eq.k) then
            mean(k,j)= mean(k,j)+y(i,j)
           else 
            continue
           end if 
 45        continue
           mean(k,j)=mean(k,j)/(lambda(k)*float(nsub))
 40        continue
 50        continue
         sum=0.d0
         do 80 i=1,nsub
         k=istak(jmem-1+i)
         do 70 j=1,nspr
 70         sum=sum+dble(y(i,j)-mean(k,j))**2
 80         continue
         sig=sngl(sum/dfloat(nspr*nsub))
         call istkrl(1)
         return
         end
      subroutine gtpar(nsub,nstim,nspr,ndimlo,ndimhi,ns,
     1    nlo,nhi,istart,iseed,nrep,icarlo,iconsc)
c
      character*80 title
c
      in=i1mach(1)
      lp=i1mach(2)
c
      read(in,10) title
 10   format(a80)
      write(lp,20) title
 20   format(' Latent Class Analysis of Distance Data'//
     1       ' Version 7.01'/
     1       ' ------------------------------------'//
     2       ' Job title: ',a80)
c
      read(in,30) nsub,nstim,nlo,nhi,ns,ndimlo,ndimhi,
     1    istart,iseed,nrep,icarlo,iconsc
 30   format(12i5)
      if (nstim .lt.10) call seterr
     1   ('gtpar - too few stimuli',23,1,2)
      if (nsub .lt. 15) call seterr
     1   ('gtpar - too few subjects',24,2,2)
      if (nlo .lt. 1 .or. nhi .gt. nsub) call seterr
     1   ('gtpar - invalid no. of classes',30,3,2)
c      if (nhi .gt. nsub/8) call seterr
c     1   ('gtpar - invalid no. of classes',30,3,2)
      if(ndimlo.lt.1) call seterr
     1   ('gtpar - too few dimensions',26,4,2) 
      if (istart .lt. 0) istart=0
      if(iseed.lt.1)iseed=326
c
      if((nrep.lt.20).or.(icarlo.eq.1)) then
      write(lp,40) nsub,nstim,nlo,nhi,ndimlo,ndimhi,
     1    iseed,ns
 40   format(/' No. of subjects:      ',i5/
     1        ' No. of stimuli:       ',i5/
     2        ' Lower No. of latent classes:',i5/
     2        ' Upper No. of Latent classes:',i5/    
     3        ' Lower No. of dimensions:    ',i5/
     3        ' Upper No. of dimensions:    ',i5/
     4        ' iseed =',i5/
     5        ' specificities parameter:    ',i5)
      else
      write(lp,50) nsub,nstim,nlo,nhi,
     1    iseed,nrep
 50   format(/' No. of subjects:      ',i5/
     1        ' No. of stimuli:       ',i5/
     2        ' Lower No. of latent classes:',i5/
     2        ' Upper No. of Latent classes:',i5/    
     4        ' iseed =',i5/
     5    ' no. of replicatins for monte carlo test on null model:',i5)
       end if
      if (istart .eq. 0) write(lp,80)
 80   format(' Rational initial parameter estimates')
      if (istart .gt. 0) write(lp,90) istart
 90   format(' Random initial parameter estimates, seed:',i6)
c
      return
      end
       subroutine gtpmod(ns,ndim,nstim,nclass,nspr,nmax,nmaxw,
     1    x,w,s,dtmp)
       real x(nmax),w(nmaxw),dtmp(nclass,nspr),s(1)
       do 10 nc=1,nclass
           m=0
           do 20 i=2,nstim
             im1=i-1
             do 30 j=1,im1
               m=m+1
               dm=0.0
               do 40 ndimc=1,ndim 
                inxi=(ndimc-1)*nstim+i
                inxj=(ndimc-1)*nstim+j
                indw=(ndimc-1)*nclass+nc
                indsi=(nc-1)*nstim+i
                indsj=(nc-1)*nstim+j
                dm=dm+w(indw)*((x(inxi)-x(inxj))**2)
   40             continue
                if(ns.gt.0) then
                 if(nclass.gt.1) then
                   if(ns.eq.1)dm=dm+w(ndim*nclass+nc)*(s(i)+s(j))
                   if(ns.eq.2)dm=dm+s(indsi)+s(indsj)
                 else
                   dm=dm+s(i)+s(j)
                 end if
                else
                  continue
                end if 
                dm=sqrt(dm)
                dtmp(nc,m)=dm
   30             continue
   20               continue
   10                 continue
             return
             end
       subroutine gtpmodr(ns,ndim,nstim,nclass,nspr,nmax,nmaxw,
     1    x,w,s,dtmp,lambda)
       real x(nmax),w(nmaxw),dtmp(nclass,nspr),s(1),lambda(nclass)
       nsd=nstim*ndim
       do 2 i=1,nsd
 2        x(i)=runif(iseed)
       do 3 i=1,nstim
 3        s(i)=runif(iseed)
       nw=ndim*nclass
       ncm1=nclass-1
       in=0
       do 6 j=1,ndim
       wc=0
       do 4 i=1,ncm1
         w(i+in)=runif(iseed)
 4       wc=wc+w(i+in)
       w(in+nclass)=nclass-wc
 6     in=in+nclass
       cuml=0.0
       do 5 i=1,ncm1
          lambda(i)=runif(iseed)
 5        cuml=cuml+lambda(i)
        lambda(nclass)=1.0-cuml
       do 10 nc=1,nclass
           m=0
           do 20 i=2,nstim
             im1=i-1
             do 30 j=1,im1
               m=m+1
               dm=0.0
               do 40 ndimc=1,ndim 
                inxi=(ndimc-1)*nstim+i
                inxj=(ndimc-1)*nstim+j
                indw=(ndimc-1)*nclass+nc
                indsi=(nc-1)*nstim+i
                indsj=(nc-1)*nstim+j
                dm=dm+w(indw)*((x(inxi)-x(inxj))**2)
   40             continue
                if(ns.gt.0) then
                 if(nclass.gt.1) then
                   if(ns.eq.1)dm=dm+w(ndim*nclass+nc)*(s(i)+s(j))
                   if(ns.eq.2)dm=dm+s(indsi)+s(indsj)
                 else
                   dm=dm+s(i)+s(j)
                 end if
                else
                  continue
                end if 
                dm=sqrt(dm)
                dtmp(nc,m)=dm
   30             continue
   20               continue
   10                 continue
             return
             end
C
C
C
        INTEGER FUNCTION IGETIJ(I,J)
C
C  RETURNS THE INDEX OF ELEMENT (I,J) OF A SYMMETRIC MATRIX WHOSE
C  LOWER TRIANGULAR PART IS STORED ROW-WISE IN A LINEAR ARRAY
C
        IF(I.GT.J) GOTO 10
        IGETIJ=(J-2)*(J-1)/2+I
        RETURN
 10     IGETIJ=(I-2)*(I-1)/2+J
        RETURN
        END
      subroutine logl(lambda,mean,sigma2,y,z,
     1   nclass,nspr,nsub,temp,loglik,loglikp,npar,aic,bic)
c
c  compute total log-likelihood and a posteriori membership probabilities
c
c
      real lambda(nclass),mean(nclass,nspr),y(nsub,nspr)
      real z(nsub,nclass)
      double precision pi,temp(nclass),loglik,loglikp
      double precision lli,tmp,twosig
      save pi
c
      data pi / 3.14159265358979323846d0 /
c
      twosig=2.d0*dble(sigma2)
      loglik=0.d0
      do 40 i=1,nsub
        index=0
        lli=0.d0
        do 20 k=1,nclass
           tmp=0.d0
           do 10 j=1,nspr
 10           tmp=tmp+ dble(y(i,j)-mean(k,j))**2
           if(twosig.eq.0.d0)print 1717
 1717      format(//"twosig is zero")
           tmp=tmp/twosig
           if (tmp .gt.1500) then
              tmp=0.d0
              if (lambda(k) .gt. 10.e-30) lli=lli+dble(lambda(k))*tmp
           else
              tmp=dexp(-tmp)
              if (lambda(k) .gt. 10.e-30) lli=lli+dble(lambda(k))*tmp
           endif
              temp(k)=tmp
 20         continue
c
c
c
c$$$ 1979       format(10x,2hi=,i2,2x,4hlli=,g20.10)
c$$$            print 1979,i,lli
c
c
c
        if(lli.le.10.e-30)go to 97
        loglik=loglik+dlog(lli)
   97        do 30 k=1,nclass
                if (temp(k) .le. 10.d-60 .or. lambda(k) .le. 10.e-30) 
     1               then
                   z(i,k)=0.0
                   if(nclass.eq.1)z(i,k)=1.0
                else
                   if(lli.eq.0.d0)print 1718
 1718              format(//"lli is zero")
                   z(i,k)=sngl(dble(lambda(k))*temp(k)/lli)
                endif
           if(z(i,k).ne.0.0)index=1
 30         continue
          if(index.eq.0) then
            do 35 k=1,nclass
   35        z(i,k)=lambda(k)
          else
            go to 40
          end if
 40      continue
c
c
c
c$$$ 1980    format(10x,5hnsub=,i2,2x,5hnspr=,i3,2x,7htwosig=,g20.10)
c$$$         print 1980,nsub,nspr,twosig
c
c
c
c
      loglik = loglik - dfloat(nsub*nspr)*dlog(twosig*pi)/2.d0
      loglikp=loglik-dfloat(npar)
      aic=-2.0*sngl(loglik)+2.0*npar
      bic=-2.0*sngl(loglik)+log(float(nspr*nsub))*npar
c
      return
      end
      subroutine zinit(nsub,nspr,z,y,lambda,mean)
      dimension z(nsub,nsub),y(nsub,nspr),mean(nsub,nspr)
      dimension lambda(nsub) 
      do 10 i=1,nsub
      do 20 j=1,nsub
        if(i.eq.j) then
           z(i,j)=1.0
        else
           z(i,j)=0.0
        end if
 20      continue
 10       continue
       do 30 i=1,nsub
 30       lambda(i)=1.0/float(nsub)
       do 40 i=1,nsub
       do 50 j=1,nspr
        mean(nsub,nspr)=y(nsub,nspr)
 50     continue
 40     continue
      return
      end
      subroutine dimstp(irklim,ssn,x,w,s,ft,mean1,
     1g,trial,abort,
     2wword,pgword,he,a1,b1,c1,dif,yy,xs,ahe,nmax,
     3inf,ghat,sp,hld,mean,nclass,hew,nmaxw)    
c   this subroutine updates the dimensions
      common/parblk/nstim,ndim,nsub,nwgt,nmodel,ndat,
     1nsd,nspr,ns,nspace
      common/cnvblk/conv,xsconv,wconv,
     1itmax,itxmax,itwmax
      logical wword,abort,pgword,xwrd
      double precision sso,ssn,slim,gnorm,eps,gnormx,
     1grmin    
      real mean1(nclass,nspr),mean(nclass,nspr),ft(nclass)
      real x(nmax),w(1),s(1)
      real dif(1),yy(1),xs(nmax)
      double precision a1(1),b1(1),c1(1),
     1he(nmax,nmax),ahe(1),g(1),hew(nmaxw,nmaxw)
      double precision ghat(1),sp(nmax,nmax),hld(nmax,nmax)
      integer inf(1)
c  needed for gradient check
c      double precision su(200),f0,fv(200),gv(200),del
c      real xu(100),ssu(100)
c      data del/0.01d0/
c  end of needed for gradient check
      data eps/0.000001d0/
      slim=0.0d0
c   set up xs
c       nw=nclass*ndim
c       nsd=nstim*ndim
c       print 1718,nsd,ns,nmodel,nmax,nspace
c 1718  format(1x,5x,"nsd,ns,nmodel,nmax,nspace   ",5(i3,1x))
      call rvecex(nsd,1,1,1,1,x,xs,1.0)
c       do 12 ii=1,nsd
c 12       xs(ii)=1.0 *sngl( x(ii))
      if (nmodel.eq.2)nss=nstim
      if(nmodel.eq.3)nss=nstim
      if(nmodel.eq.4)nss=nstim*nclass
      call rvecex(nss,1,1,nsd+1,1,s,xs,1.0)
c   compute initial ss and g and hessian
c      print 1717
c 1717 format(/1x,5x,"in dim")
   30 call dimfgh(x,w,s,ft,mean1,sso,g,abort,
     1wword,pgword,he,nmax,dif,nclass,mean)
c      print 1717
c   ****************
c   *****************
c   gradient checker dble length
c   nm dim of x;x is base test pt; diff finite diff interval;
c   s test quantities ; h numerical hessian;
c       f0=sso
c       do 5 j=1,nspace
c        do 3 l=1,nsd
c    3     xu(l)=x(l)
c        if(nmodel.gt.2) then
c          do 7 l=1,nstim
c    7           ssu(l)=s(l)
c        else
c         continue
c        end if
c       if(j.le.nsd)xu(j)=xu(j)+ sngl(del)
c       if(j.gt.nsd)ssu(j-nsd)=ssu(j-nsd)+ sngl(del)
c       call dimfgh(xu,w,ssu,ft,mean1,fv(j),gv,abort,
c     1wword,pgword,he,nmax,dif,nclass,mean)
c       do 4 i=1,nspace
c    4    he(i,j)=(gv(i)-g(i))/del
c         su(j)=he(j,j) + (1.d0/del)*(g(j)-(fv(j)-f0)/del)
c    5 continue
c      print 1050
c      do 6 i=1,nspace
c    6 print 1051,su(i),g(i)     
c 1050 format(11x,14hgradient check,"in dimstp")
c 1051 format(5x,2hs=,f20.5,5x,2hg=,f20.5)
c   *******************
c   *******************
      if(abort)go to 98
      itxs=0
       do 814 ii=1,nspace
  814      a1(ii)=dble(xs(ii))
       gnorm=0.d0
       gnormx=0.d0
c       print 1717
       do 821 ii=1,nspace
         gnorm=gnorm+g(ii)*g(ii)
         gnormx=gnormx+a1(ii)*a1(ii)
  821    continue
       if(gnorm.ge.0.d0)gnorm=dsqrt(gnorm)
       if(gnormx.ge.0.d0)gnormx=dsqrt(gnormx)
       gnormx=gnorm*gnormx/dble(nspace)
      if(pgword)print 1000,itxs,sso
      if(pgword)print 1001,gnorm,gnormx
 1001 format(1x,11x,5hgnorm,1x,f15.5,2x,6hgnormx,1x,f15.5)
 1000 format(1x,5x,6hitx.no,1x,i3,2x,2hss,f10.3)
c   loop for xs teration
      iccount=0
c      print 1717
      do 100 itxs=1,itxmax
c  compute initial reduced gradient and hessian store pos coeff in free
  101     continue
       iccount=iccount+1
       if(iccount.gt.7)go to 100
       if(itxs.gt.1)
     1call dimfgh(x,w,s,ft,mean1,sso,g,abort,wword,pgword,
     2he,nmax,dif,nclass,mean)    
   50     if(nmodel.ge.2) then
        do 201 ii=1,nspace
  201    inf(ii)=1
        do 200 in=1,nss
          ind=in+nsd
  200   if(xs(ind).le.sngl(slim+0.005d0))inf(nsd+in)=0
        if(pgword.and.(nmodel.ge.2))
     1     print 1816,(inf(ip),ip=1,nspace)
 1816 format(1x,5x,"ifree",36(1x,i2))
        nfrees=nsd
        do 202 in=1,nss
  202   nfrees=nfrees+inf(in+nsd)
      if(nfrees.eq.nspace)go to 111
        grmin=0.0
        index1=0
        do 116 in=1,nss
          if(grmin.gt.g(in).and.inf(in+nsd).eq.0) then
              grmin=g(in)
              index1=in
          else
                 go to 116
          end if
  116   continue
      if(grmin.ge.-1d-3)go to 111
      if(index1.eq.0)go to 111
      inf(index1)=1
  111  continue
        if(nfrees.eq.nspace) then
          call vecex(nspace,1,1,1,1,g,ghat,1.d0)
        else
          call reduc2(g,ghat,he,inf,
     1nfrees,nspace,sp,hld,nmax)
          gnorm=0.d0
         do 851 ii=1,nfrees
  851        gnorm=gnorm+ghat(ii)*ghat(ii)
          gnorm=dsqrt(gnorm)
          if(pgword)print 1014,gnorm
 1014     format(1x,5x,"reduced grnorm=",1x,f20.10)
        end if
      else
        nfrees=nsd+ns
        call vecex(nspace,1,1,1,1,g,ghat,1.d0)
      end if
c   compute rank limit for hessian
      irklim=nfrees-ndim-(ndim*(ndim-1))/2
	if(wword)irklim=nfrees-2*ndim
	if(wword.and.(nmodel.eq.3))irklim=irklim-1
c   compute search direction
        call vecex(nfrees,1,1,1,1,ghat,b1,-1.d0)
        call hsym(he,ahe,nfrees,nmax)
        call dmfss(ahe,nfrees,eps,irank,c1)
        call dmlss(ahe,nfrees,irank,c1,0,b1,ier)
        if(ier.ne.0)print 2001,ier
 2001 format(11x,16her from dml i1r=,i2)
   76 if(nfrees.eq.nspace)go to 75
      call vecclr(nspace,a1)
      do 45 il=1,nspace
        do 46 li=1,nfrees
   46   a1(il)=a1(il)+sp(li,il)*b1(li)
   45 continue
      call vecex(nspace,1,1,1,1,a1,b1,1.d0)
c   do linesearch
   75      xwrd=.true.
        trial=1.d0
        if(nmodel.lt.2)trial=0.5d0
        call vecclr(nspace,a1)
c
c
c
c$$$ 1979   format(10x,f20.10,1h,)
c$$$        do 1980 ij=1,nstim
c$$$           do 1981 ir=1,ndim
c$$$              print 1979,b1(ij+(ir-1)*nstim)
c$$$ 1981      continue
c$$$ 1980   continue
c$$$        print 1979,(b1(is+nsd),is=1,nstim)
c
c
c
        call sear(nspace,xs,yy,g,c1,b1,sso,ssn,ind,slim,
     1       5,trial,x,w,s,ft,mean1,abort,wword,pgword,xwrd,dif,
     2      mean,nclass,he,nmax,hew,nmaxw)
      if(pgword)print 1607,sso,ssn
 1607 format(1x,"linesearch output,sso,ssn",2(1x,f20.10))
      if(ssn.gt.sso)go to 85
      call vecex(nspace,1,1,1,1,c1,g,1.d0)
      call rvecex(nspace,1,1,1,1,yy,xs,1.0)
       gnorm=0.d0
       gnormx=0.d0
       do 817 ii=1,nspace
  817      a1(ii)=dble(xs(ii))
       do 841 ii=1,nspace
          gnorm=gnorm+g(ii)*g(ii)
          gnormx=gnormx+a1(ii)*a1(ii)
  841      continue
         gnorm=dsqrt(gnorm)
       gnormx=dsqrt(gnormx)
       gnormx=gnorm*gnormx/dble(nspace)
   85 if(pgword)print 1004,itxs,ind,ssn
      if(pgword)print 1001,gnorm,gnormx 
 1004 format(1x,1x,"itx no ",1x,i3,1x,"ind",1x,i1,
     11x,"ss",1x,f10.3)
      if(ssn.le.sso)go to 90
      ssn=sso
      if(pgword)print 1007
 1007 format(1x,5x,21hfun not imp this step)
      go to 100
   90 if(sngl((sso-ssn)/sso).ge.(xsconv))go to 100
      if(nfrees.eq.nspace)go to 110
        grmin=0.d0
        index1=0
          do 115 in=1,nstim
c
c
c
c$$$             print 1979,g(in)
c$$$ 1979        format(10x,f20.10)
c
c
c
          if(grmin.gt.g(in).and.inf(in+nsd).eq.0) then
              grmin=g(in)
              index1=in
          else
                 go to 115
          end if
  115   continue
      if(grmin.ge.-1d-3)go to 110
      if(index1.eq.0)go to 110
      inf(index1)=1
       go to 101
  100 sso=ssn
      go to 99
  110 if(pgword)print 1002
 1002 format(1x,11x,14hconfig converg)
   99 call rvecex(nsd,1,1,1,1,xs,x,1.0)
c   99   continue
      if(nmodel.ge.2)call rvecex(nss,nsd+1,1,1,1,xs,s,1.0)
   98 continue
      return
      end
c  _________________________________
      subroutine wstp(irklim,ssn,x,w,s,
     1ft,mean1,g,trial,abort,
     2wword,pgword,he,a1,b1,c1,dif,yy,xs,ahe,nmax,
     3inf,ghat,sp,hld,mean,nclass,hex,nmaxx)    
c   this subroutine updates the dimensions
      common/parblk/nstim,ndim,nsub,nwgt,nmodel,ndat,
     1nsd,nspr,ns,nspace
      common/cnvblk/conv,xsconv,wconv,
     1itmax,itxmax,itwmax
      real mean1(nclass,nspr),mean(nclass,nspr),ft(nclass)
      logical wword,abort,pgword,xwrd
      double precision sso,ssn,slim,eps,
     1grmin    
      real x(1),w(1),s(1),
     2dif(1),yy(1),xs(1)
      double precision a1(1),b1(1),c1(1),
     1he(nmax,nmax),ahe(1),g(1),hex(nmaxx,nmaxx)
      double precision ghat(1),sp(nmax,nmax),hld(nmax,nmax)
      integer inf(1)
c  needed for gradient check

c      double precision su(200),f0,fv(200),gv(200),del
c      real xu(100)
c      data del/0.001d0/
c  end of needed for gradient check
      data eps/0.000001d0/
      slim=0.0d0
      nw=nclass*ndim
       if(nmodel.eq.3)nw=nw+nclass
c   set up xs
      icount=0
      call rvecex(nw,1,1,1,1,w,xs,1.0)
c   compute initial ss and g and hessian
   30 call wfgh(x,w,s,ft,mean1,sso,g,abort,
     1wword,pgword,he,nmax,dif,nclass)
c   ****************
c   *****************
c   gradient checker dble length
c   nm dim of x; x is base test pt; diff finite diff interval;
c   s test quantities ; h numerical hessian;
c       f0=sso
c       do 5 j=1,nw
c        do 3 l=1,nw
c    3     xu(l)=x(l)
c       xu(j)=xu(j)+ del
c       call wfgh(xu,w,s,ft,mean1,fv(j),gv,abort,
c     1wword,pgword,he,nmax,dif,nclass)
c       do 4 i=1,nw
c    4    he(i,j)=(gv(i)-g(i))/del
c         su(j)=he(j,j) + (1.d0/del)*(g(j)-(fv(j)-f0)/del)
c    5 continue
c      print 1050
c      do 6 i=1,nw
c    6 print 1051,su(i),g(i)     
c 1050 format(11x,14hgradient check,"inwstp")
c 1051 format(5x,2hs=,f20.5,5x,2hg=,f20.5)
c   *******************
c   *******************
      if(abort)go to 98
      itws=0
      if(pgword)print 1000,itws,sso
 1000 format(1x,5x,6hitw.no,1x,i3,2x,2hss,f20.3)
c   loop for xs iteration
      do 100 itws=1,itwmax
c  compute initial reduced gradient and hessian store pos coeff in free
  101  icount=icount+1
        if(icount.ge.5)go to 100  
        if(itws.gt.1)
     1call wfgh(x,w,s,ft,mean1,sso,g,abort,wword,pgword,
     2he,nmax,dif,nclass)    
        do 201 ii=1,nw
  201    inf(ii)=1
        do 200 in=1,nw
          ind=in
  200   if(xs(ind).le.sngl(slim+0.005d0))inf(in)=0
        if(pgword)
     1     print 1816,(inf(ip),ip=1,nw)
 1816 format(1x,5x,"ifree",36(1x,i2))
        nfrees=0
        do 202 in=1,nw
  202   nfrees=nfrees+inf(in)
c        do 445 lq=1,ndim
c         isum=0
c         do 446 nq=1,nclass
c  446        isum=isum+inf(nclass*(lq-1)+nq)
c         if(isum.eq.0) call seterr('in w-step- too many dimensions',
c     1            30,6,2)       
c  445     continue
      if(nfrees.eq.nw)go to 111
        grmin=0.0
        index1=0
        do 116 in=1,nstim
          if(grmin.gt.g(in).and.inf(in).eq.0) then
              grmin=g(in)
              index1=in
          else
                 go to 116
          end if
  116   continue
      if(grmin.ge.-1d-3)go to 111
      inf(index1)=1
  111  continue
        if(nfrees.eq.nw) then
          call vecex(nw,1,1,1,1,g,ghat,1.d0)
        else
          call reduc2(g,ghat,he,inf,
     1nfrees,nw,sp,hld,nmax)
        end if
c   compute rank limit for hessian
      irklim=nfrees
c   compute search direction
        call vecex(nfrees,1,1,1,1,ghat,b1,-1.d0)
        call hsym(he,ahe,nfrees,nmax)
        call dmfss(ahe,nfrees,eps,irank,c1)
        call dmlss(ahe,nfrees,irank,c1,0,b1,ier)
        if(ier.ne.0)print 2001,ier
 2001 format(11x,16her from dml i2r=,i2)
   76 if(nfrees.eq.nw)go to 75
      call vecclr(nw,a1)
      do 45 il=1,nw
        do 46 li=1,nfrees
   46   a1(il)=a1(il)+sp(li,il)*b1(li)
   45 continue
      call vecex(nw,1,1,1,1,a1,b1,1.d0)
c   do linesearch
   75      xwrd=.false.
        trial=1.d0
        if(nmodel.lt.2)trial=0.5d0
        call vecclr(nw,a1)
c
c
c
c$$$ 1979   format(10x,f20.10,1h,)
c$$$        do 1980 ir=1,ndim
c$$$           do 1981 it=1,nclass
c$$$              print 1979,b1((3-it)+(ir-1)*nclass)
c$$$ 1981      continue
c$$$ 1980   continue
c$$$        print 1979,(b1((3-is)+ndim*nclass),is=1,nclass)
c
c
c
        call sear(nw,xs,yy,g,c1,b1,sso,ssn,ind,slim,
     1     5,trial,x,w,s,ft,mean1,abort,wword,pgword,xwrd,dif,
     2      mean,nclass,hex,nmaxx,he,nmax)
      if(pgword)print 1607,sso,ssn
 1607 format(1x,"linesearch output,sso,ssn",2(1x,f20.10))
      if(ssn.gt.sso)go to 85
      call vecex(nw,1,1,1,1,c1,g,1.d0)
      call rvecex(nw,1,1,1,1,yy,xs,1.0)
   85 if(pgword)print 1004,itws,ind,ssn
 1004 format(1x,1x,"itw no ",1x,i3,1x,"ind",1x,i1,
     11x,"ss",1x,f10.3)
      if(ssn.le.sso)go to 90
      ssn=sso
      if(pgword)print 1007
 1007 format(1x,5x,21hfun not imp this step)
      go to 100
   90 if(sngl((sso-ssn)/sso).ge.(wconv))go to 100
      if(nfrees.eq.nw)go to 110
        grmin=0.d0
        index1=0
          do 115 in=1,nstim
          if(grmin.gt.g(in).and.inf(in).eq.0) then
              grmin=g(in)
              index1=in
          else
                 go to 115
          end if
  115   continue
      if(grmin.ge.-1d-3)go to 110
      if(index1.eq.0)go to 110
      inf(index1)=1
       go to 101
  100 sso=ssn
      go to 99
  110 if(pgword)print 1002
 1002 format(1x,11x,14hweight converg)
   99 call rvecex(nw,1,1,1,1,xs,w,1.0)
      call wnorm(x,w,s,nstim,ndim,nclass,nmodel)
      if(pgword)print 1111
 1111 format(1x,"w after wstp")
      if(pgword)print *,(w(ik),ik=1,nw)
   98 continue
      return
      end
c  _________________________________
      subroutine onecl(lambda,mean,sigma2,y,nsub,nstim)
c
c  handle special case of 1 latent class
c
      real lambda(1),mean(1,nstim),y(nsub,nstim)
      double precision sig 
c
      lambda(1)=1.0
c
c  means are just column means of y
c
      do 20 j=1,nstim
        sum=0.0
        do 10 i=1,nsub
 10         sum=sum+y(i,j)
         mean(1,j)=sum/float(nsub)
 20      continue
c
      sig=0.d0
      do 30 i=1,nsub
      do 30 j=1,nstim
 30      sig=sig+(dble(y(i,j))-dble(mean(1,j)))**2
      sig=sig/dfloat(nsub*nstim)
      sigma2=sngl(sig)
c
      return
      end
        SUBROUTINE PRMAT(A,N,M)
C
        REAL A(N,M)
C
        LU=I1MACH(2)
C
        DO 20 K=1,M,10
           L=MIN0(M,K+9)
           WRITE(LU,30)(I,I=K,L)
           WRITE(LU,40)
           DO 10 I=1,N
 10           WRITE(LU,40)I,(A(I,J),J=K,L)
 20        CONTINUE
        RETURN
C
 30     FORMAT(/2X,10I12)
 40     FORMAT(1X,I4,10F12.4)
        END
      subroutine prtpar(lambda,mean,sigma2,nclass,nstim,loglik,
     1    aic,bic)
c
      real lambda(nclass),mean(nclass,nstim),sigma2
      double precision loglik
c
      lp=i1mach(2)
c
      write(lp,10) (lambda(i),i=1,nclass)
 10   format(/' Mixing parameters (lambda):'/(1x,10f12.4))
      write(lp,20)
 20   format(/' Class means (classes x stimuli):')
      call prmat(mean,nclass,nstim)
      write(lp,30) sigma2
 30   format(/' Variance parameter (sigma2): ',e15.7)
      write(lp,40) loglik
 40   format(/' Log-likelihood for these estimates: ',d17.9)
      write(lp,50) aic
 50   format(/' AIC statistic for these estimates:  ',f20.10)
      write(lp,60) bic
 60   format(/' BIC statistic for these estimates:  ',f20.10) 
c
      return
      end
       subroutine dimfgh(x,w,s,ft,mean1,ss,g,
     1abort,wword,pgword,he,nmax,dif,nclass,mean)
c   this is a subroutine  to compute f and g for dimensions
      common/parblk/nstim,ndim,nsub,nwgt,nmodel,ndat,
     1nsd,nspr,ns,nspace     
      real mean1(nclass,nspr),mean(nclass,nspr),ft(nclass)
      double precision ss,res,gim,gim1,
     1gw,gws,g(1),dterm,gwp,gwq,gwj
      double precision he(nmax,nmax),hterm,sterm,sterm1
      logical wword,abort,pgword
      real x(1),w(1),s(1),dif(1)
      call vecclr(nspace,g)
      do 2 il=1,nspace
      do 2 ij=1,nspace
    2 he(il,ij)=0.d0
      ss=0.d0
c
c
c
c      print 1717
c 1717 format(1x,5x,"in fgi")
c
c
c
      do 300 nc=1,nclass
      do 100 m=1,nspr
c   go thru data getting dstar and difsav(xi-xj) for each dimension
      call dcomp(wword,i,j,nc,m,x,w,s,dstar,dif,nclass)
c      if(d(l,m).lt.0.0) then
c	print 6111,i,j,l,d(l,m)
c 6111 format(11x,3hd<0,2x,4hi,j=,i3,i3,2x,2hl=,i3,2x,2hd=,f20.10)
c      else
c	continue
c      end if
c
c
c
c 6111 format(10x,7h(j,k)=(,i2,1h,,i2,1h),2x,2ht=,i1,2x,2hd=,f16.10)
c      print 6111,j-1,i-1,nc-1,dstar
c
c
c	
      dhat=dstar
      if(dstar.eq.0.0)go to 100
      mean(nc,m)=dstar
      res=dble(mean1(nc,m)-dhat)
c
c
c
c 6112 format(10x,7h(j,k)=(,i2,1h,,i2,1h),2x,2ht=,i1,2x,4hres=,f16.10)
c      print 6112,j-1,i-1,nc-1,res
c
c 6113 format(10x,7h(j,k)=(,i2,1h,,i2,1h),2x,2ht=,i1,2x,4hwre=,f16.10)
c      print 6113,j-1,i-1,nc-1,res*dble((1.0/dstar)*ft(nc))
c
c
c
      ss=ss+res*res*dble(ft(nc))
      do 15 im=1,ndim
       imj=j+nstim*(im-1)
       imi=i+nstim*(im-1)
              gim1=dble((1.0/dstar)*ft(nc))
          if(wword) then
              gw=dble(w((im-1)*nclass+nc))
          else
            gw=1.d0
          end if
        gim=-2.d0*res*gim1*gw*dble(dif(im))
        g(imi)=g(imi)+gim
   15  g(imj)=g(imj)-gim
       if(nmodel.lt.2)go to 101
       gws=1.d0
       if(nmodel.eq.3)gws=dble(w(ndim*nclass+nc))
       if(nmodel.le.3) then
         ili=nsd+i
         ilj=nsd+j
       else
         ili=nsd+i+(nc-1)*nstim
         ilj=nsd+j+(nc-1)*nstim
       end if
       g(ilj)=g(ilj)-res*gim1*gws
       g(ili)=g(ili)-res*gim1*gws
  101     dterm=dble((1.0/dstar)**2*ft(nc))
      do 20 iq=1,ndim
        inxi=i+(iq-1)*nstim
        do 30 ip=1,ndim
          inxj=j+(ip-1)*nstim
           if(wword) then
             gwp=dble(w(nclass*(ip-1)+nc))
             gwq=dble(w(nclass*(iq-1)+nc))
           else
             gwp=1.d0
             gwq=1.d0
           end if
      he(inxi,inxj)=he(inxi,inxj)-2.d0*gwp*gwq*dble(dif(ip))
     1*dble(dif(iq))*dterm
      he(inxj,inxi)=he(inxj,inxi)-2.d0*gwp*gwq*dble(dif(ip))
     1*dble(dif(iq))*dterm
   30   continue
   20 continue
      do 50 iq=1,ndim
          do 40 ij=1,ndim
            if(wword) then
               gwq=dble(w(nclass*(iq-1)+nc))
               gwj=dble(w(nclass*(ij-1)+nc))
            else
               gwq=1.d0
               gwj=1.d0
            end if
            inxi1=i+(iq-1)*nstim
            inxi2=i+(ij-1)*nstim
            inxj1=j+(iq-1)*nstim
            inxj2=j+(ij-1)*nstim
            hterm=2.d0*gwq*gwj*dterm*dble(dif(ij)*dif(iq))
            he(inxj1,inxj2)=he(inxj1,inxj2)+hterm
            he(inxi1,inxi2)=he(inxi1,inxi2)+hterm
   40     continue
   50 continue
          if(nmodel.lt.2) go to 100
           if(nmodel.eq.3) then
             gws=dble(w(nclass*ndim+nc))
           else
             gws=1.d0
           end if
          if(nmodel.le.3) then
            insi=nsd+i
            insj=nsd+j
          else
            insi=nsd+i+(nc-1)*nstim
            insj=nsd+j+(nc-1)*nstim
          end if
      do 70 iq=1,ndim
        inxi=i+(iq-1)*nstim
        inxj=j+(iq-1)*nstim
         if(wword) then
           gwq=dble(w(nclass*(iq-1)+nc))
         else
           gwq=1.d0
         end if
       sterm1=gws*gwq*dterm*dble(dif(iq))
       he(inxi,insi)=he(inxi,insi)+sterm1
       he(insi,inxi)=he(insi,inxi)+sterm1
       he(inxj,insi)=he(inxj,insi)-sterm1
       he(insi,inxj)=he(insi,inxj)-sterm1
       he(inxi,insj)=he(inxi,insj)+sterm1
       he(insj,inxi)=he(insj,inxi)+sterm1
       he(inxj,insj)=he(inxj,insj)-sterm1
       he(insj,inxj)=he(insj,inxj)-sterm1
   70 continue
       sterm=gws*gws*dterm/2.d0
       he(insj,insi)=he(insj,insi)+sterm
       he(insi,insj)=he(insi,insj)+sterm
       he(insj,insj)=he(insj,insj)+sterm
       he(insi,insi)=he(insi,insi)+sterm
  100 continue
  300   continue
      go to 98
   99 abort=.true.   
   98 continue
c      print 1717
c
c
c
c$$$ 1718 format(10x,64(f20.10,1h,,2x))
c$$$c
c$$$      do 301 ij=1,nstim
c$$$         do 302 ir1=1,ndim
c$$$            ig1=ij+(ir1-1)*nstim
c$$$            do 303 ik=1,nstim
c$$$               do 304 ir2=1,ndim
c$$$                  ig2=ik+(ir2-1)*nstim
c$$$                  print 1718,he(ig1,ig2)
c$$$ 304           continue
c$$$ 303        continue
c$$$            do 305 is2=1,nstim
c$$$               ig2=is2+nsd
c$$$               print 1718,he(ig1,ig2)
c$$$ 305        continue
c$$$ 302     continue
c$$$ 301  continue
c$$$      do 306 is1=1,nstim
c$$$         ig1=is1+nsd
c$$$         do 307 ik=1,nstim
c$$$            do 308 ir2=1,ndim
c$$$               ig2=ik+(ir2-1)*nstim
c$$$               print 1718,he(ig1,ig2)
c$$$ 308        continue
c$$$ 307     continue
c$$$         do 309 is2=1,nstim
c$$$            ig2=is2+nsd
c$$$            print 1718,he(ig1,ig2)
c$$$ 309     continue
c$$$ 306  continue
c$$$c
c$$$      print 1718,(g(ig),ig=1,nspace)
c
c
c 
      return
      end
      subroutine dmfss(a,n,eps,irank,trac)
c
c  given a symmetric positive semi-definite matrix, dmfss will
c  (1) determine the rank and linearly independent rows and columns
c  (2) factor a symmetric submatrix of maximal rank
c  (3) express nonbasic rows in terms of basic ones
c      express nonbasic columns in terms of basic ones
c      express basic variables in terms of free ones
c  subroutine dmfss may be used as a preparatory step for the calcalation of
c  the least squares solution of minimal length of a system of linear 
c  equations with symmetric positive semi-definite coefficient matrix
c
c  description of parameters
c  -------------------------
c
c     a         upper triangular part of given symmetric semi-definite
c               matrix stored columnwise in compress form.
c               on return a contains the matrix t and, if irank is
c               less than n, the matrix u and tu
c     n         dimension of given matrix a
c     eps       test value for zero affect by round-off noise
c     irank     resultant variable, containing the rank of given matrix a if a
c               is semi-definite
c               irank =  0 means a has no positive diagonal element and/or
c                          eps is not absolutely less than one
c               irank = -1 means dimension n is not positive
c               irank = -2 means complete failure, possibly due to
c                          inadequate relative tolerance eps
c     trac      vector of dimension n containing the source index of the i-th
c               pivot row and its i-th location.
c               this means that trac contains the product representation
c               of the permutation which is applied to rows and columns of a
c               in terms of transpositions
c
c  eps must be absolutely less than one. a sensible value is somewhere
c  between 10**(-4) and 10**(-6).
c  the absolute value of input parameter eps is used as relative tolerance.
c  in order to preserve symmetry only pivoting along the diagonal is built in.
c  all pivot elements must be greater than the absolute value of eps times
c  the original diagonal element, otherwise they are treated as if they were
c  zero.
c  matrix a remans unchanged if the resultant value irank equals zero.
c
      dimension a(1),trac(1)
      double precision a,trac,piv,hold,sum,eps
      if(n)36,36,1
    1 irank=0
      isub=0
      kpiv=0
      j=0
      piv=0.0
      do 3 k=1,n
      j=j+k
      trac(k)=a(j)
      if(a(j)-piv)3,3,2
    2 piv=a(j)
      ksub=j
      kpiv=k
    3 continue
      do 32 i=1,n
      isub=isub+i
      im1=i-1
    4 kmi=kpiv-i
      if(kmi)35,9,5
    5 ji=ksub-kmi
      idc=ji-isub
      jj=isub-im1
      do 6 k=jj,isub
      kk=k+idc
      hold=a(k)
      a(k)=a(kk)
    6 a(kk)=hold
      kk=ksub
      do 7 k=kpiv,n
      ii=kk-kmi
      hold=a(kk)
      a(kk)=a(ii)
      a(ii)=hold
    7 kk=kk+k
      jj=kpiv-1
      ii=isub
      do 8 k=i,jj
      hold=a(ii)
      a(ii)=a(ji)
      a(ji)=hold
      ii=ii+k
    8 ji=ji+1
    9 if(irank)22,10,10
   10 trac(kpiv)=trac(i)
      trac(i)=kpiv
      kk=im1-irank
      kmi=isub-kk
      piv=0.0
      idc=irank+1
      ji=isub-1
      jk=kmi
      jj=isub-i
      do 19 k=i,n
      sum=0.d0
      if(kk)13,13,11
   11 do 12 j=kmi,ji
      sum=sum-a(j)*a(jk)
   12 jk=jk+1
   13 jj=jj+k
      if(k-i)14,14,16
   14 sum=a(isub)+sum
      if(sum-dabs(a(isub)*eps))20,20,15
   15 a(isub)=dsqrt(sum)
      kpiv=i+1
      goto 19
   16 sum=(a(jk)+sum)/a(isub)
      a(jk)=sum
      if(a(jj))19,19,17
   17 trac(k)=trac(k)-sum*sum
      hold=trac(k)/a(jj)
      if(piv-hold)18,19,19
   18 piv=hold
      kpiv=k
      ksub=jj
   19 jk=jj+idc
      goto 32
   20 if(irank)21,21,37
   21 irank=-1
      goto 4
   22 irank=im1
      ii=isub-irank
      ji=ii
      do 26 k=1,irank
      ji=ji-1
      jk=isub-1
      jj=k-1
      do 26 j=i,n
      idc=irank
      sum=0.d0
      kmi=ji
      kk=jk
      if(jj)25,25,23
   23 do 24 l=1,jj
      idc=idc-1
      sum=sum-a(kmi)*a(kk)
      kmi=kmi-idc
   24 kk=kk-1
   25 a(kk)=(sum+a(kk))/a(kmi)
   26 jk=jk+j
      jj=isub-i
      piv=0.0
      kk=isub-1
      do 31 k=i,n
      jj=jj+k
      idc=0
      do 28 j=k,n
      sum=0.d0
      kmi=jj+idc
      do 27 l=ii,kk
      jk=l+idc
   27 sum=sum+a(l)*a(jk)
      a(kmi)=sum
   28 idc=idc+j
      a(jj)=a(jj)+1.0
      trac(k)=a(jj)
      if(piv-a(jj))29,30,30
   29 kpiv=k
      ksub=jj
      piv=a(jj)
   30 ii=ii+k
      kk=kk+k
   31 continue
      goto 4
   32 continue
   33 if(irank)35,34,35
   34 irank=n
   35 return
   36 irank=-1
      return
   37 irank=-2
      return
      end
      subroutine dmlss(a,n,irank,trac,inc,rhs,ier)
c
c  subroutine dmlss is the 2nd step in the procedure for calculating the least
c  squares solution of minimal length of a system of simultaneous linear
c  equations with symmetric positive semi-definite coefficient matrix
c
c  description of parameters
c  -------------------------
c
c     a         coefficient matrix in factor form as generated by subroutine
c               dmfss from initially given symmetric coefficient matrix a
c               stored in n*(n+1)/2 locations.
c               a remains unchanged.
c     n         dimension of coefficient matrix
c     irank     rank of coefficient matrix, calculated by means of subroutine
c               dmfss
c     trac      vector of dimension n containing the subscripts of pivot
c               rows and columns. i.e. the product representation in
c               transpositions of the permutation which as applied to rows
c               and columns of a in the factorization process.
c               trac is a resultant array of subroutine dmfss.
c     inc       input variable which should contain the value zero if the
c               system of linear equations is known to be compatible and a
c               a nonzero value otherwise
c     rhs       vector of dimension n containing the right hand side.
c               on return rhs contains the minimal length solution
c     ier       resultant error parameter
c               ier =  0 means no errors
c               ier = -1 means n and/or irank is not positive and/or irank
c                        is greater than n
c               ier = +1 means the factorization contained in a has zero
c                        divisors and/or trac contains  values outside
c                        the feasible range 1 up to n
c
c  the minimal length solution is produced in the storage locations occupied
c  by the right hand side.
c  subroutine dmlss does take care of the permutation which was applied to rows
c  and columns of a. 
c
      dimension a(1),trac(1),rhs(1)
      double precision a,trac,rhs,sum,hold
      idef=n-irank
      if(n)33,33,1
    1 if(irank)33,33,2
    2 if(idef)33,3,3
    3 ite=irank*(irank+1)/2
      ix2=irank+1
      np1=n+1
      ier=0
      jj=1
      ii=1
    4 do 6 i=1,n
      j=trac(ii)
      if(j)31,31,5
    5 hold=rhs(ii)
      rhs(ii)=rhs(j)
      rhs(j)=hold
    6 ii=ii+jj
      if(jj)32,7,7
    7 isw=1
      if(inc*idef)8,28,8
    8 ista=ite
      do 10 i=1,irank
      ista=ista+1
      jj=ista
      sum=0.d0
      do 9 j=ix2,n
      sum=sum+a(jj)*rhs(j)
    9 jj=jj+j
   10 rhs(i)=rhs(i)+sum
      goto(11,28,11),isw
   11 ista=ite
      do 15 i=ix2,n
      jj=ista
      sum=0.d0
      do 12 j=1,irank
      jj=jj+1
   12 sum=sum+a(jj)*rhs(j)
      goto(13,13,14),isw
   13 sum=-sum
   14 rhs(i)=sum
   15 ista=ista+i
      goto(16,29,30),isw
   16 ista=ix2
      iend=n
      jj=ite+ista
   17 sum=0.d0
      do 20 i=ista,iend
      if(a(jj))18,31,18
   18 rhs(i)=(rhs(i)-sum)/a(jj)
      if(i-iend)19,21,21
   19 jj=jj+ista
      sum=0.d0
      do 20 j=ista,i
      sum=sum+a(jj)*rhs(j)
   20 jj=jj+1
   21 sum=0.d0
      ii=iend
      do 24 i=ista,iend
      rhs(ii)=(rhs(ii)-sum)/a(jj)
      if(ii-ista)25,25,22
   22 kk=jj-1
      sum=0.d0
      do 23 j=ii,iend
      sum=sum+a(kk)*rhs(j)
   23 kk=kk+j
      jj=jj-ii
   24 ii=ii-1
   25 if(idef)26,30,26
   26 goto(27,11,8),isw
   27 isw=2
      goto 8
   28 ista=1
      iend=irank
      jj=1
      isw=2
      goto 17
   29 isw=3
      goto 16
   30 ii=n
      jj=-1
      goto 4
   31 ier=1
   32 return
   33 ier=-1
      return
      end
c ------------------
      subroutine reduc2(g,grhat,he,infree,nfrees,nce,
     1span,hehld,nmax)
c  this subroutine computes reduced hess and grad for model=2
      double precision g(1),grhat(1),hehld(nmax,nmax),
     1he(nmax,nmax),span(nmax,nmax)
      integer infree(1)
      call vecclr(nfrees,grhat)
      do 3 in=1,nce
      do 3 il=1,nce
    3 hehld(in,il)=0.d0
    5 m1=0
      do 20 in=1,nce
       if(infree(in).eq.1) then
         m1=m1+1
         do 21 ik=1,nfrees
            if(ik.eq.m1) then
                 span(ik,in)=1.d0
            else
                span(ik,in)=0.d0
            end if
   21    continue
       else
         do 22 ik=1,nfrees
   22    span(ik,in)=0.d0
       end if
   20 continue
      do 30 il=1,nfrees
        do 31 ik=1,nce
   31   grhat(il)=grhat(il)+span(il,ik)*g(ik)
   30 continue
      do 50 l1=1,nfrees
        do 52 l2=1,nce
          do 51 ik=1,nce
   51     hehld(l1,l2)=hehld(l1,l2)+he(ik,l2)*span(l1,ik)
   52   continue
   50 continue
      do 53 l1=1,nce
      do 53 l2=1,nce
   53 he(l1,l2)=0.0
      do 55 l1=1,nfrees
      do 55 l2=1,nfrees
        do 56 ik=1,nce
   56   he(l1,l2)=he(l1,l2)+span(l2,ik)*hehld(l1,ik)
   55 continue
   99 return
      end
c   -----------------------------------------------------------
      subroutine wfgh(x,w,s,ft,mean1,ss,g,abort,wword,pgword,
     1he,nmax,dif,nclass)    
      common/parblk/nstim,ndim,nsub,nwgt,nmodel,ndat,
     1nsd,nspr,ns,nspace
       real ft(nclass),mean1(nclass,nspr)
      double precision dterm,hterm,sterm,sterm1,gim,gim1
      double precision he(nmax,nmax),g(1),ss,res,ssum
      logical abort,wword,pgword 
      real x(1),w(1),s(1),dif(1)
      nw=nclass*ndim
      if(nmodel.eq.3)nw=nw+nclass
      call vecclr(nw,g)
      ss=0.d0
      do 2 il=1,nw
      do 2 ij=1,nw
    2 he(il,ij)=0.d0
      do 300 nc=1,nclass
      do 100 m=1,nspr
c   go thru data getting dstar and dif(xi-xj) for each dimension
      call dcomp(wword,i,j,nc,m,x,w,s,dstar,dif,nclass)
      dhat=dstar
      if (dstar.eq.0.0)go to 100
      res=dble(mean1(nc,m)-dhat)
      ss=ss+res*res*dble(ft(nc))
      do 15 im=1,ndim
      iml=nc+nclass*(im-1)
              gim1=dble((1.0/dstar)*ft(nc))
            gim=-res*(dif(im)**2)*gim1
   15 g(iml)=g(iml)+gim
      if(nmodel.ne.3)go to 116
      ssum=dble(s(i))+dble(s(j))
      ilw=nclass*ndim+nc
      g(ilw)=g(ilw)-res*ssum*gim1
  116         dterm=dble(((1.0/dstar)**2)*ft(nc))
      do 50 iq=1,ndim
          do 40 ij=1,ndim
            inw1=nc+(iq-1)*nclass
            inw2=nc+(ij-1)*nclass
            hterm=dterm*(dble(dif(ij)*dif(iq))**2)/2.d0
            he(inw1,inw2)=he(inw1,inw2)+hterm
            he(inw2,inw1)=he(inw2,inw1)+hterm
   40     continue
   50 continue
          if(nmodel.ne.3) go to 100
       inws=nc+nclass*ndim
      do 70 iq=1,ndim
        inw1=nc+(iq-1)*nclass
       sterm1=dterm*dble((dif(iq)**2)*ssum/2.d0)
       he(inw1,inws)=he(inw1,inws)+sterm1
       he(inws,inw1)=he(inws,inw1)+sterm1
   70 continue
       sterm=(dterm/2.d0)*ssum*ssum
       he(inws,inws)=he(inws,inws)+sterm
  100 continue
  200  continue
  300   continue
c
c
c
c$$$ 1718   format(10x,8(f20.10,1h,,2x))
c$$$        do 301 ir1=1,ndim
c$$$           do 302 it1=1,nclass
c$$$              ig1=(3-it1)+(ir1-1)*nclass
c$$$              do 303 ir2=1,ndim
c$$$                 do 304 it2=1,nclass
c$$$                    ig2=(3-it2)+(ir2-1)*nclass
c$$$                    print 1718,he(ig1,ig2)
c$$$ 304             continue
c$$$ 303          continue
c$$$              do 305 is2=1,nclass
c$$$                 ig2=(3-is2)+nclass*ndim
c$$$                 print 1718,he(ig1,ig2)
c$$$ 305          continue
c$$$ 302       continue
c$$$ 301    continue
c$$$        do 306 is1=1,nclass
c$$$           ig1=(3-is1)+nclass*ndim
c$$$           do 307 ir2=1,ndim
c$$$              do 308 it2=1,nclass
c$$$                 ig2=(3-it2)+(ir2-1)*nclass
c$$$                 print 1718,he(ig1,ig2)
c$$$ 308          continue
c$$$ 307       continue
c$$$           do 309 is2=1,nclass
c$$$              ig2=(3-is2)+nclass*ndim
c$$$              print 1718,he(ig1,ig2)
c$$$ 309       continue
c$$$ 306    continue
c
c        print 1718,(g(ig),ig=1,(nclass+nclass*ndim)) 
c
c
c
   98 return
      end
c   ------------------------------------------
      subroutine wnorm(x,w,s,nstim,ndim,nsub,nmodel)
      real x(1),w(1),s(1)
      double precision sum
      nend=ndim
      if(nmodel.eq.3)nend=nend+1
      do 10 i=1,nend
         ist=(i-1)*nsub+1
         iend=i*nsub
         sum=0.d0
         do 20 j=ist,iend
   20    sum=sum+dble(w(j))
         sum=sum/dfloat(nsub)
         if(sum.eq.0.d0)go to 10
         do 25 j=ist,iend
   25    w(j)=w(j)/sngl(sum)
         ijst=(i-1)*nstim+1
         ijend=i*nstim
         if(i.le.ndim) then
           do 30 ij=ijst,ijend
   30      x(ij)=x(ij)*sngl(dsqrt(sum))
         else
           if(nmodel.le.3)nss=nstim
           if(nmodel.eq.4)nss=nstim*nsub
           do 35 ij=1,nss
   35      s(ij)=s(ij)*sngl(sum)
         end if
   10 continue
      return
      end
      subroutine ranin(lambda,mean,sigma2,y,nsub,nstim,nclass,iseed)
c
c  generate random initial parameter estimates
c
      real lambda(nclass),mean(nclass,nstim),y(nsub,nstim)
c
c  generate uniform distributed random number for mixing parameters
c  make sure they add up to 1
      sum=0.0
      do 10 i=1,nclass
        lambda(i)=runif(iseed)
 10      sum=sum+lambda(i)
      do 20 i=1,nclass
 20      lambda(i)=lambda(i)/sum
c
c  sample for each stimulus j the class means from the uniform distribution
c  on [yj_min,yj_max], where yj_min (yj_max) is the min. (max.) value
c  over the subjects for stimulus j.
      do 50 j=1,nstim
        ymax=y(1,j)
        ymin=ymax
        do 30 i=1,nsub
           ymax=amax1(y(i,j),ymax)
           ymin=amin1(y(i,j),ymin)
 30         continue
        do 40 i=1,nclass
 40         mean(i,j)=ymin+(ymax-ymin)*runif(iseed)
 50      continue
c
c  use pooled variance with repect to the average class mean per
c  stimulus as initial estimate of sigma2
      sigma2=0.0
      do 80 j=1,nstim
        av=0.0
        do 60 i=1,nclass
 60        av=av+mean(i,j)
        av=av/float(nclass)
        do 70 i=1,nsub
 70         sigma2=sigma2+(y(i,j)-av)**2
 80      continue
      sigma2=sigma2/float(nsub*nstim)
c
      return
      end
      subroutine ratin(lambda,mean,sigma2,y,nsub,nstim,nclass)
c
c  compute rational initial estimates
c
      common /cstak/ dstak
c
      real lambda(nclass),mean(nclass,nstim),y(nsub,nstim)
      double precision sum
c
      double precision dstak(300000000)
      real rstak(600000000)
      integer istak(600000000)
c
      equivalence(dstak(1),rstak(1))
      equivalence(dstak(1),istak(1))
c
      lp=i1mach(2)
c
c  perform a k-means clustering on Y
      npsub=nsub*(nsub-1)/2
      jmem=istkgt(nsub,2)
      jsize=istkgt(nclass,2)
      jtmp1=istkgt(npsub,3)
      jtmp2=istkgt(nclass*nstim,3)
      call rkmean(y,istak(jmem),istak(jsize),rstak(jtmp1),
     1   rstak(jtmp2),mean,nsub,nstim,npsub,nclass)
      call istkrl(2)
c
c      write(lp,10) (istak(jmem+i-1),i=1,nsub)
 10   format(/' Results of k-kmeans clustering on data ',
     1   ' (cluster membership per subject):'//(1x,40i3))
c      write(lp,20) (istak(jsize+i-1),i=1,nclass)
 20   format(/' Number of subjects per cluster: ',20i5/
     1   (1x,20i5))
c
c  compute an initial estimate of lambda by taking the relative cluster sizes
      do 30 i=1,nclass
 30      lambda(i)=float(istak(jsize+i-1))/float(nsub)
c
c  compute initial estimate of sigma2
      sum=0.d0
      do 80 i=1,nsub
        k=istak(jmem+i-1)
        do 70 j=1,nstim
 70         sum=sum+dble(y(i,j)-mean(k,j))**2
 80      continue
      sigma2=sngl(sum/dfloat(nstim*nsub))
c
c  release memory for mem and size
      call istkrl(2)
      return
      end
      subroutine ratini(lambda,mean,y,nsub,nstim,nclass,mem)
c
c  compute rational initial estimates
c
      common /cstak/ dstak
c
      real lambda(nclass),mean(nclass,nstim),y(nsub,nstim)
c
      double precision dstak(300000000)
      real rstak(600000000)
      integer istak(600000000),mem(nsub)
c
      equivalence(dstak(1),rstak(1))
      equivalence(dstak(1),istak(1))
c
      lp=i1mach(2)
c
c  perform a k-means clustering on Y
      npsub=nsub*(nsub-1)/2
      jsize=istkgt(nclass,2)
      jtmp1=istkgt(npsub,3)
      jtmp2=istkgt(nclass*nstim,3)
      call rkmean(y,mem,istak(jsize),rstak(jtmp1),
     1   rstak(jtmp2),mean,nsub,nstim,npsub,nclass)
      call istkrl(2)
c
      write(lp,10) (mem(i),i=1,nsub)
 10   format(/' Results of k-kmeans clustering on data ',
     1   ' (cluster membership per subject):'//(1x,40i3))
      write(lp,20) (istak(jsize+i-1),i=1,nclass)
 20   format(/' Number of subjects per cluster: ',20i5/
     1   (1x,20i5))
c
c  compute an initial estimate of lambda by taking the relative cluster sizes
      do 30 i=1,nclass
 30      lambda(i)=float(istak(jsize+i-1))/float(nsub)
c
c release mem for size
      call istkrl(1)
      return
      end
      subroutine rkmean(x,mem,klus,dsq,ksum,centr,n,nvar,npair,nclus)
c
c  x(n,nvar)   contains the scores of n objects on nvar real variables
c  dsq(npair)  squared Euclidean distances between objects
c              npair = n*(n-1)/2
c  klus(nclus) used for temporary storage (contains upon return the
c              no. of objects per cluster)
c  mem(n)      mem(i) is number of cluster to which object i belongs
c  centr(nclus,nvar)   contains the cluster centroids
c  ksum(nclus,nvar)    used for temporary storage
c
      real x(n,nvar),dsq(npair),maxdsq,mindsq
      integer klus(nclus),mem(n)
      real ksum(nclus,nvar),centr(nclus,nvar)
      logical chngd
c
      if (nclus .lt. 2)
     1   call seterr('rkmean - nclus .lt. 2',21,1,2)
c
c  compute squared Euclidean distances between the rows of x
c  and find the most distant pair of objects
      ij=0
      maxdsq=0.0
      do 30 i=2,n
        im1=i-1
        do 20 j=1,im1
           ij=ij+1
           dsq(ij)=0.0
           do 10 k=1,nvar
 10            dsq(ij)=dsq(ij)+(x(i,k)-x(j,k))**2
c
c
c
c$$$ 1979          format(10x,i4,i4,i4,f20.0)
c$$$               print 1979,ij-1,j-1,i-1,dsq(ij)
c
c
c
           if (dsq(ij) .le. maxdsq) goto 20
              maxdsq=dsq(ij)
              maxi=i
              maxj=j
 20         continue
 30      continue
      if (maxdsq .eq. 0.0)
     1  call seterr('rkmean - no variance in data',28,2,2)
c
c  find representative objects for the other clusters.
c  if s and t denote the objects representing the first 2 clusters,
c  find object i such that d(is)**2 + d(it)**2 is maximal, next
c  find object j such that d(js)**2 + d(jt)**2 + d(ij)**2 is max., etc.
      klus(1)=maxi
      klus(2)=maxj
      m=2
      if (nclus .eq. 2) goto 90
      do 80 k=3,nclus
        maxdsq=0.0
        do 70 i=1,n
           do 50 l=1,m
              if (i .eq. klus(l)) goto 70
 50            continue
           sum=0.0
           do 60 j=1,m
              l=klus(j)
              il=igetij(i,l)
 60            sum=sum+dsq(il)
           if (sum .le. maxdsq) goto 70
              maxdsq=sum
              maxi=i
 70         continue
        if (maxdsq .eq. 0.0)
     1      call seterr('rkmean - cannot find seeds',26,3,2)
        klus(k)=maxi
        m=m+1
 80      continue
 90   continue
c
c
c
c$$$ 1979 format(10(i4))
c$$$      print 1979,(klus(k)-1,k=1,nclus)
c
c
c
c
c  assign all objects to one of the initial clusters
      do 110 i=1,n
        mindsq=1.e+35
        do 100 k=1,nclus
           l=klus(k)
           if (i .ne. l) goto 105
              mem(i)=k
              goto 110
 105       il=igetij(i,l)
           if (dsq(il) .ge. mindsq) goto 100
              mindsq=dsq(il)
              kmin=k
 100        continue
        mem(i)=kmin
 110     continue
c
c  compute initial cluster centroids in centr
      do 130 i=1,nclus
        klus(i)=0
        do 120 j=1,nvar
 120        ksum(i,j)=0.0
 130     continue
      do 150 i=1,n
        k=mem(i)
        do 140 j=1,nvar
 140        ksum(k,j)=ksum(k,j)+x(i,j)
        klus(k)=klus(k)+1
 150     continue
      do 170 i=1,nclus
        f=float(klus(i))
        do 160 j=1,nvar
 160        centr(i,j)=ksum(i,j)/f
 170     continue
c
c  main iteration loop starts here
c
 200  continue
c
      chngd=.FALSE.
c
c  check for each object to which centroid it is closest
      do 240 i=1,n
        dsqmin=1.e+35
        do 220 k=1,nclus
           dsqk=0.0
           do 210 j=1,nvar
 210           dsqk=dsqk+(x(i,j)-centr(k,j))**2
           if (dsqk .ge. dsqmin) goto 220
              dsqmin=dsqk
              kmin=k
 220        continue
c
c  if kmin.ne.mem(i), update sum and centroid for cluster mem(i)
c  and cluster kmin
        if (kmin .eq. mem(i)) goto 240
           k=mem(i)
           if (klus(k) .eq. 1) goto 240
              chngd=.TRUE.
              klus(k)=klus(k)-1
              klus(kmin)=klus(kmin)+1
              do 230 j=1,nvar
                 ksum(k,j)=ksum(k,j)-x(i,j)
                 ksum(kmin,j)=ksum(kmin,j)+x(i,j)
                 centr(k,j)=ksum(k,j)/float(klus(k))
                 centr(kmin,j)=(ksum(kmin,j))/float(klus(kmin))
 230              continue
              mem(i)=kmin
 240     continue
      if (chngd) goto 200
c
      return
      end
      subroutine rundvr(lambda,mean,sigma2,y,z,
     1      nclass,eps,logli,temp,npar,aic,bic,
c   stuff for mstep job
c   logical arguments
     1test,wword,abort,pgword,indscwd,
c   real arguments
     3x,w,s,
     4nmax,nmaxp,nmaxs,
     2nmaxw,nmaxpw,ndimmax,nrep,icarlo,iseed,win,iconsc)
      common /cstak/dstak
      double precision dstak(300000000)
      real rstak(600000000)
      integer istak(600000000) 
      common/parblk/nstim,ndim,nsub,nwgt,nmodel,ndat,
     1nsd,nspr,ns,nspace
      common/cnvblk/conv,xsconv,wconv,
     1itmax,itxmax,itwmax
      real x(nmax),w(nmaxw),s(nmaxs),lrmc(1000),lrobs
      logical test,wword,abort,pgword,indscwd,wword1
      real lambda(nclass),mean(nclass,nspr),y(nsub,nspr)
      real z(nsub,nclass),win(nsub,ndim)
      double precision temp(nclass),fit1,fit2,logli
      equivalence(dstak(1),rstak(1))
      equivalence(dstak(1),istak(1))
      save maxit
c
      in=i1mach(1)
      lp=i1mach(2)
      lunout=12
      if(nsub.eq.nclass)go to 5
      if((nrep.gt.19).and.(icarlo.eq.1)) go to 20
c
c   n o   m o n t e   c a r l o   s i g n i f i c a n c e   t e s t i n g
c
 5    continue
      if(iconsc.eq.1)read(in,1008)fit0
 1008 format(d20.7)
      if(iconsc.eq.1)open(lunout,file='lratio',status='unknown')
       write(lp,10)nclass
   10    format(///'Analysis with',i4,' latent classes'/
     1       '---------------------------------'/)
       call emalg(lambda,mean,sigma2,y,z,
     1      nclass,eps,logli,temp,npar,aic,bic,
c   stuff for mstep job
c   logical arguments
     1test,wword,abort,pgword,
c   real arguments
     3x,w,s,
     4nmax,nmaxp,nmaxs,
     2nmaxw,nmaxpw,ndimmax,nrep)
      if(iconsc.eq.1) then
        lrobs=sngl(2.d0*(logli-fit0))
        write(lunout,1009)lrobs
 1009   format(2x,f12.3)
        close(lunout)
      end if
      call wnorm(x,w,s,nstim,ndim,nclass,nmodel)
      return
c
c
c      m o n t e   c a r l o   s i g n i f i c a n c e   t e s t i n g
c
   20   continue
       read(in,1000)nclass1,ndim1,ns1
 1000 format(3i5)
       write(lp,1001)
 1001  format(/" Monte Carlo Hope Test for test of Model 1 vs Model 2")
       write(lp,1002)nclass,ndim
 1002  format(/' For Model 1 nclass=',i5,' ndim=',i5)
       if(ns.eq.0)write(lp,1003)
       if(ns.eq.1)write(lp,1004)
       if(ns.eq.2)write(lp,1005)
       write(lp,1006)nclass1,ndim1
 1006  format(/' For Model 2 nclass=',i5,' ndim=',i5)
       if(ns1.eq.0)write(lp,1003)
       if(ns1.eq.1)write(lp,1004)
       if(ns1.eq.2)write(lp,1005)
 1003  format('   No specificities: Euclidean Model')
 1004  format('   Specificites: Restricted Extended Euclidean Model')
 1005  format('   Different Pattern of Specificites:',
     1 ' Restricted Extended Euclidean Model')
       nclp1=nclass1
       nsd1=nstim*ndim1
       nspace1=nsd1
       nmax1=nspace1+nstim
       if(ns1.eq.1)nspace1=nspace1+nstim
       if(ns1.eq.2)nspace1=nspace1+nstim*nclass1
          if((nclass1.eq.1).and.(ns1.eq.0))nmodel1=0
          if((nclass1.eq.1).and.(ns1.ge.1))nmodel1=2
          if((nclass1.gt.1).and.(ns1.eq.0))nmodel1=1
          if((nclass1.gt.1).and.(ns1.eq.1))nmodel1=3
          if((nclass1.gt.1).and.(ns1.eq.2))nmodel1=4
           wword1=.true.
           if(nclass1.eq.1)wword1=.false.
	nw1=nclass1*ndim1
	if(nmodel1.eq.3)nw1=nw1+nclass1
        nwgt1=nw1
        nmaxw1=nw1 + ndim1 + 1
c
       if(nmodel1.ge.2)nmaxw1=nmaxw1+nclp1
       nmaxpw1=(nmaxw1*(nmaxw1+1))/2
c  compute no of active parameters for model2
c  adjust for model
       go to(11,21,31,41,51)nmodel1+1
   11    npar1=nsd1-ndim1*(ndim1+1)/2
         go to 33
   21    npar1=nsd1+ndim1*nclass1-2*ndim1
         go to 33
   31    npar1=nsd1-ndim1*(ndim1+1)/2+nstim
         go to 33
   41    npar1=nsd1+nclass1*(ndim1+1)-2*ndim1-1+nstim
         go to 33
   51      npar1=nsd1+nclass1*ndim1-2*ndim1+nstim*nclass1
   33  continue 
         npar1= npar1+nclass1
       jlam0=istkgt(nclass,3)
       jlam1=istkgt(nclp1,3)
       jmean0=istkgt(nspr*nclass,3)
       jmean1=istkgt(nspr*nclp1,3)
       jz0=istkgt(nsub*nclass,3)
       jz1=istkgt(nsub*nclp1,3)
       jtemp0=istkgt(nclass,4)
       jtemp1=istkgt(nclp1,4)
       jw0=istkgt(nmaxw,3)
       jwmean0=istkgt(ndim*nclass,3)
       jw1=istkgt(nmaxw1,3)
       jwmean1=istkgt(ndim1*nclp1,3)
       jxnew=istkgt(nmax1,3)
       jynew=istkgt(nspr*nsub,3)
       jcum1=istkgt(nclp1,3)
        if(ns1.lt.2) then
         nmaxs1=nstim
        else
         nmaxs1=nstim*nclp1
        end if 
       jsnew=istkgt(nmaxs1,3)
c  initialize xnew w0 w1
        do 16 ii=1,nmaxw
   16        rstak(jw0+ii-1)=1.0
        do 12 ii=1,nmaxw1
   12        rstak(jw1+ii-1)=1.0
        do 13 ii=1,nmax1
   13        rstak(jxnew+ii-1)=x(ii)
         if(ns1.lt.2)nup=nstim
         if (ns1.eq.2)nup=nstim*nclp1
        do 14 ii=1,nup
   14        rstak(jsnew+ii-1)=1.0
        write(lp,10)nclass
       call emalg(lambda,mean,sigma2,y,z,
     1      nclass,eps,fit1,temp,npar,aic,bic,
c   stuff for mstep job
c   logical arguments
     1test,wword,abort,pgword,
c   real arguments
     3x,w,s,
     4nmax,nmaxp,nmaxs,
     2nmaxw,nmaxpw,ndimmax,nrep)
      call wnorm(x,w,s,nstim,ndim,nclass,nmodel)
 1616   format(//'sig=',2x,e15.6)
        if(nclass1.eq.1) then
          call onecl(rstak(jlam1),rstak(jmean1),si21,y,nsub,nstim)
        else
         if (indscwd) then
          call wstart(nmodel1,rstak(jlam1),win,rstak(jw1),nsub,ndim1,
     1     nclp1,rstak(jwmean1),rstak(jmean1),nspr,nmaxw1,y,si21)
         else
          call ratin(rstak(jlam1),rstak(jmean1),si21,y,nsub,
     1      nspr,nclp1)
         end if
        end if
        call logl(rstak(jlam1),rstak(jmean1),si21,y,rstak(jz1),
     1      nclp1,nspr,nsub,dstak(jtemp1),fit2,logli,npar1,ac,bc)
       ndim0=ndim
       ndim=ndim1
       nwgt0=nwgt
       nwgt=nwgt1
       nmodel0=nmodel
       nmodel=nmodel1
       nsd0=nsd
       nsd=nsd1
       ns0=ns
       ns=ns1
       nspace0=nspace
       nspace=nspace1
c
      call emalg(rstak(jlam1),rstak(jmean1),si21,y,rstak(jz1),
     1      nclp1,eps,fit2,dstak(jtemp1),npar1,ac,bc,
c   stuff for mstep job
c   logical arguments
     1test,wword1,abort,pgword,
c   real arguments
     3rstak(jxnew),rstak(jw1),rstak(jsnew),
     4nmax1,nmaxp1,nmaxs1,
     2nmaxw1,nmaxpw1,ndimmax,nrep)
      ndim=ndim0
      nwgt=nwgt0
      nmodel=nmodel0
      nsd=nsd0
      ns=ns0
      nspace=nspace0
      lrobs=sngl(2.d0*(fit2-fit1))
      logli=fit1
      write(lp,30)lrobs
   30   format(/'Likelihood ratio statistic for testing',
     1      ' model 1 vs. model 2 :',f20.5)
c  using model parameters generate model distances for each class
       jdtmp=istkgt(nclass*nspr,3)
       call gtpmod(ns,ndim,nstim,nclass,nspr,nmax,nmaxw,
     1     x,w,s,rstak(jdtmp))
       do 40 irep=1,nrep
c       print 9898,irep
 9898  format(/1x,5x,"replication no : ",i3)
c       call gtpmodr(ns,ndim,nstim,nclass,nspr,nmax,nmaxw,
c     1 rstak(jxnew),rstak(jw0),rstak(jsnew),rstak(jdtmp),rstak(jlam0))
c       ncm1=nclass-1
c       cuml=0.0
c       do 41 i=1,ncm1
c           rstak(jlam0+i-1)=runif(iseed)
c 41        cuml=cuml+rstak(jlam0+i-1)
c       rstak(jlam0+ncm1)=1.0-cuml
c      call genery(iseed,nsub,nspr,nclass,nclp1,rstak(jlam0),rstak(jcum1)
c     1    ,rstak(jdtmp),rstak(jynew),sigma2)
      call genery(iseed,nsub,nspr,nclass,nclp1,lambda,rstak(jcum1),
     1    rstak(jdtmp),rstak(jynew),sigma2)
        if(nclass.eq.1) then
          call onecl(rstak(jlam0),rstak(jmean0),si20,rstak(jynew),
     1        nsub,nspr)
             if(ns.eq.0) then
             nmodel=0
           else
             nmodel=2
           end if
           wword=.false.
        else
          if(indscwd) then
            call wstart(nmodel,rstak(jlam0),win,rstak(jw0),nsub,ndim,
     1nclass,rstak(jwmean0),rstak(jmean0),nspr,nmaxw,rstak(jynew),si20)
          else
           call ratin(rstak(jlam0),rstak(jmean0),si20,rstak(jynew),
     1       nsub,nspr,nclass)
          end if
        end if
        call logl(rstak(jlam0),rstak(jmean0),si20,rstak(jynew),
     1      rstak(jz0),nclass,nspr,nsub,dstak(jtemp0),fit1,logli,
     2      npar,ac,bc)
      call emalg(rstak(jlam0),rstak(jmean0),si20,rstak(jynew),
     1      rstak(jz0),nclass,eps,fit1,dstak(jtemp0),npar,ac,bc,
c   stuff for mstep job
c   logical arguments
     1test,wword,abort,pgword,
c   real arguments
     3rstak(jxnew),rstak(jw0),rstak(jsnew),
     4nmax,nmaxp,nmaxs,
     2nmaxw,nmaxpw,ndimmax,nrep)
        if(nclass1.eq.1) then
          call onecl(rstak(jlam1),rstak(jmean1),si21,
     1      rstak(jynew),nsub,nstim)
        else
         if(indscwd) then
          call wstart(nmodel1,rstak(jlam1),win,rstak(jw1),nsub,ndim1,
     1nclp1,rstak(jwmean1),rstak(jmean1),nspr,nmaxw1,rstak(jynew),si21)
         else 
          call ratin(rstak(jlam1),rstak(jmean1),si21,rstak(jynew),
     1      nsub,nspr,nclp1)
         end if
        end if
       call logl(rstak(jlam1),rstak(jmean1),si21,rstak(jynew),
     1      rstak(jz1),nclp1,nspr,nsub,dstak(jtemp1),fit2,logli,
     2      npar1,ac,bc)
       ndim=ndim1
       nwgt=nwgt1
       nmodel=nmodel1
       nsd=nsd1
       ns=ns1
       nspace=nspace1
c
      call emalg(rstak(jlam1),rstak(jmean1),si21,rstak(jynew),
     1      rstak(jz1),nclp1,eps,fit2,dstak(jtemp1),npar1,ac,bc,
c   stuff for mstep job
c   logical arguments
     1test,wword1,abort,pgword,
c   real arguments
     3rstak(jxnew),rstak(jw1),rstak(jsnew),
     4nmax1,nmaxp1,nmaxs1,
     2nmaxw1,nmaxpw1,ndimmax,nrep)
c
      ndim=ndim0
      nwgt=nwgt0
      nmodel=nmodel0
      nsd=nsd0
      ns=ns0
      nspace=nspace0
      lrmc(irep)=sngl(2.d0*(fit2-fit1))    
 40      continue
c
      call rsort(lrmc,nrep)
      write(lp,45)(lrmc(i),i=1,nrep)
 45   format(/' Likelihood ratio statistic for Monte Carlo runs:'/
     1   (1x,5f20.5))
c
c  compute p-value
      nsmall=0
      nlarge=0
      do 50 irep=1,nrep
         if (lrmc(irep) .lt. lrobs) nsmall=nsmall+1
         if (lrmc(irep) .gt. lrobs) nlarge=nlarge+1
 50      continue
      p1=float(nlarge)/float(nrep+1)
      p2=1.0-float(nsmall)/float(nrep+1)
      write(lp,60) p1,p2
 60   format(/' Estimated probability of the observed value of',
     1   ' the test statistic: ',f5.2,' < p <',f5.2)
      i=ifix(0.95*float(nrep+1))
      if (nsmall .ge. i) then
             write(lp,70)
 70   format(' Null hypothesis model 1 is true',
     1   ' rejected at 0.05 significance level')
      else
             test=.true. 
             write(lp,80)
 80   format(' Null hypothesis model 1 is true',
     1   ' not rejected at 0.05 significance level')
      end if
c     
c  release memory 
      call istkrl(17)
      return
      end
      subroutine rundvrp(lambda,mean,sigma2,y,z,
     1      nclass,nclp1,logli,temp,npar,aic,bic,
     2nrep,iseed,test)
      common /cstak/dstak
      double precision dstak(300000000)
      real rstak(600000000)
      integer istak(600000000) 
      common/parblk/nstim,ndim,nsub,nwgt,nmodel,ndat,
     1nsd,nspr,ns,nspace
      common/cnvblk/conv,xsconv,wconv,
     1itmax,itxmax,itwmax
      real lrmc(1000),lrobs
      real lambda(nclass),mean(nclass,nspr),y(nsub,nspr)
      real z(nsub,nclass)
      logical test
      double precision temp(nclass),fit1,fit2,logli
      equivalence(dstak(1),rstak(1))
      equivalence(dstak(1),istak(1))
      save maxit
c
      lp=i1mach(2)
c
c
c      m o n t e   c a r l o   s i g n i f i c a n c e   t e s t i n g
c
c
       npar=nspr*nclass +nclass
       npar1=nspr*nclp1+nclp1
       jlam0=istkgt(nclass,3)
       jlam1=istkgt(nclp1,3)
       jmean0=istkgt(nspr*nclass,3)
       jmean1=istkgt(nspr*nclp1,3)
       jz0=istkgt(nsub*nclass,3)
       jz1=istkgt(nsub*nclp1,3)
       jtemp0=istkgt(nclass,4)
       jtemp1=istkgt(nclp1,4)
       jynew=istkgt(nspr*nsub,3)
       jcum1=istkgt(nclp1,3)
       if(nclass.eq.1)go to 23
       call ratin(lambda,mean,sigma2,y,nsub,nspr,nclass)
       go to 43
 23    call onecl(lambda,mean,sigma2,y,nsub,nspr)
 43    continue
       call logl(lambda,mean,sigma2,y,z,nclass,nspr,nsub,temp,
     1     fit1,logli,npar,aic,bic)
         call ratin(rstak(jlam1),rstak(jmean1),si21,y,nsub,
     1      nspr,nclp1)
        call logl(rstak(jlam1),rstak(jmean1),si21,y,rstak(jz1),
     1      nclp1,nspr,nsub,dstak(jtemp1),fit2,logli,npar1,ac,bc)
      lrobs=sngl(2.d0*(fit2-fit1))
      logli=fit1
      write(lp,30)nclass,nclp1,lrobs
   30   format(/'Likelihood ratio statistic for testing',i3,
     1      ' classes vs.',i3,' classes:',f12.3)
c  using model parameters generate model distances for each class
       jdtmp=istkgt(nclass*nspr,3)
       nprc=nclass*nspr
       do 40 irep=1,nrep
         icc=0
       do 78 ii=1,nclass
        do 76 jj=1,nspr
         amean=mean(ii,jj)
         rstak(jdtmp+icc)=amean
 76      icc=icc+1
 78      continue
c       print 9898,irep
 9898  format(/1x,5x,"replication no : ",i3)
        if(nclass.eq.1)go to 90
        nm1=nclass-1
        cuml=0.0
        do 91 ii=1,nm1
           ttlam=lambda(ii)
           rstak(jlam0+ii-1)=ttlam
           cuml=cuml+rstak(jlam0+ii-1)
 91         continue
        rstak(jlam0+nm1)=1.0-cuml
 90     continue
      var=sigma2
      if(var.le.0.0)var=sigma2
      call genery(iseed,nsub,nspr,nclass,nclp1,rstak(jlam0),rstak(jcum1)
     1   ,rstak(jdtmp),rstak(jynew),var)
        if(nclass.eq.1) then
          call onecl(rstak(jlam0),rstak(jmean0),si20,rstak(jynew),
     1        nsub,nspr)
         else
           call ratin(rstak(jlam0),rstak(jmean0),si20,rstak(jynew),
     1       nsub,nspr,nclass)
          end if
        call logl(rstak(jlam0),rstak(jmean0),si20,rstak(jynew),
     1      rstak(jz0),nclass,nspr,nsub,dstak(jtemp0),fit1,logli,
     2      npar,ac,bc)
        call ratin(rstak(jlam1),rstak(jmean1),si21,rstak(jynew),
     1      nsub,nspr,nclp1)
       call logl(rstak(jlam1),rstak(jmean1),si21,rstak(jynew),
     1      rstak(jz1),nclp1,nspr,nsub,dstak(jtemp1),fit2,logli,
     2      npar1,ac,bc)
      lrmc(irep)=sngl(2.d0*(fit2-fit1))    
 40      continue
c
      call rsort(lrmc,nrep)
      write(lp,45)(lrmc(i),i=1,nrep)
 45   format(/' Likelihood ratio statistic for Monte Carlo runs:'/
     1   (1x,10f12.3))
c
c  compute p-value
      nsmall=0
      nlarge=0
      do 50 irep=1,nrep
         if (lrmc(irep) .lt. lrobs) nsmall=nsmall+1
         if (lrmc(irep) .gt. lrobs) nlarge=nlarge+1
 50      continue
      p1=float(nlarge)/float(nrep+1)
      p2=1.0-float(nsmall)/float(nrep+1)
      write(lp,60) p1,p2
 60   format(/' Estimated probability of the observed value of',
     1   ' the test statistic: ',f5.2,' < p <',f5.2)
      i=ifix(0.95*float(nrep+1))
      if (nsmall .ge. i) then
             write(lp,70) nclass
 70   format(' Null hypothesis #class =',i3,
     1   ' rejected at 0.05 significance level')
      else
c             test=.true. 
             write(lp,80) nclass
 80   format(' Null hypothesis #class =',i3,
     1   ' not rejected at 0.05 significance level')
      end if
c     
c  release memory 
      call istkrl(11)
      return
      end
        REAL FUNCTION RUNIF(IX)
C
C  PORTABLE RANDOM NUMBER GENERATOR DUE TO L. SCHRAGE.
C
C***** RUNS ONLY ON COMPUTERS WHICH CAN REPRESENT ALL INTEGERS
C***** BETWEEN 0 AND (2**31)-1.
C
        INTEGER A,P,IX,B15,B16,XHI,XALO,LEFTLO,FHI,K
       SAVE A, B15, B16, P
        DATA A/16807/,B15/32768/,B16/65536/,P/2147483647/
C
C  GET 15 HI ORDER BITS OF IX
        XHI=IX/B16
C  GET 16 LO BITS OF IX AND FORM LO PRODUCT
        XALO=(IX-XHI*B16)*A
C  GET 15 HI ORDER BITS OF LO PRODUCT
        LEFTLO=XALO/B16
C  FORM THE 31 HIGHEST BITS OF FULL PRODUCT
        FHI=XHI*A+LEFTLO
C  GET OVERFLO PAST 31ST BIT OF FULL PRODUCT
        K=FHI/B15
C  ASSEMBLE ALL THE PARTS AND PRESUBTRACT P
C  (THE PARENTHESES ARE ESSENTIAL)
        IX=(((XALO-LEFTLO*B16)-P)+(FHI-K*B15)*B16)+K
C  ADD BACK IN IF NECESSARY
        IF(IX.LT.0) IX=IX+P
C  MULTIPLY BY 1/(2**31-1)
        RUNIF=FLOAT(IX)*4.656612875E-10
        RETURN
        END
      subroutine sear(n,x,y,g,h,d,fo,f,ind,slim,itermx,
     1trial,xd,w,s,ft,mean1,abort,wword,pgword,xwrd,dif,
     2mean,nclass,he,nmax,hew,nmaxw)    
c  this subroutine does a linesearch for dimension parameters
      common / cstak /dstak
      common/parblk/nstim,ndim,nsub,nwgt,nmodel,ndat,
     1nsd,nspr,ns,nspace
      double precision dstak(300000000)
      real rstak(600000000)
      integer istak(600000000)
      logical pgword,test1,test2,test3,test4,test5
     1,abort,wword,xwrd
c
      double precision f,fo,f0,f1,f2,f4,a0,a1,a2,b0,b1,b2,ss,
     1a3,a4,b3,b4,f3
      real mean1(nclass,nspr),mean(nclass,nspr),ft(nclass)
      double precision w1,amax,z,slim,atemp
      real x(1),w(1),dif(1),s(1),y(1),xd(1)
      double precision g(1),h(1),d(1),he(nmax,nmax),hew(nmaxw,nmaxw)
      double precision fact,aa,bb
      equivalence(dstak(1),rstak(1))
      equivalence(dstak(1),istak(1))
      data fact,aa,bb/0.1d0,2.d0,0.5d0/
      ind=0
c  ind is error parameter:0 normal return; 2 no conv of iterations;
c  3 initial slope is nonnegative
c  pgword true will print output at each iteration
c  compute slope  ... slopes are stored in b variables
      ss=0.d0
      b0=0.d0 
      do 3 ii=1,n
        ss=ss+d(ii)*d(ii)
        b0=b0+g(ii)*d(ii)
    3      continue
      ss=dsqrt(ss)
      if(b0) 30,20,20
   20 ind=3
      return
   30 if(xwrd)call vecsc(n,d,1.d0/ss)
c  steps are stored in a variables, functions in f variables
      call swap(a0,f0,b0,0.d0,fo,b0/ss)
      call swap(a1,b1,f1,0.d0,b0,fo)
      call swap(a2,b2,f2,0.d0,b0,fo)
      ips=0
c  first step set to trial
      a4=trial
      itr=0
      if(pgword) print 1000,itr,a0,b0,f0
 1000 format(1x,5x,"linesearch",2x,i3,3(1x,f12.5))
c add check on step size if model eq 2 and xwrd or wword
      if ((nmodel.lt.2).and.(xwrd))go to 37
c compute largest admissible step size
      amax=1d+30
      if(slim)37,35,35
   35 if(xwrd) then
       nst=nsd+1
      else
       nst=1
      end if
      do 50 i=nst,n
        if(d(i).ge.0.d0)go to 50
        if(amax.le.((slim-dble(x(i)))/d(i))) then
        amax=amax
         else
         amax=(slim-dble(x(i)))/d(i)
      end if
   50 continue
      if(a4.gt.amax)a4=amax
c  main loop for search
   37 do 40 itr=1,itermx
      if(a4.gt.0.d0) go to 70
      if(pgword)print 1419,a4
 1419 format(1x,5x,"stepsize 0 or neg =",1x,f10.3)
       f=fo
      call rvecex(n,1,1,1,1,x,y,1.0)
      call vecex(n,1,1,1,1,g,h,1.d0)
      ind=0
      return
   70 do 80 i=1,n
   80      y(i)=x(i)+sngl(a4*d(i))
        if(xwrd) then
         jxx=istkgt(nmax,3)
         if(ns.lt.2)nss=nstim
         if(ns.eq.2)nss=nstim*nclass
         jss=istkgt(nss,3)
        do 81 i=1,nsd
   81  rstak(jxx-1+i)=y(i)
      if(nmodel.lt.2)go to 83
       do 82 i=1,nss
       iup=i+nsd
   82  rstak(jss-1+i)=y(iup)
   83  call dimfgh(rstak(jxx),w,rstak(jss),ft,mean1,f,h,
     1abort,wword,pgword,he,nmax,dif,nclass,mean)
        call istkrl (2)
      else
c
c
c
c$$$ 1979    format(10x,10(f12.6,1h,,2x))
c$$$         print 1979,(x(iash),iash=1,(nclass+nclass*ndim))
c$$$         print 1979,(d(iash),iash=1,(nclass+nclass*ndim))
c$$$         print 1979,(y(iash),iash=1,(nclass+nclass*ndim))
c
c
c
       call wfgh(xd,y,s,ft,mean1,f,h,abort,wword,pgword,hew,
     1     nmaxw,dif,nclass)
      end if
   84 f4=f
      b4=0.d0
      do 11 ii=1,n
       b4=b4+h(ii)*d(ii)
   11    continue
      if(pgword)print 1000,itr,a4,b4,f4
      atemp=a4
c  compute new step
      test1=dabs(b4).lt.fact*dabs(b0)
      test2=f4.gt.f0
      test3=b4.gt.0d0
      test4=(test1.or..not.test3).and.test2
      test5=(b4.gt.b2.and.b2.ge.b1)
      if(.not.test4) go to 110
c  function is worse and either slope is satisfactory or negative
      ips=0
c  step is reduced by factor bb
      a4=a4*bb
      if(pgword)print 1417
 1417 format(1x,5x,"fcn is worse ..step is halved")
      call swap(a1,b1,f1,a0,b0,f0)
      call swap(a2,b2,f2,a0,b0,f0)
      go to 200
  110 if(.not.test1)go to 120
c  test1 means successful convergence
      if(pgword)print 1607,f
 1607 format(1x,5x,"linesearch convergence,fcn=",1x,f20.10)
      trial=a4
      ind=0
      return
  120 if(.not.test3)go to 180
c  current slope positive
      if(pgword)print 1418
 1418 format(1x,5x,"current slope is pos")
      ips=1
      call swap(a3,b3,f3,a4,b4,f4)
      z=(3.d0/(a4-a2))*(f2-f4)+b2+b4
      if(abs(b2+b4+2d0*z).lt.1d-5) go to 140
      w1=z*z-b2*b4
      if(w1.lt.0d0) go to 140
      w1=dsqrt(w1)
      a4=a2+(1.d0-((b4+w1-z)/(b4-b2+2.d0*w1)))*(a4-a2)
       go to 200
  140 a4=a2-b2*((a4-a2)/(b4-b2))
       go to 200
  150 a4=a4-b4*((a3-a4)/(b3-b4))
       go to 200
  160 if(.not.test5)go to 180
      z=(3.d0/(a2-a1))*(f1-f2)+b1+b2
      if(dabs(b1+b2+2.d0*z).lt.1d-5) go to 170
      w1=z*z-b1*b2
      if(w1.le.0d0)go to 170
      w1=dsqrt(w1)
      a4=a1-(1.d0-((b2+w1-z)/(b2-b1+2.d0*w1)))*(a2-a1)
      go to 200
  170 a4=a1-b1*((a2-a1)/(b2-b1))
      go to 200
  180 a4=aa*a4
c  next step now computed
  200 if((nmodel.lt.2).and.(xwrd))go to 40
c check to see if boundary was and is tight
      if(a4.lt.amax)go to 40
      a4=amax
   40 continue
c  no convergence
      ind=2
      return
      end
      subroutine emalg(lambda,mean,sigma2,y,z,

     1      nclass,eps,loglik,temp,npar,aic,bic,
c   stuff for mstep job
c   logical arguments
     1test,wword,abort,pgword,
c   real arguments
     3x,w,s,
     4nmax,nmaxp,nmaxs,
     2nmaxw,nmaxpw,ndimmax,nrep)
      common /cstak/dstak
      double precision dstak(300000000)
      real rstak(600000000)
      integer istak(600000000) 
      common/parblk/nstim,ndim,nsub,nwgt,nmodel,ndat,
     1nsd,nspr,ns,nspace
      common/cnvblk/conv,xsconv,wconv,
     1itmax,itxmax,itwmax
      real x(nmax),w(nmaxw),s(nmaxs)
      logical test,wword,abort,pgword
      real lambda(nclass),y(nsub,nspr)
      real z(nsub,nclass),mean(nclass,nspr)
      double precision oldlgl,loglik,llval(5),diff1,diff5,crit
      double precision temp(nclass),sig2,loglikp
      equivalence(dstak(1),rstak(1))
      equivalence(dstak(1),istak(1))
      save maxit
      data maxit /150 /
c
      lp=i1mach(2)
      iter=1
c
      if(nrep.lt.20)write(lp,10)
   10   format(//' History of iterations'//
     1         ' Iteration',5x,'Log-Likelihood',5x,'Improvement')
c
c  iteration starts here
c
c
c
c$$$ 1980   format(10x,5(f20.10,1h,,2x))
c$$$        do 1979 i=1,nsub
c$$$ 1979      print 1980,(z(i,t),t=1,nclass)
c
c
c
           go to 80
c
   20   continue
c
c  compute likelihood and posterior probabilities for next E step
      call logl(lambda,mean,sigma2,y,z,nclass,nspr,nsub,temp,loglik,
     1  loglikp,npar,aic,bic)
      itind=mod(iter,5)
      if (iter .gt. 1) go to 40
        if(nrep.lt.20)write(lp,30) iter,loglik
   30      format(1x,i5,5d20.7)
        go to 60
c
c  three convergence tests are used:
c     1. Max. no. of iterations
c     2. No negative improvements w.r.t. previous iteration
c     3. Relative improvement over 5 iterations.
   40   continue
        diff1=loglik-oldlgl
      if (iter .ge. 6) goto 50
        if(nrep.lt.20)write(lp,30) iter,loglik,diff1
        if (diff1 .le. 0.d0) go to 70
        go to 60
   50   diff5=loglik-llval(itind)
      crit=eps*dabs(llval(itind))
c      crit=eps*dabs(oldlgl)
      if(nrep.lt.20)write(lp,30) iter,loglik,diff1,diff5,crit
c      write(lp,30) iter,loglik,diff1
      if (diff1 .le. 0.d0 .or. diff5 .le. crit) go to 70
c      if (diff1 .le. 0.d0 .or. diff1 .le. crit) go to 70
   60   iter=iter+1
      if (iter .gt. maxit) then
        if(nrep.lt.20)write(lp,65) maxit
   65      format(/' The maximum no. of iterations (',i4,
     1      ') was performed without reaching the convergence',
     2      ' criterion')
        go to 70
      endif
      oldlgl=loglik
      llval(itind)=loglik
c
c

c  M-step
c  create memory for m-step
   80   jxs=istkgt(nmax,3) 
       jxw=istkgt(nmaxw,3)
       jdif=istkgt(ndim,3)
       jg=istkgt(nmax,4)
       jgw=istkgt(nmaxw,4)
       nmax2=nmax*nmax
       nmax2w=nmaxw*nmaxw
       jhe=istkgt(nmax2,4)
       jhew=istkgt(nmax2w,4)
       nmaxp=(nmax*(nmax+1))/2
       nmaxpw=(nmaxw*(nmaxw+1))/2
       jahe=istkgt(nmaxp,4)
       jahew=istkgt(nmaxpw,4)
       ja1=istkgt(nmax,4)
       jb1=istkgt(nmax,4)
       jc1=istkgt(nmax,4)  
       ja1w=istkgt(nmaxw,4)
        jb1w=istkgt(nmaxw,4)
        jc1w=istkgt(nmaxw,4)
       jyy=istkgt(nmax,3)
       jyyw=istkgt(nmaxw,3)
       jinf=istkgt(nmax,2)
       jinfw=istkgt(nmaxw,2)
       jghat=istkgt(nmax,4)
       jghatw=istkgt(nmaxw,4)
       jsp=istkgt(nmax2,4)
       jspw=istkgt(nmax2w,4)
       jhld=istkgt(nmax2,4) 
       jhldw=istkgt(nmax2w,4)
       jmean1=istkgt(nclass*nspr,3)
       jft=istkgt(nclass,3)
c
             call job(
c   logical arguments
     1test,wword,abort,pgword,
c   real arguments
     3x,w,s,y,
     5rstak(jdif),dstak(jg),dstak(jhe),dstak(jahe),
     7dstak(ja1),dstak(jb1),dstak(jc1),rstak(jyy),rstak(jxs),
     4nmax,nmaxp,istak(jinf),dstak(jghat),dstak(jsp),dstak(jhld),
     1rstak(jxw),dstak(jgw),dstak(jhew),dstak(jahew),
     8dstak(ja1w),dstak(jb1w),dstak(jc1w),rstak(jyyw),istak(jinfw),
     2dstak(jghatw),dstak(jspw),dstak(jhldw),nmaxw,nmaxpw,
     6z,mean,nclass,sig2,rstak(jmean1),rstak(jft),nrep)    
      call istkrl(27)
      call usigm2(lambda,mean,sigma2,y,z,nsub,nspr,nclass)
      go to 20
c
c  iteration ends here
c
   70   continue
      return
      end
      subroutine genery(iseed,nsub,nspr,nclass,nclp1,lambda,
     1    cumlam,d,y,sigma2)
      real lambda(nclass),cumlam(nclp1),d(nclass,nspr),y(nsub,nspr)
      cumlam(1)=lambda(1)
      if(nclass.lt.2) go to 20
         do 10 i=2,nclass
   10         cumlam(i)=cumlam(i-1)+lambda(i)
   20     cumlam(nclass)=1.0
        do 60 i=1,nsub
           if(nclass .gt.  1) go to 25
               k=1
               go to 40
   25            r=runif(iseed)
               do 30 j=1,nclass
                 if(r .gt. cumlam(j)) go to 30
                 k=j
                 go to 40
   30              continue
         call seterr('genery - latent class unknown',29,1,2)
   40      do 50 j=1,nspr
             ydev=gasdev(iseed)
             y(i,j)=d(k,j)+sqrt(sigma2)*ydev
   50          continue
   60            continue
             return
             end
      subroutine job(
c   logical arguments
     1test,wword,abort,pgword,
c   real arguments
     3x,w,s,d,
     5dif,g,he,ahe,a1,b1,c1,yy,xs,
     4nmax,nmaxp,infree,ghat,sp,hld,
     1xw,gw,hew,ahew,a1w,b1w,c1w,yyw,infw,
     2ghatw,spw,hldw,nmaxw,nmaxpw,
     6z,mean,nclass,flglik,mean1,ft,nrep)    
      common/parblk/nstim,ndim,nsub,nwgt,nmodel,ndat,
     1nsd,nspr,ns,nspace
      common/cnvblk/conv,xsconv,wconv,
     1itmax,itxmax,itwmax
      real x(nmax),w(nmaxw),s(1),d(nsub,nspr),
     1dif(ndim),yy(nmax),xs(nmax)
      double precision flglik,ss,ssn,ssnew,
     2gnorm,gnorm1,gnorm2,flglo,
     3g(nmax),he(nmax,nmax),ahe(nmaxp),
     4a1(nmax),b1(nmax),c1(nmax),ghat(nmax),sp(nmax,nmax),
     5hld(nmax,nmax)
      real mean(nclass,nspr),z(nsub,nclass)
      real mean1(nclass,nspr),ft(nclass)
      real xw(nmaxw),yyw(nmaxw)
      double precision gw(nmaxw),hew(nmaxw,nmaxw),ahew(nmaxpw),
     1a1w(nmaxw),b1w(nmaxw),c1w(nmaxw),ghatw(nmaxw),
     2spw(nmaxw,nmaxw),hldw(nmaxw,nmaxw)
      integer infw(nmaxw)    
      integer infree(nmax)    
      logical test,wword,abort,pgword
      nw=nmaxw
c  initialize various things
        iter=0
        do 71 nc=1,nclass
            ft(nc)=0.0
         do 72 ii=1,nsub
   72     ft(nc)=ft(nc)+z(ii,nc)
   71    continue     
          do 75 nc=1,nclass
            do 76 jk=1,nspr
             mean1(nc,jk)=0.0
             do 79 ii=1,nsub
   79        mean1(nc,jk)=mean1(nc,jk)+z(ii,nc)*d(ii,jk)
         if(ft(nc).ne.0.0)mean1(nc,jk)=mean1(nc,jk)/ft(nc)
   76    continue
   75     continue          
	flglo=-1d10
c   initialize ss and grad
 1122 format(1x,5(1x,f5.1))
       ss=0.d0
      call dimfgh(x,w,s,ft,mean1,ss,g,abort,wword,
     1pgword,he,nmax,dif,nclass,mean)
        if(abort)go to 700
c   initialize loglikelihood
   35    call lglcal(ss,flglik)
         do 814 ii=1,nspace
  814        a1(ii)=dble(x(ii))
        gnorm1=0.d0
        gnorm2=0.d0
        do 821 ii=1,nspace
         gnorm1=gnorm1+g(ii)*g(ii)
         gnorm2=gnorm2+a1(ii)*a1(ii)
  821       continue
        gnorm1=dsqrt(gnorm1)
        gnorm2=dsqrt(gnorm2)
        gnorm=(gnorm1*gnorm2)/  dble(nspace)
      if(pgword)print 1000,iter,flglik,ss,gnorm
 1000 format(//1x,5x,"iter no.",3x,i3,3(3x,3e15.7))
c   irklim is rank of configuration matrix
	irklim=nsd
      do 500 iter=1,itmax
c   update configuration
        call dimstp(irklim,ssnew,x,w,s,
     1ft,mean1,g,trial,abort,
     2wword,pgword,
     4he,a1,b1,c1,dif,yy,xs,
     5ahe,nmax,infree,ghat,sp,hld,mean,nclass,hew,nmaxw)      
      if(abort)go to 700
c   update metric
      if(wword) then
        call wstp(irklim,ssn,x,w,s,ft,mean1,gw,
     1trial,abort,wword,pgword,hew,a1w,
     2b1w,c1w,dif,yyw,xw,
     3ahew,nmaxw,infw,ghatw,spw,
     4hldw,mean,nclass,he,nmax)      
          continue
       else
         continue
       end if 
c   update ss and means
       call dimfgh(x,w,s,ft,mean1,ss,g,abort,
     1     wword,pgword,he,nmax,dif,nclass,mean)
c   update loglikelihood
   45 call lglcal(ss,flglik)
       do 810 ii=1,nspace
  810      a1(ii)=dble(x(ii))
       gnorm1=0.d0
       gnorm2=0.d0
       do 841 ii=1,nspace
          gnorm1=gnorm1+g(ii)*g(ii)
          gnorm2=gnorm2+a1(ii)*a1(ii)
  841      continue
         gnorm1=dsqrt(gnorm1)
         gnorm2=dsqrt(gnorm2)
        gnorm=(gnorm1*gnorm2)/dble(nspace)
      if(pgword)print 1000,iter,flglik,ss,gnorm
      nrcd=iter+1
c   test for convergence
      if(flglo.gt.flglik)go to 500
c      if(sngl((flglik-flglo)).ge.conv)go to 500
      if(sngl((flglik-flglo)).ge.conv*sngl(flglo))go to 500
        go to 600
  500 flglo=flglik
        if(nrep.lt.20)print 5053
 5053  format(/1x,5x,"max no. of m-step iterations")
       go to 700
  600  continue
        if(nrep.lt.20)print 5052
 5052  format(/1x,5x,"mstep convergence")
  700  call center(x,nstim,ndim)
 5050 format(1x,5x,"x after m step")
 5051 format(1x,5(1x,f15.4))
      if(pgword)print 5050
      if(pgword)print 5051,(x(ik),ik=1,nsd)
c
c
c
c$$$      print 1979,ss
c$$$ 1979 format(8x,5hSSE =,f20.10)
c
c
c
      return
      end
c   ---------------------------------------
      subroutine outprn(nmax,nmaxw,nstim,ndim,nmodel,nclass,
     1x,w,s,wword)
c  this subroutine prints the output
      real x(nmax),w(nmaxw),s(nstim)
      logical wword
      lp=i1mach(2)
      write(lp,2000)
 2000 format(//1x,5x,'x')
      do 10 i=1,nstim
   10 print * ,i,(x(nstim*(j-1)+i),j=1,ndim)
      if(wword)then
       write(lp,3000)
 3000 format(//1x,5x,'w')
	ndc=ndim
	if(nmodel.eq.3)ndc=ndc+1
	do 53 k=1,ndc
         do 54 l=1,nclass
           write(lp,4000) k,l,w(l+nclass*(k-1))
   54        continue
   53          continue
 4000      format(1x,5x,"dimension no.:",i2,3x,"class no:",i2,
     1          3x,"weight:",e15.7)         
      else
       continue
      end if 
   87 if(nmodel.lt.2)go to 20
 2001 format(5x,5(1x,f10.3))
      print 2002
 2002 format(1x,5x,"s")
      if(nmodel.eq.4)go to 20
      print 2001,(s(i),i=1,nstim)
   20   continue
       if(nmodel.eq.4) then
         do 77 j=1,nclass
        print 2003,j
 2003   format(/1x,5x,"class no.",i3)
          print 2001,(s(i+(j-1)*nstim),i=1,nstim)
   77       continue
        else
          continue
        end if
      return
      end
c  -----------------------------------------------
      subroutine lglcal(ss,flglik)
      common/parblk/nstim,ndim,nsub,nwgt,nmodel,ndat,
     1nsd,nspr,ns,nspace
c  this subroutine computes the loglikelihood from ss
      integer ndat
      double precision ss,flglik
      if(ss.gt.0.0) then
      flglik=ss/dfloat(ndat)
      else
       flglik=1e+10
      end if
      return
      end
c------------------------------------------------
      subroutine rsort(x,n)
c
c  sort x(1),...x(n) in ascending order
c
      real x(n)
c
      m=n
 10   m=m/2
      if (m .eq. 0) return
      k=n-m
      j=1
 20   i=j
 30   im=i+m
      if (x(i) .le. x(im)) goto 40
      temp=x(i)
      x(i)=x(im)
      x(im)=temp
      i=i-m
      if (i .ge. 1) goto 30
 40   j=j+1
      if (j .gt. k) goto 10
      goto 20
      end
      function gasdev(idum)
      data iset/0/
      if(iset.eq.0) then
    1  v1=2.*runif(idum)-1.
       v2=2.*runif(idum)-1.
       r=v1**2+v2**2
       if(r.ge.1)go to 1
        fac=sqrt(-2.*log(r)/r)
         gset=v1*fac
         gasdev=v2*fac
         iset=1
       else
         gasdev=gset
         iset=0
       endif
       return
       end
      subroutine usigm2(lambda,mean,sigma2,y,z,nsub,nspr,nclass)
c
c  update estimate of sigma2 and lambda
c
      real lambda(nclass),mean(nclass,nspr),y(nsub,nspr)
      real z(nsub,nclass)
      double precision sig2,sum
c
      do 50 nc=1,nclass
       lambda(nc)=0.0
      do 60 i=1,nsub
 60      lambda(nc)=lambda(nc)+z(i,nc)
       lambda(nc)=lambda(nc)/float(nsub)
 50    continue
      sig2=0.d0
      do 20 i=1,nsub
      do 20 k=1,nclass
        sum=0.d0
        do 10 j=1,nspr
 10         sum=sum+dble(y(i,j)-mean(k,j))**2
        sig2=sig2+dble(z(i,k))*sum
 20      continue
      sig2=sig2/dfloat(nsub*nspr)
      sigma2=sngl(sig2)
c
      return
      end
      SUBROUTINE ENTER(IRNEW)
C
C  THIS ROUTINE SAVES
C
C    1) THE CURRENT NUMBER OF OUTSTANDING STORAGE ALLOCATIONS, LOUT, AND
C    2) THE CURRENT RECOVERY LEVEL, LRECOV,
C
C  IN AN ENTER-BLOCK IN THE STACK.
C
C  IT ALSO SETS LRECOV = IRNEW IF IRNEW = 1 OR 2.
C  IF IRNEW = 0, THEN THE RECOVERY LEVEL IS NOT ALTERED.
C
C  SCRATCH SPACE ALLOCATED - 3 INTEGER WORDS ARE LEFT ON THE STACK.
C
C  ERROR STATES -
C
C    1 - MUST HAVE IRNEW = 0, 1 OR 2.
C
      COMMON /CSTAK/DSTACK
      DOUBLE PRECISION DSTACK(10000000)
      INTEGER ISTACK(20000000)
      EQUIVALENCE (DSTACK(1),ISTACK(1))
      EQUIVALENCE (ISTACK(1),LOUT)
C
C/6S
      IF (0.GT.IRNEW .OR. IRNEW.GT.2)
     1  CALL SETERR(35HENTER - MUST HAVE IRNEW = 0, 1 OR 2,35,1,2)
C/7S
C     IF (0.GT.IRNEW .OR. IRNEW.GT.2)
C    1  CALL SETERR('ENTER - MUST HAVE IRNEW = 0, 1 OR 2',35,1,2)
C/
C
C  ALLOCATE SPACE FOR SAVING THE ABOVE 2 ITEMS
C  AND A BACK-POINTER FOR CHAINING THE ENTER-BLOCKS TOGETHER.
C
      INOW=ISTKGT(3,2)
C
C  SAVE THE CURRENT NUMBER OF OUTSTANDING ALLOCATIONS.
C
      ISTACK(INOW)=LOUT
C
C  SAVE THE CURRENT RECOVERY LEVEL.
C
      CALL ENTSRC(ISTACK(INOW+1),IRNEW)
C
C  SAVE A BACK-POINTER TO THE START OF THE PREVIOUS ENTER-BLOCK.
C
      ISTACK(INOW+2)=I8TSEL(INOW)
C
      RETURN
C
      END
      SUBROUTINE LEAVE
C
C  THIS ROUTINE
C
C    1) DE-ALLOCATES ALL SCRATCH SPACE ALLOCATED SINCE THE LAST ENTER,
C       INCLUDING THE LAST ENTER-BLOCK.
C    2) RESTORES THE RECOVERY LEVEL TO ITS VALUE
C       AT THE TIME OF THE LAST CALL TO ENTER.
C
C  ERROR STATES -
C
C    1 - CANNOT LEAVE BEYOND THE FIRST ENTER.
C    2 - ISTACK(INOW) HAS BEEN OVERWRITTEN.
C    3 - TOO MANY ISTKRLS OR ISTACK(1 AND/OR INOW) CLOBBERED.
C    4 - ISTACK(INOW+1) HAS BEEN OVERWRITTEN.
C    5 - ISTACK(INOW+2) HAS BEEN OVERWRITTEN.
C
      COMMON /CSTAK/DSTACK
      DOUBLE PRECISION DSTACK(10000000)
      INTEGER ISTACK(20000000)
      EQUIVALENCE (DSTACK(1),ISTACK(1))
      EQUIVALENCE (ISTACK(1),LOUT)
C
C  GET THE POINTER TO THE CURRENT ENTER-BLOCK.
C
      INOW=I8TSEL(-1)
C
C/6S
      IF (INOW.EQ.0)
     1  CALL SETERR(43HLEAVE - CANNOT LEAVE BEYOND THE FIRST ENTER,43,
     2              1,2)
      IF (ISTACK(INOW).LT.1)
     1  CALL SETERR(41HLEAVE - ISTACK(INOW) HAS BEEN OVERWRITTEN,41,2,2)
      IF (LOUT.LT.ISTACK(INOW)) CALL SETERR(
     1  59HLEAVE - TOO MANY ISTKRLS OR ISTACK(1 AND/OR INOW) CLOBBERED,
     2  59,3,2)
      IF (ISTACK(INOW+1).LT.1 .OR. ISTACK(INOW+1).GT.2)
     1  CALL SETERR(43HLEAVE - ISTACK(INOW+1) HAS BEEN OVERWRITTEN,
     2              43,4,2)
      IF (ISTACK(INOW+2).GT.INOW-3 .OR. ISTACK(INOW+2).LT.0)
     1  CALL SETERR(43HLEAVE - ISTACK(INOW+2) HAS BEEN OVERWRITTEN,
     2              43,5,2)
C/7S
C     IF (INOW.EQ.0)
C    1  CALL SETERR('LEAVE - CANNOT LEAVE BEYOND THE FIRST ENTER',43,
C    2              1,2)
C     IF (ISTACK(INOW).LT.1)
C    1  CALL SETERR('LEAVE - ISTACK(INOW) HAS BEEN OVERWRITTEN',41,2,2)
C     IF (LOUT.LT.ISTACK(INOW)) CALL SETERR(
C    1  'LEAVE - TOO MANY ISTKRLS OR ISTACK(1 AND/OR INOW) CLOBBERED',
C    2  59,3,2)
C     IF (ISTACK(INOW+1).LT.1 .OR. ISTACK(INOW+1).GT.2)
C    1  CALL SETERR('LEAVE - ISTACK(INOW+1) HAS BEEN OVERWRITTEN',
C    2              43,4,2)
C     IF (ISTACK(INOW+2).GT.INOW-3 .OR. ISTACK(INOW+2).LT.0)
C    1  CALL SETERR('LEAVE - ISTACK(INOW+2) HAS BEEN OVERWRITTEN',
C    2              43,5,2)
C/
C
C  DE-ALLOCATE THE SCRATCH SPACE.
C
      CALL ISTKRL(LOUT-ISTACK(INOW)+1)
C
C  RESTORE THE RECOVERY LEVEL.
C
      CALL RETSRC(ISTACK(INOW+1))
C
C  LOWER THE BACK-POINTER.
C
      ITEMP=I8TSEL(ISTACK(INOW+2))
C
      RETURN
C
      END
      SUBROUTINE ENTSRC(IROLD,IRNEW)
C
C  THIS ROUTINE RETURNS IROLD = LRECOV AND SETS LRECOV = IRNEW.
C
C  IF THERE IS AN ACTIVE ERROR STATE, THE MESSAGE IS PRINTED
C  AND EXECUTION STOPS.
C
C  IRNEW = 0 LEAVES LRECOV UNCHANGED, WHILE
C  IRNEW = 1 GIVES RECOVERY AND
C  IRNEW = 2 TURNS RECOVERY OFF.
C
C  ERROR STATES -
C
C    1 - ILLEGAL VALUE OF IRNEW.
C    2 - CALLED WHILE IN AN ERROR STATE.
C
C/6S
      IF (IRNEW.LT.0 .OR. IRNEW.GT.2)
     1   CALL SETERR(31HENTSRC - ILLEGAL VALUE OF IRNEW,31,1,2)
C/7S
C     IF (IRNEW.LT.0 .OR. IRNEW.GT.2)
C    1   CALL SETERR('ENTSRC - ILLEGAL VALUE OF IRNEW',31,1,2)
C/
C
      IROLD=I8SAVE(2,IRNEW,IRNEW.NE.0)
C
C  IF HAVE AN ERROR STATE, STOP EXECUTION.
C
C/6S
      IF (I8SAVE(1,0,.FALSE.) .NE. 0) CALL SETERR
     1   (39HENTSRC - CALLED WHILE IN AN ERROR STATE,39,2,2)
C/7S
C     IF (I8SAVE(1,0,.FALSE.) .NE. 0) CALL SETERR
C    1   ('ENTSRC - CALLED WHILE IN AN ERROR STATE',39,2,2)
C/
C
      RETURN
C
      END
      SUBROUTINE RETSRC(IROLD)
C
C  THIS ROUTINE SETS LRECOV = IROLD.
C
C  IF THE CURRENT ERROR BECOMES UNRECOVERABLE,
C  THE MESSAGE IS PRINTED AND EXECUTION STOPS.
C
C  ERROR STATES -
C
C    1 - ILLEGAL VALUE OF IROLD.
C
C/6S
      IF (IROLD.LT.1 .OR. IROLD.GT.2)
     1  CALL SETERR(31HRETSRC - ILLEGAL VALUE OF IROLD,31,1,2)
C/7S
C     IF (IROLD.LT.1 .OR. IROLD.GT.2)
C    1  CALL SETERR('RETSRC - ILLEGAL VALUE OF IROLD',31,1,2)
C/
C
      ITEMP=I8SAVE(2,IROLD,.TRUE.)
C
C  IF THE CURRENT ERROR IS NOW UNRECOVERABLE, PRINT AND STOP.
C
      IF (IROLD.EQ.1 .OR. I8SAVE(1,0,.FALSE.).EQ.0) RETURN
C
        CALL EPRINT
        STOP
C
      END
      INTEGER FUNCTION I8TSEL(INOW)
C
C  TO RETURN I8TSEL = THE POINTER TO THE CURRENT ENTER-BLOCK AND
C  SET THE CURRENT POINTER TO INOW.
C
C  START WITH NO BACK-POINTER.
C
      DATA IENTER/0/
C
      I8TSEL=IENTER
      IF (INOW.GE.0) IENTER=INOW
C
      RETURN
C
      END
      INTEGER FUNCTION ISTKQU(ITYPE)
C
C  RETURNS THE NUMBER OF ITEMS OF TYPE ITYPE THAT REMAIN
C  TO BE ALLOCATED IN ONE REQUEST.
C
C  ERROR STATES -
C
C    1 - LNOW, LUSED, LMAX OR LBOOK OVERWRITTEN
C    2 - ITYPE .LE. 0 .OR. ITYPE .GE. 6
C
      COMMON /CSTAK/DSTAK
C
      double precision dstak(300000000)
      integer istak(600000000)
      INTEGER ISIZE(5)
C
      LOGICAL INIT
C
      EQUIVALENCE (DSTAK(1),ISTAK(1))
      EQUIVALENCE (ISTAK(2),LNOW)
      EQUIVALENCE (ISTAK(3),LUSED)
      EQUIVALENCE (ISTAK(4),LMAX)
      EQUIVALENCE (ISTAK(5),LBOOK)
      EQUIVALENCE (ISTAK(6),ISIZE(1))
C
      DATA INIT/.TRUE./
C
      IF (INIT) CALL I0TK00(INIT,500,4)
C
C/6S
      IF (LNOW.LT.LBOOK.OR.LNOW.GT.LUSED.OR.LUSED.GT.LMAX) CALL SETERR
     1   (47HISTKQU - LNOW, LUSED, LMAX OR LBOOK OVERWRITTEN,
     2    47,1,2)
C/7S
C     IF (LNOW.LT.LBOOK.OR.LNOW.GT.LUSED.OR.LUSED.GT.LMAX) CALL SETERR
C    1   ('ISTKQU - LNOW, LUSED, LMAX OR LBOOK OVERWRITTEN',
C    2    47,1,2)
C/
C
C/6S
      IF (ITYPE.LE.0.OR.ITYPE.GE.6) CALL SETERR
     1   (33HISTKQU - ITYPE.LE.0.OR.ITYPE.GE.6,33,2,2)
C/7S
C     IF (ITYPE.LE.0.OR.ITYPE.GE.6) CALL SETERR
C    1   ('ISTKQU - ITYPE.LE.0.OR.ITYPE.GE.6',33,2,2)
C/
C
      ISTKQU = MAX0( ((LMAX-2)*ISIZE(2))/ISIZE(ITYPE)
     1             - (LNOW*ISIZE(2)-1)/ISIZE(ITYPE)
     2             - 1, 0 )
C
      RETURN
C
      END
      INTEGER FUNCTION ISTKMD(NITEMS)
C
C  CHANGES THE LENGTH OF THE FRAME AT THE TOP OF THE STACK
C  TO NITEMS.
C
C  ERROR STATES -
C
C    1 - LNOW OVERWRITTEN
C    2 - ISTAK(LNOWO-1) OVERWRITTEN
C
      COMMON /CSTAK/DSTAK
C
      double precision dstak(300000000)
      integer istak(600000000)
C
      EQUIVALENCE (DSTAK(1),ISTAK(1))
      EQUIVALENCE (ISTAK(2),LNOW)
C
      LNOWO = LNOW
      CALL ISTKRL(1)
C
      ITYPE = ISTAK(LNOWO-1)
C
C/6S
      IF (ITYPE.LE.0.OR.ITYPE.GE.6) CALL SETERR
     1   (35HISTKMD - ISTAK(LNOWO-1) OVERWRITTEN,35,1,2)
C/7S
C     IF (ITYPE.LE.0.OR.ITYPE.GE.6) CALL SETERR
C    1   ('ISTKMD - ISTAK(LNOWO-1) OVERWRITTEN',35,1,2)
C/
C
      ISTKMD = ISTKGT(NITEMS,ITYPE)
C
      RETURN
C
      END
      SUBROUTINE ISTKRL(NUMBER)
C
C  DE-ALLOCATES THE LAST (NUMBER) ALLOCATIONS MADE IN THE STACK
C  BY ISTKGT.
C
C  ERROR STATES -
C
C    1 - NUMBER .LT. 0
C    2 - LNOW, LUSED, LMAX OR LBOOK OVERWRITTEN
C    3 - ATTEMPT TO DE-ALLOCATE NON-EXISTENT ALLOCATION
C    4 - THE POINTER AT ISTAK(LNOW) OVERWRITTEN
C
      COMMON /CSTAK/DSTAK
C
      double precision dstak(300000000)
      integer istak(600000000)
      LOGICAL INIT
C
      EQUIVALENCE (DSTAK(1),ISTAK(1))
      EQUIVALENCE (ISTAK(1),LOUT)
      EQUIVALENCE (ISTAK(2),LNOW)
      EQUIVALENCE (ISTAK(3),LUSED)
      EQUIVALENCE (ISTAK(4),LMAX)
      EQUIVALENCE (ISTAK(5),LBOOK)
C
      DATA INIT/.TRUE./
C
      IF (INIT) CALL I0TK00(INIT,500,4)
C
C/6S
      IF (NUMBER.LT.0) CALL SETERR(20HISTKRL - NUMBER.LT.0,20,1,2)
C/7S
C     IF (NUMBER.LT.0) CALL SETERR('ISTKRL - NUMBER.LT.0',20,1,2)
C/
C
C/6S
      IF (LNOW.LT.LBOOK.OR.LNOW.GT.LUSED.OR.LUSED.GT.LMAX) CALL SETERR
     1   (47HISTKRL - LNOW, LUSED, LMAX OR LBOOK OVERWRITTEN,
     2    47,2,2)
C/7S
C     IF (LNOW.LT.LBOOK.OR.LNOW.GT.LUSED.OR.LUSED.GT.LMAX) CALL SETERR
C    1   ('ISTKRL - LNOW, LUSED, LMAX OR LBOOK OVERWRITTEN',
C    2    47,2,2)
C/
C
      IN = NUMBER
 10      IF (IN.EQ.0) RETURN
C
C/6S
         IF (LNOW.LE.LBOOK) CALL SETERR
     1   (55HISTKRL - ATTEMPT TO DE-ALLOCATE NON-EXISTENT ALLOCATION,
     2    55,3,2)
C/7S
C        IF (LNOW.LE.LBOOK) CALL SETERR
C    1   ('ISTKRL - ATTEMPT TO DE-ALLOCATE NON-EXISTENT ALLOCATION',
C    2    55,3,2)
C/
C
C     CHECK TO MAKE SURE THE BACK POINTERS ARE MONOTONE.
C
C/6S
         IF (ISTAK(LNOW).LT.LBOOK.OR.ISTAK(LNOW).GE.LNOW-1) CALL SETERR
     1   (47HISTKRL - THE POINTER AT ISTAK(LNOW) OVERWRITTEN,
     2    47,4,2)
C/7S
C        IF (ISTAK(LNOW).LT.LBOOK.OR.ISTAK(LNOW).GE.LNOW-1) CALL SETERR
C    1   ('ISTKRL - THE POINTER AT ISTAK(LNOW) OVERWRITTEN',
C    2    47,4,2)
C/
C
         LOUT = LOUT-1
         LNOW = ISTAK(LNOW)
         IN = IN-1
         GO TO 10
C
      END
      INTEGER FUNCTION ISTKGT(NITEMS,ITYPE)
C
C  ALLOCATES SPACE OUT OF THE INTEGER ARRAY ISTAK (IN COMMON
C  BLOCK CSTAK) FOR AN ARRAY OF LENGTH NITEMS AND OF TYPE
C  DETERMINED BY ITYPE AS FOLLOWS
C
C    1 - LOGICAL
C    2 - INTEGER
C    3 - REAL
C    4 - DOUBLE PRECISION
C    5 - COMPLEX
C
C  ON RETURN, THE ARRAY WILL OCCUPY
C
C    STAK(ISTKGT), STAK(ISTKGT+1), ..., STAK(ISTKGT-NITEMS+1)
C
C  WHERE STAK IS AN ARRAY OF TYPE ITYPE EQUIVALENCED TO ISTAK.
C
C  (FOR THOSE WANTING TO MAKE MACHINE DEPENDENT MODIFICATIONS
C  TO SUPPORT OTHER TYPES, CODES 6,7,8,9,10,11 AND 12 HAVE
C  BEEN RESERVED FOR 1/4 LOGICAL, 1/2 LOGICAL, 1/4 INTEGER,
C  1/2 INTEGER, QUAD PRECISION, DOUBLE COMPLEX AND QUAD
C  COMPLEX, RESPECTIVELY.)
C
C  THE ALLOCATOR RESERVES THE FIRST TEN INTEGER WORDS OF THE STACK
C  FOR ITS OWN INTERNAL BOOK-KEEPING. THESE ARE INITIALIZED BY
C  THE INITIALIZING SUBPROGRAM I0TK00 UPON THE FIRST CALL
C  TO A SUBPROGRAM IN THE ALLOCATION PACKAGE.
C
C  THE USE OF THE FIRST FIVE WORDS IS DESCRIBED BELOW.
C
C    ISTAK( 1) - LOUT,  THE NUMBER OF CURRENT ALLOCATIONS.
C    ISTAK( 2) - LNOW,  THE CURRENT ACTIVE LENGTH OF THE STACK.
C    ISTAK( 3) - LUSED, THE MAXIMUM VALUE OF ISTAK(2) ACHIEVED.
C    ISTAK( 4) - LMAX,  THE MAXIMUM LENGTH THE STACK.
C    ISTAK( 5) - LBOOK, THE NUMBER OF WORDS USED FOR BOOKEEPING.
C
C  THE NEXT FIVE WORDS CONTAIN INTEGERS DESCRIBING THE AMOUNT
C  OF STORAGE ALLOCATED BY THE FORTRAN SYSTEM TO THE VARIOUS
C  DATA TYPES.  THE UNIT OF MEASUREMENT IS ARBITRARY AND MAY
C  BE WORDS, BYTES OR BITS OR WHATEVER IS CONVENIENT.  THE
C  VALUES CURRENTLY ASSUMED CORRESPOND TO AN ANS FORTRAN
C  ENVIRONMENT.  FOR SOME MINI-COMPUTER SYSTEMS THE VALUES MAY
C  HAVE TO BE CHANGED (SEE I0TK00).
C
C    ISTAK( 6) - THE NUMBER OF UNITS ALLOCATED TO LOGICAL
C    ISTAK( 7) - THE NUMBER OF UNITS ALLOCATED TO INTEGER
C    ISTAK( 8) - THE NUMBER OF UNITS ALLOCATED TO REAL
C    ISTAK( 9) - THE NUMBER OF UNITS ALLOCATED TO DOUBLE PRECISION
C    ISTAK(10) - THE NUMBER OF UNITS ALLOCATED TO COMPLEX
C
C  ERROR STATES -
C
C    1 - NITEMS .LT. 0
C    2 - ITYPE .LE. 0 .OR. ITYPE .GE. 6
C    3 - LNOW, LUSED, LMAX OR LBOOK OVERWRITTEN
C    4 - STACK OVERFLOW
C
      COMMON /CSTAK/DSTAK
C
      double precision dstak(300000000)
      integer istak(600000000)
      INTEGER ISIZE(5)
C
      LOGICAL INIT
C
      EQUIVALENCE (DSTAK(1),ISTAK(1))
      EQUIVALENCE (ISTAK(1),LOUT)
      EQUIVALENCE (ISTAK(2),LNOW)
      EQUIVALENCE (ISTAK(3),LUSED)
      EQUIVALENCE (ISTAK(4),LMAX)
      EQUIVALENCE (ISTAK(5),LBOOK)
      EQUIVALENCE (ISTAK(6),ISIZE(1))
C
      DATA INIT/.TRUE./
C
      IF (INIT) CALL I0TK00(INIT,500,4)
C
C/6S
      IF (NITEMS.LT.0) CALL SETERR(20HISTKGT - NITEMS.LT.0,20,1,2)
C/7S
C     IF (NITEMS.LT.0) CALL SETERR('ISTKGT - NITEMS.LT.0',20,1,2)
C/
C
C/6S
      IF (ITYPE.LE.0 .OR. ITYPE.GE.6) CALL SETERR
     1   (33HISTKGT - ITYPE.LE.0.OR.ITYPE.GE.6,33,2,2)
C/7S
C     IF (ITYPE.LE.0 .OR. ITYPE.GE.6) CALL SETERR
C    1   ('ISTKGT - ITYPE.LE.0.OR.ITYPE.GE.6',33,2,2)
C/
C
C/6S
      IF (LNOW.LT.LBOOK.OR.LNOW.GT.LUSED.OR.LUSED.GT.LMAX) CALL SETERR
     1   (47HISTKGT - LNOW, LUSED, LMAX OR LBOOK OVERWRITTEN,
     2    47,3,2)
C/7S
C     IF (LNOW.LT.LBOOK.OR.LNOW.GT.LUSED.OR.LUSED.GT.LMAX) CALL SETERR
C    1   ('ISTKGT - LNOW, LUSED, LMAX OR LBOOK OVERWRITTEN',
C    2    47,3,2)
C/
C
      ISTKGT = (LNOW*ISIZE(2)-1)/ISIZE(ITYPE) + 2
      I = ( (ISTKGT-1+NITEMS)*ISIZE(ITYPE) - 1 )/ISIZE(2) + 3
C
C  STACK OVERFLOW IS AN UNRECOVERABLE ERROR.
C
C/6S
      IF (I.GT.LMAX) CALL SETERR(69HISTKGT - STACK TOO SHORT. ENLARGE IT
     1 AND CALL ISTKIN IN MAIN PROGRAM.,69,4,2)
C/7S
C     IF (I.GT.LMAX) CALL SETERR('ISTKGT - STACK TOO SHORT. ENLARGE IT A
C    *ND CALL ISTKIN IN MAIN PROGRAM.',69,4,2)
C/
C
C  ISTAK(I-1) CONTAINS THE TYPE FOR THIS ALLOCATION.
C  ISTAK(I  ) CONTAINS A POINTER TO THE END OF THE PREVIOUS
C             ALLOCATION.
C
      ISTAK(I-1) = ITYPE
      ISTAK(I  ) = LNOW
      LOUT = LOUT+1
      LNOW = I
      LUSED = MAX0(LUSED,LNOW)
C
      RETURN
C
      END
      SUBROUTINE ISTKIN(NITEMS,ITYPE)
C
C  INITIALIZES THE STACK ALLOCATOR, SETTING THE LENGTH OF THE STACK.
C
C  ERROR STATES -
C
C    1 - NITEMS .LE. 0
C    2 - ITYPE .LE. 0 .OR. ITYPE .GE. 6
C
      LOGICAL INIT
C
      DATA INIT/.TRUE./
C
C/6S
      IF (NITEMS.LE.0) CALL SETERR(20HISTKIN - NITEMS.LE.0,20,1,2)
C/7S
C     IF (NITEMS.LE.0) CALL SETERR('ISTKIN - NITEMS.LE.0',20,1,2)
C/
C
C/6S
      IF (ITYPE.LE.0.OR.ITYPE.GE.6) CALL SETERR
     1   (33HISTKIN - ITYPE.LE.0.OR.ITYPE.GE.6,33,2,2)
C/7S
C     IF (ITYPE.LE.0.OR.ITYPE.GE.6) CALL SETERR
C    1   ('ISTKIN - ITYPE.LE.0.OR.ITYPE.GE.6',33,2,2)
C/
C
      IF (INIT) CALL I0TK00(INIT,NITEMS,ITYPE)
C
      RETURN
C
      END
      INTEGER FUNCTION ISTKST(NFACT)
C
C  RETURNS CONTROL INFORMATION AS FOLLOWS
C
C  NFACT    ITEM RETURNED
C
C    1         LOUT,  THE NUMBER OF CURRENT ALLOCATIONS
C    2         LNOW,  THE CURRENT ACTIVE LENGTH
C    3         LUSED, THE MAXIMUM USED
C    4         LMAX,  THE MAXIMUM ALLOWED
C
      COMMON /CSTAK/DSTAK
C
      double precision dstak(300000000)
      integer istak(600000000)
      INTEGER ISTATS(4)
      LOGICAL INIT
C
      EQUIVALENCE (DSTAK(1),ISTAK(1))
      EQUIVALENCE (ISTAK(1),ISTATS(1))
C
      DATA INIT/.TRUE./
C
      IF (INIT) CALL I0TK00(INIT,500,4)
C
C/6S
      IF (NFACT.LE.0.OR.NFACT.GE.5) CALL SETERR
     1   (33HISTKST - NFACT.LE.0.OR.NFACT.GE.5,33,1,2)
C/7S
C     IF (NFACT.LE.0.OR.NFACT.GE.5) CALL SETERR
C    1   ('ISTKST - NFACT.LE.0.OR.NFACT.GE.5',33,1,2)
C/
C
      ISTKST = ISTATS(NFACT)
C
      RETURN
C
      END
      REAL FUNCTION R1MACH(I)
C
C  SINGLE-PRECISION MACHINE CONSTANTS
C
C  R1MACH(1) = B**(EMIN-1), THE SMALLEST POSITIVE MAGNITUDE.
C
C  R1MACH(2) = B**EMAX*(1 - B**(-T)), THE LARGEST MAGNITUDE.
C
C  R1MACH(3) = B**(-T), THE SMALLEST RELATIVE SPACING.
C
C  R1MACH(4) = B**(1-T), THE LARGEST RELATIVE SPACING.
C
C  R1MACH(5) = LOG10(B)
C
C  TO ALTER THIS FUNCTION FOR A PARTICULAR ENVIRONMENT,
C  THE DESIRED SET OF DATA STATEMENTS SHOULD BE ACTIVATED BY
C  REMOVING THE C FROM COLUMN 1.
C
C  FOR IEEE-ARITHMETIC MACHINES (BINARY STANDARD), THE FIRST
C  SET OF CONSTANTS BELOW SHOULD BE APPROPRIATE.
C
C  WHERE POSSIBLE, OCTAL OR HEXADECIMAL CONSTANTS HAVE BEEN USED
C  TO SPECIFY THE CONSTANTS EXACTLY WHICH HAS IN SOME CASES
C  REQUIRED THE USE OF EQUIVALENT INTEGER ARRAYS.
C
      INTEGER SMALL(2)
      INTEGER LARGE(2)
      INTEGER RIGHT(2)
      INTEGER DIVER(2)
      INTEGER LOG10(2)
C
      REAL RMACH(5)
C
      EQUIVALENCE (RMACH(1),SMALL(1))
      EQUIVALENCE (RMACH(2),LARGE(1))
      EQUIVALENCE (RMACH(3),RIGHT(1))
      EQUIVALENCE (RMACH(4),DIVER(1))
      EQUIVALENCE (RMACH(5),LOG10(1))
C
C     MACHINE CONSTANTS FOR IEEE ARITHMETIC MACHINES, SUCH AS THE AT&T
C     3B SERIES, MOTOROLA 68000 BASED MACHINES (E.G. SUN 3 AND AT&T
C     PC 7300), AND 8087 BASED MICROS (E.G. IBM PC AND AT&T 6300).
C
       DATA SMALL(1) /     8388608 /
       DATA LARGE(1) /  2139095039 /
       DATA RIGHT(1) /   864026624 /
       DATA DIVER(1) /   872415232 /
       DATA LOG10(1) /  1050288283 /
C
C     MACHINE CONSTANTS FOR AMDAHL MACHINES.
C
C      DATA SMALL(1) /    1048576 /
C      DATA LARGE(1) / 2147483647 /
C      DATA RIGHT(1) /  990904320 /
C      DATA DIVER(1) / 1007681536 /
C      DATA LOG10(1) / 1091781651 /
C
C     MACHINE CONSTANTS FOR THE BURROUGHS 1700 SYSTEM.
C
C      DATA RMACH(1) / Z400800000 /
C      DATA RMACH(2) / Z5FFFFFFFF /
C      DATA RMACH(3) / Z4E9800000 /
C      DATA RMACH(4) / Z4EA800000 /
C      DATA RMACH(5) / Z500E730E8 /
C
C     MACHINE CONSTANTS FOR THE BURROUGHS 5700/6700/7700 SYSTEMS.
C
C      DATA RMACH(1) / O1771000000000000 /
C      DATA RMACH(2) / O0777777777777777 /
C      DATA RMACH(3) / O1311000000000000 /
C      DATA RMACH(4) / O1301000000000000 /
C      DATA RMACH(5) / O1157163034761675 /
C
C     MACHINE CONSTANTS FOR THE CDC 6000/7000 SERIES.
C
C      DATA RMACH(1) / 00014000000000000000B /
C      DATA RMACH(2) / 37767777777777777777B /
C      DATA RMACH(3) / 16404000000000000000B /
C      DATA RMACH(4) / 16414000000000000000B /
C      DATA RMACH(5) / 17164642023241175720B /
C
C     MACHINE CONSTANTS FOR CONVEX C-1.
C
C      DATA RMACH(1) / '00800000'X /
C      DATA RMACH(2) / '7FFFFFFF'X /
C      DATA RMACH(3) / '34800000'X /
C      DATA RMACH(4) / '35000000'X /
C      DATA RMACH(5) / '3F9A209B'X /
C
C     MACHINE CONSTANTS FOR THE CRAY 1, XMP, 2, AND 3.
C
C      DATA RMACH(1) / 200034000000000000000B /
C      DATA RMACH(2) / 577767777777777777776B /
C      DATA RMACH(3) / 377224000000000000000B /
C      DATA RMACH(4) / 377234000000000000000B /
C      DATA RMACH(5) / 377774642023241175720B /
C
C     MACHINE CONSTANTS FOR THE DATA GENERAL ECLIPSE S/200.
C
C     NOTE - IT MAY BE APPROPRIATE TO INCLUDE THE FOLLOWING LINE -
C     STATIC RMACH(5)
C
C      DATA SMALL/20K,0/,LARGE/77777K,177777K/
C      DATA RIGHT/35420K,0/,DIVER/36020K,0/
C      DATA LOG10/40423K,42023K/
C
C     MACHINE CONSTANTS FOR THE HARRIS SLASH 6 AND SLASH 7.
C
C      DATA SMALL(1),SMALL(2) / '20000000, '00000201 /
C      DATA LARGE(1),LARGE(2) / '37777777, '00000177 /
C      DATA RIGHT(1),RIGHT(2) / '20000000, '00000352 /
C      DATA DIVER(1),DIVER(2) / '20000000, '00000353 /
C      DATA LOG10(1),LOG10(2) / '23210115, '00000377 /
C
C     MACHINE CONSTANTS FOR THE HONEYWELL DPS 8/70 SERIES.
C
C      DATA RMACH(1) / O402400000000 /
C      DATA RMACH(2) / O376777777777 /
C      DATA RMACH(3) / O714400000000 /
C      DATA RMACH(4) / O716400000000 /
C      DATA RMACH(5) / O776464202324 /
C
C     MACHINE CONSTANTS FOR THE IBM 360/370 SERIES,
C     THE XEROX SIGMA 5/7/9 AND THE SEL SYSTEMS 85/86.
C
C      DATA RMACH(1) / Z00100000 /
C      DATA RMACH(2) / Z7FFFFFFF /
C      DATA RMACH(3) / Z3B100000 /
C      DATA RMACH(4) / Z3C100000 /
C      DATA RMACH(5) / Z41134413 /
C
C     MACHINE CONSTANTS FOR THE INTERDATA 8/32
C     WITH THE UNIX SYSTEM FORTRAN 77 COMPILER.
C
C     FOR THE INTERDATA FORTRAN VII COMPILER REPLACE
C     THE Z'S SPECIFYING HEX CONSTANTS WITH Y'S.
C
C      DATA RMACH(1) / Z'00100000' /
C      DATA RMACH(2) / Z'7EFFFFFF' /
C      DATA RMACH(3) / Z'3B100000' /
C      DATA RMACH(4) / Z'3C100000' /
C      DATA RMACH(5) / Z'41134413' /
C
C     MACHINE CONSTANTS FOR THE PDP-10 (KA OR KI PROCESSOR).
C
C      DATA RMACH(1) / "000400000000 /
C      DATA RMACH(2) / "377777777777 /
C      DATA RMACH(3) / "146400000000 /
C      DATA RMACH(4) / "147400000000 /
C      DATA RMACH(5) / "177464202324 /
C
C     MACHINE CONSTANTS FOR PDP-11 FORTRANS SUPPORTING
C     32-BIT INTEGERS (EXPRESSED IN INTEGER AND OCTAL).
C
C      DATA SMALL(1) /    8388608 /
C      DATA LARGE(1) / 2147483647 /
C      DATA RIGHT(1) /  880803840 /
C      DATA DIVER(1) /  889192448 /
C      DATA LOG10(1) / 1067065499 /
C
C      DATA RMACH(1) / O00040000000 /
C      DATA RMACH(2) / O17777777777 /
C      DATA RMACH(3) / O06440000000 /
C      DATA RMACH(4) / O06500000000 /
C      DATA RMACH(5) / O07746420233 /
C
C     MACHINE CONSTANTS FOR PDP-11 FORTRANS SUPPORTING
C     16-BIT INTEGERS  (EXPRESSED IN INTEGER AND OCTAL).
C
C      DATA SMALL(1),SMALL(2) /   128,     0 /
C      DATA LARGE(1),LARGE(2) / 32767,    -1 /
C      DATA RIGHT(1),RIGHT(2) / 13440,     0 /
C      DATA DIVER(1),DIVER(2) / 13568,     0 /
C      DATA LOG10(1),LOG10(2) / 16282,  8347 /
C
C      DATA SMALL(1),SMALL(2) / O000200, O000000 /
C      DATA LARGE(1),LARGE(2) / O077777, O177777 /
C      DATA RIGHT(1),RIGHT(2) / O032200, O000000 /
C      DATA DIVER(1),DIVER(2) / O032400, O000000 /
C      DATA LOG10(1),LOG10(2) / O037632, O020233 /
C
C     MACHINE CONSTANTS FOR THE SEQUENT BALANCE 8000.
C
C      DATA SMALL(1) / $00800000 /
C      DATA LARGE(1) / $7F7FFFFF /
C      DATA RIGHT(1) / $33800000 /
C      DATA DIVER(1) / $34000000 /
C      DATA LOG10(1) / $3E9A209B /
C
C     MACHINE CONSTANTS FOR THE UNIVAC 1100 SERIES.
C
C      DATA RMACH(1) / O000400000000 /
C      DATA RMACH(2) / O377777777777 /
C      DATA RMACH(3) / O146400000000 /
C      DATA RMACH(4) / O147400000000 /
C      DATA RMACH(5) / O177464202324 /
C
C     MACHINE CONSTANTS FOR THE VAX UNIX F77 COMPILER.
C
C      DATA SMALL(1) /       128 /
C      DATA LARGE(1) /    -32769 /
C      DATA RIGHT(1) /     13440 /
C      DATA DIVER(1) /     13568 /
C      DATA LOG10(1) / 547045274 /
C
C     MACHINE CONSTANTS FOR THE VAX-11 WITH
C     FORTRAN IV-PLUS COMPILER.
C
C      DATA RMACH(1) / Z00000080 /
C      DATA RMACH(2) / ZFFFF7FFF /
C      DATA RMACH(3) / Z00003480 /
C      DATA RMACH(4) / Z00003500 /
C      DATA RMACH(5) / Z209B3F9A /
C
C     MACHINE CONSTANTS FOR VAX/VMS VERSION 2.2.
C
C      DATA RMACH(1) /       '80'X /
C      DATA RMACH(2) / 'FFFF7FFF'X /
C      DATA RMACH(3) /     '3480'X /
C      DATA RMACH(4) /     '3500'X /
C      DATA RMACH(5) / '209B3F9A'X /
C
C     MACHINE CONSTANTS FOR OPUS 300PM (BASED ON CLIPPER)
C
C      DATA RMACH(1) / 1.40129846E-45 /
C      DATA RMACH(2) / 3.40282346E+38 /
C      DATA RMACH(3) / 5.9604645E-08 /
C      DATA RMACH(4) / 1.1920929E-07 /
C      DATA RMACH(5) / 3.01029996E-01 /
C/6S
      IF (I .LT. 1  .OR.  I .GT. 5)
     1   CALL SETERR(24HR1MACH - I OUT OF BOUNDS,24,1,2)
C/7S
C     IF (I .LT. 1  .OR.  I .GT. 5)
C    1   CALL SETERR('R1MACH - I OUT OF BOUNDS',24,1,2)
C/
C
      R1MACH = RMACH(I)
      RETURN
C
      END
      DOUBLE PRECISION FUNCTION D1MACH(I)
C
C  DOUBLE-PRECISION MACHINE CONSTANTS
C
C  D1MACH( 1) = B**(EMIN-1), THE SMALLEST POSITIVE MAGNITUDE.
C
C  D1MACH( 2) = B**EMAX*(1 - B**(-T)), THE LARGEST MAGNITUDE.
C
C  D1MACH( 3) = B**(-T), THE SMALLEST RELATIVE SPACING.
C
C  D1MACH( 4) = B**(1-T), THE LARGEST RELATIVE SPACING.
C
C  D1MACH( 5) = LOG10(B)
C
C  TO ALTER THIS FUNCTION FOR A PARTICULAR ENVIRONMENT,
C  THE DESIRED SET OF DATA STATEMENTS SHOULD BE ACTIVATED BY
C  REMOVING THE C FROM COLUMN 1.
C
C  FOR IEEE-ARITHMETIC MACHINES (BINARY STANDARD), ONE OF THE FIRST
C  TWO SETS OF CONSTANTS BELOW SHOULD BE APPROPRIATE.
C
C  WHERE POSSIBLE, OCTAL OR HEXADECIMAL CONSTANTS HAVE BEEN USED
C  TO SPECIFY THE CONSTANTS EXACTLY, WHICH HAS IN SOME CASES
C  REQUIRED THE USE OF EQUIVALENT INTEGER ARRAYS.
C
      INTEGER SMALL(4)
      INTEGER LARGE(4)
      INTEGER RIGHT(4)
      INTEGER DIVER(4)
      INTEGER LOG10(4)
C
      DOUBLE PRECISION DMACH(5)
C
      EQUIVALENCE (DMACH(1),SMALL(1))
      EQUIVALENCE (DMACH(2),LARGE(1))
      EQUIVALENCE (DMACH(3),RIGHT(1))
      EQUIVALENCE (DMACH(4),DIVER(1))
      EQUIVALENCE (DMACH(5),LOG10(1))
C
C     MACHINE CONSTANTS FOR IEEE ARITHMETIC MACHINES, SUCH AS THE AT&T
C     3B SERIES AND MOTOROLA 68000 BASED MACHINES (E.G. SUN 3 AND AT&T
C     PC 7300), IN WHICH THE MOST SIGNIFICANT BYTE IS STORED FIRST.
C
C      DATA SMALL(1),SMALL(2) /    1048576,          0 /
C      DATA LARGE(1),LARGE(2) / 2146435071,         -1 /
C      DATA RIGHT(1),RIGHT(2) / 1017118720,          0 /
C      DATA DIVER(1),DIVER(2) / 1018167296,          0 /
C      DATA LOG10(1),LOG10(2) / 1070810131, 1352628735 /
C
C     MACHINE CONSTANTS FOR IEEE ARITHMETIC MACHINES AND 8087-BASED
C     MICROS, SUCH AS THE IBM PC AND AT&T 6300, IN WHICH THE LEAST
C     SIGNIFICANT BYTE IS STORED FIRST.
C
       DATA SMALL(1),SMALL(2) /          0,    1048576 /
       DATA LARGE(1),LARGE(2) /         -1, 2146435071 /
       DATA RIGHT(1),RIGHT(2) /          0, 1017118720 /
       DATA DIVER(1),DIVER(2) /          0, 1018167296 /
       DATA LOG10(1),LOG10(2) / 1352628735, 1070810131 /
C
C     MACHINE CONSTANTS FOR AMDAHL MACHINES.
C
C      DATA SMALL(1),SMALL(2) /    1048576,          0 /
C      DATA LARGE(1),LARGE(2) / 2147483647,         -1 /
C      DATA RIGHT(1),RIGHT(2) /  856686592,          0 /
C      DATA DIVER(1),DIVER(2) /  873463808,          0 /
C      DATA LOG10(1),LOG10(2) / 1091781651, 1352628735 /
C
C     MACHINE CONSTANTS FOR THE BURROUGHS 1700 SYSTEM.
C
C      DATA SMALL(1) / ZC00800000 /
C      DATA SMALL(2) / Z000000000 /
C
C      DATA LARGE(1) / ZDFFFFFFFF /
C      DATA LARGE(2) / ZFFFFFFFFF /
C
C      DATA RIGHT(1) / ZCC5800000 /
C      DATA RIGHT(2) / Z000000000 /
C
C      DATA DIVER(1) / ZCC6800000 /
C      DATA DIVER(2) / Z000000000 /
C
C      DATA LOG10(1) / ZD00E730E7 /
C      DATA LOG10(2) / ZC77800DC0 /
C
C     MACHINE CONSTANTS FOR THE BURROUGHS 5700 SYSTEM.
C
C      DATA SMALL(1) / O1771000000000000 /
C      DATA SMALL(2) / O0000000000000000 /
C
C      DATA LARGE(1) / O0777777777777777 /
C      DATA LARGE(2) / O0007777777777777 /
C
C      DATA RIGHT(1) / O1461000000000000 /
C      DATA RIGHT(2) / O0000000000000000 /
C
C      DATA DIVER(1) / O1451000000000000 /
C      DATA DIVER(2) / O0000000000000000 /
C
C      DATA LOG10(1) / O1157163034761674 /
C      DATA LOG10(2) / O0006677466732724 /
C
C     MACHINE CONSTANTS FOR THE BURROUGHS 6700/7700 SYSTEMS.
C
C      DATA SMALL(1) / O1771000000000000 /
C      DATA SMALL(2) / O7770000000000000 /
C
C      DATA LARGE(1) / O0777777777777777 /
C      DATA LARGE(2) / O7777777777777777 /
C
C      DATA RIGHT(1) / O1461000000000000 /
C      DATA RIGHT(2) / O0000000000000000 /
C
C      DATA DIVER(1) / O1451000000000000 /
C      DATA DIVER(2) / O0000000000000000 /
C
C      DATA LOG10(1) / O1157163034761674 /
C      DATA LOG10(2) / O0006677466732724 /
C
C     MACHINE CONSTANTS FOR THE CDC 6000/7000 SERIES.
C
C      DATA SMALL(1) / 00604000000000000000B /
C      DATA SMALL(2) / 00000000000000000000B /
C
C      DATA LARGE(1) / 37767777777777777777B /
C      DATA LARGE(2) / 37167777777777777777B /
C
C      DATA RIGHT(1) / 15604000000000000000B /
C      DATA RIGHT(2) / 15000000000000000000B /
C
C      DATA DIVER(1) / 15614000000000000000B /
C      DATA DIVER(2) / 15010000000000000000B /
C
C      DATA LOG10(1) / 17164642023241175717B /
C      DATA LOG10(2) / 16367571421742254654B /
C
C     MACHINE CONSTANTS FOR CONVEX C-1
C
C      DATA SMALL(1),SMALL(2) / '00100000'X, '00000000'X /
C      DATA LARGE(1),LARGE(2) / '7FFFFFFF'X, 'FFFFFFFF'X /
C      DATA RIGHT(1),RIGHT(2) / '3CC00000'X, '00000000'X /
C      DATA DIVER(1),DIVER(2) / '3CD00000'X, '00000000'X /
C      DATA LOG10(1),LOG10(2) / '3FF34413'X, '509F79FF'X /
C
C     MACHINE CONSTANTS FOR THE CRAY 1, XMP, 2, AND 3.
C
C      DATA SMALL(1) / 201354000000000000000B /
C      DATA SMALL(2) / 000000000000000000000B /
C
C      DATA LARGE(1) / 577767777777777777777B /
C      DATA LARGE(2) / 000007777777777777776B /
C
C      DATA RIGHT(1) / 376434000000000000000B /
C      DATA RIGHT(2) / 000000000000000000000B /
C
C      DATA DIVER(1) / 376444000000000000000B /
C      DATA DIVER(2) / 000000000000000000000B /
C
C      DATA LOG10(1) / 377774642023241175717B /
C      DATA LOG10(2) / 000007571421742254654B /
C
C     MACHINE CONSTANTS FOR THE DATA GENERAL ECLIPSE S/200
C
C     NOTE - IT MAY BE APPROPRIATE TO INCLUDE THE FOLLOWING LINE -
C     STATIC DMACH(5)
C
C      DATA SMALL/20K,3*0/,LARGE/77777K,3*177777K/
C      DATA RIGHT/31420K,3*0/,DIVER/32020K,3*0/
C      DATA LOG10/40423K,42023K,50237K,74776K/
C
C     MACHINE CONSTANTS FOR THE HARRIS SLASH 6 AND SLASH 7
C
C      DATA SMALL(1),SMALL(2) / '20000000, '00000201 /
C      DATA LARGE(1),LARGE(2) / '37777777, '37777577 /
C      DATA RIGHT(1),RIGHT(2) / '20000000, '00000333 /
C      DATA DIVER(1),DIVER(2) / '20000000, '00000334 /
C      DATA LOG10(1),LOG10(2) / '23210115, '10237777 /
C
C     MACHINE CONSTANTS FOR THE HONEYWELL DPS 8/70 SERIES.
C
C      DATA SMALL(1),SMALL(2) / O402400000000, O000000000000 /
C      DATA LARGE(1),LARGE(2) / O376777777777, O777777777777 /
C      DATA RIGHT(1),RIGHT(2) / O604400000000, O000000000000 /
C      DATA DIVER(1),DIVER(2) / O606400000000, O000000000000 /
C      DATA LOG10(1),LOG10(2) / O776464202324, O117571775714 /
C
C     MACHINE CONSTANTS FOR THE IBM 360/370 SERIES,
C     THE XEROX SIGMA 5/7/9 AND THE SEL SYSTEMS 85/86.
C
C      DATA SMALL(1),SMALL(2) / Z00100000, Z00000000 /
C      DATA LARGE(1),LARGE(2) / Z7FFFFFFF, ZFFFFFFFF /
C      DATA RIGHT(1),RIGHT(2) / Z33100000, Z00000000 /
C      DATA DIVER(1),DIVER(2) / Z34100000, Z00000000 /
C      DATA LOG10(1),LOG10(2) / Z41134413, Z509F79FF /
C
C     MACHINE CONSTANTS FOR THE INTERDATA 8/32
C     WITH THE UNIX SYSTEM FORTRAN 77 COMPILER.
C
C     FOR THE INTERDATA FORTRAN VII COMPILER REPLACE
C     THE Z'S SPECIFYING HEX CONSTANTS WITH Y'S.
C
C      DATA SMALL(1),SMALL(2) / Z'00100000', Z'00000000' /
C      DATA LARGE(1),LARGE(2) / Z'7EFFFFFF', Z'FFFFFFFF' /
C      DATA RIGHT(1),RIGHT(2) / Z'33100000', Z'00000000' /
C      DATA DIVER(1),DIVER(2) / Z'34100000', Z'00000000' /
C      DATA LOG10(1),LOG10(2) / Z'41134413', Z'509F79FF' /
C
C     MACHINE CONSTANTS FOR THE PDP-10 (KA PROCESSOR).
C
C      DATA SMALL(1),SMALL(2) / "033400000000, "000000000000 /
C      DATA LARGE(1),LARGE(2) / "377777777777, "344777777777 /
C      DATA RIGHT(1),RIGHT(2) / "113400000000, "000000000000 /
C      DATA DIVER(1),DIVER(2) / "114400000000, "000000000000 /
C      DATA LOG10(1),LOG10(2) / "177464202324, "144117571776 /
C
C     MACHINE CONSTANTS FOR THE PDP-10 (KI PROCESSOR).
C
C      DATA SMALL(1),SMALL(2) / "000400000000, "000000000000 /
C      DATA LARGE(1),LARGE(2) / "377777777777, "377777777777 /
C      DATA RIGHT(1),RIGHT(2) / "103400000000, "000000000000 /
C      DATA DIVER(1),DIVER(2) / "104400000000, "000000000000 /
C      DATA LOG10(1),LOG10(2) / "177464202324, "047674776746 /
C
C     MACHINE CONSTANTS FOR PDP-11 FORTRANS SUPPORTING
C     32-BIT INTEGERS (EXPRESSED IN INTEGER AND OCTAL).
C
C      DATA SMALL(1),SMALL(2) /    8388608,           0 /
C      DATA LARGE(1),LARGE(2) / 2147483647,          -1 /
C      DATA RIGHT(1),RIGHT(2) /  612368384,           0 /
C      DATA DIVER(1),DIVER(2) /  620756992,           0 /
C      DATA LOG10(1),LOG10(2) / 1067065498, -2063872008 /
C
C      DATA SMALL(1),SMALL(2) / O00040000000, O00000000000 /
C      DATA LARGE(1),LARGE(2) / O17777777777, O37777777777 /
C      DATA RIGHT(1),RIGHT(2) / O04440000000, O00000000000 /
C      DATA DIVER(1),DIVER(2) / O04500000000, O00000000000 /
C      DATA LOG10(1),LOG10(2) / O07746420232, O20476747770 /
C
C     MACHINE CONSTANTS FOR PDP-11 FORTRANS SUPPORTING
C     16-BIT INTEGERS (EXPRESSED IN INTEGER AND OCTAL).
C
C      DATA SMALL(1),SMALL(2) /    128,      0 /
C      DATA SMALL(3),SMALL(4) /      0,      0 /
C
C      DATA LARGE(1),LARGE(2) /  32767,     -1 /
C      DATA LARGE(3),LARGE(4) /     -1,     -1 /
C
C      DATA RIGHT(1),RIGHT(2) /   9344,      0 /
C      DATA RIGHT(3),RIGHT(4) /      0,      0 /
C
C      DATA DIVER(1),DIVER(2) /   9472,      0 /
C      DATA DIVER(3),DIVER(4) /      0,      0 /
C
C      DATA LOG10(1),LOG10(2) /  16282,   8346 /
C      DATA LOG10(3),LOG10(4) / -31493, -12296 /
C
C      DATA SMALL(1),SMALL(2) / O000200, O000000 /
C      DATA SMALL(3),SMALL(4) / O000000, O000000 /
C
C      DATA LARGE(1),LARGE(2) / O077777, O177777 /
C      DATA LARGE(3),LARGE(4) / O177777, O177777 /
C
C      DATA RIGHT(1),RIGHT(2) / O022200, O000000 /
C      DATA RIGHT(3),RIGHT(4) / O000000, O000000 /
C
C      DATA DIVER(1),DIVER(2) / O022400, O000000 /
C      DATA DIVER(3),DIVER(4) / O000000, O000000 /
C
C      DATA LOG10(1),LOG10(2) / O037632, O020232 /
C      DATA LOG10(3),LOG10(4) / O102373, O147770 /
C
C     MACHINE CONSTANTS FOR THE PRIME 50 SERIES SYSTEMS
C     WTIH 32-BIT INTEGERS AND 64V MODE INSTRUCTIONS,
C     SUPPLIED BY IGOR BRAY.
C
C      DATA SMALL(1),SMALL(2) / :10000000000, :00000100001 /
C      DATA LARGE(1),LARGE(2) / :17777777777, :37777677775 /
C      DATA RIGHT(1),RIGHT(2) / :10000000000, :00000000122 /
C      DATA DIVER(1),DIVER(2) / :10000000000, :00000000123 /
C      DATA LOG10(1),LOG10(2) / :11504046501, :07674600177 /
C
C     MACHINE CONSTANTS FOR THE SEQUENT BALANCE 8000
C
C      DATA SMALL(1),SMALL(2) / $00000000,  $00100000 /
C      DATA LARGE(1),LARGE(2) / $FFFFFFFF,  $7FEFFFFF /
C      DATA RIGHT(1),RIGHT(2) / $00000000,  $3CA00000 /
C      DATA DIVER(1),DIVER(2) / $00000000,  $3CB00000 /
C      DATA LOG10(1),LOG10(2) / $509F79FF,  $3FD34413 /
C
C     MACHINE CONSTANTS FOR THE UNIVAC 1100 SERIES.
C
C      DATA SMALL(1),SMALL(2) / O000040000000, O000000000000 /
C      DATA LARGE(1),LARGE(2) / O377777777777, O777777777777 /
C      DATA RIGHT(1),RIGHT(2) / O170540000000, O000000000000 /
C      DATA DIVER(1),DIVER(2) / O170640000000, O000000000000 /
C      DATA LOG10(1),LOG10(2) / O177746420232, O411757177572 /
C
C     MACHINE CONSTANTS FOR THE VAX UNIX F77 COMPILER
C
C      DATA SMALL(1),SMALL(2) /        128,           0 /
C      DATA LARGE(1),LARGE(2) /     -32769,          -1 /
C      DATA RIGHT(1),RIGHT(2) /       9344,           0 /
C      DATA DIVER(1),DIVER(2) /       9472,           0 /
C      DATA LOG10(1),LOG10(2) /  546979738,  -805796613 /
C
C     MACHINE CONSTANTS FOR THE VAX-11 WITH
C     FORTRAN IV-PLUS COMPILER
C
C      DATA SMALL(1),SMALL(2) / Z00000080, Z00000000 /
C      DATA LARGE(1),LARGE(2) / ZFFFF7FFF, ZFFFFFFFF /
C      DATA RIGHT(1),RIGHT(2) / Z00002480, Z00000000 /
C      DATA DIVER(1),DIVER(2) / Z00002500, Z00000000 /
C      DATA LOG10(1),LOG10(2) / Z209A3F9A, ZCFF884FB /
C
C     MACHINE CONSTANTS FOR VAX/VMS VERSION 2.2
C
C      DATA SMALL(1),SMALL(2) /       '80'X,        '0'X /
C      DATA LARGE(1),LARGE(2) / 'FFFF7FFF'X, 'FFFFFFFF'X /
C      DATA RIGHT(1),RIGHT(2) /     '2480'X,        '0'X /
C      DATA DIVER(1),DIVER(2) /     '2500'X,        '0'X /
C      DATA LOG10(1),LOG10(2) / '209A3F9A'X, 'CFF884FB'X /
C
C     MACHINE CONSTANTS FOR OPUS 300PM (BASED ON CLIPPER)
C 
C      DATA DMACH(1) / 4.94065645841246544D-324 /
C      DATA DMACH(2) / 1.79769313486231470D+308 /
C      DATA DMACH(3) / 1.110223024625157D-16 /
C      DATA DMACH(4) / 2.220446049250318E-16 /
C      DATA DMACH(5) / 3.010299956639812E-01 /
C
C/6S
      IF (I .LT. 1  .OR.  I .GT. 5)
     1   CALL SETERR(24HD1MACH - I OUT OF BOUNDS,24,1,2)
C/7S
C     IF (I .LT. 1  .OR.  I .GT. 5)
C    1   CALL SETERR('D1MACH - I OUT OF BOUNDS',24,1,2)
C/
C
      D1MACH = DMACH(I)
      RETURN
C
      END
      SUBROUTINE N5ERR(MESSG, NMESSG, NERR, IOPT)
      INTEGER MESSG(1), NMESSG, NERR, IOPT
C  N5ERR IS A PROCEDURE USED TO REDEFINE AN ERROR STATE.
      CALL ERROFF
      CALL SETERR(MESSG, NMESSG, NERR, IOPT)
      RETURN
      END
      INTEGER FUNCTION NERROR(NERR)
C
C  RETURNS NERROR = NERR = THE VALUE OF THE ERROR FLAG LERROR.
C
      NERROR=I8SAVE(1,0,.FALSE.)
      NERR=NERROR
      RETURN
C
      END
      SUBROUTINE ERROFF
C
C  TURNS OFF THE ERROR STATE OFF BY SETTING LERROR=0.
C
      I=I8SAVE(1,0,.TRUE.)
      RETURN
C
      END
      SUBROUTINE SETERR(MESSG,NMESSG,NERR,IOPT)
C
C  SETERR SETS LERROR = NERR, OPTIONALLY PRINTS THE MESSAGE AND DUMPS
C  ACCORDING TO THE FOLLOWING RULES...
C
C    IF IOPT = 1 AND RECOVERING      - JUST REMEMBER THE ERROR.
C    IF IOPT = 1 AND NOT RECOVERING  - PRINT AND STOP.
C    IF IOPT = 2                     - PRINT, DUMP AND STOP.
C
C  INPUT
C
C    MESSG  - THE ERROR MESSAGE.
C    NMESSG - THE LENGTH OF THE MESSAGE, IN CHARACTERS.
C    NERR   - THE ERROR NUMBER. MUST HAVE NERR NON-ZERO.
C    IOPT   - THE OPTION. MUST HAVE IOPT=1 OR 2.
C
C  ERROR STATES -
C
C    1 - MESSAGE LENGTH NOT POSITIVE.
C    2 - CANNOT HAVE NERR=0.
C    3 - AN UNRECOVERED ERROR FOLLOWED BY ANOTHER ERROR.
C    4 - BAD VALUE FOR IOPT.
C
C  ONLY THE FIRST 72 CHARACTERS OF THE MESSAGE ARE PRINTED.
C
C  THE ERROR HANDLER CALLS A SUBROUTINE NAMED SDUMP TO PRODUCE A
C  SYMBOLIC DUMP.
C
C/6S
      INTEGER MESSG(1)
C/7S
C     CHARACTER*1 MESSG(NMESSG)
C/
C
C  THE UNIT FOR ERROR MESSAGES.
C
      IWUNIT=I1MACH(4)
C
      IF (NMESSG.GE.1) GO TO 10
C
C  A MESSAGE OF NON-POSITIVE LENGTH IS FATAL.
C
        WRITE(IWUNIT,9000)
 9000   FORMAT(52H1ERROR    1 IN SETERR - MESSAGE LENGTH NOT POSITIVE.)
        GO TO 60
C
C  NW IS THE NUMBER OF WORDS THE MESSAGE OCCUPIES.
C  (I1MACH(6) IS THE NUMBER OF CHARACTERS PER WORD.)
C
 10   NW=(MIN0(NMESSG,72)-1)/I1MACH(6)+1
C
      IF (NERR.NE.0) GO TO 20
C
C  CANNOT TURN THE ERROR STATE OFF USING SETERR.
C  (I8SAVE SETS A FATAL ERROR HERE.)
C
        WRITE(IWUNIT,9001)
 9001   FORMAT(42H1ERROR    2 IN SETERR - CANNOT HAVE NERR=0//
     1         34H THE CURRENT ERROR MESSAGE FOLLOWS///)
        CALL E9RINT(MESSG,NW,NERR,.TRUE.)
        ITEMP=I8SAVE(1,1,.TRUE.)
        GO TO 50
C
C  SET LERROR AND TEST FOR A PREVIOUS UNRECOVERED ERROR.
C
 20   IF (I8SAVE(1,NERR,.TRUE.).EQ.0) GO TO 30
C
        WRITE(IWUNIT,9002)
 9002   FORMAT(23H1ERROR    3 IN SETERR -,
     1         48H AN UNRECOVERED ERROR FOLLOWED BY ANOTHER ERROR.//
     2         48H THE PREVIOUS AND CURRENT ERROR MESSAGES FOLLOW.///)
        CALL EPRINT
        CALL E9RINT(MESSG,NW,NERR,.TRUE.)
        GO TO 50
C
C  SAVE THIS MESSAGE IN CASE IT IS NOT RECOVERED FROM PROPERLY.
C
 30   CALL E9RINT(MESSG,NW,NERR,.TRUE.)
C
      IF (IOPT.EQ.1 .OR. IOPT.EQ.2) GO TO 40
C
C  MUST HAVE IOPT = 1 OR 2.
C
        WRITE(IWUNIT,9003)
 9003   FORMAT(42H1ERROR    4 IN SETERR - BAD VALUE FOR IOPT//
     1         34H THE CURRENT ERROR MESSAGE FOLLOWS///)
        GO TO 50
C
C  IF THE ERROR IS FATAL, PRINT, DUMP, AND STOP
C
 40   IF (IOPT.EQ.2) GO TO 50
C
C  HERE THE ERROR IS RECOVERABLE
C
C  IF THE RECOVERY MODE IS IN EFFECT, OK, JUST RETURN
C
      IF (I8SAVE(2,0,.FALSE.).EQ.1) RETURN
C
C  OTHERWISE PRINT AND STOP
C
      CALL EPRINT
      STOP
C
 50   CALL EPRINT
 60   CALL SDUMP
      STOP
C
      END
      SUBROUTINE SDUMP
C   THIS IS THE STANDARD DUMP ROUTINE FOR THE PORT LIBRARY.
C   FIRST IT PROVIDES A FORMATTED DUMP OF THE PORT STACK.
C   THEN IT CALLS THE LOCAL (PREFERABLY SYMBOLIC) DUMP ROUTINE.
C     CALL STKDMP
C     CALL FDUMP
      RETURN
      END
      SUBROUTINE I0TK00(LARG,NITEMS,ITYPE)
C
C  INITIALIZES THE STACK TO NITEMS OF TYPE ITYPE
C
      COMMON /CSTAK/DSTAK
C
      double precision dstak(300000000)
      integer istak(600000000)
      LOGICAL LARG,INIT
      INTEGER ISIZE(5)
C
      EQUIVALENCE (DSTAK(1),ISTAK(1))
      EQUIVALENCE (ISTAK(1),LOUT)
      EQUIVALENCE (ISTAK(2),LNOW)
      EQUIVALENCE (ISTAK(3),LUSED)
      EQUIVALENCE (ISTAK(4),LMAX)
      EQUIVALENCE (ISTAK(5),LBOOK)
      EQUIVALENCE (ISTAK(6),ISIZE(1))
C
      DATA INIT/.FALSE./
C
      LARG = .FALSE.
      IF (INIT) RETURN
C
C  HERE TO INITIALIZE
C
      INIT = .TRUE.
C
C  SET DATA SIZES APPROPRIATE FOR A STANDARD CONFORMING
C  FORTRAN SYSTEM USING THE FORTRAN STORAGE UNIT AS THE
C  MEASURE OF SIZE.
C
C  LOGICAL
      ISIZE(1) = 1
C  INTEGER
      ISIZE(2) = 1
C  REAL
      ISIZE(3) = 1
C  DOUBLE PRECISION
      ISIZE(4) = 2
C  COMPLEX
      ISIZE(5) = 2
C
      LBOOK = 10
      LNOW  = LBOOK
      LUSED = LBOOK
      LMAX  = MAX0( (NITEMS*ISIZE(ITYPE))/ISIZE(2), 12 )
      LOUT  = 0
C
      RETURN
C
      END
      SUBROUTINE EPRINT
C
C  THIS SUBROUTINE PRINTS THE LAST ERROR MESSAGE, IF ANY.
C
C/6S
      INTEGER MESSG(1)
C/7S
C     CHARACTER*1 MESSG(1)
C/
C
      CALL E9RINT(MESSG,1,1,.FALSE.)
      RETURN
C
      END
      SUBROUTINE E9RINT(MESSG,NW,NERR,SAVE)
C
C  THIS ROUTINE STORES THE CURRENT ERROR MESSAGE OR PRINTS THE OLD ONE,
C  IF ANY, DEPENDING ON WHETHER OR NOT SAVE = .TRUE. .
C
C  CHANGED, BY P.FOX, MAY 18, 1983, FROM THE ORIGINAL VERSION IN ORDER
C  TO GET RID OF THE FORTRAN CARRIAGE CONTROL LINE OVERWRITE
C  CHARACTER +, WHICH HAS ALWAYS CAUSED TROUBLE.
C  FOR THE RECORD, THE PREVIOUS VERSION HAD THE FOLLOWING ARRAY
C  AND CALLS -   (WHERE CCPLUS WAS DECLARED OF TYPE INTEGER)
C
C      DATA CCPLUS  / 1H+ /
C
C      DATA FMT( 1) / 1H( /
C      DATA FMT( 2) / 1HA /
C      DATA FMT( 3) / 1H1 /
C      DATA FMT( 4) / 1H, /
C      DATA FMT( 5) / 1H1 /
C      DATA FMT( 6) / 1H4 /
C      DATA FMT( 7) / 1HX /
C      DATA FMT( 8) / 1H, /
C      DATA FMT( 9) / 1H7 /
C      DATA FMT(10) / 1H2 /
C      DATA FMT(11) / 1HA /
C      DATA FMT(12) / 1HX /
C      DATA FMT(13) / 1HX /
C      DATA FMT(14) / 1H) /
C
C        CALL S88FMT(2,I1MACH(6),FMT(12))
C        WRITE(IWUNIT,FMT) CCPLUS,(MESSGP(I),I=1,NWP)
C
C/6S
      INTEGER MESSG(NW)
C/7S
C     CHARACTER*1 MESSG(NW)
C/
      LOGICAL SAVE
C
C  MESSGP STORES AT LEAST THE FIRST 72 CHARACTERS OF THE PREVIOUS
C  MESSAGE. ITS LENGTH IS MACHINE DEPENDENT AND MUST BE AT LEAST
C
C       1 + 71/(THE NUMBER OF CHARACTERS STORED PER INTEGER WORD).
C
C/6S
      INTEGER MESSGP(36),FMT(10), FMT10(10)
      EQUIVALENCE (FMT(1),FMT10(1))
C/7S
C     CHARACTER*1 MESSGP(72),FMT(10)
C     CHARACTER*10 FMT10
C     EQUIVALENCE (FMT(1),FMT10)
C/
C
C  START WITH NO PREVIOUS MESSAGE.
C
C/6S
      DATA MESSGP(1)/1H1/, NWP/0/, NERRP/0/
C/7S
C     DATA MESSGP(1)/'1'/, NWP/0/, NERRP/0/
C/
C
C  SET UP THE FORMAT FOR PRINTING THE ERROR MESSAGE.
C  THE FORMAT IS SIMPLY (A1,14X,72AXX) WHERE XX=I1MACH(6) IS THE
C  NUMBER OF CHARACTERS STORED PER INTEGER WORD.
C
C/6S
      DATA FMT( 1) / 1H( /
      DATA FMT( 2) / 1H3 /
      DATA FMT( 3) / 1HX /
      DATA FMT( 4) / 1H, /
      DATA FMT( 5) / 1H7 /
      DATA FMT( 6) / 1H2 /
      DATA FMT( 7) / 1HA /
      DATA FMT( 8) / 1HX /
      DATA FMT( 9) / 1HX /
      DATA FMT(10) / 1H) /
C/7S
C     DATA FMT( 1) / '(' /
C     DATA FMT( 2) / '3' /
C     DATA FMT( 3) / 'X' /
C     DATA FMT( 4) / ',' /
C     DATA FMT( 5) / '7' /
C     DATA FMT( 6) / '2' /
C     DATA FMT( 7) / 'A' /
C     DATA FMT( 8) / 'X' /
C     DATA FMT( 9) / 'X' /
C     DATA FMT(10) / ')' /
C/
C
      IF (.NOT.SAVE) GO TO 20
C
C  SAVE THE MESSAGE.
C
        NWP=NW
        NERRP=NERR
        DO 10 I=1,NW
 10     MESSGP(I)=MESSG(I)
C
        GO TO 30
C
 20   IF (I8SAVE(1,0,.FALSE.).EQ.0) GO TO 30
C
C  PRINT THE MESSAGE.
C
        IWUNIT=I1MACH(4)
        WRITE(IWUNIT,9000) NERRP
 9000   FORMAT(7H ERROR ,I4,4H IN )
C
        CALL S88FMT(2,I1MACH(6),FMT( 8))
        WRITE(IWUNIT,FMT10) (MESSGP(I),I=1,NWP)
C
 30   RETURN
C
      END
      INTEGER FUNCTION I1MACH(I)
C
C  I/O UNIT NUMBERS.
C
C    I1MACH( 1) = THE STANDARD INPUT UNIT.
C
C    I1MACH( 2) = THE STANDARD OUTPUT UNIT.
C
C    I1MACH( 3) = THE STANDARD PUNCH UNIT.
C
C    I1MACH( 4) = THE STANDARD ERROR MESSAGE UNIT.
C
C  WORDS.
C
C    I1MACH( 5) = THE NUMBER OF BITS PER INTEGER STORAGE UNIT.
C
C    I1MACH( 6) = THE NUMBER OF CHARACTERS PER INTEGER STORAGE UNIT.
C
C  INTEGERS.
C
C    ASSUME INTEGERS ARE REPRESENTED IN THE S-DIGIT, BASE-A FORM
C
C               SIGN ( X(S-1)*A**(S-1) + ... + X(1)*A + X(0) )
C
C               WHERE 0 .LE. X(I) .LT. A FOR I=0,...,S-1.
C
C    I1MACH( 7) = A, THE BASE.
C
C    I1MACH( 8) = S, THE NUMBER OF BASE-A DIGITS.
C
C    I1MACH( 9) = A**S - 1, THE LARGEST MAGNITUDE.
C
C  FLOATING-POINT NUMBERS.
C
C    ASSUME FLOATING-POINT NUMBERS ARE REPRESENTED IN THE T-DIGIT,
C    BASE-B FORM
C
C               SIGN (B**E)*( (X(1)/B) + ... + (X(T)/B**T) )
C
C               WHERE 0 .LE. X(I) .LT. B FOR I=1,...,T,
C               0 .LT. X(1), AND EMIN .LE. E .LE. EMAX.
C
C    I1MACH(10) = B, THE BASE.
C
C  SINGLE-PRECISION
C
C    I1MACH(11) = T, THE NUMBER OF BASE-B DIGITS.
C
C    I1MACH(12) = EMIN, THE SMALLEST EXPONENT E.
C
C    I1MACH(13) = EMAX, THE LARGEST EXPONENT E.
C
C  DOUBLE-PRECISION
C
C    I1MACH(14) = T, THE NUMBER OF BASE-B DIGITS.
C
C    I1MACH(15) = EMIN, THE SMALLEST EXPONENT E.
C
C    I1MACH(16) = EMAX, THE LARGEST EXPONENT E.
C
C  TO ALTER THIS FUNCTION FOR A PARTICULAR ENVIRONMENT,
C  THE DESIRED SET OF DATA STATEMENTS SHOULD BE ACTIVATED BY
C  REMOVING THE C FROM COLUMN 1.  ALSO, THE VALUES OF
C  I1MACH(1) - I1MACH(4) SHOULD BE CHECKED FOR CONSISTENCY
C  WITH THE LOCAL OPERATING SYSTEM.  FOR FORTRAN 77, YOU MAY WISH
C  TO ADJUST THE DATA STATEMENT SO IMACH(6) IS SET TO 1, AND
C  THEN TO COMMENT OUT THE EXECUTABLE TEST ON I .EQ. 6 BELOW.
C
C  FOR IEEE-ARITHMETIC MACHINES (BINARY STANDARD), THE FIRST
C  SET OF CONSTANTS BELOW SHOULD BE APPROPRIATE, EXCEPT PERHAPS
C  FOR IMACH(1) - IMACH(4).
C
      INTEGER IMACH(16),OUTPUT,SANITY
C
      EQUIVALENCE (IMACH(4),OUTPUT)
C
C     MACHINE CONSTANTS FOR IEEE ARITHMETIC MACHINES, SUCH AS THE AT&T
C     3B SERIES, MOTOROLA 68000 BASED MACHINES (E.G. SUN 3 AND AT&T
C     PC 7300), AND 8087 BASED MICROS (E.G. IBM PC AND AT&T 6300).
C
       DATA IMACH( 1) /    5 /
       DATA IMACH( 2) /    6 /
       DATA IMACH( 3) /    7 /
       DATA IMACH( 4) /    6 /
       DATA IMACH( 5) /   64 /
       DATA IMACH( 6) /    4 /
       DATA IMACH( 7) /    2 /
       DATA IMACH( 8) /   31 /
       DATA IMACH( 9) / 2147483647 /
       DATA IMACH(10) /    2 /
       DATA IMACH(11) /   24 /
       DATA IMACH(12) / -125 /
       DATA IMACH(13) /  128 /
       DATA IMACH(14) /   53 /
       DATA IMACH(15) / -1021 /
       DATA IMACH(16) /  1024 /, SANITY/987/
C
C     MACHINE CONSTANTS FOR AMDAHL MACHINES.
C
C      DATA IMACH( 1) /   5 /
C      DATA IMACH( 2) /   6 /
C      DATA IMACH( 3) /   7 /
C      DATA IMACH( 4) /   6 /
C      DATA IMACH( 5) /  32 /
C      DATA IMACH( 6) /   4 /
C      DATA IMACH( 7) /   2 /
C      DATA IMACH( 8) /  31 /
C      DATA IMACH( 9) / 2147483647 /
C      DATA IMACH(10) /  16 /
C      DATA IMACH(11) /   6 /
C      DATA IMACH(12) / -64 /
C      DATA IMACH(13) /  63 /
C      DATA IMACH(14) /  14 /
C      DATA IMACH(15) / -64 /
C      DATA IMACH(16) /  63 /, SANITY/987/
C
C     MACHINE CONSTANTS FOR THE BURROUGHS 1700 SYSTEM.
C
C      DATA IMACH( 1) /    7 /
C      DATA IMACH( 2) /    2 /
C      DATA IMACH( 3) /    2 /
C      DATA IMACH( 4) /    2 /
C      DATA IMACH( 5) /   36 /
C      DATA IMACH( 6) /    4 /
C      DATA IMACH( 7) /    2 /
C      DATA IMACH( 8) /   33 /
C      DATA IMACH( 9) / Z1FFFFFFFF /
C      DATA IMACH(10) /    2 /
C      DATA IMACH(11) /   24 /
C      DATA IMACH(12) / -256 /
C      DATA IMACH(13) /  255 /
C      DATA IMACH(14) /   60 /
C      DATA IMACH(15) / -256 /
C      DATA IMACH(16) /  255 /, SANITY/987/
C
C     MACHINE CONSTANTS FOR THE BURROUGHS 5700 SYSTEM.
C
C      DATA IMACH( 1) /   5 /
C      DATA IMACH( 2) /   6 /
C      DATA IMACH( 3) /   7 /
C      DATA IMACH( 4) /   6 /
C      DATA IMACH( 5) /  48 /
C      DATA IMACH( 6) /   6 /
C      DATA IMACH( 7) /   2 /
C      DATA IMACH( 8) /  39 /
C      DATA IMACH( 9) / O0007777777777777 /
C      DATA IMACH(10) /   8 /
C      DATA IMACH(11) /  13 /
C      DATA IMACH(12) / -50 /
C      DATA IMACH(13) /  76 /
C      DATA IMACH(14) /  26 /
C      DATA IMACH(15) / -50 /
C      DATA IMACH(16) /  76 /, SANITY/987/
C
C     MACHINE CONSTANTS FOR THE BURROUGHS 6700/7700 SYSTEMS.
C
C      DATA IMACH( 1) /   5 /
C      DATA IMACH( 2) /   6 /
C      DATA IMACH( 3) /   7 /
C      DATA IMACH( 4) /   6 /
C      DATA IMACH( 5) /  48 /
C      DATA IMACH( 6) /   6 /
C      DATA IMACH( 7) /   2 /
C      DATA IMACH( 8) /  39 /
C      DATA IMACH( 9) / O0007777777777777 /
C      DATA IMACH(10) /   8 /
C      DATA IMACH(11) /  13 /
C      DATA IMACH(12) / -50 /
C      DATA IMACH(13) /  76 /
C      DATA IMACH(14) /  26 /
C      DATA IMACH(15) / -32754 /
C      DATA IMACH(16) /  32780 /, SANITY/987/
C
C     MACHINE CONSTANTS FOR THE CDC 6000/7000 SERIES.
C
C      DATA IMACH( 1) /    5 /
C      DATA IMACH( 2) /    6 /
C      DATA IMACH( 3) /    7 /
C      DATA IMACH( 4) /    6 /
C      DATA IMACH( 5) /   60 /
C      DATA IMACH( 6) /   10 /
C      DATA IMACH( 7) /    2 /
C      DATA IMACH( 8) /   48 /
C      DATA IMACH( 9) / 00007777777777777777B /
C      DATA IMACH(10) /    2 /
C      DATA IMACH(11) /   48 /
C      DATA IMACH(12) / -974 /
C      DATA IMACH(13) / 1070 /
C      DATA IMACH(14) /   96 /
C      DATA IMACH(15) / -927 /
C      DATA IMACH(16) / 1070 /, SANITY/987/
C
C     MACHINE CONSTANTS FOR CONVEX C-1.
C
C      DATA IMACH( 1) /    5 /
C      DATA IMACH( 2) /    6 /
C      DATA IMACH( 3) /    7 /
C      DATA IMACH( 4) /    6 /
C      DATA IMACH( 5) /   32 /
C      DATA IMACH( 6) /    4 /
C      DATA IMACH( 7) /    2 /
C      DATA IMACH( 8) /   31 /
C      DATA IMACH( 9) / 2147483647 /
C      DATA IMACH(10) /    2 /
C      DATA IMACH(11) /   24 /
C      DATA IMACH(12) / -128 /
C      DATA IMACH(13) /  127 /
C      DATA IMACH(14) /   53 /
C      DATA IMACH(15) /-1024 /
C      DATA IMACH(16) / 1023 /, SANITY/987/
C
C     MACHINE CONSTANTS FOR THE CRAY 1, XMP, 2, AND 3.
C
C      DATA IMACH( 1) /     5 /
C      DATA IMACH( 2) /     6 /
C      DATA IMACH( 3) /   102 /
C      DATA IMACH( 4) /     6 /
C      DATA IMACH( 5) /    64 /
C      DATA IMACH( 6) /     8 /
C      DATA IMACH( 7) /     2 /
C      DATA IMACH( 8) /    63 /
C      DATA IMACH( 9) /  777777777777777777777B /
C      DATA IMACH(10) /     2 /
C      DATA IMACH(11) /    47 /
C      DATA IMACH(12) / -8189 /
C      DATA IMACH(13) /  8190 /
C      DATA IMACH(14) /    94 /
C      DATA IMACH(15) / -8099 /
C      DATA IMACH(16) /  8190 /, SANITY/987/
C
C     MACHINE CONSTANTS FOR THE DATA GENERAL ECLIPSE S/200.
C
C      DATA IMACH( 1) /   11 /
C      DATA IMACH( 2) /   12 /
C      DATA IMACH( 3) /    8 /
C      DATA IMACH( 4) /   10 /
C      DATA IMACH( 5) /   16 /
C      DATA IMACH( 6) /    2 /
C      DATA IMACH( 7) /    2 /
C      DATA IMACH( 8) /   15 /
C      DATA IMACH( 9) /32767 /
C      DATA IMACH(10) /   16 /
C      DATA IMACH(11) /    6 /
C      DATA IMACH(12) /  -64 /
C      DATA IMACH(13) /   63 /
C      DATA IMACH(14) /   14 /
C      DATA IMACH(15) /  -64 /
C      DATA IMACH(16) /   63 /, SANITY/987/
C
C     MACHINE CONSTANTS FOR THE HARRIS SLASH 6 AND SLASH 7.
C
C      DATA IMACH( 1) /       5 /
C      DATA IMACH( 2) /       6 /
C      DATA IMACH( 3) /       0 /
C      DATA IMACH( 4) /       6 /
C      DATA IMACH( 5) /      24 /
C      DATA IMACH( 6) /       3 /
C      DATA IMACH( 7) /       2 /
C      DATA IMACH( 8) /      23 /
C      DATA IMACH( 9) / 8388607 /
C      DATA IMACH(10) /       2 /
C      DATA IMACH(11) /      23 /
C      DATA IMACH(12) /    -127 /
C      DATA IMACH(13) /     127 /
C      DATA IMACH(14) /      38 /
C      DATA IMACH(15) /    -127 /
C      DATA IMACH(16) /     127 /, SANITY/987/
C
C     MACHINE CONSTANTS FOR THE HONEYWELL DPS 8/70 SERIES.
C
C      DATA IMACH( 1) /    5 /
C      DATA IMACH( 2) /    6 /
C      DATA IMACH( 3) /   43 /
C      DATA IMACH( 4) /    6 /
C      DATA IMACH( 5) /   36 /
C      DATA IMACH( 6) /    4 /
C      DATA IMACH( 7) /    2 /
C      DATA IMACH( 8) /   35 /
C      DATA IMACH( 9) / O377777777777 /
C      DATA IMACH(10) /    2 /
C      DATA IMACH(11) /   27 /
C      DATA IMACH(12) / -127 /
C      DATA IMACH(13) /  127 /
C      DATA IMACH(14) /   63 /
C      DATA IMACH(15) / -127 /
C      DATA IMACH(16) /  127 /, SANITY/987/
C
C     MACHINE CONSTANTS FOR THE IBM 360/370 SERIES,
C     THE XEROX SIGMA 5/7/9 AND THE SEL SYSTEMS 85/86.
C
C      DATA IMACH( 1) /   5 /
C      DATA IMACH( 2) /   6 /
C      DATA IMACH( 3) /   7 /
C      DATA IMACH( 4) /   6 /
C      DATA IMACH( 5) /  32 /
C      DATA IMACH( 6) /   4 /
C      DATA IMACH( 7) /   2 /
C      DATA IMACH( 8) /  31 /
C      DATA IMACH( 9) / Z7FFFFFFF /
C      DATA IMACH(10) /  16 /
C      DATA IMACH(11) /   6 /
C      DATA IMACH(12) / -64 /
C      DATA IMACH(13) /  63 /
C      DATA IMACH(14) /  14 /
C      DATA IMACH(15) / -64 /
C      DATA IMACH(16) /  63 /, SANITY/987/
C
C     MACHINE CONSTANTS FOR THE INTERDATA 8/32
C     WITH THE UNIX SYSTEM FORTRAN 77 COMPILER.
C
C     FOR THE INTERDATA FORTRAN VII COMPILER REPLACE
C     THE Z'S SPECIFYING HEX CONSTANTS WITH Y'S.
C
C      DATA IMACH( 1) /   5 /
C      DATA IMACH( 2) /   6 /
C      DATA IMACH( 3) /   6 /
C      DATA IMACH( 4) /   6 /
C      DATA IMACH( 5) /  32 /
C      DATA IMACH( 6) /   4 /
C      DATA IMACH( 7) /   2 /
C      DATA IMACH( 8) /  31 /
C      DATA IMACH( 9) / Z'7FFFFFFF' /
C      DATA IMACH(10) /  16 /
C      DATA IMACH(11) /   6 /
C      DATA IMACH(12) / -64 /
C      DATA IMACH(13) /  62 /
C      DATA IMACH(14) /  14 /
C      DATA IMACH(15) / -64 /
C      DATA IMACH(16) /  62 /, SANITY/987/
C
C     MACHINE CONSTANTS FOR THE PDP-10 (KA PROCESSOR).
C
C      DATA IMACH( 1) /    5 /
C      DATA IMACH( 2) /    6 /
C      DATA IMACH( 3) /    7 /
C      DATA IMACH( 4) /    6 /
C      DATA IMACH( 5) /   36 /
C      DATA IMACH( 6) /    5 /
C      DATA IMACH( 7) /    2 /
C      DATA IMACH( 8) /   35 /
C      DATA IMACH( 9) / "377777777777 /
C      DATA IMACH(10) /    2 /
C      DATA IMACH(11) /   27 /
C      DATA IMACH(12) / -128 /
C      DATA IMACH(13) /  127 /
C      DATA IMACH(14) /   54 /
C      DATA IMACH(15) / -101 /
C      DATA IMACH(16) /  127 /, SANITY/987/
C
C     MACHINE CONSTANTS FOR THE PDP-10 (KI PROCESSOR).
C
C      DATA IMACH( 1) /    5 /
C      DATA IMACH( 2) /    6 /
C      DATA IMACH( 3) /    7 /
C      DATA IMACH( 4) /    6 /
C      DATA IMACH( 5) /   36 /
C      DATA IMACH( 6) /    5 /
C      DATA IMACH( 7) /    2 /
C      DATA IMACH( 8) /   35 /
C      DATA IMACH( 9) / "377777777777 /
C      DATA IMACH(10) /    2 /
C      DATA IMACH(11) /   27 /
C      DATA IMACH(12) / -128 /
C      DATA IMACH(13) /  127 /
C      DATA IMACH(14) /   62 /
C      DATA IMACH(15) / -128 /
C      DATA IMACH(16) /  127 /, SANITY/987/
C
C     MACHINE CONSTANTS FOR PDP-11 FORTRANS SUPPORTING
C     32-BIT INTEGER ARITHMETIC.
C
C      DATA IMACH( 1) /    5 /
C      DATA IMACH( 2) /    6 /
C      DATA IMACH( 3) /    7 /
C      DATA IMACH( 4) /    6 /
C      DATA IMACH( 5) /   32 /
C      DATA IMACH( 6) /    4 /
C      DATA IMACH( 7) /    2 /
C      DATA IMACH( 8) /   31 /
C      DATA IMACH( 9) / 2147483647 /
C      DATA IMACH(10) /    2 /
C      DATA IMACH(11) /   24 /
C      DATA IMACH(12) / -127 /
C      DATA IMACH(13) /  127 /
C      DATA IMACH(14) /   56 /
C      DATA IMACH(15) / -127 /
C      DATA IMACH(16) /  127 /, SANITY/987/
C
C     MACHINE CONSTANTS FOR PDP-11 FORTRANS SUPPORTING
C     16-BIT INTEGER ARITHMETIC.
C
C      DATA IMACH( 1) /    5 /
C      DATA IMACH( 2) /    6 /
C      DATA IMACH( 3) /    7 /
C      DATA IMACH( 4) /    6 /
C      DATA IMACH( 5) /   16 /
C      DATA IMACH( 6) /    2 /
C      DATA IMACH( 7) /    2 /
C      DATA IMACH( 8) /   15 /
C      DATA IMACH( 9) / 32767 /
C      DATA IMACH(10) /    2 /
C      DATA IMACH(11) /   24 /
C      DATA IMACH(12) / -127 /
C      DATA IMACH(13) /  127 /
C      DATA IMACH(14) /   56 /
C      DATA IMACH(15) / -127 /
C      DATA IMACH(16) /  127 /, SANITY/987/
C
C     MACHINE CONSTANTS FOR THE PRIME 50 SERIES SYSTEMS
C     WTIH 32-BIT INTEGERS AND 64V MODE INSTRUCTIONS,
C     SUPPLIED BY IGOR BRAY.
C
C      DATA IMACH( 1) /            1 /
C      DATA IMACH( 2) /            1 /
C      DATA IMACH( 3) /            2 /
C      DATA IMACH( 4) /            1 /
C      DATA IMACH( 5) /           32 /
C      DATA IMACH( 6) /            4 /
C      DATA IMACH( 7) /            2 /
C      DATA IMACH( 8) /           31 /
C      DATA IMACH( 9) / :17777777777 /
C      DATA IMACH(10) /            2 /
C      DATA IMACH(11) /           23 /
C      DATA IMACH(12) /         -127 /
C      DATA IMACH(13) /         +127 /
C      DATA IMACH(14) /           47 /
C      DATA IMACH(15) /       -32895 /
C      DATA IMACH(16) /       +32637 /, SANITY/987/
C
C     MACHINE CONSTANTS FOR THE SEQUENT BALANCE 8000.
C
C      DATA IMACH( 1) /     0 /
C      DATA IMACH( 2) /     0 /
C      DATA IMACH( 3) /     7 /
C      DATA IMACH( 4) /     0 /
C      DATA IMACH( 5) /    32 /
C      DATA IMACH( 6) /     1 /
C      DATA IMACH( 7) /     2 /
C      DATA IMACH( 8) /    31 /
C      DATA IMACH( 9) /  2147483647 /
C      DATA IMACH(10) /     2 /
C      DATA IMACH(11) /    24 /
C      DATA IMACH(12) /  -125 /
C      DATA IMACH(13) /   128 /
C      DATA IMACH(14) /    53 /
C      DATA IMACH(15) / -1021 /
C      DATA IMACH(16) /  1024 /, SANITY/987/
C
C     MACHINE CONSTANTS FOR THE UNIVAC 1100 SERIES.
C
C     NOTE THAT THE PUNCH UNIT, I1MACH(3), HAS BEEN SET TO 7
C     WHICH IS APPROPRIATE FOR THE UNIVAC-FOR SYSTEM.
C     IF YOU HAVE THE UNIVAC-FTN SYSTEM, SET IT TO 1.
C
C      DATA IMACH( 1) /    5 /
C      DATA IMACH( 2) /    6 /
C      DATA IMACH( 3) /    7 /
C      DATA IMACH( 4) /    6 /
C      DATA IMACH( 5) /   36 /
C      DATA IMACH( 6) /    6 /
C      DATA IMACH( 7) /    2 /
C      DATA IMACH( 8) /   35 /
C      DATA IMACH( 9) / O377777777777 /
C      DATA IMACH(10) /    2 /
C      DATA IMACH(11) /   27 /
C      DATA IMACH(12) / -128 /
C      DATA IMACH(13) /  127 /
C      DATA IMACH(14) /   60 /
C      DATA IMACH(15) /-1024 /
C      DATA IMACH(16) / 1023 /, SANITY/987/
C
C     MACHINE CONSTANTS FOR VAX.
C
C      DATA IMACH( 1) /    5 /
C      DATA IMACH( 2) /    6 /
C      DATA IMACH( 3) /    7 /
C      DATA IMACH( 4) /    6 /
C      DATA IMACH( 5) /   32 /
C      DATA IMACH( 6) /    4 /
C      DATA IMACH( 7) /    2 /
C      DATA IMACH( 8) /   31 /
C      DATA IMACH( 9) / 2147483647 /
C      DATA IMACH(10) /    2 /
C      DATA IMACH(11) /   24 /
C      DATA IMACH(12) / -127 /
C      DATA IMACH(13) /  127 /
C      DATA IMACH(14) /   56 /
C      DATA IMACH(15) / -127 /
C      DATA IMACH(16) /  127 /, SANITY/987/
C
C     MACHINE CONSTANTS FOR OPUS 300PM (BASED ON CLIPPER)
C
C      DATA IMACH( 1) / 5 /
C      DATA IMACH( 2) / 6 /
C      DATA IMACH( 3) / 8 /
C      DATA IMACH( 4) / 0 /
C      DATA IMACH( 5) / 32 /
C      DATA IMACH( 6) / 4 /
C      DATA IMACH( 7) / 2 /
C      DATA IMACH( 8) / 31 /
C      DATA IMACH( 9) / 2147483647 /
C      DATA IMACH(10) / 2 /
C      DATA IMACH(11) / 24 /
C      DATA IMACH(12) / -148 /
C      DATA IMACH(13) /  128 /
C      DATA IMACH(14) / 53 /
C      DATA IMACH(15) / -1024 /
C      DATA IMACH(16) /  1073 /, SANITY/987/
C
C  ***  ISSUE STOP 777 IF ALL DATA STATEMENTS ARE COMMENTED...
      IF (SANITY .NE. 987) STOP 777
      IF (I .LT. 1  .OR.  I .GT. 16) GO TO 10
C
      I1MACH = IMACH(I)
C/6S
C/7S
C     IF (I .EQ. 6) I1MACH = 1
C/
      RETURN
C
 10   WRITE(OUTPUT,9000)
 9000 FORMAT(39H1ERROR    1 IN I1MACH - I OUT OF BOUNDS)
C
      CALL FDUMP
C
      STOP
C
      END
      INTEGER FUNCTION I8SAVE(ISW,IVALUE,SET)
C
C  IF (ISW = 1) I8SAVE RETURNS THE CURRENT ERROR NUMBER AND
C               SETS IT TO IVALUE IF SET = .TRUE. .
C
C  IF (ISW = 2) I8SAVE RETURNS THE CURRENT RECOVERY SWITCH AND
C               SETS IT TO IVALUE IF SET = .TRUE. .
C
      LOGICAL SET
C
      INTEGER IPARAM(2)
      EQUIVALENCE (IPARAM(1),LERROR) , (IPARAM(2),LRECOV)
C
C  START EXECUTION ERROR FREE AND WITH RECOVERY TURNED OFF.
C
      DATA LERROR/0/ , LRECOV/2/
C
      I8SAVE=IPARAM(ISW)
      IF (SET) IPARAM(ISW)=IVALUE
C
      RETURN
C
      END
      SUBROUTINE S88FMT( N, W, IFMT )
C
C  S88FMT  REPLACES IFMT(1), ... , IFMT(N) WITH
C  THE CHARACTERS CORRESPONDING TO THE N LEAST SIGNIFICANT
C  DIGITS OF W.
C
      INTEGER N,W
C/6S
      INTEGER IFMT(N)
C/7S
C     CHARACTER*1 IFMT(N)
C/
C
      INTEGER NT,WT
C
C/6S
      INTEGER DIGITS(10)
      DATA DIGITS( 1) / 1H0 /
      DATA DIGITS( 2) / 1H1 /
      DATA DIGITS( 3) / 1H2 /
      DATA DIGITS( 4) / 1H3 /
      DATA DIGITS( 5) / 1H4 /
      DATA DIGITS( 6) / 1H5 /
      DATA DIGITS( 7) / 1H6 /
      DATA DIGITS( 8) / 1H7 /
      DATA DIGITS( 9) / 1H8 /
      DATA DIGITS(10) / 1H9 /
C/7S
C     CHARACTER*1 DIGITS(10)
C     DATA DIGITS( 1) / '0' /
C     DATA DIGITS( 2) / '1' /
C     DATA DIGITS( 3) / '2' /
C     DATA DIGITS( 4) / '3' /
C     DATA DIGITS( 5) / '4' /
C     DATA DIGITS( 6) / '5' /
C     DATA DIGITS( 7) / '6' /
C     DATA DIGITS( 8) / '7' /
C     DATA DIGITS( 9) / '8' /
C     DATA DIGITS(10) / '9' /
C/
C
      NT = N
      WT = W
C
 10   IF (NT .LE. 0) RETURN
        IDIGIT = MOD( WT, 10 )
        IFMT(NT) = DIGITS(IDIGIT+1)
        WT = WT/10
        NT = NT - 1
        GO TO 10
C
      END
      SUBROUTINE FDUMP
C  THIS IS A DUMMY ROUTINE TO BE SENT OUT ON
C  THE PORT SEDIT TAPE
C
      RETURN
      END
