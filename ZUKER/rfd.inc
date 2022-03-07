      implicit integer (a-z)
      parameter (maxn=1500,maxn2=3000)
      parameter (fldmax=maxn2)
      parameter (infinity=16000,sortmax=30000)
      parameter (mxbits=(maxn*(maxn+1)+31)/32)
      parameter (maxtloops=40)
      parameter (maxsiz=10000)
 
      integer*2 vst(maxn*maxn),wst1(maxn*maxn),wst2(maxn*maxn)
      integer poppen(4),maxpen
      real prelog
 
      dimension newnum(maxsiz),hstnum(maxn2),force(maxn2),numseq(maxn2),
     . work1(maxn2,0:2),work2(maxn2),
     . stack(5,5,5,5),tstk(5,5,5,5),dangle(5,5,5,2),hairpin(30)
      dimension bulge(30),inter(30),eparam(10),cntrl(10),nsave(2)
c      common /main/ newnum,hstnum,force,work1,work2,
      common /main/ newnum,hstnum,force,work1,work2,
     . stack,tstk,dangle,hairpin,bulge,inter,eparam,cntrl,nsave,n,
     . numseq,poppen,prelog,maxpen,vst,wst1,wst2
 
      character*1 seq(maxsiz)
c      character*5 inbuf
      character*10 progtitle
      character*30 seqlab
      common /seq/ seq,seqlab
      dimension list(100,4)
      common /list/ list,listsz
      common /nm/ vmin
      data progtitle/'lrna'/
 
      dimension basepr(maxn)
      common /traceback/ basepr
 
      dimension heapi(sortmax+1),heapj(sortmax+1)
      common /heap/ heapi,heapj,num

      integer*2 marks(mxbits),force2(mxbits)
      common/bits/marks,force2
 
      integer*2 tloop(maxtloops,2),numoftloops
      common/tloops/tloop,numoftloops


















