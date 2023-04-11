      implicit integer (a-z)
 
c      parameter (maxn=625,fldmax=2*maxn)   
      parameter (maxn=1500,maxn2=3000)   
      parameter (fldmax=maxn2)   
      parameter (infinity=16000,sortmax=30000)
      parameter (mxbits=(maxn*(maxn+1)+31)/32)
      parameter (maxtloops=40)
      parameter (maxsiz=10000)
 
      integer*2 vst(maxn*maxn),wst(maxn*maxn)
      integer poppen(4),maxpen
      real prelog
      dimension newnum(maxsiz),hstnum(fldmax),force(fldmax),
     . numseq(fldmax), work(fldmax,0:2),
     . stack(5,5,5,5),tstk(5,5,5,5),dangle(5,5,5,2),hairpin(30)
      dimension bulge(30),inter(30),eparam(10),cntrl(10),nsave(2)
      common /main/ vst,wst,newnum,hstnum,force,numseq,work,stack,tstk,
     . dangle,hairpin,bulge,inter,eparam,cntrl,nsave,poppen,maxpen,prelog
 
      character*1 seq(maxsiz)
c      character*5 inbuf
      character*10 progtitle
      character*30 seqlab
      common /seq/ seq,seqlab
      data progtitle/'crna'/

      dimension list(100,4)
      common /list/ list,listsz
      common /nm/ n,vmin
 
      dimension basepr(maxn)
      common /traceback/ basepr
 
      dimension heapi(sortmax+1),heapj(sortmax+1)
      common /heap/ heapi,heapj,num

      integer*2 marks(mxbits),force2(mxbits)
      common/bits/marks,force2
 
      integer*2 tloop(maxtloops,2),numoftloops
      common/tloops/tloop,numoftloops
