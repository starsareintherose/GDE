c     Energy funtion.
c     ERG is the energy of a loop closed by I,J (new numbering).
c     IP,JP is the other closing base-pair when MODE = 2 or 3.
c
c 1/2 Asym. loop correction
c Extrapolate loops with dG(n)=dG(30)+1.75*ln(n/30)
c Hairpins of 3 have no terminal stack.
c
      function erg(mode,i,j,ip,jp)
      include 'rfd.inc'
      dimension e(4)
      integer*2 tlink,tlptr
      logical fce
 
100   if (mode.eq.1) then
c     Read energy files
        call ergread
        erg = 0
        return
      endif
 
      erg = 0
c     Do not allow prohibited bases to pair.
      if (force(i).eq.1.or.force(j).eq.1) then
         erg = infinity
         return
      endif
 
      if (mode.lt.6) then
c     Add bonus energy to force base-pairs.
         if (force(i).eq.2.or.force(j).eq.2.or.fce(i,j)) then
            erg = erg + eparam(9)
            if (force(i).eq.2.and.force(j).eq.2) erg = erg + eparam(9)
         endif
      endif
 
      goto (100,200,300,400,500,600,700),mode
 
c     Nucleotide accesssibility option.
200   if (force(i).eq.3.or.force(jp).eq.3) then
         erg = infinity
         return
      endif
c     Stacking energy.
      erg = erg + stack(numseq(i),numseq(j),numseq(ip),numseq(jp))
     .          + eparam(1)
      return
 
300   size1 = ip - i - 1
      size2 = j - jp - 1
      if (size1.eq.0.or.size2.eq.0) then
c     Check for nucleotide accessibility.
          if (size1.eq.0.and.force(i).eq.3) then
             erg = infinity
             return
          endif
          if (size2.eq.0.and.force(jp).eq.3) then
             erg = infinity
             return
          endif
          size = size1+size2
c     Bulge loop.
          if (size.eq.1) then
            erg = erg + stack(numseq(i),numseq(j),numseq(ip),numseq(jp))
     .              + bulge(1) + eparam(2)
          else if (size.gt.30) then
              loginc = int(prelog*log((float(size)/30.0)))
              erg = erg + bulge(30) + loginc + eparam(2)
          else
              erg = erg + bulge(size) + eparam(2)
          endif
          return
      else
          size = size1+size2
          lopsid = abs((size1-size2))
c     Interior loop.
          if (size.gt.30) then
              loginc = int(prelog*log((float(size)/30.0)))
           erg = erg + tstk(numseq(i),numseq(j),numseq(i+1),numseq(j-1))
     .        + tstk(numseq(jp),numseq(ip),numseq(jp+1),numseq(ip-1))
     .        + inter(30) + loginc + eparam(3)
     .        + min0(maxpen,(lopsid*poppen(min0(4,size1,size2))))
          else
           erg = erg + tstk(numseq(i),numseq(j),numseq(i+1),numseq(j-1))
     .        + tstk(numseq(jp),numseq(ip),numseq(jp+1),numseq(ip-1))
     .        + inter(size) + eparam(3)
     .        + min0(maxpen,(lopsid*poppen(min0(4,size1,size2))))
          endif
          return
      endif
 
400   size = j-i-1
c     Hairpin loop.
      if ((size.eq.3).and.fce(i,j).and.seq(hstnum(i+1)).eq.' ') then
c        Closed excision
         erg = eparam(9)
         return
      endif
      if (size.gt.30) then
          loginc = int(prelog*log((float(size)/30.0)))
          erg = erg + tstk(numseq(i),numseq(j),numseq(i+1),numseq(j-1))
     .          + hairpin(30) + loginc + eparam(4)
      else if (size.lt.4) then
          erg=erg+hairpin(size)+eparam(4)
      else
c
          tlink=0
          if (size.eq.4) then
           key=((numseq(i+4)*8+numseq(i+3))*8+numseq(i+2))*8+numseq(i+1)
             tlptr=1
           do while ((tlptr.le.numoftloops).and.(tloop(tlptr,1).ne.key))
                 tlptr=tlptr+1
             enddo
             if (tlptr.le.numoftloops) tlink=tloop(tlptr,2)
          endif
          erg = erg + tstk(numseq(i),numseq(j),numseq(i+1),numseq(j-1))
     .          + hairpin(size) + eparam(4) + tlink
      endif
      return
 
c     Multi-branch loop.
500   do 501 ii = 1,4
501     e(ii) = infinity
 
      if (i.le.n-2) then
 
        ind1 = (n-1)*i
        ind2 = (n-1)*(i+1)
 
      else if (i.eq.n-1) then
 
        ind1 = (n-1)*i
        ind2 = -n
 
      else
 
        ind1 = -n
        ind2 = -1
 
      endif
 
 
      do k = i+2,j-3
c     EPARAM(6) is the energy penalty for each single-stranded base
c     in a multi-loop. EPARAM(10) is the energy penalty for each base-pair
c     closing a multi-loop.
c     No dangling ends next to the I,J base-pair.
         e(1) = min0(e(1),wst(ind1+k)+work(k+1,mod(j-1,3)))
c     I+1 dangles on the I,J base-pair.
         e(2) = min0(e(2),dangle(numseq(i),numseq(j),numseq(i+1),1)
     .                 + wst(ind2+k) + work(k+1,mod(j-1,3)) + eparam(6))
c     J-1 dangles on the I,J base-pair.
         e(3) = min0(e(3),dangle(numseq(i),numseq(j),numseq(j-1),2)
     .                 + wst(ind1+k) + work(k+1,mod(j-2,3)) + eparam(6))
c     Both I+1 and J-1 dangle on the I,J base-pair.
         e(4) = min0(e(4),dangle(numseq(i),numseq(j),numseq(i+1),1)
     .                  + dangle(numseq(i),numseq(j),numseq(j-1),2)
     .               + wst(ind2+k) + work(k+1,mod(j-2,3)) + 2*eparam(6))
      enddo
c     EPARAM(5) is the energy penalty for closing a multi-loop.
      erg = erg + eparam(5) + eparam(10) + min0(e(1),e(2),e(3),e(4))
      return
 
c     Dangling base stacking energy. IP dangles over the I,J
c     base-pair. 3' or 5' dangle if JP = 1 or 2 respectively.
600   erg = erg + dangle(numseq(i),numseq(j),numseq(ip),jp)
      return
 
700   if (force(i).eq.3.or.force(jp).eq.3) then
         erg = infinity
         return
      endif
c     Terminal stack or mismatch energy.
      erg = erg + tstk(numseq(i),numseq(j),numseq(ip),numseq(jp))
      return
      end
 
 
 
 
      subroutine fill
c     This subroutine computes the arrays of optimal energies.
      include 'rfd.inc'
      dimension inc(5,5)
      data loop/3/,inc/0,0,0,1,0,0,0,1,0,0,0,1,0,1,0,1,0,1,0,0,0,0,0,0,0
     ./
 
      vmin = infinity
      if (n.le.80) then
         pinc = 5
      elseif (n.le.100) then
         pinc = 2
      else
         pinc = 1
      endif
      pcnt = pinc
      crit = n*n*n/50
 
      do j = 1,2*n-1
c     How far along is the computation?
         if (n.gt.10) then
           if (j.le.n) then
              if (j**3.ge.pcnt*crit) then
                 write (6,1000) pcnt
                 pcnt = pcnt + pinc
              endif
           else
              if ((2*n-j)**3.le.(100-pcnt)*crit) then
                 write (6,1000) pcnt
                 pcnt = pcnt + pinc
              endif
           endif
         endif
1000     format ('+',5x,i4,'%')
 
         do i = min0(j,n),max0(1,j-n+1),-1
           vij =  infinity
           wij = infinity
           if (j-i.le.loop) goto 300
c          Test for a prohibited base-pair or a pair which cannot form
c          a base-pair.
           if (vst((n-1)*(i-1)+j).eq.1.or.inc(numseq(i),numseq(j)).eq.0)
     .       goto 200
c          Compute VIJ, the minimum energy of the fragment from I to J
c          inclusive where I and J base-pair with one another.
c          Perhaps I,J closes a hairpin loop.
           vij = min0(vij,erg(4,i,j,i,j))
           if (j-i-1.ge.loop+2) then
c          Perhaps I,J stacks over I+1,J-1.
              vij = min0(vij,erg(2,i,j,i+1,j-1)+v(i+1,j-1))
           endif
c          Search for a bulge or interior loop.
           if (j-i-1.ge.loop+3) then
              do d = j-i-3,1,-1
                 do ip = i+1,j-1-d
                    jp = d+ip
                    if (j-i-2-d.gt.eparam(7)) goto 100
                    if (abs(ip-i+jp-j).le.eparam(8)) then
                       if (ip.gt.n) then
                          vij = min0(vij,erg(3,i,j,ip,jp)+vst((n-1)*
     .                       (ip-n-1)+jp-n))
                       else
                          vij = min0(vij,erg(3,i,j,ip,jp)+vst((n-1)*
     .                       (ip-1)+jp))
                       endif
                    endif
                 enddo
              enddo
           endif
 
100        if (j-i-1.ge.2*loop+4) then
c             Perhaps I,J closes a multi-loop.
              vij=min0(vij,erg(5,i,j,i,j))
           endif
 
 
c          Compute WIJ, the minimum energy of a non-empty folding on I to
c          J inclusive. This is the circular folding program and so there
c          are no exterior bases.
200     wij = min0 ( wij, eparam(10)+vij, v(i+1,j)+eparam(6)+eparam(10)+
     .                 erg(6,j,i+1,i,2), v(i,j-1)+eparam(6)+eparam(10)+
     .             erg(6,j-1,i,j,1), v(i+1,j-1)+2*eparam(6)+eparam(10)+
     .                 erg(6,j-1,i+1,i,2)+erg(6,j-1,i+1,j,1),
     .                 w(i+1,j)+eparam(6), w(i,j-1)+eparam(6) )
 
           if (j-i-1.gt.2*loop+2) then
              index = (n-1)*(i-1)
c             Check for open bifurcation.
              do k = i,j-1
                wij = min0(wij,wst(index+k)+work(k+1,mod(j,3)))
              enddo
           endif
 
c          Store VIJ and WIJ. They can be regarded as elements V(I,J)
c          and W(I,J) in a two dimensional array. They are actually
c          stored in the one dimensional arrays VST and WST in position
c          (N-1)*(I-1) + J.
c          Columns J,J-1 and J-2 of W are stored again in the WORK array.
c          This is done to reduce virtual memory swaps.
300        vst((n-1)*(i-1)+j) = vij
           wst((n-1)*(i-1)+j) = wij
           work(i,mod(j,3)) = wij
           if (j.gt.n) then
c          VMIN is the minimum folding energy of the entire sequence.
c              vmin = min0(vmin,vst((n-1)*(i-1)+j)+vst((n-1)*(j-n-1)+i))
              vmin = min(vmin,vst((n-1)*(i-1)+j)+vst((n-1)*(j-n-1)+i))
           endif
         enddo
         if (j.ge.n) then
           do k = j+1,n+1,-1
c             Fill in some WORK array values before beginning work on the
c             next column.
              work(k,mod(j+1,3)) = wst((k-n-1)*(n-1)+j+1-n)
           enddo
         endif
      enddo
      return
      end
 
c     Used to recall values of V which are actually stored in VST.
      function v(i,j)
      include 'rfd.inc'
 
      if (i.gt.n) then
        v = vst((n-1)*(i-n-1)+j-n)
      else
        v = vst((n-1)*(i-1)+j)
      endif
      return
      end
 
 
c     Used to recall values of W which are actually stored in WST.
      function w(i,j)
      include 'rfd.inc'
 
      if (i.gt.n) then
        w = wst((n-1)*(i-n-1)+j-n)
      else
        w = wst((n-1)*(i-1)+j)
      endif
      return
      end
 
 
c     Computes an optimal structure on the subsequence II to JI where
c     II and JI must base-pair with one another. ERROR = 0 is normal
c     termination.
c     NFORCE is the number of forced base-pairs encountered in the traceback.
      subroutine trace(ii,ji,nforce,error)
      include 'rfd.inc'
      logical fce

      error = 0
 
c     Zero the appropriate region of BASEPR.
      if (ji.le.n) then
        do k=ii,ji
          basepr(k) = 0
        enddo
      else
        do k=1,ji-n
          basepr(k) = 0
        enddo
        do k = ii,n
          basepr(k) = 0
        enddo
      endif
c     Initialize the stack of outstanding base-pairs and push
c     II, JI and V(II,JI) on to the stack. The fourth stack position
c     is unused in this subroutine.
      call initst
      call push(ii,ji,v(ii,ji),0)
      nforce = 0
 
c     Pull a fragment and its expected energy from the stack.
c     End if there are no fragments left.
100   stz = pull(i,j,e,xx)
      if (stz.ne.0) return
 
c     Do I and J base-pair with one another?
      if (e.eq.v(i,j)) goto 300
 
      tst = w(i+1,j) + eparam(6)
      do while (e.eq.tst)
c       Whittle away from the 5' end.
        i = i + 1
        if (i.ge.j) goto 100
        e = w(i,j)
        tst = w(i+1,j) + eparam(6)
      enddo
 
      tst = w(i,j-1) + eparam(6)
      do while (e.eq.tst)
c       Whittle away from the 3' end.
        j = j - 1
        if (i.ge.j) goto 100
        e = w(i,j)
        tst = w(i,j-1) + eparam(6)
      enddo
 
      tst1 = v(i+1,j) + eparam(6) + eparam(10) +
     .        dangle(numseq(j),numseq(i+1),numseq(i),2)
      tst2 = v(i,j-1) + eparam(6) + eparam(10) +
     .        dangle(numseq(j-1),numseq(i),numseq(j),1)
      tst3 = v(i+1,j-1) + 2*eparam(6) + eparam(10) +
     .         dangle(numseq(j-1),numseq(i+1),numseq(i),2)
     .       + dangle(numseq(j-1),numseq(i+1),numseq(j),1)
      if (e.eq.tst1) then
c        I dangles over I+1,J.
         i = i + 1
         e = v(i,j)
      else if (e.eq.tst2) then
c        J dangles over I,J-1.
         j = j - 1
         e = v(i,j)
      else if (e.eq.tst3) then
c        Both I and J dangle over I+1,J-1.
         i = i + 1
         j = j - 1
         e = v(i,j)
      endif
c     Check for stem closing a multi-loop.
      if (e.eq.v(i,j)+eparam(10)) e = v(i,j)
 
      if (e.ne.v(i,j)) then
c        Cannot chop away at ends any more and still the ends do not
c        base-pair with one another. Structure MUST bifurcate (OPEN).
         k = i+1
200      if (k.eq.j) then
c          Structure will not split. Error
           ii = hstnum(i)
           ji = hstnum(j)
           error = 10
           return
         endif
         if (e.eq.w(i,k) + w(k+1,j)) then
c           Best structure on I,J splits into best structures on I,K and
c           K+1,J. Push these fragments on to the stack.
            call push(i,k,w(i,k),0)
            call push(k+1,j,w(k+1,j),0)
            goto 100
         else
            k = k + 1
            goto 200
         endif
      endif
 
c     Base-pair found. Base-pairs are stored in the range 1 <= I < J <= N.
300   if (j.le.n) then
         basepr(i) = j
         basepr(j) = i
      else if (i.gt.n) then
         basepr(i-n) = j-n
         basepr(j-n) = i-n
         i = i - n
         j = j - n
      else
         basepr(j-n) = i
         basepr(i) = j-n
      endif
c     Check if this is a forced base-pair.
      if (force(i).eq.2.or.force(j).eq.2.or.fce(i,j)) 
     .  nforce = nforce + 1
      if (force(i).eq.2.and.force(j).eq.2) nforce = nforce + 1
c     Perhaps I,J stacks over I+1,J-1?
      if (e.eq.erg(2,i,j,i+1,j-1)
     . + v(i+1,j-1)) then
         i = i + 1
         j = j - 1
         e = v(i,j)
         goto 300
      endif
c     Perhaps I,J closes a hairpin loop?
      if (e.eq.erg(4,i,j,i,j)) goto 100
 
c     Define E' ( EP in the program ) to be E corrected by a
c     possible bonus energy for forced base-pairing.
c
      ep = e
      if (force(i).eq.2.or.force(j).eq.2.or.fce(i,j)) 
     .  ep = ep - eparam(9)
      if (force(i).eq.2.and.force(j).eq.2) ep = ep - eparam(9)
 
      k = i+2
c     Perhaps I,J closes a multi-loop?
400   if (k.ge.j-3) goto 500
      if (ep.eq.w(i+1,k) + w(k+1,j-1) + eparam(10) + eparam(5)) then
c        Multi-loop. No dangling ends on I,J.
         call push(i+1,k,w(i+1,k),0)
         call push(k+1,j-1,w(k+1,j-1),0)
         goto 100
      else if (ep.eq.erg(6,i,j,i+1,1)+w(i+2,k)+w(k+1,j-1)+eparam(10)+
     .                     eparam(6)+eparam(5)) then
c        Multi-loop. I+1 dangles over I,J base-pair.
         call push(i+2,k,w(i+2,k),0)
         call push(k+1,j-1,w(k+1,j-1),0)
         goto 100
      else if (ep.eq.erg(6,i,j,j-1,2)+w(i+1,k)+w(k+1,j-2)+eparam(10)+
     .                     eparam(6)+eparam(5)) then
c        Multi-loop. J-1 dangles over I,J base-pair.
         call push(i+1,k,w(i+1,k),0)
         call push(k+1,j-2,w(k+1,j-2),0)
         goto 100
      else if (ep.eq.erg(6,i,j,i+1,1)+erg(6,i,j,j-1,2)+w(i+2,k)
     .     +w(k+1,j-2)+ eparam(10)+2*eparam(6)+eparam(5)) then
c        Multi-loop. Both I+1 and J-1 dangle over the I,J base-pair.
         call push(i+2,k,w(i+2,k),0)
         call push(k+1,j-2,w(k+1,j-2),0)
         goto 100
      else
         k = k + 1
         goto 400
      endif
 
 
c     None of the above work. I,J MUST close a bulge or interior loop.
500   do d = j-i-3,1,-1
        do ip = i+1,j-1-d
           jp = d+ip
           if (j-i-2-d.gt.eparam(7)) then
c             Error, bulge or interior loop not found.
              ii = hstnum(i)
              ji = hstnum(j)
              error = 11
              return
           endif
           if (abs(ip-i+jp-j).le.eparam(8)) then
              if (e.eq.erg(3,i,j,ip,jp)+v(ip,jp)) then
                 i = ip
                 j = jp
                 e = v(i,j)
                 goto 300
              endif
           endif
        enddo
      enddo
c     Error, bulge or interior loop not found.
      ii = hstnum(i)
      ji = hstnum(j)
      error = 11
      return
      end
 
 
c     Store results of a SAVE run for a CONTINUATION run.
      subroutine putcont
      include 'rfd.inc'
 
      write(30) n,nsave,vmin,listsz,seqlab
      write(30) stack,tstk,dangle,hairpin,bulge,inter,eparam
      write(30) (vst(i),i=1,n*n)
      write(30) (wst(i),i=1,n*n)
      write(30) (seq(i),i=nsave(1),nsave(2))
      write(30) ((list(i,j),i=1,100),j=1,4)
      write(30) tloop,numoftloops
      write(30) (poppen(i),i=1,4),maxpen,prelog
      return
      end
 
c     Read results of a SAVE run for a CONTINUATION run.
      subroutine getcont
      include 'rfd.inc'
 
      read(30,end=10) n,nsave,vmin,listsz,seqlab
      read(30,end=10) stack,tstk,dangle,hairpin,bulge,inter,eparam
      read(30,end=10) (vst(i),i=1,n*n)
      read(30,end=10) (wst(i),i=1,n*n)
      read(30,end=10) (seq(i),i=nsave(1),nsave(2))
      read(30,end=10) ((list(i,j),i=1,100),j=1,4)
      read(30,err=10) tloop,numoftloops
      read(30,err=10) (poppen(i),i=1,4),maxpen,prelog
      goto 11
 
10    call errmsg(40,0,0)
 
11    return
      end
 
