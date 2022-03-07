c     Energy funtion.
c     ERG is the energy of a loop closed by I,J (new numbering).
c     IP,JP is the other closing base-pair when MODE = 2 or 3.
c     The ends of the sequence cannot be contained in a hairpin, bulge
c     or interior loop. By convention, the ends of the sequence are
c     put into a special kind of multi-loop. This can be called an
c     exterior loop or an open multi-loop.
c
c 1/2 Asym. loop correction
c Extrapolate loops with dG(n)=dG(30)+1.75*ln(n/30)
c
c* Hairpins of 3 have no terminal stack.
c
      function erg(mode,i,j,ip,jp)
      include 'rfd.inc'
      dimension e(4)
      integer*2 tlink,tlptr
      logical fce
c
 
100   if (mode.eq.1) then
c       Read energy files.
        call ergread
        erg = 0
        return
      endif
 
      erg = 0
c     Do not allow prohibited base to pair.
      if (force(i).eq.1.or.force(j).eq.1) then
         erg = infinity
         return
      endif

      if (mode.lt.6) then
c        Add bonus energy to force base-pairs.
         if ((force(i).eq.2).or.(force(j).eq.2).or.fce(i,j)) 
     .    then
            erg = erg + eparam(9)
            if ((force(i).eq.2).and.(force(j).eq.2)) 
     .          erg = erg + eparam(9)
         endif


      endif

      goto (100,200,300,400,500,600,700),mode
 
c     Nucleotide accessibility option.
200   if (force(i).eq.3.or.force(jp).eq.3) then
         erg = infinity
         return
      endif
c     Molecule is not circular. N is not covalently bonded to N+1.
      if (i.eq.n.or.j.eq.n+1) then
         erg = infinity
         return
      endif
c     Stacking energy.
      erg = erg + stack(numseq(i),numseq(j),numseq(ip),numseq(jp))
     .          + eparam(1)
      return
 
300   if ((i.le.n.and.ip.gt.n).or.(jp.le.n.and.j.gt.n)) then
c          Loop is not allowed to contain the ends of the sequence.
           erg = infinity
           return
      endif
c

      size1 = ip - i - 1
      size2 = j - jp - 1
      if (size1.eq.0.or.size2.eq.0) then
c         Check for nucleotide accessibility.
          if (size1.eq.0.and.force(i).eq.3) then
             erg = infinity
             return
          endif
          if (size2.eq.0.and.force(jp).eq.3) then
             erg = infinity
             return
          endif
          size = size1+size2
c         Bulge loop energy.
          if (size.eq.1) then
            erg = erg + stack(numseq(i),numseq(j),numseq(ip),numseq(jp))
     .              + bulge(size) + eparam(2)
          elseif (size.gt.30) then
             loginc = int(prelog*log((float(size)/30.0)))
             erg = erg + bulge(30) + loginc + eparam(2)
          else
             erg = erg + bulge(size) + eparam(2)
          endif
          return
      else
          size = size1+size2
          lopsid = abs((size1-size2))
c         Interior loop.
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
 
400   if (i.le.n.and.j.gt.n) then
c         Hairpin loop must not contain the ends of the sequence.
          erg = infinity
          return
      endif
c
      size = j-i-1
      if ((size.eq.3).and.fce(i,j).and.seq(hstnum(i+1)).eq.' ') then
c        Closed excision
         erg = eparam(9)
         return
      endif
      if (size.gt.30) then
         loginc = int(prelog*log((float(size)/30.0)))
         erg = erg + tstk(numseq(i),numseq(j),numseq(i+1),numseq(j-1))
     .          + hairpin(30) + loginc + eparam(4)
      else if (size .lt. 4) then
c
c* Special case for hairpin of 3
c
         erg = erg + hairpin(size) + eparam(4)
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
 
c     Multi-branch  (or multi-) loop closed by I,J.
500   do 501 ii = 1,4
501     e(ii) = infinity

      if (i+2.gt.j-3) then
c        There are at most 3 bases between I and J. The fragment from
c        I to J inclusive contains the origen.
         e(1) = 0
         if (i.ne.n) e(2) = dangle(numseq(i),numseq(j),numseq(i+1),1)
         if (j.ne.n+1) e(3) = dangle(numseq(i),numseq(j),numseq(j-1),2)
         if (i.ne.n.and.j.ne.n+1) then
           e(4) = dangle(numseq(i),numseq(j),numseq(i+1),1) +
     .            dangle(numseq(i),numseq(j),numseq(j-1),2)
         endif
      else if (i.ge.n-1) then
c        I is at or next to the end of the sequence.
         e(1) = w2(n+1,j-1)
         if (i.ne.n) then
          e(2) = dangle(numseq(i),numseq(j),numseq(i+1),1) + w2(n+1,j-1)
           e(4) = dangle(numseq(i),numseq(j),numseq(i+1),1) +
     .           dangle(numseq(i),numseq(j),numseq(j-1),2) + w2(n+1,j-2)
         endif
         e(3) = dangle(numseq(i),numseq(j),numseq(j-1),2) + w2(n+1,j-2)
      else if (j.eq.n+1.or.j.eq.n+2) then
c        J is at or next to the end of the sequence.
         e(1) = wst2((n-1)*i+n)
         e(2) = dangle(numseq(i),numseq(j),numseq(i+1),1) + 
     .     wst2((n-1)*(i+1)+n)
         if (j.ne.n+1) then
        e(3) = dangle(numseq(i),numseq(j),numseq(j-1),2)+wst2((n-1)*i+n)
           e(4) = dangle(numseq(i),numseq(j),numseq(i+1),1) +
     .   dangle(numseq(i),numseq(j),numseq(j-1),2) + wst2((n-1)*(i+1)+n)
         endif
       else
         ind1 = (n-1)*i
         ind2 = ind1 + n - 1
         do k = i+2,j-3
           if (k.eq.n) then
c            When K = N, the structure splits into two disconnected
c            pieces. This open multi-loop ( exterior loop ) is not given the
c            usual EPARAM(5),EPARAM(6) and EPARAM(10) destabilizing energies.
c
             ind3 = -n
c            No dangling ends next to the I,J base-pair.
c             e(1) = min0(e(1),wst2(ind1+k)+wst2(ind3+j-1))
             e(1) = min(e(1),wst2(ind1+k)+wst2(ind3+j-1))
c            I+1 dangles on the I,J base-pair.
             e(2) = min0(e(2),dangle(numseq(i),numseq(j),numseq(i+1),1)
     .                    + wst2(ind2+k) + wst2(ind3+j-1))
c            J-1 dangles on the I,J base-pair.
             e(3) = min0(e(3),dangle(numseq(i),numseq(j),numseq(j-1),2)
     .                    + wst2(ind1+k) + wst2(ind3+j-2))
c            Both I+1 and J-1 dangle on the I,J base-pair.
             e(4) = min0(e(4),dangle(numseq(i),numseq(j),numseq(i+1),1)
     .                  + dangle(numseq(i),numseq(j),numseq(j-1),2)
     .                    + wst2(ind2+k) + wst2(ind3+j-2))
           else
c            When K is not N, the ends of the sequence are not in the
c            loop. This is a proper multi-loop with an energy of EPARAM(6)
c            for each single-stranded base, an energy of EPARAM(10) for
c            each closing base-pair, plus and extra energy of EPARAM(5).
c            No dangling ends next to the I,J base-pair.
           e(1) = min0(e(1),wst1(ind1+k)+work1(k+1,mod(j-1,3))+eparam(5)
     .                    +eparam(10))
c            I+1 dangles on the I,J base-pair.
             e(2) = min0(e(2),dangle(numseq(i),numseq(j),numseq(i+1),1)
     .                    + wst1(ind2+k) + work1(k+1,mod(j-1,3)) + eparam(5)
     .                    + eparam(6) + eparam(10) )
c            J-1 dangles on the I,J base-pair.
             e(3) = min0(e(3),dangle(numseq(i),numseq(j),numseq(j-1),2)
     .                    + wst1(ind1+k) + work1(k+1,mod(j-2,3)) + eparam(5)
     .                    + eparam(6) + eparam(10) )
c            Both I+1 and J-1 dangle on the I,J base-pair.
             e(4) = min0(e(4),dangle(numseq(i),numseq(j),numseq(i+1),1)
     .                  + dangle(numseq(i),numseq(j),numseq(j-1),2)
     .                    + wst1(ind2+k) + work1(k+1,mod(j-2,3)) + eparam(5)
     .                      + 2*eparam(6) + eparam(10) )
            endif
         enddo
      endif
 
      erg = erg + min0(e(1),e(2),e(3),e(4))
      return
 
c     Dangling base stacking energy. IP dangles over the I,J base-pair.
c     3' or 5' dangle if JP = 1 or 2 respectively.
600   erg = erg + dangle(numseq(i),numseq(j),numseq(ip),jp)
      return

700   if (force(i).eq.3.or.force(jp).eq.3) then
         erg = infinity
         return
      endif
c     Terminal stack or mismatch energy.
      erg = erg + tstk(numseq(i),numseq(j),numseq(ip),numseq(jp))
      return
      end function erg





      subroutine fill
c     This subroutine computes the arrays of optimal energies.
      include 'rfd.inc'
      dimension inc(5,5),e1(5),e2(5)
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
c        How far along is the computation?
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
           w1ij = infinity
           w2ij = 0
           if (j.le.n) then
              if (j-i.le.loop) goto 300
           else
              if (i.eq.n.or.j.eq.n+1) goto 100
           endif
c          Test for a prohibited base-pair or a pair which cannot form.
           if (vst((n-1)*(i-1)+j).eq.1.or.inc(numseq(i),numseq(j)).eq.0) 
     .        goto 200
c          Compute VIJ, the minimum energy of the fragment from I to J
c          inclusive where I and J base-pair with one another.
c          Perhaps I,J closes a hairpin loop.
           vij = min0(vij,erg(4,i,j,i,j))
           if (j-i-1.ge.loop+2.or.j.gt.n) then
c             Perhaps I,J stacks over I+1,J-1.
              vij = min0(vij,erg(2,i,j,i+1,j-1)+v(i+1,j-1))
           endif
c          Search for the best bulge or interior loop closed by I,J.
           if (j-i-1.ge.loop+3.or.j.gt.n) then
              do d = j-i-3,1,-1
                 do ip = i+1,j-1-d
                    jp = d+ip
                    if (j-i-2-d.gt.eparam(7)) goto 100
                    if (abs(ip-i+jp-j).le.eparam(8)) then
                       if (ip.gt.n) then
                          vij = min0(vij,erg(3,i,j,ip,jp)+vst((n-1)*
     .                         (ip-n-1)+jp-n))
                       else
                          vij = min0(vij,erg(3,i,j,ip,jp)+vst((n-1)*
     .                         (ip-1)+jp))
                       endif
                    endif
                 enddo
              enddo
           endif
 
100        if (j-i-1.ge.2*loop+4.or.j.gt.n) then
c             Search for the best multi-loop closed by I,J.
              vij=min0(vij,erg(5,i,j,i,j))
           endif
c          Compute W1IJ and W2IJ.
c          A multi-loop containing N and 1 (ie. N+1) as single-stranded
c          bases is called an exterior loop. W1IJ is the minimum folding
c          energy of a non-empty folding on I to J inclusive where an
c          exterior loop is given an energy of EPARAM(5) plus EPARAM(6)
c          per single-stranded exterior base plus EPARAM(10) per
c          double-stranded exterior base-pair in addition to
c          possible dangling base energies. W2IJ is similarly
c          defined except that the folding can be empty and that an
c          exterior loop is given no energy other than possible dangling
c          base energies.
200        do ii = 1,5
             e1(ii) = infinity
             e2(ii) = infinity
           enddo
           if (i.ne.n) then
c            Add single-stranded I to an optimal structure containing
c            the base-pair I,J.
            e1(1) = v(i+1,j) + eparam(10) + eparam(6) + erg(6,j,i+1,i,2)
             e1(4) = w1(i+1,j) + eparam(6)
             e2(1) = v(i+1,j) + erg(6,j,i+1,i,2)
             e2(4) = w2(i+1,j)
           endif
           if (j.ne.n+1) then
c            Add single-stranded J to an optimal structure containing
c            the base-pair I,J.
            e1(2) = v(i,j-1) + eparam(10) + eparam(6) + erg(6,j-1,i,j,1)
             e1(5) = w1(i,j-1) + eparam(6)
             e2(2) = v(i,j-1) + erg(6,j-1,i,j,1)
             e2(5) = w2(i,j-1)
           endif
           if (i.ne.n.and.j.ne.n+1) then
c            Add single-stranded I and J to an optimal structure containing
c            the base-pair I+1,J-1.
             e1(3) = v(i+1,j-1) + eparam(10) + 2*eparam(6) +
     .                 erg(6,j-1,i+1,i,2) + erg(6,j-1,i+1,j,1)
             e2(3) = v(i+1,j-1) + erg(6,j-1,i+1,i,2)
     .              + erg(6,j-1,i+1,j,1)
           endif
 
           w1ij = min0(eparam(10)+vij,e1(1),e1(2),e1(3),e1(4),e1(5))
           w2ij = min0(vij,w2ij,e2(1),e2(2),e2(3),e2(4),e2(5))
 
           if (j-i-1.gt.2*loop+2.or.j.gt.n) then
              index = (n-1)*(i-1)
c             Search for an open bifurcation.
              do k = i,j-1
                if (k.eq.n) then
                   w1ij = min0(w1ij,wst2(index+k)+work2(k+1))
                else
                   w1ij = min0(w1ij,wst1(index+k)+work1(k+1,mod(j,3)))
                   w2ij = min0(w2ij,wst2(index+k)+work2(k+1))
                endif
              enddo
           endif
 
c          Store VIJ, W1IJ and W2IJ. They can be regarded as elements
c          V(I,J), W1(I,J) and W2(I,J) in two dimensional arrays. They
c          are actually stored in one dimensional arrays VST, WST1 and
c          WST2 is position (N-1)*(I-1) + J.
c          Columns J,J-1 and J-2 of W1 are also stored in the work array,
c          WORK1. Column J of W2 is stored again in the work array WORK2.
c          This is done to reduce virtual memory swaps.
300        vst((n-1)*(i-1)+j) = vij
           wst1((n-1)*(i-1)+j) = w1ij
           wst2((n-1)*(i-1)+j) = w2ij
           work1(i,mod(j,3)) = w1ij
           work2(i) = w2ij
           if (j.gt.n) then
c             VMIN is the minimum folding energy of the entire sequence.
c              vmin = min0(vmin,vst((n-1)*(i-1)+j)+vst((n-1)*(j-n-1)+i))
              vmin = min(vmin,vst((n-1)*(i-1)+j)+vst((n-1)*(j-n-1)+i))
           endif
         enddo
         if (j.ge.n) then
            do k = j+1,n+1,-1
c              Fill in some work array values before beginning on
c              the next column.
               work1(k,mod(j+1,3)) = wst1((k-n-1)*(n-1)+j+1-n)
               work2(k) = wst2((k-n-1)*(n-1)+j+1-n)
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
 
c     Used to recall values of W1 which are actually stored in WST1.
      function w1(i,j)
      include 'rfd.inc'
 
      if (i.gt.n) then
        w1 = wst1((n-1)*(i-n-1)+j-n)
      else
        w1 = wst1((n-1)*(i-1)+j)
      endif
      return
      end
 
c     Used to recall values of W2 which are actually stored in WST2.
      function w2(i,j)
      include 'rfd.inc'
 
      if (i.gt.n) then
        w2 = wst2((n-1)*(i-n-1)+j-n)
      else
        w2 = wst2((n-1)*(i-1)+j)
      endif
      return
      end
 
 
 
 
 
c     Computes an optimal structure on the subsequence from II to JI
c     where II and JI must base-pair with each other. ERROR = 0
c     indicates a normal termination.
c     NFORCE is the number of forced base-pairs encountered during the
c     traceback.
c     Base-pair information is stored in the array BASEPR.
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
c     Initialize the stack of outstanding base-pairs and push II,JI,
c     V(II,JI) and 0 on to the stack.
      call initst
      call push(ii,ji,v(ii,ji),0)
      nforce = 0
 
100   i = j
      do while (i.eq.j)
c       Pull a fragment ( I to J ) and its expected energy ( E ) from
c       the stack. OPENL = 1 indicates that the free bases are part of
c       an exterior loop. OPENL = 0 (ie. closed) indicates that the
c       free bases are part of a multi-loop.
        stz = pull(i,j,e,openl)
        if (stz.ne.0) return
      enddo
c     Do I and J base-pair with one another?
      if (e.eq.v(i,j)) goto 300
 
      if (openl.eq.0) then
 
             do while (e.eq.w1(i+1,j)+eparam(6))
c            Whittle away from the 5' end.
             i = i + 1
             e = w1(i,j)
             if (i.ge.j) goto 100
         enddo
         do while (e.eq.w1(i,j-1)+eparam(6))
c            Whittle away from the 3' end.
             j = j - 1
             e = w1(i,j)
             if (i.ge.j) goto 100
         enddo
 
         if (e.eq.v(i+1,j)+eparam(10)+eparam(6)+erg(6,j,i+1,i,2)) then
c            I dangles over I+1,J.
             i = i + 1
             e = v(i,j)
       elseif (e.eq.v(i,j-1)+eparam(10)+eparam(6)+erg(6,j-1,i,j,1)) then
c            J dangles over I,J-1.
             j = j - 1
             e = v(i,j)
         elseif (e.eq.v(i+1,j-1) + eparam(10) + 2*eparam(6) +
     .                erg(6,j-1,i+1,i,2) + erg(6,j-1,i+1,j,1) ) then
c            Both I and J dangle over I+1,J-1.
             i = i + 1
             j = j - 1
             e = v(i,j)
         endif
c        Check for stem closing a multi-loop.
         if (e.eq.v(i,j)+eparam(10)) e = v(i,j)
 
      else
 
         do while (e.eq.w2(i+1,j))
c            Whittle away at the 5' end.
             i = i + 1
             if (i.ge.j) goto 100
         enddo
         do while (e.eq.w2(i,j-1))
c            Whittle away at the 3' end.
             j = j - 1
             if (i.ge.j) goto 100
         enddo
 
         if (e.eq.v(i+1,j)+erg(6,j,i+1,i,2)) then
c            I dangles over I+1,J.
             i = i + 1
             e = v(i,j)
         elseif (e.eq.v(i,j-1) + erg(6,j-1,i,j,1)) then
c            J dangles over I,J-1.
             j = j - 1
             e = v(i,j)
         elseif (e.eq.v(i+1,j-1)+erg(6,j-1,i+1,i,2)+erg(6,j-1,i+1,j,1))
     .      then
c            Bothe I and J dangle over I+1,J-1.
             i = i + 1
             j = j - 1
             e = v(i,j)
         endif
      endif
 
      if (e.ne.v(i,j)) then
c        Cannot chop away at the ends any more and still the ends do not
c        base-pair with one another. Structure MUST bifucate (OPENL).
         k = i
200      if (k.eq.j) then
c          Structure will not split. Error
           ii = hstnum(i)
           ji = hstnum(j)
           error = 10
           return
         endif
         if (openl.eq.0.and.e.eq.w1(i,k) + w1(k+1,j)) then
c           Best structure on I,J splits into best structures on I,K
c           and K+1,J. Push these fragments on to the stack. (OPENL = 0)
            call push(i,k,w1(i,k),0)
            call push(k+1,j,w1(k+1,j),0)
            goto 100
         else if (openl.eq.1.and.e.eq.w2(i,k) + w2(k+1,j)) then
c           Best structure on I,J splits into best structures on I,K
c           and K+1,J. Push these fragments on to the stack. (OPENL = 1)
            call push(i,k,w2(i,k),1)
            call push(k+1,j,w2(k+1,j),1)
            goto 100
         else
            k = k + 1
            goto 200
         endif
      endif
 
c     Base-pair found. All base-pairs are stored in the range 1 <= I < J <= N.
c     If I and J form a base-pair, then BASEPR(I) = J and BASEPR(J) = I.
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

      openl = 0

c     Check if this is a forced base-pair.
      if (force(i).eq.2.or.force(j).eq.2.or.fce(i,j)) 
     .  nforce = nforce + 1
      if (force(i).eq.2.and.force(j).eq.2) nforce = nforce + 1
c     Perhaps I,J stacks over I+1,J-1?
      if (i.ne.n.and.j.ne.n+1) then
         if (e.eq.erg(2,i,j,i+1,j-1)
     .    + v(i+1,j-1)) then
            i = i + 1
            j = j - 1
            e = v(i,j)
            goto 300
         endif
      endif
 
c     Perhaps I,J closes a hairpin loop?
      if (e.eq.erg(4,i,j,i,j)) goto 100
 
c     E' ( EP in the program ) is E corrcted for possible forced
c     base-pairs.
c
      ep = e
      if (force(i).eq.2.or.force(j).eq.2.or.fce(i,j))
     .   ep = ep - eparam(9)
      if (force(i).eq.2.and.force(j).eq.2) ep = ep - eparam(9)
 
      if (i+2.gt.j-3) then
c        Tidy up loose ends (trivial).
         if (ep.eq.0.or.(i.ne.n.and.ep.eq.erg(6,i,j,i+1,1))) then
 
             goto 100
 
         elseif (j.ne.n+1.and.ep.eq.erg(6,i,j,j-1,2)) then
 
             goto 100
 
         elseif (i.ne.n.and.j.ne.n+1.and.
     .           ep.eq.erg(6,i,j,i+1,1) + erg(6,i,j,j-1,2)) then
 
             goto 100
 
         else
             ii = hstnum(i)
             ji = hstnum(j)
             error = 12
             return
         endif
 
      else if (i.ge.n-1) then
c        Up to one base hanging on to I.
         if (ep.eq.w2(n+1,j-1)) then
           call push(n+1,j-1,w2(n+1,j-1),1)
           goto 100
         elseif (i.ne.n) then
 
           if (ep.eq.erg(6,i,j,i+1,1) + w2(n+1,j-1)) then
              call push(n+1,j-1,w2(n+1,j-1),1)
              goto 100
         else if (ep.eq.erg(6,i,j,i+1,1)+erg(6,i,j,j-1,2) + w2(n+1,j-2)) 
     .     then
              call push(n+1,j-2,w2(n+1,j-2),1)
              goto 100
           endif
 
         elseif (ep.eq.erg(6,i,j,j-1,2) + w2(n+1,j-2)) then
           call push(n+1,j-2,w2(n+1,j-2),1)
           goto 100
 
         else
           ii = hstnum(i)
           ji = hstnum(j)
           error = 12
           return
         endif
 
      else if (j.eq.n+1.or.j.eq.n+2) then
c        Up to one base hanging on to J.
         if (ep.eq.w2(i+1,n)) then
            call push(i+1,n,w2(i+1,n),1)
            goto 100
 
         elseif (ep.eq.erg(6,i,j,i+1,1) + w2(i+2,n)) then
            call push(i+2,n,w2(i+2,n),1)
            goto 100
 
         elseif (j.ne.n+1) then
 
           if (ep.eq.erg(6,i,j,j-1,2)+w2(i+1,n)) then
              call push(i+1,n,w2(i+1,n),1)
              goto 100
          elseif (ep.eq.erg(6,i,j,i+1,1) + erg(6,i,j,j-1,2) + w2(i+2,n)) 
     .      then
              call push(i+2,n,w2(i+2,n),1)
              goto 100
           endif
 
         else
           ii = hstnum(i)
           ji = hstnum(j)
           error = 12
           return
         endif
 
      else
 
         k = i+2
c        Perhaps I,J closes a multi-loop?
400      do while (k.le.j-3)
            if (k.ne.n) then
             if (ep.eq.w1(i+1,k) + w1(k+1,j-1) + eparam(10) + eparam(5))
     .          then
c                 Multi-loop. No dangling ends on I,J.
                  call push(i+1,k,w1(i+1,k),0)
                  call push(k+1,j-1,w1(k+1,j-1),0)
                  goto 100
               else if (ep.eq.erg(6,i,j,i+1,1)+w1(i+2,k)+w1(k+1,j-1) +
     .                   eparam(10) + eparam(6) + eparam(5)) then
c                 Multi-loop. I+1 dangles over the I,J base-pair.
                  call push(i+2,k,w1(i+2,k),0)
                  call push(k+1,j-1,w1(k+1,j-1),0)
                  goto 100
               else if (ep.eq.erg(6,i,j,j-1,2)+w1(i+1,k)+w1(k+1,j-2) +
     .                   eparam(10) + eparam(6) + eparam(5)) then
c                 Multi-loop. J-1 dangles over the I,J base-pair.
                  call push(i+1,k,w1(i+1,k),0)
                  call push(k+1,j-2,w1(k+1,j-2),0)
                  goto 100
              else if (ep.eq.erg(6,i,j,i+1,1)+erg(6,i,j,j-1,2)+w1(i+2,k)
     .            +w1(k+1,j-2)+eparam(10)+2*eparam(6)+eparam(5)) then
c                 Multi-loop. Both I+1 and J-1 dangle over the I,J base-pair.
                  call push(i+2,k,w1(i+2,k),0)
                  call push(k+1,j-2,w1(k+1,j-2),0)
                  goto 100
               endif
            else
               if (ep.eq.w2(i+1,k) + w2(k+1,j-1)) then
c                 Exterior loop. No ends dangling on I,J.
                  call push(i+1,k,w2(i+1,k),1)
                  call push(k+1,j-1,w2(k+1,j-1),1)
                  goto 100
             else if (ep.eq.erg(6,i,j,i+1,1)+w2(i+2,k)+w2(k+1,j-1)) then
c                 Exterior loop. I+1 dangles over the I,J base-pair.
                  call push(i+2,k,w2(i+2,k),1)
                  call push(k+1,j-1,w2(k+1,j-1),1)
                  goto 100
             else if (ep.eq.erg(6,i,j,j-1,2)+w2(i+1,k)+w2(k+1,j-2)) then
c                 Exterior loop. J-1 dangles over the I,J base-pair.
                  call push(i+1,k,w2(i+1,k),1)
                  call push(k+1,j-2,w2(k+1,j-2),1)
                  goto 100
               else if (ep.eq.erg(6,i,j,i+1,1)+
     .                erg(6,i,j,j-1,2)+w2(i+2,k)+w2(k+1,j-2)) then
c                 Exterior loop. Both I+1 and J-1 dangle over the I,J base-pair.
                  call push(i+2,k,w2(i+2,k),1)
                  call push(k+1,j-2,w2(k+1,j-2),1)
                  goto 100
               endif
            endif
            k = k + 1
         enddo
 
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
c    Error, bulge or interior loop not found.
      ii = hstnum(i)
      ji = hstnum(j)
      error = 11
      return
      end
 
 
 
 
c     Store results of a SAVE run for a continuation run.
      subroutine putcont
      include 'rfd.inc'
 
      write(30) n,nsave,vmin,listsz,seqlab
      write(30) stack,tstk,dangle,hairpin,bulge,inter,eparam
      write(30) (vst(i),i=1,n*n)
      write(30) (wst1(i),i=1,n*n)
      write(30) (wst2(i),i=1,n*n)
      write(30) (seq(i),i=nsave(1),nsave(2))
      write(30) ((list(i,j),i=1,listsz),j=1,4)
      write(30) tloop,numoftloops
      write(30) (poppen(i),i=1,4),maxpen,prelog
      return
      end
c     Read results from a SAVE run for a CONTINUATION run.
      subroutine getcont
      include 'rfd.inc'
 
      read(30,err=10) n,nsave,vmin,listsz,seqlab
      read(30,err=10) stack,tstk,dangle,hairpin,bulge,inter,eparam
      read(30,err=10) (vst(i),i=1,n*n)
      read(30,err=10) (wst1(i),i=1,n*n)
      read(30,err=10) (wst2(i),i=1,n*n)
      read(30,err=10) (seq(i),i=nsave(1),nsave(2))
      read(30,err=10) ((list(i,j),i=1,listsz),j=1,4)
      read(30,err=10) tloop,numoftloops
      read(30,err=10) (poppen(i),i=1,4),maxpen,prelog
      goto 11
 
10    call errmsg(40,0,0)
 
11    return
      end


