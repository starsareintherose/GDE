      subroutine sortout(i,j,rep,err)
      include 'rfd.inc'
      logical mark
c     The first time in (REP = 1), valid I,J base-pairs are
c     sorted by energy.
      err = 0
      if (rep.eq.1.or.cntrl(7).eq.2) then
         call build_heap
         call heap_sort
         cntr = num
      endif
c     Select the next valid unmarked base-pair
      do while (mark(heapi(cntr),heapj(cntr)))
         if (cntr.eq.1) then
            err = 30
            return
         endif
         cntr = cntr - 1
      enddo
c     The base-pair I,J will be used to create a folding.
      i = heapi(cntr)
      j = heapj(cntr)
 
      return
      end
 
 
c     Add I,J to HEAPI and HEAPJ if the best energy of a folding containing
c     I,J is no greater than a given percent ( CNTRL(8) ) of the minimum
c     folding energy.
      subroutine build_heap
      include 'rfd.inc'
 
      crit = vmin + abs(vmin)*cntrl(8)/100
 
      num = 0
      i = 1
      j = 2
      do while (i.lt.n)
         if (ene(i,j).le.crit) then
            if (num.eq.sortmax) then
               err = 31
               call errmsg(err,hstnum(i),hstnum(j))
               goto 10
            endif
            num =  num + 1
            heapi(num) = i
            heapj(num) = j
            j = j + cntrl(9) + 1
            if (j.gt.n) then
               i = i + 1
               j = i + 1
            endif
         else
            j = j +1
            if (j.gt.n) then
               i = i + 1
               j = i + 1
            endif
         endif
      enddo
 
      do i = num+1,sortmax+1
        heapi(i) = 0
        heapj(i) = 0
      enddo
 
10    do q = 2,num
         cur = q
         up = cur/2
         do while
     . (ene(heapi(cur),heapj(cur)).lt.ene(heapi(up),heapj(up)).and.
     .up.ge.1)
              call swap(heapi(cur),heapi(up))
              call swap(heapj(cur),heapj(up))
              cur = cur/2
              up = cur/2
         enddo
      enddo
 
      return
      end
 
 
c     Efficient sort of heap.
      subroutine heap_sort
      include 'rfd.inc'
 
      do ir = num-1,2,-1
         call swap(heapi(ir+1),heapi(1))
         call swap(heapj(ir+1),heapj(1))
 
         up = 1
         c = 2
         do while (c.le.ir)
           if (c.ne.ir) then
             if (ene(heapi(c+1),heapj(c+1)).lt.ene(heapi(c),heapj(c)))
     .         then
                c = c + 1
             endif
           endif
           if (ene(heapi(c),heapj(c)).lt.ene(heapi(up),heapj(up))) then
              call swap(heapi(c),heapi(up))
              call swap(heapj(c),heapj(up))
              up = c
              c = 2 * c
           else
              c = ir+1
           endif
         enddo
      enddo
      return
      end
 
c     ENE(k) is the minimum energy of a folding containing the base-pair
c     I,J at heap(k).
      function ene(i,j)
      include 'rfd.inc'
 
      ene = v(i,j)+ v(j,i+n)
 
      return
      end
 
 
 
      subroutine swap(i,j)
 
       k = i
       i = j
       j = k
 
      return
      end
