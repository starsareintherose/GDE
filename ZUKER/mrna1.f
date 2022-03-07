      subroutine process
c     Process RNA sequence to be folded.
      include 'rfd.inc'

c     Selected fragment is from NSAVE(1) to NSAVE(2) in historical
c     numbering.
      do i = nsave(1),nsave(2)
         newnum(i) = 0
      enddo
c     LIST contains information on excisions, and on forced or prohibited
c     base-pairs.
      ptr = 0
100   if (ptr.eq.listsz) goto 400
      ptr = ptr + 1
      if (list(ptr,1).eq.4) goto 200
      if (list(ptr,1).eq.5) goto 300
      goto 100
c     Closed excision beween LIST(PTR,2) and LIST(PTR,3) ( historical
c     numbering ) .
200   do i = list(ptr,2)+4,list(ptr,3)-1
         newnum(i) = 1
      enddo
      goto 100
c     Open excision beween LIST(PTR,2) and LIST(PTR,3) ( historical
c     numbering ) .
300   do i = list(ptr,2),list(ptr,3)
         newnum(i) = 1
      enddo
      goto 100
 
400   n = 0
      do k = nsave(1),nsave(2)
c     Generate new numbering of fragment ( 1 to N ).
         if (newnum(k).eq.0) then
            n = n+1
            newnum(k) = n
         else
c           An excised base gets a new numbering of 0.
            newnum(k) = 0
         endif
      enddo
 
      if (n*2.gt.fldmax) goto 700
 
c     Zero the FORCE and VST arrays.
      do i = 1,n
        force(i) = 0
        if (cntrl(1).ne.2) then
           do j = i,i+n-1
             vst((n-1)*(i-1)+j) = 0
           enddo
        endif
      enddo
 
      do k = nsave(1),nsave(2)
         i = newnum(k)
         if (i.gt.0) then
c           Non-excised bases are examined to determine their type.
c           A - type 1
c           B - an A accessible to nuclease cleavage
c           C - type 2
c           Z - a C accessible to nuclease cleavage
c           G - type 3
c           H - a G accessible to nuclease cleavage
c           U/T - type 4
c           V/W - a U/T accessible to nuclease cleavage
c           anything else - type 5
c           HSTNUM stores historical numbering
c           NUMSEQ stores nucleotide type.
            hstnum(i) = k
            numseq(i) = 5
            if (seq(k) .eq. 'A') numseq(i) = 1
            if (seq(k) .eq. 'B') then
               numseq(i) = 1
               force(i) = 3
            endif
            if (seq(k) .eq. 'C') numseq(i) = 2
            if (seq(k) .eq. 'Z') then
               numseq(i) = 2
               force(i) = 3
            endif
            if (seq(k) .eq. 'G') numseq(i) = 3
            if (seq(k) .eq. 'H') then
               numseq(i) = 3
               force(i) = 3
            endif
            if (seq(k) .eq. 'U'.or.seq(k).eq.'T') numseq(i) = 4
            if (seq(k) .eq. 'V'.or.seq(k).eq.'W') then
               numseq(i) = 4
               force(i) = 3
            endif
         endif
      enddo
 
      ptr = 0
500   if (ptr.eq.listsz) goto 600
      ptr = ptr + 1
      i = list(ptr,2)
      j = list(ptr,3)
      k = list(ptr,4)
      if (list(ptr,1).eq.2.or.list(ptr,1).eq.6) k = j
      goto (500,520,530,540,500,560,570),list(ptr,1)
c     Force bases I to I+K-1 to be double-stranded.
520   do x = i,i+k-1
          force(newnum(x)) = 2
      enddo
      goto 500
c     Force base-pairs I.J , I+1.J-1 , ... I+K-1.J-K+1.
530   do x = 0,k-1
        call sfce(newnum(i+x),newnum(j-x))
      enddo
      goto 500
c     Force the ends of a closed excision to base-pair.
540   call sfce(newnum(i),newnum(j))
      do ii = i+1,i+3
        seq(ii) = ' '
      enddo
      goto 500
c     Prohibit bases I to I+K-1 from base-pairing.
560   do ii = i,i+k-1
          force(newnum(ii)) = 1
      enddo
      goto 500
c     Prohibit the base-pairs I.J , I+1,J-1 , ... I+K-1.J-K+1.
570   if (cntrl(1).ne.2) then
         do x = 0,k-1
           vst((n-1)*(newnum(i+x)-1)+newnum(j-x)) = 1
           vst((n-1)*(newnum(j-x)-1)+newnum(i+x)+n) = 1
         enddo
      endif
      goto 500
c     Double up the sequence.
600   do i = 1,n
        hstnum(i+n) = hstnum(i)
        force(i+n) = force(i)
        numseq(i+n) = numseq(i)
      enddo
700   return
      end
 
 
c     Used in reading the energy files.
      function convt(str)
      implicit integer (a-z)
      character*5 str
      logical neg
 
      neg = .false.
      place = 0
      convt = 0
 
      do i = 5,1,-1
        if (str(i:i).eq.'-') then
          neg = .true.
        else
          if (str(i:i).ge.'0'.and.str(i:i).le.'9') then
             convt = convt + 10**place * (ichar(str(i:i)) - ichar('0'))
             place = place+1
          endif
        endif
      enddo
      if (neg) convt = convt * -1
      return
      end
 
c     Reads energy file names and open the files for reading.
      subroutine enefiles
      character*40 filen
 
10    write (6,*) 'Enter dangle energy file name (default dangle.dat)'
      read (5,100,end=1) filen
      if (filen.eq.'         ') filen = 'dangle.dat'
      open(10,file=filen,status='OLD',err=10)
 
20    write (6,*) 'Enter loop energy file name (default loop.dat)'
      read (5,100,end=1) filen
      if (filen.eq.'         ') filen = 'loop.dat'
      open(11,file=filen,status='OLD',err=20)
 
30    write (6,*) 'Enter stack energy file name (default stack.dat)'
      read (5,100,end=1) filen
      if (filen.eq.'         ') filen = 'stack.dat'
      open(12,file=filen,status='OLD',err=30)
 
40    write (6,*) 'Enter tstack energy file name (default tstack.dat)'
      read (5,100,end=1) filen
      if (filen.eq.'         ') filen = 'tstack.dat'
      open(13,file=filen,status='OLD',err=40)
 
50    write (6,*) 'Enter tloop energy file name (default tloop.dat)'
      read (5,100,end=1) filen
      if (filen.eq.'         ') filen = 'tloop.dat'
      open(29,file=filen,status='OLD',err=50)
 
60    write (6,*) 'Enter misc. loop energy file name (default 
     . miscloop.dat)'
      read (5,100,end=1) filen
      if (filen.eq.'         ') filen = 'miscloop.dat'
      open(32,file=filen,status='OLD',err=60)
 
100   format(a40)
      goto 2
1     call exit(1)
2     return
      end
 
c     Error message subroutine.
      subroutine errmsg(err,i,j)
      include 'rfd.inc'
 
      if (err.eq.10) write (6,10) i,j
      if (err.eq.11) write (6,11) i,j
      if (err.eq.12) write (6,12) i,j
      if (err.eq.20) then
         write (6,20) i,j
         stop
      endif
      if (err.eq.21) then
         write (6,21)
         err = 0
      endif
      if (err.eq.30) write (6,30) i
      if (err.eq.31) then
         write (6,31) sortmax,i,j
         err = 0
      endif
      if (err.eq.40) then
         write (6,40)
         stop
      endif
      return
 
10    format(' Open bifurcation not found between ',i4,' and ',i4)
11    format(' Bulge or interior loop closed by (',i4,',',i4,
     .') not found')
12    format(' Closed bifurcation not found between ',i4,' and ',i4)
20    format(' Base pair between ',i3,' and ',i3,' conflicts with ',
     . 'at least one other pair')
21    format(' Buffer overflow in lineout')
30    format(' End reached at traceback ',i4)
31    format(' More than ',i5,' basepairs in sort at (',i4,',',i4,')')
40    format(' Premature end of save file')
      end
c     Initialize the stack.
      subroutine initst
      implicit integer (a-z)
      dimension stk(50,4)
      common /stk/ stk,sp
 
      sp = 0
      return
      end
c     Add A,B,C,D to the bottom of the stack.
      subroutine push(a,b,c,d)
      implicit integer (a-z)
      dimension stk(50,4)
      common /stk/ stk,sp
 
      sp = sp + 1
      if (sp.gt.50) then
         write (6,*) 'ERROR - STACK OVERFLOW'
         stop
      endif
      stk(sp,1) = a
      stk(sp,2) = b
      stk(sp,3) = c
      stk(sp,4) = d
      return
      end
c     Retrieve A,B,C,D from the bottom of the stack and decrease the
c     stack size by one.
      function pull(a,b,c,d)
      implicit integer (a-z)
      dimension stk(50,4)
      common /stk/ stk,sp
 
      if (sp.eq.0) then
         pull = 1
         return
      endif
      a = stk(sp,1)
      b = stk(sp,2)
      c = stk(sp,3)
      d = stk(sp,4)
      sp = sp - 1
      pull = 0
      return
      end
 
 
c     Line printer output of a secondary structure.
      subroutine linout(n1,n2,energy,iret,jret,error)
c
      include 'rfd.inc'
      character array(6,900),dash,bl,dot
      real energy
      integer unit
c
      data dash/'-'/,bl/' '/,dot/'.'/,amax/900/

      print *,' Output of a secondary structure...'
c
c      WRITE SEQUENCE LABEL AND ENERGY
c
      unit = cntrl(4)
      if(unit.eq.0) then
         print *,' NO PRINTER OUTPUT SELECTED...'
         print *,' ...exiting subroutine outputs.'
         return
      endif
      print *,'unit is ',unit,'.'              
      hstn1 = hstnum(n1)
      hstn2 = hstnum(n2)
      write(unit,103) hstn1,hstn2,seqlab,energy
c
c      INITIALIZE TRACEBACK
c
      call initst
      call push(n1,n2,0,0)
      nstem = 0
      go to 3
c
c      OUTPUT PORTION OF STRUCTURE
c
5     do while (1.eq.1)
         write(unit,106)
         ll = countr
         if(cntrl(3).lt.ll) ll = cntrl(3)
         do k = 1,6
             if(unit.eq.6) write(unit,105) (array(k,i),i = 1,ll)
             if(unit.ne.6) write(unit,104) (array(k,i),i = 1,ll)
         enddo
         if(countr.le.cntrl(3)) go to 3
         do k = 1,5
            do j = 1,6
               array(j,k) = bl
               array(j,k+5) = bl
            enddo
            array(2,k+5) = dot
            array(5,k+5) = dot
         enddo
         k = 10
         ll = cntrl(3)+1
         do i = ll,countr
           k = k+1
           do j = 1,6
             array(j,k) = array(j,i)
           enddo
         enddo
         countr = k
      enddo

c
3     do k = 1,amax
        do j = 1,6
         array(j,k) = bl
        enddo
      enddo
c
c      FILL IN OUTPUT MATRIX
c
      nstem = pull(i,j,countr,xx)
      if (nstem.ne.0) go to 99
c
c      LOOK FOR DANGLING ENDS
c
12    ip = i
      jp = j
      do while (basepr(ip).eq.0)
         ip = ip+1
         if(ip.ge.j) go to 16
      enddo
      do while (basepr(jp).eq.0)
         jp = jp-1
      enddo
      k = ip-i
      if(j-jp.gt.k) k = j-jp
      if(k.eq.0) go to 17
      ii = ip
      jj = jp
      pos = countr+k+1
      if(pos.gt.amax) then
          error = 21
          go to 99
      endif
      do kk = 1,k
         pos = pos-1
         ii = ii-1
         jj = jj+1
         if(ii.ge.i) then
            i2 = hstnum(ii)
            array(2,pos) = seq(i2)
            if(10*(i2/10).eq.i2) call digit(1,i2,pos,amax,array)
         else
            array(2,pos) = dash
         endif
         if(jj.le.j) then
            j2 = hstnum(jj)
            array(5,pos) = seq(j2)
            if(10*(j2/10).eq.j2) call digit(6,j2,pos,amax,array)
         else
            array(5,pos) = dash
         endif
      enddo
      countr = countr+k
      go to 17
c
c      HAIRPIN LOOP
c
16    if(i.ge.j) go to 5
      half = (j-i+2)/2
      ii = i-1
      jj = j+1
      do k = 1,half
         ii = ii+1
         jj = jj-1
         countr = countr+1
         if(countr.gt.amax) then
            error = 21
            go to 99
         endif
         if(seq(hstnum(ii)).eq.' ') go to 40
         i2 = hstnum(ii)
         j2 = hstnum(jj)
         if(10*(i2/10).eq.i2.and.ii.lt.jj) 
     .       call digit(1,i2,countr,amax,array)
         if(10*(j2/10).eq.j2) call digit(6,j2,countr,amax,array)
         if(k.ne.half) then
            array(2,countr) = seq(i2)
            array(5,countr) = seq(j2)
         else
            if(ii.lt.jj) array(3,countr) = seq(i2)
            array(4,countr) = seq(j2)
         endif
22    enddo
      go to 5
c
c      'CLOSED' EXCISION FOUND
c
40    array(3,countr) = dash
      array(4,countr) = dash
      go to 5
c
c      STACKING OR BIFURCATION
c
17    i = ip
      j = jp
      if(basepr(i).eq.j) go to 24
c
c      CHECK FOR KNOT
c
      if(basepr(i).ge.basepr(j).or.i.ge.basepr(i).or.basepr(j).ge.j) 
     . then
         iret = hstnum(i)
         jret = hstnum(basepr(i))
         error = 20
         go to 99
      endif
c
c      BIFURCATION MUST OCCUR
c
      countr = countr+2
      if(countr.gt.amax) then
         error = 21
         go to 99
      endif
      call push(basepr(i)+1,j,countr,0)
      j = basepr(i)
c
c      STACKING REGION
c
24    countr = countr+1
      if(countr.gt.amax) then
         error = 21
         go to 99
      endif
      ii = hstnum(i)
      jj = hstnum(j)
      array(3,countr) = seq(ii)
      array(4,countr) = seq(jj)
      if(10*(ii/10).eq.ii) call digit(1,ii,countr,amax,array)
      if(10*(jj/10).eq.jj) call digit(6,jj,countr,amax,array)
      if(i.eq.iret.and.j.eq.jret) then
         array(2,countr) = '|'
         array(5,countr) = '^'
      end if
      i = i+1
      j = j-1
      if(basepr(i)-j) 12,24,12

99    return
 
103   format(' FOLDING BASES ',i4,' TO ',i4,' OF ',a30,/
     .' ENERGY  =  ',f8.1)
104   format(220a1)
105   format(' ',220a1)
106   format(' ')
      end
 
c     Puts the number COLUMN in row ROW and column POS of the array B.
c     The least significant digit ends up in column POS. If the number
c     is too large to fit, a period is put in column POS and row ROW.
      subroutine digit(row,column,pos,bmax,b)
      implicit integer (a-z)
      integer pos,column,bmax,d(10)
      character*1 b(6,bmax),bl,c(10),dot
      data bl/' '/,c/'0','1','2','3','4','5','6','7','8','9'/,dot/'.'/
c
      size=1
      n=column
1     p=n/10
      q=n-10*p
      d(size)=q
      if(p.eq.0) go to 2
      n=p
      size=size+1
      go to 1
2     if(pos-size.lt.0) go to 3
      do k=1,size
         q=pos-k+1
         if(b(row,q).ne.bl) go to 3
      enddo
      p=pos
      do 4 k=1,size
      q=d(k)
      b(row,p)=c(q+1)
4     p=p-1
      return
3     b(row,pos)=dot
      return
      end
 
c     Generates a region table for the Shapiro and Maizel DRAW program.
      subroutine regtab
      include 'rfd.inc'
      real r
 
      k = 1
      region = 1
      do while (k.lt.n)
         r = 0.0
         regsz = 1
         kst = k
         if (k.lt.basepr(k)) then
           do while (basepr(k+1).eq.basepr(k)-1.and.k.lt.n)
              regsz = regsz + 1
              r = r + float(erg(2,k,basepr(k),k+1,basepr(k+1))) / 10.0
              k = k + 1
           enddo
           write (22,100) region,hstnum(kst),hstnum(basepr(kst)),regsz,r
           region = region + 1
         endif
         k = k + 1
      enddo
      return
 
100   format(' (',i5,')',3x,3(i5,3x),f7.1)
      end
 
c     Generates a CT file. (Richard Feldmann)
      subroutine ct(r)
      include 'rfd.inc'
      real r

      write(21,100) n,r,seqlab
      do k = 1,n
        k1 = k+1
        if (k.eq.n) k1 = 0
        write (21,200) k, seq(hstnum(k)),k-1,k1,basepr(k),hstnum(k)
      enddo
      return
 
100   format(i5,1x,'ENERGY = ',f7.1,4x,a30)
200   format(i5,1x,a1,3x,4i5)
      end
 
 
 
c     Menu subroutine for RNA folding program.
c     Allows the user to set energy parameters and to
c     add auxiliary information.
      subroutine menu
 
      include 'rfd.inc'
c      data listsz/0/   rtm 11.II.99
      listsz=0

 
10    if (listsz.ge.100) goto 800
      write (6,900)
50    write (6,901)
      read  (5,*,end=1,err=1) choice
      if (choice.lt.1.or.choice.gt.10) goto 50
      goto (100,200,300,400,400,200,300,800,60,70),choice
 
60    call listout(6)
      goto 10
 
70    listsz = 0
      goto 10
1     call exit(1)
 
 
100   write (6,1000) (eparam(i),i=1,10)
101   write (6,1001)
      read (5,1002,end=10,err=10) parm
      if (parm.lt.1.or.parm.gt.10) goto 10
      write (6,1003)
      read (5,*,end=10,err=10) val
      eparam(parm) = val
      goto 100
1000  format(/,
     .  10x,'   Energy Parameters (10ths kcal/mole)',//,
     .  10x,' 1 Extra stack energy                        [',i5,']',/,
     .  10x,' 2 Extra bulge energy                        [',i5,']',/,
     .  10x,' 3 Extra loop energy (interior)              [',i5,']',/,
     .  10x,' 4 Extra loop energy (hairpin)               [',i5,']',/,
     .  10x,' 5 Extra loop energy (multi)                 [',i5,']',/,
     .  10x,' 6 Multi loop energy/single-stranded base    [',i5,']',/,
     .  10x,' 7 Maximum size of interior loop             [',i5,']',/,
     .  10x,' 8 Maximum lopsidedness of an interior loop  [',i5,']',/,
     .  10x,' 9 Bonus Energy                              [',i5,']',/,
     .  10x,'10 Multi loop energy/closing base-pair       [',i5,']',//)
1001   format(' Enter Parameter to be changed (<return> for main menu) '
     .,$)
1002   format(i6)
1003   format(' Enter new value  ',$)
 
200   write (6,2001)
      read (5,*,end=10,err=10) i,k
      listsz = listsz + 1
      list(listsz,1) = choice
      list(listsz,2) = i
      list(listsz,3) = k
      list(listsz,4) = -1
      goto 10
2001  format(' Enter base and length  ',$)
 
300   write (6,3001)
      read (5,*,end=10,err=10) i,j,k
      listsz = listsz + 1
      list(listsz,1) = choice
      list(listsz,2) = i
      list(listsz,3) = j
      list(listsz,4) = k
      goto 10
3001  format(' Enter base pair and length    ',$)
 
400   write (6,4001)
      read (5,*,end=10,err=10) i,j
      listsz = listsz + 1
      list(listsz,1) = choice
      list(listsz,2) = i
      list(listsz,3) = j
      list(listsz,4) = -1
      goto 10
4001  format(' Enter begining and end    ',$)
 
800   return
 
900   format(/,
     .  10x,'1  Energy Parameter           6  Single Prohibit',/,
     .  10x,'2  Single Force               7  Double Prohibit',/,
     .  10x,'3  Double Force               8  Begin Folding  ',/,
     .  10x,'4  Closed Excision            9  Show current   ',/,
     .  10x,'5  Open Excision             10  Clear current  ',//)
901   format(' Enter Choice   ',$)
      end
