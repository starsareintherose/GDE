      subroutine listout(u)
c     This subroutine lists current choices on excisions and on
c     forced or prohibited base-pairs.
      integer u
      common /list/ list,listsz
      dimension list(100,4)
      character*20 choices(7)
      data choices/'Energy Parameter   ','Single Force       ',
     .'Double Force       ','Closed Excision    ','Open Excision      ',
     .'Single Prohibit    ','Double Prohibit    '/
 
      if (listsz.eq.0) then
         write(u,*) ' No choices currently defined'
      else
         write(u,*) ' '
         write(u,*) ' Current Choices'
         do 100 i = 1,listsz
             if (list(i,1).eq.3.or.list(i,1).eq.7) then
                write(u,1000) choices(list(i,1)),(list(i,k),k = 2,4)
             else
                write(u,1001) choices(list(i,1)),(list(i,k),k = 2,3)
             endif
100      continue
         write(u,*) ' '
      endif
      return
1000  format(10x,a20,': (',i4,',',i4,') ',i4)
1001  format(10x,a20,':  ',i4,',',i4)
      end
 
c     Control subroutine for RNA folding.
      subroutine device
      include 'rfd.inc'
      character*40 sfile,str
      logical used
      used = .false.
 
c     What kind of run is this? ( regular, save or continuation )
      write (6,2000)
      read (5,2001,end=1) cntrl(1)
      write (6,*) ' '
      if (cntrl(1).lt.0.or.cntrl(1).gt.2) cntrl(1) = 0
 
      if (cntrl(1).eq.1) then
         cntrl(7) = 1
      else
c        What mode is the program to be run in?
c        dot plot, automatic sorted tracebacks of one sequence
c        fragment or suboptimal foldings of every complete
c        sequence in a file.
9        write (6,1002)
         read (5,2001,end=1) cntrl(7)
         if (cntrl(7).lt.0.or.cntrl(7).gt.2) cntrl(7) = 0
         if (cntrl(1).eq.2.and.cntrl(7).eq.2) then
            write (6,*)    
     .'Combination of continuation run and multiple foldings disallowed'
            write (6,*) ' '
            goto 9
         endif
         write (6,*) ' '
      endif
c     Folding multiple sequences is treated as a sort run.
c     Find total number of sequences to be folded in a multiple
c     sequence run.
 
         if (cntrl(7).eq.2) call mseq(cntrl(5))
 
         if (cntrl(7).eq.0) then
c           Prompt for terminal type. NOTE: NOT NEEDED ON IRIS VERSION
c$$$            write (6,1000)
c$$$            read (5,2001,end=7,err=7) cntrl(5)
c$$$7           if (cntrl(5).lt.1.or.cntrl(5).gt.3) cntrl(5) = 2
            write (6,1001)
         elseif (cntrl(1).ne.1) then
c           Prompt for controls on sort.
            write (6,1004)
            read (5,2001,end=1) cntrl(8)
            write (6,1003)
         endif
         if (cntrl(1).ne.1) then
            read (5,2001,end=1) cntrl(6)
            if (cntrl(6).lt.1) cntrl(6) = 1
 
            write (6,1005)
            read (5,2001,end=1) cntrl(9)
            write (6,*) ' '
         endif
c        Prompt for SAVE file name for a save/continuation run.
         if (cntrl(1).ne.0) then
4           write (6,3000)
            read (5,3001,end=1) sfile
            if (cntrl(1).eq.1) then
                str = 'unknown'
            else
                str = 'OLD'
            endif
            if (sfile.eq.'         ') sfile= 'fold.sav'
            open(30,err=2,file=sfile,status=str,form='UNFORMATTED')
            goto 3
2           if (cntrl(1).eq.2) goto 4
3           if (cntrl(1).eq.2) call getcont
         endif
c        Obtain sequence. Original length is N.
c        A fragment from NSAVE(1) to NSAVE(2) is selected.
c        After PROCESS, N becomes the length of the processed sequence
c        to be folded.
         if (cntrl(1).ne.2.and.cntrl(7).ne.2) then
            call formid(seqlab,seq,n,maxsiz,used)
            write(6,4000) seqlab,n
            write (6,*) 'Enter start of fragment (default 1)'
            read (5,4001,end=1) nsave(1)
            if (nsave(1).le.0) nsave(1) = 1
            write (6,4002) n
            read (5,4001,end=1) nsave(2)
            if (nsave(2).le.0) nsave(2) = n
         endif

 
1000  format(1x,'Enter terminal type',/,5x,'1  VGT100',
     .   /,5x,'2  Visual 102 (default)',/,5x,'3  Tektronics 4105')
1001  format(/,' Enter minimum vector size for plot (default 1)  ',$)
1002  format(1x,'Enter run mode',/,5x,'0  Sub-optimal plot (default)',
     .   /,5x,'1  N Best',/,5x,'2  Multiple Molecules')
1003  format(/,' Enter number of tracebacks (default 1)  ',$)
1004  format(/,' Enter percentage for sort (default 0)  ',$)
1005  format(/,' Enter window size (default 0)  ',$)
2000  format(1x,'Enter run type',/,5x,'0  Regular run (default)',
     .   /,5x,'1  Save run',/,5x,'2  Continuation run')
2001  format(i6)
3000  format(' Enter save file name (default fold.sav)')
3001  format(a30)
4000  format(/,' ',a30,5x,i5,' nucleotides',/)
4001  format(i10)
4002  format(1x,'Enter end of fragment (default ',i5,')')
      return
1     call exit(1)
      end
 
c     Obtain multiple sequences from a sequence file using MULTID.
      subroutine mseq(i)
      include 'rfd.inc'
      logical used
      data used/.false./
 
      if (.not.used) then
         call multid(seqlab,seq,n,maxsiz,used,i)
         write (6,*) ' '
      else
         call multid(seqlab,seq,n,maxsiz,used,i)
         write(6,4000) seqlab,n
      endif
      nsave(1) = 1
      nsave(2) = n
      return
4000  format(/,' ',a30,5x,i5,' nucleotides',/)
      end
 
c     Set up output units and files for RNA folding.
      subroutine outputs
      include 'rfd.inc'
      character*40 str,dstr
      character*1  in
 
      str(1:1) = ' ';

c     Examine sequence label to get default names for output files.
 
      k = 1
      do while ((seqlab(k:k).lt.'A'.or.seqlab(k:k).gt.'Z').and.
     .          (seqlab(k:k).lt.'a'.or.seqlab(k:k).gt.'z'))
         k = k + 1
      enddo
      slen = min0(30,25+k)
      do while (seqlab(slen:slen).eq.' ')
         slen = slen - 1
      enddo
      j = 0
      do i = k,slen
         j = j + 1
         if ((seqlab(i:i).ge.'A'.and.seqlab(i:i).le.'Z').or.
     .       (seqlab(i:i).ge.'a'.and.seqlab(i:i).le.'z').or.
     .       (seqlab(i:i).ge.'0'.and.seqlab(i:i).le.'9')) then
 
            dstr(j:j) = seqlab(i:i)
         else
            dstr(j:j) = '_'
         endif
      enddo
      slen = j
 
c     Line printer output. Get name and open file for write.
      cntrl(2) = 0
      write (6,5010)
      read (5,5000,end=1) in
      if (in.ne.'N'.and.in.ne.'n') then
          cntrl(2) = 1
          write (6,5011)
          read (5,5000,end=1) in
          if (in.eq.'N'.or.in.eq.'n') then
             dstr(slen+1:slen+4) = '.out'
51           write (6,5012) dstr(1:slen+4)
             read (5,5001) str
             if (str(1:1).eq.' ') str = dstr(1:slen+4)
             cntrl(4) = 20
c             open(20,file=str,recl=255,status='unknown',err=51)
             open(20,file=str,status='unknown',err=51)
c vf90: Warning, line 203: RECL with ACCESS=SEQUENTIAL could not be translated.
c vf90: Warning, line 203: Specifier removed, but may not yield same results. (RECL)
          else
             cntrl(4) = 6
          endif
          write (6,5013)
          read (5,5014,end=1) cntrl(3)
          if (cntrl(3).eq. 0) cntrl(3) = 80
      endif
 
c     CT file output. Get name and open file for write.
      write (6,5020)
      read (5,5000,end=1) in
      if (in.eq.'Y'.or.in.eq.'y') then
         cntrl(2) = 2 + 2*cntrl(2)
         dstr(slen+1:slen+3) = '.ct'
52      write (6,5021) dstr(1:slen+3)
         read (5,5001) str
         if (str(1:1).eq.' ') str = dstr(1:slen+3)
         open(21,file=str,status='unknown',err=52)
      endif
 
c     Region table output. Get name and open file for write.
      write (6,5030)
      read (5,5000,end=1) in
      if (in.eq.'Y'.or.in.eq.'y') then
         if (cntrl(2).eq.1.or.cntrl(2).eq.2) cntrl(2) = cntrl(2) + 1
         cntrl(2) = cntrl(2) + 3
         dstr(slen+1:slen+4) = '.reg'
53       write (6,5031) dstr(1:slen+4)
         read (5,5001) str
         if (str(1:1).eq.' ') str = dstr(1:slen+4)
c rtm 11.II.98 :  52 below -> 53
         open(22,file=str,status='unknown',err=53)
      endif
      write (6,*) ' '
      return 
1     call exit(1)
 
5000  format(a1)
5001  format(a40)
5010  format(' Do you want printer output? (Y,n) ',$)
5011  format(' Output to terminal? (Y,n) ',$)
5012  format(' Enter output file name (default ',a,')')
5013  format(' Enter number of columns on printer (default 80)  ',$)
5014  format(i10)
5020  format(' Do you want ct file? (y,N) ',$)
5021  format(' Enter ct file name (default ',a,')')
5030  format(' Do you want region table? (y,N) ',$)
5031  format(' Enter region table file name (default ',a,')')
 
      end
 
c     Reads energy files.
      subroutine ergread
 
      include 'rfd.inc'
      logical  endfile
      logical find
      character*80 inrec
      character*5 temp
      real a,b,c,d
      integer*2 convt
 
c     TLoop INFORMATION IN
      call gettloops
 
c     Get misc loop info
      if(find(32,3,' > ')) stop 'Premature end of MISCLOOP.DAT'
      read (32,*) prelog
      prelog=prelog*10
      endfile = find(32,3,' > ')
      read (32,*) a
      maxpen=int(a*10)
      endfile = find(32,3,' > ')
      read (32,*) a,b,c,d
      poppen(1)=int(a*10)
      poppen(2)=int(b*10)
      poppen(3)=int(c*10)
      poppen(4)=int(d*10)
      endfile = find(32,3,' > ')
c     Set default values of eparam.
      eparam(1) = 0
      eparam(2) = 0
      eparam(3) = 0
      eparam(4) = 0
      eparam(7) = 30
      eparam(8) = 30
      eparam(9) = -500
      read (32,*) a,b,c
      eparam(5)=int(a*10)
      eparam(6)=int(b*10)
      eparam(10)=int(c*10)
 
c     DANGLE IN
 
      do a = 1,5
        do b = 1,5
          do c = 1,5
            do d = 1,2
              dangle(a,b,c,d) = 0
            enddo
          enddo
        enddo
      enddo
      endfile = find(10,3,'<--')
      if (.not.endfile) then
         do var4 = 1,2
            do var1 = 1,4
               if (endfile) goto 150
               read(10,100,end=150) inrec
               do var2 = 1,4
                  do var3 = 1,4
                     j = 0
                     tstart = (var2-1)*20 + (var3-1)*5 + 1
                     temp = inrec(tstart:tstart+4)
                     do i = 2,4
                        if (temp(i-1:i+1).eq.' . ') j = infinity
                     enddo
                     if (temp(1:1).eq.'.'.or.temp(5:5).eq.'.') j = infinity
                     if (j.eq.0) j = convt(temp)
                     if(j.ne.infinity) dangle(var1,var2,var3,var4) = j
                  enddo
               enddo
               endfile = find(10,3,'<--')
            enddo
         enddo
      else
         write (6,*) 'ERROR - DANGLE ENERGY FILE NOT FOUND'
         stop
      endif
 
100   format(a80)
      goto 200
 
150   write (6,*) 'ERROR - PREMATURE END OF DANGLE ENERGY FILE'
      stop
 
 
c    INTERNAL,BULGE AND HAIRPIN IN
 
200   endfile = find(11,5,'-----')
      i = 1
201   read(11,100,end=300) inrec
      j = -1
      do ii = 1,3
         j = j + 6
         do while (inrec(j:j).eq.' ')
            j = j + 1
         enddo
         temp = inrec(j:j+4)
         k = 0
         do jj = 2,4
            if (temp(jj-1:jj+1).eq.' . ') k = infinity
         enddo
         if (temp(1:1).eq.'.'.or.temp(5:5).eq.'.') k = infinity
         if (k.eq.0) k = convt(temp)
         if (ii.eq.1) inter(i) = k
         if (ii.eq.2) bulge(i) = k
         if (ii.eq.3) hairpin(i) = k
      enddo
      i = i + 1
      if (i.le.30) goto 201
 
c     STACK IN
 
300   do a = 1,5
        do b = 1,5
          do c = 1,5
            do d = 1,5
              stack(a,b,c,d) = infinity
            enddo
          enddo
        enddo
      enddo
      endfile = find(12,3,'<--')
      if (.not.endfile) then
         do var1 = 1,4
            do var3 = 1,4
               if (endfile) goto 350
               read(12,100,end=350) inrec
               do var2 = 1,4
                  do var4 = 1,4
                     j = 0
                     tstart = (var2-1)*20 + (var4-1)*5 + 1
                     temp = inrec(tstart:tstart+4)
                     do i = 2,4
                        if (temp(i-1:i+1).eq.' . ') j = infinity
                     enddo
                     if (temp(1:1).eq.'.'.or.temp(5:5).eq.'.') j = infinity
                     if (j.eq.0) j = convt(temp)
                     stack(var1,var2,var3,var4) = j
                  enddo
               enddo
            enddo
            endfile = find(12,3,'<--')
         enddo
      else
         write (6,*) 'ERROR - STACK ENERGY FILE NOT FOUND'
         stop
      endif
      call stest(stack,'STACK ')
 
      goto 400
 
350   write (6,*) 'ERROR - PREMATURE END OF STACK ENERGY FILE'
      stop
 
400   do a = 1,5
        do b = 1,5
          do c = 1,5
            do d = 1,5
              tstk(a,b,c,d) = infinity
            enddo
          enddo
        enddo
      enddo
      endfile = find(13,3,'<--')
      if (.not.endfile) then
         do var1 = 1,4
            do var3 = 1,4
               if (endfile) goto 350
               read(13,100,end=450) inrec
               do var2 = 1,4
                  do var4 = 1,4
                     j = 0
                     tstart = (var2-1)*20 + (var4-1)*5 + 1
                     temp = inrec(tstart:tstart+4)
                     do i = 2,4
                        if (temp(i-1:i+1).eq.' . ') j = infinity
                     enddo
                     if (temp(1:1).eq.'.'.or.temp(5:5).eq.'.') j = infinity
                     if (j.eq.0) j = convt(temp)
                     tstk(var1,var2,var3,var4) = j
                  enddo
               enddo
            enddo
            endfile = find(13,3,'<--')
         enddo
      else
         write (6,*) 'ERROR - STACK ENERGY FILE NOT FOUND'
         stop
      endif
c**   CALL STEST(TSTK,'TSTACK')
 
      close(10)
      close(11)
      close(12)
      close(13)
      goto 500
 
450   write (6,*) 'ERROR - PREMATURE END OF TSTACK ENERGY FILE'
      stop
 
500   return
      end
c     Symmetry test on stacking and terminal stacking energies.
c     For all i,j,k,l between 1 and 4, STACK(i,j,k,l) MUST equal
c     STACK(l,k,j,i). If this fails at some i,j,k,l; these numbers
c     are printed out and the programs grinds to an abrupt halt!
      subroutine stest(stack,sname)
      integer stack(5,5,5,5),a,b,c,d
      character*6 sname
 
      do a = 1,4
        do b = 1,4
          do c = 1,4
            do d = 1,4
              if (stack(a,b,c,d).ne.stack(d,c,b,a)) then
                 write (6,*) 'SYMMETRY ERROR'
                 write (6,101) sname,a,b,c,d,stack(a,b,c,d)
                 write (6,101) sname,d,c,b,a,stack(d,c,b,a)
                 stop
              endif
            enddo
          enddo
        enddo
      enddo
      return
101   format(5x,a6,'(',3(i1,','),i1,') = ',i10)
      end
 
c     Writes out the numbers in the energy arrays of the folding program.
      subroutine out(u)
      include 'rfd.inc'
      integer*2 tlptr,bptr,nbase
      integer key
      character*4 tlbuf
 
c     used for testing contents of energy arrays only
c     not used in the mature program
 
      write (u,100) 'DANGLE'
      do var4 = 1,2
        do var1 = 1,4
          do var2 = 1,4
            do var3 = 1,4
               o = dangle(var1,var2,var3,var4)
               if (o.ne.infinity) then
                  write (u,101) o
               else
                  write (u,102)
               endif
            enddo
          enddo
          write (6,103)
        enddo
        write (6,104)
      enddo
 
 
      write (u,100) 'TSTACK'
      do var1 = 1,4
        do var3 = 1,4
          do var2 = 1,4
            do var4 = 1,4
               o = tstk(var1,var2,var3,var4)
               if (o.ne.infinity) then
                  write (u,101) o
               else
                  write (u,102)
               endif
            enddo
          enddo
          write (6,103)
        enddo
        write (6,104)
      enddo
 
      write (u,100) 'STACK'
      do var1 = 1,4
        do var3 = 1,4
          do var2 = 1,4
            do var4 = 1,4
               o = stack(var1,var2,var3,var4)
               if (o.ne.infinity) then
                  write (u,101) o
               else
                  write (u,102)
               endif
            enddo
          enddo
          write (6,103)
        enddo
        write (6,104)
      enddo
 
      write (u,200) 'INTER','BULGE','HAIRPIN'
      do i = 1,30
        write (u,201) i,inter(i),bulge(i),hairpin(i)
      enddo
 
      write (u,100) 'TLoops'
      do tlptr=1,numoftloops
      key=tloop(tlptr,1)
         do bptr=1,4
            nbase=mod(key, 8)
            key=int(key/8)
            if (nbase.eq.1) then
               tlbuf(bptr:bptr)='A'
            elseif (nbase.eq.2) then
               tlbuf(bptr:bptr)='C'
            elseif (nbase.eq.3) then
               tlbuf(bptr:bptr)='G'
            else
               tlbuf(bptr:bptr)='U'
           endif
         enddo
      write (u,205) tlbuf,tloop(tlptr,2)
      enddo
      return
100   format(//,a40,//)
101   format('+',i4,1x,$)
102   format('+',4('*'),1x,$)
103   format(' ')
104   format(/)
200   format(3a20,/,60('-'),/)
201   format(i4,i16,2i20)
205   format(a4,2x,i8)
      end
 
c     Used in reading the energy files.
c     Locates markers in the energy files so that data can be read
c     properly.
      function find(unit,len,str)
      implicit integer (a-z)
      logical find,flag
      character*20 str
      character*80 inrec
 
      find = .false.
      flag = .false.
      do  while(.not.flag)
         read(unit,100,end=200) inrec
         count = 1
         do 101 i = 1,80-len+1
           if (inrec(i:i).eq.str(count:count)) then
             count = count + 1
             if (count.gt.len) flag = .true.
             if (inrec(i+1:i+1).ne.str(count:count)) count = 1
           endif
101      continue
      enddo
 
      return
100   format(a80)
200   find = .true.
      return
      end
 
 
      subroutine cdump
      include 'rfd.inc'
      character*40 name
      character yn
 
      write (6,*)
     . 'Enter file name for continuation dump (return for terminal)'
      read (5,100,end=1) name
      if (name.eq.' ') then
         u = 6
      else
         u = 31
         open(31,status='unknown',file=name)
      endif
      call listout(u)
      write (u,101) 'Energy Parameters'
      write (u,1000) eparam
      write (6,*) 'Listing of energy files? (y/N)'
      read(5,102) yn
      if (yn.eq.'Y'.or.yn.eq.'y') then
         call out(u)
      endif
      return
1     call exit(1)
100   format(a30)
101   format(a20,/)
102   format(a1)
1000  format(/,
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
      end
      subroutine gettloops
c
c* Read in TLoop sequences, convert to numeric form, and
c* convert energy to an integer (*10)
c
      parameter (maxtloops=40,bufsiz=80)
      integer*2 i,ptr,tloop(maxtloops,2),nseq(4)
      integer*2 numoftloops
      integer*2 convt
      character*5 buffa
      character*80 inbuf
c
      common /tloops/tloop,numoftloops
c
      numoftloops=0
c
c* Throw out header
c
      read (29,1)
c
c* Read a line and convert to numeric sequence and energy until EOF
c
10    read (29,2,end=99) inbuf
        ptr=1
        numoftloops=numoftloops+1
        do while ((ptr.lt.bufsiz).and.
     2             (inbuf(ptr:ptr).eq.' '))
                ptr=ptr+1
        enddo
c Only take first four characters, since they're TETRAloops
        buffa(1:4)=inbuf(ptr:ptr+3)
        buffa(5:5)=' '
        call tonum(buffa,nseq)
        tloop(numoftloops,1)=((nseq(4)*8+nseq(3))*8+nseq(2))*8+nseq(1)
        ptr=ptr+4
        do while ((ptr.lt.bufsiz).and.
     2             (inbuf(ptr:ptr).eq.' '))
                ptr=ptr+1
        enddo
c Simple error czeck.
        if (inbuf(ptr+4:ptr+4).ne.' ') then
                write (*,5) inbuf
        endif
        buffa(1:4)=inbuf(ptr:ptr+3)
        buffa(5:5)=' '
        tloop(numoftloops,2)=convt(buffa)
        do i=1,ptr+4
                inbuf(i:i)=' '
        enddo
      goto 10
c
c* Normal ending
c
99    close(unit=29,status='KEEP')
      return
1     format(/)
2     format (a)
5     format (1x,'Too many characters in numeric field of this line of',/,
     1        1x,'tloop.dat file: ',a)
      end
 
      subroutine tonum(tloopseq,numeric)
c
c* Convert TLoopSeq to numeric format in Numeric.
c
      character*5 tloopseq
      integer*2 i,numeric(4)
c
      do i=1,4
        if     (tloopseq(i:i).eq.'A') then
                numeric(i)=1
        elseif (tloopseq(i:i).eq.'C') then
                numeric(i)=2
        elseif (tloopseq(i:i).eq.'G') then
                numeric(i)=3
        elseif (tloopseq(i:i).eq.'U') then
                numeric(i)=4
        elseif (tloopseq(i:i).eq.'T') then
                numeric(i)=4
        else
                write (*,1) tloopseq(i:i)
1               format (1x,'Unknown base in TLOOP file: ',a1)
                stop
        endif
      enddo
      return
      end
