      subroutine formid(seqid,seq,nseq,nmax,used)
c     RESEARCHER: M. ZUKER
c     JUN 1986
c     DS TUDHOPE
c     WRITTEN TO WORK ON VAX-11-750 UNDER VMS4.3 USING STANDARD F-77
c     ####  NOTE #### USE OF CARRIAGE CONTROL LINE FEED SUPPRESSION
c     IS USED IN THIS PROGRAM AND IS NOT STANDARD
c     FORTRAN-77
c     SEE '$' IN FORMAT STATEMENT 102
      implicit integer*2 (a-z)
      logical found,endfil,valid,used
      integer*2 sline(500),seqnum(500)
      integer*4 nseq,nmax
      character*8  stype
      character*30 seqids(500)
      character*30 choice,seqid
      character*50 filnam,fmtseq(500)
      character*80 reclin
      character*1  seq(nmax)
      data   found/.false./,
     .     endfil/.false./,
     .     valid/.false./
c     **************************************************************
c     subroutine number FS372, test main driver FP475
c     biomath program catalogued as FS372.FOR named formid.for
c     **************************************************************
c     PURPOSE:
c     subroutine used to extract sequences from various
c     format type files including
c     STANFORD,GENBANK,EMBL,PIR,and,NRC
c     **************************************************************
c     variable list table:
c     * -- sent down from main returned unchanged
c     ** - returned to main from subroutine
c     ***- sent down from main and returned changed
c     INTEGERS:IDCNT   -- to keep track of the number of sequence
c     identifiers found in the file
c     I,K     -- loop counters
c     LINE    -- to point to the record line number of a file
c     N       -- counter to extract the correct number of
c     sequence elements
c     *   NMAX    -- maximum length of sequence expected by user
c     **   NSEQ    -- length of the sequence retrieved from file
c     POINTR  -- pointer to point to sequence identifier
c     chosen from array of identifiers
c     SEQNUM  -- an ARRAY of the length of the sequences
c     for the sequence identifiers found in the
c     NRC-format type files
c     SLINE   -- an ARRAY of the record-line numbers where
c     a sequence starts in a file
c     START   -- defines what column of a record to start
c     reading the sequence in.
c     CHARACTERS:CHOICE-- character string of length 20 to read in
c     choice for sequence identifier to retrieve
c     FILNAM  -- character string of length 50 to read in
c     the filename.
c     FMTSEQ  -- an ARRAY of characters, each length 50,
c     describing how the sequence is to be read
c     RECLIN  -- a record of a file, length 80 characters
c     ** SEQ     -- an ARRAY of characters each length 1 to
c     store the sequence
c     ** SEQID   -- retrieved name of the sequence
c     SEQIDS  -- an ARRAY of retrieved names of sequences
c     STYPE   -- character string of length 8 defining the
c     format type of the file.
c     LOGICALS:FOUND   -- logical variable used in looping until
c     something retrieved
c     *** USED    -- logical variable used in determining if the
c     subroutine has been previously called
c     VALID   -- logical variable used in looping until
c     some input valid
c     
c     
c     if the subroutine has NOT been USED then must
c     input the filename and do error checking on filename
c     else branch to listing the sequences available in this file
      if (.not.used) then
 10      found = .false.
         valid = .false.
         used  = .false.
c     initialize for a new file by setting variables back
c     to zero and blanking out old names
         if (nmax.eq.0) nmax = 9999
         nseq = 0
         stype = '        '
         seqid = '                    '
         filnam= '                                                  '
         do 500 i = 1,500
            sline(i) = 0
            seqnum(i)= 0
            seqids(i)='                    '
            fmtseq(i)='                                             '
 500     continue
 111     write(6,102) 'Input sequence file name (/ to end)  '
 102     format(1x,a37,$)
         read(5,110,end=999,err=111) filnam
         if (filnam(1:1).eq.'/') then
            goto 999
         endif
 110     format(a50)
c     open the file only after a valid filename has been retrieved
c     error in filename results in prompting for the input again
         open(66,file=filnam,status='OLD',err=10)
c     find sequence file format type and the sequence identifiers
         idcnt = 0
         line = 1
c     DO WHILE (.NOT.FOUND)
 600     read(66,120,end = 410,err = 991) reclin
 120     format(a80)
c     STANFORD format, recognized by the ';' in the first column
c     of the record
         if (reclin(1:1).eq.';') then
c     DO WHILE (.NOT.ENDFIL)
 610        line = line + 1
            read(66,120,end = 410,err=991) reclin
c     to find the next sequence identifer scroll through the
c     file until the first character in the line is not ';'
            if (reclin(1:1).ne.';') then
               found = .true.
               stype ='STANFORD'
               idcnt = idcnt + 1
               sline(idcnt)= line
               seqids(idcnt) = reclin(1:30)
               line = line + 1
               read(66,120,end = 410,err=991) reclin
c     DO WHILE ((INDEX(RECLIN,'1').EQ.0).AND.
c     .                       (INDEX(RECLIN,'2').EQ.0))
c     ENDDO
 615           if (index(reclin,'1').eq.0) then
                  if (index(reclin,'2').eq.0) then
                     line = line + 1
                     read(66,120,end = 410,err=991) reclin
                     goto 615
                  endif
               endif
c     assume at least one line not a sequence identifier occurs
c     after last sequence read in to get around the CTRL-L problem.
c     therefore, after reading in the record containing the end
c     of sequence identifier, read in another record
               line = line + 1
               read(66,120,end = 410,err=991) reclin
            endif
c     ENDDO
            if (.not.endfil) goto 610
c     GENBANK format, recognized by the word LOCUS starting in the
c     first position of the record
         elseif (reclin(1:5).eq.'LOCUS') then
            found = .true.
            stype = 'GENBANK '
            idcnt = 1
            seqids(idcnt) = reclin(13:27)
c     DO WHILE (.NOT.ENDFIL)
 620        line = line + 1
            read(66,120,end = 410,err=991) reclin
c     scrolling through to find the key phrase ORIGIN because
c     the line after this key phrase occurrance is where the sequence
c     occurs in the file;
c     read through the file obtaining other sequences by scrolling
c     through to the '//' which signals the end of sequence
c     Then start looking for another sequence identifier.
            if (reclin(1:6).eq.'ORIGIN') then
               sline(idcnt) = line
c     DO WHILE (RECLIN(1:2).NE.'//')
 625           if (reclin(1:2).ne.'//') then
                  line = line + 1
                  read(66,120,end = 410,err=991) reclin
                  goto 625
c     ENDDO
               endif
            elseif (reclin(1:5).eq.'LOCUS') then
               idcnt = idcnt + 1
               seqids(idcnt) = reclin(13:27)
            endif
c     ENDDO
            if (.not.endfil) goto 620
c     PIR format, recognized by the '>' occurring in column 1 of a
c     record.  The sequence itself start 2 lines down from the record
c     containing '>'.
         elseif (reclin(1:1).eq.'>') then
            found = .true.
            stype = 'PIR     '
c     DO WHILE (.NOT.ENDFIL)
 630        if (reclin(1:1).eq.'>') then
               idcnt = idcnt + 1
               sline(idcnt) = line + 1
               seqids(idcnt) = reclin(5:25)
            endif
            line = line + 1
            read(66,120,end = 410,err=991) reclin
c     ENDDO
            if (.not.endfil) goto 630
c     EMBL (european) format, recognized by the key phrase ID in the
c     first column of the record.
         elseif (reclin(1:5).eq.'ID   ') then
            found = .true.
            stype = 'EMBL    '
c     scrolling through to find the key phrase 'SQ   Sequence' because
c     the line after this key phrase occurrance is where the sequence
c     occurs in the file;
c     DO WHILE (.NOT.ENDFIL)
 640        if (reclin(1:5).eq.'ID   ') then
               idcnt = idcnt + 1
               seqids(idcnt) = reclin(6:index(reclin,' '))
            elseif (reclin(1:13).eq.'SQ   Sequence') then
               sline(idcnt) = line
            endif
            line = line + 1
            read(66,120,end = 410,err=991) reclin
c     ENDDO
            if (.not.endfil) goto 640
c     NRC format, recognized by '(' occurring in the first column of
c     a record.  The record line following this line holds the number
c     of elements in the sequence and the name of the sequence. The
c     following line signals the beginning of the sequence.
c     The sequence itself is read in by using variable format
c     described by the '(' record line.  This format statement on the
c     record containing '(' ,must also be retrieved.
         elseif (reclin(1:1).eq.'(') then
            found = .true.
            stype = 'NRC     '
c     DO WHILE (.NOT.ENDFIL)
 650        if (reclin(1:1).eq.'(') then
               idcnt = idcnt + 1
               fmtseq(idcnt) = reclin(1:index(reclin,')'))
               line = line + 1
               read(66,140,end = 410,err = 991)
     .              seqnum(idcnt),seqids(idcnt)
               sline(idcnt) = line
            endif
            line = line + 1
            read(66,120,end = 410,err=991) reclin
c     ENDDO
 140        format(i4,5x,a30)
            if (.not.endfil) goto 650
         endif
c     keep scrolling through the file until a key phrase signalling
c     a format type is recognized or the end of file is found
         line = line + 1
c     ENDDO
         if (.not.found) goto 600
c     if a format type has not been found then this cannot be a
c     sequence file and return to main.
 410     if (.not.found) then
            write(6,105) ' No sequence identifiers found in this file '
 105        format(a50)
            return
         endif
c     if this is a valid format type proceed
         used = .true.
      endif
c     IF USED this call before start at listing identifiers
c     listing the sequence identifiers
c     allow the user to input the number requested
c     or control d to finish here and return to
c     entering a new filename
      valid = .false.
c     keep outputting the list of sequence identifiers
c     while NOT VALID choice of sequence identifier is inputted
c     DO WHILE (.NOT.VALID)
 660  write(6,150) 'Available sequences in ',filnam
 150  format(/,1x,a23,a50)
      do 525 k = 1,idcnt,2
         if (sline(k).ne.0) then
            if (k .eq. idcnt) then
               write (6,151) k,'.',seqids(k)
            else
               write (6,191) k,'.',seqids(k),k+1,'.',seqids(k+1)
            endif
 151        format(1x,i3,a1,1x,a30)
 191        format(1x,i3,a1,1x,a30,3x,i3,a1,1x,a30)
         endif
 525  continue
      choice = '                    '
      seqid  = '                    '
      pointr = 0
      write(6,152)
     .     'Choose sequence by number or name , or ? for relist;       '
      write(6,152)
     .     '      <return> defaults to the first one, / for new file.'
 152  format(1x,a60)
      read(5,153,end=10) choice
 153  format(a30)
c     error checking of the inputted CHOICE and determining
c     if this input is a  number, a default value, a relist command
c     or the name of the sequence identifier
      i = 1
      do while (choice(i:i).eq.' '.and.i.lt.6)
         i = i + 1
      enddo
      if ((choice(i:i).ge.'1').and.(choice(i:i).le.'9')) then
c     using COLLATING sequence to convert character to number value
         pointr = ichar(choice(i:i)) - ichar('0')
         if (choice(i+1:i+1).ne.' ') then
            pointr = ichar(choice(i+1:i+1)) - ichar('0') + 10 * pointr
         endif
         if (choice(i+2:i+2).ne.' ') then
            pointr = ichar(choice(i+2:i+2)) - ichar('0') + 10 * pointr
         endif
         if (pointr.gt.idcnt) then
            write(6,154) '  NUMBER CHOICE BETWEEN 1 AND ',idcnt
 154        format(a30,i2)
         else
            seqid   = seqids(pointr)
            valid = .true.
         endif
      elseif (ichar(choice(6:6)).eq.32) then
         pointr = 1
         seqid   = seqids(1)
         valid   = .true.
      elseif (choice(1:1).eq.'?') then
         continue
      else
c     find out if the inputted choice is in the list of seq identifiers.
         pointr = 0
         do 535 k = 1,idcnt
            if (seqids(k).eq.choice) then
               seqid   = seqids(k)
               pointr = k
               valid   = .true.
            endif
 535     continue
         if (pointr.eq.0) then
            write(6,105) 'Does not match any in given list'
         endif
      endif
c     ENDDO
      if (.not.valid) goto 660
c     having obtained a valid sequence choice or default of
c     one available sequence, rewind the file and retrieve the
c     sequence
      rewind(66,err = 992)
c     having retrieved the line number in the file where this
c     identifier occurs from SLINE(POINTR) scroll through the
c     file until this line is reached
c
c     SUN WARNING: brancd into block.
      do 550 i = 1,sline(pointr)
         read(66,120,end = 410,err=991) reclin
 550  continue
c     if nrc type then read sequence according to format type
      if (stype.eq.'NRC     ') then
         nseq = seqnum(pointr)
c     if the number in the NRC sequence is greater than
c     the maximum number sent down from main, then
c     truncate to NMAX and output a message to the user.
         if (nseq.gt.nmax) then
            write(6,160) '        Sequence truncated to ',nseq
 160        format(1x,a30,i5)
            nseq = nmax
         endif
         read(66,fmtseq(pointr),end = 420,err=991)
     .        (seq(k),k=1,seqnum(pointr))
c     Else if not type NRC then find the sequence by taking
c     each letter of the next records until end-of-sequence
c     indicator found or number in the sequence NMAX is reached
      else
         n = 1
         found = .false.
c     DO WHILE (.NOT.FOUND)
 670     read(66,120,end = 420,err=991) reclin
c     there is a need for the two positions to start reading the
c     sequence; to accommodate the GENBANK format and the STANFORD
c     end-of-sequence checks.
         if (stype.eq.'GENBANK ') then
            start = 10
         else
            start = 1
         endif
         do 565 i = start,80
c     if the number in the sequence is less than the desired
c     number, NMAX sent down from main, then retrieve sequence
            if (n.le.nmax) then
               if (.not.found) then
                  seq(n) = reclin(i:i)
c     check if an early end of sequence
                  if ((seq(n).eq.'1').or.(seq(n).eq.'2').or.
     .                 (reclin(1:1).eq.'/').or.(seq(n).eq.'*')) then
                     seq(n) = ' '
                     nseq = n - 1
                     found = .true.
c     if not an end-of-sequence character but have gone far
c     enough, then truncate to NMAX
                  elseif (n.eq.nmax) then
                     nseq = nmax
                     found = .true.
                     write(6,160) '        Sequence truncated to ',nseq
c     if not a end-of-sequence character, check to see if it is
c     not a blank character. Blank characters will not be added
c     to the sequence.
                  elseif (seq(n).ne.' ') then
                     n = n + 1
                  endif
               endif
c     if the number of sequence characters found is Greater than
c     NMAX then all of the sequence has been found.
            else
               found = .true.
               nseq = n
               write(6,160)
     .              '        Sequence truncated to ',nseq
            endif
 565     continue
c     ENDDO
         if (.not.found) goto 670
      endif
 420  return
c 991  stop '  ERROR IN READING FILE  '
c 992  stop '  ERROR IN REWINDING FILE'
c 999  stop '  END of SESSION...GOOD BYE'
 991  call exit(1)
 992  call exit(1)
 999  call exit(1)
      end
