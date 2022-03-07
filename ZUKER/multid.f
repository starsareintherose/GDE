       subroutine multid(seqid,seq,nseq,nmax,used,rnum)
c      REVISED VERSION OF FORMID -- TO WORK WITH MULTIPLE ALIGNMENTS
c      RESEARCHER: M. ZUKER
c      JUN 1986
c      DS TUDHOPE
c      WRITTEN TO WORK ON VAX-11-750 UNDER VMS4.3 USING STANDARD F-77
c      ####  NOTE #### USE OF CARRIAGE CONTROL LINE FEED SUPPRESSION
c                      IS USED IN THIS PROGRAM AND IS NOT STANDARD
c                      FORTRAN-77
c                      SEE '$' IN FORMAT STATEMENT 102
       implicit integer (a-z)
       logical found,endfil,valid,used
       integer sline(500),seqnum(500),nseq,nmax,rnum
       character*8  stype
       character*30 seqids(500)
c       character*30 choice,seqid
       character*30 seqid
       character*50 filnam,fmtseq(500)
       character*80 reclin
       character*1  seq(nmax)
       data   found/.false./,
     .        endfil/.false./,
     .         valid/.false./
c      **************************************************************
c      subroutine MULTID
c      **************************************************************
c      PURPOSE:
c      subroutine used to extract sequences from various
c      format type files including
c      STANFORD,GENBANK,EMBL,PIR,and,NRC
c      This revised version of FORMID is shorter.  Another
c      variable RNUM is passed first into the subroutine on
c      the first call to return the number of sequences in the
c      file.  The same variable is then used to request specific
c      sequences in do_loops of the main program.
c      This eliminates the task of having to select the sequences.
c      **************************************************************
c      variable list table:
c      * -- sent down from main returned unchanged
c      ** - returned to main from subroutine
c      ***- sent down from main and returned changed
c      INTEGERS:IDCNT   -- to keep track of the number of sequence
c                          identifiers found in the file
c               I,K     -- loop counters
c               LINE    -- to point to the record line number of a file
c               N       -- counter to extract the correct number of
c                          sequence elements
c           *   NMAX    -- maximum length of sequence expected by user
c          **   NSEQ    -- length of the sequence retrieved from file
c               POINTR  -- pointer to point to sequence identifier
c                          chosen from array of identifiers
c               RNUM    -- the number of identified sequences returned
c                          to main if first use of the subroutine,
c                          the number sent down from main, identifying
c                          a sequence if not the first time used
c               SEQNUM  -- an ARRAY of the length of the sequences
c                          for the sequence identifiers found in the
c                          NRC-format type files
c               SLINE   -- an ARRAY of the record-line numbers where
c                          a sequence starts in a file
c               START   -- defines what column of a record to start
c                          reading the sequence in.
c      CHARACTERS:CHOICE-- character string of length 20 to read in
c                          choice for sequence identifier to retrieve
c               FILNAM  -- character string of length 50 to read in
c                          the filename.
c               FMTSEQ  -- an ARRAY of characters, each length 50,
c                          describing how the sequence is to be read
c               RECLIN  -- a record of a file, length 80 characters
c            ** SEQ     -- an ARRAY of characters each length 1 to
c                          store the sequence
c            ** SEQID   -- retrieved name of the sequence
c               SEQIDS  -- an ARRAY of retrieved names of sequences
c               STYPE   -- character string of length 8 defining the
c                          format type of the file.
c      LOGICALS:FOUND   -- logical variable used in looping until
c                          something retrieved
c           *** USED    -- logical variable used in determining if the
c                          subroutine has been previously called
c               VALID   -- logical variable used in looping until
c                          some input valid
c
c
c      if the subroutine has NOT been USED then must
c      input the filename and do error checking on filename
c      else branch to listing the sequences available in this file
       if (.not.used) then
10     found = .false.
       valid = .false.
       used  = .false.
c      initialize for a new file by setting variables back
c      to zero and blanking out old names
       if (nmax.eq.0) nmax = 9999
       nseq = 0
       stype = '        '
       seqid = '                              '
       filnam= '                                                  '
       do 500 i = 1,100
          sline(i) = 0
          seqnum(i)= 0
          seqids(i)='                              '
       fmtseq(i)= '                                                  '
500    continue
111    write(6,102) 'Input sequence file name (/ to end)  '
102    format(1x,a38,$)
       read(5,110,end=999,err=111) filnam
       if (filnam(1:1).eq.'/') goto 999
110    format(a50)
c      open the file only after a valid filename has been retrieved
c      error in filename results in prompting for the input again
       open(66,file=filnam,status='OLD',err=10)
c      find sequence file format type and the sequence identifiers
       idcnt = 0
       line = 1
c      DO WHILE (.NOT.FOUND)
600       read(66,120,end = 410,err = 991) reclin
120       format(a80)
c  STANFORD format, recognized by the ';' in the first column
c  of the record
          if (reclin(1:1).eq.';') then
c            DO WHILE (.NOT.ENDFIL)
610            line = line + 1
               read(66,120,end = 410,err=991) reclin
c    to find the next sequence identifer scroll through the
c    file until the first character in the line is not ';'
               if (reclin(1:1).ne.';') then
                       found = .true.
                       stype ='STANFORD'
                       idcnt = idcnt + 1
                       sline(idcnt)= line
                   seqids(idcnt) = reclin(1:30)
                   line = line + 1
                   read(66,120,end = 410,err=991) reclin
c                   DO WHILE ((INDEX(RECLIN,'1').EQ.0).AND.
c     .                       (INDEX(RECLIN,'2').EQ.0))
c                   ENDDO
615                   if (index(reclin,'1').eq.0) then
                       if (index(reclin,'2').eq.0) then
                        line = line + 1
                        read(66,120,end = 410,err=991) reclin
                         goto 615
                       endif
                      endif
c   assume at least one line not a sequence identifier occurs
c   after last sequence read in to get around the CTRL-L problem.
c   therefore, after reading in the record containing the end
c   of sequence identifier, read in another record
                   line = line + 1
                   read(66,120,end = 410,err=991) reclin
             endif
c            ENDDO
             if (.not.endfil) goto 610
c   GENBANK format, recognized by the word LOCUS starting in the
c   first position of the record
          elseif (reclin(1:5).eq.'LOCUS') then
                 found = .true.
                 stype = 'GENBANK '
                 idcnt = 1
                 seqids(idcnt) = reclin(13:27)
c                DO WHILE (.NOT.ENDFIL)
620                  line = line + 1
                     read(66,120,end = 410,err=991) reclin
c   scrolling through to find the key phrase ORIGIN because
c   the line after this key phrase occurrance is where the sequence
c   occurs in the file;
c   read through the file obtaining other sequences by scrolling
c   through to the '//' which signals the end of sequence
c   Then start looking for another sequence identifier.
                    if (reclin(1:6).eq.'ORIGIN') then
                        sline(idcnt) = line
c                       DO WHILE (RECLIN(1:2).NE.'//')
625                     if (reclin(1:2).ne.'//') then
                           line = line + 1
                           read(66,120,end = 410,err=991) reclin
                           goto 625
c                       ENDDO
                        endif
                     elseif (reclin(1:5).eq.'LOCUS') then
                           idcnt = idcnt + 1
                           seqids(idcnt) = reclin(13:27)
                     endif
c                  ENDDO
                   if (.not.endfil) goto 620
c  PIR format, recognized by the '>' occurring in column 1 of a
c  record.  The sequence itself start 2 lines down from the record
c  containing '>'.
          elseif (reclin(1:1).eq.'>') then
                 found = .true.
                 stype = 'PIR     '
c                DO WHILE (.NOT.ENDFIL)
630                 if (reclin(1:1).eq.'>') then
                        idcnt = idcnt + 1
                        sline(idcnt) = line + 1
                        seqids(idcnt) = reclin(5:25)
                    endif
                    line = line + 1
                    read(66,120,end = 410,err=991) reclin
c                  ENDDO
                   if (.not.endfil) goto 630
c  EMBL (european) format, recognized by the key phrase ID in the
c  first column of the record.
          elseif (reclin(1:5).eq.'ID   ') then
                 found = .true.
                 stype = 'EMBL    '
c   scrolling through to find the key phrase 'SQ   Sequence' because
c   the line after this key phrase occurrance is where the sequence
c   occurs in the file;
c                 DO WHILE (.NOT.ENDFIL)
640                  if (reclin(1:5).eq.'ID   ') then
                           idcnt = idcnt + 1
                     seqids(idcnt) = reclin(6:index(reclin,' '))
                     elseif (reclin(1:13).eq.'SQ   Sequence') then
                           sline(idcnt) = line
                     endif
                     line = line + 1
                     read(66,120,end = 410,err=991) reclin
c                  ENDDO
                   if (.not.endfil) goto 640
c  NRC format, recognized by '(' occurring in the first column of
c  a record.  The record line following this line holds the number
c  of elements in the sequence and the name of the sequence. The
c  following line signals the beginning of the sequence.
c  The sequence itself is read in by using variable format
c  described by the '(' record line.  This format statement on the
c  record containing '(' ,must also be retrieved.
          elseif (reclin(1:1).eq.'(') then
                  found = .true.
                  stype = 'NRC     '
c                 DO WHILE (.NOT.ENDFIL)
650                  if (reclin(1:1).eq.'(') then
                        idcnt = idcnt + 1
                        fmtseq(idcnt) = reclin(1:index(reclin,')'))
                        line = line + 1
                        read(66,140,end = 410,err = 991)
     .                  seqnum(idcnt),seqids(idcnt)
                        sline(idcnt) = line
                     endif
                     line = line + 1
                     read(66,120,end = 410,err=991) reclin
c                  ENDDO
140              format(i4,5x,a30)
                 if (.not.endfil) goto 650
          endif
c      keep scrolling through the file until a key phrase signalling
c      a format type is recognized or the end of file is found
       line = line + 1
c      ENDDO
       if (.not.found) goto 600
c      if a format type has not been found then this cannot be a
c      sequence file and return to main.
410    if (.not.found) then
         write(6,105) ' No sequence identifiers found in this file '
105    format(a50)
           return
       endif
c      if this is a valid format type proceed
       used = .true.
c      if this is the first time used then return to main
c      with RNUM, the number of identified sequences in the file
       rnum = idcnt
       endif
c      IF USED this call before retreive a sequence using RNUM
       seqid  = '                    '
       pointr = rnum
       if ((pointr.lt.1).or.(pointr.gt.idcnt)) then
c             stop ' ERROR IN MULTID, SEQUENCE REQUESTED NOT FOUND'
          call exit(1)
       else
             seqid   = seqids(pointr)
       endif
c      having obtained a valid sequence choice or default of
c      one available sequence, rewind the file and retrieve the
c      sequence
       rewind(66,err = 992)
c      having retrieved the line number in the file where this
c      identifier occurs from SLINE(POINTR) scroll through the
c      file until this line is reached
            do 550 i = 1,sline(pointr)
              read(66,120,end = 410,err=991) reclin
550         continue
c      if nrc type then read sequence according to format type
       if (stype.eq.'NRC     ') then
              nseq = seqnum(pointr)
c      if the number in the NRC sequence is greater than
c      the maximum number sent down from main, then
c      truncate to NMAX and output a message to the user.
              if (nseq.gt.nmax) then
                 write(6,160) '        Sequence truncated to ',nseq
160              format(1x,a30,i4)
                 nseq = nmax
              endif
              read(66,fmtseq(pointr),end = 420,err=991)
     .            (seq(k),k=1,seqnum(pointr))
c      Else if not type NRC then find the sequence by taking
c      each letter of the next records until end-of-sequence
c      indicator found or number in the sequence NMAX is reached
       else
       n = 1
       found = .false.
c       DO WHILE (.NOT.FOUND)
670         read(66,120,end = 420,err=991) reclin
c      there is a need for the two positions to start reading the
c      sequence; to accommodate the GENBANK format and the STANFORD
c      end-of-sequence checks.
            if (stype.eq.'GENBANK ') then
                start = 10
            else
                start = 1
            endif
            do 565 i = start,80
c      if the number in the sequence is less than the desired
c      number, NMAX sent down from main, then retrieve sequence
               if (n.le.nmax) then
                  if (.not.found) then
                  seq(n) = reclin(i:i)
c       check if an early end of sequence
                    if ((seq(n).eq.'1').or.(seq(n).eq.'2').or.
     .              (reclin(1:1).eq.'/').or.(seq(n).eq.'*')) then
                        seq(n) = ' '
                        nseq = n - 1
                        found = .true.
c      if not an end-of-sequence character but have gone far
c      enough, then truncate to NMAX
                    elseif (n.eq.nmax) then
                        nseq = nmax
                        found = .true.
                 write(6,160) '        Sequence truncated to ',nseq
c      if not a end-of-sequence character, check to see if it is
c      not a blank character. Blank characters will not be added
c      to the sequence.
                    elseif (seq(n).ne.' ') then
                         n = n + 1
                       endif
                    endif
c      if the number of sequence characters found is Greater than
c      NMAX then all of the sequence has been found.
               else
                    found = .true.
                    nseq = n
                    write(6,160)
     .              '        Sequence truncated to ',nseq
            endif
565         continue
c          ENDDO
           if (.not.found) goto 670
          endif
420       return
c 991    stop '  ERROR IN READING FILE  '
c 992    stop '  ERROR IN REWINDING FILE'
c 999    stop '  END of SESSION...GOOD BYE'
991    call exit(1)
992    call exit(1)
999    call exit(1)
       end
