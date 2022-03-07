c      MFOLD - Prediction of RNA secondary structure by free energy
c              minimization.
c            - Version 2.0
c            - Michael Zuker and John Jaeger
c            - LRNA : folds linear RNA sequences
c            - CRNA : folds circular RNA sequences
c
c      The original version (1.0) was designed by Michael Zuker and
c      programmed by Eric Nelson in the summer of 1987 in the Division of
c      Biological Sciences at the National Research Council of Canada. John
c      Jaeger added the tetraloop bonus energy feature and created the
c      BATGEN program for batch file generation.
c
c      Version 2.0 corrects a number of small bugs from the original
c      program. These were added to version 1 and itemized in the
c      ERRATA.LIST file that was distributed along with version 1. The major
c      improvements of version 2 are : 1. During the generation of
c      suboptimal foldings, the number of new base pairs that are
c      sufficiently different from base pairs that have already been found
c      must be greater than the WINDOW parameter. This feature was added
c      during the summer of 1989, and was made part of version 1 (item 11 in
c      the ERRATA.LIST file distributed with this version). The effect is to
c      eliminate structures that contain just a few new base pairs. 2.
c      Temperature dependent folding. This was added in the fall of 1989 and
c      was never a feature of version 1.
c
c      METHOD : A dynamic programming algorithm is used to find optimal and
c      suboptimal foldings of an RNA molecule starting from linear sequence
c      data. Auxiliary information can be used to constrain the folding.
c
c      NB : Base pairs are forced by giving them a bonus energy (EPARAM(9)
c      in the program code). These energies are subtracted during the
c      traceback algorithm so that the computed sturctures have the correct
c      energies. Unfortunately, there is no way to subtract the bonus
c      energies from the energy dot plots. Moreover, each forced base pair
c      contain two bonus energies because of the nature of the algorithm.
c      For example, suppose that an optimal folding of an RNA contains 3
c      forced base pairs ( default bonus energy is 50.0 kcal per forced base
c      pair ) and that the correct folding energy is -180.0 kcal/mole.
c      Internally, the energy will be -180.0 - (3+1) x 50.0 = -380.0
c      kcal/mole. To find foldings within 10% of the correct energy, one
c      needs to compute foldings to within 18.0 kcal of -180.0 - 3 x 50.0 =
c      -330.0 kcal/mole. This comes out to -312.0 kcal/mole. The ratio of
c      -312.0 to -380.0 is 82%, so that one would request the 18% level of
c      suboptimality! This confustion only exists when base pairs are
c      forced. Each closed excision counts as one forced base pair.
c
c      Energy data from :
c      S.M. Freier et al., Proc. Natl. Acad. Sci. USA, 83, 9373-9377, 1986.
c      D.H. Turner et al., Cold Spring Harbor Symposia on Quantitative Biology,
c         52, 123-133, 1987.
c      D.H. Turner et al., Annu. Rev. Biophys. Biophys. Chem 17, 167-192 (1988).
c      This last reference has all the dangling end and terminal mismatch data.
c
c      References :
c      M. Zuker
c      On Finding All Suboptimal Foldings of an RNA Molecule.
c      Science, 244, 48-52, (1989)
c
c      J. A. Jaeger, D. H. Turner and M. Zuker
c      Improved Predictions of Secondary Structures for RNA.
c      Proc. Natl. Acad. Sci. USA, BIOCHEMISTRY, 86, 7706-7710, (1989)
c
c      J. A. Jaeger, D. H. Turner and M. Zuker
c      Predicting Optimal and Suboptimal Secondary Structure for RNA.
c      in "Molecular Evolution: Computer Analysis of Protein and
c      Nucleic Acid Sequences", R. F. Doolittle ed.
c      Methods in Enzymology, 183, 281-306 (1989)
c
      include 'rfd.inc'
      real energy
      logical flag,mark

c     Fill screen with author and reference data.
      call begin
c     Initial setup for run.
5     call device
      if (cntrl(1).ne.2) then
c        Read energy information if this is not a continuation run.
         call enefiles
         call erg(1,0,0,0,0)
      else
c        dump out information read in from continuation file
         call cdump
      endif
c     Determine output specifications if this is not a save run.
      if (cntrl(1).ne.1) call outputs
c     Call the menu if this is not a continuation run.
      if (cntrl(1).ne.2) call menu
      mrep = 1
c     CNTRL(7) = 0  -  suboptimal dot plot
c                1  -  N best sorted by energy
c                2  -  best folding for all sequences in a file
10    if (cntrl(7).eq.2) call mseq(mrep)
c     Process sequence before folding.
      do i = 1,mxbits
        marks(i) = 0
        force2(i) = 0
      enddo
      call process
      if (n*2.gt.fldmax) then
         tt = fldmax/2
c        Fragment is too long. Try again.
         write (6,*) 'Segment larger than ',tt
         stop
      endif
      if (cntrl(1).ne.2) then
c        Fill the optimal energy arrays except in a continuation run.
         call fill
      endif
 
      if (cntrl(1).eq.1) then
c        Save the results from FILL in a SAVE run and then stop.
         call putcont
         stop
      endif
 
 
      rep = 1
      jump = 1
      flag = .true.
      err = 0
 
      do while (flag)
         if (cntrl(7).eq.0) then
c           Interactive dot plot returns IRET, JRET (new numbering).
c Zuker comments out call to dotplt : do not choose this option
C If you do, the program stops dead here.
c           call dotplt(iret,jret,jump)
            if(1.eq.1) stop ' Energy dot plot disabled.'
            jump = 2
         else
c           Automatic sort returns IRET,JRET (new numbering).
            print *,'traceback'
            call sortout(iret,jret,rep,err)
            if (err.eq.30) then
                flag = .false.
                call errmsg(err,rep-1,0)
                err = 0
            endif
            rep = rep + 1
         endif
c        First traceback yields the best structure on the included fragment
c        from IRET to JRET.
         if (flag) call trace(iret,jret,nforc1,err)
         if (err.ne.0) call errmsg(err,iret,jret)
         if (flag) then
           it = iret+n
c        Second traceback yields the best structure on the excluded fragment
c        from IRET to JRET.
           call trace(jret,it,nforc2,err)
           if (err.ne.0) then
             call errmsg(err,jret,it)
           else
c            The energy of the best structure containing the base-pair IRET,
c            JRET is the sum of the energies of the optimal foldings on
c            the included and excluded fragments. A correction is made for
c            forced base-pairs.
c*           CALL EFN(ENE,1,N)
c*           WRITE (6,*) 'NEW ENERGY ',ENE
       ene = v(iret,jret) + v(jret,iret+n) - eparam(9) * (nforc1+nforc2)
             energy = float(ene) / 10.0
c            Count the number of new base pairs not within WINDOW
c            of existing base pairs.
             numbp = 0
             do k = 1,n
               if(k.lt.basepr(k)) then
                  if(.not.mark(k,basepr(k))) numbp = numbp + 1
               endif
             enddo
             do k = 1,n
               if(k.lt.basepr(k)) then
c                 Mark "traced-back" base pairs and also base-pairs
c                 which are close (within WINDOW = CNTRL(9) ).
                  call smark(k,basepr(k))
                  if(cntrl(9).gt.0) then
                     do k1 = -cntrl(9),cntrl(9)
                        do k2 = -cntrl(9),cntrl(9)
                          if(k+k1.gt.0.and.k+k1.lt.basepr(k)+k2.and.
     1                 basepr(k)+k2.le.n) call smark(k+k1,basepr(k)+k2)
                        enddo
                     enddo
                  endif
               endif
             enddo
             if(numbp.le.cntrl(9)) then
               rep = rep - 1
               go to 900
             endif
             write (6,1010) rep - 1
c  1010         format('+',i5$)
1010         format('+',i5)
             if (cntrl(2).ne.2.and.cntrl(2).ne.3.and.cntrl(2).ne.6) then
c               Line printer output.
                call linout(1,n,energy,iret,jret,err)
             endif
             if (err.ne.0) then
                call errmsg(err,iret,jret)
             else
                if (cntrl(2).ge.3.and.cntrl(2).ne.4) then
c                  Region table output.
                   call regtab
                endif
                if (mod(cntrl(2),2).eq.0.or.cntrl(2).eq.7) then
c                  CT file output.
                   call ct(energy)
                endif
             endif
           endif
         endif
         if (cntrl(7).eq.1.and.rep.gt.cntrl(6)) flag = .false.
900      continue
      enddo
c
c     Multiple sequence option (CNTRL(7) = 2)
c     If sequence number (MREP) is < total number of sequences
c     (CNTRL(5)), go get another sequence.
c
      if (cntrl(7).eq.2.and.mrep.lt.cntrl(5)) then
         mrep = mrep + 1
         goto 10
      endif
      stop
      end

c     Marks a base-pair I,J.
c     Assumes that 1 <= I <= J <= N.
c     The information is stored in a single bit in the MARKS
c     array.
c     The conversion from double dimension to single is through the
c     transformation  I,J  ==>  (J-1)*J/2 + I .
      subroutine smark(i,j)
      include 'rfd.inc'
      integer*2 bit
 
      posn = (((j-1)*j)/2) + i
      word = (posn+15) / 16
      bit = mod(posn,16)
c      marks(word) = iibset(marks(word),bit)
      marks(word) = ibset(marks(word),bit)
      return
      end
 
c     Marks a forced base-pair I,J.
c     The incoming base-pair II,JI is processed to an I,J
c     base-pair satisfying 1 <= I <= J <= N.
c     The information is stored in a single bit in the FORCE2
c     array.
c     The conversion from double dimension to single is through the
c     transformation  I,J  ==>  (J-1)*J/2 + I .
      subroutine sfce(ii,ji)
      include 'rfd.inc'
      integer*2 bit
 
      if (ii.gt.n) then
        i = ii - n
        j = ji - n
      elseif (ji.gt.n) then
        i = ji-n
        j = ii
      else
        i = ii
        j = ji
      endif
 
 
      posn = (((j-1)*j)/2) + i
      word = (posn+15) / 16
      bit = mod(posn,16)
 
c      force2(word) = iibset(force2(word),bit)
      force2(word) = ibset(force2(word),bit)
      return
      end
 
c     Retrieves information on whether or not the base-pair I,J
c     has been marked by a traceback passing through or close to
c     this pair.
      logical function mark(i,j)
      include 'rfd.inc'
      integer*2 bit,one
 
      one = 1
      posn = (((j-1)*j)/2) + i
      word = (posn+15) / 16
      bit = mod(posn,16)
 
c      set = iibits(marks(word),bit,one)
      set = ibits(marks(word),bit,one)
      mark = .false.
      if (set.ne.0) mark = .true.
      return
      end
 
 
c     Retrieves information on whether or not the base-pair I,J
c     has been forced.
      logical function fce(ii,ji)
      include 'rfd.inc'
      integer*2 bit,one
 
      if (ii.gt.n) then
        i = ii - n
        j = ji - n
      elseif (ji.gt.n) then
        i = ji-n
        j = ii
      else
        i = ii
        j = ji
      endif
 
      one = 1
      posn = (((j-1)*j)/2) + i
      word = (posn+15) / 16
      bit = mod(posn,16)
 
c      set = iibits(force2(word),bit,one)
      set = ibits(force2(word),bit,one)
      fce = .false.
      if (set.ne.0) fce = .true.
      return
      end
 
c     fills screen with author and reference information
      subroutine begin
      character*1 ans
      character*80 record
      open(3,file='begin.dat',status='old',err=5)
      write(6,1010) 
1010  format(' ')
1     read(3,1020,end=2) record
1020  format(a80)
      write(6,1030) record
c1030  format(' ',a80)
1030  format(a80)
      go to 1
2     write(6,1040) 
c     1040 format(' Press <return> to continue ...'$)
c 1040  format('  Press <return> to continue ...'$)
1040  format('  Press <return> to continue ...',$)
      read(5,1050) ans
1050  format(a1)
      return
5     write(6,1060)
1060  format(//' Author and reference file not available.'//)
c
c     C.Wang copied the next 2 lines here.  6/12/91
c     to keep consistency.
 3    write(6,1040) 
      read(5,1050) ans
      return
      end
