!> @file
!!
!!
!!   Copyright (C) 2019-2020
!!   @authors: Ioannis Deretzis, Giuseppe Fisicaro and Antonino La Magna
!!   This file is part of the MulSKIPS code.
!!   This file is distributed under the terms of the
!!   GNU General Public License, see ~/COPYING file
!!   or http://www.gnu.org/licenses/gpl-3.0.txt .
!!   For the list of contributors, see ~/AUTHORS

!!   MulSKIPS is a free software: you can redistribute it and/or modify
!!   it under the terms of the GNU General Public License as published by
!!   the Free Software Foundation, either version 3 of the License, or
!!   (at your option) any later version.

!!   MulSKIPS is distributed in the hope that it will be useful,
!!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!!   GNU General Public License for more details
!!********************************************************************
!!*     It picks the leaf of the current even in the tree and returns
!!*     index in the list and related event probability
!!********************************************************************
      SUBROUTINE PickTreeEvent(Index,Patom,Run)

!     Index - output
!     Patom - output
!     Levels,SizeTree,Tree - input, info su albero da usare
       USE Definitions
       USE DefSystem
       IMPLICIT NONE
       INTEGER :: Index
       REAL(8) :: Patom,Ptot
       CHARACTER(Len=1)  :: Run
       INTEGER Parent,i
       REAL(8) Random, currrand, numberslist
       EXTERNAL Random
       Parent = 1
       Index=1
       DO
	   IF (Parent.GT.Levels-1) EXIT
              IF(Run.EQ.'T')THEN
                  counter = counter + 1
                  currrand = numberslist(counter)
                  !write(51,'(I8,1x,e15.8)')counter,currrand
              ELSE
                  currrand = Random(idum)
                  !write(50,'(e15.8)')currrand
              END IF
	   Patom=currrand*Tree(Index)
	   Index = ISHFT(Index,1)
	   IF(Tree(Index).LT.Patom) Index=Index+1
	   Parent=Parent+1
       END DO
       Patom = Tree(Index)

       IF (Index-ISHFT(1,Levels-1)+1.GT.NumAdAtom)THEN
         Ptot=0
         DO i=ISHFT(1,Levels-1),ISHFT(1,Levels-1)+NumAdAtom-1
            Ptot=Ptot+Tree(i)
            write(*,*)i,i - ISHFT(1,Levels-1)+1,Tree(i)
         END DO
         write(*,*)'write Index> NumAdAtom',Ptot,Tree(1)
         STOP
       END IF

       Index = Index - ISHFT(1,Levels-1)+1

      END SUBROUTINE PickTreeEvent

      function numberslist(ind) result(currrandom)
       IMPLICIT NONE
       INTEGER, intent(in)  :: ind ! input
       REAL(8) :: currrandom
       INTEGER, PARAMETER :: rdim = 1000
       REAL(8), dimension(rdim) :: numbers = (/ 0.22188498E+00, 
     >  0.17518266E+00,
     >  0.40687769E-01,
     >  0.97799197E+00,
     >  0.99381448E+00,
     >  0.46786985E+00,
     >  0.17739029E-02,
     >  0.91279589E+00,
     >  0.77248372E+00,
     >  0.43659709E+00,
     >  0.90317065E+00,
     >  0.33233442E+00,
     >  0.70141988E+00,
     >  0.45583782E+00,
     >  0.34291628E+00,
     >  0.93128070E+00,
     >  0.51518693E+00,
     >  0.98950539E+00,
     >  0.69884517E+00,
     >  0.18424178E+00,
     >  0.46166453E+00,
     >  0.64453616E+00,
     >  0.46135350E+00,
     >  0.11669448E+00,
     >  0.14244262E+00,
     >  0.99587454E+00,
     >  0.19534932E+00,
     >  0.43198875E+00,
     >  0.90081907E+00,
     >  0.49299808E+00,
     >  0.93762256E+00,
     >  0.95671479E+00,
     >  0.14616436E+00,
     >  0.98751275E+00,
     >  0.70262590E+00,
     >  0.46481522E+00,
     >  0.17628363E+00,
     >  0.27303514E+00,
     >  0.25591099E+00,
     >  0.69921336E+00,
     >  0.20193184E+00,
     >  0.67403284E+00,
     >  0.58784899E+00,
     >  0.89822364E-01,
     >  0.88113567E+00,
     >  0.76797499E+00,
     >  0.13904146E+00,
     >  0.50090251E+00,
     >  0.65402841E+00,
     >  0.98326706E+00,
     >  0.62375480E-01,
     >  0.32208758E+00,
     >  0.67007289E+00,
     >  0.70181879E+00,
     >  0.87261965E+00,
     >  0.26517017E+00,
     >  0.29018298E-01,
     >  0.53174986E-01,
     >  0.27536607E+00,
     >  0.52899926E+00,
     >  0.29158622E+00,
     >  0.72873873E+00,
     >  0.65688490E+00,
     >  0.73270355E-01,
     >  0.23466525E+00,
     >  0.22913781E+00,
     >  0.74448540E+00,
     >  0.61159752E+00,
     >  0.57470212E+00,
     >  0.57494126E+00,
     >  0.79223925E+00,
     >  0.14284420E-01,
     >  0.33547698E+00,
     >  0.71557808E+00,
     >  0.12186630E+00,
     >  0.13957696E+00,
     >  0.97446324E+00,
     >  0.75953468E+00,
     >  0.24407480E+00,
     >  0.87727243E+00,
     >  0.96685624E+00,
     >  0.14217434E+00,
     >  0.15662269E+00,
     >  0.37181981E+00,
     >  0.20141185E+00,
     >  0.20888383E+00,
     >  0.29982989E+00,
     >  0.72894003E-01,
     >  0.75284750E+00,
     >  0.47348809E+00,
     >  0.72032979E+00,
     >  0.56468609E+00,
     >  0.69833299E+00,
     >  0.68096970E+00,
     >  0.90697409E+00,
     >  0.18764742E+00,
     >  0.33855586E+00,
     >  0.87227088E+00,
     >  0.96795603E+00,
     >  0.74155871E+00,
     >  0.79351172E+00,
     >  0.37950675E+00,
     >  0.25682771E+00,
     >  0.77675596E+00,
     >  0.16410818E-01,
     >  0.92020112E+00,
     >  0.16546489E+00,
     >  0.29825308E+00,
     >  0.50040693E+00,
     >  0.66373582E+00,
     >  0.96534025E+00,
     >  0.95612427E+00,
     >  0.30032746E+00,
     >  0.80187794E+00,
     >  0.80866944E+00,
     >  0.72690011E+00,
     >  0.30405739E-01,
     >  0.97591517E+00,
     >  0.16629624E+00,
     >  0.47017829E-01,
     >  0.89058192E+00,
     >  0.87221449E+00,
     >  0.64364146E+00,
     >  0.83314338E+00,
     >  0.78142951E+00,
     >  0.41273250E+00,
     >  0.75745668E+00,
     >  0.55872100E+00,
     >  0.69916726E+00,
     >  0.20166516E+00,
     >  0.97411204E+00,
     >  0.67621016E+00,
     >  0.25912775E+00,
     >  0.58033895E+00,
     >  0.91193215E+00,
     >  0.10731975E-01,
     >  0.84184685E+00,
     >  0.35474472E+00,
     >  0.56315034E+00,
     >  0.47451172E+00,
     >  0.17847809E+00,
     >  0.32391469E+00,
     >  0.90659774E+00,
     >  0.70582967E+00,
     >  0.58290615E+00,
     >  0.84811527E+00,
     >  0.92104460E+00,
     >  0.86518958E+00,
     >  0.89954016E+00,
     >  0.49424159E+00,
     >  0.43019071E+00,
     >  0.77983484E+00,
     >  0.17310362E+00,
     >  0.76629088E+00,
     >  0.76744665E+00,
     >  0.11730156E+00,
     >  0.12037900E+00,
     >  0.67648873E+00,
     >  0.86482378E+00,
     >  0.56788428E-02,
     >  0.78354264E-01,
     >  0.81072014E+00,
     >  0.73510271E+00,
     >  0.25895214E-01,
     >  0.48525773E+00,
     >  0.64142556E+00,
     >  0.49526532E-01,
     >  0.59449775E+00,
     >  0.21897180E+00,
     >  0.96055414E+00,
     >  0.80585547E+00,
     >  0.16521613E+00,
     >  0.76375008E-01,
     >  0.67205462E+00,
     >  0.61682709E+00,
     >  0.11074708E+00,
     >  0.69911087E+00,
     >  0.87735055E+00,
     >  0.65696739E-01,
     >  0.66412795E+00,
     >  0.29235350E+00,
     >  0.80967951E-01,
     >  0.69389718E+00,
     >  0.69348842E+00,
     >  0.12331089E+00,
     >  0.16339190E+00,
     >  0.94110742E+00,
     >  0.23323253E+00,
     >  0.95081221E-01,
     >  0.27050659E+00,
     >  0.96120541E+00,
     >  0.24734910E+00,
     >  0.13577292E+00,
     >  0.60259617E+00,
     >  0.66865622E+00,
     >  0.13261959E-01,
     >  0.24753968E+00,
     >  0.23454311E+00,
     >  0.89002587E-01,
     >  0.47215906E+00,
     >  0.14900440E+00,
     >  0.43694051E-01,
     >  0.79949284E+00,
     >  0.23541221E+00,
     >  0.20188809E+00,
     >  0.34922276E+00,
     >  0.85937654E-01,
     >  0.47961517E+00,
     >  0.64297998E+00,
     >  0.60405475E+00,
     >  0.17619411E+00,
     >  0.88714644E+00,
     >  0.58140751E+00,
     >  0.59431720E+00,
     >  0.44473400E-01,
     >  0.83100514E+00,
     >  0.67494722E+00,
     >  0.13250654E+00,
     >  0.35723897E+00,
     >  0.47199577E+00,
     >  0.39388587E+00,
     >  0.81498339E+00,
     >  0.50549517E+00,
     >  0.74681270E+00,
     >  0.81154975E+00,
     >  0.76216142E+00,
     >  0.36572326E+00,
     >  0.84096277E+00,
     >  0.47016654E+00,
     >  0.26760432E+00,
     >  0.24809428E-01,
     >  0.21949570E+00,
     >  0.23437057E+00,
     >  0.46164196E+00,
     >  0.48793384E+00,
     >  0.40520703E+00,
     >  0.49956041E+00,
     >  0.99579987E-01,
     >  0.64901502E+00,
     >  0.29230573E+00,
     >  0.48844465E+00,
     >  0.80860089E+00,
     >  0.87599354E+00,
     >  0.62308542E+00,
     >  0.87662068E+00,
     >  0.14622202E+00,
     >  0.74185391E+00,
     >  0.38896019E+00,
     >  0.79104639E+00,
     >  0.90649477E+00,
     >  0.64753867E+00,
     >  0.40657689E+00,
     >  0.76437655E+00,
     >  0.82139824E+00,
     >  0.44734964E+00,
     >  0.92950867E+00,
     >  0.80932345E+00,
     >  0.33785088E+00,
     >  0.74747834E+00,
     >  0.79668102E+00,
     >  0.84966232E+00,
     >  0.98635764E+00,
     >  0.83060013E+00,
     >  0.35067425E+00,
     >  0.11561010E+00,
     >  0.36759320E+00,
     >  0.11152903E-01,
     >  0.95832206E+00,
     >  0.71769648E+00,
     >  0.89825135E+00,
     >  0.89151231E-01,
     >  0.28598703E+00,
     >  0.34146011E+00,
     >  0.45074417E+00,
     >  0.82445708E+00,
     >  0.98730896E+00,
     >  0.50606847E-01,
     >  0.68409690E+00,
     >  0.29946307E+00,
     >  0.88204105E+00,
     >  0.95283794E+00,
     >  0.27872387E-01,
     >  0.93484429E-01,
     >  0.67348548E+00,
     >  0.41794197E+00,
     >  0.38451762E-01,
     >  0.38889554E+00,
     >  0.88369629E+00,
     >  0.34603186E+00,
     >  0.12034064E+00,
     >  0.39405413E+00,
     >  0.54123833E+00,
     >  0.38188348E+00,
     >  0.75076364E+00,
     >  0.20315450E+00,
     >  0.20245762E+00,
     >  0.46714077E+00,
     >  0.42524936E+00,
     >  0.79862831E+00,
     >  0.88931169E+00,
     >  0.95615176E-01,
     >  0.57757003E-01,
     >  0.89497121E-01,
     >  0.90900532E+00,
     >  0.95365679E+00,
     >  0.61966628E+00,
     >  0.31309246E+00,
     >  0.90891060E-01,
     >  0.40345626E+00,
     >  0.40889787E+00,
     >  0.54061313E+00,
     >  0.92562713E+00,
     >  0.99181899E+00,
     >  0.62713770E+00,
     >  0.40262689E+00,
     >  0.30842400E+00,
     >  0.60447416E+00,
     >  0.79836483E-01,
     >  0.14751975E+00,
     >  0.91315246E+00,
     >  0.90045240E+00,
     >  0.58590351E+00,
     >  0.15969375E+00,
     >  0.82838476E+00,
     >  0.80263617E+00,
     >  0.31394228E-01,
     >  0.19648991E+00,
     >  0.43245477E+00,
     >  0.49708735E+00,
     >  0.20479080E+00,
     >  0.67421650E+00,
     >  0.95971576E+00,
     >  0.28064064E+00,
     >  0.89056516E+00,
     >  0.34142792E+00,
     >  0.27210813E-01,
     >  0.36053370E-01,
     >  0.46634670E+00,
     >  0.27085859E+00,
     >  0.10951798E+00,
     >  0.43397757E+00,
     >  0.30905906E+00,
     >  0.73617654E+00,
     >  0.43287937E+00,
     >  0.21988822E+00,
     >  0.80815059E+00,
     >  0.38154458E+00,
     >  0.55349869E+00,
     >  0.94812744E+00,
     >  0.17176027E+00,
     >  0.59677048E-02,
     >  0.34686002E-01,
     >  0.92816198E+00,
     >  0.59383751E+00,
     >  0.21509519E+00,
     >  0.13589939E+00,
     >  0.77711634E+00,
     >  0.19893193E+00,
     >  0.56757739E+00,
     >  0.92644598E+00,
     >  0.58361291E+00,
     >  0.84674573E+00,
     >  0.82003244E+00,
     >  0.29393828E+00,
     >  0.97492028E+00,
     >  0.23155407E+00,
     >  0.18945059E+00,
     >  0.55893962E+00,
     >  0.40724948E+00,
     >  0.59447627E+00,
     >  0.92687939E+00,
     >  0.50975477E-01,
     >  0.13170901E+00,
     >  0.97575945E+00,
     >  0.90718475E+00,
     >  0.86576640E+00,
     >  0.65774150E+00,
     >  0.56585620E+00,
     >  0.61328957E+00,
     >  0.66673678E+00,
     >  0.25427786E+00,
     >  0.99755795E+00,
     >  0.86487735E+00,
     >  0.57064134E+00,
     >  0.62117786E+00,
     >  0.82747075E+00,
     >  0.13968332E+00,
     >  0.98670233E+00,
     >  0.91564486E+00,
     >  0.10987386E+00,
     >  0.83776019E+00,
     >  0.47711373E+00,
     >  0.59097220E-01,
     >  0.67638229E+00,
     >  0.18263856E+00,
     >  0.38300209E+00,
     >  0.17735005E+00,
     >  0.76041706E+00,
     >  0.52569459E+00,
     >  0.35412179E+00,
     >  0.15040909E+00,
     >  0.81568835E+00,
     >  0.94020909E+00,
     >  0.28139066E+00,
     >  0.91748238E+00,
     >  0.84097228E-02,
     >  0.16980863E+00,
     >  0.35752064E+00,
     >  0.97265962E+00,
     >  0.38762441E+00,
     >  0.99621604E+00,
     >  0.79041398E+00,
     >  0.28328704E+00,
     >  0.45770354E+00,
     >  0.88685787E-01,
     >  0.10649918E+00,
     >  0.78764850E+00,
     >  0.14365015E+00,
     >  0.11129972E+00,
     >  0.59191818E+00,
     >  0.54204016E-01,
     >  0.42903351E+00,
     >  0.33245023E-01,
     >  0.53127691E-01,
     >  0.44406718E+00,
     >  0.11119104E+00,
     >  0.11076636E+00,
     >  0.85031832E+00,
     >  0.58277077E-01,
     >  0.89877503E+00,
     >  0.69595777E+00,
     >  0.30022086E+00,
     >  0.59319656E+00,
     >  0.22566516E+00,
     >  0.67052071E+00,
     >  0.46386386E+00,
     >  0.71427091E+00,
     >  0.40717381E+00,
     >  0.48195555E+00,
     >  0.51467868E+00,
     >  0.39822247E-01,
     >  0.99603314E+00,
     >  0.87540261E+00,
     >  0.32372668E+00,
     >  0.55669840E-01,
     >  0.40872669E+00,
     >  0.44386870E+00,
     >  0.59695288E-02,
     >  0.23231511E+00,
     >  0.71447521E-01,
     >  0.27223573E+00,
     >  0.32703170E+00,
     >  0.70213998E+00,
     >  0.62691953E+00,
     >  0.65816399E+00,
     >  0.85018821E+00,
     >  0.22249179E+00,
     >  0.71454393E+00,
     >  0.61086992E+00,
     >  0.45361852E+00,
     >  0.29413878E+00,
     >  0.76263479E+00,
     >  0.87556506E+00,
     >  0.45798094E+00,
     >  0.34780216E+00,
     >  0.18290799E-03,
     >  0.91501134E+00,
     >  0.95956033E+00,
     >  0.40203370E+00,
     >  0.67995907E+00,
     >  0.66263045E+00,
     >  0.78167898E+00,
     >  0.91133502E+00,
     >  0.39852199E-01,
     >  0.31968245E+00,
     >  0.72717229E+00,
     >  0.72689350E+00,
     >  0.40632546E+00,
     >  0.39496367E+00,
     >  0.59387895E+00,
     >  0.88869922E+00,
     >  0.39622241E+00,
     >  0.23944840E+00,
     >  0.60465853E+00,
     >  0.60463625E+00,
     >  0.93332295E+00,
     >  0.42465577E+00,
     >  0.13521562E+00,
     >  0.87786297E+00,
     >  0.67033780E+00,
     >  0.54885249E+00,
     >  0.75471056E+00,
     >  0.51401099E-02,
     >  0.80199645E+00,
     >  0.85204820E+00,
     >  0.25814324E+00,
     >  0.84698120E-01,
     >  0.83555041E+00,
     >  0.40442319E-02,
     >  0.32849752E+00,
     >  0.68183316E+00,
     >  0.37543242E-01,
     >  0.61100583E+00,
     >  0.63843613E+00,
     >  0.18274827E+00,
     >  0.87601330E+00,
     >  0.87583296E-01,
     >  0.97481448E-01,
     >  0.22283282E-01,
     >  0.72484102E+00,
     >  0.42553244E+00,
     >  0.87276176E-01,
     >  0.83668093E+00,
     >  0.94053209E+00,
     >  0.90476600E+00,
     >  0.53942820E+00,
     >  0.75749468E+00,
     >  0.73568606E-01,
     >  0.60593271E+00,
     >  0.89658920E-01,
     >  0.91548476E+00,
     >  0.79460930E-01,
     >  0.95551610E+00,
     >  0.73536177E-01,
     >  0.99812589E+00,
     >  0.62508721E+00,
     >  0.17067314E+00,
     >  0.27289889E+00,
     >  0.85710390E+00,
     >  0.44366912E+00,
     >  0.63958900E+00,
     >  0.62941205E+00,
     >  0.38404218E+00,
     >  0.67012262E+00,
     >  0.16834651E+00,
     >  0.80142305E+00,
     >  0.55954145E+00,
     >  0.29891628E+00,
     >  0.69989250E+00,
     >  0.65208055E-01,
     >  0.17582827E+00,
     >  0.35108716E+00,
     >  0.52928288E+00,
     >  0.78820405E+00,
     >  0.75485301E+00,
     >  0.46939156E+00,
     >  0.79919443E+00,
     >  0.93160390E+00,
     >  0.80387054E+00,
     >  0.22696099E+00,
     >  0.87470099E-01,
     >  0.81179920E+00,
     >  0.97844648E+00,
     >  0.56037508E+00,
     >  0.68890850E+00,
     >  0.52421108E-01,
     >  0.65350103E+00,
     >  0.94088318E+00,
     >  0.47008962E+00,
     >  0.38132519E+00,
     >  0.31647185E+00,
     >  0.78866699E+00,
     >  0.39758892E+00,
     >  0.95707520E+00,
     >  0.54901275E+00,
     >  0.74445272E-01,
     >  0.55799327E+00,
     >  0.48476881E-01,
     >  0.18567908E+00,
     >  0.43537444E+00,
     >  0.74023374E+00,
     >  0.82589075E+00,
     >  0.26969804E+00,
     >  0.37897171E+00,
     >  0.21888219E-02,
     >  0.10368556E+00,
     >  0.10101442E+00,
     >  0.39514102E+00,
     >  0.38462765E+00,
     >  0.94570478E+00,
     >  0.97158615E+00,
     >  0.22978993E+00,
     >  0.80280924E+00,
     >  0.47577871E+00,
     >  0.12719727E+00,
     >  0.85092198E+00,
     >  0.23182313E+00,
     >  0.42696695E+00,
     >  0.12110987E+00,
     >  0.11323720E+00,
     >  0.26451806E+00,
     >  0.32497429E+00,
     >  0.34993749E+00,
     >  0.81389126E-01,
     >  0.15031116E+00,
     >  0.78601523E+00,
     >  0.65116746E+00,
     >  0.36837714E+00,
     >  0.40405341E+00,
     >  0.54697625E+00,
     >  0.85816573E+00,
     >  0.25537482E+00,
     >  0.85768014E+00,
     >  0.89899627E-02,
     >  0.50266777E+00,
     >  0.43317781E+00,
     >  0.83798649E+00,
     >  0.82059794E+00,
     >  0.81977331E+00,
     >  0.37618838E+00,
     >  0.13789542E+00,
     >  0.80540726E+00,
     >  0.67542979E+00,
     >  0.13307085E+00,
     >  0.63210091E+00,
     >  0.19907527E+00,
     >  0.99305612E+00,
     >  0.40768211E+00,
     >  0.26246163E+00,
     >  0.53451160E+00,
     >  0.66997298E-01,
     >  0.33618033E+00,
     >  0.27891450E+00,
     >  0.41153228E+00,
     >  0.12359689E+00,
     >  0.14450866E+00,
     >  0.94695593E-01,
     >  0.54664114E+00,
     >  0.12510683E+00,
     >  0.74505206E+00,
     >  0.41001659E+00,
     >  0.42662086E+00,
     >  0.33788329E+00,
     >  0.32178999E+00,
     >  0.17549219E+00,
     >  0.98752282E-01,
     >  0.79486601E+00,
     >  0.92203457E+00,
     >  0.10084509E+00,
     >  0.83574764E+00,
     >  0.24860294E+00,
     >  0.57872557E+00,
     >  0.19752077E+00,
     >  0.98879393E+00,
     >  0.71022990E-01,
     >  0.66985682E+00,
     >  0.55647186E+00,
     >  0.77003050E+00,
     >  0.44209021E+00,
     >  0.33511699E-03,
     >  0.73305890E+00,
     >  0.51032273E+00,
     >  0.44766354E+00,
     >  0.58236907E+00,
     >  0.16478448E+00,
     >  0.11138783E+00,
     >  0.66249430E+00,
     >  0.72184566E+00,
     >  0.43166804E+00,
     >  0.89773871E+00,
     >  0.27534328E+00,
     >  0.30214774E+00,
     >  0.55680431E+00,
     >  0.96704214E-01,
     >  0.93555006E+00,
     >  0.12805228E+00,
     >  0.32319930E+00,
     >  0.38096784E+00,
     >  0.62095503E+00,
     >  0.97803970E+00,
     >  0.29696677E+00,
     >  0.89409009E+00,
     >  0.27857938E+00,
     >  0.67847335E+00,
     >  0.61327414E+00,
     >  0.69684509E+00,
     >  0.51232649E+00,
     >  0.43356213E+00,
     >  0.85057535E+00,
     >  0.88414681E+00,
     >  0.40326114E+00,
     >  0.31338402E+00,
     >  0.51227785E+00,
     >  0.15127758E+00,
     >  0.78787977E-01,
     >  0.16320220E+00,
     >  0.15155906E+00,
     >  0.79398230E+00,
     >  0.77764577E+00,
     >  0.45477981E+00,
     >  0.62764789E+00,
     >  0.60068584E+00,
     >  0.90055396E+00,
     >  0.94703845E-01,
     >  0.79244358E+00,
     >  0.99138344E+00,
     >  0.41344011E+00,
     >  0.94466145E+00,
     >  0.44145369E-01,
     >  0.33646836E+00,
     >  0.59151484E+00,
     >  0.11618828E+00,
     >  0.32979776E+00,
     >  0.12904894E+00,
     >  0.34640216E+00,
     >  0.58370632E+00,
     >  0.55864347E+00,
     >  0.28010898E+00,
     >  0.10375641E+00,
     >  0.49769749E+00,
     >  0.84736791E+00,
     >  0.92915640E+00,
     >  0.49601834E+00,
     >  0.34996094E-01,
     >  0.54860310E+00,
     >  0.33560867E+00,
     >  0.33181583E+00,
     >  0.96752770E+00,
     >  0.67629355E+00,
     >  0.93389434E+00,
     >  0.96049838E+00,
     >  0.30257525E+00,
     >  0.76145940E+00,
     >  0.81235005E-01,
     >  0.30451320E+00,
     >  0.50417319E+00,
     >  0.30044049E+00,
     >  0.84461765E+00,
     >  0.33275042E-01,
     >  0.40852144E+00,
     >  0.65358007E+00,
     >  0.18836761E+00,
     >  0.83582922E+00,
     >  0.58276960E+00,
     >  0.12820610E+00,
     >  0.60295593E+00,
     >  0.45837363E+00,
     >  0.44582994E+00,
     >  0.48725208E+00,
     >  0.95135431E+00,
     >  0.66679148E+00,
     >  0.64270785E+00,
     >  0.99710465E+00,
     >  0.18320205E+00,
     >  0.96291034E+00,
     >  0.31955167E-01,
     >  0.87341646E-01,
     >  0.81574776E+00,
     >  0.48518008E+00,
     >  0.16366367E+00,
     >  0.52686422E+00,
     >  0.77751140E+00,
     >  0.94068130E+00,
     >  0.51057291E+00,
     >  0.93671597E-03,
     >  0.43043736E+00,
     >  0.67715302E+00,
     >  0.64538275E+00,
     >  0.51867547E-01,
     >  0.36011583E+00,
     >  0.82922684E+00,
     >  0.94940514E-01,
     >  0.75647451E+00,
     >  0.70555616E+00,
     >  0.68910795E+00,
     >  0.97042302E+00,
     >  0.49309149E+00,
     >  0.97098397E+00,
     >  0.92854322E+00,
     >  0.21523361E+00,
     >  0.34664331E+00,
     >  0.86349548E+00,
     >  0.25267176E+00,
     >  0.23459518E+00,
     >  0.30372358E+00,
     >  0.36383187E+00,
     >  0.99360025E+00,
     >  0.29950377E+00,
     >  0.41418028E+00,
     >  0.60171252E+00,
     >  0.82825176E+00,
     >  0.85802713E+00,
     >  0.75354274E+00,
     >  0.33265588E-01,
     >  0.84648139E+00,
     >  0.75281744E+00,
     >  0.75672196E+00,
     >  0.51682903E+00,
     >  0.45826282E+00,
     >  0.69580748E+00,
     >  0.11512335E-01,
     >  0.57689496E+00,
     >  0.28340917E+00,
     >  0.77921234E+00,
     >  0.74443289E+00,
     >  0.94860684E+00,
     >  0.65918676E+00,
     >  0.66812327E+00,
     >  0.51624399E+00,
     >  0.70999796E-01,
     >  0.80754165E+00,
     >  0.76372553E+00,
     >  0.17579888E+00,
     >  0.11242954E+00,
     >  0.65254576E+00,
     >  0.24739395E+00,
     >  0.39717178E+00,
     >  0.83067160E+00,
     >  0.89256528E+00,
     >  0.29514556E+00,
     >  0.84328677E+00,
     >  0.51953924E+00,
     >  0.13341936E+00,
     >  0.83428180E-01,
     >  0.17957955E+00,
     >  0.42214699E+00,
     >  0.90989559E+00,
     >  0.31179722E+00,
     >  0.26041995E+00,
     >  0.12149224E+00,
     >  0.83039930E+00,
     >  0.79249569E+00,
     >  0.44513009E+00,
     >  0.47086962E+00,
     >  0.12792470E+00,
     >  0.25140233E+00,
     >  0.34105450E+00,
     >  0.52109825E-01,
     >  0.17008509E-01,
     >  0.52545036E+00,
     >  0.87057335E+00,
     >  0.30656696E+00,
     >  0.98496496E+00,
     >  0.33848788E+00,
     >  0.62012338E+00,
     >  0.94983738E+00,
     >  0.33067045E+00,
     >  0.84682635E+00,
     >  0.29083890E+00,
     >  0.91377817E+00,
     >  0.38401026E+00,
     >  0.75109236E+00,
     >  0.45540272E+00,
     >  0.45300984E+00,
     >  0.98671662E+00,
     >  0.29930280E+00,
     >  0.47773722E+00,
     >  0.53126205E+00,
     >  0.41672094E+00,
     >  0.75268684E+00,
     >  0.46413416E+00,
     >  0.53991287E-01,
     >  0.28209128E+00,
     >  0.89315215E+00,
     >  0.86923189E+00,
     >  0.62727054E+00,
     >  0.44733437E+00,
     >  0.16376976E+00,
     >  0.56189483E+00,
     >  0.44831918E+00,
     >  0.55244787E+00,
     >  0.60576105E+00,
     >  0.74940907E+00,
     >  0.33233579E+00,
     >  0.72417680E+00,
     >  0.96913712E+00,
     >  0.92317894E+00,
     >  0.92668730E+00,
     >  0.66747399E-01,
     >  0.78053514E+00,
     >  0.84369899E+00,
     >  0.36880537E+00,
     >  0.36626513E+00,
     >  0.73850440E+00,
     >  0.12393779E+00,
     >  0.26996623E-01,
     >  0.42483925E+00,
     >  0.56967411E+00,
     >  0.36168060E+00,
     >  0.30867852E+00,
     >  0.85824775E+00,
     >  0.43251709E+00,
     >  0.73272681E+00,
     >  0.87071428E+00,
     >  0.61750159E+00,
     >  0.94272501E+00,
     >  0.36153331E+00,
     >  0.92364738E+00,
     >  0.36415157E+00,
     >  0.84703077E+00,
     >  0.60347510E+00,
     >  0.90739334E+00,
     >  0.86597348E-01,
     >  0.24821222E+00,
     >  0.13626400E+00,
     >  0.90001975E+00,
     >  0.27256926E+00,
     >  0.29278315E+00,
     >  0.72569022E+00,
     >  0.39294908E-01,
     >  0.48431715E+00,
     >  0.92041065E+00,
     >  0.58447363E+00,
     >  0.10984134E-01,
     >  0.69494742E+00,
     >  0.58133104E+00,
     >  0.75655623E+00,
     >  0.82983275E+00,
     >  0.22104472E+00,
     >  0.20036152E+00,
     >  0.52467177E+00,
     >  0.18829630E+00,
     >  0.42494242E+00,
     >  0.63757946E+00,
     >  0.88239241E+00,
     >  0.67496672E+00,
     >  0.79042330E+00,
     >  0.16672762E+00,
     >  0.50796587E+00,
     >  0.55091584E+00,
     >  0.64311513E+00,
     >  0.32697023E+00,
     >  0.25418725E+00,
     >  0.24262813E+00,
     >  0.99324379E+00,
     >  0.24770865E+00,
     >  0.42899033E+00,
     >  0.44566556E+00,
     >  0.66828300E+00,
     >  0.73984133E+00,
     >  0.14063588E+00,
     >  0.33357598E+00,
     >  0.24422080E+00,
     >  0.97399653E+00,
     >  0.72478031E+00,
     >  0.19255917E+00,
     >  0.30514555E+00,
     >  0.47914087E+00,
     >  0.24868066E+00,
     >  0.57372824E+00,
     >  0.68030315E+00,
     >  0.95509225E-01,
     >  0.35647751E+00,
     >  0.44348219E+00,
     >  0.75977446E+00,
     >  0.99402493E+00,
     >  0.89363584E+00,
     >  0.90677594E+00,
     >  0.24860608E-01,
     >  0.86379279E+00,
     >  0.74447579E+00,
     >  0.77977477E+00,
     >  0.47615662E+00,
     >  0.67740813E+00,
     >  0.45072663E+00,
     >  0.60733448E+00,
     >  0.31775919E-01,
     >  0.63727359E+00,
     >  0.91589914E+00,
     >  0.72122062E+00,
     >  0.27599111E+00,
     >  0.61456803E+00,
     >  0.78427102E-01,
     >  0.50424746E-01,
     >  0.68464913E-01,
     >  0.19409726E+00,
     >  0.12261795E+00,
     >  0.68094176E+00,
     >  0.89678743E+00,
     >  0.36309046E+00,
     >  0.95595831E+00,
     >  0.50971144E+00,
     >  0.46285332E+00,
     >  0.51708717E+00,
     >  0.57030050E+00,
     >  0.97826368E+00,
     >  0.83833105E+00,
     >  0.63650708E+00,
     >  0.10256775E+00,
     >  0.22473670E+00,
     >  0.38709636E+00,
     >  0.57584865E-01,
     >  0.62965274E+00,
     >  0.89556943E+00,
     >  0.67435557E+00,
     >  0.12409425E+00,
     >  0.11104829E+00,
     >  0.35652292E+00,
     >  0.42035149E+00,
     >  0.61240393E+00,
     >  0.66935446E+00,
     >  0.80391729E-01,
     >  0.80381611E+00,
     >  0.48431349E+00,
     >  0.43078252E+00,
     >  0.38968877E+00,
     >  0.45456008E+00,
     >  0.88552908E+00,
     >  0.44169359E+00,
     >  0.73450478E+00 /)

        currrandom = numbers(mod(ind-1,rdim)+1) 

      end function numberslist

