      Program scfEnergyTerms
!
!     This program reads a series of matrices from a Gaussian matrix file
!     generated during an SCF calculation and then computes 1-electron and
!     2-electron energy contributions. The program also reads atomic numbers and
!     Cartesian coordinates from the matrix file to evaluate the nuclear-nuclear
!     repulsion energy. Using these component values, this program ends by
!     reporting the total SCF energy.
!
!     -H. P. Hratchian, 2020.
!
!
!     USE Connections
!
      use mqc_general
      use mqc_molecule
      use mqc_gaussian
      use mqc_algebra2
      use iso_fortran_env
!
!     Variable Declarations
!
      implicit none
      integer,parameter::IOut=6
      integer(kind=int64)::nCommands,i,j,k1,k2,nAtoms,nAt3
      integer(kind=int64),dimension(:),allocatable::atomicNumbers
      real(kind=real64)::Vnn,Escf
      real(kind=real64),dimension(3)::tmp3Vec
      real(kind=real64),dimension(:),allocatable::cartesians
      real(kind=real64),dimension(:,:),allocatable::distanceMatrix
      character(len=512)::matrixFilename,tmpString
      type(mqc_gaussian_unformatted_matrix_file)::GMatrixFile
      type(MQC_Variable)::tmpMQCvar
      type(MQC_Variable)::nEalpha,nEbeta,nEtot,KEnergy,VEnergy,OneElEnergy,  &
        TwoElEnergy,scfEnergy
      type(MQC_Variable)::SMatrixAO,TMatrixAO,VMatrixAO,HCoreMatrixAO,  &
        FMatrixAlpha,FMatrixBeta,PMatrixAlpha,PMatrixBeta,PMatrixTotal
!
!     Format Statements
!
 1000 Format(1x,'Enter Test Program scfEnergyTerms.')
 1010 Format(3x,'Matrix File: ',A,/)
 1100 Format(1x,'nAtoms=',I4)
 1200 Format(1x,'Atomic Coordinates (Angstrom)')
 1210 Format(3x,I3,2x,A2,5x,F7.4,3x,F7.4,3x,F7.4)
 1300 Format(1x,'Nuclear Repulsion Energy = ',F20.6)
 8999 Format(/,1x,'END OF TEST PROGRAM scfEnergyTerms.')
!
!
      write(IOut,1000)
!
!     Open the Gaussian matrix file and load the number of atomic centers.

      nCommands = command_argument_count()
      if(nCommands.eq.0)  &
        call mqc_error('No command line arguments provided. The input Gaussian matrix file name is required.')
      call get_command_argument(1,matrixFilename)
      call GMatrixFile%load(matrixFilename)
      write(IOut,1010) TRIM(matrixFilename)
      nAtoms = GMatrixFile%getVal('nAtoms')
      write(IOut,1100) nAtoms
!
!     Figure out nAt3, then allocate memory for key arrays.
!
      nAt3 = 3*nAtoms
      Allocate(cartesians(NAt3),atomicNumbers(NAtoms))
!
!     Load up a few matrices from the matrix file.
!
      call GMatrixFile%getArray('OVERLAP',mqcVarOut=SMatrixAO)
      call GMatrixFile%getArray('KINETIC ENERGY',mqcVarOut=TMatrixAO)
      call GMatrixFile%getArray('CORE HAMILTONIAN ALPHA',mqcVarOut=HCoreMatrixAO)
      call GMatrixFile%getArray('ALPHA FOCK MATRIX',mqcVarOut=FMatrixAlpha)
      if(GMatrixFile%isUnrestricted()) then
        call GMatrixFile%getArray('BETA FOCK MATRIX',mqcVarOut=FMatrixBeta)
      else
        FMatrixBeta  = FMatrixAlpha
      endIf
      call GMatrixFile%getArray('ALPHA DENSITY MATRIX',mqcVarOut=PMatrixAlpha)
      if(GMatrixFile%isUnrestricted()) then
        call GMatrixFile%getArray('BETA DENSITY MATRIX',mqcVarOut=PMatrixBeta)
      else
        PMatrixBeta  = PMatrixAlpha
      endIf
      PMatrixTotal = PMatrixAlpha+PMatrixBeta
      VMatrixAO = HCoreMatrixAO-TMatrixAO
!
!     Calculate the number of electrons using <PS>.
!
      nEalpha = Contraction(PMatrixAlpha,SMatrixAO)
      nEbeta  = Contraction(PMatrixBeta,SMatrixAO)
      nEtot   = Contraction(PMatrixTotal,SMatrixAO)
      call nEalpha%print(IOut,' <P(Alpha)S>=')
      call nEbeta%print(IOut,' <P(Beta )S>=')
      call nEtot%print(IOut,' <P(Total)S>=')
!
!     Calculate the 1-electron energy and component pieces of the 1-electron
!     energy. Also, calculate the 2-electron energy.
!
      KEnergy     = Contraction(PMatrixTotal,TMatrixAO)
      VEnergy     = Contraction(PMatrixTotal,VMatrixAO)
      OneElEnergy = Contraction(PMatrixTotal,HCoreMatrixAO)
      TwoElEnergy = Contraction(PMatrixTotal,FMatrixAlpha)
      if(GMatrixFile%isUnrestricted()) then
        tmpMQCvar = Contraction(PMatrixTotal,FMatrixBeta)
        TwoElEnergy = TwoElEnergy + tmpMQCvar
      endIf
      TwoElEnergy = TwoElEnergy - OneElEnergy
      TwoElEnergy = MQC(0.5)*TwoElEnergy
      call KEnergy%print(IOut,' <P.K> = ')
      call VEnergy%print(IOut,' <P.V> =')
      call OneElEnergy%print(IOut,' <P.H> = ')
      call TwoElEnergy%print(IOut,' <P.F>-<P.H> = ')
!
!     Load the atommic numbers and Cartesian coordinates into our intrinsic
!     arrays.
!
      atomicNumbers = GMatrixFile%getAtomicNumbers()
      cartesians = GMatrixFile%getAtomCarts()
      cartesians = cartesians*angPBohr
!
!     Print out the atomic numbers and Cartesian coordiantes for each atomic
!     center.
!
      write(IOut,1200)
      do i = 1,NAtoms
        j = 3*(i-1)
        write(IOut,1210) i,mqc_element_symbol(atomicNumbers(i)),  &
          cartesians(j+1),cartesians(j+2),cartesians(j+3)
      endDo
!
!     Form the distance matrix between atomic centers.
!
      Allocate(distanceMatrix(nAtoms,nAtoms))
      do i = 1,nAtoms-1
        distanceMatrix(i,i) = float(0)
        k1 = 3*(i-1)+1
        do j = i+1,NAtoms
          k2 = 3*(j-1)+1
          tmp3Vec = cartesians(k1:k1+2)-cartesians(k2:k2+2)
          distanceMatrix(i,j) = sqrt(dot_product(tmp3Vec,tmp3Vec))
          distanceMatrix(j,i) = distanceMatrix(i,j)
        endDo
      endDo
!
!     Calculate the nuclear-nuclear repulsion energy.
!
      distanceMatrix = distanceMatrix/angPBohr
      Vnn = float(0)
      do i = 1,NAtoms-1
        do j = i+1,NAtoms
          Vnn = Vnn + float(atomicNumbers(i)*atomicNumbers(j))/distanceMatrix(i,j)
        endDo
      endDo
      write(iOut,1300) Vnn
!
!     Put things together and report the SCF energy.
!
      scfEnergy = oneElEnergy + twoElEnergy
      scfEnergy = scfEnergy + MQC(Vnn)
      call scfEnergy%print(IOut,' SCF Energy = ')
!
  999 Continue
      write(iOut,8999)
      end program scfEnergyTerms
