

! Declaration of all variables used throughout the text
MODULE GlobalVariables
  IMPLICIT NONE
  INTEGER :: ISEED
!  real(kind=16), dimension(1001,1001) :: SmallWorld1d
END MODULE GlobalVariables


! Set of all secondary routines used in the program
MODULE SecondaryRoutines
  USE GlobalVariables
  IMPLICIT NONE
  CONTAINS
  ! Generator of random numbers. It is called with a root given in the main program section (Iseed =1313) to have a different string of random nubers, change this root 
  FUNCTION RandomGenerator(Seed)
    INTEGER, INTENT(INOUT) :: Seed
    REAL(KIND=16) :: RandomGenerator
    ISEED = mod(8121*seed+28411, 134456)
    RandomGenerator = real(ISEED)/134456.
    RETURN
  END FUNCTION RandomGenerator
  ! The coupling function used in the text f(X_t^j). 
  FUNCTION CouplingFunction(x)
    REAL(KIND=16), INTENT(IN) :: x
    REAL(KIND=16) :: n
    REAL(KIND=16) :: CouplingFunction
    

    !The default is Linear
    n = 1
    CouplingFunction = x
    ! If one wishes for n-th potency, comment the lines aboves, specify n and uncoment the lines below
    
    !n=3
    !CouplingFunction = x**n
    
    !If one wishes for exponential or more exotic probabilities, ignore n, and uncoment the line below
    
    !CouplingFunction = exp(x)
    RETURN
  END FUNCTION CouplingFunction
  ! This is the GNUplot script routine that acompanies the routine TorusImages
  SUBROUTINE GnuScriptTorusImages
    INTEGER :: i
    CHARACTER(LEN=100) :: Arq1
    write(Arq1,"(A19)")"TorusImagesFilm.dat"
    open(1,FILE=TRIM(Arq1))
    write(1,*)"reset"
    write(1,*)"set terminal qt"
    write(1,*)"set xtics 1,200"
    write(1,*)"set ytics 1,200"
    write(1,*)"unset key"
    write(1,*)"set xrange [1:200]"
    write(1,*)"set yrange [1:200]"
    write(1,*)"set cbrange [0:1]"
    write(1,*)"set pm3d map"
    do i = 1,200
      write(1,*)"set terminal png"
      write(1,"(A22,I3.3,A5)")"set output 'TorusImage",i,".png'"
      write(1,"(A17,I3.3,A21)")"splot 'TorusImage",i,".dat' u 1:2:3 w image"
    end do
    write(1,*)"reset"
    write(1,*)"set terminal qt"
    close(1)
    RETURN
  END SUBROUTINE GnuScriptTorusImages
  ! This is the GNUplot script routine that acompanies the routine TorusImages
  SUBROUTINE GnuScriptTorusTriangularCellImages
    INTEGER :: i
    CHARACTER(LEN=100) :: Arq1
    write(Arq1,"(A28)")"TorusImagesTriangularCellFilm.dat"
    open(1,FILE=TRIM(Arq1))
    write(1,*)"reset"
    write(1,*)"set terminal qt"
    write(1,*)"set xtics 1,200"
    write(1,*)"set ytics 1,200"
    write(1,*)"unset key"
    write(1,*)"set xrange [1:200]"
    write(1,*)"set yrange [1:200]"
    write(1,*)"set cbrange [0:1]"
    write(1,*)"set pm3d map"
    do i = 1,200
      write(1,*)"set terminal png"
      write(1,"(A36,I3.3,A5)")"set output 'TorusTriangularCellImage",i,".png'"
      write(1,"(A31,I3.3,A21)")"splot 'TorusTriangularCellImage",i,".dat' u 1:2:3 w image"
    end do
    write(1,*)"reset"
    write(1,*)"set terminal qt"
    close(1)
    RETURN
  END SUBROUTINE GnuScriptTorusTriangularCellImages
  
END MODULE SecondaryRoutines


! Set of all secondary routines used in the program
MODULE MainRoutines
  USE GlobalVariables
  USE SecondaryRoutines
  IMPLICIT NONE
  CONTAINS
  !The main dificulty in implementing the different geometries is in the boundary conditions so I have separated the implementations in geometric groups and not number of neighbors
  !We begin with ring geometries. This routine produces a diagra saptio-temporal of the evolution of a ring geometry that may have one or two first neighbor connections in principal, but keep in mind that this routine can easely be adapted to include any number of interacting neighbors. It does not need to be the first ones. 
  SUBROUTINE RingImages
    IMPLICIT NONE
    integer :: i, j, k, NodesNumber, ResolutionTime, MeasuringTime
    real(kind=16) :: p0, CouplingStrengh, eps
    real(kind=16), dimension(:),allocatable :: ProbabilitiesTable, IntrinsicalProbabilitiesTable, ProbabilitySum 
    character(len=100) :: Arq1
    write(Arq1,"(A13)")"RingImage.dat"
    open(1,FILE=TRIM(Arq1)) 
    
    !Inicialization of the size of the ring and the time one weight the calculations of the probabilities. See the article for more comprehensive explanation. 
    NodesNumber = 200
    ResolutionTime = 1
    MeasuringTime = 1000
    
    !Initialization of the average probability, the coupling strenght and the size of the noise implicit in the initialization of the system. See how the probabilities are initialized below. Also the initialization of measurement time, do not measure for much time as it may complicate the printing. is recomended that if one wishes to know the long time behavior, just duplicate the main loop without measuring and them start another shorter one measuring. 
    p0 =  0.5
    CouplingStrengh = .51
    eps = .1
    
    Allocate(IntrinsicalProbabilitiesTable(0:NodesNumber+1))
    Allocate(ProbabilitiesTable(NodesNumber))
    Allocate(ProbabilitySum(0:NodesNumber+1))
    
    ! Initialize the intrinsic probabilities with a random noize around p0, but if you wish to make another variation you can include it in the loop below
    do i = 1,NodesNumber
      IntrinsicalProbabilitiesTable(i) = p0 + eps*(2*RandomGenerator(ISEED)-1.)
      if(IntrinsicalProbabilitiesTable(i)>1.)IntrinsicalProbabilitiesTable(i) = 1.
      if(IntrinsicalProbabilitiesTable(i)<0.)IntrinsicalProbabilitiesTable(i) = 0.
    end do
    ! periodic boundary conditions  
    IntrinsicalProbabilitiesTable(NodesNumber+1) = IntrinsicalProbabilitiesTable(1)
    IntrinsicalProbabilitiesTable(0) = IntrinsicalProbabilitiesTable(NodesNumber)
    ! p_0 = p_(t=0)
    ProbabilitiesTable = IntrinsicalProbabilitiesTable
    

    do i = 1,MeasuringTime
      ProbabilitySum = 0.0
      !Sum over the first neighbors
      do j = 1,ResolutionTime
        do k =1,NodesNumber
          if(RandomGenerator(ISEED)<=ProbabilitiesTable(k))then
            ProbabilitySum(k) = ProbabilitySum(k) + 1.
          else 
            ProbabilitySum(k) = ProbabilitySum(k) - 1.
          end if
        end do
      end do
      !boundary conditions
      ProbabilitySum(NodesNumber+1) =ProbabilitySum(1) 
      ProbabilitySum(0)=ProbabilitySum(NodesNumber)
      !Average
      ProbabilitySum = ProbabilitySum/ResolutionTime
      !probability evaluation, notice you can change the number of neighbors considered here
      do k=1,NodesNumber
        ProbabilitiesTable(k)=IntrinsicalProbabilitiesTable(k)+CouplingStrengh*&
        &CouplingFunction((ProbabilitySum(k+1)+ProbabilitySum(k-1))/2)
        if(ProbabilitiesTable(k)>1.)ProbabilitiesTable(k)=1.
        if(ProbabilitiesTable(k)<0.)ProbabilitiesTable(k)=0.
      end do
      !The results are printed in a matrix like file
      write(1,*)ProbabilitiesTable
      !It can be printed in GNUplot using the following comands
      !set pm3d map
      !set xrange [0:NodesNumber]
      !set yrange [0:MeasuringTime]
      !splot "RingImage.dat" matrix
      !Notice that you is also able to pdrocuce time series of any of the nodes using the comands
      !reset 
      !plot "RingImage.dat" u 0:NodeYouWishToSee w lp
    end do
    close(1)
    DeAllocate(IntrinsicalProbabilitiesTable)
    DeAllocate(ProbabilitiesTable)
    DeAllocate(ProbabilitySum)
  END SUBROUTINE RingImages
  !The same as the previus routine, but the geometry is a torus and not a ring 
  SUBROUTINE TorusImages
    IMPLICIT NONE
    integer :: i, j, k, l, NodesNumber, ResolutionTime, MeasuringTime
    real(kind=16) :: p0, CouplingStrengh, eps
    real(kind=16), dimension(:,:),allocatable :: ProbabilitiesTable, IntrinsicalProbabilitiesTable, ProbabilitySum 
    character(len=100) :: Arq1
    
    
    !Inicialization of the size of the torus (linear size, the overall size is NodesNumber^2) and the time one weight the calculations of the probabilities.
    NodesNumber = 200
    ResolutionTime = 1
    
    p0 =  0.5
    CouplingStrengh = .51
    eps = .1
    MeasuringTime = 100
    
    Allocate(IntrinsicalProbabilitiesTable(0:NodesNumber+1,0:NodesNumber+1))
    Allocate(ProbabilitiesTable(NodesNumber,NodesNumber))
    Allocate(ProbabilitySum(0:NodesNumber+1,0:NodesNumber+1))
    
    Do i=1,NodesNumber
      Do j=1,NodesNumber
        IntrinsicalProbabilitiesTable(i,j) = p0 + eps*(2*RandomGenerator(ISEED)-1.)
        if(IntrinsicalProbabilitiesTable(i,j)>1.)IntrinsicalProbabilitiesTable(i,j) = 1.
        if(IntrinsicalProbabilitiesTable(i,j)<0.)IntrinsicalProbabilitiesTable(i,j) = 0.
      End do  
    End do 
    !periodic contour conditions
    IntrinsicalProbabilitiesTable(NodesNumber+1,NodesNumber+1) = IntrinsicalProbabilitiesTable(1,1)
    IntrinsicalProbabilitiesTable(0,0) = IntrinsicalProbabilitiesTable(NodesNumber,NodesNumber)
    IntrinsicalProbabilitiesTable(NodesNumber+1,0) = IntrinsicalProbabilitiesTable(1,NodesNumber)
    IntrinsicalProbabilitiesTable(0,NodesNumber+1) = IntrinsicalProbabilitiesTable(NodesNumber,1)  
    ! p_0 = p_(t=0)  
    ProbabilitiesTable = IntrinsicalProbabilitiesTable
    

    do i = 1,MeasuringTime
      ProbabilitySum = 0.0
      !Sum over the first neighbors
      do j = 1,ResolutionTime
        do k =1,NodesNumber
          do l =1,NodesNumber
            if(RandomGenerator(ISEED)<=ProbabilitiesTable(k,l))then
              ProbabilitySum(k,l) = ProbabilitySum(k,l) + 1.
            else 
              ProbabilitySum(k,l) = ProbabilitySum(k,l) - 1.
            end if
          end do 
        end do
      end do
      !boundary conditions
      ProbabilitySum(NodesNumber+1,NodesNumber+1) = ProbabilitySum(1,1)
      ProbabilitySum(0,0) = ProbabilitySum(NodesNumber,NodesNumber)
      ProbabilitySum(NodesNumber+1,0) = ProbabilitySum(1,NodesNumber)
      ProbabilitySum(0,NodesNumber+1) = ProbabilitySum(NodesNumber,1)  
      !Average
      ProbabilitySum = ProbabilitySum/ResolutionTime
      !probability evaluation, notice you can change the number of neighbors considered here
      do j=1,NodesNumber 
        do k=1,NodesNumber
          ProbabilitiesTable(j,k)=IntrinsicalProbabilitiesTable(j,k)+CouplingStrengh*&
          &CouplingFunction((ProbabilitySum(j,k+1)+ProbabilitySum(j,k-1)+ProbabilitySum(j+1,k)+&
          &ProbabilitySum(j-1,k))/4)
          if(ProbabilitiesTable(j,k)>1.)ProbabilitiesTable(j,k)=1.
          if(ProbabilitiesTable(j,k)<0.)ProbabilitiesTable(j,k)=0.
        end do
      end do
      !The results are printed in files containing tree columns, the first two are the coordinates and the last th probability at the time of measement. Each file represents a snapshot of the system. In the case below, we assume that  MeasuringTime has 3 decimal numbers, change the comand below and teh GNUplot script routine that acompany this routine.
      write(Arq1,"(A10,I3.3,A4)")"TorusImage",i,".dat"
      open(1,FILE=TRIM(Arq1)) 
      Do j=1,NodesNumber
        Do k=1,NodesNumber
          write(1,*)1.0*j,1.0*k,ProbabilitiesTable(j,k)
        End Do
      ENd DO  
      close(1)
      !To build the images of this data, after calling this routine, call the routine GnuScriptTorusImages and type in the GNUplot terminal 
      !load "TorusImagesFilm.dat"
      ! This will produce MeasuringTime snapshots that you can use to built a film in here https://imgflip.com/gif-maker, for examle.
      ! The figure 03 a) in the article was made using this routine at parameters 
      !ResolutionTime = 1000
      !eps = 0.05
      !NodesNumber = 100
      !CouplingStrengh = 0.7
    end do
    DeAllocate(IntrinsicalProbabilitiesTable)
    DeAllocate(ProbabilitiesTable)
    DeAllocate(ProbabilitySum)
  END SUBROUTINE TorusImages
  ! Still in the torus geometry, it is possible to consider many more sum configurations, at least, in geometries which space can be completely filled by symmetric Wigner-Seitz cells. Here we give an example of a Wigner-Seitz cell in the shape of a equilateral triangles
  SUBROUTINE TorusTriangularCellImages
    IMPLICIT NONE
    integer :: i, j, k, l, NodesNumber, ResolutionTime, MeasuringTime, Aux1, Aux2
    real(kind=16) :: p0, CouplingStrengh, eps
    real(kind=16), dimension(:,:),allocatable :: ProbabilitiesTable, IntrinsicalProbabilitiesTable, ProbabilitySum 
    character(len=100) :: Arq1
    
    
    NodesNumber = 200
    ResolutionTime = 1
    
    p0 =  0.5
    CouplingStrengh = .51
    eps = .1
    MeasuringTime = 100
    
    Allocate(IntrinsicalProbabilitiesTable(0:NodesNumber+1,0:NodesNumber+1))
    Allocate(ProbabilitiesTable(NodesNumber,NodesNumber))
    Allocate(ProbabilitySum(0:NodesNumber+1,0:NodesNumber+1))
    
    Do i=1,NodesNumber
      Do j=1,NodesNumber
        IntrinsicalProbabilitiesTable(i,j) = p0 + eps*(2*RandomGenerator(ISEED)-1.)
        if(IntrinsicalProbabilitiesTable(i,j)>1.)IntrinsicalProbabilitiesTable(i,j) = 1.
        if(IntrinsicalProbabilitiesTable(i,j)<0.)IntrinsicalProbabilitiesTable(i,j) = 0.
      End do  
    End do 
    !periodic contour conditions
    IntrinsicalProbabilitiesTable(NodesNumber+1,NodesNumber+1) = IntrinsicalProbabilitiesTable(1,1)
    IntrinsicalProbabilitiesTable(0,0) = IntrinsicalProbabilitiesTable(NodesNumber,NodesNumber)
    IntrinsicalProbabilitiesTable(NodesNumber+1,0) = IntrinsicalProbabilitiesTable(1,NodesNumber)
    IntrinsicalProbabilitiesTable(0,NodesNumber+1) = IntrinsicalProbabilitiesTable(NodesNumber,1)  
    ! p_0 = p_(t=0)  
    ProbabilitiesTable = IntrinsicalProbabilitiesTable
    

    do i = 1,MeasuringTime
      ProbabilitySum = 0.0
      !Sum over the first neighbors
      do j = 1,ResolutionTime
        do k =1,NodesNumber
          do l =1,NodesNumber
            if(RandomGenerator(ISEED)<=ProbabilitiesTable(k,l))then
              ProbabilitySum(k,l) = ProbabilitySum(k,l) + 1.
            else 
              ProbabilitySum(k,l) = ProbabilitySum(k,l) - 1.
            end if
          end do 
        end do
      end do
      !boundary conditions
      ProbabilitySum(NodesNumber+1,NodesNumber+1) = ProbabilitySum(1,1)
      ProbabilitySum(0,0) = ProbabilitySum(NodesNumber,NodesNumber)
      ProbabilitySum(NodesNumber+1,0) = ProbabilitySum(1,NodesNumber)
      ProbabilitySum(0,NodesNumber+1) = ProbabilitySum(NodesNumber,1)  
      !Average
      ProbabilitySum = ProbabilitySum/ResolutionTime
      !probability evaluation, here we must have a logical clause to sum the right neighbors. In a triangular lattice, consider that we are in a coordinated system such that we sum both horizontal neighbors, we sum the upper neighbor to complete tre, but sum the lower neighbor in the next site to complete the tree. Is this way is easy to see that this construction completes the space
      do j=1,NodesNumber
        Aux1 = MOD(j,2)
        SELECT CASE (Aux1)
        CASE(1) ! Odd number 
          do k=1,NodesNumber
            Aux2 = MOD(k,2)
            SELECT CASE (Aux2)
            CASE(0)
              ProbabilitiesTable(j,k)=IntrinsicalProbabilitiesTable(j,k)+CouplingStrengh*&
              &CouplingFunction((ProbabilitySum(j,k+1)+ProbabilitySum(j,k-1)+ProbabilitySum(j+1,k))/3)
            CASE(1)  
              ProbabilitiesTable(j,k)=IntrinsicalProbabilitiesTable(j,k)+CouplingStrengh*&
              &CouplingFunction((ProbabilitySum(j,k+1)+ProbabilitySum(j,k-1)+ProbabilitySum(j-1,k))/3)
            END SELECT  
            if(ProbabilitiesTable(j,k)>1.)ProbabilitiesTable(j,k)=1.
            if(ProbabilitiesTable(j,k)<0.)ProbabilitiesTable(j,k)=0.
          end do
        CASE(0)! Even number 
          do k=1,NodesNumber
            Aux2 = MOD(k,2)
            SELECT CASE (Aux2)
            CASE(0)
              ProbabilitiesTable(j,k)=IntrinsicalProbabilitiesTable(j,k)+CouplingStrengh*&
              &CouplingFunction((ProbabilitySum(j,k+1)+ProbabilitySum(j,k-1)+ProbabilitySum(j-1,k))/3)
            CASE(1)  
              ProbabilitiesTable(j,k)=IntrinsicalProbabilitiesTable(j,k)+CouplingStrengh*&
              &CouplingFunction((ProbabilitySum(j,k+1)+ProbabilitySum(j,k-1)+ProbabilitySum(j+1,k))/3)
            END SELECT  
            if(ProbabilitiesTable(j,k)>1.)ProbabilitiesTable(j,k)=1.
            if(ProbabilitiesTable(j,k)<0.)ProbabilitiesTable(j,k)=0.
          end do
        END SELECT
      end do
      
      
    
      !The results are printed in files containing tree columns, the first two are the coordinates and the last th probability at the time of measement. Each file represents a snapshot of the system. In the case below, we assume that  MeasuringTime has 3 decimal numbers, change the comand below and teh GNUplot script routine that acompany this routine.
      
      !To plot it in the right geometry, more logical loops are necessary
      write(Arq1,"(A24,I3.3,A4)")"TorusTriangularCellImage",i,".dat"
      open(1,FILE=TRIM(Arq1)) 
      Do j=1,NodesNumber
        Aux1 = MOD(j,1)
        SELECT CASE (Aux1)
        CASE(0)
          Do k=1,NodesNumber
            Aux2 = MOD(j,1)
            SELECT CASE (Aux2)
            CASE(0)
              write(1,*)1.0*j + 4.0/81 - SQRT(3.)/72, k + SQRT(3.)/72,ProbabilitiesTable(j,k)
            CASE(1)
              write(1,*)1.0*j + 4.0/81 - SQRT(3.)/72, k - SQRT(3.)/72,ProbabilitiesTable(j,k)
            END SELECT
          End do
        CASE(1)
          Do k=1,NodesNumber
            Aux2 = MOD(j,1)
            SELECT CASE (Aux2)
            CASE(0)
              write(1,*)1.0*j + 4.0/81 - SQRT(3.)/72,k - SQRT(3.)/72,ProbabilitiesTable(j,k) 
            CASE(1)
              write(1,*)1.0*j + 4.0/81 - SQRT(3.)/72,k + SQRT(3.)/72,ProbabilitiesTable(j,k) 
            END SELECT
          End do
        END SELECT
      end do
      close(1)
      !To build the images of this data, after calling this routine, call the routine GnuScriptTorusImages and type in the GNUplot terminal 
      !load "TorusImagesTriangularCellFilm.dat"
      ! This will produce MeasuringTime snapshots that you can use to built a film in here https://imgflip.com/gif-maker, for examle.
      
      ! The figure 03 b) in the article was made using this routine at parameters 
      !ResolutionTime = 1000
      !eps = 0.05
      !NodesNumber = 100
      !CouplingStrengh = 0.7
    end do
    DeAllocate(IntrinsicalProbabilitiesTable)
    DeAllocate(ProbabilitiesTable)
    DeAllocate(ProbabilitySum)
  END SUBROUTINE TorusTriangularCellImages
  
  ! Now we turn to more quantitative measurements starting from the phases diagram
  
  ! This routine prints the Kxp_m of a population in a ring geometry 
  SUBROUTINE RingPhasesDiagram
    IMPLICIT NONE
    integer :: i, j, k, l, m, NodesNumber, ResolutionTime, FinalInstant
    real(kind=16) :: p0, CouplingStrengh, eps, Aux
    real(kind=16), dimension(:),allocatable :: ProbabilitiesTable, IntrinsicalProbabilitiesTable, ProbabilitySum 
    character(len=100) :: Arq1
    write(Arq1,"(A21)")"RingPhasesDiagram.dat"
    open(1,FILE=TRIM(Arq1)) 
    
    !Inicialization of the size of the ring, the time one weight the calculations of the probabilities and the instant the measurements start 
    NodesNumber = 100
    ResolutionTime = 1
    FinalInstant = 1000
    eps = .01
    Allocate(IntrinsicalProbabilitiesTable(0:NodesNumber+1))
    Allocate(ProbabilitiesTable(NodesNumber))
    Allocate(ProbabilitySum(0:NodesNumber+1))
    
    ! The p0 and the CouplingStrengh are variable parameters
    Do l = 1,100
      p0 =  0.01*l
      ! Initialize the intrinsic probabilities with a random noize around p0, but if you wish to make another variation you can include it in the loop below
      do i = 1,NodesNumber
        IntrinsicalProbabilitiesTable(i) = p0 + eps*(2*RandomGenerator(ISEED)-1.)
        if(IntrinsicalProbabilitiesTable(i)>1.)IntrinsicalProbabilitiesTable(i) = 1.
        if(IntrinsicalProbabilitiesTable(i)<0.)IntrinsicalProbabilitiesTable(i) = 0.
      end do
      ! periodic boundary conditions  
      IntrinsicalProbabilitiesTable(NodesNumber+1) = IntrinsicalProbabilitiesTable(1)
      IntrinsicalProbabilitiesTable(0) = IntrinsicalProbabilitiesTable(NodesNumber)
      ! p_0 = p_(t=0)
      ProbabilitiesTable = IntrinsicalProbabilitiesTable
      Do m = 1,100
        CouplingStrengh = 0.01*m
        do i = 1,FinalInstant
          ProbabilitySum = 0.0
          !Sum over the first neighbors
          do j = 1,ResolutionTime
            do k =1,NodesNumber
              if(RandomGenerator(ISEED)<=ProbabilitiesTable(k))then
                ProbabilitySum(k) = ProbabilitySum(k) + 1.
              else 
                ProbabilitySum(k) = ProbabilitySum(k) - 1.
              end if
            end do
          end do
          !boundary conditions
          ProbabilitySum(NodesNumber+1) =ProbabilitySum(1) 
          ProbabilitySum(0)=ProbabilitySum(NodesNumber)
          !Average
          ProbabilitySum = ProbabilitySum/ResolutionTime
          !probability evaluation, notice you can change the number of neighbors considered here
          do k=1,NodesNumber
            ! To alter the probability dependency, one must alter the CouplingFunction AND the line below
            ProbabilitiesTable(k)=IntrinsicalProbabilitiesTable(k)+CouplingStrengh*&
            &CouplingFunction((ProbabilitySum(k+1)+ProbabilitySum(k-1))/2)
            if(ProbabilitiesTable(k)>1.)ProbabilitiesTable(k)=1.
            if(ProbabilitiesTable(k)<0.)ProbabilitiesTable(k)=0.
          end do
        end do
        ! At the instant FinalInstant we measure the system
        Aux = 0.0
        Do k=1,NodesNumber
          Aux=Aux+ProbabilitiesTable(k) 
        End do
        Aux=Aux/NodesNumber
        write(1,*)CouplingStrengh,p0,Aux
      end do
    end do  
    !The figures 01 a) and 01 c) were taken using this routine with eps = 0, FinalInstant = 10.000 and NodesNumber = 100
    close(1)
    DeAllocate(IntrinsicalProbabilitiesTable)
    DeAllocate(ProbabilitiesTable)
    DeAllocate(ProbabilitySum)
  END SUBROUTINE RingPhasesDiagram
  !The same as above, but in a Torus geometry
  SUBROUTINE TorusPhasesDiagram
    IMPLICIT NONE
    integer :: i, j, k, l, m, n, NodesNumber, ResolutionTime, FinalInstant
    real(kind=16) :: p0, CouplingStrengh, eps, Aux
    real(kind=16), dimension(:,:),allocatable :: ProbabilitiesTable, IntrinsicalProbabilitiesTable, ProbabilitySum 
    character(len=100) :: Arq1
    write(Arq1,"(A22)")"TorusPhasesDiagram.dat"
    open(1,FILE=TRIM(Arq1)) 
    
    NodesNumber = 10
    ResolutionTime = 1
    eps = .01
    FinalInstant = 1000
    Allocate(IntrinsicalProbabilitiesTable(0:NodesNumber+1,0:NodesNumber+1))
    Allocate(ProbabilitiesTable(NodesNumber,NodesNumber))
    Allocate(ProbabilitySum(0:NodesNumber+1,0:NodesNumber+1))
    Do m = 0,100
      p0 = 0.01*m
      Do i=1,NodesNumber
        Do j=1,NodesNumber
          IntrinsicalProbabilitiesTable(i,j) = p0 + eps*(2*RandomGenerator(ISEED)-1.)
          if(IntrinsicalProbabilitiesTable(i,j)>1.)IntrinsicalProbabilitiesTable(i,j) = 1.
          if(IntrinsicalProbabilitiesTable(i,j)<0.)IntrinsicalProbabilitiesTable(i,j) = 0.
        End do  
      End do 
      IntrinsicalProbabilitiesTable(NodesNumber+1,NodesNumber+1) = IntrinsicalProbabilitiesTable(1,1)
      IntrinsicalProbabilitiesTable(0,0) = IntrinsicalProbabilitiesTable(NodesNumber,NodesNumber)
      IntrinsicalProbabilitiesTable(NodesNumber+1,0) = IntrinsicalProbabilitiesTable(1,NodesNumber)
      IntrinsicalProbabilitiesTable(0,NodesNumber+1) = IntrinsicalProbabilitiesTable(NodesNumber,1)  
      ProbabilitiesTable = IntrinsicalProbabilitiesTable
      Do n = 0,100
        CouplingStrengh = 0.01*n
        do i = 1,FinalInstant
          ProbabilitySum = 0.0
          do j = 1,ResolutionTime
            do k =1,NodesNumber
              do l =1,NodesNumber
                if(RandomGenerator(ISEED)<=ProbabilitiesTable(k,l))then
                  ProbabilitySum(k,l) = ProbabilitySum(k,l) + 1.
                else 
                  ProbabilitySum(k,l) = ProbabilitySum(k,l) - 1.
                end if
              end do 
            end do
          end do
          ProbabilitySum(NodesNumber+1,NodesNumber+1) = ProbabilitySum(1,1)
          ProbabilitySum(0,0) = ProbabilitySum(NodesNumber,NodesNumber)
          ProbabilitySum(NodesNumber+1,0) = ProbabilitySum(1,NodesNumber)
          ProbabilitySum(0,NodesNumber+1) = ProbabilitySum(NodesNumber,1)  
          ProbabilitySum = ProbabilitySum/ResolutionTime
          do j=1,NodesNumber 
            do k=1,NodesNumber
              ProbabilitiesTable(j,k)=IntrinsicalProbabilitiesTable(j,k)+CouplingStrengh*&
              &CouplingFunction((ProbabilitySum(j,k+1)+ProbabilitySum(j,k-1)+ProbabilitySum(j+1,k)+&
              &ProbabilitySum(j-1,k))/4)
              if(ProbabilitiesTable(j,k)>1.)ProbabilitiesTable(j,k)=1.
              if(ProbabilitiesTable(j,k)<0.)ProbabilitiesTable(j,k)=0.
            end do
          end do
        end do
        Aux = 0.0
        Do j=1,NodesNumber
          Do k=1,NodesNumber
            Aux=Aux+ProbabilitiesTable(j,k) 
          End do
        End do
        Aux=Aux/(NodesNumber*NodesNumber)
        write(1,*)CouplingStrengh,p0,Aux
      End do   
    End do
    close(1)
    DeAllocate(IntrinsicalProbabilitiesTable)
    DeAllocate(ProbabilitiesTable)
    DeAllocate(ProbabilitySum)
  END SUBROUTINE TorusPhasesDiagram
  
  
  !We direct our attention now to the perturbation analysis
  !The perturbation of the second type, heterogeneous intrinsic probabilities, were already considered in the  RingImages routine. Now we construct the first type, the Dirac perturbation on the central node 
  SUBROUTINE RingDiracPerturbation
    IMPLICIT NONE
    integer :: i, j, k, NodesNumber, ResolutionTime, MeasuringTime
    real(kind=16) :: p0, CouplingStrengh, eps
    real(kind=16), dimension(:),allocatable :: ProbabilitiesTable, IntrinsicalProbabilitiesTable, ProbabilitySum 
    character(len=100) :: Arq1
    write(Arq1,"(A25)")"RingDiracPerturbation.dat"
    open(1,FILE=TRIM(Arq1)) 
    
    !Inicialization of the size of the ring and the time one weight the calculations of the probabilities. See the article for more comprehensive explanation. 
    NodesNumber = 100
    ResolutionTime = 1000
    MeasuringTime = 1000
    
    !Initialization of the average probability, the coupling strenght and the size of the noise implicit in the initialization of the system. See how the probabilities are initialized below. Also the initialization of measurement time, do not measure for much time as it may complicate the printing. is recomended that if one wishes to know the long time behavior, just duplicate the main loop without measuring and them start another shorter one measuring. 
    p0 =  0.5
    CouplingStrengh = .51
    eps = .0
    
    Allocate(IntrinsicalProbabilitiesTable(0:NodesNumber+1))
    Allocate(ProbabilitiesTable(NodesNumber))
    Allocate(ProbabilitySum(0:NodesNumber+1))
    
    ! Initialize the intrinsic probabilities with a random noize around p0, but if you wish to make another variation you can include it in the loop below
    do i = 1,NodesNumber
      IntrinsicalProbabilitiesTable(i) = p0 + eps*(2*RandomGenerator(ISEED)-1.)
      if(IntrinsicalProbabilitiesTable(i)>1.)IntrinsicalProbabilitiesTable(i) = 1.
      if(IntrinsicalProbabilitiesTable(i)<0.)IntrinsicalProbabilitiesTable(i) = 0.
    end do
    ! periodic boundary conditions  
    IntrinsicalProbabilitiesTable(NodesNumber+1) = IntrinsicalProbabilitiesTable(1)
    IntrinsicalProbabilitiesTable(0) = IntrinsicalProbabilitiesTable(NodesNumber)
    ! p_0 = p_(t=0)
    ProbabilitiesTable = IntrinsicalProbabilitiesTable
    

    do i = 1,MeasuringTime
      ProbabilitySum = 0.0
      !Sum over the first neighbors
      do j = 1,ResolutionTime
        do k =1,NodesNumber
          if(RandomGenerator(ISEED)<=ProbabilitiesTable(k))then
            ProbabilitySum(k) = ProbabilitySum(k) + 1.
          else 
            ProbabilitySum(k) = ProbabilitySum(k) - 1.
          end if
        end do
      end do
      !boundary conditions
      ProbabilitySum(NodesNumber+1) =ProbabilitySum(1) 
      ProbabilitySum(0)=ProbabilitySum(NodesNumber)
      !Average
      ProbabilitySum = ProbabilitySum/ResolutionTime
      !probability evaluation, notice you can change the number of neighbors considered here
      do k=1,NodesNumber
        ProbabilitiesTable(k)=IntrinsicalProbabilitiesTable(k)+CouplingStrengh*&
        &CouplingFunction((ProbabilitySum(k+1)+ProbabilitySum(k-1))/2)
        if(ProbabilitiesTable(k)>1.)ProbabilitiesTable(k)=1.
        if(ProbabilitiesTable(k)<0.)ProbabilitiesTable(k)=0.
      end do
      ! We introduce the perturbation at the instant t with a logical clause
      If(i==100)then
        ProbabilitiesTable(NodesNumber/2)=ProbabilitiesTable(NodesNumber/2)+eps
        if(ProbabilitiesTable(NodesNumber/2)>1.)ProbabilitiesTable(NodesNumber/2)=1.
        if(ProbabilitiesTable(NodesNumber/2)<0.)ProbabilitiesTable(NodesNumber/2)=0.
      End if
      !The results are printed in a matrix like file
      write(1,*)ProbabilitiesTable
      !The time series of any of interest can be printed in GNUplot typing 
      !plot "RingDiracPerturbation.dat" u 0:NodeYouWishToSee 
    end do
    close(1)
    DeAllocate(IntrinsicalProbabilitiesTable)
    DeAllocate(ProbabilitiesTable)
    DeAllocate(ProbabilitySum)
  END SUBROUTINE RingDiracPerturbation
  ! The same as above but in a toroidal geometry
  SUBROUTINE TorusDiracPerturbation
    IMPLICIT NONE
    integer :: i, j, k, l, NodesNumber, ResolutionTime, MeasuringTime
    real(kind=16) :: p0, CouplingStrengh, eps
    real(kind=16), dimension(:,:),allocatable :: ProbabilitiesTable, IntrinsicalProbabilitiesTable, ProbabilitySum 
    character(len=100) :: Arq1
    write(Arq1,"(A26)")"TorusDiracPerturbation.dat"
    open(1,FILE=TRIM(Arq1)) 
    
    !Inicialization of the size of the torus (linear size, the overall size is NodesNumber^2) and the time one weight the calculations of the probabilities.
    NodesNumber = 50
    ResolutionTime = 1
    eps = .0
    MeasuringTime = 100
    p0 = 0.5
    CouplingStrengh = 0.51
    Allocate(IntrinsicalProbabilitiesTable(0:NodesNumber+1,0:NodesNumber+1))
    Allocate(ProbabilitiesTable(NodesNumber,NodesNumber))
    Allocate(ProbabilitySum(0:NodesNumber+1,0:NodesNumber+1))
    
    Do i=1,NodesNumber
      Do j=1,NodesNumber
        IntrinsicalProbabilitiesTable(i,j) = p0 + eps*(2*RandomGenerator(ISEED)-1.)
        if(IntrinsicalProbabilitiesTable(i,j)>1.)IntrinsicalProbabilitiesTable(i,j) = 1.
        if(IntrinsicalProbabilitiesTable(i,j)<0.)IntrinsicalProbabilitiesTable(i,j) = 0.
      End do  
    End do 
    IntrinsicalProbabilitiesTable(NodesNumber+1,NodesNumber+1) = IntrinsicalProbabilitiesTable(1,1)
    IntrinsicalProbabilitiesTable(0,0) = IntrinsicalProbabilitiesTable(NodesNumber,NodesNumber)
    IntrinsicalProbabilitiesTable(NodesNumber+1,0) = IntrinsicalProbabilitiesTable(1,NodesNumber)
    IntrinsicalProbabilitiesTable(0,NodesNumber+1) = IntrinsicalProbabilitiesTable(NodesNumber,1)  
    ProbabilitiesTable = IntrinsicalProbabilitiesTable
    Do i = 1,MeasuringTime
      ProbabilitySum = 0.0
      Do j = 1,ResolutionTime
        Do k =1,NodesNumber
          Do l =1,NodesNumber
            If(RandomGenerator(ISEED)<=ProbabilitiesTable(k,l))then
              ProbabilitySum(k,l) = ProbabilitySum(k,l) + 1.
            Else 
              ProbabilitySum(k,l) = ProbabilitySum(k,l) - 1.
            End if
          End do 
        End do
      End do
      ProbabilitySum(NodesNumber+1,NodesNumber+1) = ProbabilitySum(1,1)
      ProbabilitySum(0,0) = ProbabilitySum(NodesNumber,NodesNumber)
      ProbabilitySum(NodesNumber+1,0) = ProbabilitySum(1,NodesNumber)
      ProbabilitySum(0,NodesNumber+1) = ProbabilitySum(NodesNumber,1)  
      ProbabilitySum = ProbabilitySum/ResolutionTime
      Do j=1,NodesNumber 
        Do k=1,NodesNumber
          ProbabilitiesTable(j,k)=IntrinsicalProbabilitiesTable(j,k)+CouplingStrengh*&
          &CouplingFunction((ProbabilitySum(j,k+1)+ProbabilitySum(j,k-1)+ProbabilitySum(j+1,k)+&
          &ProbabilitySum(j-1,k))/4)
          If(ProbabilitiesTable(j,k)>1.)ProbabilitiesTable(j,k)=1.
          If(ProbabilitiesTable(j,k)<0.)ProbabilitiesTable(j,k)=0.
        End do
      End do
      If(i==100)then
        ProbabilitiesTable(NodesNumber/2,NodesNumber/2)=ProbabilitiesTable(NodesNumber/2,NodesNumber/2)+eps
        If(ProbabilitiesTable(NodesNumber/2,NodesNumber/2)>1.)ProbabilitiesTable(NodesNumber/2,NodesNumber/2)=1.
        If(ProbabilitiesTable(NodesNumber/2,NodesNumber/2)<0.)ProbabilitiesTable(NodesNumber/2,NodesNumber/2)=0.
      End if
      Do j=1,NodesNumber
        Do k=1,NodesNumber
          write(1,"(f9.2,A1)",ADVANCE='NO')ProbabilitiesTable(j,k)," " 
        End Do
      End do
      write(1,*)" "
      !The results of ALL the nodes are pinted in the same line so that the node of interes Node(i,j) can be printed typing 
      !!plot "TorusDiracPerturbation.dat" u 0:(NodesNumber*j+k)
    End do
    close(1)
    DeAllocate(IntrinsicalProbabilitiesTable)
    DeAllocate(ProbabilitiesTable)
    DeAllocate(ProbabilitySum)
  END SUBROUTINE TorusDiracPerturbation

  
  
END MODULE MainRoutines





! This module present simple tests of the routines implemented so far

MODULE Tests
  USE GlobalVariables
  USE SecondaryRoutines
  USE MainRoutines
  IMPLICIT NONE
  CONTAINS
  
  
  ! Test of the random generator using histograms
  ! This routine divides the frequency of the random numbers in AmostralSize discreet steps and 
  ! measure their frequency over ExperimentTime time steps. It presents the results in a histogram 
  ! saved in the file RandoGeneratorHistogram.dat
  SUBROUTINE RandomGeneratorTest
    IMPLICIT NONE
    integer :: i, indice, AmostralSize, ExperimentTime
    real(kind=16), dimension(:),allocatable :: Histogram
    character(len=100) :: Arq1
    write(Arq1,"(A27)")"RandoGeneratorHistogram.dat"
    open(1,FILE=TRIM(Arq1)) 
    
    AmostralSize = 10000
    ExperimentTime = 100000*AmostralSize
    Allocate(Histogram(AmostralSize))
    Histogram = 0.0
    
    Do i = 1,ExperimentTime
      indice = INT(AmostralSize*RandomGenerator(ISEED))+1
      Histogram(indice)=Histogram(indice)+1.0
    End Do 
    Histogram = Histogram/ExperimentTime
    Do i = 1,AmostralSize
      write(1,*)i,Histogram(i)
    End do  
    close(1)
    DeAllocate(Histogram)
    !You can verify the results with Gnuplot typing in the console: plot "RandoGeneratorHistogram.dat" w boxes
  END SUBROUTINE RandomGeneratorTest

  
END MODULE Tests


PROGRAM MAIN
  USE GlobalVariables
  USE SecondaryRoutines
  USE MainRoutines
  USE Tests
  IMPLICIT NONE
  Iseed =1313
  
  !Call RandomGeneratorTest
  !Call RingImages
  
  !Call TorusImages
  !Call GnuScriptTorusImages
  
  !Call TorusTriangularCellImages
  !Call GnuScriptTorusTriangularCellImages
  
  !Call RingPhasesDiagram
  
  !Call TorusPhasesDiagram
 
  !Call RingDiracPerturbation
  
  Call TorusDiracPerturbation
END PROGRAM MAIN
