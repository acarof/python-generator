program analytic_rabi_oscillation

   INTEGER, PARAMETER :: dp=KIND(1.0D0)
   REAL(KIND=dp)      :: amplitude, frequence, atomic_unit_time, DeltaT
   REAL(KIND=dp)      :: Hab, Ea, Eb
   INTEGER            :: number_steps
 
   atomic_unit_time = 2.418884326505E-17/1.0E-15

   print*, "Write: Hab, Eb, Ea, DeltaT (fs), number_steps"
   read(*,*) Hab, Eb, Ea, DeltaT, number_steps

   DeltaT = DeltaT/atomic_unit_time
   amplitude = (4*Hab**2)/ ( (Eb-Ea)**2 + 4*Hab**2)
   frequence = SQRT( (Eb-Ea)**2 + 4*Hab**2 )

   print*, "amplitude =", amplitude
   print*, "frequence =", frequence
 
   OPEN(unit=101, file="analytic_rabi_oscillation.dat")
   DO i=1,number_steps
      write(101,*)  amplitude*( SIN( frequence*(i-1)*DeltaT/2) )**2
   ENDDO
   CLOSE(101)
      
END 
