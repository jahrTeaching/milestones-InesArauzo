real, target :: a
real, pointer :: pa

a = 4
pa=> a

a=3

write(*,*) "a=", a, "pa=", pa



################################################

real, target :: a(3)
real, pointer :: pa(:) #es un pointer a vectores, porque hasta que no apunte algo no se el tamaño

a = [1, 2, 3] 
pa=> a

a=a+1

write(*,*) "a=", a, "pa=", pa# sale pa=a=[2, 3, 4] 

####################################################

real, target :: a(4)
real, pointer :: pa(:, :) #es un pointer a vectores, porque hasta que no apunte algo no se el tamaño

a = [1, 2, 3, 4] 
pa(1:2 , 1:2)=> a # te mete en la primera col de la matriz, 1, 2 y en la segunda COL 3, 4. OJO QUE PYTHON Y C LO HACEN POR FILAS

write(*,*) "a=", a, "pa=", pa#  
a= 0
